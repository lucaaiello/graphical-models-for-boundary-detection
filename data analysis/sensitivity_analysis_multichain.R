rm(list = ls())

# SEER sensitivity analysis from adaptive-tempering chains
# =============================================================================
# This script does not run MCMC. It compares the completed six-chain
# Adjacency, Mean--Adjacency, and Mean analyses. Only the cold posterior chains
# from the adaptive-tempering workflow are accepted as input.
#
# Run from the repository root with:
#   source("data analysis/sensitivity_analysis_multichain.R")
#
# Required inputs for each of adj, meanadj, and mean:
#   data analysis/multiple_chains_tempering_output/
#     <specification>/hotter_final/ladders/ladder_[1-6]_cold_chain.rds
#   data analysis/multiple_chains_tempering_output/
#     <specification>/hotter_final/postprocessing/
#       tempering_fdr_boundary_results.rds
#
# All generated artifacts are written below data analysis/sensitivity/.

options(stringsAsFactors = FALSE)

input_root <- Sys.getenv(
  "SEER_MULTICHAIN_INPUT_ROOT",
  unset = file.path(
    "data analysis", "multiple_chains_tempering_output"
  )
)
run_mode <- Sys.getenv("SEER_MULTICHAIN_RUN_MODE", unset = "hotter_final")
output_dir <- Sys.getenv(
  "SEER_MULTICHAIN_SENSITIVITY_OUTPUT",
  unset = file.path("data analysis", "sensitivity")
)
figures_dir <- file.path(output_dir, "figures")
tables_dir <- file.path(output_dir, "tables")
diagnostics_dir <- file.path(output_dir, "diagnostics")
for (directory in c(output_dir, figures_dir, tables_dir, diagnostics_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

required_packages <- c("Matrix", "ggplot2", "sf")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Run install.R first.", call. = FALSE
  )
}

specifications <- list(
  adj = list(label = "Adjacency", order = 1L),
  meanadj = list(label = "Mean--Adjacency", order = 2L),
  mean = list(label = "Mean", order = 3L)
)
disease_names <- c("Lung", "Esophageal", "Larynx", "Colorectal")
n_chains <- 6L
fdr_alpha <- 0.05

stopf <- function(...) stop(sprintf(...), call. = FALSE)
title_case_county <- function(x) tools::toTitleCase(as.character(x))

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_text <- function(x, path) {
  writeLines(x, path, useBytes = TRUE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

md5_file <- function(path) unname(tools::md5sum(path))

# Reconstruct the audited data ordering and California geometry. Evaluation is
# deliberately stopped before any sampler is compiled or run.
load_seer_setup <- function() {
  driver <- file.path("data analysis", "main_multiple_chains.R")
  if (!file.exists(driver)) stopf("Missing model-setup driver: %s", driver)
  source_lines <- readLines(driver, warn = FALSE)
  marker <- grep(
    "^# Baseline multiple-chain MCMC workflow", source_lines
  )
  if (length(marker) != 1L || marker <= 2L) {
    stop("Could not locate the model-setup boundary in the model driver.",
         call. = FALSE)
  }
  environment <- new.env(parent = .GlobalEnv)
  suppressPackageStartupMessages(eval(
    parse(text = source_lines[2L:(marker - 1L)]), envir = environment
  ))
  environment
}

build_design_matrix <- function(setup, specification) {
  X1 <- get("X1", setup)
  X2 <- get("X2", setup)
  X3 <- get("X3", setup)
  X4 <- get("X4", setup)
  if (specification %in% c("meanadj", "mean")) {
    as.matrix(Matrix::bdiag(
      Matrix::bdiag(X1[, c(1, 2, 4, 6)], X2[, c(1, 2, 4, 6)]),
      Matrix::bdiag(X3[, c(1, 2, 4, 6)], X4[, c(1, 2, 4, 6)])
    ))
  } else {
    as.matrix(Matrix::bdiag(
      Matrix::bdiag(X1[, 1], X2[, 1]),
      Matrix::bdiag(X3[, 1], X4[, 1])
    ))
  }
}

setup <- load_seer_setup()
reference <- list(
  final_perm = as.integer(get("final_perm", setup)),
  edges = unname(as.matrix(get("neighbor_list0", setup))),
  edge_names = as.character(get("neighbor_name", setup)),
  county_names = as.character(get("county.ID", setup)),
  county_map = get("ca.poly", setup),
  Y = as.numeric(get("Y", setup)),
  E = as.numeric(get("E", setup))
)
design_matrices <- lapply(names(specifications), function(specification) {
  build_design_matrix(setup, specification)
})
names(design_matrices) <- names(specifications)
rm(setup)
invisible(gc(FALSE))

n_counties <- length(reference$county_names)
n_diseases <- length(disease_names)
n_edges <- nrow(reference$edges)
if (
  n_counties != 58L || n_diseases != 4L || n_edges != 139L ||
    length(reference$Y) != n_counties * n_diseases ||
    length(reference$E) != n_counties * n_diseases ||
    !identical(sort(reference$final_perm), seq_len(n_counties))
) {
  stop("The reconstructed SEER data dimensions/order are not as expected.",
       call. = FALSE)
}

reorder_vector_to_geographic <- function(x) {
  inverse_permutation <- order(reference$final_perm)
  unlist(lapply(seq_len(n_diseases), function(disease) {
    indices <- ((disease - 1L) * n_counties + 1L):(disease * n_counties)
    x[indices][inverse_permutation]
  }), use.names = FALSE)
}

validate_chain <- function(fit, specification, chain_file) {
  required <- c("beta", "phi", "settings", "seed", "runtime")
  missing <- setdiff(required, names(fit))
  if (length(missing)) {
    stopf("%s is missing: %s", chain_file, paste(missing, collapse = ", "))
  }
  saved_specification <- fit$settings$model_specification
  if (is.null(saved_specification)) saved_specification <- "adj"
  if (!identical(saved_specification, specification)) {
    stopf(
      "%s contains specification '%s', not '%s'.",
      chain_file, saved_specification, specification
    )
  }
  if (
    !identical(fit$settings$run_mode, run_mode) ||
      is.null(fit$settings$sampler) ||
      !grepl("tempering", fit$settings$sampler, fixed = TRUE) ||
      is.null(fit$settings$n_temperatures) ||
      fit$settings$n_temperatures < 2L
  ) {
    stopf(
      "%s is not a completed adaptive-tempering '%s' chain.",
      chain_file, run_mode
    )
  }
  X <- design_matrices[[specification]]
  if (
    !is.matrix(fit$beta) || !is.matrix(fit$phi) ||
      nrow(fit$beta) < 2L || nrow(fit$phi) != nrow(fit$beta) ||
      ncol(fit$beta) != ncol(X) ||
      ncol(fit$phi) != n_counties * n_diseases
  ) {
    stopf("%s has incompatible beta/phi draws.", chain_file)
  }
  if (any(!is.finite(fit$beta)) || any(!is.finite(fit$phi))) {
    stopf("%s contains non-finite beta/phi draws.", chain_file)
  }
  invisible(TRUE)
}

load_specification <- function(specification) {
  input_dir <- file.path(input_root, specification, run_mode)
  chain_files <- file.path(
    input_dir, "ladders",
    sprintf("ladder_%d_cold_chain.rds", seq_len(n_chains))
  )
  boundary_file <- file.path(
    input_dir, "postprocessing",
    "tempering_fdr_boundary_results.rds"
  )
  missing <- c(chain_files, boundary_file)[
    !file.exists(c(chain_files, boundary_file))
  ]
  if (length(missing)) {
    stop(
      "The completed multiple-chain input for '", specification,
      "' is incomplete. Run its hotter_final analysis and post-processing.\n",
      paste(" -", missing, collapse = "\n"), call. = FALSE
    )
  }

  boundary <- readRDS(boundary_file)
  if (!identical(boundary$model_specification, specification)) {
    stopf("Boundary summary %s has the wrong specification.", boundary_file)
  }
  if (!isTRUE(all.equal(boundary$fdr_alpha, fdr_alpha))) {
    stopf("Boundary summary %s does not use FDR alpha %.3f.",
          boundary_file, fdr_alpha)
  }
  probabilities <- boundary$pooled_boundary_probabilities
  if (!is.matrix(probabilities) ||
      !all(dim(probabilities) == c(n_edges, n_diseases))) {
    stopf("Boundary probabilities in %s have incompatible dimensions.",
          boundary_file)
  }
  if (!is.null(rownames(probabilities))) {
    if (!setequal(rownames(probabilities), reference$edge_names)) {
      stopf("Boundary names in %s do not match the SEER graph.", boundary_file)
    }
    probabilities <- probabilities[reference$edge_names, , drop = FALSE]
  }
  if (!is.null(colnames(probabilities))) {
    if (!setequal(colnames(probabilities), disease_names)) {
      stopf("Disease names in %s are incompatible.", boundary_file)
    }
    probabilities <- probabilities[, disease_names, drop = FALSE]
  }
  dimnames(probabilities) <- list(reference$edge_names, disease_names)

  selected_sets <- boundary$pooled_selected_boundary_sets
  if (!is.list(selected_sets) || !setequal(names(selected_sets), disease_names)) {
    stopf("Selected boundary sets in %s are incomplete.", boundary_file)
  }
  selections <- vapply(disease_names, function(disease) {
    reference$edge_names %in% selected_sets[[disease]]
  }, logical(n_edges))
  dimnames(selections) <- list(reference$edge_names, disease_names)

  X <- design_matrices[[specification]]
  fitted_risk_sum <- numeric(n_counties * n_diseases)
  total_draws <- 0L
  metadata <- vector("list", n_chains)
  for (chain in seq_len(n_chains)) {
    fit <- readRDS(chain_files[[chain]])
    validate_chain(fit, specification, chain_files[[chain]])
    relative_risk <- exp(fit$beta %*% t(X) + fit$phi)
    if (any(!is.finite(relative_risk))) {
      stopf("%s produced non-finite fitted relative risks.",
            chain_files[[chain]])
    }
    fitted_risk_sum <- fitted_risk_sum + colSums(relative_risk)
    total_draws <- total_draws + nrow(relative_risk)
    metadata[[chain]] <- data.frame(
      specification = specification,
      specification_label = specifications[[specification]]$label,
      chain = chain,
      seed = fit$seed,
      retained_draws = nrow(fit$beta),
      warmup_iterations = fit$settings$warmup_iterations,
      sampling_thin = fit$settings$sampling_thin,
      n_temperatures = fit$settings$n_temperatures,
      beta_hot = fit$settings$learned_beta_hot,
      runtime_hours = fit$runtime$elapsed_seconds / 3600,
      chain_file = normalizePath(chain_files[[chain]], winslash = "/"),
      chain_md5 = md5_file(chain_files[[chain]]),
      stringsAsFactors = FALSE
    )
    rm(fit, relative_risk)
    invisible(gc(FALSE))
  }
  fitted_risk <- reorder_vector_to_geographic(fitted_risk_sum / total_draws)
  fitted_risk <- matrix(
    fitted_risk, nrow = n_counties, ncol = n_diseases,
    dimnames = list(title_case_county(reference$county_names), disease_names)
  )

  list(
    probabilities = probabilities,
    selections = selections,
    counts = colSums(selections),
    fitted_risk = fitted_risk,
    total_draws = total_draws,
    metadata = do.call(rbind, metadata),
    boundary_file = normalizePath(boundary_file, winslash = "/"),
    boundary_md5 = md5_file(boundary_file)
  )
}

cat("Multiple-chain SEER sensitivity analysis\n")
cat("Started:", format(Sys.time()), "\n")
cat("Input root:", normalizePath(input_root, winslash = "/"), "\n")
cat("Run mode:", run_mode, "\n")
cat("Output:", normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    "\n\n")

results <- lapply(names(specifications), load_specification)
names(results) <- names(specifications)
chain_metadata <- do.call(rbind, lapply(results, `[[`, "metadata"))
if (anyDuplicated(chain_metadata[, c("specification", "seed")])) {
  stop("Duplicate chain seeds were detected within a specification.",
       call. = FALSE)
}

# RESULT: complete edge-level probabilities and FDR selections ----------------
edge_comparison <- data.frame(
  edge_id = rep(seq_len(n_edges), times = n_diseases),
  county_1 = rep(
    title_case_county(reference$county_names[reference$edges[, 1]]),
    times = n_diseases
  ),
  county_2 = rep(
    title_case_county(reference$county_names[reference$edges[, 2]]),
    times = n_diseases
  ),
  cancer = rep(disease_names, each = n_edges),
  stringsAsFactors = FALSE
)
for (specification in names(specifications)) {
  label <- specifications[[specification]]$label
  prefix <- gsub("--", "_", tolower(label), fixed = TRUE)
  edge_comparison[[paste0(prefix, "_boundary_probability")]] <-
    as.vector(results[[specification]]$probabilities)
  edge_comparison[[paste0(prefix, "_fdr_selected")]] <-
    as.vector(results[[specification]]$selections)
}

classify_boundary <- function(adjacency, mean_adjacency, mean) {
  ifelse(
    adjacency & mean_adjacency & mean, "Selected under all three",
    ifelse(
      adjacency & mean_adjacency, "Adjacency and Mean--Adjacency only",
      ifelse(
        adjacency & mean, "Adjacency and Mean only",
        ifelse(
          mean_adjacency & mean, "Mean--Adjacency and Mean only",
          ifelse(
            adjacency, "Adjacency only",
            ifelse(mean_adjacency, "Mean--Adjacency only",
                   ifelse(mean, "Mean only", "Not selected"))
          )
        )
      )
    )
  )
}
edge_comparison$selection_classification <- classify_boundary(
  edge_comparison$adjacency_fdr_selected,
  edge_comparison$mean_adjacency_fdr_selected,
  edge_comparison$mean_fdr_selected
)
edge_file <- file.path(
  tables_dir, "edge_boundary_probabilities_multichain_three_specifications.csv"
)
write_csv(edge_comparison, edge_file)

classification_levels <- c(
  "Selected under all three", "Adjacency and Mean--Adjacency only",
  "Adjacency and Mean only", "Mean--Adjacency and Mean only",
  "Adjacency only", "Mean--Adjacency only", "Mean only", "Not selected"
)
classification_summary <- do.call(rbind, lapply(disease_names, function(cancer) {
  counts <- table(factor(
    edge_comparison$selection_classification[
      edge_comparison$cancer == cancer
    ],
    levels = classification_levels
  ))
  data.frame(
    cancer = cancer,
    classification = names(counts),
    edge_count = as.integer(counts),
    proportion_of_edges = as.integer(counts) / n_edges,
    stringsAsFactors = FALSE
  )
}))
classification_file <- file.path(
  tables_dir,
  "boundary_selection_classification_multichain_three_specifications.csv"
)
write_csv(classification_summary, classification_file)

pair_definitions <- list(
  list(A = "adj", B = "meanadj", label_A = "Adjacency",
       label_B = "Mean--Adjacency", comparison = "Mean adjustment"),
  list(A = "meanadj", B = "mean", label_A = "Mean--Adjacency",
       label_B = "Mean", comparison = "Adjacency learning"),
  list(A = "adj", B = "mean", label_A = "Adjacency",
       label_B = "Mean", comparison = "Combined specification difference")
)

overlap_summary <- do.call(rbind, lapply(pair_definitions, function(pair) {
  do.call(rbind, lapply(seq_len(n_diseases), function(disease) {
    selected_A <- results[[pair$A]]$selections[, disease]
    selected_B <- results[[pair$B]]$selections[, disease]
    probability_A <- results[[pair$A]]$probabilities[, disease]
    probability_B <- results[[pair$B]]$probabilities[, disease]
    intersection_size <- sum(selected_A & selected_B)
    union_size <- sum(selected_A | selected_B)
    data.frame(
      cancer = disease_names[[disease]],
      specification_A = pair$label_A,
      specification_B = pair$label_B,
      comparison = pair$comparison,
      count_A = sum(selected_A),
      count_B = sum(selected_B),
      intersection_size = intersection_size,
      union_size = union_size,
      jaccard_similarity = if (union_size) intersection_size / union_size else NA_real_,
      A_only = sum(selected_A & !selected_B),
      B_only = sum(!selected_A & selected_B),
      pearson_probability_correlation = stats::cor(probability_A, probability_B),
      spearman_probability_correlation = stats::cor(
        probability_A, probability_B, method = "spearman"
      ),
      mean_absolute_probability_difference = mean(abs(
        probability_A - probability_B
      )),
      stringsAsFactors = FALSE
    )
  }))
}))
overlap_file <- file.path(
  tables_dir, "boundary_overlap_multichain_three_specifications.csv"
)
write_csv(overlap_summary, overlap_file)

edge_changes <- edge_comparison
edge_changes$mean_adjacency_minus_mean_probability_change <-
  edge_changes$mean_adjacency_boundary_probability -
  edge_changes$mean_boundary_probability
edge_changes$absolute_probability_change <- abs(
  edge_changes$mean_adjacency_minus_mean_probability_change
)
edge_changes <- edge_changes[order(
  edge_changes$cancer, -edge_changes$absolute_probability_change
),]
edge_changes_file <- file.path(
  tables_dir, "edge_changes_meanadj_vs_mean_multichain.csv"
)
write_csv(edge_changes, edge_changes_file)

# RESULT: main sensitivity table ----------------------------------------------
clean_overlap <- overlap_summary[
  overlap_summary$specification_A == "Mean--Adjacency" &
    overlap_summary$specification_B == "Mean",
]
sensitivity_table <- data.frame(
  Cancer = disease_names,
  Adjacency = as.integer(results$adj$counts[disease_names]),
  `Mean--Adjacency` = as.integer(results$meanadj$counts[disease_names]),
  Mean = as.integer(results$mean$counts[disease_names]),
  Jaccard = clean_overlap$jaccard_similarity[
    match(disease_names, clean_overlap$cancer)
  ],
  `Mean--Adjacency only` = clean_overlap$A_only[
    match(disease_names, clean_overlap$cancer)
  ],
  `Mean only` = clean_overlap$B_only[
    match(disease_names, clean_overlap$cancer)
  ],
  check.names = FALSE
)
sensitivity_csv <- file.path(
  tables_dir,
  "boundary_sensitivity_table_multichain_three_specifications.csv"
)
write_csv(sensitivity_table, sensitivity_csv)

latex_rows <- sprintf(
  "%s & %d & %d & %d & %.2f & %d & %d \\\\",
  sensitivity_table$Cancer,
  sensitivity_table$Adjacency,
  sensitivity_table[["Mean--Adjacency"]],
  sensitivity_table$Mean,
  sensitivity_table$Jaccard,
  sensitivity_table[["Mean--Adjacency only"]],
  sensitivity_table[["Mean only"]]
)
sensitivity_tex <- file.path(
  tables_dir,
  "boundary_sensitivity_table_multichain_three_specifications.tex"
)
write_text(c(
  "\\begin{table}[!ht]",
  "\\centering",
  paste0(
    "\\caption{Disease-specific sensitivity of pooled posterior-FDR ",
    "boundary selections based on 150,000 retained cold-posterior draws ",
    "from six adaptive-tempering runs per specification.}"
  ),
  "\\label{tab:seer_boundary_sensitivity}",
  "\\begin{tabular}{lrrrrrr}",
  "\\hline",
  paste0(
    "Cancer & Adjacency & Mean--Adjacency & Mean & Jaccard & ",
    "Mean--Adjacency only & Mean only \\\\"
  ),
  "\\hline", latex_rows, "\\hline", "\\end{tabular}", "\\end{table}"
), sensitivity_tex)

# RESULT: fitted relative-risk sensitivity ------------------------------------
fitted_risk_table <- do.call(rbind, lapply(names(specifications), function(s) {
  data.frame(
    specification = specifications[[s]]$label,
    county = rep(rownames(results[[s]]$fitted_risk), times = n_diseases),
    cancer = rep(disease_names, each = n_counties),
    posterior_mean_relative_risk = as.vector(results[[s]]$fitted_risk),
    stringsAsFactors = FALSE
  )
}))
fitted_risk_file <- file.path(
  tables_dir,
  "fitted_relative_risk_means_multichain_three_specifications.csv"
)
write_csv(fitted_risk_table, fitted_risk_file)

fitted_risk_summary <- do.call(rbind, lapply(pair_definitions, function(pair) {
  do.call(rbind, lapply(seq_len(n_diseases), function(disease) {
    A <- results[[pair$A]]$fitted_risk[, disease]
    B <- results[[pair$B]]$fitted_risk[, disease]
    data.frame(
      cancer = disease_names[[disease]],
      specification_A = pair$label_A,
      specification_B = pair$label_B,
      comparison = pair$comparison,
      pearson_correlation = stats::cor(A, B),
      spearman_correlation = stats::cor(A, B, method = "spearman"),
      mean_absolute_difference = mean(abs(A - B)),
      maximum_absolute_difference = max(abs(A - B)),
      stringsAsFactors = FALSE
    )
  }))
}))
fitted_risk_summary_file <- file.path(
  tables_dir, "fitted_risk_sensitivity_multichain_three_specifications.csv"
)
write_csv(fitted_risk_summary, fitted_risk_summary_file)

# RESULT: boundary-probability agreement figure -------------------------------
agreement_pairs <- pair_definitions[1:2]
agreement_data <- do.call(rbind, lapply(agreement_pairs, function(pair) {
  do.call(rbind, lapply(seq_len(n_diseases), function(disease) {
    data.frame(
      probability_A = results[[pair$A]]$probabilities[, disease],
      probability_B = results[[pair$B]]$probabilities[, disease],
      cancer = disease_names[[disease]],
      comparison = paste(pair$comparison, pair$label_A, "vs", pair$label_B,
                         sep = "\n"),
      stringsAsFactors = FALSE
    )
  }))
}))
agreement_annotations <- do.call(rbind, lapply(split(
  agreement_data, interaction(agreement_data$cancer, agreement_data$comparison)
), function(x) {
  data.frame(
    probability_A = 0.03, probability_B = 0.97,
    cancer = x$cancer[[1]], comparison = x$comparison[[1]],
    label = sprintf("Spearman = %.2f", stats::cor(
      x$probability_A, x$probability_B, method = "spearman"
    ))
  )
}))

# FIGURE: California county names used in the supplementary discussion --------
# st_point_on_surface keeps labels inside multipart county geometries more
# reliably than raw centroids. The map rows use the same geographic ordering as
# the adjacency edge list and the fitted-risk summaries.
county_label_map <- suppressWarnings(sf::st_point_on_surface(
  reference$county_map
))
county_label_map$county_label <- title_case_county(reference$county_names)
named_map_plot <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = reference$county_map, fill = "grey97", color = "grey45",
    linewidth = 0.25
  ) +
  ggplot2::geom_sf_text(
    data = county_label_map, ggplot2::aes(label = county_label),
    size = 2.25, color = "black", lineheight = 0.9
  ) +
  ggplot2::coord_sf(datum = NA) +
  ggplot2::theme_void()
named_map_pdf <- file.path(
  figures_dir, "named_map_multichain_california_counties.pdf"
)
named_map_png <- sub("\\.pdf$", ".png", named_map_pdf)
ggplot2::ggsave(named_map_pdf, named_map_plot, width = 8.5, height = 10.5)
ggplot2::ggsave(
  named_map_png, named_map_plot, width = 8.5, height = 10.5,
  dpi = 300, bg = "white"
)

agreement_plot <- ggplot2::ggplot(
  agreement_data,
  ggplot2::aes(x = probability_A, y = probability_B)
) +
  ggplot2::geom_abline(slope = 1, color = "grey50", linetype = "dashed") +
  ggplot2::geom_point(alpha = 0.65, size = 1.5, color = "#2C3E50") +
  ggplot2::geom_text(
    data = agreement_annotations, ggplot2::aes(label = label),
    hjust = 0, vjust = 1, size = 3
  ) +
  ggplot2::facet_grid(cancer ~ comparison) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    x = "Boundary probability under the first named model",
    y = "Boundary probability under the second named model"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 9),
    panel.grid.minor = ggplot2::element_blank()
  )
agreement_pdf <- file.path(
  figures_dir,
  "boundary_probability_agreement_multichain_three_specifications.pdf"
)
agreement_png <- sub("\\.pdf$", ".png", agreement_pdf)
ggplot2::ggsave(agreement_pdf, agreement_plot, width = 12, height = 11)
ggplot2::ggsave(agreement_png, agreement_plot, width = 12, height = 11,
                dpi = 300)

# RESULT: Mean--Adjacency versus Mean FDR boundary map -------------------------
# Intersect county boundaries rather than the polygon interiors. Some polygons
# in the maps database are not valid GEOS polygon rings, although their
# boundaries are valid and were used successfully by the main postprocessor.
# Preserve point-only Queen adjacencies as points for plotting with crosses.
county_geometry <- sf::st_geometry(reference$county_map)
edge_geometries <- vector("list", n_edges)
valid_geometry <- logical(n_edges)
for (edge in seq_len(n_edges)) {
  counties <- reference$edges[edge,]
  shared_boundary <- suppressWarnings(sf::st_intersection(
    sf::st_boundary(county_geometry[counties[[1]]]),
    sf::st_boundary(county_geometry[counties[[2]]])
  ))
  if (!length(shared_boundary) || all(sf::st_is_empty(shared_boundary))) {
    next
  }
  dimensions <- sf::st_dimension(shared_boundary)
  linear_part <- shared_boundary[dimensions == 1L]
  point_part <- shared_boundary[dimensions == 0L]
  if (length(linear_part) && !all(sf::st_is_empty(linear_part))) {
    edge_geometries[[edge]] <- sf::st_union(linear_part)[[1L]]
    valid_geometry[[edge]] <- TRUE
  } else if (length(point_part) && !all(sf::st_is_empty(point_part))) {
    edge_geometries[[edge]] <- sf::st_union(point_part)[[1L]]
    valid_geometry[[edge]] <- TRUE
  }
}
if (!all(valid_geometry)) {
  warning(
    "No drawable geometry was recovered for ", sum(!valid_geometry),
    " neighboring county pair(s); those edges are omitted from the map."
  )
}
edge_geometry_sf <- sf::st_sf(
  edge_id = which(valid_geometry),
  geometry = sf::st_sfc(
    edge_geometries[valid_geometry], crs = sf::st_crs(reference$county_map)
  )
)

map_rows <- do.call(rbind, lapply(seq_len(n_diseases), function(disease) {
  adaptive <- results$meanadj$selections[, disease]
  fixed <- results$mean$selections[, disease]
  category <- ifelse(
    adaptive & fixed, "Common",
    ifelse(adaptive, "Mean--Adjacency only", ifelse(fixed, "Mean only", NA))
  )
  data.frame(
    edge_id = which(!is.na(category)),
    cancer = disease_names[[disease]],
    category = category[!is.na(category)],
    stringsAsFactors = FALSE
  )
}))
map_edges <- merge(map_rows, edge_geometry_sf, by = "edge_id", sort = FALSE)
map_edges <- sf::st_as_sf(map_edges)
map_edges$category <- factor(
  map_edges$category,
  levels = c("Common", "Mean--Adjacency only", "Mean only")
)
map_edges$cancer <- factor(map_edges$cancer, levels = disease_names)
geometry_types <- as.character(sf::st_geometry_type(map_edges))
line_edges <- map_edges[grepl("LINE", geometry_types),]
point_edges <- map_edges[grepl("POINT", geometry_types),]

boundary_map_plot <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = reference$county_map, fill = "grey97", color = "grey55",
    linewidth = 0.25
  ) +
  ggplot2::geom_sf(
    data = line_edges, ggplot2::aes(color = category), linewidth = 1.1,
    lineend = "round"
  ) +
  ggplot2::geom_sf(
    data = point_edges, ggplot2::aes(color = category), shape = 4,
    size = 2.4, stroke = 1.1
  ) +
  ggplot2::facet_wrap(~ cancer, ncol = 2) +
  ggplot2::scale_color_manual(values = c(
    "Common" = "#3B3B3B", "Mean--Adjacency only" = "#009E73",
    "Mean only" = "#0072B2"
  ), drop = FALSE) +
  ggplot2::coord_sf(datum = NA) +
  ggplot2::labs(
    title = "Mean--Adjacency versus Mean",
    color = "Pooled FDR boundary class",
    caption = paste(
      "Both models use the same covariate-adjusted mean;",
      "only adjacency learning differs."
    )
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 12),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )
boundary_map_pdf <- file.path(
  figures_dir, "boundary_maps_multichain_meanadj_vs_mean.pdf"
)
boundary_map_png <- sub("\\.pdf$", ".png", boundary_map_pdf)
ggplot2::ggsave(boundary_map_pdf, boundary_map_plot, width = 10, height = 9)
ggplot2::ggsave(boundary_map_png, boundary_map_plot, width = 10, height = 9,
                dpi = 300, bg = "white")

# Reproducibility records ------------------------------------------------------
metadata_file <- file.path(
  diagnostics_dir, "chain_metadata_multichain_three_specifications.csv"
)
write_csv(chain_metadata, metadata_file)
input_manifest <- rbind(
  chain_metadata[, c(
    "specification", "chain", "chain_file", "chain_md5",
    "retained_draws", "seed"
  )],
  do.call(rbind, lapply(names(results), function(specification) {
    data.frame(
      specification = specification,
      chain = NA_integer_,
      chain_file = results[[specification]]$boundary_file,
      chain_md5 = results[[specification]]$boundary_md5,
      retained_draws = results[[specification]]$total_draws,
      seed = NA_integer_,
      stringsAsFactors = FALSE
    )
  }))
)
names(input_manifest)[names(input_manifest) == "chain_file"] <- "input_file"
names(input_manifest)[names(input_manifest) == "chain_md5"] <- "md5"
input_manifest_file <- file.path(
  diagnostics_dir, "input_manifest_sensitivity_multichain.csv"
)
write_csv(input_manifest, input_manifest_file)

artifact_paths <- c(
  edge_file, classification_file, overlap_file, edge_changes_file,
  sensitivity_csv, sensitivity_tex,
  fitted_risk_file, fitted_risk_summary_file,
  named_map_pdf, named_map_png,
  agreement_pdf, agreement_png, boundary_map_pdf, boundary_map_png,
  metadata_file, input_manifest_file
)
artifact_manifest <- data.frame(
  artifact = basename(artifact_paths),
  type = ifelse(
    dirname(artifact_paths) == figures_dir, "figure",
    ifelse(dirname(artifact_paths) == tables_dir, "table", "diagnostic")
  ),
  path = normalizePath(artifact_paths, winslash = "/", mustWork = TRUE),
  md5 = vapply(artifact_paths, md5_file, character(1L)),
  stringsAsFactors = FALSE
)
artifact_manifest_file <- file.path(
  output_dir, "artifact_manifest_sensitivity_multichain.csv"
)
write_csv(artifact_manifest, artifact_manifest_file)

saveRDS(
  list(
    generated_at = Sys.time(),
    run_mode = run_mode,
    fdr_alpha = fdr_alpha,
    boundary_counts = lapply(results, `[[`, "counts"),
    boundary_overlap = overlap_summary,
    sensitivity_table = sensitivity_table,
    fitted_risk_sensitivity = fitted_risk_summary,
    input_manifest = input_manifest,
    artifact_manifest = artifact_manifest
  ),
  file.path(output_dir, "sensitivity_results_multichain.rds")
)

cat("\nSensitivity analysis complete.\n")
cat("Artifacts:", normalizePath(output_dir, winslash = "/"), "\n")
