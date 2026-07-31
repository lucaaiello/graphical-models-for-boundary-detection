# Posterior predictive checks and three-specification SEER sensitivity analysis
#
# This script uses the three saved posterior workspaces that generated the
# manuscript and supplementary results. It never calls MADAGAR(), sourceCpp(),
# or any MCMC sampler.
#
# Object orientation and model construction were traced from main.R:
#   * mcmc_samples$beta: retained draws x regression coefficients
#   * mcmc_samples$phi: retained draws x (58 counties x 4 cancers)
#   * X, Y, and E: model order, with cancer blocks ordered as
#       Lung, Esophageal, Larynx, Colorectal
#   * final_perm: maps model-order counties to the original map order
#   * X1[, c(1, 2, 4, 6)] etc.: intercept, smoking, age 65+, poverty
#   * neighbor_list0: the canonical 139-edge ordering used in the figures
#   * path: plotting coordinates for those same 139 edges
#
# The FDR implementation below reproduces main.R exactly: it considers the
# first 135 ordered thresholds, estimates FDR using all probabilities at or
# above each threshold (including ties), and uses zeta = 0.05.

rm(list = ls())

# Configuration ----------------------------------------------------------------

set.seed(20260724L)

project_root_markers <- c("README.md", "data analysis")
if (!all(file.exists(project_root_markers))) {
  stop(
    "Run this script from the repository root. Expected to find README.md ",
    "and the 'data analysis' directory."
  )
}

required_packages <- c("ggplot2", "sf", "maps", "matrixStats")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Run install.R before this analysis."
  )
}

output_dir <- file.path(
  "data analysis",
  "posterior_predictive_sensitivity_output"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cancer_names <- c("Lung", "Esophageal", "Larynx", "Colorectal")
n_cancers <- length(cancer_names)
n_counties <- 58L
fdr_zeta <- 0.05
maximum_fdr_rank <- 135L
unusual_support_threshold <- 0.02
ppc_seed <- 20260724L

specification_definitions <- list(
  primary_adaptive = list(
    key = "primary_adaptive",
    file = file.path("data analysis", "results_adj.RData"),
    label = "Adjacency",
    short_label = "Adjacency",
    mean_covariate_columns = 1L,
    comparison_role = paste(
      "Adjacency analysis; comparison with Mean--Adjacency",
      "isolates mean adjustment while retaining adaptive adjacency."
    )
  ),
  adjusted_adaptive = list(
    key = "adjusted_adaptive",
    file = file.path("data analysis", "results_meanadj.RData"),
    label = "Mean--Adjacency",
    short_label = "Mean--Adjacency",
    mean_covariate_columns = c(1L, 2L, 4L, 6L),
    comparison_role = paste(
      "Mean--Adjacency specification; comparison with Mean is the clean",
      "adjacency-learning comparison."
    )
  ),
  adjusted_fixed = list(
    key = "adjusted_fixed",
    file = file.path("data analysis", "results_mean.RData"),
    label = "Mean",
    short_label = "Mean",
    mean_covariate_columns = c(1L, 2L, 4L, 6L),
    comparison_role = paste(
      "Mean specification with fixed adjacency; it has the same mean",
      "covariates as Mean--Adjacency."
    )
  )
)

input_files <- vapply(
  specification_definitions,
  function(definition) definition$file,
  character(1)
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing saved posterior workspace(s): ",
    paste(missing_files, collapse = ", ")
  )
}

# General utilities -------------------------------------------------------------

stopf <- function(format_string, ...) {
  stop(sprintf(format_string, ...), call. = FALSE)
}

write_output_csv <- function(object, filename) {
  path <- file.path(output_dir, filename)
  utils::write.csv(object, path, row.names = FALSE, na = "")
  path
}

write_output_text <- function(lines, filename) {
  path <- file.path(output_dir, filename)
  writeLines(lines, path, useBytes = TRUE)
  path
}

format_number <- function(x, digits = 3L) {
  formatC(x, digits = digits, format = "f")
}

title_case_county <- function(x) {
  tools::toTitleCase(as.character(x))
}

row_standard_deviation <- function(x) {
  if (!is.matrix(x) || ncol(x) < 2L) {
    stop("row_standard_deviation() requires a matrix with at least two columns.")
  }
  row_sum <- rowSums(x)
  variance <- (
    rowSums(x * x) - (row_sum * row_sum) / ncol(x)
  ) / (ncol(x) - 1L)
  sqrt(pmax(variance, 0))
}

reorder_draws_to_geographic <- function(draw_matrix, final_perm) {
  expected_columns <- n_counties * n_cancers
  if (!is.matrix(draw_matrix) || ncol(draw_matrix) != expected_columns) {
    stopf(
      "Expected a draw matrix with %d columns; observed %s.",
      expected_columns,
      paste(dim(draw_matrix), collapse = " x ")
    )
  }
  inverse_permutation <- order(final_perm)
  do.call(
    cbind,
    lapply(seq_len(n_cancers), function(disease) {
      columns <- (
        (disease - 1L) * n_counties + 1L
      ):(disease * n_counties)
      draw_matrix[, columns, drop = FALSE][
        ,
        inverse_permutation,
        drop = FALSE
      ]
    })
  )
}

reorder_vector_to_geographic <- function(vector, final_perm) {
  if (length(vector) != n_counties * n_cancers) {
    stopf(
      "Expected a vector of length %d; observed %d.",
      n_counties * n_cancers,
      length(vector)
    )
  }
  inverse_permutation <- order(final_perm)
  unlist(
    lapply(seq_len(n_cancers), function(disease) {
      indices <- (
        (disease - 1L) * n_counties + 1L
      ):(disease * n_counties)
      vector[indices][inverse_permutation]
    }),
    use.names = FALSE
  )
}

extract_disease_block <- function(draw_matrix, disease) {
  columns <- (
    (disease - 1L) * n_counties + 1L
  ):(disease * n_counties)
  draw_matrix[, columns, drop = FALSE]
}

extract_disease_vector <- function(vector, disease) {
  indices <- (
    (disease - 1L) * n_counties + 1L
  ):(disease * n_counties)
  vector[indices]
}

build_expected_design_matrix <- function(X_by_disease, active_columns) {
  columns_per_disease <- length(active_columns)
  result <- matrix(
    0,
    nrow = n_counties * n_cancers,
    ncol = columns_per_disease * n_cancers
  )
  for (disease in seq_len(n_cancers)) {
    rows <- (
      (disease - 1L) * n_counties + 1L
    ):(disease * n_counties)
    columns <- (
      (disease - 1L) * columns_per_disease + 1L
    ):(disease * columns_per_disease)
    result[rows, columns] <- X_by_disease[[disease]][
      ,
      active_columns,
      drop = FALSE
    ]
  }
  result
}

normalise_stored_selection <- function(selection, number_of_edges, object_name) {
  if (length(selection) == 0L) {
    return(integer(number_of_edges))
  }
  if (length(selection) != number_of_edges) {
    stopf(
      "Stored object %s has length %d; expected either 0 or %d.",
      object_name,
      length(selection),
      number_of_edges
    )
  }
  if (!all(selection %in% c(0, 1))) {
    stopf("Stored object %s is not a binary selection vector.", object_name)
  }
  as.integer(selection)
}

# Strict saved-workspace extraction --------------------------------------------

load_saved_specification <- function(definition) {
  workspace <- new.env(parent = baseenv())
  loaded_names <- load(definition$file, envir = workspace)

  required_objects <- c(
    "mcmc_samples", "X", "X1", "X2", "X3", "X4",
    "Y", "Y1", "Y2", "Y3", "Y4",
    "E", "E1", "E2", "E3", "E4",
    "final_perm", "county.ID", "Adj", "neighbor_list0",
    "neighbor_name",
    "rate_lung", "rate_esophagus", "rate_larynx", "rate_colrect"
  )
  if (definition$key == "primary_adaptive") {
    required_objects <- c(required_objects, "path")
  }
  missing_objects <- setdiff(required_objects, loaded_names)
  if (length(missing_objects) > 0L) {
    stopf(
      "%s is missing required object(s): %s",
      definition$file,
      paste(missing_objects, collapse = ", ")
    )
  }
  stored_boundary_objects <- c(
    paste0("pvij", seq_len(n_cancers)),
    paste0("est_diff", seq_len(n_cancers))
  )
  stored_boundary_objects_present <- stored_boundary_objects %in% loaded_names
  if (any(stored_boundary_objects_present) &&
      !all(stored_boundary_objects_present)) {
    stopf(
      "%s contains only some stored boundary summaries: %s",
      definition$file,
      paste(
        stored_boundary_objects[stored_boundary_objects_present],
        collapse = ", "
      )
    )
  }
  has_stored_boundary_summaries <- all(stored_boundary_objects_present)

  mcmc_samples <- get("mcmc_samples", envir = workspace)
  required_components <- c("beta", "phi")
  if (
    !is.list(mcmc_samples) ||
      !all(required_components %in% names(mcmc_samples))
  ) {
    stopf(
      "%s: mcmc_samples must be a named list containing beta and phi.",
      definition$file
    )
  }

  beta <- mcmc_samples[["beta"]]
  phi <- mcmc_samples[["phi"]]
  X <- get("X", envir = workspace)
  Y <- get("Y", envir = workspace)
  E <- get("E", envir = workspace)
  final_perm <- as.integer(get("final_perm", envir = workspace))
  county_names <- as.character(get("county.ID", envir = workspace))
  adjacency <- as.matrix(get("Adj", envir = workspace))
  edges <- as.matrix(get("neighbor_list0", envir = workspace))
  edge_names <- as.character(get("neighbor_name", envir = workspace))

  if (!is.matrix(beta) || !is.matrix(phi) || !is.matrix(X)) {
    stopf("%s: beta, phi, and X must all be matrices.", definition$file)
  }
  if (nrow(beta) != nrow(phi)) {
    stopf(
      "%s: beta has %d draws but phi has %d; orientation is incompatible.",
      definition$file,
      nrow(beta),
      nrow(phi)
    )
  }
  if (ncol(beta) != ncol(X)) {
    stopf(
      "%s: beta has %d columns but X has %d columns.",
      definition$file,
      ncol(beta),
      ncol(X)
    )
  }
  if (
    nrow(X) != n_counties * n_cancers ||
      ncol(phi) != n_counties * n_cancers ||
      length(Y) != n_counties * n_cancers ||
      length(E) != n_counties * n_cancers
  ) {
    stopf(
      "%s: expected 232 observations in X, phi, Y, and E.",
      definition$file
    )
  }
  if (
    any(!is.finite(beta)) ||
      any(!is.finite(phi)) ||
      any(!is.finite(X)) ||
      any(!is.finite(Y)) ||
      any(!is.finite(E))
  ) {
    stopf("%s contains non-finite posterior or model values.", definition$file)
  }
  if (any(E <= 0) || any(Y < 0) || any(Y != floor(Y))) {
    stopf("%s contains invalid observed or expected counts.", definition$file)
  }
  if (
    length(final_perm) != n_counties ||
      !identical(sort(final_perm), seq_len(n_counties))
  ) {
    stopf("%s: final_perm is not a permutation of 1:58.", definition$file)
  }
  if (
    length(county_names) != n_counties ||
      anyNA(county_names) ||
      anyDuplicated(county_names)
  ) {
    stopf("%s: county.ID must contain 58 unique county names.", definition$file)
  }
  if (
    !identical(dim(adjacency), c(n_counties, n_counties)) ||
      !isTRUE(all.equal(
        unname(adjacency),
        unname(t(adjacency)),
        tolerance = 0
      )) ||
      any(diag(adjacency) != 0) ||
      !all(adjacency %in% c(0, 1))
  ) {
    stopf("%s: Adj is not a valid symmetric binary adjacency matrix.", definition$file)
  }
  if (
    !is.null(colnames(adjacency)) &&
      !identical(as.character(colnames(adjacency)), county_names)
  ) {
    stopf("%s: Adj column names and county.ID are inconsistent.", definition$file)
  }
  if (
    ncol(edges) != 2L ||
      any(edges[, 1] >= edges[, 2]) ||
      any(edges < 1L) ||
      any(edges > n_counties) ||
      any(adjacency[edges] != 1)
  ) {
    stopf("%s: neighbor_list0 is not a valid canonical edge list.", definition$file)
  }
  edge_keys <- paste(edges[, 1], edges[, 2], sep = "-")
  if (anyDuplicated(edge_keys) || length(edge_names) != nrow(edges)) {
    stopf("%s: edge ordering contains duplicates or missing edge names.", definition$file)
  }

  Y_by_disease <- lapply(
    seq_len(n_cancers),
    function(disease) get(paste0("Y", disease), envir = workspace)
  )
  E_by_disease <- lapply(
    seq_len(n_cancers),
    function(disease) get(paste0("E", disease), envir = workspace)
  )
  if (
    !isTRUE(all.equal(Y, unlist(Y_by_disease), tolerance = 0)) ||
      !isTRUE(all.equal(E, unlist(E_by_disease), tolerance = 0))
  ) {
    stopf("%s: the cancer-block ordering of Y or E is inconsistent.", definition$file)
  }

  rate_objects <- c(
    "rate_lung", "rate_esophagus", "rate_larynx", "rate_colrect"
  )
  expected_sites <- c(
    "Lung and Bronchus", "Esophagus", "Larynx", "Colon and Rectum"
  )
  observed_sites <- vapply(
    rate_objects,
    function(object_name) {
      rate_data <- get(object_name, envir = workspace)
      if (!"Site.code" %in% names(rate_data)) {
        stopf("%s lacks Site.code.", object_name)
      }
      sites <- unique(as.character(rate_data[["Site.code"]]))
      if (length(sites) != 1L) {
        stopf("%s does not identify exactly one cancer.", object_name)
      }
      sites
    },
    character(1)
  )
  if (!identical(unname(observed_sites), expected_sites)) {
    stopf(
      "%s: cancer ordering is inconsistent. Observed: %s",
      definition$file,
      paste(observed_sites, collapse = ", ")
    )
  }

  X_by_disease <- lapply(
    seq_len(n_cancers),
    function(disease) get(paste0("X", disease), envir = workspace)
  )
  expected_X <- build_expected_design_matrix(
    X_by_disease,
    definition$mean_covariate_columns
  )
  if (!isTRUE(all.equal(X, expected_X, check.attributes = FALSE))) {
    stopf(
      "%s: X does not match the design construction traced from main.R.",
      definition$file
    )
  }

  stored_probabilities <- NULL
  stored_selections <- NULL
  if (has_stored_boundary_summaries) {
    stored_probabilities <- do.call(
      cbind,
      lapply(
        seq_len(n_cancers),
        function(disease) {
          probability <- as.numeric(
            get(paste0("pvij", disease), envir = workspace)
          )
          if (length(probability) != nrow(edges)) {
            stopf(
              "%s: pvij%d has length %d instead of %d.",
              definition$file,
              disease,
              length(probability),
              nrow(edges)
            )
          }
          probability
        }
      )
    )
    stored_selections <- do.call(
      cbind,
      lapply(
        seq_len(n_cancers),
        function(disease) {
          object_name <- paste0("est_diff", disease)
          normalise_stored_selection(
            get(object_name, envir = workspace),
            nrow(edges),
            object_name
          )
        }
      )
    )
    colnames(stored_probabilities) <- cancer_names
    colnames(stored_selections) <- cancer_names
  }

  path_coordinates <- NULL
  if (definition$key == "primary_adaptive") {
    path_coordinates <- get("path", envir = workspace)
    if (
      !is.list(path_coordinates) ||
        length(path_coordinates) != nrow(edges)
    ) {
      stopf("%s: path does not match neighbor_list0.", definition$file)
    }
    path_rows <- vapply(path_coordinates, nrow, integer(1))
    if (any(path_rows < 1L)) {
      stopf("%s: at least one stored edge path is empty.", definition$file)
    }
  }

  result <- list(
    key = definition$key,
    file = definition$file,
    label = definition$label,
    short_label = definition$short_label,
    comparison_role = definition$comparison_role,
    beta = unname(beta),
    phi = unname(phi),
    X = unname(X),
    Y = as.numeric(Y),
    E = as.numeric(E),
    final_perm = final_perm,
    county_names = county_names,
    adjacency = unname(adjacency),
    edges = unname(edges),
    edge_names = edge_names,
    stored_probabilities = stored_probabilities,
    stored_selections = stored_selections,
    has_stored_boundary_summaries = has_stored_boundary_summaries,
    path_coordinates = path_coordinates,
    retained_draws = nrow(beta),
    beta_dimensions = dim(beta),
    phi_dimensions = dim(phi),
    X_dimensions = dim(X),
    objects_used = c(
      "mcmc_samples$beta", "mcmc_samples$phi", "X", "X1", "X2",
      "X3", "X4", "Y", "Y1", "Y2", "Y3", "Y4", "E", "E1",
      "E2", "E3", "E4", "final_perm", "county.ID", "Adj",
      "neighbor_list0", "neighbor_name",
      if (has_stored_boundary_summaries) "pvij1:pvij4" else NULL,
      if (has_stored_boundary_summaries) "est_diff1:est_diff4" else NULL,
      if (definition$key == "primary_adaptive") "path" else NULL
    )
  )

  rm(list = loaded_names, envir = workspace)
  rm(workspace, mcmc_samples, beta, phi, X, expected_X)
  invisible(gc())
  result
}

cat("Loading and validating saved posterior workspaces...\n")
saved_specifications <- lapply(
  specification_definitions,
  load_saved_specification
)

# Common-data consistency checks ------------------------------------------------

reference <- saved_specifications[["primary_adaptive"]]
common_components <- c(
  "Y", "E", "final_perm", "county_names", "adjacency", "edges", "edge_names"
)
for (key in setdiff(names(saved_specifications), "primary_adaptive")) {
  candidate <- saved_specifications[[key]]
  for (component in common_components) {
    if (!isTRUE(all.equal(
      reference[[component]],
      candidate[[component]],
      tolerance = 0,
      check.attributes = TRUE
    ))) {
      stopf(
        "%s and %s disagree on common object %s.",
        reference$file,
        candidate$file,
        component
      )
    }
  }
}

number_of_edges <- nrow(reference$edges)
island_count <- sum(rowSums(reference$adjacency) == 0)

# Boundary probabilities and exact manuscript FDR ------------------------------

compute_boundary_probabilities <- function(phi, final_perm, edges) {
  phi_geographic <- reorder_draws_to_geographic(phi, final_perm)
  probabilities <- matrix(
    NA_real_,
    nrow = nrow(edges),
    ncol = n_cancers,
    dimnames = list(NULL, cancer_names)
  )
  for (disease in seq_len(n_cancers)) {
    disease_phi <- extract_disease_block(phi_geographic, disease)
    probabilities[, disease] <- colMeans(
      disease_phi[, edges[, 1], drop = FALSE] !=
        disease_phi[, edges[, 2], drop = FALSE]
    )
  }
  probabilities
}

apply_manuscript_fdr <- function(probabilities, zeta = 0.05) {
  if (
    any(!is.finite(probabilities)) ||
      any(probabilities < 0) ||
      any(probabilities > 1)
  ) {
    stop("Boundary probabilities must be finite and in [0, 1].")
  }
  ranks <- seq_len(min(maximum_fdr_rank, length(probabilities)))
  thresholds <- sort(probabilities, decreasing = TRUE)[ranks]
  estimated_fdr <- vapply(
    thresholds,
    function(threshold) {
      selected <- probabilities >= threshold
      mean(1 - probabilities[selected])
    },
    numeric(1)
  )
  accepted_threshold_count <- sum(estimated_fdr <= zeta)
  if (accepted_threshold_count == 0L) {
    selection <- integer(length(probabilities))
    selected_threshold <- NA_real_
  } else {
    selected_threshold <- thresholds[[accepted_threshold_count]]
    selection <- as.integer(probabilities >= selected_threshold)
  }
  list(
    selection = selection,
    count = sum(selection),
    threshold = selected_threshold,
    ranks = ranks,
    estimated_fdr = estimated_fdr,
    accepted_threshold_count = accepted_threshold_count
  )
}

boundary_results <- lapply(
  saved_specifications,
  function(specification) {
    probabilities <- compute_boundary_probabilities(
      specification$phi,
      specification$final_perm,
      specification$edges
    )
    selections <- matrix(
      0L,
      nrow = number_of_edges,
      ncol = n_cancers,
      dimnames = list(NULL, cancer_names)
    )
    counts <- integer(n_cancers)
    fdr_details <- vector("list", n_cancers)
    for (disease in seq_len(n_cancers)) {
      fdr <- apply_manuscript_fdr(
        probabilities[, disease],
        zeta = fdr_zeta
      )
      selections[, disease] <- fdr$selection
      counts[[disease]] <- fdr$count
      fdr_details[[disease]] <- fdr
    }
    names(counts) <- cancer_names

    probability_agreement <- NA
    selection_agreement <- NA
    if (specification$has_stored_boundary_summaries) {
      probability_agreement <- max(
        abs(probabilities - specification$stored_probabilities)
      ) <= 1e-12
      selection_agreement <- identical(
        unname(selections),
        unname(specification$stored_selections)
      )
    }
    if (specification$has_stored_boundary_summaries &&
        (!probability_agreement || !selection_agreement)) {
      stopf(
        "%s: recomputed boundaries do not reproduce the stored analysis.",
        specification$file
      )
    }
    list(
      probabilities = probabilities,
      selections = selections,
      counts = counts,
      fdr_details = fdr_details,
      probability_agreement = probability_agreement,
      selection_agreement = selection_agreement
    )
  }
)

required_primary_counts <- c(
  Lung = 100L,
  Esophageal = 21L,
  Larynx = 33L,
  Colorectal = 63L
)
primary_counts_reproduced <- identical(
  as.integer(boundary_results[["primary_adaptive"]]$counts),
  as.integer(required_primary_counts)
)
if (!primary_counts_reproduced) {
  stopf(
    paste0(
      "Adjacency boundary validation failed. Expected %s; reproduced %s. ",
      "Check draw orientation, edge order, and the FDR rule."
    ),
    paste(required_primary_counts, collapse = ", "),
    paste(
      boundary_results[["primary_adaptive"]]$counts,
      collapse = ", "
    )
  )
}

# Validation report -------------------------------------------------------------

validation_lines <- c(
  "SEER saved-run validation report",
  "================================",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste("Cancer ordering:", paste(cancer_names, collapse = ", ")),
  paste("Counties:", n_counties),
  paste("Geographic edges:", number_of_edges),
  paste("Adjacency weights: symmetric binary; islands:", island_count),
  paste("Bayesian FDR zeta:", fdr_zeta),
  paste("Maximum ordered FDR rank used in main.R:", maximum_fdr_rank),
  ""
)
for (key in names(saved_specifications)) {
  specification <- saved_specifications[[key]]
  boundaries <- boundary_results[[key]]
  validation_lines <- c(
    validation_lines,
    paste0(specification$label, ":"),
    paste("  File:", specification$file),
    paste(
      "  Objects used:",
      paste(specification$objects_used, collapse = ", ")
    ),
    paste("  Retained posterior draws:", specification$retained_draws),
    paste(
      "  beta dimensions:",
      paste(specification$beta_dimensions, collapse = " x ")
    ),
    paste(
      "  phi dimensions:",
      paste(specification$phi_dimensions, collapse = " x ")
    ),
    paste(
      "  X dimensions:",
      paste(specification$X_dimensions, collapse = " x ")
    ),
    paste(
      "  Reproduced FDR boundary counts:",
      paste(
        paste(cancer_names, boundaries$counts, sep = "="),
        collapse = ", "
      )
    ),
    if (specification$has_stored_boundary_summaries) {
      paste(
        "  Agreement with stored probabilities/selections:",
        boundaries$probability_agreement &&
          boundaries$selection_agreement
      )
    } else {
      paste(
        "  Stored probabilities/selections:",
        "not present; recomputed from the saved posterior draws"
      )
    },
    ""
  )
}
validation_lines <- c(
  validation_lines,
  paste(
    "Adjacency manuscript counts 100/21/33/63 reproduced:",
    primary_counts_reproduced
  )
)
validation_report_path <- write_output_text(
  validation_lines,
  "validation_report.txt"
)
cat(paste(validation_lines, collapse = "\n"), "\n\n")

# Posterior predictive diagnostics ---------------------------------------------

moran_rows <- function(values, adjacency) {
  if (ncol(values) != nrow(adjacency)) {
    stop("Moran's I input and adjacency dimensions are incompatible.")
  }
  centered <- values - rowMeans(values)
  denominator <- rowSums(centered * centered)
  numerator <- rowSums((centered %*% adjacency) * centered)
  result <- ncol(values) / sum(adjacency) * numerator / denominator
  result[denominator == 0] <- NA_real_
  result
}

edge_contrast_rows <- function(counts, expected, edges) {
  # The 95th percentile of absolute neighboring stabilized log-SIR
  # differences is used instead of the maximum raw SIR difference. This is
  # less dominated by a single zero/small count for rare cancers while still
  # targeting unusually large local disparities.
  log_sir <- log(
    sweep(counts + 0.5, 2, expected, FUN = "/")
  )
  edge_differences <- abs(
    log_sir[, edges[, 1], drop = FALSE] -
      log_sir[, edges[, 2], drop = FALSE]
  )
  matrixStats::rowQuantiles(
    edge_differences,
    probs = 0.95,
    drop = TRUE
  )
}

summarise_predictive_statistic <- function(
  observed,
  replicated,
  specification,
  cancer,
  diagnostic,
  diagnostic_label
) {
  replicated <- replicated[is.finite(replicated)]
  if (length(replicated) == 0L || !is.finite(observed)) {
    stopf(
      "No finite posterior predictive values for %s, %s, %s.",
      specification,
      cancer,
      diagnostic
    )
  }
  lower_probability <- mean(replicated <= observed)
  upper_probability <- mean(replicated >= observed)
  interval <- stats::quantile(
    replicated,
    probs = c(0.025, 0.5, 0.975),
    names = FALSE,
    type = 8
  )
  data.frame(
    specification = specification,
    cancer = cancer,
    diagnostic = diagnostic,
    diagnostic_label = diagnostic_label,
    observed = observed,
    predictive_mean = mean(replicated),
    predictive_median = interval[[2]],
    predictive_lower_95 = interval[[1]],
    predictive_upper_95 = interval[[3]],
    lower_tail_probability = lower_probability,
    upper_tail_probability = upper_probability,
    two_sided_ppp = min(
      1,
      2 * min(lower_probability, upper_probability)
    ),
    observed_within_95 = (
      observed >= interval[[1]] && observed <= interval[[3]]
    ),
    stringsAsFactors = FALSE
  )
}

pointwise_predictive_summary <- function(
  observed,
  expected,
  replicated,
  counties,
  cancer,
  specification
) {
  predictive_mean <- colMeans(replicated)
  predictive_sd <- matrixStats::colSds(replicated)
  predictive_interval <- matrixStats::colQuantiles(
    replicated,
    probs = c(0.025, 0.5, 0.975),
    drop = FALSE
  )

  # Mid-tail probabilities account for discreteness by assigning half of the
  # posterior predictive mass at the observed count to each tail.
  lower_tail <- vapply(
    seq_along(observed),
    function(index) {
      mean(replicated[, index] < observed[[index]]) +
        0.5 * mean(replicated[, index] == observed[[index]])
    },
    numeric(1)
  )
  upper_tail <- vapply(
    seq_along(observed),
    function(index) {
      mean(replicated[, index] > observed[[index]]) +
        0.5 * mean(replicated[, index] == observed[[index]])
    },
    numeric(1)
  )
  two_sided <- pmin(1, 2 * pmin(lower_tail, upper_tail))
  predictive_residual <- ifelse(
    predictive_sd > 0,
    (observed - predictive_mean) / predictive_sd,
    NA_real_
  )

  data.frame(
    county = title_case_county(counties),
    cancer = cancer,
    specification = specification,
    observed_count = observed,
    expected_count = expected,
    predictive_mean = predictive_mean,
    predictive_median = predictive_interval[, 2],
    predictive_lower_95 = predictive_interval[, 1],
    predictive_upper_95 = predictive_interval[, 3],
    pointwise_95_coverage = (
      observed >= predictive_interval[, 1] &
        observed <= predictive_interval[, 3]
    ),
    predictive_residual = predictive_residual,
    lower_tail_probability = lower_tail,
    upper_tail_probability = upper_tail,
    two_sided_tail_probability = two_sided,
    unusual_observation = two_sided < unusual_support_threshold,
    descriptive_threshold = unusual_support_threshold,
    stringsAsFactors = FALSE
  )
}

run_posterior_predictive_analysis <- function(specification, seed) {
  set.seed(seed)
  number_of_draws <- specification$retained_draws

  linear_predictor <- tcrossprod(
    specification$beta,
    specification$X
  ) + specification$phi
  relative_risk <- exp(linear_predictor)
  if (any(!is.finite(relative_risk))) {
    stopf("%s produced non-finite fitted relative risks.", specification$file)
  }
  poisson_mean <- sweep(
    relative_risk,
    2,
    specification$E,
    FUN = "*"
  )
  if (any(!is.finite(poisson_mean)) || any(poisson_mean < 0)) {
    stopf("%s produced invalid posterior predictive means.", specification$file)
  }
  replicated_model_order <- matrix(
    stats::rpois(length(poisson_mean), lambda = as.vector(poisson_mean)),
    nrow = number_of_draws,
    ncol = n_counties * n_cancers
  )
  fitted_risk_mean_model_order <- colMeans(relative_risk)

  rm(linear_predictor, relative_risk, poisson_mean)
  invisible(gc())

  replicated <- reorder_draws_to_geographic(
    replicated_model_order,
    specification$final_perm
  )
  observed <- reorder_vector_to_geographic(
    specification$Y,
    specification$final_perm
  )
  expected <- reorder_vector_to_geographic(
    specification$E,
    specification$final_perm
  )
  fitted_risk_mean <- reorder_vector_to_geographic(
    fitted_risk_mean_model_order,
    specification$final_perm
  )
  rm(replicated_model_order, fitted_risk_mean_model_order)
  invisible(gc())

  ppc_rows <- list()
  unusual_rows <- list()
  fitted_risk_matrix <- matrix(
    fitted_risk_mean,
    nrow = n_counties,
    ncol = n_cancers,
    dimnames = list(specification$county_names, cancer_names)
  )

  row_counter <- 1L
  for (disease in seq_len(n_cancers)) {
    replicated_disease <- extract_disease_block(replicated, disease)
    observed_disease <- extract_disease_vector(observed, disease)
    expected_disease <- extract_disease_vector(expected, disease)
    sir_replicated <- sweep(
      replicated_disease,
      2,
      expected_disease,
      FUN = "/"
    )
    sir_observed <- observed_disease / expected_disease

    replicated_statistics <- list(
      total = rowSums(replicated_disease),
      sir_sd = row_standard_deviation(sir_replicated),
      moran_i = moran_rows(
        sir_replicated,
        specification$adjacency
      ),
      edge_log_sir_q95 = edge_contrast_rows(
        replicated_disease,
        expected_disease,
        specification$edges
      )
    )
    observed_statistics <- c(
      total = sum(observed_disease),
      sir_sd = stats::sd(sir_observed),
      moran_i = moran_rows(
        matrix(sir_observed, nrow = 1L),
        specification$adjacency
      )[[1]],
      edge_log_sir_q95 = edge_contrast_rows(
        matrix(observed_disease, nrow = 1L),
        expected_disease,
        specification$edges
      )[[1]]
    )
    diagnostic_labels <- c(
      total = "Total disease burden",
      sir_sd = "Across-county SIR SD",
      moran_i = "Moran's I of county SIR",
      edge_log_sir_q95 = "95th percentile neighboring log-SIR contrast"
    )

    for (diagnostic in names(replicated_statistics)) {
      ppc_rows[[row_counter]] <- summarise_predictive_statistic(
        observed = observed_statistics[[diagnostic]],
        replicated = replicated_statistics[[diagnostic]],
        specification = specification$label,
        cancer = cancer_names[[disease]],
        diagnostic = diagnostic,
        diagnostic_label = diagnostic_labels[[diagnostic]]
      )
      row_counter <- row_counter + 1L
    }

    unusual_rows[[disease]] <- pointwise_predictive_summary(
      observed = observed_disease,
      expected = expected_disease,
      replicated = replicated_disease,
      counties = specification$county_names,
      cancer = cancer_names[[disease]],
      specification = specification$label
    )
  }

  rm(replicated)
  invisible(gc())
  list(
    ppc_summary = do.call(rbind, ppc_rows),
    unusual = do.call(rbind, unusual_rows),
    fitted_risk_mean = fitted_risk_matrix,
    draws_used = number_of_draws
  )
}

cat("Generating posterior predictive counts from saved draws...\n")
ppc_results <- Map(
  function(specification, index) {
    run_posterior_predictive_analysis(
      specification,
      seed = ppc_seed + index
    )
  },
  saved_specifications,
  seq_along(saved_specifications)
)

ppc_summary <- do.call(
  rbind,
  lapply(ppc_results, function(result) result$ppc_summary)
)
row.names(ppc_summary) <- NULL
ppc_summary_path <- write_output_csv(ppc_summary, "ppc_summary.csv")

unusual_observations <- do.call(
  rbind,
  lapply(ppc_results, function(result) result$unusual)
)
row.names(unusual_observations) <- NULL
unusual_observations <- unusual_observations[
  order(
    unusual_observations$two_sided_tail_probability,
    -abs(unusual_observations$predictive_residual)
  ),
]
unusual_path <- write_output_csv(
  unusual_observations,
  "ppc_unusual_county_cancer.csv"
)

primary_label <- specification_definitions$primary_adaptive$label
primary_unusual <- unusual_observations[
  unusual_observations$specification == primary_label,
]
primary_unusual <- head(primary_unusual, 12L)
primary_unusual_lines <- c(
  paste0(
    "Most unusual county-cancer observations under Adjacency"
  ),
  "==========================================================================",
  paste(
    sum(primary_unusual$unusual_observation),
    "of the observations shown crossed the descriptive mid-tail threshold <",
    unusual_support_threshold,
    "; the table is ordered from lowest posterior predictive support."
  ),
  capture.output(
    print(
      primary_unusual[
        ,
        c(
          "county", "cancer", "observed_count", "predictive_mean",
          "predictive_lower_95", "predictive_upper_95",
          "predictive_residual", "two_sided_tail_probability",
          "unusual_observation"
        )
      ],
      row.names = FALSE,
      digits = 4
    )
  )
)
primary_unusual_path <- write_output_text(
  primary_unusual_lines,
  "primary_unusual_observations.txt"
)
cat(paste(primary_unusual_lines, collapse = "\n"), "\n\n")

# Publication-ready PPC figure --------------------------------------------------

specification_levels <- vapply(
  specification_definitions,
  function(definition) definition$label,
  character(1)
)
ppc_plot_data <- ppc_summary
ppc_plot_data$specification <- factor(
  ppc_plot_data$specification,
  levels = rev(specification_levels)
)
ppc_plot_data$cancer <- factor(
  ppc_plot_data$cancer,
  levels = cancer_names
)
ppc_plot_data$diagnostic_label <- factor(
  ppc_plot_data$diagnostic_label,
  levels = c(
    "Total disease burden",
    "Across-county SIR SD",
    "Moran's I of county SIR",
    "95th percentile neighboring log-SIR contrast"
  )
)
ppc_plot_data$panel <- interaction(
  ppc_plot_data$diagnostic_label,
  ppc_plot_data$cancer,
  sep = " | ",
  lex.order = TRUE
)

ppc_plot <- ggplot2::ggplot(
  ppc_plot_data,
  ggplot2::aes(y = specification)
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = predictive_lower_95,
      xend = predictive_upper_95,
      yend = specification,
      color = specification
    ),
    linewidth = 1.1
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = predictive_median,
      color = specification,
      shape = "Predictive median"
    ),
    size = 2.5
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = observed,
      shape = "Observed"
    ),
    color = "black",
    size = 2.8,
    stroke = 1
  ) +
  ggplot2::facet_wrap(
    stats::as.formula("~ diagnostic_label + cancer"),
    scales = "free_x",
    ncol = 4,
    labeller = ggplot2::label_wrap_gen(width = 28)
  ) +
  ggplot2::scale_shape_manual(
    name = NULL,
    values = c("Predictive median" = 16, "Observed" = 4)
  ) +
  ggplot2::guides(color = "none") +
  ggplot2::labs(
    x = "Diagnostic value",
    y = NULL,
    caption = paste(
      "Lines show 95% posterior predictive intervals.",
      "Moran's I uses the common binary geographic adjacency matrix;",
      "local contrast is the 95th percentile of absolute neighboring",
      "stabilized log-SIR differences."
    )
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

ppc_pdf_path <- file.path(output_dir, "ppc_three_specifications.pdf")
ppc_png_path <- file.path(output_dir, "ppc_three_specifications.png")
ggplot2::ggsave(ppc_pdf_path, ppc_plot, width = 14, height = 10)
ggplot2::ggsave(
  ppc_png_path,
  ppc_plot,
  width = 14,
  height = 10,
  dpi = 300
)

# Edge-level sensitivity comparison --------------------------------------------

boundary_probability_comparison <- data.frame(
  edge_id = rep(seq_len(number_of_edges), times = n_cancers),
  county_1 = rep(
    title_case_county(
      reference$county_names[reference$edges[, 1]]
    ),
    times = n_cancers
  ),
  county_2 = rep(
    title_case_county(
      reference$county_names[reference$edges[, 2]]
    ),
    times = n_cancers
  ),
  cancer = rep(cancer_names, each = number_of_edges),
  stringsAsFactors = FALSE
)

for (key in names(boundary_results)) {
  probability_name <- paste0("probability_", key)
  selection_name <- paste0("selected_", key)
  boundary_probability_comparison[[probability_name]] <- as.vector(
    boundary_results[[key]]$probabilities
  )
  boundary_probability_comparison[[selection_name]] <- as.vector(
    boundary_results[[key]]$selections
  )
}

classify_boundary <- function(primary, adjusted_adaptive, adjusted_fixed) {
  ifelse(
    primary == 1L & adjusted_adaptive == 1L & adjusted_fixed == 1L,
    "Selected under all three",
    ifelse(
      primary == 1L & adjusted_adaptive == 1L & adjusted_fixed == 0L,
      "Adjacency and Mean--Adjacency only",
      ifelse(
        primary == 1L & adjusted_adaptive == 0L & adjusted_fixed == 1L,
        "Adjacency and Mean only",
        ifelse(
          primary == 0L & adjusted_adaptive == 1L & adjusted_fixed == 1L,
          "Mean--Adjacency and Mean only",
          ifelse(
            primary == 1L,
            "Adjacency only",
            ifelse(
              adjusted_adaptive == 1L,
              "Mean--Adjacency only",
              ifelse(
                adjusted_fixed == 1L,
                "Mean only",
                "Not selected"
              )
            )
          )
        )
      )
    )
  )
}

boundary_probability_comparison$selection_classification <-
  classify_boundary(
    boundary_probability_comparison$selected_primary_adaptive,
    boundary_probability_comparison$selected_adjusted_adaptive,
    boundary_probability_comparison$selected_adjusted_fixed
  )

boundary_probability_output <- boundary_probability_comparison
names(boundary_probability_output)[
  match(
    c(
      "probability_primary_adaptive", "selected_primary_adaptive",
      "probability_adjusted_adaptive", "selected_adjusted_adaptive",
      "probability_adjusted_fixed", "selected_adjusted_fixed"
    ),
    names(boundary_probability_output)
  )
] <- c(
  "Adjacency_probability", "Adjacency_selected",
  "Mean--Adjacency_probability", "Mean--Adjacency_selected",
  "Mean_probability", "Mean_selected"
)
boundary_probability_path <- write_output_csv(
  boundary_probability_output,
  "boundary_probability_comparison.csv"
)

classification_levels <- c(
  "Selected under all three",
  "Adjacency and Mean--Adjacency only",
  "Adjacency and Mean only",
  "Mean--Adjacency and Mean only",
  "Adjacency only",
  "Mean--Adjacency only",
  "Mean only",
  "Not selected"
)
boundary_sensitivity_summary <- do.call(
  rbind,
  lapply(cancer_names, function(cancer) {
    classifications <- factor(
      boundary_probability_comparison$selection_classification[
        boundary_probability_comparison$cancer == cancer
      ],
      levels = classification_levels
    )
    counts <- table(classifications)
    data.frame(
      cancer = cancer,
      classification = names(counts),
      edge_count = as.integer(counts),
      proportion_of_edges = as.integer(counts) / number_of_edges,
      stringsAsFactors = FALSE
    )
  })
)
boundary_sensitivity_path <- write_output_csv(
  boundary_sensitivity_summary,
  "boundary_sensitivity_summary.csv"
)

pair_definitions <- list(
  list(
    A = "primary_adaptive",
    B = "adjusted_adaptive",
    label_A = "Adjacency",
    label_B = "Mean--Adjacency",
    comparison = "Mean adjustment effect (both adaptive adjacency)"
  ),
  list(
    A = "adjusted_adaptive",
    B = "adjusted_fixed",
    label_A = "Mean--Adjacency",
    label_B = "Mean",
    comparison = paste(
      "Adjacency-learning effect",
      "(same covariate-adjusted mean)"
    )
  ),
  list(
    A = "primary_adaptive",
    B = "adjusted_fixed",
    label_A = "Adjacency",
    label_B = "Mean",
    comparison = paste(
      "Combined specification difference",
      "(not a pure adjacency effect)"
    )
  )
)

boundary_overlap_summary <- do.call(
  rbind,
  lapply(pair_definitions, function(pair) {
    do.call(
      rbind,
      lapply(seq_len(n_cancers), function(disease) {
        selected_A <- boundary_results[[pair$A]]$selections[, disease] == 1L
        selected_B <- boundary_results[[pair$B]]$selections[, disease] == 1L
        probability_A <- boundary_results[[pair$A]]$probabilities[, disease]
        probability_B <- boundary_results[[pair$B]]$probabilities[, disease]
        intersection_size <- sum(selected_A & selected_B)
        union_size <- sum(selected_A | selected_B)
        count_A <- sum(selected_A)
        count_B <- sum(selected_B)
        data.frame(
          cancer = cancer_names[[disease]],
          specification_A = pair$label_A,
          specification_B = pair$label_B,
          comparison = pair$comparison,
          count_A = count_A,
          count_B = count_B,
          intersection_size = intersection_size,
          union_size = union_size,
          jaccard_similarity = if (
            union_size > 0L
          ) intersection_size / union_size else NA_real_,
          proportion_A_selected_in_B = if (
            count_A > 0L
          ) intersection_size / count_A else NA_real_,
          proportion_B_selected_in_A = if (
            count_B > 0L
          ) intersection_size / count_B else NA_real_,
          pearson_probability_correlation = stats::cor(
            probability_A,
            probability_B,
            method = "pearson"
          ),
          spearman_probability_correlation = stats::cor(
            probability_A,
            probability_B,
            method = "spearman"
          ),
          mean_absolute_probability_difference = mean(
            abs(probability_A - probability_B)
          ),
          root_mean_squared_probability_difference = sqrt(
            mean((probability_A - probability_B)^2)
          ),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
boundary_overlap_path <- write_output_csv(
  boundary_overlap_summary,
  "boundary_overlap_summary.csv"
)

fixed_edge_changes <- boundary_probability_comparison
fixed_edge_changes$probability_change_adjusted_adaptive_minus_fixed <-
  fixed_edge_changes$probability_adjusted_adaptive -
  fixed_edge_changes$probability_adjusted_fixed
fixed_edge_changes$absolute_probability_change <- abs(
  fixed_edge_changes$probability_change_adjusted_adaptive_minus_fixed
)
fixed_edge_changes$clean_comparison <- paste(
  "Mean--Adjacency versus Mean:",
  "same covariate-adjusted mean; adjacency learning differs."
)
fixed_edge_changes <- fixed_edge_changes[
  order(
    fixed_edge_changes$cancer,
    -fixed_edge_changes$absolute_probability_change
  ),
]
fixed_edge_changes_output <- fixed_edge_changes
names(fixed_edge_changes_output)[
  match(
    c(
      "probability_primary_adaptive", "selected_primary_adaptive",
      "probability_adjusted_adaptive", "selected_adjusted_adaptive",
      "probability_adjusted_fixed", "selected_adjusted_fixed",
      "probability_change_adjusted_adaptive_minus_fixed"
    ),
    names(fixed_edge_changes_output)
  )
] <- c(
  "Adjacency_probability", "Adjacency_selected",
  "Mean--Adjacency_probability", "Mean--Adjacency_selected",
  "Mean_probability", "Mean_selected",
  "Mean--Adjacency_minus_Mean_probability_change"
)
fixed_edge_changes_path <- write_output_csv(
  fixed_edge_changes_output,
  "fixed_adjacency_edge_changes.csv"
)

# Main sensitivity table --------------------------------------------------------

clean_overlap <- boundary_overlap_summary[
  boundary_overlap_summary$specification_A == "Mean--Adjacency" &
    boundary_overlap_summary$specification_B == "Mean",
]
main_sensitivity_table <- data.frame(
  cancer = cancer_names,
  primary_adaptive_count = as.integer(
    boundary_results$primary_adaptive$counts[cancer_names]
  ),
  adjusted_adaptive_count = as.integer(
    boundary_results$adjusted_adaptive$counts[cancer_names]
  ),
  adjusted_fixed_count = as.integer(
    boundary_results$adjusted_fixed$counts[cancer_names]
  ),
  adjusted_adaptive_vs_fixed_jaccard =
    clean_overlap$jaccard_similarity[
      match(cancer_names, clean_overlap$cancer)
    ],
  adjusted_adaptive_only = (
    clean_overlap$count_A - clean_overlap$intersection_size
  )[match(cancer_names, clean_overlap$cancer)],
  adjusted_fixed_only = (
    clean_overlap$count_B - clean_overlap$intersection_size
  )[match(cancer_names, clean_overlap$cancer)],
  stringsAsFactors = FALSE
)
main_sensitivity_output <- main_sensitivity_table
names(main_sensitivity_output) <- c(
  "cancer", "Adjacency_count", "Mean--Adjacency_count", "Mean_count",
  "Mean--Adjacency_vs_Mean_jaccard", "Mean--Adjacency_only", "Mean_only"
)
main_sensitivity_csv_path <- write_output_csv(
  main_sensitivity_output,
  "main_sensitivity_table.csv"
)

latex_rows <- vapply(
  seq_len(nrow(main_sensitivity_table)),
  function(row) {
    paste0(
      main_sensitivity_table$cancer[[row]], " & ",
      main_sensitivity_table$primary_adaptive_count[[row]], " & ",
      main_sensitivity_table$adjusted_adaptive_count[[row]], " & ",
      main_sensitivity_table$adjusted_fixed_count[[row]], " & ",
      format_number(
        main_sensitivity_table$adjusted_adaptive_vs_fixed_jaccard[[row]],
        2L
      ),
      " & ",
      main_sensitivity_table$adjusted_adaptive_only[[row]], " & ",
      main_sensitivity_table$adjusted_fixed_only[[row]],
      " \\\\"
    )
  },
  character(1)
)
latex_table <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\caption{Disease-specific sensitivity of FDR-selected boundaries. The",
  "Adjacency and Mean--Adjacency assess mean adjustment while",
  "retaining adaptive adjacency. The Jaccard and model-only columns compare",
  "Mean--Adjacency and Mean and therefore isolate adjacency learning.}",
  "\\label{tab:seer-boundary-sensitivity}",
  "\\begin{tabular}{lrrrrrr}",
  "\\hline",
  paste0(
    "Cancer & Adjacency & Mean--Adjacency & Mean & Jaccard & ",
    "Mean--Adjacency only & Mean only \\\\"
  ),
  "\\hline",
  latex_rows,
  "\\hline",
  "\\end{tabular}",
  "\\end{table}"
)
main_sensitivity_tex_path <- write_output_text(
  latex_table,
  "main_sensitivity_table.tex"
)

# Boundary-probability agreement figure ----------------------------------------

agreement_plot_rows <- list()
agreement_annotation_rows <- list()
agreement_pairs <- pair_definitions[1:2]
counter <- 1L
for (pair in agreement_pairs) {
  for (disease in seq_len(n_cancers)) {
    x <- boundary_results[[pair$A]]$probabilities[, disease]
    y <- boundary_results[[pair$B]]$probabilities[, disease]
    comparison_label <- if (pair$A == "primary_adaptive") {
      paste0(
        "Mean adjustment effect\n",
        "Adjacency vs Mean--Adjacency"
      )
    } else {
      paste0(
        "Adjacency-learning effect\n",
        "Mean--Adjacency vs Mean"
      )
    }
    agreement_plot_rows[[counter]] <- data.frame(
      probability_A = x,
      probability_B = y,
      cancer = cancer_names[[disease]],
      comparison = comparison_label,
      stringsAsFactors = FALSE
    )
    agreement_annotation_rows[[counter]] <- data.frame(
      probability_A = 0.03,
      probability_B = 0.97,
      cancer = cancer_names[[disease]],
      comparison = comparison_label,
      label = paste0(
        "Spearman = ",
        format_number(stats::cor(x, y, method = "spearman"), 2L)
      ),
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
  }
}
agreement_plot_data <- do.call(rbind, agreement_plot_rows)
agreement_annotations <- do.call(rbind, agreement_annotation_rows)
agreement_comparison_levels <- c(
  paste0(
    "Mean adjustment effect\n",
    "Adjacency vs Mean--Adjacency"
  ),
  paste0(
    "Adjacency-learning effect\n",
    "Mean--Adjacency vs Mean"
  )
)
agreement_plot_data$comparison <- factor(
  agreement_plot_data$comparison,
  levels = agreement_comparison_levels
)
agreement_annotations$comparison <- factor(
  agreement_annotations$comparison,
  levels = agreement_comparison_levels
)
agreement_plot_data$cancer <- factor(
  agreement_plot_data$cancer,
  levels = cancer_names
)
agreement_annotations$cancer <- factor(
  agreement_annotations$cancer,
  levels = cancer_names
)

agreement_plot <- ggplot2::ggplot(
  agreement_plot_data,
  ggplot2::aes(x = probability_A, y = probability_B)
) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    color = "grey50",
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    alpha = 0.65,
    size = 1.5,
    color = "#2C3E50"
  ) +
  ggplot2::geom_text(
    data = agreement_annotations,
    ggplot2::aes(label = label),
    hjust = 0,
    vjust = 1,
    size = 3,
    inherit.aes = TRUE
  ) +
  ggplot2::facet_grid(
    cancer ~ comparison,
    labeller = ggplot2::label_wrap_gen(width = 32)
  ) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    x = "Boundary probability under the first named model",
    y = "Boundary probability under the second named model",
    caption = paste0(
      "Left: covariates are added to the mean while retaining adaptive adjacency.\n",
      "Right: the covariate-adjusted mean is held fixed and only adjacency learning changes."
    )
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 9),
    plot.caption = ggplot2::element_text(hjust = 0, size = 9),
    panel.grid.minor = ggplot2::element_blank()
  )

agreement_pdf_path <- file.path(
  output_dir,
  "boundary_probability_agreement.pdf"
)
agreement_png_path <- file.path(
  output_dir,
  "boundary_probability_agreement.png"
)
ggplot2::ggsave(
  agreement_pdf_path,
  agreement_plot,
  width = 12,
  height = 11
)
ggplot2::ggsave(
  agreement_png_path,
  agreement_plot,
  width = 12,
  height = 11,
  dpi = 300
)

# Adaptive-versus-fixed boundary classification maps ---------------------------

california_map <- sf::st_as_sf(
  maps::map("county", "california", fill = TRUE, plot = FALSE)
)
sf::st_crs(california_map) <- NA
map_counties <- sub("^.*,", "", california_map$ID)
if (!identical(map_counties, reference$county_names)) {
  stop(
    "The maps::map California county order does not match county.ID; ",
    "the classification map cannot be drawn safely."
  )
}

path_coordinates <- reference$path_coordinates
path_data <- do.call(
  rbind,
  lapply(seq_along(path_coordinates), function(edge_id) {
    coordinates <- path_coordinates[[edge_id]]
    if (ncol(coordinates) < 2L || any(!is.finite(as.matrix(coordinates[, 1:2])))) {
      stopf("Stored path for edge %d is invalid.", edge_id)
    }
    data.frame(
      edge_id = edge_id,
      longitude = as.numeric(coordinates[[1]]),
      latitude = as.numeric(coordinates[[2]]),
      point_touch = nrow(coordinates) == 1L,
      path_row = seq_len(nrow(coordinates)),
      stringsAsFactors = FALSE
    )
  })
)

map_edge_rows <- list()
counter <- 1L
for (disease in seq_len(n_cancers)) {
  adaptive <- boundary_results$adjusted_adaptive$selections[, disease]
  fixed <- boundary_results$adjusted_fixed$selections[, disease]
  category <- ifelse(
    adaptive == 1L & fixed == 1L,
    "Common",
    ifelse(
      adaptive == 1L,
      "Mean--Adjacency only",
      ifelse(fixed == 1L, "Mean only", NA_character_)
    )
  )
  selected_edges <- which(!is.na(category))
  disease_paths <- path_data[path_data$edge_id %in% selected_edges,]
  disease_paths$category <- category[disease_paths$edge_id]
  disease_paths$cancer <- cancer_names[[disease]]
  map_edge_rows[[counter]] <- disease_paths
  counter <- counter + 1L
}
map_edge_data <- do.call(rbind, map_edge_rows)
map_edge_data$category <- factor(
  map_edge_data$category,
  levels = c("Common", "Mean--Adjacency only", "Mean only")
)
map_edge_data$cancer <- factor(
  map_edge_data$cancer,
  levels = cancer_names
)
line_edge_data <- map_edge_data[!map_edge_data$point_touch,]
point_edge_data <- map_edge_data[map_edge_data$point_touch,]

boundary_map_plot <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = california_map,
    fill = "grey97",
    color = "grey55",
    linewidth = 0.25
  ) +
  ggplot2::geom_path(
    data = line_edge_data,
    ggplot2::aes(
      x = longitude,
      y = latitude,
      group = interaction(cancer, edge_id),
      color = category
    ),
    linewidth = 1.1,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    data = point_edge_data,
    ggplot2::aes(
      x = longitude,
      y = latitude,
      color = category
    ),
    shape = 4,
    size = 2.4,
    stroke = 1.1
  ) +
  ggplot2::facet_wrap(
    stats::as.formula("~ cancer"),
    ncol = 2
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Common" = "#3B3B3B",
      "Mean--Adjacency only" = "#009E73",
      "Mean only" = "#0072B2"
    ),
    drop = FALSE
  ) +
  ggplot2::coord_sf(datum = NA) +
  ggplot2::labs(
    title = "Mean--Adjacency versus Mean",
    color = "FDR boundary class",
    caption = paste0(
      "Both models use the same covariate-adjusted mean; only adjacency learning differs.\n",
      "Crosses denote queen-neighbor county pairs that meet at a point rather than along a line segment."
    )
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 12,
      color = "black"
    ),
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    legend.text = ggplot2::element_text(color = "black"),
    legend.title = ggplot2::element_text(color = "black"),
    legend.position = "bottom",
    plot.caption = ggplot2::element_text(hjust = 0, color = "black")
  )

boundary_map_pdf_path <- file.path(
  output_dir,
  "adaptive_vs_fixed_boundary_maps.pdf"
)
boundary_map_png_path <- file.path(
  output_dir,
  "adaptive_vs_fixed_boundary_maps.png"
)
ggplot2::ggsave(
  boundary_map_pdf_path,
  boundary_map_plot,
  width = 10,
  height = 9
)
ggplot2::ggsave(
  boundary_map_png_path,
  boundary_map_plot,
  width = 10,
  height = 9,
  dpi = 300,
  bg = "white"
)

# Fitted relative-risk sensitivity ---------------------------------------------

fitted_risk_means <- lapply(
  ppc_results,
  function(result) result$fitted_risk_mean
)
fitted_risk_sensitivity <- do.call(
  rbind,
  lapply(pair_definitions, function(pair) {
    do.call(
      rbind,
      lapply(seq_len(n_cancers), function(disease) {
        risk_A <- fitted_risk_means[[pair$A]][, disease]
        risk_B <- fitted_risk_means[[pair$B]][, disease]
        difference <- risk_A - risk_B
        data.frame(
          cancer = cancer_names[[disease]],
          specification_A = pair$label_A,
          specification_B = pair$label_B,
          comparison = pair$comparison,
          pearson_correlation = stats::cor(
            risk_A,
            risk_B,
            method = "pearson"
          ),
          spearman_correlation = stats::cor(
            risk_A,
            risk_B,
            method = "spearman"
          ),
          mean_absolute_difference = mean(abs(difference)),
          maximum_absolute_difference = max(abs(difference)),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
fitted_risk_path <- write_output_csv(
  fitted_risk_sensitivity,
  "fitted_risk_sensitivity.csv"
)

# Data-driven manuscript paragraphs --------------------------------------------

ppc_minimum <- ppc_summary[
  which.min(ppc_summary$two_sided_ppp),
]
ppc_inside_count <- sum(ppc_summary$observed_within_95)
ppc_total_count <- nrow(ppc_summary)
primary_flagged <- unusual_observations[
  unusual_observations$specification == primary_label &
    unusual_observations$unusual_observation,
]
primary_flagged_counties <- length(unique(primary_flagged$county))
localisation_description <- if (primary_flagged_counties <= 10L) {
  "localized to a relatively small set of counties"
} else {
  "distributed across several counties rather than confined to one location"
}

primary_unusual_statement <- if (nrow(primary_flagged) == 0L) {
  paste0(
    "No county-cancer observation under Adjacency had a ",
    "descriptive two-sided mid-tail probability below ",
    unusual_support_threshold, "."
  )
} else {
  paste(
    "Under Adjacency,",
    nrow(primary_flagged),
    "county-cancer observations had descriptive two-sided mid-tail",
    "probabilities below", unusual_support_threshold, "and were",
    paste0(localisation_description, ".")
  )
}

posterior_predictive_paragraph <- paste(
  "Posterior predictive checks based on all",
  saved_specifications$primary_adaptive$retained_draws,
  "retained draws assessed total burden, across-county SIR dispersion,",
  "Moran's I, and the upper tail of neighboring stabilized log-SIR",
  "contrasts for each cancer and specification.",
  ppc_inside_count, "of", ppc_total_count,
  "observed summaries fell within their corresponding 95% posterior",
  "predictive intervals. The least well reproduced summary was",
  paste0(ppc_minimum$diagnostic_label, " for ", ppc_minimum$cancer),
  "under the", ppc_minimum$specification,
  paste0(
    "(two-sided posterior predictive p = ",
    format_number(ppc_minimum$two_sided_ppp, 3L),
    ")."
  ),
  primary_unusual_statement,
  "These checks identify unusual predictive behavior but do not establish",
  "that the model is correct or that flagged observations are formal outliers."
)

all_three_counts <- boundary_sensitivity_summary[
  boundary_sensitivity_summary$classification == "Selected under all three",
  c("cancer", "edge_count")
]
mean_effect_overlap <- boundary_overlap_summary[
  boundary_overlap_summary$specification_A == "Adjacency" &
    boundary_overlap_summary$specification_B == "Mean--Adjacency",
]
clean_overlap_ordered <- clean_overlap[
  match(cancer_names, clean_overlap$cancer),
]
most_mean_sensitive <- mean_effect_overlap$cancer[
  which.min(mean_effect_overlap$jaccard_similarity)
]
least_mean_sensitive <- mean_effect_overlap$cancer[
  which.max(mean_effect_overlap$jaccard_similarity)
]
most_adjacency_sensitive <- clean_overlap_ordered$cancer[
  which.min(clean_overlap_ordered$jaccard_similarity)
]

three_specification_paragraph <- paste(
  "The Adjacency, Mean--Adjacency, and Mean specifications selected",
  paste(
    paste0(
      cancer_names, "=",
      boundary_results$primary_adaptive$counts[cancer_names], "/",
      boundary_results$adjusted_adaptive$counts[cancer_names], "/",
      boundary_results$adjusted_fixed$counts[cancer_names]
    ),
    collapse = ", "
  ),
  "boundaries, respectively. Comparing Adjacency and Mean--Adjacency isolates the",
  "effect of adding smoking, age 65+, and poverty to the mean; this comparison",
  "was most sensitive for", most_mean_sensitive, "boundaries.",
  least_mean_sensitive, "was least sensitive to mean adjustment",
  paste0(
    "(Jaccard ",
    format_number(
      max(mean_effect_overlap$jaccard_similarity, na.rm = TRUE),
      2L
    ),
    ")."
  ),
  "Comparing Mean--Adjacency and Mean instead directly isolates whether",
  "adjacency is learned or fixed; its smallest Jaccard similarity occurred",
  "for", most_adjacency_sensitive,
  paste0(
    "(Jaccard ",
    format_number(
      min(clean_overlap_ordered$jaccard_similarity, na.rm = TRUE),
      2L
    ),
    ")."
  ),
  "Across the three models,",
  paste(
    paste0(
      all_three_counts$cancer, "=",
      all_three_counts$edge_count
    ),
    collapse = ", "
  ),
  "boundaries were selected consistently. Thus mean adjustment and adjacency",
  "learning have distinguishable effects, and the direct Adjacency-versus-Mean",
  "comparison should not be interpreted as a pure adjacency effect. These set",
  "overlaps are distinct from agreement in fitted county-level risk."
)

clean_risk <- fitted_risk_sensitivity[
  fitted_risk_sensitivity$specification_A == "Mean--Adjacency" &
    fitted_risk_sensitivity$specification_B == "Mean",
]
clean_risk <- clean_risk[match(cancer_names, clean_risk$cancer),]
largest_clean_change <- fixed_edge_changes[
  which.max(fixed_edge_changes$absolute_probability_change),
]
fixed_comparison_paragraph <- paste(
  "The cleanest fixed-adjacency comparison holds the mean structure constant:",
  "both models adjust for smoking, age 65+, and poverty, while one learns",
  "adjacency and the other fixes geographic adjacency. Their disease-specific",
  "Jaccard similarities were",
  paste0(
    paste(
      paste0(
        clean_overlap_ordered$cancer, "=",
        format_number(clean_overlap_ordered$jaccard_similarity, 2L)
      ),
      collapse = ", "
    ),
    "."
  ),
  "Mean--Adjacency uniquely selected",
  paste(
    paste0(
      clean_overlap_ordered$cancer, "=",
      clean_overlap_ordered$count_A -
        clean_overlap_ordered$intersection_size
    ),
    collapse = ", "
  ),
  "boundaries, whereas Mean uniquely selected",
  paste0(
    paste(
      paste0(
        clean_overlap_ordered$cancer, "=",
        clean_overlap_ordered$count_B -
          clean_overlap_ordered$intersection_size
      ),
      collapse = ", "
    ),
    "."
  ),
  paste0(
    "The largest edge-probability change was ",
    largest_clean_change$county_1, "--", largest_clean_change$county_2,
    " for ", largest_clean_change$cancer, " (absolute change ",
    format_number(largest_clean_change$absolute_probability_change, 2L),
    ")."
  ),
  "Posterior mean risk correlations ranged from",
  format_number(min(clean_risk$pearson_correlation), 2L), "to",
  paste0(format_number(max(clean_risk$pearson_correlation), 2L), "."),
  "Learning adjacency therefore changes which local discontinuities receive",
  "posterior support even when the broad county-level risk surfaces remain",
  "quite strongly aligned overall.",
  "This shows that boundary sensitivity is distinct from broad risk-surface",
  "agreement and does not imply uniform superiority of either adjacency model."
)

manuscript_summary_lines <- c(
  "Posterior predictive paragraph",
  "------------------------------",
  posterior_predictive_paragraph,
  "",
  "Three-specification sensitivity paragraph",
  "-----------------------------------------",
  three_specification_paragraph,
  "",
  "Fixed-adjacency comparison paragraph",
  "------------------------------------",
  fixed_comparison_paragraph
)
manuscript_summary_path <- write_output_text(
  manuscript_summary_lines,
  "manuscript_summary.txt"
)

# Final report ------------------------------------------------------------------

boundary_count_report <- vapply(
  names(boundary_results),
  function(key) {
    paste0(
      specification_definitions[[key]]$short_label,
      ": ",
      paste(
        paste0(
          cancer_names,
          "=",
          boundary_results[[key]]$counts[cancer_names]
        ),
        collapse = ", "
      )
    )
  },
  character(1)
)

principal_jaccards <- clean_overlap_ordered[
  ,
  c("cancer", "jaccard_similarity")
]
number_flagged <- sum(unusual_observations$unusual_observation)
draws_used <- unique(vapply(ppc_results, function(x) x$draws_used, integer(1)))
if (length(draws_used) != 1L) {
  stop("The three specifications used different numbers of PPC draws.")
}

generated_files <- sort(
  basename(
    c(
      validation_report_path,
      ppc_summary_path,
      unusual_path,
      primary_unusual_path,
      ppc_pdf_path,
      ppc_png_path,
      boundary_probability_path,
      boundary_sensitivity_path,
      boundary_overlap_path,
      fixed_edge_changes_path,
      main_sensitivity_csv_path,
      main_sensitivity_tex_path,
      agreement_pdf_path,
      agreement_png_path,
      boundary_map_pdf_path,
      boundary_map_png_path,
      fitted_risk_path,
      manuscript_summary_path,
      file.path(output_dir, "analysis_report.txt")
    )
  )
)

analysis_report_lines <- c(
  "Posterior predictive and sensitivity analysis report",
  "====================================================",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Posterior predictive draws used per specification:", draws_used),
  paste(
    "Adjacency manuscript boundary counts reproduced:",
    primary_counts_reproduced
  ),
  paste(
    "Posterior predictive checks:",
    paste(
      c(
        "total disease burden",
        "across-county SIR SD",
        "Moran's I on SIR using binary geographic adjacency",
        "95th percentile neighboring stabilized log-SIR contrast"
      ),
      collapse = "; "
    )
  ),
  paste(
    "Unusual county-cancer observations flagged across all specifications:",
    number_flagged,
    paste0(
      "(descriptive two-sided mid-tail probability < ",
      unusual_support_threshold,
      ")"
    )
  ),
  "",
  "FDR-selected boundary counts:",
  paste0("  ", boundary_count_report),
  "",
  paste(
    "Mean--Adjacency versus Mean Jaccard similarities:",
    paste(
      paste0(
        principal_jaccards$cancer,
        "=",
        format_number(principal_jaccards$jaccard_similarity, 3L)
      ),
      collapse = ", "
    )
  ),
  "",
  "Comparison interpretation:",
  paste0(
    "  Mean adjustment effect: results_adj.RData versus ",
    "results_meanadj.RData (both adaptive adjacency)."
  ),
  paste0(
    "  Adjacency-learning effect: results_meanadj.RData versus ",
    "results_mean.RData (same covariate-adjusted mean)."
  ),
  paste0(
    "  results_adj.RData versus results_mean.RData differs in both the mean ",
    "and adjacency specifications and is not a pure adjacency comparison."
  ),
  "",
  "Warnings and interpretation limits:",
  paste0(
    "  The unusual-observation threshold is descriptive and unadjusted; ",
    "flags are not formal outlier declarations."
  ),
  paste0(
    "  Queen-neighbor pairs meeting only at a point are shown as crosses ",
    "in the boundary classification maps."
  ),
  "  Posterior predictive p-values are diagnostics, not proof of model correctness.",
  "",
  "Generated files:",
  paste0("  ", generated_files),
  "",
  paste("Manuscript-ready summaries:", manuscript_summary_path)
)
analysis_report_path <- write_output_text(
  analysis_report_lines,
  "analysis_report.txt"
)

cat(paste(analysis_report_lines, collapse = "\n"), "\n")
cat("\nAnalysis report saved to:", analysis_report_path, "\n")
