# Standardize reader-facing specification labels in completed PPC/sensitivity
# outputs. This script changes labels and redraws figures from existing CSV
# summaries; it performs no posterior simulation and calls no MCMC sampler.

rm(list = ls())

if (!all(file.exists(c("README.md", "data analysis")))) {
  stop("Run this script from the repository root.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required.", call. = FALSE)
}
if (!requireNamespace("sf", quietly = TRUE)) {
  stop("The sf package is required to redraw the completed boundary map.", call. = FALSE)
}
if (!requireNamespace("maps", quietly = TRUE)) {
  stop("The maps package is required to redraw the completed boundary map.", call. = FALSE)
}

output_dir <- file.path(
  "data analysis",
  "posterior_predictive_sensitivity_output"
)
csv_paths <- list.files(
  output_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)
numeric_csv_snapshot <- setNames(
  lapply(csv_paths, function(path) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    unname(x[vapply(x, is.numeric, logical(1))])
  }),
  basename(csv_paths)
)

label_replacements <- c(
  "Adaptive adjacency, intercept-only mean" = "Adjacency",
  "Primary adaptive" = "Adjacency",
  "Adaptive adjacency, covariate-adjusted mean" = "Mean--Adjacency",
  "Adjusted adaptive" = "Mean--Adjacency",
  "Fixed adjacency, covariate-adjusted mean" = "Mean",
  "Adjusted fixed" = "Mean"
)
classification_replacements <- c(
  "Both adaptive only" = "Adjacency and Mean--Adjacency only",
  "Primary and fixed only" = "Adjacency and Mean only",
  "Both covariate-adjusted only" = "Mean--Adjacency and Mean only",
  "Primary only" = "Adjacency only",
  "Adjusted adaptive only" = "Mean--Adjacency only",
  "Adjusted fixed only" = "Mean only"
)

replace_exact <- function(x, replacements) {
  x <- as.character(x)
  matched <- match(x, names(replacements))
  x[!is.na(matched)] <- unname(replacements[matched[!is.na(matched)]])
  x
}

audit_lines <- c(
  "Specification label audit",
  "=========================",
  "",
  "Canonical labels: Adjacency; Mean--Adjacency; Mean.",
  "Presentation-only update: no MCMC, posterior predictive simulation,",
  "posterior probability, FDR selection, fitted-risk, or diagnostic statistic",
  "was recomputed.",
  ""
)

update_csv_columns <- function(filename, columns, replacements) {
  path <- file.path(output_dir, filename)
  if (!file.exists(path)) {
    stop("Missing generated summary: ", path, call. = FALSE)
  }
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  changed_expressions <- character()
  for (column in columns) {
    if (!column %in% names(x)) {
      stop(filename, " lacks column ", column, ".", call. = FALSE)
    }
    before <- as.character(x[[column]])
    after <- replace_exact(before, replacements)
    if (!identical(before, after)) {
      matched_old <- unique(before[before != after])
      changed_expressions <- c(
        changed_expressions,
        paste0(
          column,
          ": ",
          matched_old,
          " -> ",
          replace_exact(matched_old, replacements)
        )
      )
      x[[column]] <- after
    }
  }
  utils::write.csv(x, path, row.names = FALSE, na = "")
  if (length(changed_expressions) > 0L) {
    audit_lines <<- c(
      audit_lines,
      paste0(filename, ":"),
      paste0("  ", changed_expressions),
      ""
    )
  }
  invisible(x)
}

rename_csv_columns <- function(filename, replacements) {
  path <- file.path(output_dir, filename)
  if (!file.exists(path)) {
    stop("Missing generated summary: ", path, call. = FALSE)
  }
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  before_names <- names(x)
  numeric_before <- x[vapply(x, is.numeric, logical(1))]
  for (old in names(replacements)) {
    if (old %in% names(x)) {
      new <- unname(replacements[[old]])
      if (new %in% names(x)) {
        stop(filename, " contains both old and canonical column names: ", old,
             " and ", new, ".", call. = FALSE)
      }
      names(x)[names(x) == old] <- new
    }
  }
  if (!identical(before_names, names(x))) {
    utils::write.csv(x, path, row.names = FALSE, na = "")
    reread <- utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    numeric_after <- reread[vapply(reread, is.numeric, logical(1))]
    if (!isTRUE(all.equal(
      unname(as.data.frame(numeric_before)),
      unname(as.data.frame(numeric_after)),
      tolerance = 0,
      check.attributes = FALSE
    ))) {
      stop("Numerical values changed while renaming columns in ", filename, ".",
           call. = FALSE)
    }
    audit_lines <<- c(
      audit_lines,
      paste0(filename, ":"),
      "  reader-facing column names standardized; numerical columns verified unchanged.",
      ""
    )
  }
  invisible(x)
}

edge_column_replacements <- c(
  "probability_primary_adaptive" = "Adjacency_probability",
  "selected_primary_adaptive" = "Adjacency_selected",
  "probability_adjusted_adaptive" = "Mean--Adjacency_probability",
  "selected_adjusted_adaptive" = "Mean--Adjacency_selected",
  "probability_adjusted_fixed" = "Mean_probability",
  "selected_adjusted_fixed" = "Mean_selected"
)
rename_csv_columns(
  "boundary_probability_comparison.csv",
  edge_column_replacements
)
rename_csv_columns(
  "fixed_adjacency_edge_changes.csv",
  c(
    edge_column_replacements,
    "probability_change_adjusted_adaptive_minus_fixed" =
      "Mean--Adjacency_minus_Mean_probability_change"
  )
)

ppc_summary <- update_csv_columns(
  "ppc_summary.csv",
  "specification",
  label_replacements
)
update_csv_columns(
  "ppc_unusual_county_cancer.csv",
  "specification",
  label_replacements
)
boundary_overlap <- update_csv_columns(
  "boundary_overlap_summary.csv",
  c("specification_A", "specification_B"),
  label_replacements
)
update_csv_columns(
  "fitted_risk_sensitivity.csv",
  c("specification_A", "specification_B"),
  label_replacements
)
update_csv_columns(
  "boundary_probability_comparison.csv",
  "selection_classification",
  classification_replacements
)
update_csv_columns(
  "boundary_sensitivity_summary.csv",
  "classification",
  classification_replacements
)
fixed_changes <- update_csv_columns(
  "fixed_adjacency_edge_changes.csv",
  "selection_classification",
  classification_replacements
)
fixed_path <- file.path(output_dir, "fixed_adjacency_edge_changes.csv")
old_clean_comparison <- paste(
  "Adjusted adaptive versus adjusted fixed:",
  "same covariate-adjusted mean; adjacency learning differs."
)
clean_comparison_changed <-
  fixed_changes$clean_comparison == old_clean_comparison
if (any(clean_comparison_changed)) {
  fixed_changes$clean_comparison[clean_comparison_changed] <- paste(
    "Mean--Adjacency versus Mean:",
    "same covariate-adjusted mean; adjacency learning differs."
  )
}
utils::write.csv(fixed_changes, fixed_path, row.names = FALSE, na = "")
if (any(clean_comparison_changed)) {
  audit_lines <- c(
    audit_lines,
    "fixed_adjacency_edge_changes.csv:",
    paste0(
      "  clean_comparison: Adjusted adaptive versus adjusted fixed -> ",
      "Mean--Adjacency versus Mean"
    ),
    ""
  )
}

main_table_path <- file.path(output_dir, "main_sensitivity_table.csv")
main_table <- utils::read.csv(
  main_table_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
expected_main_names <- c(
  "cancer", "primary_adaptive_count", "adjusted_adaptive_count",
  "adjusted_fixed_count", "adjusted_adaptive_vs_fixed_jaccard",
  "adjusted_adaptive_only", "adjusted_fixed_only"
)
if (identical(names(main_table), expected_main_names)) {
  names(main_table) <- c(
    "cancer", "Adjacency_count", "Mean--Adjacency_count", "Mean_count",
    "Mean--Adjacency_vs_Mean_jaccard", "Mean--Adjacency_only", "Mean_only"
  )
  utils::write.csv(main_table, main_table_path, row.names = FALSE, na = "")
  audit_lines <- c(
    audit_lines,
    "main_sensitivity_table.csv:",
    "  reader-facing column labels standardized to Adjacency, Mean--Adjacency, and Mean.",
    ""
  )
}

text_replacements <- c(
  "covariate-Mean--Adjacency" = "Mean--Adjacency",
  "covariate-Mean" = "Mean",
  label_replacements,
  "covariate-adjusted adaptive model" = "Mean--Adjacency specification",
  "covariate-adjusted fixed model" = "Mean specification",
  "primary specification" = "Adjacency specification",
  "primary model" = "Adjacency specification",
  "primary adaptive" = "Adjacency",
  "adjusted adaptive" = "Mean--Adjacency",
  "adjusted fixed" = "Mean",
  "primary and adjusted-adaptive columns" =
    "Adjacency and Mean--Adjacency",
  "Primary manuscript boundary counts" = "Adjacency manuscript boundary counts",
  "Primary manuscript counts" = "Adjacency manuscript counts",
  "the two adaptive models" = "Adjacency and Mean--Adjacency",
  "the two covariate-adjusted models" = "Mean--Adjacency and Mean",
  "direct primary-versus-fixed" = "direct Adjacency-versus-Mean",
  "The adaptive model uniquely selected" = "Mean--Adjacency uniquely selected",
  "the fixed model uniquely selected" = "Mean uniquely selected",
  "Mean--Adj. only" = "Mean--Adjacency only",
  "Mean--Adj." = "Mean--Adjacency",
  "Cancer & Primary & Adj. adapt. & Adj. fixed & Jaccard & Adapt. only & Fixed only" =
    "Cancer & Adjacency & Mean--Adj. & Mean & Jaccard & Mean--Adj. only & Mean only"
)
text_files <- c(
  "analysis_report.txt",
  "manuscript_summary.txt",
  "primary_unusual_observations.txt",
  "validation_report.txt",
  "main_sensitivity_table.tex",
  "main_manuscript_insertions.tex",
  "supplement_manuscript_insertions.tex"
)
for (filename in text_files) {
  path <- file.path(output_dir, filename)
  if (!file.exists(path)) {
    next
  }
  lines <- readLines(path, warn = FALSE)
  before <- lines
  changed_pairs <- character()
  for (old in names(text_replacements)) {
    if (any(grepl(old, lines, fixed = TRUE))) {
      lines <- gsub(old, text_replacements[[old]], lines, fixed = TRUE)
      changed_pairs <- c(
        changed_pairs,
        paste0(old, " -> ", text_replacements[[old]])
      )
    }
  }
  if (!identical(before, lines)) {
    writeLines(lines, path, useBytes = TRUE)
    audit_lines <- c(
      audit_lines,
      paste0(filename, ":"),
      paste0("  ", unique(changed_pairs)),
      ""
    )
  }
}

# Redraw the global PPC figure from its completed 48-row summary. No replicated
# counts or diagnostic statistics are recomputed.
specification_levels <- c("Adjacency", "Mean--Adjacency", "Mean")
cancer_levels <- c("Lung", "Esophageal", "Larynx", "Colorectal")
ppc_plot_data <- ppc_summary
ppc_plot_data$specification <- factor(
  ppc_plot_data$specification,
  levels = rev(specification_levels)
)
ppc_plot_data$cancer <- factor(ppc_plot_data$cancer, levels = cancer_levels)
ppc_plot_data$diagnostic_label <- factor(
  ppc_plot_data$diagnostic_label,
  levels = c(
    "Total disease burden",
    "Across-county SIR SD",
    "Moran's I of county SIR",
    "95th percentile neighboring log-SIR contrast"
  )
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
    ggplot2::aes(x = observed, shape = "Observed"),
    color = "black",
    size = 2.8,
    stroke = 1
  ) +
  ggplot2::facet_wrap(
    ~ diagnostic_label + cancer,
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
ggplot2::ggsave(
  file.path(output_dir, "ppc_three_specifications.pdf"),
  ppc_plot,
  width = 14,
  height = 10
)
ggplot2::ggsave(
  file.path(output_dir, "ppc_three_specifications.png"),
  ppc_plot,
  width = 14,
  height = 10,
  dpi = 300
)
audit_lines <- c(
  audit_lines,
  "ppc_three_specifications.pdf/.png:",
  "  figure labels redrawn from ppc_summary.csv using the canonical labels.",
  ""
)

# Redraw edge-probability agreement from the existing edge summary.
edge_data <- utils::read.csv(
  file.path(output_dir, "boundary_probability_comparison.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
agreement_rows <- list()
annotation_rows <- list()
counter <- 1L
comparison_definitions <- list(
  list(
    x = "Adjacency_probability",
    y = "Mean--Adjacency_probability",
    specification_A = "Adjacency",
    specification_B = "Mean--Adjacency",
    label = "Mean adjustment effect\nAdjacency vs Mean--Adjacency"
  ),
  list(
    x = "Mean--Adjacency_probability",
    y = "Mean_probability",
    specification_A = "Mean--Adjacency",
    specification_B = "Mean",
    label = "Adjacency-learning effect\nMean--Adjacency vs Mean"
  )
)
for (comparison in comparison_definitions) {
  for (cancer in cancer_levels) {
    x <- edge_data[edge_data$cancer == cancer, comparison$x]
    y <- edge_data[edge_data$cancer == cancer, comparison$y]
    saved_summary <- boundary_overlap[
      boundary_overlap$cancer == cancer &
        boundary_overlap$specification_A == comparison$specification_A &
        boundary_overlap$specification_B == comparison$specification_B,
    ]
    if (nrow(saved_summary) != 1L) {
      stop(
        "Missing saved agreement summary for ", cancer, ": ",
        comparison$specification_A, " versus ", comparison$specification_B,
        call. = FALSE
      )
    }
    agreement_rows[[counter]] <- data.frame(
      probability_A = x,
      probability_B = y,
      cancer = cancer,
      comparison = comparison$label,
      stringsAsFactors = FALSE
    )
    annotation_rows[[counter]] <- data.frame(
      probability_A = 0.03,
      probability_B = 0.97,
      cancer = cancer,
      comparison = comparison$label,
      label = sprintf(
        "Spearman = %.2f",
        saved_summary$spearman_probability_correlation
      ),
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
  }
}
agreement_data <- do.call(rbind, agreement_rows)
agreement_annotations <- do.call(rbind, annotation_rows)
comparison_levels <- vapply(
  comparison_definitions,
  function(x) x[["label"]],
  character(1)
)
agreement_data$comparison <- factor(
  agreement_data$comparison,
  levels = comparison_levels
)
agreement_annotations$comparison <- factor(
  agreement_annotations$comparison,
  levels = comparison_levels
)
agreement_data$cancer <- factor(agreement_data$cancer, levels = cancer_levels)
agreement_annotations$cancer <- factor(
  agreement_annotations$cancer,
  levels = cancer_levels
)
agreement_plot <- ggplot2::ggplot(
  agreement_data,
  ggplot2::aes(x = probability_A, y = probability_B)
) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    color = "grey50",
    linetype = "dashed"
  ) +
  ggplot2::geom_point(alpha = 0.65, size = 1.5, color = "#2C3E50") +
  ggplot2::geom_text(
    data = agreement_annotations,
    ggplot2::aes(label = label),
    hjust = 0,
    vjust = 1,
    size = 3
  ) +
  ggplot2::facet_grid(
    cancer ~ comparison,
    labeller = ggplot2::label_wrap_gen(width = 32)
  ) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    x = "Boundary probability under the first named specification",
    y = "Boundary probability under the second named specification"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 9),
    panel.grid.minor = ggplot2::element_blank()
  )
ggplot2::ggsave(
  file.path(output_dir, "boundary_probability_agreement.pdf"),
  agreement_plot,
  width = 8.5,
  height = 11
)
ggplot2::ggsave(
  file.path(output_dir, "boundary_probability_agreement.png"),
  agreement_plot,
  width = 8.5,
  height = 11,
  dpi = 300
)
audit_lines <- c(
  audit_lines,
  "boundary_probability_agreement.pdf/.png:",
  "  panel labels redrawn from boundary_probability_comparison.csv using the canonical labels.",
  "",
  "main_manuscript_insertions.tex:",
  "  formal specification names standardized to Adjacency, Mean--Adjacency, and Mean.",
  "",
  "supplement_manuscript_insertions.tex:",
  "  formal specification names, table entries, and figure captions standardized.",
  ""
)

# Redraw the Mean--Adjacency-versus-Mean boundary map from the completed edge
# selections and stored county-border paths. This changes presentation only:
# posterior probabilities and FDR selections are read from existing outputs.
map_source_path <- file.path("data analysis", "results_adj.RData")
if (!file.exists(map_source_path)) {
  stop("Missing stored map geometry source: ", map_source_path, call. = FALSE)
}
map_environment <- new.env(parent = emptyenv())
load(map_source_path, envir = map_environment)
required_map_objects <- c("path", "county.ID")
missing_map_objects <- required_map_objects[
  !vapply(required_map_objects, exists, logical(1), envir = map_environment)
]
if (length(missing_map_objects) > 0L) {
  stop(
    "The stored Adjacency result lacks map object(s): ",
    paste(missing_map_objects, collapse = ", "),
    call. = FALSE
  )
}
path_coordinates <- get("path", envir = map_environment)
county_names <- as.character(get("county.ID", envir = map_environment))
california_map <- sf::st_as_sf(
  maps::map("county", "california", fill = TRUE, plot = FALSE)
)
sf::st_crs(california_map) <- NA
map_counties <- sub("^.*,", "", california_map$ID)
if (!identical(map_counties, county_names)) {
  stop(
    "The maps::map California county order does not match county.ID; ",
    "the completed boundary selections cannot be mapped safely.",
    call. = FALSE
  )
}
if (length(path_coordinates) != length(unique(edge_data$edge_id))) {
  stop("Stored path and completed edge-selection counts differ.", call. = FALSE)
}
path_data <- do.call(
  rbind,
  lapply(seq_along(path_coordinates), function(edge_id) {
    coordinates <- path_coordinates[[edge_id]]
    if (ncol(coordinates) < 2L ||
        any(!is.finite(as.matrix(coordinates[, 1:2])))) {
      stop("Stored path for edge ", edge_id, " is invalid.", call. = FALSE)
    }
    data.frame(
      edge_id = edge_id,
      longitude = as.numeric(coordinates[[1]]),
      latitude = as.numeric(coordinates[[2]]),
      point_touch = nrow(coordinates) == 1L,
      stringsAsFactors = FALSE
    )
  })
)
map_edge_rows <- list()
for (cancer in cancer_levels) {
  cancer_edges <- edge_data[edge_data$cancer == cancer, ]
  cancer_edges <- cancer_edges[order(cancer_edges$edge_id), ]
  if (!identical(as.integer(cancer_edges$edge_id), seq_along(path_coordinates))) {
    stop("Edge ordering is incomplete for ", cancer, ".", call. = FALSE)
  }
  mean_adjacency <- as.integer(cancer_edges[["Mean--Adjacency_selected"]])
  mean_fixed <- as.integer(cancer_edges[["Mean_selected"]])
  category <- ifelse(
    mean_adjacency == 1L & mean_fixed == 1L,
    "Common",
    ifelse(
      mean_adjacency == 1L,
      "Mean--Adjacency only",
      ifelse(mean_fixed == 1L, "Mean only", NA_character_)
    )
  )
  selected_edges <- which(!is.na(category))
  disease_paths <- path_data[path_data$edge_id %in% selected_edges, ]
  disease_paths$category <- category[disease_paths$edge_id]
  disease_paths$cancer <- cancer
  map_edge_rows[[cancer]] <- disease_paths
}
map_edge_data <- do.call(rbind, map_edge_rows)
map_edge_data$category <- factor(
  map_edge_data$category,
  levels = c("Common", "Mean--Adjacency only", "Mean only")
)
map_edge_data$cancer <- factor(map_edge_data$cancer, levels = cancer_levels)
line_edge_data <- map_edge_data[!map_edge_data$point_touch, ]
point_edge_data <- map_edge_data[map_edge_data$point_touch, ]
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
    ggplot2::aes(x = longitude, y = latitude, color = category),
    shape = 4,
    size = 2.4,
    stroke = 1.1
  ) +
  ggplot2::facet_wrap(stats::as.formula("~ cancer"), ncol = 2) +
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
      "Both specifications use the same covariate-adjusted mean; only adjacency learning differs.\n",
      "Crosses denote queen-neighbor county pairs that meet at a point rather than along a line segment."
    )
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.text = ggplot2::element_text(face = "bold", size = 12),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    plot.caption = ggplot2::element_text(hjust = 0)
  )
ggplot2::ggsave(
  file.path(output_dir, "adaptive_vs_fixed_boundary_maps.pdf"),
  boundary_map_plot,
  width = 10,
  height = 9
)
ggplot2::ggsave(
  file.path(output_dir, "adaptive_vs_fixed_boundary_maps.png"),
  boundary_map_plot,
  width = 10,
  height = 9,
  dpi = 300,
  bg = "white"
)
audit_lines <- c(
  audit_lines,
  "adaptive_vs_fixed_boundary_maps.pdf/.png:",
  paste(
    "  title and legend redrawn from completed selections using",
    "Mean--Adjacency versus Mean; Common; Mean--Adjacency only; Mean only."
  ),
  ""
)

for (path in csv_paths) {
  filename <- basename(path)
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  numeric_after <- unname(x[vapply(x, is.numeric, logical(1))])
  if (!isTRUE(all.equal(
    numeric_csv_snapshot[[filename]],
    numeric_after,
    tolerance = 0,
    check.attributes = FALSE
  ))) {
    stop(
      "Numerical CSV content changed during presentation update: ",
      filename,
      call. = FALSE
    )
  }
}
audit_lines <- c(
  audit_lines,
  "Numerical preservation check:",
  paste0(
    "  PASS: all numeric columns in ", length(csv_paths),
    " existing CSV outputs were unchanged during this presentation-only run."
  ),
  ""
)

writeLines(
  audit_lines,
  file.path(output_dir, "specification_label_audit.txt"),
  useBytes = TRUE
)
message(
  "Specification labels standardized without rerunning MCMC or posterior ",
  "predictive simulation."
)
