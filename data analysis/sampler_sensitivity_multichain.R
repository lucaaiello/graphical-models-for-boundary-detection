rm(list = ls())

# SEER sampler sensitivity across three model specifications
# =============================================================================
# This script does not run MCMC. For each of Adjacency, Mean--Adjacency, and
# Mean, it compares the completed adaptive-tempering analysis with
# the completed independent non-tempered multiple-chain analysis.
#
# Run from the repository root, after post-processing both workflows for all
# three specifications, with:
#   source("data analysis/sampler_sensitivity_multichain.R")
#
# Inputs are located automatically under:
#   data analysis/multiple_chains_tempering_output/
#     <specification>/hotter_final/postprocessing/
#   data analysis/multiple_chains_output/
#     <specification>/production/postprocessing/
#
# All generated artifacts are written below:
#   data analysis/sampler_sensitivity/

options(stringsAsFactors = FALSE)

tempering_root <- Sys.getenv(
  "SEER_SAMPLER_SENSITIVITY_TEMPERING_ROOT",
  unset = file.path(
    "data analysis", "multiple_chains_tempering_output"
  )
)
non_tempered_root <- Sys.getenv(
  "SEER_SAMPLER_SENSITIVITY_NON_TEMPERED_ROOT",
  unset = file.path("data analysis", "multiple_chains_output")
)
tempering_mode <- Sys.getenv(
  "SEER_SAMPLER_SENSITIVITY_TEMPERING_MODE", unset = "hotter_final"
)
non_tempered_mode <- Sys.getenv(
  "SEER_SAMPLER_SENSITIVITY_NON_TEMPERED_MODE", unset = "production"
)
output_dir <- Sys.getenv(
  "SEER_SAMPLER_SENSITIVITY_OUTPUT",
  unset = file.path("data analysis", "sampler_sensitivity")
)
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
for (directory in c(output_dir, tables_dir, figures_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Run install.R first.", call. = FALSE
  )
}

specifications <- data.frame(
  specification = c("adj", "meanadj", "mean"),
  label = c("Adjacency", "Mean--Adjacency", "Mean"),
  order = seq_len(3L),
  stringsAsFactors = FALSE
)
disease_names <- c("Lung", "Esophageal", "Larynx", "Colorectal")

stopf <- function(...) stop(sprintf(...), call. = FALSE)
write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
write_tex <- function(lines, path) {
  writeLines(lines, path, useBytes = TRUE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
normalized_path <- function(path) {
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
jaccard_index <- function(a, b) {
  denominator <- length(union(a, b))
  if (!denominator) return(1)
  length(intersect(a, b)) / denominator
}

input_paths <- list()
for (index in seq_len(nrow(specifications))) {
  specification <- specifications$specification[[index]]
  tempering_postprocessing <- file.path(
    tempering_root, specification, tempering_mode, "postprocessing"
  )
  non_tempered_postprocessing <- file.path(
    non_tempered_root, specification, non_tempered_mode, "postprocessing"
  )
  input_paths[[specification]] <- list(
    tempering_diagnostic_overview = file.path(
      tempering_postprocessing, "diagnostics",
      paste0("diagnostic_overview_", specification, ".csv")
    ),
    non_tempered_diagnostic_overview = file.path(
      non_tempered_postprocessing, "diagnostics",
      paste0("diagnostic_overview_", specification, ".csv")
    ),
    tempering_diagnostics = file.path(
      tempering_postprocessing,
      "tempering_convergence_diagnostics.rds"
    ),
    non_tempered_diagnostics = file.path(
      non_tempered_postprocessing,
      "multiple_chains_diagnostics.rds"
    ),
    tempering_boundaries = file.path(
      tempering_postprocessing,
      "tempering_fdr_boundary_results.rds"
    ),
    non_tempered_boundaries = file.path(
      non_tempered_postprocessing,
      "multiple_chains_boundary_summaries.rds"
    )
  )
}

input_manifest <- do.call(rbind, lapply(
  names(input_paths),
  function(specification) {
    paths <- input_paths[[specification]]
    do.call(rbind, lapply(names(paths), function(input_name) {
      path <- paths[[input_name]]
      data.frame(
        specification = specification,
        input = input_name,
        file = path,
        exists = file.exists(path),
        bytes = if (file.exists(path)) file.info(path)$size else NA_real_,
        md5 = if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_,
        stringsAsFactors = FALSE
      )
    }))
  }
))
missing_inputs <- input_manifest$file[!input_manifest$exists]
if (length(missing_inputs)) {
  stop(
    "Sampler sensitivity requires completed post-processing for both ",
    "workflows and all three specifications. Missing:\n",
    paste(" -", missing_inputs, collapse = "\n"),
    "\nRun the missing sampler/post-processing workflow and try again.",
    call. = FALSE
  )
}
input_manifest$file <- vapply(
  input_manifest$file, normalized_path, character(1L)
)
diagnostic_columns <- c(
  "group", "variables", "median_rhat", "max_rhat", "rhat_above_1.01",
  "rhat_above_1.05", "rhat_above_1.10", "min_bulk_ess",
  "median_bulk_ess", "min_tail_ess", "median_tail_ess",
  "max_mcse_mean", "max_relative_mcse_mean"
)
quantity_labels <- c(
  "Core parameters" = "Core parameters",
  "Centered regression coefficients" = "Centered regression",
  "Raw covariance-factor coordinates" = "Raw covariance coordinates",
  "Disease covariance and correlation" = "Covariance/correlation",
  "Adjacency cardinality" = "Adjacency cardinality",
  "Draw-level boundary counts" = "Boundary count",
  "Centered spatial effects" = "Centered spatial effects",
  "Fitted log-risks (equivalently, relative risks)" = "Fitted log-risks",
  "Latent Gaussian r variables" = "Latent Gaussian variables"
)

read_diagnostic_overview <- function(path, specification, workflow) {
  result <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  missing <- setdiff(diagnostic_columns, names(result))
  if (length(missing)) {
    stopf("%s is missing diagnostic columns: %s", path,
          paste(missing, collapse = ", "))
  }
  order <- match(names(quantity_labels), result$group)
  if (anyNA(order)) {
    stopf("%s does not contain every required diagnostic group.", path)
  }
  result <- result[order, diagnostic_columns, drop = FALSE]
  # Mean uses a fixed geographic adjacency matrix. Its adjacency cardinality
  # is therefore constant by construction, so Rhat and ESS are not defined and
  # should not appear as infinite diagnostic values in the comparison table.
  if (identical(specification, "mean")) {
    result <- result[
      result$group != "Adjacency cardinality", , drop = FALSE
    ]
  }
  result$quantity <- unname(quantity_labels[result$group])
  result$specification <- specification
  result$specification_label <- specifications$label[
    match(specification, specifications$specification)
  ]
  result$workflow <- workflow
  result[, c(
    "specification", "specification_label", "workflow", "quantity",
    diagnostic_columns
  )]
}

convergence_comparison <- do.call(rbind, lapply(
  seq_len(nrow(specifications)),
  function(index) {
    specification <- specifications$specification[[index]]
    paths <- input_paths[[specification]]
    rbind(
      read_diagnostic_overview(
        paths$tempering_diagnostic_overview,
        specification, "Adaptive tempering"
      ),
      read_diagnostic_overview(
        paths$non_tempered_diagnostic_overview,
        specification, "Non-tempered"
      )
    )
  }
))
row.names(convergence_comparison) <- NULL

validate_result_specification <- function(x, specification, path) {
  saved <- x$model_specification
  if (is.null(saved)) saved <- "adj"
  if (!identical(saved, specification)) {
    stopf("%s contains specification '%s', not '%s'.",
          path, saved, specification)
  }
  invisible(x)
}

boundary_comparison <- list()
edge_probability_comparison <- list()
fitted_log_risk_comparison <- list()
fitted_log_risk_values <- list()

for (index in seq_len(nrow(specifications))) {
  specification <- specifications$specification[[index]]
  specification_label <- specifications$label[[index]]
  paths <- input_paths[[specification]]

  tempering_boundaries <- readRDS(paths$tempering_boundaries)
  non_tempered_boundaries <- readRDS(paths$non_tempered_boundaries)
  validate_result_specification(
    tempering_boundaries, specification, paths$tempering_boundaries
  )
  validate_result_specification(
    non_tempered_boundaries, specification, paths$non_tempered_boundaries
  )
  tempering_probability <-
    tempering_boundaries$pooled_boundary_probabilities
  non_tempered_probability <-
    non_tempered_boundaries$pooled_boundary_probabilities
  if (is.null(rownames(tempering_probability)) ||
      is.null(rownames(non_tempered_probability))) {
    stopf("Boundary probability matrices for %s require edge names.",
          specification)
  }
  non_tempered_order <- match(
    rownames(tempering_probability), rownames(non_tempered_probability)
  )
  if (anyNA(non_tempered_order)) {
    stopf("Could not align all county borders for %s.", specification)
  }
  non_tempered_probability <- non_tempered_probability[non_tempered_order, , drop = FALSE]
  if (!identical(dim(tempering_probability), dim(non_tempered_probability))) {
    stopf("Boundary probability dimensions differ for %s.", specification)
  }

  tempering_sets <- tempering_boundaries$pooled_selected_boundary_sets
  non_tempered_sets <- non_tempered_boundaries$pooled_selected_boundary_sets
  if (length(tempering_sets) != length(disease_names) ||
      length(non_tempered_sets) != length(disease_names)) {
    stopf("Selected boundary sets are incomplete for %s.", specification)
  }

  for (disease in seq_along(disease_names)) {
    probability_difference <- non_tempered_probability[, disease] -
      tempering_probability[, disease]
    boundary_comparison[[length(boundary_comparison) + 1L]] <- data.frame(
      specification = specification,
      specification_label = specification_label,
      cancer = disease_names[[disease]],
      adaptive_tempering_count = length(tempering_sets[[disease]]),
      non_tempered_count = length(non_tempered_sets[[disease]]),
      jaccard_similarity = jaccard_index(
        tempering_sets[[disease]], non_tempered_sets[[disease]]
      ),
      pearson_probability_correlation = stats::cor(
        tempering_probability[, disease], non_tempered_probability[, disease],
        method = "pearson"
      ),
      spearman_probability_correlation = stats::cor(
        tempering_probability[, disease], non_tempered_probability[, disease],
        method = "spearman"
      ),
      mean_absolute_probability_difference = mean(abs(probability_difference)),
      maximum_absolute_probability_difference = max(abs(probability_difference)),
      stringsAsFactors = FALSE
    )
    edge_rows <- data.frame(
      specification = specification,
      specification_label = specification_label,
      cancer = disease_names[[disease]],
      edge = rownames(tempering_probability),
      adaptive_tempering_probability = tempering_probability[, disease],
      non_tempered_probability = non_tempered_probability[, disease],
      difference_non_tempered_minus_tempering = probability_difference,
      absolute_difference = abs(probability_difference),
      stringsAsFactors = FALSE
    )
    edge_probability_comparison[[length(edge_probability_comparison) + 1L]] <-
      edge_rows
  }

  tempering_diagnostics <- readRDS(paths$tempering_diagnostics)
  non_tempered_diagnostics <- readRDS(paths$non_tempered_diagnostics)
  validate_result_specification(
    tempering_diagnostics, specification, paths$tempering_diagnostics
  )
  validate_result_specification(
    non_tempered_diagnostics, specification, paths$non_tempered_diagnostics
  )
  extract_pooled_log_risk <- function(x, path) {
    summary <- x$diagnostic_results$log_risk$posterior_summary
    if (is.null(summary)) stopf("%s lacks fitted log-risk summaries.", path)
    summary <- summary[summary$chain == "pooled", , drop = FALSE]
    if (!nrow(summary) || anyDuplicated(summary$variable)) {
      stopf("%s has invalid pooled fitted log-risk summaries.", path)
    }
    summary
  }
  tempering_log_risk <- extract_pooled_log_risk(
    tempering_diagnostics, paths$tempering_diagnostics
  )
  non_tempered_log_risk <- extract_pooled_log_risk(
    non_tempered_diagnostics, paths$non_tempered_diagnostics
  )
  log_risk_order <- match(
    tempering_log_risk$variable, non_tempered_log_risk$variable
  )
  if (anyNA(log_risk_order)) {
    stopf("Could not align fitted log-risks for %s.", specification)
  }
  non_tempered_log_risk <- non_tempered_log_risk[log_risk_order, , drop = FALSE]
  parsed_cancer <- sub(
    "^log_risk\\.([^.]+)\\.county_[0-9]+$", "\\1",
    tempering_log_risk$variable
  )
  if (!all(parsed_cancer %in% disease_names)) {
    stopf("Could not parse fitted log-risk cancer labels for %s.", specification)
  }
  log_risk_rows <- data.frame(
    specification = specification,
    specification_label = specification_label,
    cancer = parsed_cancer,
    variable = tempering_log_risk$variable,
    adaptive_tempering_mean_log_risk = tempering_log_risk$mean,
    non_tempered_mean_log_risk = non_tempered_log_risk$mean,
    difference_non_tempered_minus_tempering =
      non_tempered_log_risk$mean - tempering_log_risk$mean,
    stringsAsFactors = FALSE
  )
  fitted_log_risk_values[[length(fitted_log_risk_values) + 1L]] <- log_risk_rows
  fitted_log_risk_comparison[[length(fitted_log_risk_comparison) + 1L]] <-
    do.call(rbind, lapply(disease_names, function(cancer) {
      rows <- log_risk_rows[log_risk_rows$cancer == cancer, , drop = FALSE]
      difference <- rows$difference_non_tempered_minus_tempering
      data.frame(
        specification = specification,
        specification_label = specification_label,
        cancer = cancer,
        pearson_correlation = stats::cor(
          rows$adaptive_tempering_mean_log_risk,
          rows$non_tempered_mean_log_risk, method = "pearson"
        ),
        spearman_correlation = stats::cor(
          rows$adaptive_tempering_mean_log_risk,
          rows$non_tempered_mean_log_risk, method = "spearman"
        ),
        mean_absolute_log_risk_difference = mean(abs(difference)),
        maximum_absolute_log_risk_difference = max(abs(difference)),
        stringsAsFactors = FALSE
      )
    }))
  rm(
    tempering_boundaries, non_tempered_boundaries, tempering_diagnostics,
    non_tempered_diagnostics
  )
  invisible(gc(FALSE))
}

boundary_comparison <- do.call(rbind, boundary_comparison)
edge_probability_comparison <- do.call(rbind, edge_probability_comparison)
fitted_log_risk_comparison <- do.call(rbind, fitted_log_risk_comparison)
fitted_log_risk_values <- do.call(rbind, fitted_log_risk_values)
row.names(boundary_comparison) <- NULL
row.names(edge_probability_comparison) <- NULL
row.names(fitted_log_risk_comparison) <- NULL
row.names(fitted_log_risk_values) <- NULL

# RESULT: Supplementary sampler-sensitivity CSV tables ------------------------
convergence_csv <- write_csv(
  convergence_comparison,
  file.path(
    tables_dir,
    "sampler_convergence_comparison_multichain_three_specifications.csv"
  )
)
boundary_csv <- write_csv(
  boundary_comparison,
  file.path(
    tables_dir,
    "sampler_boundary_robustness_multichain_three_specifications.csv"
  )
)
log_risk_csv <- write_csv(
  fitted_log_risk_comparison,
  file.path(
    tables_dir,
    "sampler_fitted_log_risk_agreement_multichain_three_specifications.csv"
  )
)

# RESULT: Supplementary sampler convergence table -----------------------------
convergence_tex_rows <- vapply(seq_len(nrow(convergence_comparison)), function(i) {
  x <- convergence_comparison[i, ]
  paste0(
    x$specification_label, " & ", x$workflow, " & ", x$quantity, " & ",
    sprintf("%.3f", x$median_rhat), " & ",
    sprintf("%.3f", x$max_rhat), " & ",
    sprintf("%.1f", x$min_bulk_ess), " & ",
    sprintf("%.1f", x$median_bulk_ess), " & ",
    sprintf("%.1f", x$min_tail_ess), " & ",
    sprintf("%.3f", x$max_relative_mcse_mean), " \\\\"
  )
}, character(1L))
convergence_tex <- write_tex(
  c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\scriptsize",
    "\\caption{Multiple-chain convergence diagnostics comparing adaptive tempering and independent non-tempered sampling for all three SEER specifications. Statistics summarize the quantities within each row. Relative MCSE is the MCSE for the posterior mean divided by the pooled posterior standard deviation.}",
    "\\label{tab:seer_sampler_convergence_three_specifications}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{lllrrrrrr}",
    "\\hline",
    "Specification & Workflow & Quantity & Median $\\widehat R$ & Max. $\\widehat R$ & Min. bulk ESS & Median bulk ESS & Min. tail ESS & Max. rel. MCSE \\\\",
    "\\hline",
    convergence_tex_rows,
    "\\hline",
    "\\end{tabular}",
    "}",
    "\\end{table}"
  ),
  file.path(
    tables_dir,
    "sampler_convergence_comparison_multichain_three_specifications.tex"
  )
)

# RESULT: Supplementary sampler boundary-robustness table ----------------------
boundary_tex_rows <- vapply(seq_len(nrow(boundary_comparison)), function(i) {
  x <- boundary_comparison[i, ]
  paste0(
    x$specification_label, " & ", x$cancer, " & ",
    x$adaptive_tempering_count, " & ", x$non_tempered_count, " & ",
    sprintf("%.3f", x$jaccard_similarity), " & ",
    sprintf("%.3f", x$pearson_probability_correlation), " & ",
    sprintf("%.3f", x$mean_absolute_probability_difference), " \\\\"
  )
}, character(1L))
boundary_tex <- write_tex(
  c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\small",
    "\\caption{Sensitivity of pooled posterior-FDR boundary conclusions to adaptive-tempering versus independent non-tempered sampling for all three SEER specifications.}",
    "\\label{tab:seer_sampler_boundary_robustness_three_specifications}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{llrrrrr}",
    "\\hline",
    "Specification & Cancer & Tempering & Non-tempered & Jaccard & Probability correlation & Mean absolute difference \\\\",
    "\\hline",
    boundary_tex_rows,
    "\\hline",
    "\\end{tabular}",
    "}",
    "\\end{table}"
  ),
  file.path(
    tables_dir,
    "sampler_boundary_robustness_multichain_three_specifications.tex"
  )
)

# RESULT: Supplementary fitted log-risk agreement table -----------------------
log_risk_tex_rows <- vapply(
  seq_len(nrow(fitted_log_risk_comparison)),
  function(i) {
    x <- fitted_log_risk_comparison[i, ]
    paste0(
      x$specification_label, " & ", x$cancer, " & ",
      sprintf("%.3f", x$pearson_correlation), " & ",
      sprintf("%.3f", x$spearman_correlation), " & ",
      sprintf("%.3f", x$mean_absolute_log_risk_difference), " & ",
      sprintf("%.3f", x$maximum_absolute_log_risk_difference), " \\\\"
    )
  },
  character(1L)
)
log_risk_tex <- write_tex(
  c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\small",
    "\\caption{Agreement between pooled posterior mean fitted log-risks under adaptive-tempering and independent non-tempered sampling for all three SEER specifications.}",
    "\\label{tab:seer_sampler_log_risk_agreement_three_specifications}",
    "\\begin{tabular}{llrrrr}",
    "\\hline",
    "Specification & Cancer & Pearson & Spearman & Mean absolute difference & Maximum absolute difference \\\\",
    "\\hline",
    log_risk_tex_rows,
    "\\hline",
    "\\end{tabular}",
    "\\end{table}"
  ),
  file.path(
    tables_dir,
    "sampler_fitted_log_risk_agreement_multichain_three_specifications.tex"
  )
)

# RESULT: Supplementary sampler boundary-probability agreement figure ----------
specification_levels <- specifications$label
edge_probability_comparison$specification_label <- factor(
  edge_probability_comparison$specification_label,
  levels = specification_levels
)
edge_probability_comparison$cancer <- factor(
  edge_probability_comparison$cancer, levels = disease_names
)
boundary_plot <- ggplot2::ggplot(
  edge_probability_comparison,
  ggplot2::aes(
    x = adaptive_tempering_probability,
    y = non_tempered_probability
  )
) +
  ggplot2::geom_abline(
    intercept = 0, slope = 1, linetype = "dashed", color = "grey45"
  ) +
  ggplot2::geom_point(alpha = 0.55, size = 0.7, color = "#2166AC") +
  ggplot2::facet_grid(cancer ~ specification_label) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    x = "Adaptive-tempering posterior boundary probability",
    y = "Non-tempered posterior boundary probability"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
boundary_plot_pdf <- file.path(
  figures_dir,
  "sampler_boundary_probability_agreement_multichain_three_specifications.pdf"
)
ggplot2::ggsave(boundary_plot_pdf, boundary_plot, width = 9, height = 10)

# RESULT: Supplementary sampler fitted-log-risk agreement figure ---------------
fitted_log_risk_values$specification_label <- factor(
  fitted_log_risk_values$specification_label, levels = specification_levels
)
fitted_log_risk_values$cancer <- factor(
  fitted_log_risk_values$cancer, levels = disease_names
)
log_risk_plot <- ggplot2::ggplot(
  fitted_log_risk_values,
  ggplot2::aes(
    x = adaptive_tempering_mean_log_risk,
    y = non_tempered_mean_log_risk
  )
) +
  ggplot2::geom_abline(
    intercept = 0, slope = 1, linetype = "dashed", color = "grey45"
  ) +
  ggplot2::geom_point(alpha = 0.70, size = 0.9, color = "#B2182B") +
  ggplot2::facet_grid(cancer ~ specification_label, scales = "free") +
  ggplot2::labs(
    x = "Adaptive-tempering posterior mean fitted log-risk",
    y = "Non-tempered posterior mean fitted log-risk"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
log_risk_plot_pdf <- file.path(
  figures_dir,
  "sampler_fitted_log_risk_agreement_multichain_three_specifications.pdf"
)
ggplot2::ggsave(log_risk_plot_pdf, log_risk_plot, width = 9, height = 10)

artifacts <- c(
  convergence_csv = convergence_csv,
  convergence_tex = convergence_tex,
  boundary_csv = boundary_csv,
  boundary_tex = boundary_tex,
  fitted_log_risk_csv = log_risk_csv,
  fitted_log_risk_tex = log_risk_tex,
  boundary_probability_figure_pdf = normalized_path(boundary_plot_pdf),
  fitted_log_risk_figure_pdf = normalized_path(log_risk_plot_pdf)
)
if (!all(file.exists(unname(artifacts)))) {
  stop("One or more sampler-sensitivity artifacts were not created.",
       call. = FALSE)
}

cat("Sampler sensitivity analysis complete:\n")
cat("  Output:", normalized_path(output_dir), "\n")
cat("  Reported tables and figures:\n")
cat(paste0("    - ", unname(artifacts), collapse = "\n"), "\n")
