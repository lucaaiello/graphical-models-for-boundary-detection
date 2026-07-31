# Complete the observation-level SEER posterior predictive summaries.
#
# This script does not run MCMC and does not recompute the completed global
# posterior predictive checks. The existing global and pointwise CSV files are
# treated as validation references. Replicated county-level counts are
# regenerated only because the stored pointwise file lacks predictive medians;
# every previously stored pointwise quantity is checked before new outputs are
# written.
#
# Published result index (line numbers are also documented in README.md):
#   Main Table 2: lines 413--476
#   Supplementary Table S19: lines 477--550
#   Supplementary Figure S28: lines 551--665
# Run this complete script after posterior_predictive_sensitivity.R: the result
# blocks use the validated pointwise posterior predictive data built above.

rm(list = ls())

project_markers <- c("README.md", "data analysis")
if (!all(file.exists(project_markers))) {
  stop("Run this script from the repository root.", call. = FALSE)
}

required_packages <- c("ggplot2", "matrixStats", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Run install.R before this script.",
    call. = FALSE
  )
}

output_dir <- file.path(
  "data analysis",
  "posterior_predictive_sensitivity_output"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_counties <- 58L
cancer_names <- c("Lung", "Esophageal", "Larynx", "Colorectal")
n_cancers <- length(cancer_names)
unusual_threshold <- 0.02

specifications <- list(
  list(
    file = file.path("data analysis", "results_adj.RData"),
    label = "Adjacency",
    seed = 20260725L
  ),
  list(
    file = file.path("data analysis", "results_meanadj.RData"),
    label = "Mean--Adjacency",
    seed = 20260726L
  ),
  list(
    file = file.path("data analysis", "results_mean.RData"),
    label = "Mean",
    seed = 20260727L
  )
)
specification_levels <- vapply(specifications, `[[`, character(1), "label")

required_inputs <- c(
  vapply(specifications, `[[`, character(1), "file"),
  file.path(output_dir, "ppc_summary.csv"),
  file.path(output_dir, "ppc_unusual_county_cancer.csv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input file(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

standardise_specification <- function(x) {
  replacements <- c(
    "Adaptive adjacency, intercept-only mean" = "Adjacency",
    "Primary adaptive" = "Adjacency",
    "Adaptive adjacency, covariate-adjusted mean" = "Mean--Adjacency",
    "Adjusted adaptive" = "Mean--Adjacency",
    "Fixed adjacency, covariate-adjusted mean" = "Mean",
    "Adjusted fixed" = "Mean"
  )
  x <- as.character(x)
  matched <- match(x, names(replacements))
  x[!is.na(matched)] <- unname(replacements[matched[!is.na(matched)]])
  x
}

write_csv <- function(x, filename) {
  path <- file.path(output_dir, filename)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}

reorder_draws_to_geographic <- function(draw_matrix, final_perm) {
  inverse_permutation <- order(final_perm)
  do.call(
    cbind,
    lapply(seq_len(n_cancers), function(disease) {
      columns <- ((disease - 1L) * n_counties + 1L):(disease * n_counties)
      draw_matrix[, columns, drop = FALSE][
        , inverse_permutation, drop = FALSE
      ]
    })
  )
}

reorder_vector_to_geographic <- function(x, final_perm) {
  inverse_permutation <- order(final_perm)
  unlist(
    lapply(seq_len(n_cancers), function(disease) {
      indices <- ((disease - 1L) * n_counties + 1L):(disease * n_counties)
      x[indices][inverse_permutation]
    }),
    use.names = FALSE
  )
}

pointwise_summary <- function(
  observed,
  expected,
  replicated,
  counties,
  cancer,
  specification
) {
  predictive_mean <- colMeans(replicated)
  predictive_sd <- matrixStats::colSds(replicated)
  predictive_quantiles <- matrixStats::colQuantiles(
    replicated,
    probs = c(0.025, 0.5, 0.975),
    drop = FALSE
  )
  lower_mid_tail <- vapply(
    seq_along(observed),
    function(index) {
      mean(replicated[, index] < observed[[index]]) +
        0.5 * mean(replicated[, index] == observed[[index]])
    },
    numeric(1)
  )
  upper_mid_tail <- vapply(
    seq_along(observed),
    function(index) {
      mean(replicated[, index] > observed[[index]]) +
        0.5 * mean(replicated[, index] == observed[[index]])
    },
    numeric(1)
  )
  two_sided_mid_tail <- pmin(
    1,
    2 * pmin(lower_mid_tail, upper_mid_tail)
  )

  data.frame(
    specification = specification,
    county = tools::toTitleCase(as.character(counties)),
    cancer = cancer,
    observed_count = as.numeric(observed),
    expected_count = as.numeric(expected),
    predictive_mean = predictive_mean,
    predictive_median = predictive_quantiles[, 2],
    predictive_lower_95 = predictive_quantiles[, 1],
    predictive_upper_95 = predictive_quantiles[, 3],
    pointwise_95_coverage = (
      observed >= predictive_quantiles[, 1] &
        observed <= predictive_quantiles[, 3]
    ),
    standardized_predictive_residual = ifelse(
      predictive_sd > 0,
      (observed - predictive_mean) / predictive_sd,
      NA_real_
    ),
    lower_mid_tail_probability = lower_mid_tail,
    upper_mid_tail_probability = upper_mid_tail,
    two_sided_mid_tail_probability = two_sided_mid_tail,
    low_predictive_support = two_sided_mid_tail < unusual_threshold,
    stringsAsFactors = FALSE
  )
}

regenerate_one_specification <- function(definition) {
  message(
    "Regenerating only missing observation-level summaries for ",
    definition$label,
    " (seed ",
    definition$seed,
    ")..."
  )
  workspace <- new.env(parent = baseenv())
  loaded <- load(definition$file, envir = workspace)
  required <- c(
    "mcmc_samples", "X", "Y", "E", "final_perm", "county.ID"
  )
  missing <- setdiff(required, loaded)
  if (length(missing) > 0L) {
    stop(
      definition$file,
      " is missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  mcmc_samples <- get("mcmc_samples", envir = workspace)
  beta <- unname(mcmc_samples[["beta"]])
  phi <- unname(mcmc_samples[["phi"]])
  X <- unname(get("X", envir = workspace))
  Y <- as.numeric(get("Y", envir = workspace))
  E <- as.numeric(get("E", envir = workspace))
  final_perm <- as.integer(get("final_perm", envir = workspace))
  counties <- as.character(get("county.ID", envir = workspace))

  if (
    !is.matrix(beta) || !is.matrix(phi) || !is.matrix(X) ||
      nrow(beta) != 10000L || nrow(phi) != nrow(beta) ||
      ncol(phi) != n_counties * n_cancers ||
      nrow(X) != n_counties * n_cancers || ncol(beta) != ncol(X) ||
      length(Y) != n_counties * n_cancers ||
      length(E) != n_counties * n_cancers ||
      !identical(sort(final_perm), seq_len(n_counties)) ||
      length(counties) != n_counties
  ) {
    stop(definition$file, " has incompatible saved objects.", call. = FALSE)
  }

  set.seed(definition$seed)
  poisson_mean <- sweep(
    exp(tcrossprod(beta, X) + phi),
    2,
    E,
    FUN = "*"
  )
  replicated_model_order <- matrix(
    stats::rpois(length(poisson_mean), lambda = as.vector(poisson_mean)),
    nrow = nrow(beta),
    ncol = n_counties * n_cancers
  )
  rm(poisson_mean, beta, phi, X, mcmc_samples)
  invisible(gc())

  replicated <- reorder_draws_to_geographic(
    replicated_model_order,
    final_perm
  )
  observed <- reorder_vector_to_geographic(Y, final_perm)
  expected <- reorder_vector_to_geographic(E, final_perm)
  rm(replicated_model_order, workspace)
  invisible(gc())

  rows <- lapply(seq_len(n_cancers), function(disease) {
    indices <- ((disease - 1L) * n_counties + 1L):(disease * n_counties)
    pointwise_summary(
      observed = observed[indices],
      expected = expected[indices],
      replicated = replicated[, indices, drop = FALSE],
      counties = counties,
      cancer = cancer_names[[disease]],
      specification = definition$label
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  rm(replicated)
  invisible(gc())
  result
}

# Reuse and validate the completed global PPC summary --------------------------

global_ppc <- utils::read.csv(
  file.path(output_dir, "ppc_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
global_ppc$specification <- standardise_specification(
  global_ppc$specification
)
if (
  nrow(global_ppc) != 48L ||
    sum(global_ppc$observed_within_95) != 46L ||
    !setequal(global_ppc$specification, specification_levels)
) {
  stop("The stored global PPC summary does not verify 46 of 48 checks.")
}

global_exceptions <- global_ppc[!global_ppc$observed_within_95, ]
global_exceptions$direction <- ifelse(
  global_exceptions$observed > global_exceptions$predictive_upper_95,
  "Observed above the 97.5% predictive quantile",
  "Observed below the 2.5% predictive quantile"
)
global_exceptions <- global_exceptions[
  ,
  c(
    "specification", "cancer", "diagnostic_label", "observed",
    "predictive_mean", "predictive_median", "predictive_lower_95",
    "predictive_upper_95", "two_sided_ppp", "direction"
  )
]
if (nrow(global_exceptions) != 2L) {
  stop("Expected exactly two stored global PPC exceptions.")
}
write_csv(global_exceptions, "global_ppc_exceptions.csv")
writeLines(
  c(
    "Global posterior predictive exceptions",
    "======================================",
    "The stored global PPC summary contains 48 checks; 46 are covered.",
    "",
    capture.output(print(global_exceptions, row.names = FALSE, digits = 6))
  ),
  file.path(output_dir, "global_ppc_exceptions.txt"),
  useBytes = TRUE
)

# Regenerate the missing medians, then verify all existing pointwise values ----

observation_ppc <- do.call(
  rbind,
  lapply(specifications, regenerate_one_specification)
)
rownames(observation_ppc) <- NULL

if (
  nrow(observation_ppc) != 696L ||
    any(table(observation_ppc$specification) != 232L) ||
    any(table(observation_ppc$specification, observation_ppc$cancer) != 58L)
) {
  stop("The regenerated observation-level PPC does not contain 696 rows.")
}

stored_pointwise <- utils::read.csv(
  file.path(output_dir, "ppc_unusual_county_cancer.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stored_pointwise$specification <- standardise_specification(
  stored_pointwise$specification
)
key <- function(x) paste(x$specification, x$cancer, x$county, sep = "||")
stored_order <- match(key(observation_ppc), key(stored_pointwise))
if (anyNA(stored_order) || anyDuplicated(key(stored_pointwise))) {
  stop("Stored and regenerated pointwise PPC keys do not agree.")
}
stored_pointwise <- stored_pointwise[stored_order, ]

numeric_checks <- c(
  "observed_count" = "observed_count",
  "expected_count" = "expected_count",
  "predictive_mean" = "predictive_mean",
  "predictive_lower_95" = "predictive_lower_95",
  "predictive_upper_95" = "predictive_upper_95",
  "predictive_residual" = "standardized_predictive_residual",
  "lower_tail_probability" = "lower_mid_tail_probability",
  "upper_tail_probability" = "upper_mid_tail_probability",
  "two_sided_tail_probability" = "two_sided_mid_tail_probability"
)
maximum_difference <- vapply(
  names(numeric_checks),
  function(stored_name) {
    regenerated_name <- unname(numeric_checks[[stored_name]])
    max(
      abs(
        stored_pointwise[[stored_name]] -
          observation_ppc[[regenerated_name]]
      ),
      na.rm = TRUE
    )
  },
  numeric(1)
)
if (any(maximum_difference > 1e-10)) {
  stop(
    "Regenerated pointwise values disagree with stored values: ",
    paste(
      names(maximum_difference),
      format(maximum_difference, scientific = TRUE),
      collapse = "; "
    ),
    call. = FALSE
  )
}
if (!identical(
  as.logical(stored_pointwise$unusual_observation),
  observation_ppc$low_predictive_support
)) {
  stop("Stored and regenerated low-support indicators disagree.")
}

observation_ppc$specification <- factor(
  observation_ppc$specification,
  levels = specification_levels
)
observation_ppc$cancer <- factor(
  observation_ppc$cancer,
  levels = cancer_names
)
observation_ppc <- observation_ppc[
  order(observation_ppc$specification, observation_ppc$cancer),
]
observation_ppc$specification <- as.character(observation_ppc$specification)
observation_ppc$cancer <- as.character(observation_ppc$cancer)
write_csv(observation_ppc, "observation_level_ppc.csv")

# RESULT: Main Table 2 (lines 413--476) ----------------------------------------
# Writes observation_ppc_summary.csv and observation_ppc_summary.tex.

observation_summary <- do.call(
  rbind,
  lapply(specification_levels, function(specification) {
    x <- observation_ppc[observation_ppc$specification == specification, ]
    data.frame(
      Specification = specification,
      Observed = mean(x$observed_count),
      Predicted = mean(x$predictive_mean),
      Lower = mean(x$predictive_lower_95),
      Upper = mean(x$predictive_upper_95),
      Coverage = mean(x$pointwise_95_coverage),
      `Low support` = sum(x$low_predictive_support),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
)
if (length(unique(observation_summary$Observed)) != 1L) {
  stop("The observed mean differs across specifications.")
}
write_csv(observation_summary, "observation_ppc_summary.csv")

summary_tex_rows <- sprintf(
  "%s & %.1f & %.1f & %.1f & %.1f & %.3f & %d \\\\",
  observation_summary$Specification,
  observation_summary$Observed,
  observation_summary$Predicted,
  observation_summary$Lower,
  observation_summary$Upper,
  observation_summary$Coverage,
  observation_summary[["Low support"]]
)
summary_tex <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\small",
  paste0(
    "\\caption{Observation-level posterior predictive summaries across ",
    "the three specifications.}"
  ),
  "\\label{tab:observation_ppc}",
  "\\begin{tabular}{lrrrrrr}",
  "\\hline",
  "Specification & Observed & Predicted & Lower & Upper & Coverage & Low support \\\\",
  "\\hline",
  summary_tex_rows,
  "\\hline",
  "\\end{tabular}",
  "\\end{table}",
  "",
  paste0(
    "Lower and Upper are averages of the pointwise 2.5\\% and 97.5\\% ",
    "posterior predictive quantiles across the 232 county--cancer observations."
  )
)
writeLines(
  summary_tex,
  file.path(output_dir, "observation_ppc_summary.tex"),
  useBytes = TRUE
)

# RESULT: Supplementary Table S19 (lines 477--550) -----------------------------
# Writes unusual_observation_ppc.csv and unusual_observation_ppc.tex.

unusual <- observation_ppc[
  observation_ppc$two_sided_mid_tail_probability < unusual_threshold,
]
if (
  nrow(unusual) != 2L ||
    !all(unusual$county == "Mono") ||
    !all(unusual$cancer == "Lung") ||
    !setequal(unusual$specification, c("Mean--Adjacency", "Mean"))
) {
  stop(
    "The regenerated pointwise PPC does not reproduce the expected two ",
    "Mono County lung checks.",
    call. = FALSE
  )
}
unusual_output <- data.frame(
  Specification = unusual$specification,
  County = unusual$county,
  Cancer = unusual$cancer,
  Observed = unusual$observed_count,
  `Predictive mean` = unusual$predictive_mean,
  `Predictive median` = unusual$predictive_median,
  `Predictive lower 95` = unusual$predictive_lower_95,
  `Predictive upper 95` = unusual$predictive_upper_95,
  `Standardized residual` = unusual$standardized_predictive_residual,
  `Mid-tail probability` = unusual$two_sided_mid_tail_probability,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_csv(unusual_output, "unusual_observation_ppc.csv")

unusual_tex_rows <- sprintf(
  "%s & %s & %s & %d & %.1f & %.1f & [%g, %g] & %.2f & %.4f \\\\",
  unusual_output$Specification,
  unusual_output$County,
  unusual_output$Cancer,
  unusual_output$Observed,
  unusual_output[["Predictive mean"]],
  unusual_output[["Predictive median"]],
  unusual_output[["Predictive lower 95"]],
  unusual_output[["Predictive upper 95"]],
  unusual_output[["Standardized residual"]],
  unusual_output[["Mid-tail probability"]]
)
unusual_tex <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\small",
  paste0(
    "\\caption{County--cancer observations with two-sided predictive ",
    "mid-tail probability below $0.02$.}"
  ),
  "\\label{tab:unusual_observation_ppc}",
  "\\begin{tabular}{lllrrrcrr}",
  "\\hline",
  paste0(
    "Specification & County & Cancer & Observed & Pred. mean & Pred. median & ",
    "95\\% predictive interval & Std. residual & Mid-tail $p$ \\\\"
  ),
  "\\hline",
  unusual_tex_rows,
  "\\hline",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(
  unusual_tex,
  file.path(output_dir, "unusual_observation_ppc.tex"),
  useBytes = TRUE
)

# RESULT: Supplementary Figure S28 (lines 551--665) ----------------------------
# Builds and saves observed_vs_fitted_ppc.pdf and its PNG counterpart.

plot_data <- observation_ppc
plot_data$specification <- factor(
  plot_data$specification,
  levels = specification_levels
)
plot_data$cancer <- factor(plot_data$cancer, levels = cancer_names)
plot_data$panel <- factor(
  paste(plot_data$cancer, plot_data$specification, sep = "\n"),
  levels = unlist(
    lapply(
      cancer_names,
      function(cancer) paste(cancer, specification_levels, sep = "\n")
    )
  )
)
plot_data$flag_label <- ifelse(
  plot_data$low_predictive_support & plot_data$county == "Mono",
  "Mono",
  ""
)

observed_vs_fitted_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = observed_count, y = predictive_mean)
) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2,
    color = "grey45",
    linewidth = 0.45
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      xend = observed_count,
      y = predictive_lower_95,
      yend = predictive_upper_95,
      color = low_predictive_support
    ),
    linewidth = 0.35,
    alpha = 0.7
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      color = low_predictive_support,
      shape = low_predictive_support
    ),
    size = 1.35,
    alpha = 0.9
  ) +
  ggplot2::geom_label(
    data = plot_data[plot_data$flag_label != "", ],
    ggplot2::aes(label = flag_label),
    hjust = -0.15,
    vjust = 0.5,
    size = 2.6,
    color = "#B2182B",
    fill = "white",
    linewidth = 0.15,
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(
    ~ panel,
    ncol = 3,
    scales = "free"
  ) +
  ggplot2::scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 1),
    breaks = c(0, 10, 100, 1000, 10000),
    labels = scales::label_comma(),
    expand = ggplot2::expansion(mult = c(0.02, 0.12))
  ) +
  ggplot2::scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 1),
    breaks = c(0, 10, 100, 1000, 10000),
    labels = scales::label_comma(),
    expand = ggplot2::expansion(mult = c(0.02, 0.12))
  ) +
  ggplot2::scale_color_manual(
    values = c(`FALSE` = "#2166AC", `TRUE` = "#B2182B"),
    labels = c(`FALSE` = "Other", `TRUE` = "Mid-tail p < 0.02"),
    name = NULL
  ) +
  ggplot2::scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 17),
    guide = "none"
  ) +
  ggplot2::labs(
    x = "Observed county-level count (pseudo-log scale)",
    y = "Posterior predictive mean count (pseudo-log scale)"
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(output_dir, "observed_vs_fitted_ppc.pdf"),
  observed_vs_fitted_plot,
  width = 10,
  height = 13
)
ggplot2::ggsave(
  file.path(output_dir, "observed_vs_fitted_ppc.png"),
  observed_vs_fitted_plot,
  width = 10,
  height = 13,
  dpi = 300
)

# Machine-readable completion report ------------------------------------------

completion_lines <- c(
  "Observation-level PPC completion report",
  "=======================================",
  "No MCMC sampler was called.",
  "The completed global PPC was read from ppc_summary.csv, not rerun.",
  paste(
    "Observation-level replicated counts were regenerated only to recover",
    "the missing predictive medians."
  ),
  paste(
    "Maximum absolute discrepancy versus the stored pointwise output:",
    format(max(maximum_difference), scientific = TRUE)
  ),
  "Rows: 696 total; 232 per specification; 58 per cancer and specification.",
  "Global coverage: 46 of 48; exceptions: 2.",
  paste(
    "Observation-level low-support checks:",
    nrow(unusual),
    "(the same Mono County lung observation under two specifications)."
  )
)
writeLines(
  completion_lines,
  file.path(output_dir, "observation_ppc_completion_report.txt"),
  useBytes = TRUE
)
message(paste(completion_lines, collapse = "\n"))
