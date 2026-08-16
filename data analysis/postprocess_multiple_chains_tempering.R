rm(list = ls())

# Shared post-processing for SEER multiple-chain analyses
# =============================================================================
# This script does not run MCMC. It reads six completed posterior chains from
# either the baseline or adaptive-tempering workflow and
# produces:
#
#   * rank-normalized split-Rhat, bulk ESS, tail ESS, and MCSE diagnostics;
#   * posterior-summary and between-chain agreement tables;
#   * proposal acceptance and, for tempering, temperature-traffic diagnostics;
#   * chain-specific and pooled posterior-FDR boundary selections;
#   * edge-probability and selected-boundary stability diagnostics;
#   * pooled and chain-specific FDR curves and boundary maps;
#   * paper-style random-effect, shared-boundary, mutual-boundary, and
#     non-adjacency maps;
#   * pooled posterior caterpillar plots for each tempering specification
#     (Supplementary Figures S3--S6 for Adjacency).
#
# Run from the repository root with:
#   source("data analysis/postprocess_multiple_chains_tempering.R")
#
# Select the workflow with MCMC_POSTPROCESS_WORKFLOW = "tempering" or
# "non_tempered". Both workflows infer their run directory from specification and
# mode variables. Input/output locations can still be overridden with the
# corresponding MCMC_TEMPERING_POSTPROCESS_* or
# MCMC_MULTICHAIN_POSTPROCESS_* variables.
# For tempering, the default input follows the same hierarchy as the sampler:
#   multiple_chains_tempering_output/<specification>/<mode>/
# For non-tempered chains, the default input is:
#   multiple_chains_output/<specification>/<mode>/
# Derived outputs are written to <input>/postprocessing/ for both workflows.

options(stringsAsFactors = FALSE)

postprocess_workflow <- tolower(Sys.getenv(
  "MCMC_POSTPROCESS_WORKFLOW", unset = "tempering"
))
if (!postprocess_workflow %in% c("tempering", "non_tempered")) {
  stop(
    "MCMC_POSTPROCESS_WORKFLOW must be 'tempering' or 'non_tempered'.",
    call. = FALSE
  )
}
is_tempering_workflow <- identical(postprocess_workflow, "tempering")
workflow_label <- if (is_tempering_workflow) {
  "adaptive-tempering"
} else {
  "baseline multiple-chain"
}
output_prefix <- if (is_tempering_workflow) {
  "tempering"
} else {
  "multiple_chains"
}

tempering_specification <- tolower(Sys.getenv(
  "MCMC_TEMPERING_POSTPROCESS_SPECIFICATION",
  unset = Sys.getenv(
    "MCMC_TEMPERING_SPECIFICATION", unset = "adj"
  )
))
if (!tempering_specification %in% c("adj", "meanadj", "mean")) {
  stop(
    "The tempering specification must be 'adj', 'meanadj', or 'mean'.",
    call. = FALSE
  )
}
tempering_mode <- Sys.getenv(
  "MCMC_TEMPERING_POSTPROCESS_MODE",
  unset = Sys.getenv(
    "MCMC_TEMPERING_MODE", unset = "hotter_final"
  )
)
tempering_output_root <- Sys.getenv(
  "MCMC_TEMPERING_OUTPUT_ROOT",
  unset = file.path(
    "data analysis", "multiple_chains_tempering_output"
  )
)

non_tempered_specification <- tolower(Sys.getenv(
  "MCMC_MULTICHAIN_POSTPROCESS_SPECIFICATION",
  unset = Sys.getenv("MCMC_MULTICHAIN_SPECIFICATION", unset = "adj")
))
if (!non_tempered_specification %in% c("adj", "meanadj", "mean")) {
  stop(
    "The non-tempered specification must be 'adj', 'meanadj', or 'mean'.",
    call. = FALSE
  )
}
non_tempered_mode <- Sys.getenv(
  "MCMC_MULTICHAIN_POSTPROCESS_MODE",
  unset = Sys.getenv("MCMC_MULTICHAIN_MODE", unset = "production")
)
if (!non_tempered_mode %in% c("smoke", "quick", "pilot", "production")) {
  stop(
    "The non-tempered mode must be 'smoke', 'quick', 'pilot', or 'production'.",
    call. = FALSE
  )
}
non_tempered_output_root <- Sys.getenv(
  "MCMC_MULTICHAIN_OUTPUT_ROOT",
  unset = file.path("data analysis", "multiple_chains_output")
)
requested_specification <- if (is_tempering_workflow) {
  tempering_specification
} else {
  non_tempered_specification
}
requested_mode <- if (is_tempering_workflow) {
  tempering_mode
} else {
  non_tempered_mode
}

required_packages <- c(
  "posterior", "matrixStats", "ggplot2", "sf", "RColorBrewer", "ggpubr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "),
    ". Run install.R first.", call. = FALSE
  )
}

input_environment_variable <- if (is_tempering_workflow) {
  "MCMC_TEMPERING_POSTPROCESS_INPUT"
} else {
  "MCMC_MULTICHAIN_POSTPROCESS_INPUT"
}
output_environment_variable <- if (is_tempering_workflow) {
  "MCMC_TEMPERING_POSTPROCESS_OUTPUT"
} else {
  "MCMC_MULTICHAIN_POSTPROCESS_OUTPUT"
}
default_input_dir <- if (is_tempering_workflow) {
  file.path(
    tempering_output_root, tempering_specification, tempering_mode
  )
} else {
  file.path(
    non_tempered_output_root, non_tempered_specification, non_tempered_mode
  )
}
input_dir <- Sys.getenv(
  input_environment_variable, unset = default_input_dir
)
output_dir <- Sys.getenv(
  output_environment_variable,
  unset = file.path(input_dir, "postprocessing")
)
diagnostics_dir <- file.path(output_dir, "diagnostics")
figures_dir <- file.path(output_dir, "figures")
supplement_dir <- file.path(output_dir, "supplement")
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

n_chains <- 6L
chain_names <- paste0("chain_", seq_len(n_chains))
chain_files <- if (is_tempering_workflow) {
  file.path(
    input_dir, "ladders",
    sprintf("ladder_%d_cold_chain.rds", seq_len(n_chains))
  )
} else {
  file.path(
    input_dir, "chains",
    sprintf("chain_%d_results.rds", seq_len(n_chains))
  )
}
ladder_diagnostic_files <- if (is_tempering_workflow) {
  file.path(
    input_dir, "ladders",
    sprintf("ladder_%d_diagnostics.rds", seq_len(n_chains))
  )
} else {
  character()
}
required_files <- c(chain_files, ladder_diagnostic_files)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Missing completed ", postprocess_workflow, " output file(s):\n",
    paste(" -", missing_files, collapse = "\n"),
    call. = FALSE
  )
}

cat("", workflow_label, " post-processing\n", sep = "")
cat("Started:", format(Sys.time()), "\n")
cat("Input:", normalizePath(input_dir, winslash = "/"), "\n")
cat("Output:", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n")
cat("SEER specification:", requested_specification, "\n")
cat("Run mode:", requested_mode, "\n")
cat("Chains:", n_chains, "\n")
cat("Diagnostic backend: posterior",
    as.character(utils::packageVersion("posterior")), "\n\n")

# Reconstruct the audited SEER data ordering and California geometry without
# loading a manuscript-result .RData file. Evaluation stops before sampling.
model_driver <- file.path("data analysis", "main_multiple_chains.R")
if (!file.exists(model_driver)) {
  stop("Missing model driver: ", model_driver, call. = FALSE)
}
model_source <- readLines(model_driver, warn = FALSE)
setup_marker <- grep(
  "^# Baseline multiple-chain MCMC workflow",
  model_source
)
if (length(setup_marker) != 1L || setup_marker <= 2L) {
  stop("Could not locate the model driver's model-setup boundary.")
}
setup_environment <- new.env(parent = .GlobalEnv)
suppressPackageStartupMessages(eval(
  parse(text = model_source[2L:(setup_marker - 1L)]),
  envir = setup_environment
))
rm(model_source, setup_marker)

{
  assign("cvrts", requested_specification, envir = setup_environment)
  X1 <- get("X1", setup_environment)
  X2 <- get("X2", setup_environment)
  X3 <- get("X3", setup_environment)
  X4 <- get("X4", setup_environment)
  X <- if (requested_specification %in% c("mean", "meanadj")) {
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
  assign("X", X, envir = setup_environment)
  rm(X1, X2, X3, X4, X)
}

final_perm <- get("final_perm", setup_environment)
neighbor_list0 <- get("neighbor_list0", setup_environment)
neighbor_name <- get("neighbor_name", setup_environment)
Minc <- get("Minc", setup_environment)
X <- get("X", setup_environment)
ca.poly <- get("ca.poly", setup_environment)
cvrts <- get("cvrts", setup_environment)
rm(setup_environment)
invisible(gc(FALSE))

disease_names <- c("Lung", "Esophageal", "Larynx", "Colorectal")
n_diseases <- length(disease_names)
n_counties <- nrow(Minc)
n_edges <- nrow(neighbor_list0)
fdr_alpha <- 0.05

read_chain <- function(chain) {
  result <- readRDS(chain_files[[chain]])
  required_fields <- c(
    "beta", "phi", "theta", "rho", "eta", "V", "r", "A", "tau",
    "W_cardinality", "W_mean", "acceptance", "state_consistency", "settings"
  )
  missing_fields <- setdiff(required_fields, names(result))
  if (length(missing_fields)) {
    stop(
      "Cold chain ", chain, " lacks: ",
      paste(missing_fields, collapse = ", "), call. = FALSE
    )
  }
  result
}

first_chain <- read_chain(1L)
saved_specification <- first_chain$settings$model_specification
if (is.null(saved_specification)) {
  saved_specification <- first_chain$settings$cvrts
}
if (is.null(saved_specification)) saved_specification <- "adj"
if (!identical(saved_specification, requested_specification)) {
  stop(
    "Requested specification '", requested_specification,
    "' does not match the saved chains ('", saved_specification, "').",
    call. = FALSE
  )
}
retained_draws <- nrow(first_chain$beta)
parameters_per_disease <- ncol(first_chain$beta) / n_diseases
if (!isTRUE(all.equal(parameters_per_disease, as.integer(parameters_per_disease)))) {
  stop("Regression coefficients cannot be divided evenly among diseases.")
}
parameters_per_disease <- as.integer(parameters_per_disease)
intercept_columns <- seq(1L, ncol(first_chain$beta), by = parameters_per_disease)
rm(first_chain)

chain_metadata <- do.call(rbind, lapply(seq_len(n_chains), function(chain) {
  fit <- read_chain(chain)
  warmup_iterations <- if (is_tempering_workflow) {
    fit$settings$warmup_iterations
  } else {
    fit$settings$burn
  }
  sampling_thin <- if (is_tempering_workflow) {
    fit$settings$sampling_thin
  } else {
    fit$settings$thin
  }
  data.frame(
    chain = chain_names[[chain]],
    model_specification = if (!is.null(
      fit$settings$model_specification
    )) {
      fit$settings$model_specification
    } else if (!is.null(fit$settings$cvrts)) {
      fit$settings$cvrts
    } else {
      "adj"
    },
    seed = fit$seed,
    initialization_regime = fit$settings$initialization_regime,
    retained_draws = nrow(fit$beta),
    warmup_iterations = warmup_iterations,
    sampling_thin = sampling_thin,
    n_temperatures = if (is_tempering_workflow) {
      fit$settings$n_temperatures
    } else {
      1L
    },
    learned_beta_hot = if (is_tempering_workflow) {
      fit$settings$learned_beta_hot
    } else {
      NA_real_
    },
    runtime_hours = fit$runtime$elapsed_seconds / 3600,
    stringsAsFactors = FALSE
  )
}))
if (length(unique(chain_metadata$retained_draws)) != 1L) {
  stop("Cold chains have different numbers of retained draws.")
}

# Parameter transformations ----------------------------------------------------

beta_names <- if (parameters_per_disease == 1L) {
  paste0("centered_intercept.", disease_names)
} else {
  as.vector(t(outer(
    disease_names,
    c("intercept", "smoking", "old", "poverty")[
      seq_len(parameters_per_disease)
    ],
    paste,
    sep = "."
  )))
}

spatial_names <- paste0(
  rep(disease_names, each = n_counties),
  ".county_", rep(seq_len(n_counties), times = n_diseases)
)

center_chain <- function(fit) {
  centered_phi <- fit$phi
  phi_means <- matrix(0, nrow(fit$phi), n_diseases)
  for (disease in seq_len(n_diseases)) {
    indices <- ((disease - 1L) * n_counties + 1L):(disease * n_counties)
    phi_means[, disease] <- rowMeans(fit$phi[, indices, drop = FALSE])
    centered_phi[, indices] <-
      fit$phi[, indices, drop = FALSE] - phi_means[, disease]
  }
  centered_beta <- fit$beta
  centered_beta[, intercept_columns] <-
    fit$beta[, intercept_columns, drop = FALSE] + phi_means
  list(beta = centered_beta, phi = centered_phi)
}

build_core_parameters <- function(fit) {
  result <- cbind(tau = fit$tau, fit$rho)
  colnames(result) <- c("tau", paste0("rho.", disease_names))
  if (cvrts != "mean") {
    result <- cbind(result, fit$eta)
    colnames(result)[
      (ncol(result) - ncol(fit$eta) + 1L):ncol(result)
    ] <- paste0(
      "eta.",
      rep(c("smoking", "old", "poverty"), each = n_diseases),
      ".", rep(disease_names, times = 3L)
    )
  }
  result
}

build_centered_beta <- function(fit) {
  result <- center_chain(fit)$beta
  colnames(result) <- beta_names
  result
}

build_centered_phi <- function(fit) {
  result <- center_chain(fit)$phi
  colnames(result) <- paste0("centered_phi.", spatial_names)
  result
}

build_log_risk <- function(fit) {
  result <- fit$beta %*% t(X) + fit$phi
  colnames(result) <- paste0("log_risk.", spatial_names)
  result
}

build_latent_r <- function(fit) {
  result <- fit$r
  colnames(result) <- paste0("latent_r.", spatial_names)
  result
}

build_raw_A <- function(fit) {
  result <- fit$A
  colnames(result) <- paste0("A_coordinate.", seq_len(ncol(result)))
  result
}

A_vector_to_matrix <- function(A_vector) {
  A_matrix <- matrix(0, n_diseases, n_diseases)
  index <- 1L
  for (row in seq_len(n_diseases)) {
    for (column in seq_len(row)) {
      A_matrix[row, column] <- if (row == column) {
        exp(A_vector[[index]])
      } else {
        A_vector[[index]]
      }
      index <- index + 1L
    }
  }
  A_matrix
}

covariance_indices <- which(
  lower.tri(matrix(0, n_diseases, n_diseases), diag = TRUE),
  arr.ind = TRUE
)
correlation_indices <- which(
  lower.tri(matrix(0, n_diseases, n_diseases), diag = FALSE),
  arr.ind = TRUE
)

build_dependence_parameters <- function(fit) {
  result <- matrix(
    NA_real_, nrow(fit$A),
    nrow(covariance_indices) + nrow(correlation_indices)
  )
  for (draw in seq_len(nrow(fit$A))) {
    A_matrix <- A_vector_to_matrix(fit$A[draw, ])
    Sigma <- A_matrix %*% t(A_matrix)
    correlation <- stats::cov2cor(Sigma)
    result[draw, ] <- c(
      Sigma[covariance_indices], correlation[correlation_indices]
    )
  }
  colnames(result) <- c(
    paste0(
      "Sigma.", covariance_indices[, "row"], ".",
      covariance_indices[, "col"]
    ),
    paste0(
      "correlation.", correlation_indices[, "row"], ".",
      correlation_indices[, "col"]
    )
  )
  result
}

build_adjacency_cardinality <- function(fit) {
  result <- fit$W_cardinality
  colnames(result) <- paste0("adjacency_edges.", disease_names)
  result
}

build_boundary_count <- function(fit) {
  result <- matrix(NA_integer_, nrow(fit$phi), n_diseases)
  for (disease in seq_len(n_diseases)) {
    indices <- ((disease - 1L) * n_counties + 1L):(disease * n_counties)
    phi <- fit$phi[, indices, drop = FALSE]
    phi <- phi[, order(final_perm), drop = FALSE]
    result[, disease] <- rowSums(
      phi[, neighbor_list0[, 1L], drop = FALSE] !=
        phi[, neighbor_list0[, 2L], drop = FALSE]
    )
  }
  colnames(result) <- paste0("boundary_count.", disease_names)
  result
}

# Multiple-chain convergence diagnostics --------------------------------------

summarise_draw_matrix <- function(draw_matrix, chain_label) {
  quantiles <- matrixStats::colQuantiles(
    draw_matrix, probs = c(0.025, 0.50, 0.975), drop = FALSE
  )
  data.frame(
    chain = chain_label,
    variable = colnames(draw_matrix),
    mean = colMeans(draw_matrix),
    sd = matrixStats::colSds(draw_matrix),
    median = quantiles[, 2L],
    q025 = quantiles[, 1L],
    q975 = quantiles[, 3L],
    row.names = NULL
  )
}

summarise_pooled_array <- function(draw_array) {
  variable_names <- dimnames(draw_array)[[3L]]
  statistics <- vapply(seq_along(variable_names), function(variable) {
    values <- as.vector(draw_array[, , variable])
    quantiles <- stats::quantile(
      values, c(0.025, 0.50, 0.975), names = FALSE
    )
    c(
      mean = mean(values), sd = stats::sd(values),
      median = quantiles[[2L]], q025 = quantiles[[1L]],
      q975 = quantiles[[3L]]
    )
  }, numeric(5L))
  data.frame(
    chain = "pooled", variable = variable_names,
    mean = statistics[1L, ], sd = statistics[2L, ],
    median = statistics[3L, ], q025 = statistics[4L, ],
    q975 = statistics[5L, ], row.names = NULL
  )
}

summarise_posterior_agreement <- function(chain_summaries, pooled_summary) {
  do.call(rbind, lapply(seq_len(nrow(pooled_summary)), function(index) {
    variable <- pooled_summary$variable[[index]]
    summaries <- chain_summaries[
      chain_summaries$variable == variable, , drop = FALSE
    ]
    pooled_sd <- pooled_summary$sd[[index]]
    mean_range <- diff(range(summaries$mean))
    overlap_lower <- max(summaries$q025)
    overlap_upper <- min(summaries$q975)
    data.frame(
      variable = variable,
      min_chain_mean = min(summaries$mean),
      max_chain_mean = max(summaries$mean),
      between_chain_mean_range = mean_range,
      standardized_mean_range = if (is.finite(pooled_sd) && pooled_sd > 0) {
        mean_range / pooled_sd
      } else {
        NA_real_
      },
      between_chain_median_range = diff(range(summaries$median)),
      credible_intervals_all_overlap = overlap_lower <= overlap_upper,
      common_interval_lower = overlap_lower,
      common_interval_upper = overlap_upper,
      row.names = NULL
    )
  }))
}

diagnose_draw_group <- function(group_key, group_label, builder) {
  cat("Diagnostics:", group_label, "... ")
  first_fit <- read_chain(1L)
  first_matrix <- as.matrix(builder(first_fit))
  rm(first_fit)
  variable_names <- colnames(first_matrix)
  draw_array <- array(
    NA_real_,
    dim = c(nrow(first_matrix), n_chains, ncol(first_matrix)),
    dimnames = list(
      iteration = NULL, chain = chain_names, variable = variable_names
    )
  )
  draw_array[, 1L, ] <- first_matrix
  rm(first_matrix)
  if (n_chains > 1L) {
    for (chain in 2:n_chains) {
      fit <- read_chain(chain)
      draw_array[, chain, ] <- builder(fit)
      rm(fit)
      invisible(gc(FALSE))
    }
  }

  chain_summaries <- do.call(rbind, lapply(seq_len(n_chains), function(chain) {
    chain_matrix <- matrix(
      draw_array[, chain, ], nrow = dim(draw_array)[1L],
      ncol = dim(draw_array)[3L],
      dimnames = list(NULL, variable_names)
    )
    summarise_draw_matrix(chain_matrix, chain_names[[chain]])
  }))
  pooled_summary <- summarise_pooled_array(draw_array)
  posterior_summary <- rbind(chain_summaries, pooled_summary)
  posterior_summary$group <- group_label
  posterior_summary <- posterior_summary[
    , c("group", "chain", "variable", "mean", "sd", "median", "q025", "q975")
  ]
  posterior_agreement <- summarise_posterior_agreement(
    chain_summaries, pooled_summary
  )
  posterior_agreement$group <- group_label
  posterior_agreement <- posterior_agreement[
    , c("group", setdiff(names(posterior_agreement), "group"))
  ]

  diagnostic_summary <- as.data.frame(posterior::summarise_draws(
    posterior::as_draws_array(draw_array),
    "rhat", "ess_bulk", "ess_tail", "mcse_mean"
  ))
  diagnostic_summary$relative_mcse_mean <-
    diagnostic_summary$mcse_mean /
    pooled_summary$sd[
      match(diagnostic_summary$variable, pooled_summary$variable)
    ]
  diagnostic_summary$relative_mcse_mean[
    !is.finite(diagnostic_summary$relative_mcse_mean)
  ] <- NA_real_

  finite_rhat <- is.finite(diagnostic_summary$rhat)
  finite_bulk <- is.finite(diagnostic_summary$ess_bulk)
  finite_tail <- is.finite(diagnostic_summary$ess_tail)
  finite_mcse <- is.finite(diagnostic_summary$mcse_mean)
  finite_relative <- is.finite(diagnostic_summary$relative_mcse_mean)
  overview <- data.frame(
    group = group_label,
    variables = sum(finite_rhat),
    median_rhat = median(diagnostic_summary$rhat[finite_rhat]),
    max_rhat = max(diagnostic_summary$rhat[finite_rhat]),
    rhat_above_1.01 = sum(diagnostic_summary$rhat[finite_rhat] > 1.01),
    rhat_above_1.05 = sum(diagnostic_summary$rhat[finite_rhat] > 1.05),
    rhat_above_1.10 = sum(diagnostic_summary$rhat[finite_rhat] > 1.10),
    min_bulk_ess = min(diagnostic_summary$ess_bulk[finite_bulk]),
    median_bulk_ess = median(diagnostic_summary$ess_bulk[finite_bulk]),
    min_tail_ess = min(diagnostic_summary$ess_tail[finite_tail]),
    median_tail_ess = median(diagnostic_summary$ess_tail[finite_tail]),
    max_mcse_mean = max(diagnostic_summary$mcse_mean[finite_mcse]),
    max_relative_mcse_mean = max(
      diagnostic_summary$relative_mcse_mean[finite_relative]
    ),
    row.names = NULL
  )
  worst <- head(
    diagnostic_summary[
      order(diagnostic_summary$rhat, decreasing = TRUE, na.last = TRUE),
      c(
        "variable", "rhat", "ess_bulk", "ess_tail",
        "mcse_mean", "relative_mcse_mean"
      )
    ],
    20L
  )
  agreement_worst <- head(
    posterior_agreement[
      order(
        posterior_agreement$standardized_mean_range,
        decreasing = TRUE, na.last = TRUE
      ),
    ],
    20L
  )

  utils::write.csv(
    diagnostic_summary,
    file.path(
      diagnostics_dir,
      paste0("rhat_ess_mcse_", group_key, "_", cvrts, ".csv")
    ),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    worst,
    file.path(
      diagnostics_dir,
      paste0("worst_rhat_", group_key, "_", cvrts, ".csv")
    ),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    posterior_summary,
    file.path(
      diagnostics_dir,
      paste0("posterior_summary_", group_key, "_", cvrts, ".csv")
    ),
    row.names = FALSE, na = ""
  )
  utils::write.csv(
    posterior_agreement,
    file.path(
      diagnostics_dir,
      paste0("posterior_agreement_", group_key, "_", cvrts, ".csv")
    ),
    row.names = FALSE, na = ""
  )

  rm(draw_array)
  invisible(gc(FALSE))
  cat("done\n")
  list(
    overview = overview,
    summary = diagnostic_summary,
    worst = worst,
    posterior_summary = posterior_summary,
    posterior_agreement = posterior_agreement,
    posterior_agreement_worst = agreement_worst
  )
}

diagnostic_specifications <- list(
  core = list("Core parameters", build_core_parameters),
  centered_beta = list("Centered regression coefficients", build_centered_beta),
  raw_A = list("Raw covariance-factor coordinates", build_raw_A),
  dependence = list(
    "Disease covariance and correlation", build_dependence_parameters
  ),
  adjacency_cardinality = list(
    "Adjacency cardinality", build_adjacency_cardinality
  ),
  boundary_count = list("Draw-level boundary counts", build_boundary_count),
  centered_phi = list("Centered spatial effects", build_centered_phi),
  log_risk = list(
    "Fitted log-risks (equivalently, relative risks)", build_log_risk
  ),
  latent_r = list("Latent Gaussian r variables", build_latent_r)
)

diagnostic_results <- setNames(vector("list", length(diagnostic_specifications)),
                               names(diagnostic_specifications))
for (group_key in names(diagnostic_specifications)) {
  specification <- diagnostic_specifications[[group_key]]
  diagnostic_results[[group_key]] <- diagnose_draw_group(
    group_key, specification[[1L]], specification[[2L]]
  )
}

diagnostic_overview <- do.call(
  rbind, lapply(diagnostic_results, `[[`, "overview")
)
posterior_summary_all_groups <- do.call(
  rbind, lapply(diagnostic_results, `[[`, "posterior_summary")
)
posterior_agreement_all_groups <- do.call(
  rbind, lapply(diagnostic_results, `[[`, "posterior_agreement")
)
posterior_agreement_overview <- do.call(rbind, lapply(split(
  posterior_agreement_all_groups,
  posterior_agreement_all_groups$group
), function(group) {
  finite_range <- is.finite(group$standardized_mean_range)
  data.frame(
    group = group$group[[1L]],
    variables = nrow(group),
    max_standardized_mean_range = if (any(finite_range)) {
      max(group$standardized_mean_range[finite_range])
    } else {
      NA_real_
    },
    credible_intervals_all_overlap = sum(
      group$credible_intervals_all_overlap, na.rm = TRUE
    ),
    credible_intervals_not_all_overlap = sum(
      !group$credible_intervals_all_overlap, na.rm = TRUE
    ),
    row.names = NULL
  )
}))

# Proposal acceptance, state integrity, and optional tempering traffic ---------

acceptance_summary <- do.call(rbind, lapply(seq_len(n_chains), function(chain) {
  fit <- read_chain(chain)
  acceptance <- fit$acceptance
  data.frame(
    chain = chain_names[[chain]],
    beta = acceptance$beta_rate,
    theta = acceptance$theta_rate,
    r = acceptance$r_rate,
    V = acceptance$V_rate,
    rho = acceptance$rho_rate,
    eta = acceptance$eta_rate,
    A = acceptance$A_rate,
    row.names = NULL
  )
}))

state_consistency_summary <- do.call(
  rbind, lapply(seq_len(n_chains), function(chain) {
    fit <- read_chain(chain)
    checks <- fit$state_consistency
    data.frame(
      chain = chain_names[[chain]],
      checks = checks$checks,
      tolerance = checks$tolerance,
      maximum_F_gap_seen = checks$maximum_F_gap_seen,
      all_checks_passed = checks$all_checks_passed,
      row.names = NULL
    )
  })
)

replica_traffic_summary <- if (is_tempering_workflow) {
  do.call(
    rbind, lapply(seq_len(n_chains), function(chain) {
      diagnostics <- readRDS(ladder_diagnostic_files[[chain]])
      occupancy <- diagnostics$replica_temperature_occupancy
      levels_visited <- rowSums(occupancy > 0L)
      data.frame(
        chain = chain_names[[chain]],
        replicas_visiting_cold = sum(occupancy[, 1L] > 0L),
        replicas_visiting_hot = sum(occupancy[, ncol(occupancy)] > 0L),
        replicas_visiting_both = sum(
          occupancy[, 1L] > 0L & occupancy[, ncol(occupancy)] > 0L
        ),
        replicas_visiting_all_temperatures = sum(
          levels_visited == ncol(occupancy)
        ),
        minimum_levels_visited = min(levels_visited),
        median_levels_visited = median(levels_visited),
        maximum_levels_visited = max(levels_visited),
        total_round_trips = diagnostics$total_round_trips,
        minimum_sampling_swap = min(diagnostics$sampling_swap_acceptance),
        mean_sampling_swap = mean(diagnostics$sampling_swap_acceptance),
        maximum_sampling_swap = max(diagnostics$sampling_swap_acceptance),
        row.names = NULL
      )
    })
  )
} else {
  NULL
}

# Posterior boundary probabilities and FDR selection ---------------------------

cat("Computing edge-level boundary probabilities and FDR selections...\n")
boundary_probabilities <- array(
  NA_real_,
  dim = c(n_edges, n_diseases, n_chains),
  dimnames = list(
    boundary = neighbor_name, disease = disease_names, chain = chain_names
  )
)
pair_indices <- t(utils::combn(seq_len(n_diseases), 2L))
pair_names <- apply(pair_indices, 1L, function(pair) {
  paste(disease_names[pair], collapse = ", ")
})
n_pairs <- nrow(pair_indices)
shared_boundary_probability_sum <- matrix(
  0, nrow = n_edges, ncol = n_pairs,
  dimnames = list(boundary = neighbor_name, pair = pair_names)
)
mutual_boundary_probability_sum <- matrix(
  0, nrow = n_edges, ncol = n_pairs,
  dimnames = list(boundary = neighbor_name, pair = pair_names)
)
pooled_phi_sum <- matrix(
  0, nrow = n_counties, ncol = n_diseases,
  dimnames = list(county = seq_len(n_counties), disease = disease_names)
)
chain_phi_means <- array(
  NA_real_, dim = c(n_counties, n_diseases, n_chains),
  dimnames = list(
    county = seq_len(n_counties), disease = disease_names, chain = chain_names
  )
)
pooled_W_sum <- lapply(seq_len(n_diseases), function(disease) {
  matrix(0, nrow = n_counties, ncol = n_counties)
})
total_posterior_draws <- 0
cardinality_summary_by_chain <- list()
cardinality_row <- 1L
for (chain in seq_len(n_chains)) {
  fit <- read_chain(chain)
  chain_draws <- nrow(fit$phi)
  phi_by_disease <- vector("list", n_diseases)
  boundary_draws_by_disease <- vector("list", n_diseases)
  for (disease in seq_len(n_diseases)) {
    indices <- ((disease - 1L) * n_counties + 1L):(disease * n_counties)
    phi <- fit$phi[, indices, drop = FALSE]
    phi <- phi[, order(final_perm), drop = FALSE]
    boundary_draws <-
      phi[, neighbor_list0[, 1L], drop = FALSE] !=
      phi[, neighbor_list0[, 2L], drop = FALSE]
    phi_by_disease[[disease]] <- phi
    boundary_draws_by_disease[[disease]] <- boundary_draws
    boundary_probabilities[, disease, chain] <- colMeans(boundary_draws)
    chain_phi_means[, disease, chain] <- colMeans(phi)
    pooled_phi_sum[, disease] <- pooled_phi_sum[, disease] + colSums(phi)
    pooled_W_sum[[disease]] <- pooled_W_sum[[disease]] +
      chain_draws * fit$W_mean[[disease]]
    cardinality <- fit$W_cardinality[, disease]
    cardinality_quantiles <- stats::quantile(
      cardinality, c(0.05, 0.50, 0.95), names = FALSE
    )
    cardinality_summary_by_chain[[cardinality_row]] <- data.frame(
      chain = chain_names[[chain]],
      disease = disease_names[[disease]],
      mean_adjacency_edges = mean(cardinality),
      sd_adjacency_edges = stats::sd(cardinality),
      q05_adjacency_edges = cardinality_quantiles[[1L]],
      median_adjacency_edges = cardinality_quantiles[[2L]],
      q95_adjacency_edges = cardinality_quantiles[[3L]],
      row.names = NULL
    )
    cardinality_row <- cardinality_row + 1L
  }

  for (pair_index in seq_len(n_pairs)) {
    first_disease <- pair_indices[pair_index, 1L]
    second_disease <- pair_indices[pair_index, 2L]

    shared_boundary_probability_sum[, pair_index] <-
      shared_boundary_probability_sum[, pair_index] + colSums(
        boundary_draws_by_disease[[first_disease]] &
          boundary_draws_by_disease[[second_disease]]
      )

    first_phi <- phi_by_disease[[first_disease]]
    second_phi <- phi_by_disease[[second_disease]]
    mutual_boundary_probability_sum[, pair_index] <-
      mutual_boundary_probability_sum[, pair_index] + colSums(
        first_phi[, neighbor_list0[, 1L], drop = FALSE] !=
          second_phi[, neighbor_list0[, 2L], drop = FALSE] &
          second_phi[, neighbor_list0[, 1L], drop = FALSE] !=
          first_phi[, neighbor_list0[, 2L], drop = FALSE]
      )
  }

  total_posterior_draws <- total_posterior_draws + chain_draws
  rm(
    fit, phi_by_disease, boundary_draws_by_disease, phi, boundary_draws,
    first_phi, second_phi
  )
  invisible(gc(FALSE))
}
cardinality_summary_by_chain <- do.call(rbind, cardinality_summary_by_chain)

pooled_phi_mean <- pooled_phi_sum / total_posterior_draws
shared_boundary_probabilities <-
  shared_boundary_probability_sum / total_posterior_draws
mutual_boundary_probabilities <-
  mutual_boundary_probability_sum / total_posterior_draws
pooled_W_mean <- lapply(pooled_W_sum, function(W_sum) {
  W_mean <- W_sum / total_posterior_draws
  W_mean[order(final_perm), order(final_perm), drop = FALSE]
})
names(pooled_W_mean) <- disease_names
non_adjacency_probabilities <- vapply(
  seq_len(n_diseases),
  function(disease) {
    probabilities <- 1 - pooled_W_mean[[disease]][
      cbind(neighbor_list0[, 1L], neighbor_list0[, 2L])
    ]
    pmin(1, pmax(0, probabilities))
  },
  numeric(n_edges)
)
dimnames(non_adjacency_probabilities) <- list(
  boundary = neighbor_name, disease = disease_names
)

boundary_probability_mean <- apply(boundary_probabilities, c(1L, 2L), mean)
boundary_probability_range <- apply(
  boundary_probabilities, c(1L, 2L), function(x) diff(range(x))
)
boundary_probability_mean_pairwise_absolute_difference <- apply(
  boundary_probabilities, c(1L, 2L), function(probabilities) {
    differences <- abs(outer(probabilities, probabilities, "-"))
    mean(differences[upper.tri(differences)])
  }
)

fdr_curve <- function(probabilities) {
  thresholds <- sort(unique(probabilities[is.finite(probabilities)]),
                     decreasing = TRUE)
  data.frame(
    threshold = thresholds,
    selected_edges = vapply(
      thresholds, function(threshold) sum(probabilities >= threshold),
      integer(1L)
    ),
    estimated_fdr = vapply(
      thresholds,
      function(threshold) mean(1 - probabilities[probabilities >= threshold]),
      numeric(1L)
    )
  )
}

boundary_fdr_selection <- function(probabilities, alpha = fdr_alpha) {
  probabilities <- probabilities[is.finite(probabilities)]
  curve <- fdr_curve(probabilities)
  valid <- which(curve$estimated_fdr <= alpha)
  if (!length(valid)) {
    return(list(
      selected_edges = 0L, threshold = NA_real_, estimated_fdr = NA_real_,
      cutoff_ties = 0L, mean_selected_probability = NA_real_,
      min_selected_probability = NA_real_, max_probability = max(probabilities),
      selected_boundaries = character(0), curve = curve
    ))
  }
  selected_index <- valid[which.max(curve$selected_edges[valid])]
  threshold <- curve$threshold[[selected_index]]
  selected_probabilities <- sort(
    probabilities[probabilities >= threshold], decreasing = TRUE
  )
  list(
    selected_edges = length(selected_probabilities),
    threshold = threshold,
    estimated_fdr = curve$estimated_fdr[[selected_index]],
    cutoff_ties = sum(probabilities == threshold),
    mean_selected_probability = mean(selected_probabilities),
    min_selected_probability = min(selected_probabilities),
    max_probability = max(probabilities),
    selected_boundaries = names(selected_probabilities),
    curve = curve
  )
}

selected_boundary_sets <- setNames(
  lapply(seq_len(n_diseases), function(x) {
    setNames(vector("list", n_chains), chain_names)
  }),
  disease_names
)
chain_boundary_rows <- list()
fdr_curve_rows <- list()
row_index <- 1L
curve_index <- 1L
for (chain in seq_len(n_chains)) {
  for (disease in seq_len(n_diseases)) {
    probabilities <- boundary_probabilities[, disease, chain]
    names(probabilities) <- neighbor_name
    selection <- boundary_fdr_selection(probabilities)
    selected_boundary_sets[[disease_names[[disease]]]][[
      chain_names[[chain]]
    ]] <- selection$selected_boundaries
    chain_boundary_rows[[row_index]] <- data.frame(
      chain = chain_names[[chain]],
      disease = disease_names[[disease]],
      selected_edges = selection$selected_edges,
      threshold = selection$threshold,
      estimated_fdr = selection$estimated_fdr,
      cutoff_ties = selection$cutoff_ties,
      mean_selected_probability = selection$mean_selected_probability,
      min_selected_probability = selection$min_selected_probability,
      max_probability = selection$max_probability,
      boundaries_prob_ge_0.95 = sum(probabilities >= 0.95),
      boundaries_prob_ge_0.90 = sum(probabilities >= 0.90),
      row.names = NULL
    )
    curve <- selection$curve
    curve$chain <- chain_names[[chain]]
    curve$disease <- disease_names[[disease]]
    curve$pooled <- FALSE
    fdr_curve_rows[[curve_index]] <- curve
    row_index <- row_index + 1L
    curve_index <- curve_index + 1L
  }
}
boundary_count_summary <- do.call(rbind, chain_boundary_rows)

pooled_selected_boundary_sets <- setNames(
  vector("list", n_diseases), disease_names
)
pooled_boundary_rows <- list()
for (disease in seq_len(n_diseases)) {
  probabilities <- boundary_probability_mean[, disease]
  names(probabilities) <- neighbor_name
  selection <- boundary_fdr_selection(probabilities)
  pooled_selected_boundary_sets[[disease_names[[disease]]]] <-
    selection$selected_boundaries
  pooled_boundary_rows[[disease]] <- data.frame(
    chain = "pooled_mean",
    disease = disease_names[[disease]],
    selected_edges = selection$selected_edges,
    threshold = selection$threshold,
    estimated_fdr = selection$estimated_fdr,
    cutoff_ties = selection$cutoff_ties,
    mean_selected_probability = selection$mean_selected_probability,
    min_selected_probability = selection$min_selected_probability,
    max_probability = selection$max_probability,
    boundaries_prob_ge_0.95 = sum(probabilities >= 0.95),
    boundaries_prob_ge_0.90 = sum(probabilities >= 0.90),
    row.names = NULL
  )
  curve <- selection$curve
  curve$chain <- "pooled_mean"
  curve$disease <- disease_names[[disease]]
  curve$pooled <- TRUE
  fdr_curve_rows[[curve_index]] <- curve
  curve_index <- curve_index + 1L
}
pooled_boundary_count_summary <- do.call(rbind, pooled_boundary_rows)
fdr_curves <- do.call(rbind, fdr_curve_rows)
fdr_curves$disease <- factor(fdr_curves$disease, levels = disease_names)

summarize_pair_boundary_selection <- function(probability_matrix, boundary_type) {
  selected_sets <- setNames(vector("list", n_pairs), pair_names)
  count_rows <- vector("list", n_pairs)
  curve_rows <- vector("list", n_pairs)
  for (pair_index in seq_len(n_pairs)) {
    probabilities <- probability_matrix[, pair_index]
    names(probabilities) <- neighbor_name
    selection <- boundary_fdr_selection(probabilities)
    selected_sets[[pair_names[[pair_index]]]] <- selection$selected_boundaries
    count_rows[[pair_index]] <- data.frame(
      boundary_type = boundary_type,
      cancer_pair = pair_names[[pair_index]],
      selected_edges = selection$selected_edges,
      threshold = selection$threshold,
      estimated_fdr = selection$estimated_fdr,
      cutoff_ties = selection$cutoff_ties,
      mean_selected_probability = selection$mean_selected_probability,
      min_selected_probability = selection$min_selected_probability,
      max_probability = selection$max_probability,
      row.names = NULL
    )
    curve <- selection$curve
    curve$boundary_type <- boundary_type
    curve$cancer_pair <- pair_names[[pair_index]]
    curve_rows[[pair_index]] <- curve
  }
  list(
    selected_sets = selected_sets,
    count_summary = do.call(rbind, count_rows),
    fdr_curves = do.call(rbind, curve_rows)
  )
}

shared_boundary_results <- summarize_pair_boundary_selection(
  shared_boundary_probabilities, "shared"
)
mutual_boundary_results <- summarize_pair_boundary_selection(
  mutual_boundary_probabilities, "mutual"
)
pair_fdr_curves <- rbind(
  shared_boundary_results$fdr_curves,
  mutual_boundary_results$fdr_curves
)

boundary_count_ranges <- do.call(rbind, lapply(disease_names, function(disease) {
  counts <- boundary_count_summary$selected_edges[
    boundary_count_summary$disease == disease
  ]
  pooled_count <- pooled_boundary_count_summary$selected_edges[
    pooled_boundary_count_summary$disease == disease
  ]
  data.frame(
    disease = disease,
    min_chain_selected_edges = min(counts),
    median_chain_selected_edges = median(counts),
    max_chain_selected_edges = max(counts),
    pooled_selected_edges = pooled_count,
    max_minus_min_chain = max(counts) - min(counts),
    row.names = NULL
  )
}))

boundary_jaccard <- function(first, second) {
  union_size <- length(union(first, second))
  if (!union_size) return(1)
  length(intersect(first, second)) / union_size
}

boundary_pairwise_jaccard <- do.call(rbind, lapply(
  disease_names, function(disease) {
    do.call(rbind, utils::combn(chain_names, 2L, function(pair) {
      first <- selected_boundary_sets[[disease]][[pair[[1L]]]]
      second <- selected_boundary_sets[[disease]][[pair[[2L]]]]
      data.frame(
        disease = disease, chain_1 = pair[[1L]], chain_2 = pair[[2L]],
        jaccard = boundary_jaccard(first, second),
        common_selected_edges = length(intersect(first, second)),
        union_selected_edges = length(union(first, second)),
        row.names = NULL
      )
    }, simplify = FALSE))
  }
))
boundary_pairwise_jaccard_summary <- do.call(
  rbind, lapply(disease_names, function(disease) {
    values <- boundary_pairwise_jaccard$jaccard[
      boundary_pairwise_jaccard$disease == disease
    ]
    data.frame(
      disease = disease,
      min_pairwise_jaccard = min(values),
      median_pairwise_jaccard = median(values),
      mean_pairwise_jaccard = mean(values),
      max_pairwise_jaccard = max(values),
      row.names = NULL
    )
  })
)

safe_correlation <- function(first, second) {
  if (stats::sd(first) == 0 || stats::sd(second) == 0) {
    return(if (isTRUE(all.equal(first, second))) 1 else NA_real_)
  }
  stats::cor(first, second)
}

boundary_probability_pairwise <- do.call(rbind, lapply(
  seq_len(n_diseases), function(disease) {
    do.call(rbind, utils::combn(seq_len(n_chains), 2L, function(pair) {
      first <- boundary_probabilities[, disease, pair[[1L]]]
      second <- boundary_probabilities[, disease, pair[[2L]]]
      difference <- first - second
      data.frame(
        disease = disease_names[[disease]],
        chain_1 = chain_names[[pair[[1L]]]],
        chain_2 = chain_names[[pair[[2L]]]],
        correlation = safe_correlation(first, second),
        mean_absolute_difference = mean(abs(difference)),
        root_mean_squared_difference = sqrt(mean(difference^2)),
        max_absolute_difference = max(abs(difference)),
        row.names = NULL
      )
    }, simplify = FALSE))
  }
))
boundary_probability_pairwise_summary <- do.call(
  rbind, lapply(disease_names, function(disease) {
    values <- boundary_probability_pairwise[
      boundary_probability_pairwise$disease == disease, , drop = FALSE
    ]
    data.frame(
      disease = disease,
      min_pairwise_correlation = min(values$correlation, na.rm = TRUE),
      median_pairwise_correlation = median(values$correlation, na.rm = TRUE),
      mean_pairwise_correlation = mean(values$correlation, na.rm = TRUE),
      max_pairwise_mean_absolute_difference = max(
        values$mean_absolute_difference
      ),
      max_pairwise_edge_difference = max(values$max_absolute_difference),
      row.names = NULL
    )
  })
)

boundary_agreement <- data.frame(
  boundary = rep(neighbor_name, times = n_diseases),
  disease = rep(disease_names, each = n_edges),
  pooled_probability = as.vector(boundary_probability_mean),
  between_chain_range = as.vector(boundary_probability_range),
  mean_pairwise_absolute_difference = as.vector(
    boundary_probability_mean_pairwise_absolute_difference
  ),
  row.names = NULL
)
for (chain in seq_len(n_chains)) {
  boundary_agreement[[chain_names[[chain]]]] <- mapply(
    function(boundary, disease) {
      boundary_probabilities[boundary, disease, chain]
    },
    boundary_agreement$boundary, boundary_agreement$disease
  )
}
boundary_agreement <- boundary_agreement[
  order(boundary_agreement$between_chain_range, decreasing = TRUE),
]

boundary_probability_range_summary <- do.call(
  rbind, lapply(seq_len(n_diseases), function(disease) {
    ranges <- boundary_probability_range[, disease]
    pairwise_differences <-
      boundary_probability_mean_pairwise_absolute_difference[, disease]
    data.frame(
      disease = disease_names[[disease]],
      median_between_chain_range = median(ranges),
      mean_between_chain_range = mean(ranges),
      q90_between_chain_range = stats::quantile(ranges, 0.90),
      max_between_chain_range = max(ranges),
      median_pairwise_absolute_difference = median(pairwise_differences),
      mean_pairwise_absolute_difference = mean(pairwise_differences),
      boundaries_range_ge_0.25 = sum(ranges >= 0.25),
      boundaries_range_ge_0.50 = sum(ranges >= 0.50),
      boundaries_range_ge_0.75 = sum(ranges >= 0.75),
      row.names = NULL
    )
  })
)

boundary_selection_incidence <- do.call(rbind, lapply(
  disease_names, function(disease) {
    incidence <- vapply(chain_names, function(chain) {
      neighbor_name %in% selected_boundary_sets[[disease]][[chain]]
    }, logical(n_edges))
    colnames(incidence) <- chain_names
    selected_by_chains <- rowSums(incidence)
    disease_index <- match(disease, disease_names)
    data.frame(
      boundary = neighbor_name,
      disease = disease,
      incidence,
      selected_by_chains = selected_by_chains,
      selected_by_all_chains = selected_by_chains == n_chains,
      selected_by_majority = selected_by_chains >= floor(n_chains / 2L) + 1L,
      selected_by_one_chain = selected_by_chains == 1L,
      pooled_probability = boundary_probability_mean[, disease_index],
      pooled_selected = neighbor_name %in%
        pooled_selected_boundary_sets[[disease]],
      row.names = NULL,
      check.names = FALSE
    )
  }
))
boundary_selection_frequency_summary <- do.call(
  rbind, lapply(disease_names, function(disease) {
    incidence <- boundary_selection_incidence[
      boundary_selection_incidence$disease == disease, , drop = FALSE
    ]
    data.frame(
      disease = disease,
      boundaries_selected_by_all_chains = sum(
        incidence$selected_by_all_chains
      ),
      boundaries_selected_by_majority = sum(incidence$selected_by_majority),
      boundaries_selected_by_one_chain = sum(incidence$selected_by_one_chain),
      boundaries_selected_by_no_chain = sum(incidence$selected_by_chains == 0L),
      pooled_boundaries_selected = sum(incidence$pooled_selected),
      row.names = NULL
    )
  })
)

boundary_cutoff_sensitivity <- do.call(
  rbind, lapply(disease_names, function(disease) {
    pooled_row <- pooled_boundary_count_summary[
      pooled_boundary_count_summary$disease == disease, , drop = FALSE
    ]
    incidence <- boundary_selection_incidence[
      boundary_selection_incidence$disease == disease, , drop = FALSE
    ]
    incidence$pooled_fdr_threshold <- pooled_row$threshold[[1L]]
    incidence$absolute_distance_from_threshold <- abs(
      incidence$pooled_probability - pooled_row$threshold[[1L]]
    )
    head(
      incidence[order(incidence$absolute_distance_from_threshold), ],
      20L
    )
  })
)

# RESULT: Supplementary Figures S2, S12, and S21 - pooled FDR curves ----------
# The specification determines the manuscript number: Adjacency = S2,
# Mean--Adjacency = S12, and Mean = S21. The chain-specific curves saved below
# are diagnostics and are not manuscript figures.

pooled_fdr_curves <- fdr_curves[fdr_curves$pooled, , drop = FALSE]
pooled_fdr_plot <- ggplot2::ggplot(
  pooled_fdr_curves,
  ggplot2::aes(x = selected_edges, y = estimated_fdr, linetype = disease)
) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_hline(
    yintercept = fdr_alpha, color = "#B2182B", linewidth = 0.8
  ) +
  ggplot2::scale_linetype_manual(
    values = c("solid", "dashed", "dotted", "dotdash"), name = "Cancer"
  ) +
  ggplot2::labs(
    x = "Number of edges selected",
    y = "Estimated FDR"
  ) +
  ggplot2::theme_bw(base_size = 14)
pooled_fdr_file <- file.path(
  figures_dir, paste0("pooled_fdr_curves_", cvrts, ".pdf")
)
ggplot2::ggsave(pooled_fdr_file, pooled_fdr_plot, width = 9, height = 6.5)

chain_fdr_plot <- ggplot2::ggplot(
  fdr_curves,
  ggplot2::aes(
    x = selected_edges, y = estimated_fdr,
    group = chain, color = chain,
    linewidth = pooled, alpha = pooled
  )
) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(
    yintercept = fdr_alpha, color = "#B2182B", linewidth = 0.7
  ) +
  ggplot2::facet_wrap(~disease, nrow = 2L, scales = "free_x") +
  ggplot2::scale_linewidth_manual(values = c(`FALSE` = 0.45, `TRUE` = 1.2)) +
  ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.65, `TRUE` = 1)) +
  ggplot2::labs(
    x = "Number of selected boundaries",
    y = "Estimated posterior FDR",
    color = "Analysis",
    title = "Chain-specific and pooled posterior-FDR curves"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::guides(linewidth = "none", alpha = "none")
chain_fdr_file <- file.path(
  figures_dir, paste0("chain_specific_fdr_curves_", cvrts, ".pdf")
)
ggplot2::ggsave(chain_fdr_file, chain_fdr_plot, width = 11, height = 8)

# RESULTS: Supplementary posterior interval figures ---------------------------
# These are generated for all adaptive-tempering specifications. The
# exact manuscript mapping is:
#   beta/theta:       S3 (Adjacency), S8 (Mean--Adjacency), S17 (Mean)
#   gamma:            S4 (Adjacency), S9 (Mean--Adjacency), S18 (Mean)
#   V/rho:            S5 (Adjacency), S10 (Mean--Adjacency), S19 (Mean)
#   eta/covariance:   S6 (Adjacency), S11 (Mean--Adjacency)
#   covariance only:  S20 (Mean)
# The plotted intervals pool the retained cold-posterior draws from all six
# independently initialized ladders. In the sampler output, r is the latent
# Gaussian vector denoted by gamma in the manuscript.

supplement_figure_files <- character(0)
supplement_table_files <- character(0)
if (is_tempering_workflow) {
  cat(
    "Generating pooled supplementary posterior plots for ",
    tempering_specification, "...\n", sep = ""
  )
  dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)

  supplement_fdr_file <- file.path(
    supplement_dir,
    paste0("pooled_fdr_curves_", tempering_specification, ".pdf")
  )
  ggplot2::ggsave(
    supplement_fdr_file, pooled_fdr_plot, width = 9, height = 6.5
  )

  supplement_fields <- c(
    "beta", "theta", "tau", "r", "V", "rho",
    if (cvrts != "mean") "eta" else NULL,
    "A"
  )
  supplement_blocks <- setNames(
    lapply(
      supplement_fields,
      function(field) vector("list", length(chain_files))
    ),
    supplement_fields
  )
  for (chain in seq_len(n_chains)) {
    fit <- read_chain(chain)
    for (field in supplement_fields) {
      supplement_blocks[[field]][[chain]] <- as.matrix(fit[[field]])
    }
    rm(fit)
    invisible(gc(FALSE))
  }
  supplement_draws <- lapply(
    supplement_blocks,
    function(field_blocks) do.call(rbind, field_blocks)
  )

  supplement_interval_summary <- function(x) {
    x <- as.matrix(x)
    quantiles <- matrixStats::colQuantiles(
      x,
      probs = c(0.025, 0.05, 0.50, 0.95, 0.975),
      drop = FALSE
    )
    data.frame(
      q025 = quantiles[, 1L],
      q05 = quantiles[, 2L],
      median = quantiles[, 3L],
      q95 = quantiles[, 4L],
      q975 = quantiles[, 5L]
    )
  }

  plot_supplement_interval_panel <- function(
      x, labels, main = "", xlab = "Posterior interval",
      order_by_median = TRUE) {
    summary <- supplement_interval_summary(x)
    if (length(labels) != nrow(summary)) {
      stop("The number of labels does not match the number of parameters.")
    }
    ordering <- if (order_by_median) {
      order(summary$median)
    } else {
      seq_len(nrow(summary))
    }
    summary <- summary[ordering, , drop = FALSE]
    labels <- labels[ordering]
    y <- seq_len(nrow(summary))
    limits <- range(summary$q025, summary$q975, 0, finite = TRUE)
    padding <- max(diff(limits) * 0.04, .Machine$double.eps)

    graphics::plot(
      NA_real_, NA_real_,
      xlim = limits + c(-padding, padding),
      ylim = c(0.5, nrow(summary) + 0.5),
      yaxt = "n", ylab = "", xlab = xlab, main = main, bty = "l"
    )
    graphics::abline(v = 0, col = "grey80", lty = 2L)
    graphics::segments(
      summary$q025, y, summary$q975, y, col = "grey45", lwd = 1
    )
    graphics::segments(
      summary$q05, y, summary$q95, y, col = "black", lwd = 3
    )
    graphics::points(summary$median, y, pch = 16, cex = 0.65)
    graphics::axis(
      2L, at = y, labels = labels, las = 1L, cex.axis = 0.72
    )
    invisible(summary)
  }

  # RESULT: Supplementary Tables S16--S18 - pooled posterior summaries --------
  # Adjacency = S16, Mean--Adjacency = S17, and Mean = S18.
  table_variable_names <- c(
    paste0("beta[", seq_len(ncol(supplement_draws$beta)), "]"),
    paste0("theta[", seq_len(ncol(supplement_draws$theta)), "]"),
    "tau"
  )
  table_latex_names <- c(
    paste0("$\\beta_{", seq_len(ncol(supplement_draws$beta)), "}$"),
    paste0("$\\theta_{", seq_len(ncol(supplement_draws$theta)), "}$"),
    "$\\tau$"
  )
  table_draw_array <- array(
    NA_real_,
    dim = c(retained_draws, n_chains, length(table_variable_names)),
    dimnames = list(
      iteration = NULL,
      chain = chain_names,
      variable = table_variable_names
    )
  )
  for (chain in seq_len(n_chains)) {
    table_draw_array[, chain, ] <- cbind(
      supplement_blocks$beta[[chain]],
      supplement_blocks$theta[[chain]],
      supplement_blocks$tau[[chain]]
    )
  }
  table_summary <- as.data.frame(posterior::summarise_draws(
    posterior::as_draws_array(table_draw_array),
    "mean", "sd", "mcse_mean"
  ))
  table_summary <- table_summary[
    match(table_variable_names, table_summary$variable),
    c("variable", "mean", "sd", "mcse_mean"),
    drop = FALSE
  ]
  table_summary$latex_variable <- table_latex_names

  table_number <- switch(
    tempering_specification,
    adj = "S16",
    meanadj = "S17",
    mean = "S18"
  )
  table_stem <- paste0(
    "pooled_table_", table_number, "_", tempering_specification
  )
  table_csv_file <- file.path(supplement_dir, paste0(table_stem, ".csv"))
  table_tex_file <- file.path(supplement_dir, paste0(table_stem, ".tex"))
  utils::write.csv(
    table_summary[, c("variable", "mean", "sd", "mcse_mean")],
    table_csv_file,
    row.names = FALSE
  )

  specification_label <- switch(
    tempering_specification,
    adj = "Adjacency",
    meanadj = "Mean--Adjacency",
    mean = "Mean"
  )
  table_label <- switch(
    tempering_specification,
    adj = "tab:post_inf",
    meanadj = "tab:post_inf_meanadj",
    mean = "tab:post_inf_mean"
  )
  table_rows <- vapply(seq_len(nrow(table_summary)), function(row) {
    paste0(
      table_summary$latex_variable[[row]], " & ",
      sprintf("%.3f", table_summary$mean[[row]]), " & ",
      sprintf("%.3f", table_summary$sd[[row]]), " & ",
      sprintf("%.3f", table_summary$mcse_mean[[row]]), " \\\\"
    )
  }, character(1L))
  table_caption <- paste0(
    "Posterior estimates, posterior standard deviations, and multiple-chain ",
    "Monte Carlo standard errors for the regression coefficients, atoms, ",
    "and their precision under ", specification_label, ". Summaries pool ",
    format(total_posterior_draws, big.mark = ",", scientific = FALSE),
    " retained cold-posterior draws from six adaptive-tempering runs."
  )
  writeLines(
    c(
      "\\begin{table}[!htbp]",
      "\\centering",
      "\\small",
      paste0("\\caption{", table_caption, "}"),
      paste0("\\label{", table_label, "}"),
      "\\begin{tabular}{lrrr}",
      "\\hline",
      "Variable & Estimate & Standard deviation & MCSE \\\\",
      "\\hline",
      table_rows,
      "\\hline",
      "\\end{tabular}",
      "\\end{table}"
    ),
    table_tex_file,
    useBytes = TRUE
  )
  supplement_table_files <- c(
    posterior_table_csv = table_csv_file,
    posterior_table_tex = table_tex_file
  )

  # RESULT: Supplementary Figures S3, S8, and S17 - beta and theta ------------
  figure_S3_file <- file.path(
    supplement_dir,
    paste0("pooled_beta_theta_", tempering_specification, ".pdf")
  )
  grDevices::pdf(figure_S3_file, width = 10, height = 6, onefile = FALSE)
  graphics::par(mfrow = c(1L, 2L), mar = c(4.2, 4.8, 2.0, 1.0))
  plot_supplement_interval_panel(
    supplement_draws$beta,
    paste0("beta[", seq_len(ncol(supplement_draws$beta)), "]"),
    "beta"
  )
  plot_supplement_interval_panel(
    supplement_draws$theta,
    paste0("theta[", seq_len(ncol(supplement_draws$theta)), "]"),
    "theta"
  )
  grDevices::dev.off()

  # RESULT: Supplementary Figures S4, S9, and S18 - latent gamma --------------
  figure_S4_file <- file.path(
    supplement_dir,
    paste0("pooled_gamma_", tempering_specification, ".pdf")
  )
  grDevices::pdf(figure_S4_file, width = 10, height = 10, onefile = FALSE)
  graphics::par(mfrow = c(2L, 2L), mar = c(4.2, 4.3, 2.1, 0.8))
  latent_counties <- ncol(supplement_draws$r) %/% n_diseases
  for (disease in seq_len(n_diseases)) {
    columns <- ((disease - 1L) * latent_counties + 1L):
      (disease * latent_counties)
    plot_supplement_interval_panel(
      supplement_draws$r[, columns, drop = FALSE],
      paste0("gamma[", seq_len(latent_counties), "]"),
      disease_names[[disease]]
    )
  }
  grDevices::dev.off()

  # RESULT: Supplementary Figures S5, S10, and S19 - V and rho ----------------
  figure_S5_file <- file.path(
    supplement_dir,
    paste0("pooled_V_rho_", tempering_specification, ".pdf")
  )
  grDevices::pdf(figure_S5_file, width = 10, height = 6, onefile = FALSE)
  graphics::par(mfrow = c(1L, 2L), mar = c(4.2, 4.8, 2.0, 1.0))
  plot_supplement_interval_panel(
    supplement_draws$V,
    paste0("V[", seq_len(ncol(supplement_draws$V)), "]"),
    "V"
  )
  plot_supplement_interval_panel(
    supplement_draws$rho, disease_names, "rho"
  )
  grDevices::dev.off()

  # RESULT: Supplementary Figures S6, S11, and S20 - covariance summaries -----
  # Adjacency and Mean--Adjacency also display eta; Mean displays only
  # (A A^T)_{ij}, because its adjacency matrix is fixed.
  covariance_draws <- matrix(
    NA_real_,
    nrow = nrow(supplement_draws$A),
    ncol = nrow(covariance_indices)
  )
  for (draw in seq_len(nrow(supplement_draws$A))) {
    A_matrix <- A_vector_to_matrix(supplement_draws$A[draw, ])
    AA_t <- A_matrix %*% t(A_matrix)
    covariance_draws[draw, ] <- AA_t[covariance_indices]
  }
  figure_S6_file <- file.path(
    supplement_dir,
    if (cvrts != "mean") {
      paste0("pooled_eta_covariance_", tempering_specification, ".pdf")
    } else {
      paste0("pooled_covariance_", tempering_specification, ".pdf")
    }
  )
  has_adjacency_coefficients <- cvrts != "mean"
  grDevices::pdf(
    figure_S6_file,
    width = if (has_adjacency_coefficients) 12 else 7,
    height = 8,
    onefile = FALSE
  )
  graphics::par(
    mfrow = if (has_adjacency_coefficients) c(1L, 2L) else c(1L, 1L),
    mar = c(4.2, 5.0, 2.0, 1.0)
  )
  if (has_adjacency_coefficients) {
    plot_supplement_interval_panel(
      supplement_draws$eta,
      paste0("eta[", seq_len(ncol(supplement_draws$eta)), "]"),
      "eta"
    )
  }
  AA_t_labels <- sprintf(
    "(AA^T)[%d,%d]",
    covariance_indices[, "row"], covariance_indices[, "col"]
  )
  plot_supplement_interval_panel(
    covariance_draws, AA_t_labels, expression(A * A^T)
  )
  grDevices::dev.off()

  supplement_figure_files <- c(
    pooled_fdr = supplement_fdr_file,
    figure_S3 = figure_S3_file,
    figure_S4 = figure_S4_file,
    figure_S5 = figure_S5_file,
    figure_S6 = figure_S6_file
  )
  rm(
    supplement_draws, supplement_blocks, supplement_fields,
    covariance_draws, A_matrix, AA_t, AA_t_labels, table_draw_array,
    table_summary, table_rows
  )
  invisible(gc(FALSE))
}

# Boundary maps ----------------------------------------------------------------

cat("Constructing California boundary geometries and paper-style maps...\n")
county_geometry <- sf::st_geometry(ca.poly)
edge_geometries <- vector("list", n_edges)
valid_edge_geometry <- logical(n_edges)
for (edge in seq_len(n_edges)) {
  first_county <- neighbor_list0[edge, 1L]
  second_county <- neighbor_list0[edge, 2L]
  shared_border <- suppressWarnings(sf::st_intersection(
    sf::st_boundary(county_geometry[first_county]),
    sf::st_boundary(county_geometry[second_county])
  ))
  if (length(shared_border) && !all(sf::st_is_empty(shared_border))) {
    linear_part <- shared_border[sf::st_dimension(shared_border) == 1L]
    if (length(linear_part) && !all(sf::st_is_empty(linear_part))) {
      edge_geometries[[edge]] <- sf::st_union(linear_part)[[1L]]
      valid_edge_geometry[[edge]] <- TRUE
    }
  }
}
if (!all(valid_edge_geometry)) {
  warning(
    "No line geometry was found for ", sum(!valid_edge_geometry),
    " neighboring county pair(s); those edges cannot be drawn."
  )
}
edge_sf <- sf::st_sf(
  edge_index = which(valid_edge_geometry),
  boundary = neighbor_name[valid_edge_geometry],
  geometry = sf::st_sfc(
    edge_geometries[valid_edge_geometry], crs = sf::st_crs(ca.poly)
  )
)

scale_selected_widths <- function(probabilities, minimum = 0.5, maximum = 1) {
  if (!length(probabilities)) return(numeric())
  probability_range <- range(probabilities, finite = TRUE)
  if (!is.finite(diff(probability_range)) || diff(probability_range) == 0) {
    return(rep(mean(c(minimum, maximum)), length(probabilities)))
  }
  minimum + (probabilities - probability_range[[1L]]) /
    diff(probability_range) * (maximum - minimum)
}

selected_edge_sf <- function(selected_names, probabilities) {
  selected_indices <- match(selected_names, neighbor_name)
  selected_indices <- selected_indices[
    !is.na(selected_indices) & valid_edge_geometry[selected_indices]
  ]
  if (!length(selected_indices)) return(edge_sf[FALSE, , drop = FALSE])
  selected <- edge_sf[
    match(selected_indices, edge_sf$edge_index), , drop = FALSE
  ]
  selected$posterior_probability <- probabilities[selected_indices]
  selected$line_width <- scale_selected_widths(
    selected$posterior_probability, 0.5, 1
  )
  selected
}

paper_map_theme <- function(title_size = 13) {
  ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = title_size, face = "bold", hjust = 0.5
      ),
      legend.title = ggplot2::element_blank()
    )
}

build_random_effect_boundary_plot <- function(
  disease, random_effect_mean, probabilities, selected_names
) {
  map_data <- ca.poly
  map_data$posterior_mean_random_effect <- random_effect_mean
  selected <- selected_edge_sf(selected_names, probabilities)
  plot <- ggplot2::ggplot(map_data) +
    ggplot2::geom_sf(
      ggplot2::aes(fill = posterior_mean_random_effect),
      color = "black", alpha = 0.6, linewidth = 0.25
    ) +
    ggplot2::scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(9, "YlOrRd")
    ) +
    ggplot2::ggtitle(paste0(disease, " (", length(selected_names), ")")) +
    paper_map_theme(13)
  if (nrow(selected)) {
    plot <- plot +
      ggplot2::geom_sf(
        data = selected,
        ggplot2::aes(linewidth = line_width),
        color = "blue", inherit.aes = FALSE
      ) +
      ggplot2::scale_linewidth_identity()
  }
  plot
}

# RESULT: Main Figure 5 / Supplementary Figures S13 and S22 -------------------
# Disease-specific FDR boundary maps: Adjacency = Main Figure 5,
# Mean--Adjacency = S13, and Mean = S22.
pooled_disease_map_plots <- lapply(seq_len(n_diseases), function(disease) {
  build_random_effect_boundary_plot(
    disease_names[[disease]],
    pooled_phi_mean[, disease],
    boundary_probability_mean[, disease],
    pooled_selected_boundary_sets[[disease_names[[disease]]]]
  )
})
pooled_disease_map_page <- ggpubr::ggarrange(
  plotlist = pooled_disease_map_plots, nrow = 1L, ncol = n_diseases
)
pooled_map_pdf <- file.path(
  figures_dir, paste0("pooled_fdr_boundary_maps_", cvrts, ".pdf")
)
ggplot2::ggsave(
  pooled_map_pdf, pooled_disease_map_page, width = 16, height = 5.2
)

# Chain-specific maps are retained as diagnostics. Their county colors use each
# chain's posterior mean random effects; final inference uses the pooled map.
chain_map_pdf <- file.path(
  figures_dir, paste0("chain_specific_fdr_boundary_maps_", cvrts, ".pdf")
)
grDevices::pdf(chain_map_pdf, width = 16, height = 5.5, onefile = TRUE)
tryCatch(
  {
    for (chain in seq_len(n_chains)) {
      chain_plots <- lapply(seq_len(n_diseases), function(disease) {
        build_random_effect_boundary_plot(
          disease_names[[disease]],
          chain_phi_means[, disease, chain],
          boundary_probabilities[, disease, chain],
          selected_boundary_sets[[disease_names[[disease]]]][[
            chain_names[[chain]]
          ]]
        )
      })
      page <- ggpubr::ggarrange(
        plotlist = chain_plots, nrow = 1L, ncol = n_diseases
      )
      page <- ggpubr::annotate_figure(
        page,
        top = ggpubr::text_grob(
          paste0(chain_names[[chain]], " posterior boundary selection"),
          face = "bold", size = 13
        )
      )
      print(page)
    }
  },
  finally = grDevices::dev.off()
)

build_pair_boundary_plot <- function(
  pair_name, probabilities, selected_names, edge_color = "red"
) {
  selected <- selected_edge_sf(selected_names, probabilities)
  plot <- ggplot2::ggplot(ca.poly) +
    ggplot2::geom_sf(
      fill = "white", color = "black", alpha = 0.6, linewidth = 0.25
    ) +
    ggplot2::ggtitle(paste0(pair_name, " (", length(selected_names), ")")) +
    paper_map_theme(12)
  if (nrow(selected)) {
    plot <- plot +
      ggplot2::geom_sf(
        data = selected,
        ggplot2::aes(linewidth = line_width),
        color = edge_color, inherit.aes = FALSE
      ) +
      ggplot2::scale_linewidth_identity()
  }
  plot
}

# RESULT: Main Figure 6 / Supplementary Figures S14 and S23 -------------------
# Shared FDR boundary maps: Adjacency = Main Figure 6,
# Mean--Adjacency = S14, and Mean = S23.
shared_map_plots <- lapply(seq_len(n_pairs), function(pair_index) {
  build_pair_boundary_plot(
    pair_names[[pair_index]],
    shared_boundary_probabilities[, pair_index],
    shared_boundary_results$selected_sets[[pair_names[[pair_index]]]],
    edge_color = "red"
  )
})
shared_map_page <- ggpubr::ggarrange(
  plotlist = shared_map_plots, nrow = 2L, ncol = 3L
)
shared_map_pdf <- file.path(
  figures_dir, paste0("shared_fdr_boundary_maps_", cvrts, ".pdf")
)
ggplot2::ggsave(shared_map_pdf, shared_map_page, width = 12, height = 9)

# RESULT: Main Figure 7 / Supplementary Figures S15 and S24 -------------------
# Mutual cross-difference maps: Adjacency = Main Figure 7,
# Mean--Adjacency = S15, and Mean = S24.
mutual_map_plots <- lapply(seq_len(n_pairs), function(pair_index) {
  build_pair_boundary_plot(
    pair_names[[pair_index]],
    mutual_boundary_probabilities[, pair_index],
    mutual_boundary_results$selected_sets[[pair_names[[pair_index]]]],
    edge_color = "red"
  )
})
mutual_map_page <- ggpubr::ggarrange(
  plotlist = mutual_map_plots, nrow = 2L, ncol = 3L
)
mutual_map_pdf <- file.path(
  figures_dir, paste0("mutual_fdr_boundary_maps_", cvrts, ".pdf")
)
ggplot2::ggsave(mutual_map_pdf, mutual_map_page, width = 12, height = 9)

build_non_adjacency_plot <- function(disease, probabilities) {
  drawn_indices <- which(probabilities > 0 & valid_edge_geometry)
  drawn <- edge_sf[
    match(drawn_indices, edge_sf$edge_index), , drop = FALSE
  ]
  drawn$non_adjacency_probability <- probabilities[drawn_indices]
  drawn$line_width <- 2.5 * drawn$non_adjacency_probability
  plot <- ggplot2::ggplot(ca.poly) +
    ggplot2::geom_sf(
      fill = "white", color = "black", alpha = 0.6, linewidth = 0.25
    ) +
    ggplot2::ggtitle(paste(disease, "non-adjacencies")) +
    paper_map_theme(13)
  if (nrow(drawn)) {
    plot <- plot +
      ggplot2::geom_sf(
        data = drawn,
        ggplot2::aes(linewidth = line_width),
        color = "blue", inherit.aes = FALSE
      ) +
      ggplot2::scale_linewidth_identity()
  }
  plot
}

# RESULT: Main Figure 8 / Supplementary Figure S16 ----------------------------
# Learned non-adjacency maps: Adjacency = Main Figure 8 and
# Mean--Adjacency = S16. Mean uses fixed adjacency and has no such figure.
non_adjacency_map_plots <- lapply(seq_len(n_diseases), function(disease) {
  build_non_adjacency_plot(
    disease_names[[disease]], non_adjacency_probabilities[, disease]
  )
})
non_adjacency_map_page <- ggpubr::ggarrange(
  plotlist = non_adjacency_map_plots, nrow = 1L, ncol = n_diseases
)
non_adjacency_map_pdf <- file.path(
  figures_dir, paste0("pooled_non_adjacency_maps_", cvrts, ".pdf")
)
ggplot2::ggsave(
  non_adjacency_map_pdf, non_adjacency_map_page, width = 16, height = 5.2
)

# Persist tables and binary objects --------------------------------------------

boundary_probabilities_long <- as.data.frame.table(
  boundary_probabilities,
  responseName = "boundary_probability",
  stringsAsFactors = FALSE
)
shared_boundary_probabilities_long <- as.data.frame.table(
  shared_boundary_probabilities,
  responseName = "shared_boundary_probability",
  stringsAsFactors = FALSE
)
mutual_boundary_probabilities_long <- as.data.frame.table(
  mutual_boundary_probabilities,
  responseName = "mutual_boundary_probability",
  stringsAsFactors = FALSE
)
non_adjacency_probabilities_long <- as.data.frame.table(
  non_adjacency_probabilities,
  responseName = "non_adjacency_probability",
  stringsAsFactors = FALSE
)
county_labels <- if ("NAME" %in% names(ca.poly)) {
  as.character(ca.poly$NAME)
} else {
  as.character(seq_len(n_counties))
}
pooled_random_effect_means <- data.frame(
  county = rep(county_labels, times = n_diseases),
  disease = rep(disease_names, each = n_counties),
  posterior_mean_random_effect = as.vector(pooled_phi_mean),
  row.names = NULL
)

selected_boundaries_to_data_frame <- function(boundary_sets, pooled = FALSE) {
  rows <- list()
  index <- 1L
  for (disease in names(boundary_sets)) {
    disease_sets <- if (pooled) {
      list(pooled_mean = boundary_sets[[disease]])
    } else {
      boundary_sets[[disease]]
    }
    for (chain in names(disease_sets)) {
      boundaries <- disease_sets[[chain]]
      if (length(boundaries)) {
        rows[[index]] <- data.frame(
          disease = disease, chain = chain, boundary = boundaries,
          row.names = NULL
        )
        index <- index + 1L
      }
    }
  }
  if (!length(rows)) {
    return(data.frame(
      disease = character(), chain = character(), boundary = character()
    ))
  }
  do.call(rbind, rows)
}

selected_pair_boundaries_to_data_frame <- function(results) {
  rows <- lapply(names(results$selected_sets), function(pair_name) {
    boundaries <- results$selected_sets[[pair_name]]
    if (!length(boundaries)) return(NULL)
    data.frame(
      boundary_type = results$count_summary$boundary_type[[1L]],
      cancer_pair = pair_name,
      boundary = boundaries,
      row.names = NULL
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(
      boundary_type = character(), cancer_pair = character(),
      boundary = character()
    ))
  }
  do.call(rbind, rows)
}

# RESULT SOURCES: manuscript numerical summaries and diagnostics --------------
# These CSV files are the machine-readable sources for reported boundary counts,
# posterior summaries, convergence diagnostics, map colors, and FDR selections.
diagnostic_tables <- list(
  chain_metadata = chain_metadata,
  state_consistency = state_consistency_summary,
  cold_proposal_acceptance = acceptance_summary,
  diagnostic_overview = diagnostic_overview,
  posterior_summary_all_groups = posterior_summary_all_groups,
  posterior_agreement_all_groups = posterior_agreement_all_groups,
  posterior_agreement_overview = posterior_agreement_overview,
  adjacency_cardinality_by_chain = cardinality_summary_by_chain,
  boundary_counts_by_chain = boundary_count_summary,
  pooled_boundary_counts = pooled_boundary_count_summary,
  boundary_count_ranges = boundary_count_ranges,
  boundary_probability_range_summary = boundary_probability_range_summary,
  boundary_probability_pairwise = boundary_probability_pairwise,
  boundary_probability_pairwise_summary =
    boundary_probability_pairwise_summary,
  boundary_pairwise_jaccard = boundary_pairwise_jaccard,
  boundary_pairwise_jaccard_summary = boundary_pairwise_jaccard_summary,
  boundary_agreement_by_chain = boundary_agreement,
  boundary_selection_incidence = boundary_selection_incidence,
  boundary_selection_frequency = boundary_selection_frequency_summary,
  boundary_cutoff_sensitivity = boundary_cutoff_sensitivity,
  edge_boundary_probabilities_by_chain = boundary_probabilities_long,
  pooled_random_effect_means = pooled_random_effect_means,
  pooled_shared_boundary_probabilities = shared_boundary_probabilities_long,
  pooled_mutual_boundary_probabilities = mutual_boundary_probabilities_long,
  pooled_non_adjacency_probabilities = non_adjacency_probabilities_long,
  pooled_shared_boundary_counts = shared_boundary_results$count_summary,
  pooled_mutual_boundary_counts = mutual_boundary_results$count_summary,
  pooled_shared_selected_boundaries = selected_pair_boundaries_to_data_frame(
    shared_boundary_results
  ),
  pooled_mutual_selected_boundaries = selected_pair_boundaries_to_data_frame(
    mutual_boundary_results
  ),
  selected_boundaries_by_chain = selected_boundaries_to_data_frame(
    selected_boundary_sets
  ),
  selected_boundaries_pooled = selected_boundaries_to_data_frame(
    pooled_selected_boundary_sets, pooled = TRUE
  ),
  fdr_curves = fdr_curves,
  pair_fdr_curves = pair_fdr_curves
)
if (is_tempering_workflow) {
  diagnostic_tables$replica_temperature_traffic <- replica_traffic_summary
} else {
  # Main-driver-compatible alias. Writing into the main diagnostics directory
  # intentionally replaces the earlier built-in acceptance table.
  diagnostic_tables$acceptance_rates <- acceptance_summary
}

diagnostic_csv_manifest <- do.call(
  rbind,
  lapply(names(diagnostic_tables), function(table_name) {
    table_value <- diagnostic_tables[[table_name]]
    csv_file <- file.path(
      diagnostics_dir, paste0(table_name, "_", cvrts, ".csv")
    )
    utils::write.csv(
      table_value, csv_file, row.names = FALSE, na = ""
    )
    data.frame(
      diagnostic = table_name,
      rows = nrow(table_value),
      columns = ncol(table_value),
      file = normalizePath(csv_file, winslash = "/", mustWork = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
utils::write.csv(
  diagnostic_csv_manifest,
  file.path(
    diagnostics_dir, paste0("diagnostics_manifest_", cvrts, ".csv")
  ),
  row.names = FALSE
)

diagnostics <- list(
  generated_at = Sys.time(),
  workflow = postprocess_workflow,
  model_specification = cvrts,
  input_dir = normalizePath(input_dir, winslash = "/"),
  diagnostic_backend = paste0(
    "posterior ", as.character(utils::packageVersion("posterior"))
  ),
  diagnostic_csv_manifest = diagnostic_csv_manifest,
  chain_metadata = chain_metadata,
  state_consistency_summary = state_consistency_summary,
  acceptance_summary = acceptance_summary,
  replica_traffic_summary = replica_traffic_summary,
  diagnostic_overview = diagnostic_overview,
  diagnostic_results = diagnostic_results,
  posterior_summary_all_groups = posterior_summary_all_groups,
  posterior_agreement_all_groups = posterior_agreement_all_groups,
  posterior_agreement_overview = posterior_agreement_overview
)
boundary_results <- list(
  model_specification = cvrts,
  fdr_alpha = fdr_alpha,
  boundary_probabilities = boundary_probabilities,
  pooled_boundary_probabilities = boundary_probability_mean,
  pooled_random_effect_means = pooled_phi_mean,
  chain_random_effect_means = chain_phi_means,
  pooled_W_means = pooled_W_mean,
  pooled_non_adjacency_probabilities = non_adjacency_probabilities,
  pooled_shared_boundary_probabilities = shared_boundary_probabilities,
  pooled_mutual_boundary_probabilities = mutual_boundary_probabilities,
  boundary_probability_range = boundary_probability_range,
  boundary_probability_mean_pairwise_absolute_difference =
    boundary_probability_mean_pairwise_absolute_difference,
  boundary_count_summary = boundary_count_summary,
  pooled_boundary_count_summary = pooled_boundary_count_summary,
  boundary_count_ranges = boundary_count_ranges,
  selected_boundary_sets = selected_boundary_sets,
  pooled_selected_boundary_sets = pooled_selected_boundary_sets,
  pooled_shared_boundary_results = shared_boundary_results,
  pooled_mutual_boundary_results = mutual_boundary_results,
  boundary_pairwise_jaccard = boundary_pairwise_jaccard,
  boundary_pairwise_jaccard_summary = boundary_pairwise_jaccard_summary,
  boundary_probability_pairwise = boundary_probability_pairwise,
  boundary_probability_pairwise_summary =
    boundary_probability_pairwise_summary,
  boundary_agreement = boundary_agreement,
  boundary_selection_incidence = boundary_selection_incidence,
  boundary_selection_frequency_summary =
    boundary_selection_frequency_summary,
  boundary_cutoff_sensitivity = boundary_cutoff_sensitivity,
  fdr_curves = fdr_curves,
  pair_fdr_curves = pair_fdr_curves,
  supplement_table_files = supplement_table_files,
  figure_files = c(
    pooled_fdr_curves = pooled_fdr_file,
    chain_fdr_curves = chain_fdr_file,
    pooled_boundary_maps = pooled_map_pdf,
    chain_boundary_maps = chain_map_pdf,
    shared_boundary_maps = shared_map_pdf,
    mutual_boundary_maps = mutual_map_pdf,
    non_adjacency_maps = non_adjacency_map_pdf,
    supplement_figure_files
  )
)

# A single manifest provides the one-to-one inventory of every generated table
# and figure for this specification. It is also stored in both summary RDS
# objects so downstream manuscript assembly can validate its inputs.
artifact_path <- function(path) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  if (!file.exists(path)) return(path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
figure_artifacts <- boundary_results$figure_files
table_artifacts <- supplement_table_files
empty_artifact_rows <- function() {
  data.frame(
    artifact_type = character(),
    artifact = character(),
    model_specification = character(),
    file = character(),
    stringsAsFactors = FALSE
  )
}
supplement_artifact_rows <- if (length(table_artifacts)) {
  data.frame(
    artifact_type = "supplement_table",
    artifact = names(table_artifacts),
    model_specification = cvrts,
    file = vapply(table_artifacts, artifact_path, character(1L)),
    stringsAsFactors = FALSE
  )
} else {
  empty_artifact_rows()
}
figure_artifact_rows <- if (length(figure_artifacts)) {
  data.frame(
    artifact_type = "figure",
    artifact = names(figure_artifacts),
    model_specification = cvrts,
    file = vapply(figure_artifacts, artifact_path, character(1L)),
    stringsAsFactors = FALSE
  )
} else {
  empty_artifact_rows()
}
artifact_manifest <- rbind(
  data.frame(
    artifact_type = "diagnostic_table",
    artifact = diagnostic_csv_manifest$diagnostic,
    model_specification = cvrts,
    file = diagnostic_csv_manifest$file,
    stringsAsFactors = FALSE
  ),
  supplement_artifact_rows,
  figure_artifact_rows
)
artifact_manifest$status <- ifelse(
  is.na(artifact_manifest$file) | !nzchar(artifact_manifest$file),
  "not_applicable",
  ifelse(file.exists(artifact_manifest$file), "available", "missing")
)
artifact_manifest_file <- file.path(
  output_dir, paste0("artifact_manifest_", cvrts, ".csv")
)
utils::write.csv(
  artifact_manifest,
  artifact_manifest_file,
  row.names = FALSE,
  na = ""
)
diagnostics$artifact_manifest <- artifact_manifest
diagnostics$artifact_manifest_file <- artifact_manifest_file
boundary_results$artifact_manifest <- artifact_manifest
boundary_results$artifact_manifest_file <- artifact_manifest_file
diagnostics_rds_file <- if (is_tempering_workflow) {
  file.path(
    output_dir,
    paste0(output_prefix, "_convergence_diagnostics.rds")
  )
} else {
  file.path(output_dir, "multiple_chains_diagnostics.rds")
}
boundary_rds_file <- if (is_tempering_workflow) {
  file.path(
    output_dir,
    paste0(output_prefix, "_fdr_boundary_results.rds")
  )
} else {
  file.path(output_dir, "multiple_chains_boundary_summaries.rds")
}
saveRDS(diagnostics, diagnostics_rds_file)
saveRDS(boundary_results, boundary_rds_file)
report_file <- file.path(output_dir, "postprocessing_summary.txt")
report <- capture.output({
  cat(
    "", toupper(workflow_label),
    " POST-PROCESSING\n", sep = ""
  )
    cat("Generated:", format(Sys.time()), "\n")
    cat("Input:", normalizePath(input_dir, winslash = "/"), "\n")
    cat("SEER specification:", cvrts, "\n")
    cat("Chains:", n_chains, "\n")
  cat("Retained draws per chain:", retained_draws, "\n\n")
  cat("Rhat, ESS, and MCSE overview\n")
  print(diagnostic_overview, row.names = FALSE)
  cat("\nCold proposal acceptance\n")
  print(acceptance_summary, row.names = FALSE)
  if (is_tempering_workflow) {
    cat("\nTemperature traffic\n")
    print(replica_traffic_summary, row.names = FALSE)
  }
  cat("\nChain-specific FDR-selected boundary counts\n")
  print(boundary_count_summary, row.names = FALSE)
  cat("\nPooled FDR-selected boundary counts\n")
  print(pooled_boundary_count_summary, row.names = FALSE)
  cat("\nPooled shared FDR-selected boundary counts\n")
  print(shared_boundary_results$count_summary, row.names = FALSE)
  cat("\nPooled mutual FDR-selected boundary counts\n")
  print(mutual_boundary_results$count_summary, row.names = FALSE)
  cat("\nPairwise selected-boundary Jaccard agreement\n")
  print(boundary_pairwise_jaccard_summary, row.names = FALSE)
  cat("\nEdge-probability agreement\n")
  print(boundary_probability_pairwise_summary, row.names = FALSE)
})
writeLines(report, report_file)

cat("\nPost-processing complete:", format(Sys.time()), "\n")
cat("Summary:", normalizePath(report_file, winslash = "/"), "\n")
cat("Diagnostics:", normalizePath(diagnostics_dir, winslash = "/"), "\n")
cat("Figures:", normalizePath(figures_dir, winslash = "/"), "\n")
