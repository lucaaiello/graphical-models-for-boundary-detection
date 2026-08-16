rm(list = ls())

# SEER posterior predictive checks from adaptive-tempering chains
# =============================================================================
# This script does not run MCMC. By default, it uses all 150,000 retained cold
# posterior draws (25,000 from each of six chains) and generates one replicated
# count dataset per draw for each of the Adjacency, Mean--Adjacency, and Mean
# specifications. No single-chain manuscript workspace is read or accepted as
# a fallback. Set SEER_MULTICHAIN_PPC_DRAWS to override the default.
#
# Run from the repository root with:
#   source("data analysis/posterior_predictive_checks_multichain.R")
#
# All generated artifacts are written below data analysis/ppc/.

options(stringsAsFactors = FALSE)

input_root <- Sys.getenv(
  "SEER_MULTICHAIN_INPUT_ROOT",
  unset = file.path(
    "data analysis", "multiple_chains_tempering_output"
  )
)
run_mode <- Sys.getenv("SEER_MULTICHAIN_RUN_MODE", unset = "hotter_final")
output_dir <- Sys.getenv(
  "SEER_MULTICHAIN_PPC_OUTPUT",
  unset = file.path("data analysis", "ppc")
)
figures_dir <- file.path(output_dir, "figures")
tables_dir <- file.path(output_dir, "tables")
diagnostics_dir <- file.path(output_dir, "diagnostics")
for (directory in c(output_dir, figures_dir, tables_dir, diagnostics_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

required_packages <- c("Matrix", "matrixStats", "ggplot2", "scales")
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
  adj = list(label = "Adjacency", seed_offset = 1000L),
  meanadj = list(label = "Mean--Adjacency", seed_offset = 2000L),
  mean = list(label = "Mean", seed_offset = 3000L)
)
disease_names <- c("Lung", "Esophageal", "Larynx", "Colorectal")
n_chains <- 6L
ppc_seed <- as.integer(Sys.getenv("SEER_MULTICHAIN_PPC_SEED", unset = "84731"))
ppc_draws_per_specification <- as.integer(Sys.getenv(
  "SEER_MULTICHAIN_PPC_DRAWS", unset = "150000"
))
unusual_threshold <- as.numeric(Sys.getenv(
  "SEER_MULTICHAIN_PPC_MIDTAIL_THRESHOLD", unset = "0.02"
))
if (
  is.na(ppc_seed) || is.na(ppc_draws_per_specification) ||
    ppc_draws_per_specification < n_chains ||
    !is.finite(unusual_threshold) || unusual_threshold <= 0 ||
    unusual_threshold >= 1
) {
  stop("Invalid PPC seed, draw count, or mid-tail threshold.", call. = FALSE)
}

stopf <- function(...) stop(sprintf(...), call. = FALSE)
title_case_county <- function(x) tools::toTitleCase(as.character(x))
md5_file <- function(path) unname(tools::md5sum(path))

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_text <- function(x, path) {
  writeLines(x, path, useBytes = TRUE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# Reconstruct data and graph ordering, stopping before sampler compilation.
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
  adjacency = unname(as.matrix(get("Adj", setup))),
  county_names = as.character(get("county.ID", setup)),
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
if (
  n_counties != 58L || n_diseases != 4L || nrow(reference$edges) != 139L ||
    !all(dim(reference$adjacency) == c(n_counties, n_counties)) ||
    length(reference$Y) != n_counties * n_diseases ||
    length(reference$E) != n_counties * n_diseases ||
    any(reference$E <= 0) ||
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

reorder_draws_to_geographic <- function(draws) {
  inverse_permutation <- order(reference$final_perm)
  do.call(cbind, lapply(seq_len(n_diseases), function(disease) {
    indices <- ((disease - 1L) * n_counties + 1L):(disease * n_counties)
    draws[, indices, drop = FALSE][, inverse_permutation, drop = FALSE]
  }))
}

chain_paths <- function(specification) {
  file.path(
    input_root, specification, run_mode, "ladders",
    sprintf("ladder_%d_cold_chain.rds", seq_len(n_chains))
  )
}

all_chain_files <- unlist(lapply(names(specifications), chain_paths),
                          use.names = FALSE)
missing_chain_files <- all_chain_files[!file.exists(all_chain_files)]
if (length(missing_chain_files)) {
  stop(
    "The multiple-chain PPC inputs are incomplete. Complete the ",
    "hotter_final run for all three specifications.\n",
    paste(" -", missing_chain_files, collapse = "\n"), call. = FALSE
  )
}

validate_chain <- function(fit, specification, chain_file) {
  required <- c("beta", "phi", "settings", "seed", "runtime")
  missing <- setdiff(required, names(fit))
  if (length(missing)) {
    stopf("%s is missing: %s", chain_file, paste(missing, collapse = ", "))
  }
  saved_specification <- fit$settings$model_specification
  if (is.null(saved_specification)) saved_specification <- "adj"
  X <- design_matrices[[specification]]
  if (!identical(saved_specification, specification)) {
    stopf("%s contains specification '%s', not '%s'.", chain_file,
          saved_specification, specification)
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
  if (
    !is.matrix(fit$beta) || !is.matrix(fit$phi) ||
      nrow(fit$beta) < 2L || nrow(fit$phi) != nrow(fit$beta) ||
      ncol(fit$beta) != ncol(X) ||
      ncol(fit$phi) != n_counties * n_diseases ||
      any(!is.finite(fit$beta)) || any(!is.finite(fit$phi))
  ) {
    stopf("%s has incompatible or non-finite beta/phi draws.", chain_file)
  }
  invisible(TRUE)
}

allocate_draws <- function(total, chains) {
  allocation <- rep(total %/% chains, chains)
  allocation[seq_len(total %% chains)] <-
    allocation[seq_len(total %% chains)] + 1L
  allocation
}

select_evenly_spaced <- function(number_available, number_requested) {
  if (number_requested > number_available) {
    stopf("Requested %d draws from a chain containing only %d.",
          number_requested, number_available)
  }
  if (number_requested == 1L) return(as.integer(ceiling(number_available / 2)))
  indices <- unique(as.integer(round(seq(
    1, number_available, length.out = number_requested
  ))))
  if (length(indices) != number_requested) {
    stop("Could not construct the requested deterministic draw subset.",
         call. = FALSE)
  }
  indices
}

moran_rows <- function(values, adjacency) {
  centered <- values - rowMeans(values)
  denominator <- rowSums(centered * centered)
  numerator <- rowSums((centered %*% adjacency) * centered)
  result <- ncol(values) / sum(adjacency) * numerator / denominator
  result[denominator == 0] <- NA_real_
  result
}

row_standard_deviation <- function(x) {
  matrixStats::rowSds(x)
}

edge_contrast_rows <- function(counts, expected, edges) {
  log_sir <- log(sweep(counts + 0.5, 2, expected, FUN = "/"))
  edge_differences <- abs(
    log_sir[, edges[, 1], drop = FALSE] -
      log_sir[, edges[, 2], drop = FALSE]
  )
  matrixStats::rowQuantiles(edge_differences, probs = 0.95, drop = TRUE)
}

summarise_predictive_statistic <- function(
  observed, replicated, specification, cancer, diagnostic, diagnostic_label
) {
  replicated <- replicated[is.finite(replicated)]
  if (!length(replicated) || !is.finite(observed)) {
    stopf("No finite PPC values for %s, %s, %s.", specification, cancer,
          diagnostic)
  }
  interval <- stats::quantile(
    replicated, c(0.025, 0.5, 0.975), names = FALSE, type = 8
  )
  lower <- mean(replicated <= observed)
  upper <- mean(replicated >= observed)
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
    lower_tail_probability = lower,
    upper_tail_probability = upper,
    two_sided_ppp = min(1, 2 * min(lower, upper)),
    observed_within_95 = observed >= interval[[1]] && observed <= interval[[3]],
    stringsAsFactors = FALSE
  )
}

pointwise_predictive_summary <- function(
  observed, expected, replicated, counties, cancer, specification
) {
  predictive_mean <- colMeans(replicated)
  predictive_sd <- matrixStats::colSds(replicated)
  predictive_interval <- matrixStats::colQuantiles(
    replicated, probs = c(0.025, 0.5, 0.975), drop = FALSE
  )
  lower_mid_tail <- vapply(seq_along(observed), function(index) {
    mean(replicated[, index] < observed[[index]]) +
      0.5 * mean(replicated[, index] == observed[[index]])
  }, numeric(1L))
  upper_mid_tail <- vapply(seq_along(observed), function(index) {
    mean(replicated[, index] > observed[[index]]) +
      0.5 * mean(replicated[, index] == observed[[index]])
  }, numeric(1L))
  two_sided <- pmin(1, 2 * pmin(lower_mid_tail, upper_mid_tail))
  data.frame(
    specification = specification,
    county = title_case_county(counties),
    cancer = cancer,
    observed_count = observed,
    expected_count = expected,
    predictive_mean = predictive_mean,
    predictive_median = predictive_interval[, 2],
    predictive_lower_95 = predictive_interval[, 1],
    predictive_upper_95 = predictive_interval[, 3],
    pointwise_95_coverage = observed >= predictive_interval[, 1] &
      observed <= predictive_interval[, 3],
    standardized_predictive_residual = ifelse(
      predictive_sd > 0, (observed - predictive_mean) / predictive_sd, NA_real_
    ),
    lower_mid_tail_probability = lower_mid_tail,
    upper_mid_tail_probability = upper_mid_tail,
    two_sided_mid_tail_probability = two_sided,
    low_predictive_support = two_sided < unusual_threshold,
    descriptive_threshold = unusual_threshold,
    stringsAsFactors = FALSE
  )
}

generate_replicates <- function(specification, specification_index) {
  files <- chain_paths(specification)
  allocation <- allocate_draws(ppc_draws_per_specification, n_chains)
  X <- design_matrices[[specification]]
  replicated_parts <- vector("list", n_chains)
  fitted_risk_sum <- numeric(n_counties * n_diseases)
  selected_rows <- vector("list", n_chains)
  metadata <- vector("list", n_chains)

  for (chain in seq_len(n_chains)) {
    cat(sprintf(
      "  %s: chain %d/%d (%d replicated draws)\n",
      specifications[[specification]]$label, chain, n_chains,
      allocation[[chain]]
    ))
    fit <- readRDS(files[[chain]])
    validate_chain(fit, specification, files[[chain]])
    indices <- select_evenly_spaced(nrow(fit$beta), allocation[[chain]])
    linear_predictor <-
      fit$beta[indices, , drop = FALSE] %*% t(X) +
      fit$phi[indices, , drop = FALSE]
    relative_risk <- exp(linear_predictor)
    poisson_mean <- sweep(relative_risk, 2, reference$E, FUN = "*")
    if (any(!is.finite(poisson_mean)) || any(poisson_mean < 0)) {
      stopf("%s produced invalid posterior predictive means.", files[[chain]])
    }
    set.seed(ppc_seed + specifications[[specification]]$seed_offset + chain)
    replicated_model_order <- matrix(
      stats::rpois(length(poisson_mean), lambda = as.vector(poisson_mean)),
      nrow = length(indices), ncol = n_counties * n_diseases
    )
    replicated_parts[[chain]] <-
      reorder_draws_to_geographic(replicated_model_order)
    fitted_risk_sum <- fitted_risk_sum + colSums(relative_risk)
    selected_rows[[chain]] <- data.frame(
      specification = specification,
      chain = chain,
      source_draw = indices,
      ppc_seed = ppc_seed + specifications[[specification]]$seed_offset + chain,
      stringsAsFactors = FALSE
    )
    metadata[[chain]] <- data.frame(
      specification = specification,
      specification_label = specifications[[specification]]$label,
      chain = chain,
      seed = fit$seed,
      retained_draws_available = nrow(fit$beta),
      ppc_draws_used = length(indices),
      warmup_iterations = fit$settings$warmup_iterations,
      sampling_thin = fit$settings$sampling_thin,
      n_temperatures = fit$settings$n_temperatures,
      beta_hot = fit$settings$learned_beta_hot,
      chain_file = normalizePath(files[[chain]], winslash = "/"),
      md5 = md5_file(files[[chain]]),
      stringsAsFactors = FALSE
    )
    rm(fit, linear_predictor, relative_risk, poisson_mean,
       replicated_model_order)
    invisible(gc(FALSE))
  }

  replicated <- do.call(rbind, replicated_parts)
  if (nrow(replicated) != ppc_draws_per_specification) {
    stopf("%s produced %d rather than %d PPC draws.", specification,
          nrow(replicated), ppc_draws_per_specification)
  }
  list(
    replicated = replicated,
    fitted_risk_mean = reorder_vector_to_geographic(
      fitted_risk_sum / ppc_draws_per_specification
    ),
    selected_rows = do.call(rbind, selected_rows),
    metadata = do.call(rbind, metadata)
  )
}

analyse_replicates <- function(generated, specification) {
  replicated <- generated$replicated
  observed <- reorder_vector_to_geographic(reference$Y)
  expected <- reorder_vector_to_geographic(reference$E)
  global_rows <- list()
  pointwise_rows <- list()
  diagnostic_labels <- c(
    total = "Total disease burden",
    sir_sd = "Across-county SIR SD",
    moran_i = "Moran's I of county SIR",
    edge_log_sir_q95 = "95th percentile neighboring log-SIR contrast"
  )

  for (disease in seq_len(n_diseases)) {
    columns <- ((disease - 1L) * n_counties + 1L):(disease * n_counties)
    replicated_disease <- replicated[, columns, drop = FALSE]
    observed_disease <- observed[columns]
    expected_disease <- expected[columns]
    replicated_sir <- sweep(replicated_disease, 2, expected_disease, FUN = "/")
    observed_sir <- observed_disease / expected_disease
    replicated_statistics <- list(
      total = rowSums(replicated_disease),
      sir_sd = row_standard_deviation(replicated_sir),
      moran_i = moran_rows(replicated_sir, reference$adjacency),
      edge_log_sir_q95 = edge_contrast_rows(
        replicated_disease, expected_disease, reference$edges
      )
    )
    observed_statistics <- c(
      total = sum(observed_disease),
      sir_sd = stats::sd(observed_sir),
      moran_i = moran_rows(matrix(observed_sir, nrow = 1L),
                           reference$adjacency)[[1]],
      edge_log_sir_q95 = edge_contrast_rows(
        matrix(observed_disease, nrow = 1L), expected_disease,
        reference$edges
      )[[1]]
    )
    global_rows[[disease]] <- do.call(rbind, lapply(
      names(replicated_statistics), function(diagnostic) {
        summarise_predictive_statistic(
          observed_statistics[[diagnostic]],
          replicated_statistics[[diagnostic]],
          specifications[[specification]]$label,
          disease_names[[disease]], diagnostic,
          diagnostic_labels[[diagnostic]]
        )
      }
    ))
    pointwise_rows[[disease]] <- pointwise_predictive_summary(
      observed_disease, expected_disease, replicated_disease,
      reference$county_names, disease_names[[disease]],
      specifications[[specification]]$label
    )
  }
  list(
    global = do.call(rbind, global_rows),
    pointwise = do.call(rbind, pointwise_rows),
    fitted_risk_mean = generated$fitted_risk_mean,
    selected_rows = generated$selected_rows,
    metadata = generated$metadata
  )
}

cat("Multiple-chain SEER posterior predictive checks\n")
cat("Started:", format(Sys.time()), "\n")
cat("Input root:", normalizePath(input_root, winslash = "/"), "\n")
cat("Run mode:", run_mode, "\n")
cat("Replicated draws per specification:", ppc_draws_per_specification, "\n")
cat("Base PPC seed:", ppc_seed, "\n")
cat("Output:", normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    "\n\n")

ppc_results <- lapply(seq_along(specifications), function(index) {
  specification <- names(specifications)[[index]]
  cat("Generating", specifications[[specification]]$label, "PPC...\n")
  generated <- generate_replicates(specification, index)
  analysed <- analyse_replicates(generated, specification)
  rm(generated)
  invisible(gc(FALSE))
  analysed
})
names(ppc_results) <- names(specifications)

global_ppc <- do.call(rbind, lapply(ppc_results, `[[`, "global"))
pointwise_ppc <- do.call(rbind, lapply(ppc_results, `[[`, "pointwise"))
selected_rows <- do.call(rbind, lapply(ppc_results, `[[`, "selected_rows"))
chain_metadata <- do.call(rbind, lapply(ppc_results, `[[`, "metadata"))
row.names(global_ppc) <- row.names(pointwise_ppc) <- NULL

if (
  nrow(global_ppc) != length(specifications) * n_diseases * 4L ||
    nrow(pointwise_ppc) != length(specifications) * n_diseases * n_counties ||
    any(table(pointwise_ppc$specification, pointwise_ppc$cancer) != n_counties)
) {
  stop("The completed PPC summaries have unexpected dimensions.",
       call. = FALSE)
}

# RESULT: complete machine-readable global and observation-level PPC tables ----
global_ppc_file <- file.path(
  diagnostics_dir, "global_ppc_multichain_three_specifications.csv"
)
pointwise_ppc_file <- file.path(
  diagnostics_dir, "observation_level_ppc_multichain_three_specifications.csv"
)
write_csv(global_ppc, global_ppc_file)
write_csv(pointwise_ppc, pointwise_ppc_file)

global_exceptions <- global_ppc[!global_ppc$observed_within_95,]
if (nrow(global_exceptions)) {
  global_exceptions$direction <- ifelse(
    global_exceptions$observed > global_exceptions$predictive_upper_95,
    "Observed above the 97.5% predictive quantile",
    "Observed below the 2.5% predictive quantile"
  )
}
global_exceptions_file <- file.path(
  tables_dir, "global_ppc_exceptions_multichain_three_specifications.csv"
)
write_csv(global_exceptions, global_exceptions_file)

# RESULT: observation-level PPC summary table ---------------------------------
specification_labels <- vapply(specifications, `[[`, character(1L), "label")
observation_summary <- do.call(rbind, lapply(
  specification_labels, function(specification) {
    x <- pointwise_ppc[pointwise_ppc$specification == specification,]
    data.frame(
      Specification = specification,
      Observed = mean(x$observed_count),
      Predicted = mean(x$predictive_mean),
      Lower = mean(x$predictive_lower_95),
      Upper = mean(x$predictive_upper_95),
      Coverage = mean(x$pointwise_95_coverage),
      `Low support` = sum(x$low_predictive_support),
      check.names = FALSE
    )
  }
))
observation_summary_csv <- file.path(
  tables_dir,
  "observation_ppc_summary_multichain_three_specifications.csv"
)
write_csv(observation_summary, observation_summary_csv)
summary_rows <- sprintf(
  "%s & %.1f & %.1f & %.1f & %.1f & %.3f & %d \\\\",
  observation_summary$Specification, observation_summary$Observed,
  observation_summary$Predicted, observation_summary$Lower,
  observation_summary$Upper, observation_summary$Coverage,
  observation_summary[["Low support"]]
)
observation_summary_tex <- file.path(
  tables_dir,
  "observation_ppc_summary_multichain_three_specifications.tex"
)
write_text(c(
  "\\begin{table}[!ht]", "\\centering", "\\small",
  paste0(
    "\\caption{Observation-level posterior predictive summaries from ",
    "150,000 retained cold-posterior draws pooled across six adaptive-",
    "tempering runs per specification. Observed counts, predictive means, ",
    "and predictive limits are averages over the 232 county--cancer ",
    "observations; coverage and low support are evaluated pointwise.}"
  ),
  "\\label{tab:observation_ppc}", "\\begin{tabular}{lrrrrrr}",
  "\\hline",
  paste0(
    "Specification & Mean observed & Mean predicted & Mean lower & Mean upper & ",
    "Coverage & Low support \\\\"
  ),
  "\\hline", summary_rows, "\\hline", "\\end{tabular}",
  "\\end{table}"
), observation_summary_tex)

# RESULT: low-predictive-support observations table ---------------------------
unusual <- pointwise_ppc[pointwise_ppc$low_predictive_support,]
unusual <- unusual[order(
  unusual$two_sided_mid_tail_probability,
  -abs(unusual$standardized_predictive_residual)
),]
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
  check.names = FALSE
)
unusual_csv <- file.path(
  tables_dir,
  "unusual_observations_ppc_multichain_three_specifications.csv"
)
write_csv(unusual_output, unusual_csv)
if (nrow(unusual_output)) {
  unusual_rows <- sprintf(
    "%s & %s & %s & %d & %.1f & %.1f & [%g, %g] & %.2f & %.4f \\\\",
    unusual_output$Specification, unusual_output$County,
    unusual_output$Cancer, unusual_output$Observed,
    unusual_output[["Predictive mean"]],
    unusual_output[["Predictive median"]],
    unusual_output[["Predictive lower 95"]],
    unusual_output[["Predictive upper 95"]],
    unusual_output[["Standardized residual"]],
    unusual_output[["Mid-tail probability"]]
  )
} else {
  unusual_rows <- paste0(
    "\\multicolumn{9}{c}{No observation crossed the descriptive threshold.}",
    " \\\\"
  )
}
unusual_tex <- file.path(
  tables_dir,
  "unusual_observations_ppc_multichain_three_specifications.tex"
)
write_text(c(
  "\\begin{table}[!ht]", "\\centering", "\\small",
  paste0(
    "\\caption{County--cancer observations with two-sided predictive ",
    "mid-tail probability below $", unusual_threshold, "$ in the pooled ",
    "six-chain posterior predictive analysis.}"
  ),
  "\\label{tab:unusual_observation_ppc}",
  "\\begin{tabular}{lllrrrcrr}", "\\hline",
  paste0(
    "Specification & County & Cancer & Observed & Pred. mean & Pred. median & ",
    "95\\% predictive interval & Std. residual & Mid-tail $p$ \\\\"
  ),
  "\\hline", unusual_rows, "\\hline", "\\end{tabular}", "\\end{table}"
), unusual_tex)

# RESULT: global PPC figure ----------------------------------------------------
global_plot_data <- global_ppc
global_plot_data$specification <- factor(
  global_plot_data$specification, levels = rev(specification_labels)
)
global_plot_data$cancer <- factor(global_plot_data$cancer,
                                  levels = disease_names)
global_ppc_plot <- ggplot2::ggplot(
  global_plot_data, ggplot2::aes(y = specification)
) +
  ggplot2::geom_segment(ggplot2::aes(
    x = predictive_lower_95, xend = predictive_upper_95,
    yend = specification, color = specification
  ), linewidth = 1.1) +
  ggplot2::geom_point(ggplot2::aes(
    x = predictive_median, color = specification, shape = "Predictive median"
  ), size = 2.5) +
  ggplot2::geom_point(ggplot2::aes(x = observed, shape = "Observed"),
                      color = "black", size = 2.8, stroke = 1) +
  ggplot2::facet_wrap(
    ~ diagnostic_label + cancer, scales = "free_x", ncol = 4,
    labeller = ggplot2::label_wrap_gen(width = 28)
  ) +
  ggplot2::scale_shape_manual(
    name = NULL, values = c("Predictive median" = 16, "Observed" = 4)
  ) +
  ggplot2::guides(color = "none") +
  ggplot2::labs(
    x = "Diagnostic value", y = NULL,
    caption = paste(
      "Lines show 95% posterior predictive intervals.",
      "The local contrast is the 95th percentile of absolute neighboring",
      "stabilized log-SIR differences."
    )
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    legend.position = "bottom", panel.grid.minor = ggplot2::element_blank()
  )
global_ppc_pdf <- file.path(
  figures_dir, "ppc_multichain_three_specifications.pdf"
)
global_ppc_png <- sub("\\.pdf$", ".png", global_ppc_pdf)
ggplot2::ggsave(global_ppc_pdf, global_ppc_plot, width = 14, height = 10)
ggplot2::ggsave(global_ppc_png, global_ppc_plot, width = 14, height = 10,
                dpi = 300)

# RESULT: observed-versus-fitted PPC figure -----------------------------------
observed_plot_data <- pointwise_ppc
observed_plot_data$specification <- factor(
  observed_plot_data$specification, levels = specification_labels
)
observed_plot_data$cancer <- factor(observed_plot_data$cancer,
                                    levels = disease_names)
observed_plot_data$panel <- factor(
  paste(observed_plot_data$cancer, observed_plot_data$specification, sep = "\n"),
  levels = unlist(lapply(disease_names, function(cancer) {
    paste(cancer, specification_labels, sep = "\n")
  }))
)
observed_plot_data$flag_label <- ifelse(
  observed_plot_data$low_predictive_support,
  observed_plot_data$county, ""
)
observed_vs_fitted_plot <- ggplot2::ggplot(
  observed_plot_data,
  ggplot2::aes(x = observed_count, y = predictive_mean)
) +
  ggplot2::geom_abline(
    intercept = 0, slope = 1, linetype = 2, color = "grey45",
    linewidth = 0.45
  ) +
  ggplot2::geom_segment(ggplot2::aes(
    xend = observed_count, y = predictive_lower_95,
    yend = predictive_upper_95, color = low_predictive_support
  ), linewidth = 0.35, alpha = 0.7) +
  ggplot2::geom_point(ggplot2::aes(
    color = low_predictive_support, shape = low_predictive_support
  ), size = 1.35, alpha = 0.9) +
  ggplot2::geom_label(
    data = observed_plot_data[observed_plot_data$flag_label != "",],
    ggplot2::aes(label = flag_label), hjust = -0.15, vjust = 0.5,
    size = 2.6, color = "#B2182B", fill = "white", linewidth = 0.15,
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~ panel, ncol = 3, scales = "free") +
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
    labels = c(`FALSE` = "Other", `TRUE` = paste0(
      "Mid-tail p < ", unusual_threshold
    )), name = NULL
  ) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                              guide = "none") +
  ggplot2::labs(
    x = "Observed county-level count (pseudo-log scale)",
    y = "Posterior predictive mean count (pseudo-log scale)"
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold", size = 8),
    panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom"
  )
observed_vs_fitted_pdf <- file.path(
  figures_dir, "observed_vs_fitted_ppc_multichain_three_specifications.pdf"
)
observed_vs_fitted_png <- sub("\\.pdf$", ".png", observed_vs_fitted_pdf)
ggplot2::ggsave(observed_vs_fitted_pdf, observed_vs_fitted_plot,
                width = 10, height = 13)
ggplot2::ggsave(observed_vs_fitted_png, observed_vs_fitted_plot,
                width = 10, height = 13, dpi = 300)

# Reproducibility records ------------------------------------------------------
selected_rows_file <- file.path(
  diagnostics_dir, "selected_posterior_draws_ppc_multichain.csv"
)
metadata_file <- file.path(
  diagnostics_dir, "chain_metadata_ppc_multichain.csv"
)
write_csv(selected_rows, selected_rows_file)
write_csv(chain_metadata, metadata_file)

artifact_paths <- c(
  global_ppc_file, pointwise_ppc_file, global_exceptions_file,
  observation_summary_csv, observation_summary_tex, unusual_csv, unusual_tex,
  global_ppc_pdf, global_ppc_png, observed_vs_fitted_pdf,
  observed_vs_fitted_png, selected_rows_file, metadata_file
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
  output_dir, "artifact_manifest_ppc_multichain.csv"
)
write_csv(artifact_manifest, artifact_manifest_file)

saveRDS(list(
  generated_at = Sys.time(),
  run_mode = run_mode,
  ppc_seed = ppc_seed,
  ppc_draws_per_specification = ppc_draws_per_specification,
  unusual_threshold = unusual_threshold,
  global_ppc = global_ppc,
  pointwise_ppc = pointwise_ppc,
  observation_summary = observation_summary,
  unusual_observations = unusual_output,
  fitted_risk_means = lapply(ppc_results, `[[`, "fitted_risk_mean"),
  selected_draws = selected_rows,
  chain_metadata = chain_metadata,
  artifact_manifest = artifact_manifest
), file.path(output_dir, "ppc_results_multichain.rds"))

cat("\nPosterior predictive checks complete.\n")
cat("Global checks inside 95% intervals:",
    sum(global_ppc$observed_within_95), "of", nrow(global_ppc), "\n")
cat("Low-support county-cancer checks:", nrow(unusual_output), "\n")
cat("Artifacts:", normalizePath(output_dir, winslash = "/"), "\n")
