rm(list = ls())

# Adaptive parallel tempering for the SEER application
# =============================================================================
# This workflow is built from sampler_multiple_chains.cpp, not from
# the deleted experimental tempering code. The number of replicas is explicit
# and defaults to 16. During warm-up, the interior temperatures are adapted to remove local
# exchange bottlenecks. The cold endpoint beta=1 and the selected hot endpoint
# (0.05 by default and 0.01 for the hotter modes) are fixed because swap
# acceptance alone cannot determine how much heating is needed to cross
# posterior modes. The complete ladder and all local proposal scales are frozen
# and validated before any retained draw is generated.
#
# Run from the repository root, for example in PowerShell:
#   $env:MCMC_TEMPERING_SPECIFICATION = "adj"
#   $env:MCMC_TEMPERING_MODE = "smoke"
#   Rscript "data analysis/main_multiple_chains_tempering.R"
# Specifications: adj, meanadj, mean.
# Modes: smoke, quick, calibration, rhat_smoke, rhat, hotter_smoke,
# hotter_pilot, hotter_final, production.
# Results are stored under:
#   data analysis/multiple_chains_tempering_output/
#     <specification>/<mode>/

options(stringsAsFactors = FALSE)

model_specification <- tolower(Sys.getenv(
  "MCMC_TEMPERING_SPECIFICATION", unset = "adj"
))
if (!model_specification %in% c("adj", "meanadj", "mean")) {
  stop(
    "MCMC_TEMPERING_SPECIFICATION must be one of ",
    "'adj', 'meanadj', or 'mean'.", call. = FALSE
  )
}

# Reuse the audited SEER data-construction section from the model driver.
# Evaluation stops before that file's sampling workflow starts.
model_driver <- file.path(
  "data analysis", "main_multiple_chains.R"
)
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
eval(
  parse(text = model_source[2L:(setup_marker - 1L)]),
  envir = .GlobalEnv
)
rm(model_source, setup_marker)

# The shared setup defaults to Adjacency. Override that choice explicitly and
# reconstruct the design matrix so all three SEER specifications use exactly
# the common SEER data construction used by the multiple-chain workflows.
cvrts <- model_specification
if (cvrts %in% c("mean", "meanadj")) {
  X <- as.matrix(Matrix::bdiag(
    Matrix::bdiag(X1[, c(1, 2, 4, 6)], X2[, c(1, 2, 4, 6)]),
    Matrix::bdiag(X3[, c(1, 2, 4, 6)], X4[, c(1, 2, 4, 6)])
  ))
} else {
  X <- as.matrix(Matrix::bdiag(
    Matrix::bdiag(X1[, 1], X2[, 1]),
    Matrix::bdiag(X3[, 1], X4[, 1])
  ))
}

# Controls --------------------------------------------------------------------

run_mode <- Sys.getenv(
  "MCMC_TEMPERING_MODE", unset = "smoke"
)
stopifnot(run_mode %in% c(
  "smoke", "quick", "calibration", "rhat_smoke", "rhat",
  "hotter_smoke", "hotter_pilot", "hotter_final", "production"
))

n_temperatures <- as.integer(Sys.getenv(
  "MCMC_TEMPERING_N_TEMPERATURES", unset = "16"
))
beta_cold <- 1.0
initial_beta_hot <- as.numeric(Sys.getenv(
  "MCMC_TEMPERING_BETA_HOT", unset = "0.05"
))
endpoint_strategy <- tolower(Sys.getenv(
  "MCMC_TEMPERING_ENDPOINT_STRATEGY", unset = "fixed"
))
stopifnot(endpoint_strategy %in% c("fixed", "adaptive"))
target_swap_acceptance <- 0.25
minimum_acceptable_swap_acceptance <- 0.10
adaptation_exponent <- 0.60
adaptation_gain <- 1.0
# Broad numerical safeguards. Hitting either bound is a failed calibration,
# not a silently accepted final ladder.
minimum_learned_beta_hot <- 0.005
maximum_learned_beta_hot <- 0.50
stabilization_fraction <- 0.20
joint_adaptation_fraction <- 0.40
temperature_only_adaptation_fraction <- 0.30
state_consistency_tolerance <- 1e-12
reuse_completed_ladders <- identical(
  tolower(Sys.getenv("MCMC_TEMPERING_REUSE", unset = "true")),
  "true"
)
complete_exchange_cycle <- run_mode %in% c(
  "hotter_smoke", "hotter_pilot", "hotter_final"
)

if (run_mode == "hotter_smoke") {
  n_ladders <- 1L
  warmup_iterations <- 20L
  retained_draws <- 4L
  sampling_thin <- 2L
  local_sweeps_per_swap <- 1L
  sampling_sweeps_per_swap <- 1L
  retained_draws_per_swap <- 1L
  sampling_chunk_draws <- 2L
  adaptation_window_rounds <- 2L
} else if (run_mode == "rhat_smoke") {
  n_ladders <- 6L
  warmup_iterations <- 20L
  retained_draws <- 4L
  sampling_thin <- 2L
  local_sweeps_per_swap <- 1L
  sampling_sweeps_per_swap <- 1L
  retained_draws_per_swap <- 1L
  sampling_chunk_draws <- 2L
  adaptation_window_rounds <- 2L
} else if (run_mode == "smoke") {
  n_ladders <- 1L
  warmup_iterations <- 20L
  retained_draws <- 4L
  sampling_thin <- 1L
  local_sweeps_per_swap <- 2L
  sampling_sweeps_per_swap <- sampling_thin
  retained_draws_per_swap <- 1L
  sampling_chunk_draws <- retained_draws
  adaptation_window_rounds <- 2L
} else if (run_mode == "quick") {
  n_ladders <- 1L
  warmup_iterations <- 2000L
  retained_draws <- 200L
  sampling_thin <- 5L
  local_sweeps_per_swap <- 5L
  sampling_sweeps_per_swap <- sampling_thin
  retained_draws_per_swap <- 1L
  sampling_chunk_draws <- retained_draws
  adaptation_window_rounds <- 20L
} else if (run_mode == "calibration") {
  n_ladders <- 1L
  warmup_iterations <- 30000L
  retained_draws <- 5000L
  sampling_thin <- 10L
  local_sweeps_per_swap <- 10L
  sampling_sweeps_per_swap <- sampling_thin
  retained_draws_per_swap <- 1L
  sampling_chunk_draws <- retained_draws
  adaptation_window_rounds <- 50L
} else if (run_mode == "rhat") {
  n_ladders <- 6L
  warmup_iterations <- 50000L
  retained_draws <- 25000L
  sampling_thin <- 10L
  local_sweeps_per_swap <- 5L
  sampling_sweeps_per_swap <- 5L
  retained_draws_per_swap <- 1L
  sampling_chunk_draws <- 500L
  adaptation_window_rounds <- 50L
} else if (run_mode == "hotter_pilot") {
  n_ladders <- 1L
  warmup_iterations <- 30000L
  retained_draws <- 5000L
  sampling_thin <- 10L
  local_sweeps_per_swap <- 5L
  sampling_sweeps_per_swap <- 5L
  retained_draws_per_swap <- 1L
  sampling_chunk_draws <- 500L
  adaptation_window_rounds <- 50L
} else if (run_mode == "hotter_final") {
  n_ladders <- 6L
  warmup_iterations <- 50000L
  retained_draws <- 25000L
  sampling_thin <- 10L
  local_sweeps_per_swap <- 5L
  sampling_sweeps_per_swap <- 5L
  retained_draws_per_swap <- 1L
  sampling_chunk_draws <- 500L
  adaptation_window_rounds <- 50L
} else {
  n_ladders <- 6L
  warmup_iterations <- 300000L
  retained_draws <- 50000L
  sampling_thin <- 10L
  local_sweeps_per_swap <- 10L
  sampling_sweeps_per_swap <- sampling_thin
  retained_draws_per_swap <- 1L
  sampling_chunk_draws <- 1000L
  adaptation_window_rounds <- 100L
}

stopifnot(
  is.finite(n_temperatures), n_temperatures >= 3L,
  is.finite(initial_beta_hot), initial_beta_hot > 0,
  initial_beta_hot < beta_cold,
  target_swap_acceptance > 0, target_swap_acceptance < 1,
  minimum_acceptable_swap_acceptance > 0,
  minimum_acceptable_swap_acceptance < target_swap_acceptance,
  stabilization_fraction >= 0,
  joint_adaptation_fraction > 0,
  temperature_only_adaptation_fraction > 0,
  stabilization_fraction + joint_adaptation_fraction +
    temperature_only_adaptation_fraction < 1,
  warmup_iterations %% local_sweeps_per_swap == 0L,
  retained_draws %% retained_draws_per_swap == 0L,
  sampling_thin %% sampling_sweeps_per_swap == 0L,
  sampling_chunk_draws > 0L
)
if (run_mode %in% c("hotter_smoke", "hotter_pilot", "hotter_final") &&
    abs(initial_beta_hot - 0.01) > 1e-12) {
  stop(
    run_mode, " requires MCMC_TEMPERING_BETA_HOT=0.01.",
    call. = FALSE
  )
}

ladder_seeds <- c(12345L, 23456L, 34567L, 45678L, 56789L, 67890L)[
  seq_len(n_ladders)
]
initialization_regimes <- c(
  "original", "dense_adjacency", "sparse_adjacency",
  "random_prior_like", "strong_cross_dependence",
  "weak_cross_dependence"
)[seq_len(n_ladders)]

output_root <- Sys.getenv(
  "MCMC_TEMPERING_OUTPUT_ROOT",
  unset = file.path(
    "data analysis", "multiple_chains_tempering_output"
  )
)
output_dir <- file.path(output_root, model_specification, run_mode)
ladder_dir <- file.path(output_dir, "ladders")
checkpoint_dir <- file.path(output_dir, "checkpoints")
diagnostics_dir <- file.path(output_dir, "diagnostics")
dir.create(ladder_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
diagnostic_table_file <- function(stem) {
  file.path(
    diagnostics_dir, paste0(stem, "_", model_specification, ".csv")
  )
}

sampler_path <- normalizePath(
  file.path(
    "data analysis", "sampler_multiple_chains_tempering.cpp"
  ),
  mustWork = TRUE
)
sampler_md5 <- unname(tools::md5sum(sampler_path))

run_manifest <- list(
  model_specification = model_specification,
  run_mode = run_mode,
  n_ladders = n_ladders,
  n_temperatures = n_temperatures,
  endpoint_strategy = endpoint_strategy,
  initial_beta_hot = initial_beta_hot,
  warmup_iterations = warmup_iterations,
  retained_draws = retained_draws,
  sampling_thin = sampling_thin,
  local_sweeps_per_swap = local_sweeps_per_swap,
  sampling_sweeps_per_swap = sampling_sweeps_per_swap,
  complete_exchange_cycle = complete_exchange_cycle,
  sampler_md5 = sampler_md5
)
manifest_file <- file.path(output_dir, "run_manifest.rds")
if (file.exists(manifest_file)) {
  saved_manifest <- readRDS(manifest_file)
  manifest_comparison <- all.equal(
    saved_manifest, run_manifest, check.attributes = FALSE
  )
  if (!isTRUE(manifest_comparison)) {
    stop(
      "The existing output directory belongs to a different run ",
      "configuration:\n", output_dir, "\n", paste(manifest_comparison,
      collapse = "\n"), "\nUse a different mode or output root.",
      call. = FALSE
    )
  }
} else {
  existing_ladder_files <- list.files(
    ladder_dir, pattern = "^ladder_[0-9]+_cold_chain\\.rds$",
    full.names = TRUE
  )
  if (length(existing_ladder_files)) {
    existing_chain <- readRDS(existing_ladder_files[[1L]])
    existing_specification <- existing_chain$settings$model_specification
    if (is.null(existing_specification)) existing_specification <- "adj"
    existing_configuration <- list(
      model_specification = existing_specification,
      run_mode = existing_chain$settings$run_mode,
      n_temperatures = existing_chain$settings$n_temperatures,
      endpoint_strategy = existing_chain$settings$endpoint_strategy,
      initial_beta_hot = existing_chain$settings$initial_beta_hot,
      warmup_iterations = existing_chain$settings$warmup_iterations,
      retained_draws = existing_chain$settings$retained_draws,
      sampling_thin = existing_chain$settings$sampling_thin,
      sampling_sweeps_per_swap =
        existing_chain$settings$sampling_sweeps_per_swap,
      complete_exchange_cycle =
        existing_chain$settings$complete_exchange_cycle,
      sampler_md5 = existing_chain$sampler_md5
    )
    expected_configuration <- run_manifest[names(existing_configuration)]
    configuration_comparison <- all.equal(
      existing_configuration, expected_configuration, check.attributes = FALSE
    )
    if (!isTRUE(configuration_comparison)) {
      stop(
        "Existing ladder files do not match the requested run:\n",
        output_dir, "\n", paste(configuration_comparison, collapse = "\n"),
        call. = FALSE
      )
    }
  }
  saveRDS(run_manifest, manifest_file)
}

model_arguments <- list(
  y = Y, X = X, Z1 = Z1, Z2 = Z2, Z3 = Z3, E = E,
  cvrts = cvrts, q = 4L, Winc = Winc, Minc = Minc,
  alpha = 1, n_atoms = 15L
)

# Ladder adaptation -----------------------------------------------------------

initial_inverse_temperatures <- function() {
  if (run_mode %in% c("rhat_smoke", "rhat") && n_temperatures == 16L &&
      endpoint_strategy == "fixed" && abs(initial_beta_hot - 0.05) < 1e-12) {
    return(c(
      1.00000000, 0.87621658, 0.75968419, 0.65203311,
      0.54580492, 0.46603103, 0.38759378, 0.32345111,
      0.26991076, 0.22110959, 0.17807059, 0.14000361,
      0.11317276, 0.08519953, 0.06585721, 0.05000000
    ))
  }
  exp(seq(
    log(beta_cold), log(initial_beta_hot), length.out = n_temperatures
  ))
}

adapt_inverse_temperatures <- function(
    inverse_temperatures, accepted, proposed, adaptation_index) {
  if (any(proposed <= 0L)) return(inverse_temperatures)
  log_temperature <- -log(inverse_temperatures)
  gaps <- diff(log_temperature)
  rates <- accepted / proposed
  gain <- adaptation_gain / (adaptation_index + 10)^adaptation_exponent
  adaptation_center <- if (endpoint_strategy == "adaptive") {
    target_swap_acceptance
  } else {
    # With both endpoints fixed, only relative gaps can change. Centering by
    # the current mean makes the objective explicit: equalize pairwise rates.
    mean(rates)
  }
  increment <- gain * (rates - adaptation_center)
  # Clipping guards against one very noisy adaptation window.
  increment <- pmax(-0.50, pmin(0.50, increment))
  log_gaps <- log(gaps) + increment
  gaps <- exp(log_gaps)
  # Preserve ordering and impose only broad numerical bounds on the learned
  # hot endpoint. A calibration that reaches a bound is reported as failed.
  minimum_gap <- 1e-5
  gaps <- pmax(gaps, minimum_gap)
  total_span <- sum(gaps)
  if (endpoint_strategy == "fixed") {
    fixed_span <- log(beta_cold / initial_beta_hot)
    gaps <- gaps * fixed_span / total_span
  } else {
    minimum_span <- log(beta_cold / maximum_learned_beta_hot)
    maximum_span <- log(beta_cold / minimum_learned_beta_hot)
    if (total_span < minimum_span) {
      gaps <- gaps * minimum_span / total_span
    } else if (total_span > maximum_span) {
      gaps <- gaps * maximum_span / total_span
    }
  }
  exp(-c(0, cumsum(gaps)))
}

swap_pairs_for_round <- function(round_index) {
  if (round_index %% 2L == 1L) {
    seq.int(1L, n_temperatures - 1L, by = 2L)
  } else {
    seq.int(2L, n_temperatures - 1L, by = 2L)
  }
}

perform_swaps <- function(
    states, log_likelihoods, labels, inverse_temperatures,
    round_index, proposed, accepted) {
  exchange_passes <- if (complete_exchange_cycle) {
    c(round_index, round_index + 1L)
  } else {
    round_index
  }
  for (exchange_pass in exchange_passes) {
    for (lower in swap_pairs_for_round(exchange_pass)) {
      upper <- lower + 1L
      proposed[[lower]] <- proposed[[lower]] + 1L
      log_ratio <-
        (inverse_temperatures[[lower]] - inverse_temperatures[[upper]]) *
        (log_likelihoods[[upper]] - log_likelihoods[[lower]])
      if (log(runif(1L)) < min(0, log_ratio)) {
        temporary_state <- states[[lower]]
        states[[lower]] <- states[[upper]]
        states[[upper]] <- temporary_state
        temporary_log_likelihood <- log_likelihoods[[lower]]
        log_likelihoods[[lower]] <- log_likelihoods[[upper]]
        log_likelihoods[[upper]] <- temporary_log_likelihood
        temporary_label <- labels[[lower]]
        labels[[lower]] <- labels[[upper]]
        labels[[upper]] <- temporary_label
        accepted[[lower]] <- accepted[[lower]] + 1L
      }
    }
  }
  list(
    states = states,
    log_likelihoods = log_likelihoods,
    labels = labels,
    proposed = proposed,
    accepted = accepted
  )
}

update_round_trips <- function(labels, traffic_phase, round_trips) {
  # Per-replica phases: 0 = has not visited cold; 1 = visited cold and is
  # travelling toward hot; 2 = subsequently visited hot and is returning.
  # A round trip is counted only after the complete cold -> hot -> cold path.
  cold_label <- labels[[1L]]
  if (traffic_phase[[cold_label]] == 2L) {
    round_trips[[cold_label]] <- round_trips[[cold_label]] + 1L
    traffic_phase[[cold_label]] <- 1L
  } else if (traffic_phase[[cold_label]] == 0L) {
    traffic_phase[[cold_label]] <- 1L
  }

  hot_label <- labels[[n_temperatures]]
  if (traffic_phase[[hot_label]] == 1L) {
    traffic_phase[[hot_label]] <- 2L
  }
  list(traffic_phase = traffic_phase, round_trips = round_trips)
}

# Segment and cold-draw helpers -----------------------------------------------

run_segment <- function(
    state, proposal, inverse_temperature, seed, chain_id,
    iterations, retain, thin, adaptation_offset, adapt_proposals,
    initialization_regime) {
  do.call(
    MADAGAR_tempered_segment,
    c(
      model_arguments,
      list(
        runs = as.integer(retain),
        burn = as.integer(iterations - retain * thin),
        thin = as.integer(thin),
        chain_seed = as.integer(seed),
        chain_id = as.integer(chain_id),
        worker_pid = Sys.getpid(),
        initialization_regime = initialization_regime,
        progress_log_path = "",
        print_progress = FALSE,
        state_consistency_tolerance = state_consistency_tolerance,
        inverse_temperature = inverse_temperature,
        supplied_state = state,
        supplied_proposal = proposal,
        adaptation_offset = as.integer(adaptation_offset),
        adapt_proposals = adapt_proposals
      )
    )
  )
}

matrix_fields <- c(
  "beta", "phi", "theta", "u", "rho", "V", "r", "F_r",
  "eta", "W_cardinality", "A"
)

strip_cold_segment <- function(segment) {
  output <- lapply(matrix_fields, function(field) segment[[field]])
  names(output) <- matrix_fields
  output$tau <- segment$tau
  output
}

combine_cold_segments <- function(segments, W_sum, total_draws) {
  output <- lapply(matrix_fields, function(field) {
    do.call(rbind, lapply(segments, `[[`, field))
  })
  names(output) <- matrix_fields
  output$tau <- unlist(lapply(segments, `[[`, "tau"), use.names = FALSE)
  output$W_mean <- lapply(W_sum, `/`, total_draws)
  output
}

accumulate_acceptance <- function(total, current) {
  count_names <- c(
    "beta_proposed", "beta_accepted", "theta_proposed", "theta_accepted",
    "r_proposed", "r_accepted", "V_proposed", "V_accepted",
    "rho_proposed", "rho_accepted", "eta_proposed", "eta_accepted",
    "A_proposed", "A_accepted"
  )
  if (is.null(total)) {
    total <- setNames(as.list(rep(0, length(count_names))), count_names)
    if (!is.null(current$r_accepted_by_coordinate)) {
      total$r_accepted_by_coordinate <-
        rep(0, length(current$r_accepted_by_coordinate))
      total$r_proposed_per_coordinate <- 0
    }
  }
  for (name in count_names) total[[name]] <- total[[name]] + current[[name]]
  if (!is.null(current$r_accepted_by_coordinate)) {
    total$r_accepted_by_coordinate <-
      total$r_accepted_by_coordinate + current$r_accepted_by_coordinate
    total$r_proposed_per_coordinate <-
      total$r_proposed_per_coordinate +
      current$r_proposed / length(current$r_accepted_by_coordinate)
  }
  total
}

finalize_acceptance <- function(counts) {
  output <- counts
  for (block in c("beta", "theta", "r", "V", "rho", "eta", "A")) {
    output[[paste0(block, "_rate")]] <-
      counts[[paste0(block, "_accepted")]] /
      counts[[paste0(block, "_proposed")]]
  }
  if (!is.null(counts$r_accepted_by_coordinate)) {
    output$r_rate_by_coordinate <-
      counts$r_accepted_by_coordinate / counts$r_proposed_per_coordinate
  }
  output
}

# One complete independent ladder --------------------------------------------

run_one_ladder <- function(ladder_id) {
  ladder_seed <- ladder_seeds[[ladder_id]]
  regime <- initialization_regimes[[ladder_id]]
  ladder_file <- file.path(
    ladder_dir, paste0("ladder_", ladder_id, "_cold_chain.rds")
  )
  diagnostics_file <- file.path(
    ladder_dir, paste0("ladder_", ladder_id, "_diagnostics.rds")
  )
  final_tempering_state_file <- file.path(
    ladder_dir,
    paste0("ladder_", ladder_id, "_final_tempering_state.rds")
  )
  sampling_chunk_dir <- file.path(
    checkpoint_dir, paste0("ladder_", ladder_id, "_sampling_chunks")
  )

  if (reuse_completed_ladders && file.exists(ladder_file) &&
      file.exists(diagnostics_file) &&
      file.exists(final_tempering_state_file)) {
    existing <- tryCatch(readRDS(ladder_file), error = identity)
    if (!inherits(existing, "error") &&
        identical(existing$sampler_md5, sampler_md5) &&
        identical(existing$settings$run_mode, run_mode) &&
        identical(
          if (is.null(existing$settings$model_specification)) {
            "adj"
          } else {
            existing$settings$model_specification
          },
          model_specification
        ) &&
        identical(existing$settings$n_temperatures, n_temperatures) &&
        identical(existing$settings$warmup_iterations, warmup_iterations) &&
        identical(existing$settings$retained_draws, retained_draws) &&
        identical(existing$settings$sampling_thin, sampling_thin)) {
      return(list(status = "reused", ladder = ladder_id,
                  chain_file = ladder_file,
                  diagnostics_file = diagnostics_file,
                  final_tempering_state_file = final_tempering_state_file,
                  elapsed_seconds = NA_real_))
    }
  }

  dir.create(sampling_chunk_dir, recursive = TRUE, showWarnings = FALSE)
  stale_chunks <- list.files(
    sampling_chunk_dir,
    pattern = "^chunk_[0-9]+\\.rds$",
    full.names = TRUE
  )
  if (length(stale_chunks)) unlink(stale_chunks)

  set.seed(ladder_seed)
  started <- Sys.time()
  inverse_temperatures <- initial_inverse_temperatures()
  states <- vector("list", n_temperatures)
  proposals <- vector("list", n_temperatures)
  log_likelihoods <- rep(NA_real_, n_temperatures)
  adaptation_offsets <- integer(n_temperatures)
  labels <- seq_len(n_temperatures)
  warmup_phase_names <- c(
    "stabilization", "joint_adaptation",
    "temperature_adaptation", "validation"
  )
  warmup_phase_proposed <- setNames(lapply(
    warmup_phase_names,
    function(x) integer(n_temperatures - 1L)
  ), warmup_phase_names)
  warmup_phase_accepted <- setNames(lapply(
    warmup_phase_names,
    function(x) integer(n_temperatures - 1L)
  ), warmup_phase_names)
  warmup_phase_occupancy <- setNames(lapply(
    warmup_phase_names,
    function(x) matrix(0L, n_temperatures, n_temperatures)
  ), warmup_phase_names)
  warmup_phase_traffic_state <- setNames(lapply(
    warmup_phase_names,
    function(x) integer(n_temperatures)
  ), warmup_phase_names)
  warmup_phase_round_trips <- setNames(lapply(
    warmup_phase_names,
    function(x) integer(n_temperatures)
  ), warmup_phase_names)
  window_proposed <- integer(n_temperatures - 1L)
  window_accepted <- integer(n_temperatures - 1L)
  adaptation_history <- list()
  initial_values <- NULL
  consistency_checks <- 0L
  maximum_F_gap <- 0

  warmup_rounds <- warmup_iterations %/% local_sweeps_per_swap
  stabilization_rounds <- max(
    1L, as.integer(floor(warmup_rounds * stabilization_fraction))
  )
  joint_adaptation_rounds <- max(
    1L, as.integer(floor(warmup_rounds * joint_adaptation_fraction))
  )
  temperature_only_adaptation_rounds <- max(
    1L,
    as.integer(floor(
      warmup_rounds * temperature_only_adaptation_fraction
    ))
  )
  allocated_adaptation_rounds <-
    stabilization_rounds + joint_adaptation_rounds +
    temperature_only_adaptation_rounds
  if (allocated_adaptation_rounds >= warmup_rounds) {
    temperature_only_adaptation_rounds <-
      warmup_rounds - stabilization_rounds -
      joint_adaptation_rounds - 1L
  }
  validation_rounds <-
    warmup_rounds - stabilization_rounds - joint_adaptation_rounds -
    temperature_only_adaptation_rounds
  adaptation_first_round <- stabilization_rounds + 1L
  joint_adaptation_last_round <-
    stabilization_rounds + joint_adaptation_rounds
  adaptation_last_round <-
    joint_adaptation_last_round + temperature_only_adaptation_rounds
  adaptation_step <- 0L
  adaptation_gain_index <- 0L
  for (round_index in seq_len(warmup_rounds)) {
    phase <- if (round_index <= stabilization_rounds) {
      "stabilization"
    } else if (round_index <= joint_adaptation_last_round) {
      "joint_adaptation"
    } else if (round_index <= adaptation_last_round) {
      "temperature_adaptation"
    } else {
      "validation"
    }
    if (round_index == joint_adaptation_last_round + 1L) {
      # Proposal scales are now fixed, so restart the temperature-adaptation
      # gain to respond to the new frozen-kernel exchange landscape.
      adaptation_gain_index <- 0L
    }
    proposal_adaptation_active <- phase %in% c(
      "stabilization", "joint_adaptation"
    )
    for (temperature in seq_len(n_temperatures)) {
      segment <- run_segment(
        states[[temperature]], proposals[[temperature]],
        inverse_temperatures[[temperature]],
        ladder_seed + 100000L + round_index * 20L + temperature,
        ladder_id * 100L + temperature,
        iterations = local_sweeps_per_swap,
        retain = 1L,
        thin = 1L,
        adaptation_offset = adaptation_offsets[[temperature]],
        adapt_proposals = proposal_adaptation_active,
        initialization_regime = regime
      )
      if (is.null(initial_values) && temperature == 1L) {
        initial_values <- segment$initial_values
      }
      states[[temperature]] <- segment$final_state
      proposals[[temperature]] <- segment$proposal_state
      log_likelihoods[[temperature]] <- segment$log_likelihood
      adaptation_offsets[[temperature]] <-
        segment$proposal_state$adaptation_iterations
      consistency_checks <- consistency_checks +
        segment$state_consistency$checks
      maximum_F_gap <- max(
        maximum_F_gap,
        segment$state_consistency$maximum_F_gap_seen
      )
      if (!isTRUE(segment$state_consistency$all_checks_passed)) {
        stop("A warm-up segment failed state consistency.")
      }
    }

    before_proposed <- warmup_phase_proposed[[phase]]
    before_accepted <- warmup_phase_accepted[[phase]]
    swapped <- perform_swaps(
      states, log_likelihoods, labels, inverse_temperatures,
      round_index,
      warmup_phase_proposed[[phase]],
      warmup_phase_accepted[[phase]]
    )
    states <- swapped$states
    log_likelihoods <- swapped$log_likelihoods
    labels <- swapped$labels
    warmup_phase_proposed[[phase]] <- swapped$proposed
    warmup_phase_accepted[[phase]] <- swapped$accepted
    temperature_adaptation_active <- phase %in% c(
      "joint_adaptation", "temperature_adaptation"
    )
    if (temperature_adaptation_active) {
      window_proposed <- window_proposed +
        warmup_phase_proposed[[phase]] - before_proposed
      window_accepted <- window_accepted +
        warmup_phase_accepted[[phase]] - before_accepted
    }
    warmup_phase_occupancy[[phase]][
      cbind(labels, seq_len(n_temperatures))
    ] <- warmup_phase_occupancy[[phase]][
      cbind(labels, seq_len(n_temperatures))
    ] + 1L
    traffic <- update_round_trips(
      labels,
      warmup_phase_traffic_state[[phase]],
      warmup_phase_round_trips[[phase]]
    )
    warmup_phase_traffic_state[[phase]] <- traffic$traffic_phase
    warmup_phase_round_trips[[phase]] <- traffic$round_trips

    adaptation_phase_round <- round_index - stabilization_rounds
    adaptation_window_complete <-
      temperature_adaptation_active &&
      (adaptation_phase_round %% adaptation_window_rounds == 0L ||
         round_index == adaptation_last_round)
    if (adaptation_window_complete) {
      adaptation_step <- adaptation_step + 1L
      adaptation_gain_index <- adaptation_gain_index + 1L
      old_temperatures <- inverse_temperatures
      inverse_temperatures <- adapt_inverse_temperatures(
        inverse_temperatures,
        window_accepted,
        window_proposed,
        adaptation_gain_index
      )
      adaptation_history[[adaptation_step]] <- data.frame(
        adaptation = adaptation_step,
        adaptation_phase = phase,
        phase_adaptation = adaptation_gain_index,
        round = round_index,
        pair = seq_len(n_temperatures - 1L),
        beta_lower = old_temperatures[-n_temperatures],
        beta_upper = old_temperatures[-1L],
        proposed = window_proposed,
        accepted = window_accepted,
        acceptance = window_accepted / pmax(1L, window_proposed)
      )
      window_proposed[] <- 0L
      window_accepted[] <- 0L
      saveRDS(
        list(
          phase = "warmup", round = round_index,
          warmup_subphase = phase,
          inverse_temperatures = inverse_temperatures,
          states = states, proposals = proposals,
          log_likelihoods = log_likelihoods,
          labels = labels,
          round_trips = warmup_phase_round_trips,
          sampler_md5 = sampler_md5
        ),
        file.path(
          checkpoint_dir,
          paste0("ladder_", ladder_id, "_latest_checkpoint.rds")
        )
      )
    }

    if (round_index %% max(1L, warmup_rounds %/% 10L) == 0L ||
        round_index %in% c(stabilization_rounds, adaptation_last_round)) {
      phase_rates <- warmup_phase_accepted[[phase]] /
        pmax(1L, warmup_phase_proposed[[phase]])
      cat(
        "### Ladder", ladder_id, "warmup", phase,
        round_index, "/", warmup_rounds,
        "| beta:", paste(round(inverse_temperatures, 4), collapse = ","),
        "| phase swaps:", paste(round(phase_rates, 2), collapse = ","),
        "\n"
      )
    }
  }

  frozen_inverse_temperatures <- inverse_temperatures
  warmup_proposed <- Reduce(`+`, warmup_phase_proposed)
  warmup_accepted <- Reduce(`+`, warmup_phase_accepted)
  validation_swap_acceptance <-
    warmup_phase_accepted$validation /
    pmax(1L, warmup_phase_proposed$validation)
  learned_beta_hot <- tail(frozen_inverse_temperatures, 1L)
  hot_endpoint_hit_bound <-
    endpoint_strategy == "adaptive" &&
    (learned_beta_hot <= minimum_learned_beta_hot * (1 + 1e-8) ||
       learned_beta_hot >= maximum_learned_beta_hot * (1 - 1e-8))
  if (hot_endpoint_hit_bound) {
    warning(
      "The learned hot endpoint reached a numerical safeguard. ",
      "Treat this ladder calibration as failed."
    )
  }
  sampling_proposed <- integer(n_temperatures - 1L)
  sampling_accepted <- integer(n_temperatures - 1L)
  sampling_occupancy <- matrix(0L, n_temperatures, n_temperatures)
  sampling_traffic_state <- integer(n_temperatures)
  sampling_traffic_state[[labels[[1L]]]] <- 1L
  sampling_round_trips <- integer(n_temperatures)
  temperature_cardinality_sum <- matrix(0, n_temperatures, model_arguments$q)
  temperature_cardinality_sum_squares <- matrix(
    0, n_temperatures, model_arguments$q
  )
  temperature_cardinality_min <- matrix(
    Inf, n_temperatures, model_arguments$q
  )
  temperature_cardinality_max <- matrix(
    -Inf, n_temperatures, model_arguments$q
  )
  temperature_cardinality_changes <- matrix(
    0L, n_temperatures, model_arguments$q
  )
  temperature_cardinality_last <- vector("list", n_temperatures)
  temperature_cardinality_draws <- integer(n_temperatures)
  sampling_rounds <- as.integer(
    retained_draws * sampling_thin / sampling_sweeps_per_swap
  )
  cold_segments <- list()
  cold_draw_index <- 0L
  cold_chunk_index <- 0L
  cold_chunk_files <- character()
  cold_W_sum <- lapply(
    seq_len(model_arguments$q),
    function(disease) matrix(0, nrow(model_arguments$Minc), ncol(model_arguments$Minc))
  )
  cold_W_draws <- 0L
  cold_acceptance <- NULL

  for (round_index in seq_len(sampling_rounds)) {
    segment_fits <- vector("list", n_temperatures)
    for (temperature in seq_len(n_temperatures)) {
      segment <- run_segment(
        states[[temperature]], proposals[[temperature]],
        frozen_inverse_temperatures[[temperature]],
        ladder_seed + 100000000L + round_index * 20L + temperature,
        ladder_id * 100L + temperature,
        iterations = sampling_sweeps_per_swap,
        retain = 1L,
        thin = sampling_sweeps_per_swap,
        adaptation_offset = adaptation_offsets[[temperature]],
        adapt_proposals = FALSE,
        initialization_regime = regime
      )
      states[[temperature]] <- segment$final_state
      proposals[[temperature]] <- segment$proposal_state
      log_likelihoods[[temperature]] <- segment$log_likelihood
      segment_fits[[temperature]] <- segment
      cardinality <- segment$W_cardinality
      temperature_cardinality_sum[temperature, ] <-
        temperature_cardinality_sum[temperature, ] + colSums(cardinality)
      temperature_cardinality_sum_squares[temperature, ] <-
        temperature_cardinality_sum_squares[temperature, ] +
        colSums(cardinality^2)
      temperature_cardinality_min[temperature, ] <- pmin(
        temperature_cardinality_min[temperature, ],
        apply(cardinality, 2L, min)
      )
      temperature_cardinality_max[temperature, ] <- pmax(
        temperature_cardinality_max[temperature, ],
        apply(cardinality, 2L, max)
      )
      for (disease in seq_len(model_arguments$q)) {
        trajectory <- cardinality[, disease]
        if (!is.null(temperature_cardinality_last[[temperature]])) {
          trajectory <- c(
            temperature_cardinality_last[[temperature]][[disease]],
            trajectory
          )
        }
        temperature_cardinality_changes[temperature, disease] <-
          temperature_cardinality_changes[temperature, disease] +
          sum(diff(trajectory) != 0)
      }
      temperature_cardinality_last[[temperature]] <-
        cardinality[nrow(cardinality), ]
      temperature_cardinality_draws[[temperature]] <-
        temperature_cardinality_draws[[temperature]] + nrow(cardinality)
      if (!isTRUE(segment$state_consistency$all_checks_passed)) {
        stop("A sampling segment failed state consistency.")
      }
    }
    cold_segment <- segment_fits[[1L]]
    cold_acceptance <- accumulate_acceptance(
      cold_acceptance, cold_segment$acceptance
    )
    retain_cold_draw <-
      (round_index * sampling_sweeps_per_swap) %% sampling_thin == 0L
    if (retain_cold_draw) {
      cold_draw_index <- cold_draw_index + 1L
      cold_segments[[length(cold_segments) + 1L]] <-
        strip_cold_segment(cold_segment)
      segment_draws <- nrow(cold_segment$beta)
      for (disease in seq_len(model_arguments$q)) {
        cold_W_sum[[disease]] <- cold_W_sum[[disease]] +
          cold_segment$W_mean[[disease]] * segment_draws
      }
      cold_W_draws <- cold_W_draws + segment_draws
      flush_cold_chunk <-
        length(cold_segments) >= sampling_chunk_draws ||
        cold_draw_index == retained_draws
      if (flush_cold_chunk) {
        cold_chunk_index <- cold_chunk_index + 1L
        chunk <- lapply(matrix_fields, function(field) {
          do.call(rbind, lapply(cold_segments, `[[`, field))
        })
        names(chunk) <- matrix_fields
        chunk$tau <- unlist(
          lapply(cold_segments, `[[`, "tau"), use.names = FALSE
        )
        chunk$last_retained_draw <- cold_draw_index
        chunk$first_retained_draw <-
          cold_draw_index - length(cold_segments) + 1L
        chunk_file <- file.path(
          sampling_chunk_dir,
          sprintf("chunk_%04d.rds", cold_chunk_index)
        )
        saveRDS(chunk, chunk_file)
        cold_chunk_files <- c(cold_chunk_files, chunk_file)
        cold_segments <- list()
        rm(chunk)
      }
    }

    swapped <- perform_swaps(
      states, log_likelihoods, labels, frozen_inverse_temperatures,
      round_index, sampling_proposed, sampling_accepted
    )
    states <- swapped$states
    log_likelihoods <- swapped$log_likelihoods
    labels <- swapped$labels
    sampling_proposed <- swapped$proposed
    sampling_accepted <- swapped$accepted
    sampling_occupancy[cbind(labels, seq_len(n_temperatures))] <-
      sampling_occupancy[cbind(labels, seq_len(n_temperatures))] + 1L
    traffic <- update_round_trips(
      labels, sampling_traffic_state, sampling_round_trips
    )
    sampling_traffic_state <- traffic$traffic_phase
    sampling_round_trips <- traffic$round_trips

    if (round_index %% max(1L, sampling_rounds %/% 20L) == 0L) {
      cat(
        "### Ladder", ladder_id, "sampling",
        round_index, "/", sampling_rounds,
        "| swaps:",
        paste(round(sampling_accepted / pmax(1L, sampling_proposed), 2), collapse = ","),
        "| round trips:", sum(sampling_round_trips), "\n"
      )
    }
  }

  if (cold_draw_index != retained_draws) {
    stop(
      "The number of retained cold draws does not match the requested total."
    )
  }

  cold_chunks <- lapply(cold_chunk_files, readRDS)
  cold <- lapply(matrix_fields, function(field) {
    do.call(rbind, lapply(cold_chunks, `[[`, field))
  })
  names(cold) <- matrix_fields
  cold$tau <- unlist(
    lapply(cold_chunks, `[[`, "tau"), use.names = FALSE
  )
  cold$W_mean <- lapply(cold_W_sum, `/`, cold_W_draws)
  rm(cold_chunks)
  gc(verbose = FALSE)
  temperature_cardinality_mean <-
    temperature_cardinality_sum / temperature_cardinality_draws
  temperature_cardinality_variance <-
    temperature_cardinality_sum_squares / temperature_cardinality_draws -
    temperature_cardinality_mean^2
  temperature_cardinality_variance[
    temperature_cardinality_variance < 0
  ] <- 0
  temperature_cardinality_summary <- do.call(rbind, lapply(
    seq_len(n_temperatures),
    function(temperature) {
      data.frame(
        temperature = temperature,
        inverse_temperature = frozen_inverse_temperatures[[temperature]],
        disease = c("Lung", "Esophageal", "Larynx", "Colorectal"),
        mean = temperature_cardinality_mean[temperature, ],
        sd = sqrt(temperature_cardinality_variance[temperature, ]),
        minimum = temperature_cardinality_min[temperature, ],
        maximum = temperature_cardinality_max[temperature, ],
        changes = temperature_cardinality_changes[temperature, ],
        change_rate = temperature_cardinality_changes[temperature, ] /
          pmax(1L, temperature_cardinality_draws[[temperature]] - 1L)
      )
    }
  ))
  cold_chain <- c(
    list(
      chain_id = ladder_id,
      seed = ladder_seed,
      initial_values = initial_values,
      acceptance = finalize_acceptance(cold_acceptance),
      state_consistency = list(
        checks = consistency_checks,
        tolerance = state_consistency_tolerance,
        maximum_F_gap_seen = maximum_F_gap,
        maximum_u_mismatches_seen = 0L,
        maximum_W_eta_mismatches_seen = 0L,
        all_checks_passed = maximum_F_gap <= state_consistency_tolerance
      ),
      final_state = states[[1L]],
      log_likelihood = log_likelihoods[[1L]],
      settings = list(
        model_specification = model_specification,
        run_mode = run_mode,
        sampler = "adaptive_parallel_tempering",
        endpoint_strategy = endpoint_strategy,
        initialization_regime = regime,
        n_temperatures = n_temperatures,
        inverse_temperatures = frozen_inverse_temperatures,
        initial_beta_hot = initial_beta_hot,
        learned_beta_hot = learned_beta_hot,
        target_swap_acceptance = target_swap_acceptance,
        warmup_iterations = warmup_iterations,
        stabilization_rounds = stabilization_rounds,
        joint_adaptation_rounds = joint_adaptation_rounds,
        temperature_only_adaptation_rounds =
          temperature_only_adaptation_rounds,
        validation_rounds = validation_rounds,
        retained_draws = retained_draws,
        sampling_thin = sampling_thin,
        sampling_sweeps_per_swap = sampling_sweeps_per_swap,
        complete_exchange_cycle = complete_exchange_cycle,
        sampling_chunk_draws = sampling_chunk_draws,
        sampling_chunk_files = cold_chunk_files
      )
    ),
    cold
  )
  cold_chain$sampler_md5 <- sampler_md5
  cold_chain$runtime <- list(
    start = started, finish = Sys.time(),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )

  ladder_diagnostics <- list(
    ladder = ladder_id,
    model_specification = model_specification,
    initial_inverse_temperatures = initial_inverse_temperatures(),
    frozen_inverse_temperatures = frozen_inverse_temperatures,
    adaptation_history = if (length(adaptation_history)) {
      do.call(rbind, adaptation_history)
    } else data.frame(),
    warmup_swap_proposed = warmup_proposed,
    warmup_swap_accepted = warmup_accepted,
    warmup_swap_acceptance = warmup_accepted / pmax(1L, warmup_proposed),
    warmup_phase_swap_proposed = warmup_phase_proposed,
    warmup_phase_swap_accepted = warmup_phase_accepted,
    warmup_phase_swap_acceptance = Map(
      function(accepted, proposed) accepted / pmax(1L, proposed),
      warmup_phase_accepted,
      warmup_phase_proposed
    ),
    validation_swap_acceptance = validation_swap_acceptance,
    warmup_phase_replica_temperature_occupancy = warmup_phase_occupancy,
    warmup_phase_round_trips_by_replica = warmup_phase_round_trips,
    sampling_swap_proposed = sampling_proposed,
    sampling_swap_accepted = sampling_accepted,
    sampling_swap_acceptance = sampling_accepted / pmax(1L, sampling_proposed),
    replica_temperature_occupancy = sampling_occupancy,
    round_trips_by_replica = sampling_round_trips,
    total_round_trips = sum(sampling_round_trips),
    learned_beta_hot = learned_beta_hot,
    hot_endpoint_hit_bound = hot_endpoint_hit_bound,
    temperature_cardinality_summary = temperature_cardinality_summary
  )
  final_tempering_state <- list(
    ladder = ladder_id,
    seed = ladder_seed,
    model_specification = model_specification,
    run_mode = run_mode,
    sampler_md5 = sampler_md5,
    inverse_temperatures = frozen_inverse_temperatures,
    states = states,
    proposals = proposals,
    log_likelihoods = log_likelihoods,
    labels = labels,
    adaptation_offsets = adaptation_offsets,
    rng_state = if (exists(".Random.seed", inherits = FALSE)) {
      .Random.seed
    } else {
      NULL
    },
    completed_sampling_rounds = sampling_rounds,
    retained_draws = retained_draws,
    sampling_thin = sampling_thin,
    sampling_sweeps_per_swap = sampling_sweeps_per_swap,
    complete_exchange_cycle = complete_exchange_cycle,
    saved_at = Sys.time()
  )
  saveRDS(cold_chain, ladder_file)
  saveRDS(ladder_diagnostics, diagnostics_file)
  saveRDS(final_tempering_state, final_tempering_state_file)
  list(
    status = "ok", ladder = ladder_id,
    chain_file = ladder_file, diagnostics_file = diagnostics_file,
    final_tempering_state_file = final_tempering_state_file,
    elapsed_seconds = cold_chain$runtime$elapsed_seconds
  )
}

# Execute ---------------------------------------------------------------------

cat("\nAdaptive parallel-tempering workflow\n")
cat("Start:", format(Sys.time()), "\n")
cat("SEER specification:", model_specification, "\n")
cat("Mode:", run_mode, "\n")
cat("Ladders:", n_ladders, "\n")
cat("Temperatures per ladder:", n_temperatures, "\n")
cat("Initial beta ladder:", paste(round(initial_inverse_temperatures(), 5), collapse = ", "), "\n")
cat("Fixed cold endpoint:", beta_cold, "\n")
cat("Endpoint strategy:", endpoint_strategy, "\n")
if (endpoint_strategy == "fixed") {
  cat("Fixed hot endpoint:", initial_beta_hot, "\n")
  cat("Ladder objective: equalize adjacent-swap acceptance\n")
} else {
  cat("Initial (adaptive) hot endpoint:", initial_beta_hot, "\n")
  cat("Target adjacent-swap acceptance:", target_swap_acceptance, "\n")
}
cat("Warm-up sweeps per replica:", warmup_iterations, "\n")
cat(
  paste0(
    "Warm-up fractions ",
    "(stabilize/joint-adapt/temperature-only/validate):"
  ),
  stabilization_fraction, "/", joint_adaptation_fraction, "/",
  temperature_only_adaptation_fraction, "/",
  1 - stabilization_fraction - joint_adaptation_fraction -
    temperature_only_adaptation_fraction,
  "\n"
)
cat("Retained cold draws per ladder:", retained_draws, "\n")
cat("Sampling thinning:", sampling_thin, "\n")
cat("Sampling sweeps per swap:", sampling_sweeps_per_swap, "\n")
cat("Complete odd-even exchange cycle:", complete_exchange_cycle, "\n")
cat("Retained draws per disk chunk:", sampling_chunk_draws, "\n")
cat("Output:", output_dir, "\n\n")

Rcpp::sourceCpp(sampler_path)

if (n_ladders == 1L) {
  run_status <- list(run_one_ladder(1L))
} else {
  workers <- min(n_ladders, parallel::detectCores(logical = FALSE))
  cl <- parallel::makeCluster(workers, outfile = "")
  export_names <- c(
    "model_arguments", "sampler_path", "sampler_md5",
    "model_specification", "run_mode", "n_temperatures", "beta_cold",
    "initial_beta_hot",
    "endpoint_strategy", "complete_exchange_cycle",
    "target_swap_acceptance", "minimum_acceptable_swap_acceptance",
    "adaptation_exponent", "adaptation_gain",
    "minimum_learned_beta_hot", "maximum_learned_beta_hot",
    "stabilization_fraction", "joint_adaptation_fraction",
    "temperature_only_adaptation_fraction",
    "state_consistency_tolerance", "reuse_completed_ladders",
    "warmup_iterations", "retained_draws", "sampling_thin",
    "local_sweeps_per_swap", "sampling_sweeps_per_swap",
    "retained_draws_per_swap", "sampling_chunk_draws",
    "adaptation_window_rounds", "ladder_seeds", "initialization_regimes",
    "ladder_dir", "checkpoint_dir", "matrix_fields",
    "initial_inverse_temperatures", "adapt_inverse_temperatures",
    "swap_pairs_for_round", "perform_swaps", "update_round_trips",
    "run_segment", "strip_cold_segment", "combine_cold_segments",
    "accumulate_acceptance",
    "finalize_acceptance", "run_one_ladder"
  )
  parallel::clusterExport(cl, export_names, envir = environment())
  # Avoid simultaneous Rtools writes on Windows.
  for (worker in seq_len(workers)) {
    one <- cl[worker]
    class(one) <- class(cl)
    parallel::clusterCall(one, function(path) {
      library(Rcpp)
      library(RcppArmadillo)
      Rcpp::sourceCpp(path)
      TRUE
    }, sampler_path)
  }
  run_status <- tryCatch(
    parallel::parLapply(cl, seq_len(n_ladders), run_one_ladder),
    finally = parallel::stopCluster(cl)
  )
}

status_table <- do.call(rbind, lapply(run_status, as.data.frame))
write.csv(
  status_table,
  diagnostic_table_file("ladder_run_status"),
  row.names = FALSE
)
if (any(!status_table$status %in% c("ok", "reused"))) {
  stop("At least one ladder failed; completed ladders remain saved.")
}

cold_chains <- lapply(status_table$chain_file, readRDS)
ladder_diagnostics <- lapply(status_table$diagnostics_file, readRDS)
final_tempering_states <- lapply(
  status_table$final_tempering_state_file,
  readRDS
)

learned_ladders <- do.call(rbind, lapply(ladder_diagnostics, function(x) {
  data.frame(
    ladder = x$ladder,
    temperature = seq_len(n_temperatures),
    inverse_temperature = x$frozen_inverse_temperatures
  )
}))
swap_summary <- do.call(rbind, lapply(ladder_diagnostics, function(x) {
  data.frame(
    ladder = x$ladder,
    pair = seq_len(n_temperatures - 1L),
    beta_lower = x$frozen_inverse_temperatures[-n_temperatures],
    beta_upper = x$frozen_inverse_temperatures[-1L],
    proposed = x$sampling_swap_proposed,
    accepted = x$sampling_swap_accepted,
    acceptance = x$sampling_swap_acceptance
  )
}))
validation_swap_summary <- do.call(rbind, lapply(
  ladder_diagnostics,
  function(x) {
    data.frame(
      ladder = x$ladder,
      pair = seq_len(n_temperatures - 1L),
      beta_lower = x$frozen_inverse_temperatures[-n_temperatures],
      beta_upper = x$frozen_inverse_temperatures[-1L],
      proposed = x$warmup_phase_swap_proposed$validation,
      accepted = x$warmup_phase_swap_accepted$validation,
      acceptance = x$validation_swap_acceptance
    )
  }
))
round_trip_summary <- do.call(rbind, lapply(ladder_diagnostics, function(x) {
  data.frame(
    ladder = x$ladder,
    replica = seq_len(n_temperatures),
    round_trips = x$round_trips_by_replica
  )
}))
write.csv(learned_ladders, diagnostic_table_file("learned_temperature_ladders"), row.names = FALSE)
write.csv(validation_swap_summary, diagnostic_table_file("validation_swap_acceptance"), row.names = FALSE)
write.csv(swap_summary, diagnostic_table_file("sampling_swap_acceptance"), row.names = FALSE)
write.csv(round_trip_summary, diagnostic_table_file("replica_round_trips"), row.names = FALSE)

ladder_quality <- do.call(rbind, lapply(ladder_diagnostics, function(x) {
  validation_reference <- if (endpoint_strategy == "adaptive") {
    target_swap_acceptance
  } else {
    mean(x$validation_swap_acceptance)
  }
  data.frame(
    ladder = x$ladder,
    learned_beta_hot = x$learned_beta_hot,
    hot_endpoint_hit_bound = x$hot_endpoint_hit_bound,
    minimum_validation_swap = min(x$validation_swap_acceptance),
    mean_validation_swap = mean(x$validation_swap_acceptance),
    maximum_validation_rate_deviation = max(abs(
      x$validation_swap_acceptance - validation_reference
    )),
    minimum_sampling_swap = min(x$sampling_swap_acceptance),
    mean_sampling_swap = mean(x$sampling_swap_acceptance),
    sampling_round_trips = x$total_round_trips,
    no_severe_validation_bottleneck =
      min(x$validation_swap_acceptance) >=
      minimum_acceptable_swap_acceptance,
    no_severe_sampling_bottleneck =
      min(x$sampling_swap_acceptance) >=
      minimum_acceptable_swap_acceptance
  )
}))
write.csv(
  ladder_quality,
  diagnostic_table_file("ladder_quality_gate"),
  row.names = FALSE
)
temperature_cardinality_summary <- do.call(rbind, lapply(
  ladder_diagnostics,
  function(x) {
    result <- x$temperature_cardinality_summary
    result$ladder <- x$ladder
    result[, c("ladder", setdiff(names(result), "ladder"))]
  }
))
write.csv(
  temperature_cardinality_summary,
  diagnostic_table_file("cardinality_movement_by_temperature"),
  row.names = FALSE
)

cardinality_summary <- do.call(rbind, lapply(seq_along(cold_chains), function(chain) {
  data.frame(
    ladder = chain,
    disease = c("Lung", "Esophageal", "Larynx", "Colorectal"),
    mean = colMeans(cold_chains[[chain]]$W_cardinality),
    sd = apply(cold_chains[[chain]]$W_cardinality, 2L, sd)
  )
}))
write.csv(cardinality_summary, diagnostic_table_file("cold_adjacency_cardinality"), row.names = FALSE)

# Scalar-r adaptation diagnostics.  Proposal states remain attached to
# temperatures during swaps, so these scales are indexed by ladder and
# temperature rather than by replica label.
r_proposal_summary <- do.call(rbind, lapply(
  seq_along(final_tempering_states),
  function(ladder) {
    do.call(rbind, lapply(seq_len(n_temperatures), function(temperature) {
      proposal <- final_tempering_states[[ladder]]$proposals[[temperature]]
      if (is.null(proposal$s_r) || is.null(proposal$R_r)) {
        stop("A final tempering state lacks scalar-r adaptation fields.")
      }
      data.frame(
        ladder = ladder,
        temperature = temperature,
        coordinate = seq_along(proposal$s_r),
        s_r = proposal$s_r,
        R_r = proposal$R_r,
        proposal_sd = sqrt(exp(proposal$s_r) * proposal$R_r)
      )
    }))
  }
))
write.csv(
  r_proposal_summary,
  diagnostic_table_file("r_proposal_scales_by_temperature"),
  row.names = FALSE
)

cold_r_acceptance_by_coordinate <- do.call(rbind, lapply(
  seq_along(cold_chains),
  function(ladder) {
    rates <- cold_chains[[ladder]]$acceptance$r_rate_by_coordinate
    if (is.null(rates)) {
      stop("A cold chain lacks coordinate-specific scalar-r acceptance.")
    }
    data.frame(
      ladder = ladder,
      coordinate = seq_along(rates),
      acceptance = rates
    )
  }
))
write.csv(
  cold_r_acceptance_by_coordinate,
  diagnostic_table_file("cold_r_acceptance_by_coordinate"),
  row.names = FALSE
)

saveRDS(
  list(
    settings = list(
      model_specification = model_specification,
      run_mode = run_mode, n_ladders = n_ladders,
      n_temperatures = n_temperatures,
      endpoint_strategy = endpoint_strategy,
      initial_beta_hot = initial_beta_hot,
      target_swap_acceptance = target_swap_acceptance,
      warmup_iterations = warmup_iterations,
      retained_draws = retained_draws, sampling_thin = sampling_thin,
      sampling_sweeps_per_swap = sampling_sweeps_per_swap,
      complete_exchange_cycle = complete_exchange_cycle,
      sampling_chunk_draws = sampling_chunk_draws,
      sampler_md5 = sampler_md5
    ),
    cold_chains = cold_chains,
    ladder_diagnostics = ladder_diagnostics,
    learned_ladders = learned_ladders,
    validation_swap_summary = validation_swap_summary,
    swap_summary = swap_summary,
    round_trip_summary = round_trip_summary,
    ladder_quality = ladder_quality,
    temperature_cardinality_summary = temperature_cardinality_summary,
    cardinality_summary = cardinality_summary,
    r_proposal_summary = r_proposal_summary,
    cold_r_acceptance_by_coordinate = cold_r_acceptance_by_coordinate
  ),
  file.path(output_dir, "tempering_results.rds")
)

cat("\nWorkflow complete:", format(Sys.time()), "\n")
cat("Learned inverse temperatures:\n")
print(learned_ladders)
cat("Frozen-ladder validation swap acceptance:\n")
print(validation_swap_summary)
cat("Sampling swap acceptance:\n")
print(swap_summary)
cat("Round trips:\n")
print(round_trip_summary)
cat("Ladder quality gate:\n")
print(ladder_quality)
cat("Saved in:", output_dir, "\n")
