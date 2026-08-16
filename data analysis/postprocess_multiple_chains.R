# Post-processing entry point for the non-tempered SEER chains
# =============================================================================
# This script does not run MCMC. It reads the six chain files produced by
# main_multiple_chains.R and writes directly comparable convergence,
# FDR, boundary-stability, and paper-style map outputs.
#
# Automatic input:
#   data analysis/multiple_chains_output/<specification>/<mode>
#
# Automatic output:
#   data analysis/multiple_chains_output/<specification>/<mode>/
#   postprocessing
#
# Derived diagnostics, figures, and summary objects are written under the
# run-specific postprocessing directory. Completed chain files, run status,
# preliminary run diagnostics, and multiple_chains_all_results.rds
# are never replaced.
#
# The specification and mode are read from the same variables as the sampler:
#   MCMC_MULTICHAIN_SPECIFICATION = adj, meanadj, or mean
#   MCMC_MULTICHAIN_MODE = smoke, quick, pilot, or production
# Direct input/output overrides remain available for exceptional cases, but
# they are not needed for ordinary use.
#
# Run from the repository root with either:
#   source("data analysis/postprocess_multiple_chains.R")
# or:
#   Rscript "data analysis/postprocess_multiple_chains.R"

run_non_tempered_postprocessing <- function() {
  previous_workflow <- Sys.getenv(
    "MCMC_POSTPROCESS_WORKFLOW", unset = NA_character_
  )
  on.exit({
    if (is.na(previous_workflow)) {
      Sys.unsetenv("MCMC_POSTPROCESS_WORKFLOW")
    } else {
      Sys.setenv(MCMC_POSTPROCESS_WORKFLOW = previous_workflow)
    }
  }, add = TRUE)

  Sys.setenv(MCMC_POSTPROCESS_WORKFLOW = "non_tempered")
  source(
    file.path(
      "data analysis",
      "postprocess_multiple_chains_tempering.R"
    ),
    local = new.env(parent = globalenv()),
    chdir = FALSE
  )
}

run_non_tempered_postprocessing()
rm(run_non_tempered_postprocessing)
