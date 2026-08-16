# Validate and inventory every SEER table and figure used by the manuscript.
# This script never runs MCMC or recomputes posterior summaries.

rm(list = ls())

if (!all(file.exists(c("README.md", "data analysis")))) {
  stop("Run this script from the repository root.", call. = FALSE)
}

tempering_root <- file.path(
  "data analysis", "multiple_chains_tempering_output"
)
validation_dir <- file.path("data analysis", "validation")
sensitivity_dir <- file.path("data analysis", "sensitivity")
ppc_dir <- file.path("data analysis", "ppc")
sampler_sensitivity_dir <- file.path(
  "data analysis", "sampler_sensitivity"
)
exploratory_dir <- file.path("data analysis", "exploratory_figures")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

artifact_row <- function(
  artifact,
  artifact_type,
  specification,
  source_file,
  manuscript_destination,
  generator,
  reported_in,
  required_now = TRUE
) {
  data.frame(
    artifact = artifact,
    artifact_type = artifact_type,
    specification = specification,
    source_file = source_file,
    manuscript_destination = manuscript_destination,
    generator = generator,
    reported_in = reported_in,
    required_now = required_now,
    stringsAsFactors = FALSE
  )
}

postprocessing_dir <- function(specification) {
  file.path(
    tempering_root, specification, "hotter_final", "postprocessing"
  )
}

specification_artifacts <- function(specification, table_number) {
  root <- postprocessing_dir(specification)
  final_run_exists <- dir.exists(file.path(
    tempering_root, specification, "hotter_final"
  ))
  required_now <- specification == "adj" || final_run_exists
  posterior_covariance_file <- if (specification == "mean") {
    paste0("pooled_covariance_", specification, ".pdf")
  } else {
    paste0("pooled_eta_covariance_", specification, ".pdf")
  }
  rows <- list(
    artifact_row(
      paste0("pooled_fdr_boundary_maps_", specification),
      "figure", specification,
      file.path(
        root, "figures",
        paste0("pooled_fdr_boundary_maps_", specification, ".pdf")
      ),
      file.path(
        "Images", paste0("pooled_fdr_boundary_maps_", specification, ".pdf")
      ),
      "postprocess_multiple_chains_tempering.R",
      if (specification == "adj") "main" else "supplement",
      required_now
    ),
    artifact_row(
      paste0("shared_fdr_boundary_maps_", specification),
      "figure", specification,
      file.path(
        root, "figures",
        paste0("shared_fdr_boundary_maps_", specification, ".pdf")
      ),
      file.path(
        "Images", paste0("shared_fdr_boundary_maps_", specification, ".pdf")
      ),
      "postprocess_multiple_chains_tempering.R",
      if (specification == "adj") "main" else "supplement",
      required_now
    ),
    artifact_row(
      paste0("mutual_fdr_boundary_maps_", specification),
      "figure", specification,
      file.path(
        root, "figures",
        paste0("mutual_fdr_boundary_maps_", specification, ".pdf")
      ),
      file.path(
        "Images", paste0("mutual_fdr_boundary_maps_", specification, ".pdf")
      ),
      "postprocess_multiple_chains_tempering.R",
      if (specification == "adj") "main" else "supplement",
      required_now
    ),
    artifact_row(
      paste0("pooled_fdr_curves_", specification),
      "figure", specification,
      file.path(
        root, "supplement",
        paste0("pooled_fdr_curves_", specification, ".pdf")
      ),
      file.path(
        "Images", paste0("pooled_fdr_curves_", specification, ".pdf")
      ),
      "postprocess_multiple_chains_tempering.R",
      "supplement", required_now
    ),
    artifact_row(
      paste0("pooled_beta_theta_", specification),
      "figure", specification,
      file.path(
        root, "supplement",
        paste0("pooled_beta_theta_", specification, ".pdf")
      ),
      file.path(
        "Images", paste0("pooled_beta_theta_", specification, ".pdf")
      ),
      "postprocess_multiple_chains_tempering.R",
      "supplement", required_now
    ),
    artifact_row(
      paste0("pooled_gamma_", specification),
      "figure", specification,
      file.path(
        root, "supplement", paste0("pooled_gamma_", specification, ".pdf")
      ),
      file.path(
        "Images", paste0("pooled_gamma_", specification, ".pdf")
      ),
      "postprocess_multiple_chains_tempering.R",
      "supplement", required_now
    ),
    artifact_row(
      paste0("pooled_V_rho_", specification),
      "figure", specification,
      file.path(
        root, "supplement", paste0("pooled_V_rho_", specification, ".pdf")
      ),
      file.path(
        "Images", paste0("pooled_V_rho_", specification, ".pdf")
      ),
      "postprocess_multiple_chains_tempering.R",
      "supplement", required_now
    ),
    artifact_row(
      sub("\\.pdf$", "", posterior_covariance_file),
      "figure", specification,
      file.path(root, "supplement", posterior_covariance_file),
      file.path("Images", posterior_covariance_file),
      "postprocess_multiple_chains_tempering.R",
      "supplement", required_now
    ),
    artifact_row(
      paste0("pooled_table_", table_number, "_", specification, "_tex"),
      "table", specification,
      file.path(
        root, "supplement",
        paste0("pooled_table_", table_number, "_", specification, ".tex")
      ),
      file.path(
        "Tables", paste0("pooled_table_", table_number, "_", specification, ".tex")
      ),
      "postprocess_multiple_chains_tempering.R",
      "supplement", required_now
    ),
    artifact_row(
      paste0("pooled_table_", table_number, "_", specification, "_csv"),
      "table_source", specification,
      file.path(
        root, "supplement",
        paste0("pooled_table_", table_number, "_", specification, ".csv")
      ),
      NA_character_,
      "postprocess_multiple_chains_tempering.R",
      "source_data", required_now
    )
  )
  if (specification != "mean") {
    rows <- append(
      rows,
      list(artifact_row(
        paste0("pooled_non_adjacency_maps_", specification),
        "figure", specification,
        file.path(
          root, "figures",
          paste0("pooled_non_adjacency_maps_", specification, ".pdf")
        ),
        file.path(
          "Images", paste0("pooled_non_adjacency_maps_", specification, ".pdf")
        ),
        "postprocess_multiple_chains_tempering.R",
        if (specification == "adj") "main" else "supplement",
        required_now
      )),
      after = 3L
    )
  }
  do.call(rbind, rows)
}

exploratory_artifacts <- rbind(
  artifact_row(
    "figure_1_seer_cancer_sir_maps", "figure", "observed_data",
    file.path(exploratory_dir, "figure_1_seer_cancer_sir_maps.pdf"),
    file.path("Images", "figure_1_seer_cancer_sir_maps.pdf"),
    "exploratory_figures_seer.R", "main"
  ),
  artifact_row(
    "figure_2_seer_morans_i_by_distance_band", "figure", "observed_data",
    file.path(
      exploratory_dir, "figure_2_seer_morans_i_by_distance_band.pdf"
    ),
    file.path("Images", "figure_2_seer_morans_i_by_distance_band.pdf"),
    "exploratory_figures_seer.R", "main"
  ),
  artifact_row(
    "figure_3_seer_covariate_maps", "figure", "observed_data",
    file.path(exploratory_dir, "figure_3_seer_covariate_maps.pdf"),
    file.path("Images", "figure_3_seer_covariate_maps.pdf"),
    "exploratory_figures_seer.R", "main"
  )
)

artifacts <- rbind(
  exploratory_artifacts,
  specification_artifacts("adj", "S16"),
  specification_artifacts("meanadj", "S17"),
  specification_artifacts("mean", "S18")
)

sensitivity_artifacts <- rbind(
  artifact_row(
    "named_map_multichain_california_counties", "figure",
    "three_specifications",
    file.path(
      sensitivity_dir, "figures",
      "named_map_multichain_california_counties.pdf"
    ),
    file.path("Images", "named_map_multichain_california_counties.pdf"),
    "sensitivity_analysis_multichain.R", "supplement"
  ),
  artifact_row(
    "boundary_probability_agreement_multichain_three_specifications",
    "figure", "three_specifications",
    file.path(
      sensitivity_dir, "figures",
      "boundary_probability_agreement_multichain_three_specifications.pdf"
    ),
    file.path(
      "Images",
      "boundary_probability_agreement_multichain_three_specifications.pdf"
    ),
    "sensitivity_analysis_multichain.R", "supplement"
  ),
  artifact_row(
    "boundary_maps_multichain_meanadj_vs_mean", "figure",
    "meanadj_vs_mean",
    file.path(
      sensitivity_dir, "figures",
      "boundary_maps_multichain_meanadj_vs_mean.pdf"
    ),
    file.path("Images", "boundary_maps_multichain_meanadj_vs_mean.pdf"),
    "sensitivity_analysis_multichain.R", "supplement"
  ),
  artifact_row(
    "boundary_sensitivity_table_multichain_three_specifications_tex",
    "table", "three_specifications",
    file.path(
      sensitivity_dir, "tables",
      "boundary_sensitivity_table_multichain_three_specifications.tex"
    ),
    file.path(
      "Tables",
      "boundary_sensitivity_table_multichain_three_specifications.tex"
    ),
    "sensitivity_analysis_multichain.R", "main_and_supplement"
  ),
  artifact_row(
    "boundary_sensitivity_table_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      sensitivity_dir, "tables",
      "boundary_sensitivity_table_multichain_three_specifications.csv"
    ),
    NA_character_, "sensitivity_analysis_multichain.R", "source_data"
  )
)

sampler_sensitivity_artifacts <- rbind(
  artifact_row(
    "sampler_convergence_comparison_multichain_three_specifications_tex",
    "table", "three_specifications",
    file.path(
      sampler_sensitivity_dir, "tables",
      "sampler_convergence_comparison_multichain_three_specifications.tex"
    ),
    file.path(
      "Tables",
      "sampler_convergence_comparison_multichain_three_specifications.tex"
    ),
    "sampler_sensitivity_multichain.R", "supplement"
  ),
  artifact_row(
    "sampler_convergence_comparison_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      sampler_sensitivity_dir, "tables",
      "sampler_convergence_comparison_multichain_three_specifications.csv"
    ),
    NA_character_, "sampler_sensitivity_multichain.R", "source_data"
  ),
  artifact_row(
    "sampler_boundary_robustness_multichain_three_specifications_tex",
    "table", "three_specifications",
    file.path(
      sampler_sensitivity_dir, "tables",
      "sampler_boundary_robustness_multichain_three_specifications.tex"
    ),
    file.path(
      "Tables",
      "sampler_boundary_robustness_multichain_three_specifications.tex"
    ),
    "sampler_sensitivity_multichain.R", "supplement"
  ),
  artifact_row(
    "sampler_boundary_robustness_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      sampler_sensitivity_dir, "tables",
      "sampler_boundary_robustness_multichain_three_specifications.csv"
    ),
    NA_character_, "sampler_sensitivity_multichain.R", "source_data"
  ),
  artifact_row(
    "sampler_fitted_log_risk_agreement_multichain_three_specifications_tex",
    "table", "three_specifications",
    file.path(
      sampler_sensitivity_dir, "tables",
      "sampler_fitted_log_risk_agreement_multichain_three_specifications.tex"
    ),
    file.path(
      "Tables",
      "sampler_fitted_log_risk_agreement_multichain_three_specifications.tex"
    ),
    "sampler_sensitivity_multichain.R", "supplement"
  ),
  artifact_row(
    "sampler_fitted_log_risk_agreement_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      sampler_sensitivity_dir, "tables",
      "sampler_fitted_log_risk_agreement_multichain_three_specifications.csv"
    ),
    NA_character_, "sampler_sensitivity_multichain.R", "source_data"
  ),
  artifact_row(
    "sampler_boundary_probability_agreement_multichain_three_specifications",
    "figure", "three_specifications",
    file.path(
      sampler_sensitivity_dir, "figures",
      "sampler_boundary_probability_agreement_multichain_three_specifications.pdf"
    ),
    file.path(
      "Images",
      "sampler_boundary_probability_agreement_multichain_three_specifications.pdf"
    ),
    "sampler_sensitivity_multichain.R", "supplement"
  ),
  artifact_row(
    "sampler_fitted_log_risk_agreement_multichain_three_specifications",
    "figure", "three_specifications",
    file.path(
      sampler_sensitivity_dir, "figures",
      "sampler_fitted_log_risk_agreement_multichain_three_specifications.pdf"
    ),
    file.path(
      "Images",
      "sampler_fitted_log_risk_agreement_multichain_three_specifications.pdf"
    ),
    "sampler_sensitivity_multichain.R", "supplement"
  )
)

ppc_artifacts <- rbind(
  artifact_row(
    "ppc_multichain_three_specifications", "figure", "three_specifications",
    file.path(
      ppc_dir, "figures", "ppc_multichain_three_specifications.pdf"
    ),
    file.path("Images", "ppc_multichain_three_specifications.pdf"),
    "posterior_predictive_checks_multichain.R", "supplement"
  ),
  artifact_row(
    "observed_vs_fitted_ppc_multichain_three_specifications", "figure",
    "three_specifications",
    file.path(
      ppc_dir, "figures",
      "observed_vs_fitted_ppc_multichain_three_specifications.pdf"
    ),
    file.path(
      "Images", "observed_vs_fitted_ppc_multichain_three_specifications.pdf"
    ),
    "posterior_predictive_checks_multichain.R", "supplement"
  ),
  artifact_row(
    "observation_ppc_summary_multichain_three_specifications_tex", "table",
    "three_specifications",
    file.path(
      ppc_dir, "tables",
      "observation_ppc_summary_multichain_three_specifications.tex"
    ),
    file.path(
      "Tables", "observation_ppc_summary_multichain_three_specifications.tex"
    ),
    "posterior_predictive_checks_multichain.R", "main_and_supplement"
  ),
  artifact_row(
    "observation_ppc_summary_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      ppc_dir, "tables",
      "observation_ppc_summary_multichain_three_specifications.csv"
    ),
    NA_character_, "posterior_predictive_checks_multichain.R", "source_data"
  ),
  artifact_row(
    "unusual_observations_ppc_multichain_three_specifications_tex", "table",
    "three_specifications",
    file.path(
      ppc_dir, "tables",
      "unusual_observations_ppc_multichain_three_specifications.tex"
    ),
    file.path(
      "Tables", "unusual_observations_ppc_multichain_three_specifications.tex"
    ),
    "posterior_predictive_checks_multichain.R", "supplement"
  ),
  artifact_row(
    "unusual_observations_ppc_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      ppc_dir, "tables",
      "unusual_observations_ppc_multichain_three_specifications.csv"
    ),
    NA_character_, "posterior_predictive_checks_multichain.R", "source_data"
  ),
  artifact_row(
    "global_ppc_multichain_three_specifications_csv", "table_source",
    "three_specifications",
    file.path(
      ppc_dir, "diagnostics", "global_ppc_multichain_three_specifications.csv"
    ),
    NA_character_, "posterior_predictive_checks_multichain.R", "source_data"
  ),
  artifact_row(
    "observation_level_ppc_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      ppc_dir, "diagnostics",
      "observation_level_ppc_multichain_three_specifications.csv"
    ),
    NA_character_, "posterior_predictive_checks_multichain.R", "source_data"
  )
)

artifacts <- rbind(
  artifacts,
  sensitivity_artifacts,
  sampler_sensitivity_artifacts,
  ppc_artifacts
)
artifacts$source_file <- vapply(
  artifacts$source_file,
  function(path) {
    if (!file.exists(path)) return(path)
    normalizePath(path, winslash = "/", mustWork = TRUE)
  },
  character(1L)
)
artifacts$status <- ifelse(
  file.exists(artifacts$source_file),
  "available",
  ifelse(artifacts$required_now, "missing", "pending_final_run")
)

manifest_file <- file.path(
  validation_dir, "seer_manuscript_artifact_manifest.csv"
)
utils::write.csv(artifacts, manifest_file, row.names = FALSE, na = "")

missing_required <- artifacts[
  artifacts$required_now & artifacts$status != "available",
]
if (nrow(missing_required) > 0L) {
  stop(
    "Missing required SEER manuscript artifacts:\n",
    paste(" -", missing_required$source_file, collapse = "\n"),
    call. = FALSE
  )
}

cat("SEER manuscript artifact validation passed.\n")
cat("Available:", sum(artifacts$status == "available"), "\n")
cat("Pending final runs:", sum(artifacts$status == "pending_final_run"), "\n")
cat("Manifest:", normalizePath(manifest_file, winslash = "/"), "\n")
