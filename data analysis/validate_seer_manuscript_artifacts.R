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
  numerical_source_stems <- c(
    "diagnostic_overview",
    "posterior_summary_core",
    "posterior_summary_dependence",
    "pooled_boundary_counts",
    "pooled_random_effect_means",
    "pooled_shared_boundary_counts",
    "pooled_mutual_boundary_counts",
    "selected_boundaries_pooled",
    "fdr_curves"
  )
  rows <- c(
    rows,
    lapply(numerical_source_stems, function(stem) {
      artifact_row(
        paste0(stem, "_", specification),
        "table_source", specification,
        file.path(root, "diagnostics", paste0(stem, "_", specification, ".csv")),
        NA_character_, "postprocess_multiple_chains_tempering.R",
        "source_data", required_now
      )
    })
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
  ),
  artifact_row(
    "figure_4_disease_graph_schematics", "figure", "conceptual_model",
    file.path(exploratory_dir, "figure_4_disease_graph_schematics.pdf"),
    file.path("Images", "figure_4_disease_graph_schematics.pdf"),
    "disease_graph_schematics.R", "main"
  ),
  artifact_row(
    "named_map_multichain_california_counties", "figure",
    "observed_data",
    file.path(
      exploratory_dir, "named_map_multichain_california_counties.pdf"
    ),
    file.path("Images", "named_map_multichain_california_counties.pdf"),
    "exploratory_figures_seer.R", "supplement"
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
  ),
  artifact_row(
    "edge_boundary_probabilities_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      sensitivity_dir, "tables",
      "edge_boundary_probabilities_multichain_three_specifications.csv"
    ),
    NA_character_, "sensitivity_analysis_multichain.R", "source_data"
  ),
  artifact_row(
    "boundary_overlap_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      sensitivity_dir, "tables",
      "boundary_overlap_multichain_three_specifications.csv"
    ),
    NA_character_, "sensitivity_analysis_multichain.R", "source_data"
  ),
  artifact_row(
    "boundary_selection_classification_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      sensitivity_dir, "tables",
      "boundary_selection_classification_multichain_three_specifications.csv"
    ),
    NA_character_, "sensitivity_analysis_multichain.R", "source_data"
  ),
  artifact_row(
    "edge_changes_meanadj_vs_mean_multichain_csv",
    "table_source", "meanadj_vs_mean",
    file.path(
      sensitivity_dir, "tables", "edge_changes_meanadj_vs_mean_multichain.csv"
    ),
    NA_character_, "sensitivity_analysis_multichain.R", "source_data"
  ),
  artifact_row(
    "fitted_relative_risk_means_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      sensitivity_dir, "tables",
      "fitted_relative_risk_means_multichain_three_specifications.csv"
    ),
    NA_character_, "sensitivity_analysis_multichain.R", "source_data"
  ),
  artifact_row(
    "fitted_risk_sensitivity_multichain_three_specifications_csv",
    "table_source", "three_specifications",
    file.path(
      sensitivity_dir, "tables",
      "fitted_risk_sensitivity_multichain_three_specifications.csv"
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

artifact_code_marker <- function(artifact, generator) {
  if (generator == "disease_graph_schematics.R") {
    return("# RESULT: Main Figure 4")
  }
  if (generator == "exploratory_figures_seer.R") {
    return(switch(
      artifact,
      figure_1_seer_cancer_sir_maps = "# RESULT: Main Figure 1",
      figure_2_seer_morans_i_by_distance_band = "# RESULT: Main Figure 2",
      figure_3_seer_covariate_maps = "# RESULT: Main Figure 3",
      named_map_multichain_california_counties =
        "# RESULT: Supplementary Figure S7"
    ))
  }
  if (generator == "postprocess_multiple_chains_tempering.R") {
    if (grepl("^pooled_fdr_boundary_maps_", artifact)) {
      return("# RESULT: Main Figure 5 / Supplementary Figures S13 and S22")
    }
    if (grepl("^shared_fdr_boundary_maps_", artifact)) {
      return("# RESULT: Main Figure 6 / Supplementary Figures S14 and S23")
    }
    if (grepl("^mutual_fdr_boundary_maps_", artifact)) {
      return("# RESULT: Main Figure 7 / Supplementary Figures S15 and S24")
    }
    if (grepl("^pooled_non_adjacency_maps_", artifact)) {
      return("# RESULT: Main Figure 8 / Supplementary Figure S16")
    }
    if (grepl("^pooled_fdr_curves_", artifact)) {
      return("# RESULT: Supplementary Figures S2, S12, and S21")
    }
    if (grepl("^pooled_beta_theta_", artifact)) {
      return("# RESULT: Supplementary Figures S3, S8, and S17")
    }
    if (grepl("^pooled_gamma_", artifact)) {
      return("# RESULT: Supplementary Figures S4, S9, and S18")
    }
    if (grepl("^pooled_V_rho_", artifact)) {
      return("# RESULT: Supplementary Figures S5, S10, and S19")
    }
    if (grepl("^pooled_(eta_)?covariance_", artifact)) {
      return("# RESULT: Supplementary Figures S6, S11, and S20")
    }
    if (grepl("^pooled_table_S(16|17|18)_", artifact)) {
      return("# RESULT: Supplementary Tables S16--S18")
    }
    return("# RESULT SOURCES: manuscript numerical summaries and diagnostics")
  }
  if (generator == "sensitivity_analysis_multichain.R") {
    if (grepl("^boundary_probability_agreement_", artifact)) {
      return("# RESULT: boundary-probability agreement figure")
    }
    if (grepl("^boundary_maps_", artifact)) {
      return("# RESULT: Mean--Adjacency versus Mean FDR boundary map")
    }
    if (grepl("^fitted_", artifact)) {
      return("# RESULT: fitted relative-risk sensitivity")
    }
    if (grepl("^edge_boundary_probabilities_", artifact)) {
      return("# RESULT: complete edge-level probabilities and FDR selections")
    }
    return("# RESULT: main sensitivity table")
  }
  if (generator == "sampler_sensitivity_multichain.R") {
    if (grepl("convergence", artifact)) {
      return("# RESULT: Supplementary sampler convergence table")
    }
    if (grepl("fitted_log_risk", artifact)) {
      return("# RESULT: Supplementary fitted log-risk agreement")
    }
    if (grepl("boundary_probability_agreement", artifact)) {
      return("# RESULT: Supplementary sampler boundary-probability agreement figure")
    }
    return("# RESULT: Supplementary sampler boundary-robustness table")
  }
  if (generator == "posterior_predictive_checks_multichain.R") {
    if (grepl("^ppc_multichain", artifact) || grepl("^global_ppc", artifact)) {
      return("# RESULT: global PPC figure")
    }
    if (grepl("^observed_vs_fitted", artifact)) {
      return("# RESULT: observed-versus-fitted PPC figure")
    }
    if (grepl("^observation_ppc_summary", artifact)) {
      return("# RESULT: observation-level PPC summary table")
    }
    if (grepl("^unusual_observations", artifact)) {
      return("# RESULT: low-predictive-support observations table")
    }
    return("# RESULT: complete machine-readable global and observation-level PPC tables")
  }
  NA_character_
}

artifacts$code_marker <- mapply(
  artifact_code_marker, artifacts$artifact, artifacts$generator,
  USE.NAMES = FALSE
)
artifacts$generator_file <- file.path("data analysis", artifacts$generator)
artifacts$generator_line <- mapply(
  function(generator_file, code_marker) {
    if (!file.exists(generator_file) || is.na(code_marker)) return(NA_integer_)
    hits <- grep(code_marker, readLines(generator_file, warn = FALSE), fixed = TRUE)
    if (!length(hits)) return(NA_integer_)
    hits[[1L]]
  },
  artifacts$generator_file, artifacts$code_marker,
  USE.NAMES = FALSE
)

missing_code_markers <- artifacts[
  is.na(artifacts$code_marker) | is.na(artifacts$generator_line),
]
if (nrow(missing_code_markers) > 0L) {
  stop(
    "Missing code marker for registered SEER artifact(s):\n",
    paste(
      " -", missing_code_markers$artifact,
      "in", missing_code_markers$generator,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

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
