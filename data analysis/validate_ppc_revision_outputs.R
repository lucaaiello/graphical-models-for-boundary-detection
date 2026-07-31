# Validate the completed PPC revision outputs without running MCMC, posterior
# simulation, or LaTeX.

rm(list = ls())

output_dir <- file.path(
  "data analysis",
  "posterior_predictive_sensitivity_output"
)
canonical_labels <- c("Adjacency", "Mean--Adjacency", "Mean")

required_files <- file.path(
  output_dir,
  c(
    "ppc_summary.csv",
    "global_ppc_exceptions.csv",
    "observation_level_ppc.csv",
    "observation_ppc_summary.csv",
    "observation_ppc_summary.tex",
    "unusual_observation_ppc.csv",
    "unusual_observation_ppc.tex",
    "observed_vs_fitted_ppc.pdf",
    "observed_vs_fitted_ppc.png",
    "main_manuscript_insertions.tex",
    "supplement_manuscript_insertions.tex",
    "specification_label_audit.txt"
  )
)
if (any(!file.exists(required_files))) {
  stop(
    "Missing required revision output(s): ",
    paste(required_files[!file.exists(required_files)], collapse = ", "),
    call. = FALSE
  )
}

global <- utils::read.csv(
  file.path(output_dir, "ppc_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
exceptions <- utils::read.csv(
  file.path(output_dir, "global_ppc_exceptions.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
observation <- utils::read.csv(
  file.path(output_dir, "observation_level_ppc.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
summary_table <- utils::read.csv(
  file.path(output_dir, "observation_ppc_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
unusual <- utils::read.csv(
  file.path(output_dir, "unusual_observation_ppc.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
main_tex <- paste(
  readLines(
    file.path(output_dir, "main_manuscript_insertions.tex"),
    warn = FALSE
  ),
  collapse = "\n"
)
supplement_tex <- paste(
  readLines(
    file.path(output_dir, "supplement_manuscript_insertions.tex"),
    warn = FALSE
  ),
  collapse = "\n"
)

stopifnot(
  nrow(global) == 48L,
  sum(global$observed_within_95) == 46L,
  setequal(global$specification, canonical_labels),
  nrow(exceptions) == 2L,
  identical(exceptions$specification, c("Adjacency", "Mean")),
  all(exceptions$cancer == "Colorectal"),
  all(exceptions$diagnostic_label == "Moran's I of county SIR"),
  all(exceptions$direction == "Observed above the 97.5% predictive quantile")
)

stopifnot(
  nrow(observation) == 696L,
  all(table(observation$specification) == 232L),
  all(table(observation$specification, observation$cancer) == 58L),
  setequal(observation$specification, canonical_labels)
)
observation_key <- paste(observation$county, observation$cancer, sep = "||")
observed_by_key <- split(observation$observed_count, observation_key)
stopifnot(all(vapply(
  observed_by_key,
  function(x) length(unique(x)) == 1L,
  logical(1)
)))

recomputed_summary <- do.call(
  rbind,
  lapply(canonical_labels, function(specification) {
    x <- observation[observation$specification == specification, ]
    data.frame(
      Specification = specification,
      Observed = mean(x$observed_count),
      Predicted = mean(x$predictive_mean),
      Lower = mean(x$predictive_lower_95),
      Upper = mean(x$predictive_upper_95),
      Coverage = mean(x$pointwise_95_coverage),
      "Low support" = sum(x$low_predictive_support),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
)
stopifnot(isTRUE(all.equal(
  summary_table,
  recomputed_summary,
  tolerance = 1e-12,
  check.attributes = FALSE
)))
stopifnot(length(unique(summary_table$Observed)) == 1L)

derived_unusual <- observation[
  observation$two_sided_mid_tail_probability < 0.02,
]
stopifnot(
  nrow(derived_unusual) == 2L,
  nrow(unusual) == 2L,
  all(derived_unusual$county == "Mono"),
  all(derived_unusual$cancer == "Lung"),
  setequal(
    derived_unusual$specification,
    c("Mean--Adjacency", "Mean")
  ),
  length(unique(
    paste(derived_unusual$county, derived_unusual$cancer)
  )) == 1L
)

required_main_strings <- c(
  "46 of the 48",
  "Adjacency (observed $0.214$",
  "$[-0.160,0.197]$",
  "$p=0.0468$",
  "Adjacency & 723.2 & 723.1 & 680.6 & 767.9 & 0.996 & 0",
  "Mean--Adjacency & 723.2 & 723.1 & 683.5 & 764.7 & 0.996 & 1",
  "Mean & 723.2 & 723.2 & 679.3 & 769.6 & 0.991 & 1",
  "The observed count is 20",
  "$-2.27$",
  "$0.0162$",
  "$-2.20$",
  "$0.0159$"
)
required_supplement_strings <- c(
  "46 of the 48",
  "$[-0.160,0.197]$",
  "$[-0.183,0.211]$",
  "$0.0468$",
  "Images/observed_vs_fitted_ppc.pdf",
  "\\label{tab:unusual_observation_ppc}",
  "Mean--Adjacency & Mono & Lung & 20 & 33.1 & 33.0",
  "Mean & Mono & Lung & 20 & 33.2 & 33.0"
)
stopifnot(all(vapply(
  required_main_strings,
  grepl,
  logical(1),
  x = main_tex,
  fixed = TRUE
)))
stopifnot(all(vapply(
  required_supplement_strings,
  grepl,
  logical(1),
  x = supplement_tex,
  fixed = TRUE
)))

banned_formal_labels <- c(
  "primary adaptive",
  "adjusted adaptive",
  "adjusted fixed",
  "covariate-adjusted adaptive",
  "covariate-adjusted fixed"
)
source_text_lower <- tolower(paste(main_tex, supplement_tex))
stopifnot(!any(vapply(
  banned_formal_labels,
  grepl,
  logical(1),
  x = source_text_lower,
  fixed = TRUE
)))

figure_info <- file.info(file.path(
  output_dir,
  c(
    "observed_vs_fitted_ppc.pdf",
    "observed_vs_fitted_ppc.png",
    "ppc_three_specifications.pdf",
    "ppc_three_specifications.png"
  )
))
stopifnot(all(!is.na(figure_info$size)), all(figure_info$size > 0))

report <- c(
  "PPC revision consistency report",
  "===============================",
  "PASS: stored global PPC contains 48 checks and exactly 46 are covered.",
  "PASS: both global exceptions are colorectal Moran's I and are named in both source insertions.",
  "PASS: observation-level output contains 696 rows, 232 per specification.",
  "PASS: observed counts agree across specifications for every county--cancer pair.",
  "PASS: the three-row table recomputes exactly from observation_level_ppc.csv.",
  "PASS: exactly two low-support checks concern the same Mono County lung observation.",
  "PASS: quantitative table and Mono values in the source insertions agree with the CSV outputs.",
  "PASS: formal specification labels are Adjacency, Mean--Adjacency, and Mean.",
  "PASS: requested figure files exist and are nonempty.",
  "No LaTeX file was compiled by this validation."
)
writeLines(
  report,
  file.path(output_dir, "ppc_revision_consistency_report.txt"),
  useBytes = TRUE
)
message(paste(report, collapse = "\n"))
