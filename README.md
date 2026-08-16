# Detecting Spatial Health Disparities Using Disease Maps

This repository contains all code used in the paper  
**“Detecting Spatial Health Disparities Using Disease Maps.”**

The repository is organized into two main parts: data analysis and simulation studies.

---

## Installation and setup

Install the current CRAN/R package dependencies from the repository root with:

```r
source("install.R")
```

Open `graphical-models-for-boundary-detection.Rproj` in RStudio before running
the analyses. The project sets the repository root as the working directory, and
all scripts use paths relative to that root. From a terminal, run scripts from
the repository root, for example:

```sh
Rscript "data analysis/exploratory_figures_seer.R"
Rscript "data analysis/disease_graph_schematics.R"
```

The sampler source files are compiled with `Rcpp::sourceCpp()`, so a working C/C++
toolchain is required. On Windows this usually means Rtools; on macOS this usually
means Xcode command line tools; on Linux this usually means the system compiler
toolchain.

The C++ sampler files use `RcppArmadillo`.

## User-facing scripts and functions

This repository is organized as a collection of reproducibility scripts rather
than as an installed R package. The recommended user interface is therefore to
run the R scripts from the repository root. Functions defined inside those
scripts, such as `sd_diff_mat()` and `bd()`, are internal helpers and are not
intended as a stable public API.

### Script inputs and outputs

| Script type | Required inputs | Anticipated outputs |
|---|---|---|
| [`install.R`](install.R) | An R installation, internet access, and a C/C++ toolchain | Installs the required R packages |
| [`data analysis/exploratory_figures_seer.R`](<data analysis/exploratory_figures_seer.R>) | `SIR_adjusted.csv`, `covariates.csv`, and `smoking.csv` in `data analysis/data/` | Creates the observed-data PDFs used as main Figures 1--3 and the California county-name map used as Figure S7 under `data analysis/exploratory_figures/` |
| [`data analysis/disease_graph_schematics.R`](<data analysis/disease_graph_schematics.R>) | No data inputs | Deterministically creates the unstructured, directed-acyclic, and undirected disease-graph panels used as main Figure 4 |
| [`data analysis/main_multiple_chains_tempering.R`](<data analysis/main_multiple_chains_tempering.R>) and [`sampler_multiple_chains_tempering.cpp`](<data analysis/sampler_multiple_chains_tempering.cpp>) | The SEER data and one of `adj`, `meanadj`, or `mean` | Runs the reported six-chain adaptive-tempering analysis and saves each completed cold posterior stream |
| [`data analysis/postprocess_multiple_chains_tempering.R`](<data analysis/postprocess_multiple_chains_tempering.R>) | A completed adaptive-tempering run | Writes specification-specific diagnostics, posterior tables, FDR summaries, and manuscript figures |
| [`data analysis/main_multiple_chains.R`](<data analysis/main_multiple_chains.R>) and [`sampler_multiple_chains.cpp`](<data analysis/sampler_multiple_chains.cpp>) | The same SEER inputs and specification | Runs the independent non-tempered six-chain replication used only for sampler sensitivity |
| [`data analysis/postprocess_multiple_chains.R`](<data analysis/postprocess_multiple_chains.R>) | A completed non-tempered run | Writes the matching diagnostics and posterior summaries automatically under that run's `postprocessing/` directory |
| [`data analysis/sensitivity_analysis_multichain.R`](<data analysis/sensitivity_analysis_multichain.R>) | The three post-processed adaptive-tempering fits | Writes the scientific model-specification sensitivity tables and figures under `data analysis/sensitivity/` |
| [`data analysis/posterior_predictive_checks_multichain.R`](<data analysis/posterior_predictive_checks_multichain.R>) | The three completed adaptive-tempering fits | Writes global and observation-level PPC tables and figures under `data analysis/ppc/` |
| [`data analysis/sampler_sensitivity_multichain.R`](<data analysis/sampler_sensitivity_multichain.R>) | All three post-processed adaptive-tempering and non-tempered fits | Writes the transition-kernel convergence, boundary, and fitted-risk comparisons under `data analysis/sampler_sensitivity/` |
| [`data analysis/validate_seer_manuscript_artifacts.R`](<data analysis/validate_seer_manuscript_artifacts.R>) | All completed post-processing outputs | Verifies every SEER table and figure referenced by the revised manuscript and writes an artifact manifest |
| `sim_*.R` simulation generators | The scenario settings embedded in the script, the corresponding C++ sampler, and an existing matching `runs_*` directory | One `mcmc_samples_<seed>.rds` file per replicate, plus the simulated truth and design objects saved in the matching `RE_generation_*` directory |
| `assessment*.R` | The matching `RE_generation_*` objects and `mcmc_samples_<seed>.rds` files | Prints sensitivity, specificity, and adjacency-recovery summaries; non-FDR assessment scripts also save one WAIC vector |
| `assessment*_FDR.R` | The matching generated truth and saved MCMC replicates | Prints FDR-based sensitivity and specificity summaries |
| `rmse*.R` | The matching generated truth and saved MCMC replicates | Saves one compact RMSE vector in `Misspecified models/rmse/` |
| [`WAIC_assessment.R`](<Misspecified models/WAIC/WAIC_assessment.R>) and [`rmse_assessment.R`](<Misspecified models/rmse/rmse_assessment.R>) | The compact WAIC or RMSE `.rds` files created by the preceding scripts | Prints the final table values |
| WAIC comparison scripts | The corresponding WAIC vectors | Draws the WAIC density plots used in Figure S1 |
| Comparison folders | Their `sim_*.R` or comparison generator followed by the matching `assessment.R` | Produces the Lee & Mitchell and Li et al. comparison summaries |

The two standalone figure scripts save Figures 1--4 and Figure S7 directly. The reported
model fits and every downstream table and figure are likewise saved in their
workflow-specific output directories; those generated outputs are excluded
from Git and are recreated by the scripts above.

### Low-level sampler interface

The C++ samplers are exposed to R by `Rcpp::sourceCpp()`, but the supported
interface is through the supplied driver scripts because they construct and
validate the county ordering, disease stacking, design matrices, adjacency
matrices, seeds, and initialization regimes. The two empirical entry points are:

- `MADAGAR_tempered_segment()` in
  `sampler_multiple_chains_tempering.cpp`, called by
  `main_multiple_chains_tempering.R` for the reported analysis;
- `MADAGAR_chain()` in `sampler_multiple_chains.cpp`, called
  by `main_multiple_chains.R` for the non-tempered sampler-sensitivity
  replication.

Both return posterior draws for the regression, spatial, mixture, dependence,
and adjacency components together with acceptance, initialization, runtime, and
state-consistency information used by the matching post-processing scripts.
Users should set the documented environment variables and run the R drivers
rather than call these low-level functions directly.

The simulation samplers use the same basic data layout but a Gaussian outcome
model. Their public function names are:

- `MADAGAR()` for the multivariate DAGAR and CAR simulation scripts;
- `ADAGAR()` for the univariate Lee--Mitchell comparison;
- `MCAR_indep()` for the independent-CAR comparison.

Simulation calls use one dissimilarity matrix, `Z1`. Directed and undirected
disease-graph samplers additionally require `W_dis`, the binary `q x q`
inter-disease graph. The simulation `MADAGAR()` and `ADAGAR()` functions return
`beta`, `taud`, `phi`, `theta`, `u`, `rho`, `v`, `r`, `F_r`, `eta`, `tau`,
and `W1`--`W4`, followed by `A`, the directed-graph coefficients, or the
undirected disease-correlation parameter as appropriate. `MCAR_indep()` returns
`beta`, `taud`, `phi`, `theta`, `u`, `rho`, `v`, `r`, `F_r`, `tau`, and
`W1`--`W4`.

### Runtime notes

The installation script only installs packages; it does not run the analyses.

The following approximate wall-clock times are intended to distinguish quick
post-processing tasks from substantial model fitting. They are not formal
benchmarks. The reference machine used for the timings reported in the
supplement was an Intel i7-10750H CPU at 2.60 GHz. Runtime can differ
substantially with processor speed, available memory, compiler, operating
system, and whether replicates are run in parallel. The model-fitting estimates
below are calculated from the per-iteration timings reported in the supplement
and the iteration counts in the scripts; compilation, file writing, and
post-processing add some overhead.

| Task | Approximate runtime | Notes |
|---|---:|---|
| Install R dependencies with `install.R` | 5--30 minutes | Usually required only once; network speed and compilation dominate |
| Compile one C++ sampler | Less than 1--3 minutes | Performed by `Rcpp::sourceCpp()` when a fitting script starts |
| Generate observed-data Figures 1--3 and Figure S7 | Usually less than 1 minute after the Census geometry is cached | The first run may require downloading the fixed 2016 county geometry |
| Post-process one completed empirical fit | Seconds to a few minutes | Run the matching standalone post-processing script after all six chains finish |
| One six-chain adaptive-tempering SEER specification | Approximately 15--20 hours on the reference machine | Six parallel workers, 16 temperatures per worker, 50,000 warm-up iterations and 25,000 retained cold draws per worker |
| One six-chain non-tempered SEER replication | Approximately 1.3--1.8 hours on the reference machine | Six parallel workers, 300,000 warm-up and 500,000 sampling iterations per worker, retaining every tenth draw |
| One 100-replicate unstructured MDAGAR simulation script | Approximately 3.3 hours | 0.006 seconds x 20,000 iterations x 100 replicates |
| One 100-replicate undirected MDAGAR simulation script | Approximately 5.0 hours | 0.009 seconds x 20,000 iterations x 100 replicates |
| One 100-replicate directed MDAGAR simulation script | Approximately 7.2 hours | 0.013 seconds x 20,000 iterations x 100 replicates |
| Complete 3 x 3 MDAGAR simulation grid | Approximately 46.7 hours of sampler time | Three scripts for each disease-graph specification; plan for roughly 2--3 days sequentially after allowing for overhead |
| One 100-replicate multivariate CAR simulation script | Approximately 23.9--27.2 hours | 0.043--0.049 seconds x 20,000 iterations x 100 replicates |
| Complete three-script CAR simulation set | Approximately 71.7--81.7 hours | Unstructured, directed, and undirected disease-graph specifications |
| Lee--Mitchell univariate comparison | Approximately 2--3 hours | Includes 100 DAGAR and 100 CARBayes fits |
| Li et al. independent-CAR comparison | Approximately 23.9--27.2 hours | 100 multivariate independent-CAR fits |
| Proposed-method side of the Li comparison | Approximately 3--8 hours | 100 multivariate DAGAR fits |
| Assessment or RMSE script over 100 saved runs | Tens of minutes to several hours | Depends on whether it computes boundary probabilities, adjacency recovery, or WAIC |
| Final WAIC/RMSE aggregation and density plots | Seconds to a few minutes | Assumes the small summary `.rds` files have already been generated |

Most simulation-generating scripts use `burn = 10000`, `runs = 10000`, and
`thin = 1` for each of 100 seeds, corresponding to 20,000 C++ iterations per
replicate. The simulation scripts call `tic()` and `toc()` around each fit, so
users can obtain machine-specific per-replicate timings from the console before
committing to all 100 replicates.

### Simulation output and assessment workflow

Each simulation-generating script performs the model fitting for a sequence of
simulation seeds. For every completed seed, it saves the corresponding MCMC
output as a separate file named `mcmc_samples_<seed>.rds` in the `runs_*`
directory belonging to that simulation scenario and fitted model. For example:

```text
sim_gaussian_DAGAR.R
    -> runs_DAGAR/mcmc_samples_1.rds
    -> runs_DAGAR/mcmc_samples_2.rds
    ...
    -> runs_DAGAR/mcmc_samples_100.rds
```

The generator and its output directory must remain paired. Depending on the
model, the output directories are named `runs_DAGAR`,
`runs_directed_DAGAR`, `runs_undirected_DAGAR`, `runs_CAR`,
`runs_directed_CAR`, `runs_undirected_CAR`, `runs_CARBayes`, or
`runs_CAR_indep`. The scripts expect these directories to exist before the
simulation begins.

The assessment scripts do not refit the models. Instead, they iterate over the
saved `mcmc_samples_<seed>.rds` files in the matching `runs_*` directory and
combine results across replicates to calculate boundary-detection sensitivity
and specificity, adjacency recovery, WAIC, or other reported summaries. The
scenario-specific `rmse*.R` scripts similarly read the saved MCMC files and
write compact RMSE summaries to `Misspecified models/rmse/`. The final
aggregation scripts then read the compact WAIC or RMSE summary files to produce
the published tables and plots.

Consequently, the workflow for a fresh simulation is:

1. Create the required `runs_*` output directory if it is absent.
2. Run the corresponding `sim_*.R` generator and allow all required seeds to
   finish.
3. Confirm that the expected `mcmc_samples_<seed>.rds` files are present.
4. Run the matching assessment and, where applicable, RMSE script.
5. Run the WAIC/RMSE aggregation or plotting script.

The complete set of existing `runs_*` outputs contains 1,600 `.rds` files and
occupies approximately 63.4 GB. For a quick test, reduce the seed loop and MCMC
settings before launching the full analysis; such reduced runs are for checking
the software workflow only and do not reproduce the published numerical
results.

The scripts do not call `setwd()` or contain machine-specific absolute paths.
Input, output, and C++ source paths are relative to the repository root.

---

## Data Analysis

The `data analysis/` folder contains the data and scripts used to produce the results for the empirical analysis settings considered in the paper.

### Data

The `data/` subfolder includes the following files:

- **`SIR_adjusted.csv`**  
  Observed and expected counts for four cancers, together with standardized incidence ratios (SIRs), for each California county.

- **`covariates.csv`**  
  County-level demographic and socioeconomic covariates for U.S. counties, including:
  - percentage of population under 18,
  - percentage over 65,
  - percentage with high school education,
  - percentage of families below the poverty threshold,
  - unemployment rate.

- **`insurance.csv`**  
  County-level percentages of insurance status (e.g., uninsured, Medicaid, insured).

- **`race.csv`**  
  County-level race composition percentages (White, Black, Other).

- **`sex.csv`**  
  County-level sex composition percentages (Male, Female).

- **`smoking.csv`**  
  Smoking prevalence for each California county.

### Model Fitting

For all three empirical settings, an **unstructured disease graph** is used.

- **`main_multiple_chains_tempering.R`** prepares the SEER inputs and
  runs six independently initialized adaptive-tempering samplers.

- **`sampler_multiple_chains_tempering.cpp`** implements the
  Metropolis-within-Gibbs transition kernel and temperature exchanges in
  **RcppArmadillo**.

Set `MCMC_TEMPERING_SPECIFICATION` to control the role of covariates:

- `"adj"`: covariates enter only the adjacency (boundary) model;
- `"mean"`: covariates enter only the mean structure;
- `"meanadj"`: covariates enter both the mean and adjacency components.

The `"mean"` setting keeps the geographic adjacency matrices fixed throughout
the MCMC run. It is therefore the fixed-adjacency MARDP/MDAGAR special case most
closely corresponding to Gao et al. (2023), whereas `"adj"` and `"meanadj"`
activate the adjacency-modeling extension developed in this paper.

To reproduce an empirical setting, set the corresponding environment variable,
run the adaptive-tempering main script, and then run its
post-processing script. The resulting setting-specific outputs are:

- `cvrts = "adj"` produces main Figures 5--8, Table S16, and Figures S2--S6;
- `cvrts = "meanadj"` produces Table S17 and Figures S8--S16;
- `cvrts = "mean"` produces Table S18 and Figures S17--S24.

Figures 1--3 and the California county-name map in Figure S7 are produced independently by
`data analysis/exploratory_figures_seer.R` from the observed data and do not
depend on the fitted model variant. The conceptual disease-graph panels in
Figure 4 are produced by `data analysis/disease_graph_schematics.R` and require
no data or posterior samples.

### Three-specification sensitivity and posterior predictive workflow

The additional SEER analyses compare three specifications while preserving the
same cancer data, county ordering, likelihood, priors, MCMC settings, and FDR
rule:

| Manuscript label | `cvrts` | Mean model | Adjacency model | Purpose in the comparison |
|---|---|---|---|---|
| Adjacency | `"adj"` | Cancer-specific intercepts | Learned using smoking, population over 65, and poverty | Primary proposed specification |
| Mean--Adjacency | `"meanadj"` | The same three covariates enter the mean | Adjacency is still learned using the covariates | Comparison with Adjacency assesses sensitivity to mean adjustment |
| Mean | `"mean"` | Same covariate-adjusted mean as Mean--Adjacency | Geographic adjacency is fixed | Comparison with Mean--Adjacency isolates the effect of adjacency learning and provides the Gao et al. comparison |

Run `main_multiple_chains_tempering.R` from the repository root for
`"adj"`, `"meanadj"`, and `"mean"`, using the reported `hotter_final`
configuration, and then run `postprocess_multiple_chains_tempering.R`
for each completed specification. For example:

```r
Sys.setenv(
  MCMC_TEMPERING_SPECIFICATION = "adj",
  MCMC_TEMPERING_MODE = "hotter_final",
  MCMC_TEMPERING_N_TEMPERATURES = "16",
  MCMC_TEMPERING_BETA_HOT = "0.01",
  MCMC_TEMPERING_ENDPOINT_STRATEGY = "fixed"
)
source("data analysis/main_multiple_chains_tempering.R")
source("data analysis/postprocess_multiple_chains_tempering.R")
```

Replace `"adj"` by `"meanadj"` or `"mean"` for the other specifications.
Results are kept under
`data analysis/multiple_chains_tempering_output/<specification>/hotter_final/`,
with all derived diagnostics, tables, and figures in its `postprocessing/`
subdirectory.

The independent non-tempered multiple-chain workflow uses the same automatic
specification/mode hierarchy. For example, the complete Mean--Adjacency run
and its post-processing are executed from the repository root with:

```r
Sys.setenv(
  MCMC_MULTICHAIN_SPECIFICATION = "meanadj",
  MCMC_MULTICHAIN_MODE = "production"
)
source("data analysis/main_multiple_chains.R")
source("data analysis/postprocess_multiple_chains.R")
```

No input or output directory needs to be supplied. The sampler writes to
`data analysis/multiple_chains_output/meanadj/production/`, and the
second script automatically reads those six completed chains and writes all
derived diagnostics, tables, figures, and summary objects to its
`postprocessing/` subdirectory. Replacing `"meanadj"` by `"adj"` or `"mean"`
selects the corresponding specification.

After the primary adaptive-tempering workflow has been post-processed for all
three specifications, generate the model-specification sensitivity and
posterior predictive results with:

```r
source("data analysis/sensitivity_analysis_multichain.R")
source("data analysis/posterior_predictive_checks_multichain.R")
```

The first script writes the boundary and fitted-risk sensitivity results under
`data analysis/sensitivity/`. The second uses all 150,000 retained cold draws
per specification by default and writes the global and observation-level
posterior predictive results under `data analysis/ppc/`.

After both the adaptive-tempering and independent non-tempered workflows have
been post-processed for all three specifications, run the transition-kernel
sensitivity comparison with:

```r
source("data analysis/sampler_sensitivity_multichain.R")
```

This script automatically compares adaptive tempering with independent
non-tempered sampling within Adjacency, Mean--Adjacency, and Mean. It writes
the combined convergence, boundary-selection, edge-probability, and fitted
log-risk agreement tables and figures under `data analysis/sampler_sensitivity/`.
This sampler comparison is intentionally separate from
`sensitivity_analysis_multichain.R`, which compares the three scientific model
specifications using the primary adaptive-tempering results. Finally, run:

```r
source("data analysis/validate_seer_manuscript_artifacts.R")
```

`validate_seer_manuscript_artifacts.R` writes
`data analysis/validation/seer_manuscript_artifact_manifest.csv`, verifies
every currently reportable SEER table, figure, and machine-readable numerical
source, and records the producing script, searchable `# RESULT:` marker, and
current source line. It stops with an explicit list if an artifact or code
marker is missing. The validation directory is created automatically.

#### Published multichain result index

The relevant blocks are marked by searchable `# RESULT:` comments. Run each
complete script because the published blocks use validated objects constructed
earlier in that file.

| Published result | Script | Searchable code marker | Generated file |
|---|---|---:|---|
| Main Figures 1--3 | [`exploratory_figures_seer.R`](<data analysis/exploratory_figures_seer.R>) | Search `# RESULT: Main Figure 1`, `2`, or `3` | `exploratory_figures/figure_1_seer_cancer_sir_maps.pdf`, `figure_2_seer_morans_i_by_distance_band.pdf`, and `figure_3_seer_covariate_maps.pdf` |
| Main Figure 4 | [`disease_graph_schematics.R`](<data analysis/disease_graph_schematics.R>) | Search `# RESULT: Main Figure 4` | `exploratory_figures/figure_4_disease_graph_schematics.pdf` |
| Main Figure 5 / Figures S13 and S22 | [`postprocess_multiple_chains_tempering.R`](<data analysis/postprocess_multiple_chains_tempering.R>) | Search `# RESULT: Main Figure 5 / Supplementary Figures S13 and S22` | `pooled_fdr_boundary_maps_<specification>.pdf` |
| Main Figure 6 / Figures S14 and S23 | Same script | Search `# RESULT: Main Figure 6 / Supplementary Figures S14 and S23` | `shared_fdr_boundary_maps_<specification>.pdf` |
| Main Figure 7 / Figures S15 and S24 | Same script | Search `# RESULT: Main Figure 7 / Supplementary Figures S15 and S24` | `mutual_fdr_boundary_maps_<specification>.pdf` |
| Main Figure 8 / Figure S16 | Same script | Search `# RESULT: Main Figure 8 / Supplementary Figure S16` | `pooled_non_adjacency_maps_adj.pdf` and `pooled_non_adjacency_maps_meanadj.pdf` |
| Figures S2, S12, and S21 | Same script | Search `# RESULT: Supplementary Figures S2, S12, and S21` | `pooled_fdr_curves_<specification>.pdf` |
| Figures S3, S8, and S17 | Same script | Search `# RESULT: Supplementary Figures S3, S8, and S17` | `pooled_beta_theta_<specification>.pdf` |
| Figures S4, S9, and S18 | Same script | Search `# RESULT: Supplementary Figures S4, S9, and S18` | `pooled_gamma_<specification>.pdf` |
| Figures S5, S10, and S19 | Same script | Search `# RESULT: Supplementary Figures S5, S10, and S19` | `pooled_V_rho_<specification>.pdf` |
| Figures S6, S11, and S20 | Same script | Search `# RESULT: Supplementary Figures S6, S11, and S20` | `pooled_eta_covariance_adj.pdf`, `pooled_eta_covariance_meanadj.pdf`, and `pooled_covariance_mean.pdf` |
| Tables S16--S18 | Same script | Search `# RESULT: Supplementary Tables S16--S18` | `pooled_table_S16_adj`, `pooled_table_S17_meanadj`, and `pooled_table_S18_mean`, each as `.csv` and `.tex` |
| Supplementary Figure S7 | [`exploratory_figures_seer.R`](<data analysis/exploratory_figures_seer.R>) | Search `# RESULT: Supplementary Figure S7` | `exploratory_figures/named_map_multichain_california_counties.pdf` and `.png` |
| Main Table 1: boundary sensitivity | [`sensitivity_analysis_multichain.R`](<data analysis/sensitivity_analysis_multichain.R>) | Search `# RESULT: main sensitivity table` | `sensitivity/tables/boundary_sensitivity_table_multichain_three_specifications.csv` and `.tex` |
| Supplementary boundary-probability agreement | [`sensitivity_analysis_multichain.R`](<data analysis/sensitivity_analysis_multichain.R>) | Search `# RESULT: boundary-probability agreement figure` | `sensitivity/figures/boundary_probability_agreement_multichain_three_specifications.pdf` and `.png` |
| Supplementary Mean--Adjacency versus Mean boundary map | [`sensitivity_analysis_multichain.R`](<data analysis/sensitivity_analysis_multichain.R>) | Search `# RESULT: Mean--Adjacency versus Mean FDR boundary map` | `sensitivity/figures/boundary_maps_multichain_meanadj_vs_mean.pdf` and `.png` |
| Supplementary global posterior predictive checks | [`posterior_predictive_checks_multichain.R`](<data analysis/posterior_predictive_checks_multichain.R>) | Search `# RESULT: global PPC figure` | `ppc/figures/ppc_multichain_three_specifications.pdf` and `.png` |
| Main Table 2: observation-level PPC summary | [`posterior_predictive_checks_multichain.R`](<data analysis/posterior_predictive_checks_multichain.R>) | Search `# RESULT: observation-level PPC summary table` | `ppc/tables/observation_ppc_summary_multichain_three_specifications.csv` and `.tex` |
| Supplementary low-support county--cancer checks | [`posterior_predictive_checks_multichain.R`](<data analysis/posterior_predictive_checks_multichain.R>) | Search `# RESULT: low-predictive-support observations table` | `ppc/tables/unusual_observations_ppc_multichain_three_specifications.csv` and `.tex` |
| Supplementary observed versus predicted counts | [`posterior_predictive_checks_multichain.R`](<data analysis/posterior_predictive_checks_multichain.R>) | Search `# RESULT: observed-versus-fitted PPC figure` | `ppc/figures/observed_vs_fitted_ppc_multichain_three_specifications.pdf` and `.png` |
| Supplementary sampler convergence comparison | [`sampler_sensitivity_multichain.R`](<data analysis/sampler_sensitivity_multichain.R>) | Search `# RESULT: Supplementary sampler convergence table` | `sampler_sensitivity/tables/sampler_convergence_comparison_multichain_three_specifications.csv` and `.tex` |
| Supplementary sampler boundary robustness | [`sampler_sensitivity_multichain.R`](<data analysis/sampler_sensitivity_multichain.R>) | Search `# RESULT: Supplementary sampler boundary-robustness table` | `sampler_sensitivity/tables/sampler_boundary_robustness_multichain_three_specifications.csv` and `.tex` |
| Supplementary fitted log-risk sampler agreement | [`sampler_sensitivity_multichain.R`](<data analysis/sampler_sensitivity_multichain.R>) | Search `# RESULT: Supplementary fitted log-risk agreement table` | `sampler_sensitivity/tables/sampler_fitted_log_risk_agreement_multichain_three_specifications.csv` and `.tex` |

---

## Simulation Studies

Simulation code is contained in the `Misspecified models/` and `CAR/` folders.

### Misspecified Models

The `Misspecified models/` folder contains all code used to produce the simulation results presented in the manuscript. The folder structure is nested as follows:

1. The top level contains three folders:
   - `Directed`
   - `Undirected`
   - `Unstructured`  
   These correspond to the **data-generating process (DGP)**.

2. Each of these folders contains three subfolders:
   - `Directed`
   - `Undirected`
   - `Unstructured`  
   These correspond to the **fitted models**, including the two misspecified models and the correctly specified one.

Each subfolder includes:
- **`RE_generation_X_DAGAR/`**  
  Generated random effects (synthetic data) for the corresponding scenario.

- **`assessment_X.R`**  
  Scripts for evaluating model performance. These scripts expect results to be stored in a folder named `runs_X/`, which must be created manually.

- **`sim_gaussian_X_DAGAR.R`**  
  Script used to generate simulation runs. The outputs of this script should be saved in the corresponding `runs_X/` folder.

### Model Comparison: WAIC and RMSE

The `WAIC/` and `rmse/` folder contains scripts to compute and visualize:
- Widely Applicable Information Criterion (WAIC),
- Root Mean Squared Error (RMSE).

These scripts generate the tables and figures reported in the manuscript.

---

## CAR Model Comparison

The `CAR/` folder contains code to reproduce the results in the **Supplementary Materials**, where the proposed approach is compared with the **MCAR (Multivariate Conditional Autoregressive)** specification.

---

## Comparison with Competing Models

The repository includes two folders, `Lee_Mitchell_comparison/` and `Multivariate_Li_comparison/`, which contain the code required to reproduce comparisons with existing boundary-detection methods proposed in:

- **Lee, D., & Mitchell, R. (2012)**  
  *Boundary detection in disease mapping studies.* Biostatistics, 13(3), 415–426.

- **Li, P., Banerjee, S., Hanson, T. A., & McBean, A. M. (2015)**  
  *Bayesian models for detecting difference boundaries in areal data.* Statistica Sinica, 25(1), 385–399.

### Multivariate Li et al. (2015) Comparison

The `Multivariate_Li_comparison/` folder is further organized into two subfolders:

- **`Li_et_al_2015/`**  
  Implements the model proposed by Li et al. (2015).

- **`Our_method/`**  
  Implements the proposed method in this paper.

These folders allow for a controlled comparison under two complementary scenarios:
1. Data are generated under the proposed model and fitted using the independent CAR model with fixed adjacency.
2. Data are generated under the independent CAR model with fixed adjacency and fitted using the proposed method.

This setup facilitates a direct assessment of model robustness and performance under correct specification and misspecification.

---

## Reproducibility: figures and tables

This mapping corresponds to the current main manuscript (Figures 1--8 and
Tables 1--2) and supplementary material (Tables S1--S19 and Figures S1--S28).

### Main manuscript

| Result | Description | Repository files |
|---|---|---|
| Figure 1 | Cancer SIR maps | [`data analysis/exploratory_figures_seer.R`](<data analysis/exploratory_figures_seer.R>); see `figure_1_seer_cancer_sir_maps.pdf` |
| Figure 2 | Moran's I by distance band | Same script; see `figure_2_seer_morans_i_by_distance_band.pdf` |
| Figure 3 | Smoking, elderly-population, and poverty maps | Same script; see `figure_3_seer_covariate_maps.pdf` |
| Figure 4 | Three disease-graph schematics | [`data analysis/disease_graph_schematics.R`](<data analysis/disease_graph_schematics.R>); see `figure_4_disease_graph_schematics.pdf` |
| Figure 5 | Disease-specific difference boundaries | The Adjacency `hotter_final` run followed by [`postprocess_multiple_chains_tempering.R`](<data analysis/postprocess_multiple_chains_tempering.R>); see `pooled_fdr_boundary_maps_adj.pdf` |
| Figure 6 | Shared difference boundaries | Same workflow; see `shared_fdr_boundary_maps_adj.pdf` |
| Figure 7 | Mutual cross-difference boundaries | Same workflow; see `mutual_fdr_boundary_maps_adj.pdf` |
| Figure 8 | Estimated non-adjacencies | Same workflow; see `pooled_non_adjacency_maps_adj.pdf` |
| Table 1 | Disease-specific boundary sensitivity across the three specifications | [`sensitivity_analysis_multichain.R`](<data analysis/sensitivity_analysis_multichain.R>); see `boundary_sensitivity_table_multichain_three_specifications.csv` and `.tex` |
| Table 2 | Observation-level posterior predictive summaries | [`posterior_predictive_checks_multichain.R`](<data analysis/posterior_predictive_checks_multichain.R>); see `observation_ppc_summary_multichain_three_specifications.csv` and `.tex` |

The reported empirical workflow is based on the six-chain
adaptive-tempering sampler. The specification-specific post-processing script
creates the primary maps, posterior intervals, FDR curves, and posterior tables.
The sensitivity and posterior predictive scripts then use the three completed
fits to generate the two main-manuscript tables and the additional
supplementary results.

### Supplementary tables

| Result | Repository files |
|---|---|
| Table S1 | [`WAIC_assessment.R`](<Misspecified models/WAIC/WAIC_assessment.R>), using the nine WAIC files produced by the nine non-FDR assessment scripts in the 3 x 3 `Misspecified models` grid |
| Table S2 | Unstructured fit: [`Unstructured/Unstructured/assessment_FDR.R`](<Misspecified models/Unstructured/Unstructured/assessment_FDR.R>), [`Directed/Unstructured/assessment_FDR.R`](<Misspecified models/Directed/Unstructured/assessment_FDR.R>), and [`Undirected/Unstructured/assessment_FDR.R`](<Misspecified models/Undirected/Unstructured/assessment_FDR.R>) |
| Table S3 | Unstructured fit: [`Unstructured/Unstructured/assessment.R`](<Misspecified models/Unstructured/Unstructured/assessment.R>), [`Directed/Unstructured/assessment.R`](<Misspecified models/Directed/Unstructured/assessment.R>), and [`Undirected/Unstructured/assessment.R`](<Misspecified models/Undirected/Unstructured/assessment.R>) |
| Table S4 | The adjacency-recovery sections of the same three non-FDR assessment scripts used for Table S3 |
| Table S5 | [`rmse_assessment.R`](<Misspecified models/rmse/rmse_assessment.R>), using the nine `rmse*.R` scripts in the 3 x 3 `Misspecified models` grid |
| Table S6 | Correctly specified DAGAR fits: [`Unstructured/Unstructured/assessment_FDR.R`](<Misspecified models/Unstructured/Unstructured/assessment_FDR.R>), [`Directed/Directed/assessment_directed_FDR.R`](<Misspecified models/Directed/Directed/assessment_directed_FDR.R>), and [`Undirected/Undirected/assessment_undirected_FDR.R`](<Misspecified models/Undirected/Undirected/assessment_undirected_FDR.R>) |
| Table S7 | Correctly specified DAGAR fits: [`Unstructured/Unstructured/assessment.R`](<Misspecified models/Unstructured/Unstructured/assessment.R>), [`Directed/Directed/assessment_directed.R`](<Misspecified models/Directed/Directed/assessment_directed.R>), and [`Undirected/Undirected/assessment_undirected.R`](<Misspecified models/Undirected/Undirected/assessment_undirected.R>) |
| Table S8 | [`CAR/assessment_FDR.R`](<CAR/assessment_FDR.R>), [`CAR/assessment_directed_FDR.R`](<CAR/assessment_directed_FDR.R>), and [`CAR/assessment_undirected_FDR.R`](<CAR/assessment_undirected_FDR.R>) |
| Table S9 | [`CAR/assessment.R`](<CAR/assessment.R>), [`CAR/assessment_directed.R`](<CAR/assessment_directed.R>), and [`CAR/assessment_undirected.R`](<CAR/assessment_undirected.R>) |
| Table S10 | The adjacency-recovery sections of the three correctly specified DAGAR assessment scripts used for Table S7 |
| Table S11 | The adjacency-recovery sections of the three CAR assessment scripts used for Table S9 |
| Tables S12--S13 | [`Lee_Mitchell_comparison/assessment.R`](<Lee_Mitchell_comparison/assessment.R>); the simulations are generated by [`comparison_univariate.R`](<Lee_Mitchell_comparison/comparison_univariate.R>) and [`sampler_sim_gaussian_DAGAR_fastest_univariate.cpp`](<Lee_Mitchell_comparison/sampler_sim_gaussian_DAGAR_fastest_univariate.cpp>) |
| Table S14 | [`Multivariate_Li_comparison/Li_et_al_2015/assessment.R`](<Multivariate_Li_comparison/Li_et_al_2015/assessment.R>); upstream files are [`sim_gaussian_DAGAR.R`](<Multivariate_Li_comparison/Li_et_al_2015/sim_gaussian_DAGAR.R>) and [`sampler_sim_gaussian_CAR_fastest_indep.cpp`](<Multivariate_Li_comparison/Li_et_al_2015/sampler_sim_gaussian_CAR_fastest_indep.cpp>) |
| Table S15 | [`Multivariate_Li_comparison/Our_method/assessment.R`](<Multivariate_Li_comparison/Our_method/assessment.R>); upstream files are [`sim_gaussian_DAGAR.R`](<Multivariate_Li_comparison/Our_method/sim_gaussian_DAGAR.R>) and [`sampler_sim_gaussian_DAGAR_fastest.cpp`](<Multivariate_Li_comparison/Our_method/sampler_sim_gaussian_DAGAR_fastest.cpp>) |
| Table S16 | Run [`postprocess_multiple_chains_tempering.R`](<data analysis/postprocess_multiple_chains_tempering.R>) for `adj`; it writes `pooled_table_S16_adj.csv` and `.tex` |
| Table S17 | Run [`postprocess_multiple_chains_tempering.R`](<data analysis/postprocess_multiple_chains_tempering.R>) for `meanadj`; it writes `pooled_table_S17_meanadj.csv` and `.tex` |
| Table S18 | Run [`postprocess_multiple_chains_tempering.R`](<data analysis/postprocess_multiple_chains_tempering.R>) for `mean`; it writes `pooled_table_S18_mean.csv` and `.tex` |
| Observation-level low-support table | Search `# RESULT: low-predictive-support observations table` in [`posterior_predictive_checks_multichain.R`](<data analysis/posterior_predictive_checks_multichain.R>); see `unusual_observations_ppc_multichain_three_specifications.csv` and `.tex` |
| Sampler convergence table | [`sampler_sensitivity_multichain.R`](<data analysis/sampler_sensitivity_multichain.R>); see `sampler_convergence_comparison_multichain_three_specifications.csv` and `.tex` |
| Sampler boundary-robustness table | Same script; see `sampler_boundary_robustness_multichain_three_specifications.csv` and `.tex` |
| Sampler fitted-log-risk agreement table | Same script; see `sampler_fitted_log_risk_agreement_multichain_three_specifications.csv` and `.tex` |

For Tables S1 and S5, every cell of the 3 x 3 true-graph/fitted-graph
comparison is required. Consequently, all nine simulation scripts, all nine
non-FDR assessment scripts, and all nine RMSE scripts under
`Misspecified models/` contribute to published results.

The complete 3 x 3 source grid is:

| True graph | Fitted unstructured | Fitted directed | Fitted undirected |
|---|---|---|---|
| Unstructured | [`simulation`](<Misspecified models/Unstructured/Unstructured/sim_gaussian_DAGAR.R>), [`assessment`](<Misspecified models/Unstructured/Unstructured/assessment.R>), [`FDR`](<Misspecified models/Unstructured/Unstructured/assessment_FDR.R>), [`RMSE`](<Misspecified models/Unstructured/Unstructured/rmse.R>) | [`simulation`](<Misspecified models/Unstructured/Directed/sim_gaussian_directed_DAGAR.R>), [`assessment`](<Misspecified models/Unstructured/Directed/assessment_directed.R>), [`RMSE`](<Misspecified models/Unstructured/Directed/rmse_directed.R>) | [`simulation`](<Misspecified models/Unstructured/Undirected/sim_gaussian_undirected_DAGAR.R>), [`assessment`](<Misspecified models/Unstructured/Undirected/assessment_undirected.R>), [`RMSE`](<Misspecified models/Unstructured/Undirected/rmse_undirected.R>) |
| Directed | [`simulation`](<Misspecified models/Directed/Unstructured/sim_gaussian_DAGAR.R>), [`assessment`](<Misspecified models/Directed/Unstructured/assessment.R>), [`FDR`](<Misspecified models/Directed/Unstructured/assessment_FDR.R>), [`RMSE`](<Misspecified models/Directed/Unstructured/rmse.R>) | [`simulation`](<Misspecified models/Directed/Directed/sim_gaussian_directed_DAGAR.R>), [`assessment`](<Misspecified models/Directed/Directed/assessment_directed.R>), [`FDR`](<Misspecified models/Directed/Directed/assessment_directed_FDR.R>), [`RMSE`](<Misspecified models/Directed/Directed/rmse_directed.R>) | [`simulation`](<Misspecified models/Directed/Undirected/sim_gaussian_undirected_DAGAR.R>), [`assessment`](<Misspecified models/Directed/Undirected/assessment_undirected.R>), [`RMSE`](<Misspecified models/Directed/Undirected/rmse_undirected.R>) |
| Undirected | [`simulation`](<Misspecified models/Undirected/Unstructured/sim_gaussian_DAGAR.R>), [`assessment`](<Misspecified models/Undirected/Unstructured/assessment.R>), [`FDR`](<Misspecified models/Undirected/Unstructured/assessment_FDR.R>), [`RMSE`](<Misspecified models/Undirected/Unstructured/rmse.R>) | [`simulation`](<Misspecified models/Undirected/Directed/sim_gaussian_directed_DAGAR.R>), [`assessment`](<Misspecified models/Undirected/Directed/assessment_directed.R>), [`RMSE`](<Misspecified models/Undirected/Directed/rmse_directed.R>) | [`simulation`](<Misspecified models/Undirected/Undirected/sim_gaussian_undirected_DAGAR.R>), [`assessment`](<Misspecified models/Undirected/Undirected/assessment_undirected.R>), [`FDR`](<Misspecified models/Undirected/Undirected/assessment_undirected_FDR.R>), [`RMSE`](<Misspecified models/Undirected/Undirected/rmse_undirected.R>) |

The grid's simulation scripts compile one of three shared samplers:
[`sampler_sim_gaussian_DAGAR_fastest.cpp`](<Misspecified models/sampler_sim_gaussian_DAGAR_fastest.cpp>),
[`sampler_sim_gaussian_directed_DAGAR_fastest.cpp`](<Misspecified models/sampler_sim_gaussian_directed_DAGAR_fastest.cpp>),
or [`sampler_sim_gaussian_undirected_DAGAR_fastest.cpp`](<Misspecified models/sampler_sim_gaussian_undirected_DAGAR_fastest.cpp>).

The CAR results similarly require all three generators
([`unstructured`](<CAR/sim_gaussian_CAR.R>),
[`directed`](<CAR/sim_gaussian_directed_CAR.R>), and
[`undirected`](<CAR/sim_gaussian_undirected_CAR.R>)) and their corresponding
C++ samplers
([`unstructured`](<CAR/sampler_sim_gaussian_CAR_fastest.cpp>),
[`directed`](<CAR/sampler_sim_gaussian_directed_CAR_fastest.cpp>), and
[`undirected`](<CAR/sampler_sim_gaussian_undirected_CAR_fastest.cpp>)).

### Supplementary figures

| Result | Repository files |
|---|---|
| Figure S1 | The DAGAR panel is generated by [`WAIC_comparison.R`](<Misspecified models/WAIC/WAIC_comparison.R>); the CAR panel is generated by [`WAIC_comparison_CAR.R`](<CAR/WAIC/WAIC_comparison_CAR.R>) |
| Figure S2 | Adjacency pooled FDR curves from [`postprocess_multiple_chains_tempering.R`](<data analysis/postprocess_multiple_chains_tempering.R>) |
| Figures S3--S6 | Adjacency pooled posterior interval figures from the same post-processing script |
| Figure S7 | California county-name map from [`exploratory_figures_seer.R`](<data analysis/exploratory_figures_seer.R>) |
| Figures S8--S16 | Mean--Adjacency posterior intervals, FDR curves, boundary maps, and non-adjacency map from [`postprocess_multiple_chains_tempering.R`](<data analysis/postprocess_multiple_chains_tempering.R>) |
| Figures S17--S24 | Mean posterior intervals, FDR curves, and boundary maps from the same post-processing script |
| Model-sensitivity figures | Edge-level boundary-probability agreement and Mean--Adjacency versus Mean boundary maps from [`sensitivity_analysis_multichain.R`](<data analysis/sensitivity_analysis_multichain.R>) |
| Sampler-sensitivity figures | Boundary-probability and fitted-log-risk agreement from [`sampler_sensitivity_multichain.R`](<data analysis/sampler_sensitivity_multichain.R>) |
| Posterior-predictive figures | Global checks and observed versus fitted counts from [`posterior_predictive_checks_multichain.R`](<data analysis/posterior_predictive_checks_multichain.R>) |

### Known reproducibility limitation outside the data-analysis workflow

- Figure S1 is assembled from two separate simulation plotting scripts rather
  than produced by a single figure script. This will be addressed when the
  simulation experiments are reorganized.
