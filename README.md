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
Rscript "data analysis/main.R"
```

The sampler source files are compiled with `Rcpp::sourceCpp()`, so a working C/C++
toolchain is required. On Windows this usually means Rtools; on macOS this usually
means Xcode command line tools; on Linux this usually means the system compiler
toolchain.

The C++ sampler files use `RcppArmadillo`.

California county polygons are constructed with `sf`, replacing the archived
`maptools` package. The geometries remain as `sf` objects for centroid extraction,
adjacency construction, and boundary plotting while preserving county ordering and
the original adjacency outputs. The archived `maptools`, `rgeos`, and `gpclib`
packages are no longer required.

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
| [`data analysis/main.R`](<data analysis/main.R>) | The six CSV files in `data analysis/data/`, [`sampler.cpp`](<data analysis/sampler.cpp>), and `results_adj.RData` in the current configuration | Creates the California exploratory maps, posterior summaries, diagnostics, FDR curves, and boundary maps in the active R graphics device; the fitted samples are held in `mcmc_samples` |
| `sim_*.R` simulation generators | The scenario settings embedded in the script, the corresponding C++ sampler, and an existing matching `runs_*` directory | One `mcmc_samples_<seed>.rds` file per replicate, plus the simulated truth and design objects saved in the matching `RE_generation_*` directory |
| `assessment*.R` | The matching `RE_generation_*` objects and `mcmc_samples_<seed>.rds` files | Prints sensitivity, specificity, and adjacency-recovery summaries; non-FDR assessment scripts also save one WAIC vector |
| `assessment*_FDR.R` | The matching generated truth and saved MCMC replicates | Prints FDR-based sensitivity and specificity summaries |
| `rmse*.R` | The matching generated truth and saved MCMC replicates | Saves one compact RMSE vector in `Misspecified models/rmse/` |
| [`WAIC_assessment.R`](<Misspecified models/WAIC/WAIC_assessment.R>) and [`rmse_assessment.R`](<Misspecified models/rmse/rmse_assessment.R>) | The compact WAIC or RMSE `.rds` files created by the preceding scripts | Prints the final table values |
| WAIC comparison scripts | The corresponding WAIC vectors | Draws the WAIC density plots used in Figure S1 |
| Comparison folders | Their `sim_*.R` or comparison generator followed by the matching `assessment.R` | Produces the Lee--Mitchell and Li et al. comparison summaries |

The empirical script displays figures and tables interactively; it does not
automatically write them to image or table files. They can be saved from the R
graphics device or by wrapping the relevant plotting calls in `pdf()`,
`png()`, or `ggsave()`.

### Low-level sampler interface

The C++ samplers are exposed to R by `Rcpp::sourceCpp()`. Direct use is possible,
but the calling user is responsible for constructing the same county ordering,
disease stacking, and matrix dimensions used by the supplied R scripts.

The main empirical function is:

```r
Rcpp::sourceCpp("data analysis/sampler.cpp")

fit <- MADAGAR(
  y = Y, X = X,
  Z1 = Z1, Z2 = Z2, Z3 = Z3,
  E = E, cvrts = cvrts,
  q = 4, Winc = Winc, Minc = Minc,
  alpha = 1, n_atoms = 15,
  runs = 10000, burn = 50000, thin = 5
)
```

Its arguments are:

| Argument | Meaning and required structure |
|---|---|
| `y` | Numeric outcome vector of length `n * q`, stacked by disease after applying the script's county ordering |
| `X` | Design matrix with `n * q` rows; the columns correspond to the regression coefficients being fitted |
| `Z1`, `Z2`, `Z3` | Symmetric `n x n` matrices of standardized covariate dissimilarities between geographic neighbors |
| `E` | Expected-count vector of length `n * q` for the empirical Poisson model |
| `cvrts` | One of `"adj"`, `"mean"`, or `"meanadj"`, controlling whether covariates enter the adjacency model, mean model, or both |
| `q` | Number of diseases; `q = 4` in the paper |
| `Winc` | Lower-triangular geographic adjacency matrix supplied to the empirical sampler |
| `Minc` | Symmetric binary `n x n` geographic adjacency matrix in the reordered county indexing |
| `alpha` | Dirichlet-process concentration parameter |
| `n_atoms` | Truncation level for the discrete mixture |
| `runs` | Number of posterior draws retained |
| `burn` | Number of initial MCMC iterations discarded |
| `thin` | Number of MCMC iterations between retained draws |

The returned list contains posterior samples for `beta`, `phi`, `theta`, `u`,
`rho`, `V`, `r`, `F_r`, `eta`, `tau`, `W1`--`W4`, and `A`. Matrix-valued
parameters have one row per retained draw. Each `Wd` component is a list of
`runs` sampled `n x n` adjacency matrices for disease `d`.

| Output component | Contents |
|---|---|
| `beta` | Regression coefficients |
| `taud` | Gaussian observation precisions; present in the simulation samplers but not the empirical Poisson sampler |
| `phi` | Discrete spatial random effects |
| `theta` | Values of the mixture atoms |
| `u` | Atom-allocation indices for the region-disease effects |
| `rho` | Disease-specific spatial autocorrelation parameters |
| `V` or `v` | Stick-breaking variables |
| `r` | Latent Gaussian spatial variables |
| `F_r` | Probability-integral transforms of the latent variables |
| `eta` | Coefficients controlling covariate-informed geographic adjacency |
| `tau` | Precision parameter for the mixture atoms |
| `W1`--`W4` | Posterior samples of the disease-specific geographic adjacency matrices |
| `A`, directed coefficients, or undirected correlation | Parameters governing dependence between diseases, according to the selected disease graph |

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
| Produce exploratory maps and summaries from existing empirical results | Seconds to a few minutes | Does not include fitting the MCMC model |
| Fit one empirical SEER model in `data analysis/main.R` | Approximately 10 minutes of sampler time | The script requests 100,000 C++ iterations: `burn = 50000`, `runs = 10000`, `thin = 5`; allow additional time for compilation and post-processing |
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

- **`main.R`**  
  Data preprocessing, model fitting, and posterior analysis.

- **`sampler.cpp`**  
  Metropolis-within-Gibbs MCMC sampler implemented in **RcppArmadillo**.

Both `main.R` and `sampler.cpp` include a parameter `cvrts`, which controls the role of covariates in the model:
- `"adj"`: covariates enter only the adjacency (boundary) model;
- `"mean"`: covariates enter only the mean structure;
- `"meanadj"`: covariates enter both the mean and adjacency components.

The `"mean"` setting keeps the geographic adjacency matrices fixed throughout
the MCMC run. It is therefore the fixed-adjacency MARDP/MDAGAR special case most
closely corresponding to Gao et al. (2023), whereas `"adj"` and `"meanadj"`
activate the adjacency-modeling extension developed in this paper. The saved
`results_mean.RData` fit is used for Table S18 and Figures S18--S25.

### Result index for `data analysis/main.R`

The following index identifies the current line ranges that produce each
empirical table and figure. The file also contains searchable comments beginning
with `# RESULT:` or `# RESULTS:` so that the relevant blocks remain easy to find
if later edits shift the line numbers.

Shared data preparation, county ordering, model-matrix construction, and model
fitting occur in lines 1--349. The posterior result file selected around lines
332--336 determines which supplementary analysis is produced:

- `results_adj.RData`: main Figures 5--9, Table S16, and Figures S2--S6;
- `results_meanadj.RData`: Table S17 and Figures S8--S17;
- `results_mean.RData`: Table S18 and Figures S18--S25.

#### Main manuscript

| Result | Current lines in `data analysis/main.R` |
|---|---:|
| Figure 1: cancer SIR maps | 499--568 |
| Figure 2: Moran's I | 137--183 |
| Figure 3: smoking, elderly-population, and poverty maps | 570--614 |
| Figure 4: disease-graph schematic | Not generated by `main.R` |
| Figure 5: estimated FDR curves | 624--710 |
| Figure 6: disease-specific difference-boundary maps | 624--876 |
| Figure 7: shared difference-boundary maps | 951--1213 |
| Figure 8: mutual cross-difference-boundary maps | 1215--1446 |
| Figure 9: non-adjacency maps | 1447--1575 |

#### Supplementary material

| Result | Model result file | Current lines in `data analysis/main.R` |
|---|---|---:|
| Figure S2 | `results_adj.RData` | 1577--1623 |
| Table S16 | `results_adj.RData` | 350--369 |
| Figure S3 | `results_adj.RData` | 383--407 |
| Figure S4 | `results_adj.RData` | 408--431 |
| Figure S5 | `results_adj.RData` | 433--441 |
| Figure S6 | `results_adj.RData` | 443--453 |
| Figure S7 | -- | Not generated by `main.R` |
| Figure S8 | `results_meanadj.RData` | 1577--1623 |
| Table S17 | `results_meanadj.RData` | 350--369 |
| Figure S9 | `results_meanadj.RData` | 383--407 |
| Figure S10 | `results_meanadj.RData` | 408--431 |
| Figure S11 | `results_meanadj.RData` | 433--441 |
| Figure S12 | `results_meanadj.RData` | 443--453 |
| Figure S13 | `results_meanadj.RData` | 624--710 |
| Figure S14 | `results_meanadj.RData` | 624--876 |
| Figure S15 | `results_meanadj.RData` | 951--1213 |
| Figure S16 | `results_meanadj.RData` | 1215--1446 |
| Figure S17 | `results_meanadj.RData` | 1447--1575 |
| Table S18 | `results_mean.RData` | 350--369 |
| Figure S18 | `results_mean.RData` | 383--407 |
| Figure S19 | `results_mean.RData` | 408--431 |
| Figure S20 | `results_mean.RData` | 433--441 |
| Figure S21 | `results_mean.RData` | 455--488 |
| Figure S22 | `results_mean.RData` | 624--710 |
| Figure S23 | `results_mean.RData` | 624--876 |
| Figure S24 | `results_mean.RData` | 951--1213 |
| Figure S25 | `results_mean.RData` | 1215--1446 |

Figures S2 and S8 are adjacency-cardinality traceplots and therefore apply only
to the two models that estimate adjacency. Similarly, the non-adjacency block is
not used by the fixed-adjacency `results_mean.RData` analysis.

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

This mapping corresponds to the current main manuscript (Figures 1--9 and no
numbered tables) and supplementary material (Tables S1--S18 and Figures
S1--S25).

### Main manuscript

| Result | Description | Repository files |
|---|---|---|
| Figure 1 | Cancer SIR maps | [`data analysis/main.R`](<data analysis/main.R>), using [`SIR_adjusted.csv`](<data analysis/data/SIR_adjusted.csv>) |
| Figure 2 | Moran's I by neighbor order | [`data analysis/main.R`](<data analysis/main.R>), using [`SIR_adjusted.csv`](<data analysis/data/SIR_adjusted.csv>) |
| Figure 3 | Smoking, elderly-population, and poverty maps | [`data analysis/main.R`](<data analysis/main.R>), using [`smoking.csv`](<data analysis/data/smoking.csv>) and [`covariates.csv`](<data analysis/data/covariates.csv>) |
| Figure 4 | Three disease-graph schematics | No generating script or source graphic is currently present in the repository |
| Figure 5 | Estimated FDR curves | [`data analysis/main.R`](<data analysis/main.R>) with [`results_adj.RData`](<data analysis/results_adj.RData>) |
| Figure 6 | Disease-specific difference boundaries | [`data analysis/main.R`](<data analysis/main.R>) with [`results_adj.RData`](<data analysis/results_adj.RData>) |
| Figure 7 | Shared difference boundaries | [`data analysis/main.R`](<data analysis/main.R>) with [`results_adj.RData`](<data analysis/results_adj.RData>) |
| Figure 8 | Mutual cross-difference boundaries | [`data analysis/main.R`](<data analysis/main.R>) with [`results_adj.RData`](<data analysis/results_adj.RData>) |
| Figure 9 | Estimated non-adjacencies | [`data analysis/main.R`](<data analysis/main.R>) with [`results_adj.RData`](<data analysis/results_adj.RData>) |

The empirical result files are generated by the model in
[`data analysis/main.R`](<data analysis/main.R>), whose MCMC implementation is
[`data analysis/sampler.cpp`](<data analysis/sampler.cpp>).

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
| Table S16 | [`data analysis/main.R`](<data analysis/main.R>) with [`results_adj.RData`](<data analysis/results_adj.RData>) |
| Table S17 | [`data analysis/main.R`](<data analysis/main.R>) with [`results_meanadj.RData`](<data analysis/results_meanadj.RData>) |
| Table S18 | [`data analysis/main.R`](<data analysis/main.R>) with [`results_mean.RData`](<data analysis/results_mean.RData>) |

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
| Figure S2 | [`data analysis/main.R`](<data analysis/main.R>) with [`results_adj.RData`](<data analysis/results_adj.RData>) |
| Figures S3--S6 | HPD plots from [`data analysis/main.R`](<data analysis/main.R>) with [`results_adj.RData`](<data analysis/results_adj.RData>) |
| Figure S7 | California county-name map; no code that adds the county labels is currently present in the repository |
| Figure S8 | [`data analysis/main.R`](<data analysis/main.R>) with [`results_meanadj.RData`](<data analysis/results_meanadj.RData>) |
| Figures S9--S12 | HPD plots from [`data analysis/main.R`](<data analysis/main.R>) with [`results_meanadj.RData`](<data analysis/results_meanadj.RData>) |
| Figure S13 | Estimated FDR curves from [`data analysis/main.R`](<data analysis/main.R>) with [`results_meanadj.RData`](<data analysis/results_meanadj.RData>) |
| Figures S14--S17 | Boundary and non-adjacency maps from [`data analysis/main.R`](<data analysis/main.R>) with [`results_meanadj.RData`](<data analysis/results_meanadj.RData>) |
| Figures S18--S21 | HPD plots from [`data analysis/main.R`](<data analysis/main.R>) with [`results_mean.RData`](<data analysis/results_mean.RData>) |
| Figure S22 | Estimated FDR curves from [`data analysis/main.R`](<data analysis/main.R>) with [`results_mean.RData`](<data analysis/results_mean.RData>) |
| Figures S23--S25 | Boundary maps from [`data analysis/main.R`](<data analysis/main.R>) with [`results_mean.RData`](<data analysis/results_mean.RData>) |

### Known reproducibility gaps

- `data analysis/main.R` currently hard-codes `cvrts = "adj"` in the sampler
  call and explicitly loads only `results_adj.RData`. Reproducing Table S17,
  Table S18, and Figures S8--S25 therefore currently requires manually changing
  the model variant and loaded result file.
- Figure 4 and Figure S7 do not have complete generating code or source graphics
  in the repository.
- Figure S1 is assembled from two separate plotting scripts rather than produced
  by a single figure script.








