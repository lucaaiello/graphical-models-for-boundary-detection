# Graphical Models for Boundary Detection

This repository contains all code used in the paper  
**“Detecting Spatial Health Disparities Using Disease Maps.”**

The repository is organized into three main components: data analysis, simulation studies, and model comparison.

---

## Data Analysis

The `data_analysis/` folder contains the data and scripts used to produce the results for the three empirical analysis settings considered in the paper.

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

---

### Model Fitting

For all three empirical settings, an **unstructured disease graph** is used.

- **`main.R`**  
  Data preprocessing, model fitting, and posterior analysis.

- **`sampler.rcpp`**  
  Metropolis-within-Gibbs MCMC sampler implemented in **RcppArmadillo**.

Both `main.R` and `sampler.rcpp` include a parameter `cvrts`, which controls the role of covariates in the model:
- `"adj"`: covariates enter only the adjacency (boundary) model;
- `"mean"`: covariates enter only the mean structure;
- `"meanadj"`: covariates enter both the mean and adjacency components.

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

---

## Model Comparison: WAIC and RMSE

The `WAIC and rmse/` folder contains scripts to compute and visualize:
- Widely Applicable Information Criterion (WAIC),
- Root Mean Squared Error (RMSE).

These scripts generate the tables and figures reported in the manuscript.

---

## CAR Model Comparison

The `CAR/` folder contains code to reproduce the results in the **Supplementary Materials**, where the proposed approach is compared with the **MCAR (Multivariate Conditional Autoregressive)** specification.

Here a brief guideline on where to find code to generate the plots and tables regarding the results in the manuscript:

Figure 1, 2, 3, 5, 6, 7, 8, 9: data analysis/main.R

Table 1: Misspecified models/WAIC_assessment.R

Table 2, 3: 
- Misspecified models/Unstructured/Unstructured/assessment_FDR.R (top row),
- Misspecified models/Directed/Unstructured/assessment_FDR.R (middle row),
- Misspecified models/Undirected/Unstructured/assessment_FDR.R (bottom row)

Table S1: 
- Misspecified models/Unstructured/Unstructured/assessment.R (top rows),
- Misspecified models/Directed/Unstructured/assessment.R (middle rows),
- Misspecified models/Undirected/Unstructured/assessment.R (bottom rows)

Table S2: Misspecified models/rmse_assessment.R

Table S3:
- Misspecified models/Unstructured/Unstructured/assessment_FDR.R (top row),
- Misspecified models/Directed/Directed/assessment_directed_FDR.R (middle row),
- Misspecified models/Undirected/Undirected/assessment_undirected_FDR.R (bottom row)

Table S4, S7:
- Misspecified models/Unstructured/Unstructured/assessment.R (top rows),
- Misspecified models/Directed/Directed/assessment_directed.R (middle rows),
- Misspecified models/Undirected/Undirected/assessment_undirected.R (bottom rows)

Table S5:
- CAR/assessment_FDR.R (top row),
- CAR/assessment_directed_FDR.R (middle row),
- CAR/assessment_undirected_FDR.R (bottom row)

Table S6, S8:
- CAR/assessment.R (top rows),
- CAR/assessment_directed.R (middle rows),
- CAR/assessment_undirected.R (bottom rows)
  
Figure S1: Misspecified models/WAIC_comparison.R
Figure S2: CAR/WAIC/WAIC_comparison_CAR.R
Table S6, S7, S8: data analysis/main.R
Figure S3-S6/S8-S24: data analysis/main.R







