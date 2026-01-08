# graphical-models-for-boundary-detection
This is a repository containing all the codes used in "Detecting Spatial Health Disparities Using Disease Maps"

The data analysis folder contains the data along with the script through which obtaining the results regarding the three settings under which the analysis was conducted. The data folder contains the following files: 
- SIR_adjusted.csv: observed and expected counts of the 4 cancers along with the standardized ratio for each California county;
- covariates.csv: for each county in the US there are reported the percentages of people under-18, over-65, with high school education, the percentage of families below the poverty threshold and the percentage of unemployment;
- insurance.csv: for each county in the US there are reported the percentages of different insurance status (e.g. uninsured, medicalaid, insured);
- race.csv: for each county in the US there are reported race percentages (white, black, other);
- sex.csv: for each county in the US there are reported sex percentages (male, female);
- smoking: smoking rates for each California county.

For all of the three cases it was used the unstructured disease graph:
- main.R: data preparation, posterior analysis
- sampler.rcpp: metropolis within gibbs algorithm written with RcppArmadillo

Note that into main.R and sampler.rcpp there is a parameter cvrts that can take 3 values: "adj" for covariates only into the adjacency model, "mean" only in the mean structure and "meanadj" both in the mean and in the adjacency.
 
The simulation codes are contained in the "Misspecified models" and "CAR" folders. The Misspecified models folder contains the codes used to produce the results shown in the manuscript. To help the reader navigate the folder the folder are nested in the following way: 

1. once you access "Misspecified models" you find three folders called "Directed", "Undirected" and "Unstrucured". They refer to the data generating process.
2. each of these three folder contains other three folders named "Directed", "Undirected" and "Unstrucured", referring to the two misspecified models and the one fitted under the true one. They contain codes for generating results for all the scenarios.

Note that: in each subfolder the one called "RE_generation_X_DAGAR" contains the generated data; in each file called assessment_X.R there are some files to be called referring to the folder "runs_X", they are the results that are obtained by running the codes in "sim_gaussian_X_DAGAR.R". 

The WAIC and rmse folder contain codes to generate tables and figures regarding WAIC and rmse. 

The CAR folder contains codes to reproduce the results concerning the comparison in the supplementary materials with the MCAR specification. 


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







