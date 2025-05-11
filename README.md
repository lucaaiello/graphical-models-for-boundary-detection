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
  
The simulation folder contains the codes used to produce the results shown in the manuscript generating data through the unstructured, directed and undirected disease graph fitted under the unstructured model. In particular there are the data simulation, the simulation, the assessment for each disease graphical model.


