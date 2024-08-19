# graphical-models-for-boundary-detection
This is a repository containing all the codes used in "Detecting Spatial Health Disparities Using Disease Maps - Aiello &amp; Banerjee"

The data analysis folder contains the data along with the script through which obtaining the results regarding the three settings under which the analysis was conducted. For all of the three cases it was used the unstructured disease graph:
- main.R: data preparation, posterior analysis
- sampler.rcpp: metropolis within gibbs algorithm written with RcppArmadillo
  
The simulation folder contains the codes used to produce the results shown in the manuscript for the unstructured, directed and undirected disease graph under their true model. In particular there are the data simulation, the simulation, the assessment for each disease graphical model.


