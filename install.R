cran_packages <- c(
  "AICcmodavg",
  "blockmatrix",
  "boot",
  "CARBayes",
  "caret",
  "classInt",
  "coda",
  "cowplot",
  "dplyr",
  "fields",
  "generics",
  "ggmcmc",
  "ggplot2",
  "ggpubr",
  "gtools",
  "hesim",
  "Hmisc",
  "invgamma",
  "LaplacesDemon",
  "lattice",
  "loo",
  "magic",
  "mapproj",
  "maps",
  "MASS",
  "Matrix",
  "matrixStats",
  "mcmc",
  "mcmcse",
  "monomvn",
  "msos",
  "mvtnorm",
  "plyr",
  "RColorBrewer",
  "Rcpp",
  "RcppArmadillo",
  "readr",
  "rmapshaper",
  "sf",
  "spdep",
  "stringr",
  "tictoc",
  "tigris",
  "tidyr",
  "tidyverse"
)

installed <- rownames(installed.packages())
missing <- setdiff(cran_packages, installed)

if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

remaining <- setdiff(cran_packages, rownames(installed.packages()))

if (length(remaining) > 0) {
  warning(
    "Some packages were not installed: ",
    paste(remaining, collapse = ", "),
    call. = FALSE
  )
}

message(
  "Note: the archived maptools, rgeos, and gpclib packages have been replaced ",
  "with sf-based polygon construction, centroid extraction, and boundary plotting."
)

message(
  "Note: RcppGSL and RcppEigen are not installed by default because the ",
  "current sampler code does not use GSL or Eigen APIs."
)
