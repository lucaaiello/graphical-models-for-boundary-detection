rm(list = ls())

# Libraries --------------------------------------------------------------------

library(maps)
#Import California map
ca.county <- maps::map("county", "california", fill = TRUE, plot = FALSE)
library(readr)
library(dplyr)
library(spdep)
library(sf)
library(classInt)
library(RColorBrewer)
library(tidyr)
library(MASS)
library(Matrix)
library(Hmisc)
library(mapproj)
library(lattice)
library(mvtnorm)
library(matrixStats)
library(fields)
library(boot)
library(blockmatrix)
library(ggmcmc)
library(mcmc)
library(magic)
library(msos)
library(AICcmodavg)
library(coda)
library(invgamma)
library(mcmcse)
library(LaplacesDemon)
library(gtools)
library(ggpubr)

# Import covariates ------------------------------------------------------------

required_input_files <- file.path(
  "data analysis",
  "data",
  c(
    "SIR_adjusted.csv",
    "covariates.csv",
    "race.csv",
    "sex.csv",
    "insurance.csv",
    "smoking.csv"
  )
)
missing_input_files <- required_input_files[!file.exists(required_input_files)]
if (length(missing_input_files) > 0L) {
  stop(
    "Missing required input data file(s): ",
    paste(missing_input_files, collapse = ", ")
  )
}

rate_5y<- read_csv("data analysis/data/SIR_adjusted.csv")

covariates <- read_csv("data analysis/data/covariates.csv")
race <- read_csv("data analysis/data/race.csv")
sex <- read_csv("data analysis/data/sex.csv")
insurance <- read_csv("data analysis/data/insurance.csv")
smoking <- read.csv("data analysis/data/smoking.csv")
smoking$smoking <- as.numeric(substr(smoking$Cigarette.Smoking.Rate, 1,4))

# Import incidence data for 4 cancers in California ----------------------------

rate_lung <- rate_5y %>% filter(Site.code == "Lung and Bronchus")
rate_esophagus <- rate_5y %>% filter(Site.code == "Esophagus")
rate_larynx <- rate_5y %>% filter(Site.code == "Larynx")
rate_colrect <- rate_5y %>% filter(Site.code == "Colon and Rectum")

ca.poly <- st_as_sf(ca.county)
st_crs(ca.poly) <- NA
county.ID <- sub("^.*,", "", ca.poly$ID)
ca.poly$county <- county.ID
ca.poly$rate_lung <- rate_lung$standard_ratio
ca.poly$rate_esophagus <- rate_esophagus$standard_ratio
ca.poly$rate_larynx <- rate_larynx$standard_ratio
ca.poly$rate_colrect <- rate_colrect$standard_ratio
ca.poly$smoking <- smoking$smoking

ca.coords <- st_coordinates(st_centroid(st_geometry(ca.poly)))

# Exploratory analysis for Pearson correlation ---------------------------------

cancer_ratio <- cbind(rate_lung$standard_ratio, rate_esophagus$standard_ratio,
                      rate_larynx$standard_ratio, rate_colrect$standard_ratio)
colnames(cancer_ratio) <- c("Lung and Bronchus", "Esophageal",
                            "Larynx", "Colon and Rectum")
cor_matrix <- rcorr(cancer_ratio)
cor_matrix$r
cor_matrix$P

# Data construction ------------------------------------------------------------

county_attribute <- covariates[substr(covariates$State_county,1,2) == "CA",]
county_attribute$state <- extract_numeric(county_attribute$State_county)

county_attribute1 <- data.frame(county_attribute[order(county_attribute$state),])
county_attribute1$V_Persons_age_18_ACS_2012_2016 <- as.numeric(county_attribute1$V_Persons_age_18_ACS_2012_2016)/100
county_attribute1$V_Persons_age_65_ACS_2012_2016 <- as.numeric(county_attribute1$V_Persons_age_65_ACS_2012_2016)/100
county_attribute1$VHighschooleducationACS2012201 <- as.numeric(county_attribute1$VHighschooleducationACS2012201)/100
county_attribute1$VFamiliesbelowpovertyACS201220 <- as.numeric(county_attribute1$VFamiliesbelowpovertyACS201220)/100
county_attribute1$V_Unemployed_ACS_2012_2016 <- as.numeric(county_attribute1$V_Unemployed_ACS_2012_2016)/100

race1 <- race[substr(race$State_county,1,2) == "CA"&race$Race_recode_White_Black_Other=="Black",]
sex1 <- sex[substr(sex$State_county,1,2) == "CA"&sex$Sex=="Male",]
insurance1 <- insurance[substr(insurance$State_county,1,2) == "CA"&insurance$Insurance_Recode_2007=="Uninsured",]

rate_lung1 <- cbind(rate_lung, smoking$smoking, county_attribute1[,2:6],
                    race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
rate_lung1 <- rate_lung1[,-1]
colnames(rate_lung1) <- c("county", "O_count", "E_count", "standard_ratio",
                          "site", "smoking", "young","old", "highschool",
                          "poverty", "unemployed", "black", "male",
                          "uninsured")

rate_esophagus1 <- cbind(rate_esophagus, smoking$smoking, county_attribute1[,2:6],
                         race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
rate_esophagus1 <- rate_esophagus1[,-1]
colnames(rate_esophagus1) <- c("county", "O_count", "E_count", "standard_ratio",
                               "site", "smoking", "young","old", "highschool",
                               "poverty", "unemployed", "black", "male",
                               "uninsured")

rate_larynx1 <- cbind(rate_larynx, smoking$smoking, county_attribute1[,2:6],
                      race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
rate_larynx1 <- rate_larynx1[,-1]
colnames(rate_larynx1) <- c("county", "O_count", "E_count", "standard_ratio",
                            "site", "smoking", "young","old", "highschool",
                            "poverty", "unemployed", "black", "male",
                            "uninsured")

rate_colrect1 <- cbind(rate_colrect, smoking$smoking, county_attribute1[,2:6],
                       race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
rate_colrect1 <- rate_colrect1[,-1]
colnames(rate_colrect1) <- c("county", "O_count", "E_count", "standard_ratio",
                             "site", "smoking", "young","old", "highschool",
                             "poverty", "unemployed", "black", "male",
                             "uninsured")

# Adjacency matrix and neighbor info -------------------------------------------

ca.neighbors <- poly2nb(ca.poly, row.names = county.ID)
n <- length(ca.neighbors)

Adj <- sapply(ca.neighbors, function(x,n) { v <- rep(0,n); v[x] <- 1; v}, n)
colnames(Adj) <- county.ID

ca.coord <- st_coordinates(st_centroid(st_geometry(ca.poly)))
ca.latrange <- round(quantile(ca.coord[,2],c(0.25,0.75)))
ca.albersproj <- mapproject(ca.coord[,1],ca.coord[,2],
                            projection = "albers",param=ca.latrange)

# Data preparation: Moran's I by neighbor order -------------------------------
# This calculation is retained for the SEER setup. The reported Main Figure 2
# is generated and saved by exploratory_figures_seer.R.

projmat <- cbind(ca.albersproj$x,ca.albersproj$y)
dmat <- as.matrix(dist(projmat))

moranI <- function(y, A){
  n <- length(y)
  nom_sum <- 0
  den_sum <- 0
  for(i in 1:n){
    den_sum <- den_sum + (y[i]-mean(y))^2
    for(j in 1:n){
      nom_sum <- nom_sum + A[i,j]*(y[i]-mean(y))*(y[j]-mean(y))
    }
  }
  return(n*nom_sum/sum(A)/den_sum)
}

lung_moran <- c()
esophagus_moran <- c()
larynx_moran <- c()
colrect_moran <- c()
for(lag in 1:11){
  A_1 <- as.matrix((dmat <= 0.01*lag & dmat > 0.01*(lag-1))*1)
  diag(A_1) <- 0
  lung_moran[lag] <- as.numeric(moranI(rate_lung$standard_ratio, A_1))
  esophagus_moran[lag] <- as.numeric(moranI(rate_esophagus$standard_ratio, A_1))
  larynx_moran[lag] <- as.numeric(moranI(rate_larynx$standard_ratio, A_1))
  colrect_moran[lag] <- as.numeric(moranI(rate_colrect$standard_ratio, A_1))
}

moran_value <- c(lung_moran, esophagus_moran, larynx_moran, colrect_moran)
cancer <- c(rep("Lung", 11), rep("Esophageal", 11), rep("Larynx", 11), rep("Colorectum", 11))

df <- data.frame(cancer)
df$moran_value <- moran_value
df$r <- rep(1:11, 4)
df$cancer <- factor(df$cancer, levels = c("Lung", "Esophageal", "Larynx", "Colorectum"))

if (FALSE) {
  ggplot(df, aes(r, moran_value)) + geom_point(size = 3) +
    ylab("Moran's I") +
    facet_wrap(~cancer, nrow = 1, ncol = 4) +
    theme_bw() +
    scale_x_continuous(breaks = 1:11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 15),
      strip.text = element_text(size = 15)
    )
}

# Neighboring counties information

neighborvec0 <- NULL
neighbor_list0 <- NULL
neighbor_name <- NULL
for(i in 1:(n-1)){
  for(j in (i+1):n){
    if(Adj[i,j] == 1){
      neighborvec0 <- c(neighborvec0, paste(i, ",", j, sep=""))
      neighbor_list0 <- rbind(neighbor_list0, c(i,j))
      neighbor_name <- c(neighbor_name, paste(colnames(Adj)[i], ", ",
                                              colnames(Adj)[j], sep=""))
    }
  }
}

# Specify the order and reorder the map with new neighbor info -----------------

perm <- order(ca.albersproj$x-ca.albersproj$y)
Adj_new <- Adj[perm,perm]

n <- nrow(Adj_new)
ni <- rowSums(Adj_new)
maxn <- max(ni)
neimat <- matrix(0,n,maxn)
neighbors <- lapply(1:n,function(x) which(Adj_new[x,]==1))
#N(i): 2:n
dneighbors <- sapply(2:n,function(i) intersect(neighbors[[i]],1:(i-1)))
#n<i: 2:n
dni <- sapply(dneighbors,length)
original_perm <- 1:58
index2 <- c(1,which(dni==0)+1)

final_perm <- c(original_perm[perm][index2],
                original_perm[perm][-index2])
final_perm[order(final_perm)]

Minc <- Adj[final_perm,final_perm]

Winc <- Minc
Winc[upper.tri(Winc)] <- 0

n <- nrow(Minc)
ni <- rowSums(Minc)
maxn <- max(ni)
neimat <- matrix(0,n,maxn)
neighbors <- lapply(1:n,function(x) which(Minc[x,]==1))
#N(i): 2:n
dneighbors <- sapply(2:n,function(i) intersect(neighbors[[i]],1:(i-1)))
#n<i: 2:n

dni <- sapply(dneighbors,length)
nmax <- max(dni)
cni <- cumsum(dni)
dneimat <- sapply(dneighbors, function(nei,nmax,n) c(nei,rep(n+1,nmax+1-length(nei))),nmax,n)
udnei <- unlist(dneighbors)

ni_wo <- sapply(neighbors,length)
cni_wo <- cumsum(ni_wo)
udnei_wo <- unlist(neighbors)
cn <- c(0, cni)
ns <- dni

region <- seq(1:n)
cn <- c(0, cni)
ns <- dni
index <- list()
for(i in 1:(n-2)){
  index[[i]] <- region[-(udnei[(cn[i+1] + 1):(cn[i+1] + ns[i+1])])]
}
index1 <- unlist(index)

# Function to calculate percentage differences based on min-max range
sd_diff_mat <- function(vector,Minc) {
  n <- length(vector)

  matrix_sd_differences <- matrix(NA, nrow = n, ncol = n)

  indices <- which(Minc == 1, arr.ind = TRUE)

  for (i in 1:nrow(indices)) {
    row_idx <- indices[i, "row"]
    col_idx <- indices[i, "col"]

    abs_difference <- abs(vector[row_idx] - vector[col_idx])

    matrix_sd_differences[row_idx, col_idx] <- abs_difference

  }

  return(matrix_sd_differences/sd(matrix_sd_differences[!is.na(matrix_sd_differences)]))

}

# Data in model ----------------------------------------------------------------

Y1 <- rate_lung1$O_count[final_perm]
Y2 <- rate_esophagus1$O_count[final_perm]
Y3 <- rate_larynx1$O_count[final_perm]
Y4 <- rate_colrect1$O_count[final_perm]

E1 <- rate_lung1$E_count[final_perm]
E2 <- rate_esophagus1$E_count[final_perm]
E3 <- rate_larynx1$E_count[final_perm]
E4 <- rate_colrect1$E_count[final_perm]

X1 <- as.matrix(cbind(1,rate_lung1[,6:14]))[final_perm,]
X2 <- as.matrix(cbind(1,rate_esophagus1[,6:14]))[final_perm,]
X3 <- as.matrix(cbind(1,rate_larynx1[,6:14]))[final_perm,]
X4 <- as.matrix(cbind(1,rate_colrect1[,6:14]))[final_perm,]

Z1 <- sd_diff_mat(X1[,2],Minc)
Z2 <- sd_diff_mat(X2[,4],Minc)
Z3 <- sd_diff_mat(X3[,6],Minc)

Y <- c(Y1,Y2,Y3,Y4)
E <- c(E1, E2, E3, E4)

cvrts <- tolower(Sys.getenv("MCMC_MULTICHAIN_SPECIFICATION", unset = "adj"))
# "adj" for covariates only in the adjacency model;
# "mean" for covariates only in the mean structure, reproducing the
# Gao et al. (2023)-style specification; and
# "meanadj" for covariates in both the mean and adjacency models.

stopifnot(cvrts %in% c("adj", "mean", "meanadj"))

if (cvrts == "mean" | cvrts == "meanadj"){
  X <- as.matrix(bdiag(bdiag(X1[,c(1,2,4,6)], X2[,c(1,2,4,6)]),
                       bdiag(X3[,c(1,2,4,6)], X4[,c(1,2,4,6)])))
} else if (cvrts == "adj") {
  X <- as.matrix(bdiag(bdiag(X1[,1], X2[,1]),
                       bdiag(X3[,1], X4[,1])))
}

library(Rcpp)
library(RcppArmadillo)
library(parallel)

library(tictoc)

# Baseline multiple-chain MCMC workflow ------------------------------
# Six independent cold posterior chains are run with the original transition
# structure: no tempering, targeted cancer-specific moves, latent-r block moves,
# or covariance transport. The correction is limited to posterior consistency:
# whenever rho, eta, or A changes the marginal Gaussian-copula transform, F_r
# and u are updated in the same MH proposal and the induced likelihood change
# is included. Proposal adaptation is restricted to warm-up.

# Start with "quick" to verify compilation and data flow. Use "pilot" before
# committing to "production". The "smoke" mode exists only for a minimal
# engineering test. The environment variable permits non-interactive selection,
# for example MCMC_MULTICHAIN_MODE=pilot.
run_mode <- Sys.getenv("MCMC_MULTICHAIN_MODE", unset = "production")
run_sampling <- TRUE
# Reuse a completed chain only when its seed, MCMC settings, initialization, and
# C++ source checksum all match the current configuration.
reuse_completed_chains <- TRUE
# The map is always based on pooled posterior probabilities from all chains.
write_pooled_boundary_map <- TRUE

n_chains <- 6L
chain_seeds <- c(
  12345L, 23456L, 34567L,
  45678L, 56789L, 67890L
)
initialization_regimes <- if (identical(cvrts, "mean")) {
  c(
    "original",
    "high_spatial_dependence",
    "low_spatial_dependence",
    "random_prior_like",
    "strong_cross_dependence",
    "weak_cross_dependence"
  )
} else {
  c(
    "original",
    "dense_adjacency",
    "sparse_adjacency",
    "random_prior_like",
    "strong_cross_dependence",
    "weak_cross_dependence"
  )
}

stopifnot(run_mode %in% c("smoke", "quick", "pilot", "production"))
smoke_test <- identical(run_mode, "smoke")
quick_test <- identical(run_mode, "quick")
medium_pilot <- identical(run_mode, "pilot")

if (smoke_test) {
  warmup_iterations <- 5L
  retained_draws <- 2L
  sampling_thin <- 1L
} else if (quick_test) {
  warmup_iterations <- 1000L
  retained_draws <- 200L
  sampling_thin <- 5L
} else if (medium_pilot) {
  warmup_iterations <- 30000L
  retained_draws <- 5000L
  sampling_thin <- 10L
} else {
  warmup_iterations <- 300000L
  retained_draws <- 50000L
  sampling_thin <- 10L
}

state_consistency_tolerance <- 1e-12

sampling_iterations <- retained_draws * sampling_thin
n_workers <- min(6L, n_chains)

output_root <- Sys.getenv(
  "MCMC_MULTICHAIN_OUTPUT_ROOT",
  unset = file.path("data analysis", "multiple_chains_output")
)
output_dir <- file.path(output_root, cvrts, run_mode)
chain_output_dir <- file.path(output_dir, "chains")
log_dir <- file.path(output_dir, "logs")
diagnostics_dir <- file.path(output_dir, "diagnostics")
figures_dir <- file.path(output_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chain_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Long diagnostic tables are written to diagnostics_dir instead of flooding
# the console. Run/progress messages remain visible.
print_diagnostics_to_console <- FALSE
diagnostic_cat <- function(...) {
  if (print_diagnostics_to_console) {
    cat(...)
  }
}
diagnostic_print <- function(object, ...) {
  if (print_diagnostics_to_console) {
    print(object, ...)
  }
  invisible(object)
}

sampler_path <- normalizePath(
  file.path("data analysis", "sampler_multiple_chains.cpp"),
  mustWork = TRUE
)

stopifnot(
  n_chains == length(chain_seeds),
  n_chains == length(initialization_regimes),
  length(unique(chain_seeds)) == n_chains,
  all(chain_seeds > 0),
  warmup_iterations >= 0,
  retained_draws > 0,
  sampling_thin > 0,
  is.finite(state_consistency_tolerance),
  state_consistency_tolerance > 0
)

sampler_md5 <- unname(tools::md5sum(sampler_path))

model_arguments <- list(
  y = Y,
  X = X,
  Z1 = Z1,
  Z2 = Z2,
  Z3 = Z3,
  E = E,
  cvrts = cvrts,
  q = 4L,
  Winc = Winc,
  Minc = Minc,
  alpha = 1,
  n_atoms = 15L
)

cat("\nBaseline multiple-chain SEER run\n")
cat("Start time:", format(Sys.time()), "\n")
cat("Run mode:", run_mode, "\n")
cat("Specification:", cvrts, "\n")
cat("Run new sampling:", run_sampling, "\n")
cat("Reuse matching completed chains:", reuse_completed_chains, "\n")
cat("Write pooled boundary map:", write_pooled_boundary_map, "\n")
cat("Number of chains:", n_chains, "\n")
cat("Number of workers:", n_workers, "\n")
cat("Warm-up iterations per chain:", warmup_iterations, "\n")
cat("Sampling iterations per chain:", sampling_iterations, "\n")
cat("Retained post-warmup draws per chain:", retained_draws, "\n")
cat("Thinning:", sampling_thin, "\n")
cat("Transition structure: non-tempered baseline blocks only\n")
cat("State-consistency tolerance:", state_consistency_tolerance, "\n")
cat("Chain seeds:", paste(chain_seeds, collapse = ", "), "\n")
cat("Initialization regimes:\n")
print(data.frame(
  chain = paste0("chain_", seq_len(n_chains)),
  seed = chain_seeds,
  initialization_regime = initialization_regimes
), row.names = FALSE)
cat("Output directory:", output_dir, "\n\n")
cat("Diagnostics directory:", diagnostics_dir, "\n\n")

if (run_sampling) {
  cat("Compiling baseline C++ sampler on the main R session...\n")
  Rcpp::sourceCpp(sampler_path)
  cat("C++ compile check passed.\n")
} else {
  cat(
    "Diagnostics-only mode: existing completed chain files will be loaded.\n"
  )
}

run_one_chain <- function(chain_index) {
  chain_id <- as.integer(chain_index)
  seed <- chain_seeds[[chain_id]]
  regime <- initialization_regimes[[chain_id]]
  worker_pid <- Sys.getpid()
  chain_label <- paste0("chain_", chain_id)
  chain_file <- file.path(
    chain_output_dir,
    paste0(chain_label, "_results.rds")
  )
  error_file <- file.path(
    chain_output_dir,
    paste0(chain_label, "_error.rds")
  )
  progress_log_path <- file.path(
    log_dir,
    paste0(chain_label, "_progress.log")
  )

  if (reuse_completed_chains && file.exists(chain_file)) {
    existing_fit <- tryCatch(readRDS(chain_file), error = identity)
    existing_matches <- !inherits(existing_fit, "error") &&
      identical(existing_fit$chain_id, chain_id) &&
      identical(existing_fit$seed, seed) &&
      identical(existing_fit$settings$initialization_regime, regime) &&
      identical(existing_fit$settings$cvrts, cvrts) &&
      identical(existing_fit$settings$burn, warmup_iterations) &&
      identical(existing_fit$settings$runs, retained_draws) &&
      identical(existing_fit$settings$thin, sampling_thin) &&
      identical(existing_fit$sampler_md5, sampler_md5) &&
      isTRUE(existing_fit$state_consistency$all_checks_passed)

    if (existing_matches) {
      cat("###", chain_label, "| reusing matching completed chain\n")
      return(list(
        status = "reused",
        chain = chain_id,
        seed = seed,
        initialization_regime = regime,
        output_file = chain_file,
        progress_log_path = progress_log_path,
        elapsed_seconds = existing_fit$runtime$elapsed_seconds
      ))
    }
  }

  if (file.exists(progress_log_path)) {
    unlink(progress_log_path)
  }

  start_time <- Sys.time()
  cat(
    "\n###", chain_label,
    "| PID", worker_pid,
    "| seed", seed,
    "| regime", regime,
    "| started", format(start_time), "\n"
  )

  result <- tryCatch(
    {
      fit <- do.call(
        MADAGAR_chain,
        c(
          model_arguments,
          list(
            runs = retained_draws,
            burn = warmup_iterations,
            thin = sampling_thin,
            chain_seed = seed,
            chain_id = chain_id,
            worker_pid = worker_pid,
            initialization_regime = regime,
            progress_log_path = progress_log_path,
            print_progress = TRUE,
            state_consistency_tolerance = state_consistency_tolerance
          )
        )
      )

      finish_time <- Sys.time()
      runtime_seconds <- as.numeric(
        difftime(finish_time, start_time, units = "secs")
      )
      fit$runtime <- list(
        start_time = start_time,
        finish_time = finish_time,
        elapsed_seconds = runtime_seconds
      )
      fit$output_file <- chain_file
      fit$progress_log_path <- progress_log_path
      fit$quick_test <- quick_test
      fit$medium_pilot <- medium_pilot
      fit$run_mode <- run_mode
      fit$sampler_md5 <- sampler_md5

      saveRDS(fit, chain_file)
      if (file.exists(error_file)) {
        unlink(error_file)
      }

      cat(
        "###", chain_label,
        "| finished", format(finish_time),
        "| elapsed seconds", round(runtime_seconds), "\n"
      )

      list(
        status = "ok",
        chain = chain_id,
        seed = seed,
        initialization_regime = regime,
        output_file = chain_file,
        progress_log_path = progress_log_path,
        elapsed_seconds = runtime_seconds
      )
    },
    error = function(error) {
      finish_time <- Sys.time()
      error_result <- list(
        status = "error",
        chain = chain_id,
        seed = seed,
        initialization_regime = regime,
        message = conditionMessage(error),
        start_time = start_time,
        finish_time = finish_time,
        progress_log_path = progress_log_path
      )
      saveRDS(error_result, error_file)
      cat(
        "###", chain_label,
        "| ERROR:", conditionMessage(error), "\n"
      )
      error_result$error_file <- error_file
      error_result
    }
  )

  result
}

if (run_sampling) {
  cl <- parallel::makeCluster(n_workers, outfile = "")
  parallel::clusterSetRNGStream(cl, iseed = 20260721L)

  chain_run_results <- tryCatch(
    {
      parallel::clusterExport(
        cl,
        c(
          "sampler_path",
          "model_arguments",
          "chain_seeds",
          "initialization_regimes",
          "warmup_iterations",
          "retained_draws",
          "sampling_thin",
          "state_consistency_tolerance",
          "sampler_md5",
          "reuse_completed_chains",
          "run_mode",
          "cvrts",
          "quick_test",
          "medium_pilot",
          "chain_output_dir",
          "log_dir",
          "run_one_chain"
        ),
        envir = environment()
      )

      parallel::clusterEvalQ(cl, {
        library(Rcpp)
        library(RcppArmadillo)
        Rcpp::sourceCpp(sampler_path)
        NULL
      })

      parallel::parLapply(cl, seq_len(n_chains), run_one_chain)
    },
    finally = parallel::stopCluster(cl)
  )

  chain_run_status <- do.call(
    rbind,
    lapply(chain_run_results, function(result) {
      data.frame(
        status = result$status,
        chain = result$chain,
        seed = result$seed,
        initialization_regime = result$initialization_regime,
        output_file = if (!is.null(result$output_file)) {
          result$output_file
        } else {
          NA_character_
        },
        error_file = if (!is.null(result$error_file)) {
          result$error_file
        } else {
          NA_character_
        },
        message = if (!is.null(result$message)) {
          result$message
        } else {
          NA_character_
        },
        elapsed_seconds = if (!is.null(result$elapsed_seconds)) {
          result$elapsed_seconds
        } else {
          NA_real_
        }
      )
    })
  )
  print(chain_run_status, row.names = FALSE)

  saveRDS(
    list(
      settings = list(
        run_mode = run_mode,
        smoke_test = smoke_test,
        quick_test = quick_test,
        medium_pilot = medium_pilot,
        run_sampling = run_sampling,
        reuse_completed_chains = reuse_completed_chains,
        sampler_md5 = sampler_md5,
        n_chains = n_chains,
        chain_seeds = chain_seeds,
        initialization_regimes = initialization_regimes,
        warmup_iterations = warmup_iterations,
        retained_draws = retained_draws,
        sampling_thin = sampling_thin,
        sampling_iterations = sampling_iterations,
        transition_blocks = "non-tempered baseline blocks only",
        state_consistency_tolerance = state_consistency_tolerance,
        cvrts = cvrts
      ),
      chain_run_status = chain_run_status
    ),
    file.path(output_dir, "multiple_chains_run_status.rds")
  )

  if (any(!chain_run_status$status %in% c("ok", "reused"))) {
    stop(
      "At least one baseline chain failed. Completed chains were saved in ",
      chain_output_dir,
      "; inspect multiple_chains_run_status.rds and *_error.rds."
    )
  }

  chain_files <- chain_run_status$output_file[
    order(chain_run_status$chain)
  ]
} else {
  chain_files <- file.path(
    chain_output_dir,
    paste0("chain_", seq_len(n_chains), "_results.rds")
  )
  missing_chain_files <- chain_files[!file.exists(chain_files)]
  if (length(missing_chain_files) > 0L) {
    stop(
      "Diagnostics-only mode requires all completed chain files. Missing: ",
      paste(missing_chain_files, collapse = ", ")
    )
  }
  chain_run_status <- data.frame(
    status = "loaded",
    chain = seq_len(n_chains),
    seed = chain_seeds,
    initialization_regime = initialization_regimes,
    output_file = chain_files,
    error_file = NA_character_,
    message = NA_character_,
    elapsed_seconds = NA_real_
  )
  print(chain_run_status, row.names = FALSE)
}

chain_fits <- lapply(chain_files, readRDS)

stored_chain_ids <- vapply(chain_fits, `[[`, integer(1L), "chain_id")
stored_seeds <- vapply(chain_fits, `[[`, integer(1L), "seed")
stored_regimes <- vapply(
  chain_fits,
  function(fit) fit$settings$initialization_regime,
  character(1L)
)
retained_draw_counts <- vapply(
  chain_fits,
  function(fit) nrow(fit$beta),
  integer(1L)
)
if (!identical(stored_chain_ids, seq_len(n_chains)) ||
    !identical(stored_seeds, chain_seeds) ||
    !identical(stored_regimes, initialization_regimes)) {
  stop(
    "Saved chain metadata do not match the configured chain ids, seeds, ",
    "or initialization regimes."
  )
}
if (length(unique(retained_draw_counts)) != 1L) {
  stop("All chains must contain the same number of retained draws.")
}

state_consistency_summary <- do.call(
  rbind,
  lapply(seq_along(chain_fits), function(chain) {
    consistency <- chain_fits[[chain]]$state_consistency
    data.frame(
      chain = paste0("chain_", chain),
      checks = consistency$checks,
      tolerance = consistency$tolerance,
      maximum_F_gap_seen = consistency$maximum_F_gap_seen,
      maximum_u_mismatches_seen = consistency$maximum_u_mismatches_seen,
      maximum_W_eta_mismatches_seen = if (
        is.null(consistency$maximum_W_eta_mismatches_seen)
      ) {
        NA_integer_
      } else {
        consistency$maximum_W_eta_mismatches_seen
      },
      all_checks_passed = isTRUE(consistency$all_checks_passed)
    )
  })
)
diagnostic_cat("\nLatent-state consistency checks\n")
diagnostic_print(state_consistency_summary, row.names = FALSE)
if (!all(state_consistency_summary$all_checks_passed)) {
  stop("At least one chain failed a saved latent-state consistency check.")
}

if (!run_sampling) {
  stored_burn <- vapply(
    chain_fits,
    function(fit) fit$settings$burn,
    integer(1L)
  )
  stored_thin <- vapply(
    chain_fits,
    function(fit) fit$settings$thin,
    integer(1L)
  )
  if (length(unique(stored_burn)) != 1L ||
      length(unique(stored_thin)) != 1L) {
    stop("Saved chains do not share common warm-up and thinning settings.")
  }
  warmup_iterations <- stored_burn[[1L]]
  sampling_thin <- stored_thin[[1L]]
  retained_draws <- retained_draw_counts[[1L]]
  sampling_iterations <- retained_draws * sampling_thin
  cat(
    "Loaded", n_chains, "chains with", retained_draws,
    "retained draws per chain.\n"
  )
}

# Multiple-chain diagnostics ---------------------------------------------------

if (!requireNamespace("posterior", quietly = TRUE)) {
  stop(
    "The 'posterior' package is required for rank-normalized split-Rhat, ",
    "bulk ESS, tail ESS, and MCSE diagnostics. Run install.R first."
  )
}
diagnostic_backend <- paste0(
  "posterior ", as.character(utils::packageVersion("posterior"))
)

disease_names <- c("Lung", "Esophageal", "Larynx", "Colorectal")
n_diseases <- length(disease_names)
n_counties <- nrow(Minc)
parameters_per_disease <- ncol(chain_fits[[1]]$beta) / n_diseases
intercept_columns <- seq(
  1L,
  ncol(chain_fits[[1]]$beta),
  by = parameters_per_disease
)

beta_names <- if (parameters_per_disease == 1L) {
  paste0("centered_intercept.", disease_names)
} else {
  as.vector(t(outer(
    disease_names,
    c("intercept", "smoking", "old", "poverty"),
    paste,
    sep = "."
  )))
}

spatial_names <- paste0(
  rep(disease_names, each = n_counties),
  ".county_",
  rep(seq_len(n_counties), times = n_diseases)
)

center_chain <- function(fit) {
  centered_phi <- fit$phi
  phi_means <- matrix(
    0,
    nrow = nrow(fit$phi),
    ncol = n_diseases
  )

  for (disease in seq_len(n_diseases)) {
    indices <- ((disease - 1L) * n_counties + 1L):
      (disease * n_counties)
    phi_means[, disease] <- rowMeans(fit$phi[, indices, drop = FALSE])
    centered_phi[, indices] <-
      fit$phi[, indices, drop = FALSE] - phi_means[, disease]
  }

  centered_beta <- fit$beta
  centered_beta[, intercept_columns] <-
    fit$beta[, intercept_columns, drop = FALSE] + phi_means

  list(beta = centered_beta, phi = centered_phi)
}

build_centered_beta <- function(fit) {
  result <- center_chain(fit)$beta
  colnames(result) <- beta_names
  result
}

build_centered_phi <- function(fit) {
  result <- center_chain(fit)$phi
  colnames(result) <- paste0("centered_phi.", spatial_names)
  result
}

build_log_risk <- function(fit) {
  result <- fit$beta %*% t(X) + fit$phi
  colnames(result) <- paste0("log_risk.", spatial_names)
  result
}

build_core_parameters <- function(fit) {
  result <- cbind(tau = fit$tau, fit$rho)
  colnames(result) <- c("tau", paste0("rho.", disease_names))

  if (cvrts != "mean") {
    result <- cbind(result, fit$eta)
    colnames(result)[
      (ncol(result) - ncol(fit$eta) + 1L):ncol(result)
    ] <- paste0(
      "eta.",
      rep(c("smoking", "old", "poverty"), each = n_diseases),
      ".",
      rep(disease_names, times = 3L)
    )
  }

  result
}

A_vector_to_matrix <- function(A_vector) {
  A_matrix <- matrix(0, n_diseases, n_diseases)
  index <- 1L

  for (row in seq_len(n_diseases)) {
    for (column in seq_len(row)) {
      A_matrix[row, column] <- if (row == column) {
        exp(A_vector[index])
      } else {
        A_vector[index]
      }
      index <- index + 1L
    }
  }

  A_matrix
}

covariance_indices <- which(
  lower.tri(matrix(0, n_diseases, n_diseases), diag = TRUE),
  arr.ind = TRUE
)
correlation_indices <- which(
  lower.tri(matrix(0, n_diseases, n_diseases), diag = FALSE),
  arr.ind = TRUE
)

build_dependence_parameters <- function(fit) {
  result <- matrix(
    NA_real_,
    nrow = nrow(fit$A),
    ncol = nrow(covariance_indices) + nrow(correlation_indices)
  )

  for (draw in seq_len(nrow(fit$A))) {
    A_matrix <- A_vector_to_matrix(fit$A[draw, ])
    Sigma <- A_matrix %*% t(A_matrix)
    correlation <- cov2cor(Sigma)
    result[draw, ] <- c(
      Sigma[covariance_indices],
      correlation[correlation_indices]
    )
  }

  colnames(result) <- c(
    paste0(
      "Sigma.",
      covariance_indices[, "row"],
      ".",
      covariance_indices[, "col"]
    ),
    paste0(
      "correlation.",
      correlation_indices[, "row"],
      ".",
      correlation_indices[, "col"]
    )
  )
  result
}

build_adjacency_cardinality <- function(fit) {
  result <- fit$W_cardinality
  colnames(result) <- paste0("adjacency_edges.", disease_names)
  result
}

build_boundary_count <- function(fit) {
  result <- matrix(
    NA_integer_,
    nrow = nrow(fit$phi),
    ncol = n_diseases
  )

  for (disease in seq_len(n_diseases)) {
    indices <- ((disease - 1L) * n_counties + 1L):
      (disease * n_counties)
    phi <- fit$phi[, indices, drop = FALSE]
    phi <- phi[, order(final_perm), drop = FALSE]
    result[, disease] <- rowSums(
      phi[, neighbor_list0[, 1L], drop = FALSE] !=
        phi[, neighbor_list0[, 2L], drop = FALSE]
    )
  }

  colnames(result) <- paste0("boundary_count.", disease_names)
  result
}

summarise_draw_array <- function(draw_array) {
  as.data.frame(
    posterior::summarise_draws(
      posterior::as_draws_array(draw_array),
      "rhat",
      "ess_bulk",
      "ess_tail",
      "mcse_mean"
    )
  )
}

summarise_draw_matrix <- function(draw_matrix, chain_label) {
  draw_matrix <- as.matrix(draw_matrix)
  variable_names <- colnames(draw_matrix)
  quantiles <- matrixStats::colQuantiles(
    draw_matrix,
    probs = c(0.025, 0.50, 0.975),
    drop = FALSE
  )

  data.frame(
    chain = chain_label,
    variable = variable_names,
    mean = colMeans(draw_matrix),
    sd = matrixStats::colSds(draw_matrix),
    median = quantiles[, 2L],
    q025 = quantiles[, 1L],
    q975 = quantiles[, 3L],
    row.names = NULL
  )
}

summarise_pooled_draw_array <- function(draw_array) {
  variable_names <- dimnames(draw_array)[[3L]]
  pooled_statistics <- vapply(
    seq_along(variable_names),
    function(variable) {
      values <- as.vector(draw_array[, , variable])
      quantiles <- as.numeric(
        stats::quantile(values, c(0.025, 0.50, 0.975), names = FALSE)
      )
      c(
        mean = mean(values),
        sd = stats::sd(values),
        median = quantiles[[2L]],
        q025 = quantiles[[1L]],
        q975 = quantiles[[3L]]
      )
    },
    numeric(5L)
  )

  data.frame(
    chain = "pooled",
    variable = variable_names,
    mean = pooled_statistics[1L, ],
    sd = pooled_statistics[2L, ],
    median = pooled_statistics[3L, ],
    q025 = pooled_statistics[4L, ],
    q975 = pooled_statistics[5L, ],
    row.names = NULL
  )
}

summarise_posterior_agreement <- function(chain_summaries, pooled_summary) {
  do.call(
    rbind,
    lapply(seq_len(nrow(pooled_summary)), function(index) {
      variable <- pooled_summary$variable[[index]]
      summaries <- chain_summaries[
        chain_summaries$variable == variable,
        ,
        drop = FALSE
      ]
      pooled_sd <- pooled_summary$sd[[index]]
      mean_range <- diff(range(summaries$mean))
      interval_overlap_lower <- max(summaries$q025)
      interval_overlap_upper <- min(summaries$q975)

      data.frame(
        variable = variable,
        min_chain_mean = min(summaries$mean),
        max_chain_mean = max(summaries$mean),
        between_chain_mean_range = mean_range,
        standardized_mean_range = if (is.finite(pooled_sd) && pooled_sd > 0) {
          mean_range / pooled_sd
        } else {
          NA_real_
        },
        between_chain_median_range = diff(range(summaries$median)),
        between_chain_q025_range = diff(range(summaries$q025)),
        between_chain_q975_range = diff(range(summaries$q975)),
        max_abs_chain_mean_minus_pooled = max(
          abs(summaries$mean - pooled_summary$mean[[index]])
        ),
        credible_intervals_all_overlap =
          interval_overlap_lower <= interval_overlap_upper,
        common_interval_lower = interval_overlap_lower,
        common_interval_upper = interval_overlap_upper,
        row.names = NULL
      )
    })
  )
}

diagnose_draw_group <- function(group_name, builder) {
  first_chain <- builder(chain_fits[[1L]])
  variable_names <- colnames(first_chain)
  draw_array <- array(
    NA_real_,
    dim = c(
      nrow(first_chain),
      n_chains,
      ncol(first_chain)
    ),
    dimnames = list(
      iteration = NULL,
      chain = paste0("chain_", seq_len(n_chains)),
      variable = variable_names
    )
  )
  draw_array[, 1L, ] <- first_chain

  if (n_chains > 1L) {
    for (chain in 2:n_chains) {
      draw_array[, chain, ] <- builder(chain_fits[[chain]])
    }
  }

  chain_posterior_summaries <- do.call(
    rbind,
    lapply(seq_len(n_chains), function(chain) {
      chain_matrix <- matrix(
        draw_array[, chain, ],
        nrow = dim(draw_array)[1L],
        ncol = dim(draw_array)[3L],
        dimnames = list(NULL, variable_names)
      )
      summarise_draw_matrix(chain_matrix, paste0("chain_", chain))
    })
  )
  pooled_posterior_summary <- summarise_pooled_draw_array(draw_array)
  posterior_summary <- rbind(
    chain_posterior_summaries,
    pooled_posterior_summary
  )
  posterior_summary$group <- group_name
  posterior_summary <- posterior_summary[
    ,
    c("group", "chain", "variable", "mean", "sd", "median", "q025", "q975")
  ]
  posterior_agreement <- summarise_posterior_agreement(
    chain_posterior_summaries,
    pooled_posterior_summary
  )
  posterior_agreement$group <- group_name
  posterior_agreement <- posterior_agreement[
    ,
    c("group", setdiff(names(posterior_agreement), "group"))
  ]

  diagnostic_summary <- summarise_draw_array(draw_array)
  diagnostic_summary$relative_mcse_mean <-
    diagnostic_summary$mcse_mean /
    pooled_posterior_summary$sd[
      match(diagnostic_summary$variable, pooled_posterior_summary$variable)
    ]
  diagnostic_summary$relative_mcse_mean[
    !is.finite(diagnostic_summary$relative_mcse_mean)
  ] <- NA_real_
  finite_rhat <- is.finite(diagnostic_summary$rhat)
  finite_ess <- is.finite(diagnostic_summary$ess_bulk)
  finite_tail <- is.finite(diagnostic_summary$ess_tail)
  finite_mcse <- is.finite(diagnostic_summary$mcse_mean)
  finite_relative_mcse <- is.finite(diagnostic_summary$relative_mcse_mean)

  overview <- data.frame(
    group = group_name,
    variables = sum(finite_rhat),
    median_rhat = if (any(finite_rhat)) {
      median(diagnostic_summary$rhat[finite_rhat])
    } else {
      NA_real_
    },
    max_rhat = if (any(finite_rhat)) {
      max(diagnostic_summary$rhat[finite_rhat])
    } else {
      NA_real_
    },
    rhat_above_1.01 = sum(
      diagnostic_summary$rhat[finite_rhat] > 1.01
    ),
    rhat_above_1.05 = sum(
      diagnostic_summary$rhat[finite_rhat] > 1.05
    ),
    min_bulk_ess = if (any(finite_ess)) {
      min(diagnostic_summary$ess_bulk[finite_ess])
    } else {
      NA_real_
    },
    min_tail_ess = if (any(finite_tail)) {
      min(diagnostic_summary$ess_tail[finite_tail])
    } else {
      NA_real_
    },
    max_mcse_mean = if (any(finite_mcse)) {
      max(diagnostic_summary$mcse_mean[finite_mcse])
    } else {
      NA_real_
    },
    max_relative_mcse_mean = if (any(finite_relative_mcse)) {
      max(diagnostic_summary$relative_mcse_mean[finite_relative_mcse])
    } else {
      NA_real_
    }
  )

  worst <- diagnostic_summary[
    order(diagnostic_summary$rhat, decreasing = TRUE),
    c(
      "variable", "rhat", "ess_bulk", "ess_tail",
      "mcse_mean", "relative_mcse_mean"
    )
  ]
  worst <- head(worst[is.finite(worst$rhat), ], 10L)

  agreement_worst <- head(
    posterior_agreement[
      order(
        posterior_agreement$standardized_mean_range,
        decreasing = TRUE,
        na.last = TRUE
      ),
    ],
    10L
  )

  diagnostic_cat("\n", group_name, "\n", sep = "")
  diagnostic_print(overview, row.names = FALSE)
  diagnostic_print(worst, row.names = FALSE)
  diagnostic_cat("Largest standardized between-chain posterior-mean ranges\n")
  diagnostic_print(agreement_worst, row.names = FALSE)

  rm(draw_array, first_chain)
  invisible(gc(FALSE))

  list(
    overview = overview,
    summary = diagnostic_summary,
    worst = worst,
    posterior_summary = posterior_summary,
    posterior_agreement = posterior_agreement,
    posterior_agreement_worst = agreement_worst
  )
}

diagnostic_results <- list(
  core = diagnose_draw_group(
    "Core parameters",
    build_core_parameters
  ),
  centered_beta = diagnose_draw_group(
    "Centered regression coefficients",
    build_centered_beta
  ),
  centered_phi = diagnose_draw_group(
    "Centered spatial effects",
    build_centered_phi
  ),
  log_risk = diagnose_draw_group(
    "Fitted log-risks (equivalently, relative risks)",
    build_log_risk
  ),
  dependence = diagnose_draw_group(
    "Disease covariance and correlation",
    build_dependence_parameters
  ),
  adjacency_cardinality = diagnose_draw_group(
    "Adjacency cardinality",
    build_adjacency_cardinality
  ),
  boundary_count = diagnose_draw_group(
    "Draw-level boundary counts",
    build_boundary_count
  )
)

diagnostic_overview <- do.call(
  rbind,
  lapply(diagnostic_results, `[[`, "overview")
)
diagnostic_print(diagnostic_overview, row.names = FALSE)

posterior_summary_all_groups <- do.call(
  rbind,
  lapply(diagnostic_results, `[[`, "posterior_summary")
)
posterior_agreement_all_groups <- do.call(
  rbind,
  lapply(diagnostic_results, `[[`, "posterior_agreement")
)
cardinality_posterior_summary <-
  diagnostic_results$adjacency_cardinality$posterior_summary
draw_boundary_count_posterior_summary <-
  diagnostic_results$boundary_count$posterior_summary

diagnostic_cat("\nCross-chain posterior-summary agreement overview\n")
posterior_agreement_overview <- do.call(
  rbind,
  lapply(split(
    posterior_agreement_all_groups,
    posterior_agreement_all_groups$group
  ), function(group) {
    data.frame(
      group = group$group[[1L]],
      variables = nrow(group),
      max_standardized_mean_range = if (any(
        is.finite(group$standardized_mean_range)
      )) {
        max(
          group$standardized_mean_range[
            is.finite(group$standardized_mean_range)
          ]
        )
      } else {
        NA_real_
      },
      credible_intervals_all_overlap = sum(
        group$credible_intervals_all_overlap,
        na.rm = TRUE
      ),
      credible_intervals_not_all_overlap = sum(
        !group$credible_intervals_all_overlap,
        na.rm = TRUE
      )
    )
  })
)
diagnostic_print(posterior_agreement_overview, row.names = FALSE)

fitted_risk_pairwise_agreement <- do.call(
  rbind,
  lapply(disease_names, function(disease) {
    do.call(rbind, combn(
      paste0("chain_", seq_len(n_chains)),
      2L,
      function(pair) {
        first <- posterior_summary_all_groups[
          posterior_summary_all_groups$group ==
            "Fitted log-risks (equivalently, relative risks)" &
            posterior_summary_all_groups$chain == pair[[1L]] &
            startsWith(
              posterior_summary_all_groups$variable,
              paste0("log_risk.", disease, ".")
            ),
          c("variable", "mean"),
          drop = FALSE
        ]
        second <- posterior_summary_all_groups[
          posterior_summary_all_groups$group ==
            "Fitted log-risks (equivalently, relative risks)" &
            posterior_summary_all_groups$chain == pair[[2L]] &
            startsWith(
              posterior_summary_all_groups$variable,
              paste0("log_risk.", disease, ".")
            ),
          c("variable", "mean"),
          drop = FALSE
        ]
        second <- second[match(first$variable, second$variable), ]
        difference <- first$mean - second$mean
        data.frame(
          disease = disease,
          chain_1 = pair[[1L]],
          chain_2 = pair[[2L]],
          correlation = if (stats::sd(first$mean) > 0 &&
                            stats::sd(second$mean) > 0) {
            stats::cor(first$mean, second$mean)
          } else {
            NA_real_
          },
          mean_absolute_difference = mean(abs(difference)),
          root_mean_squared_difference = sqrt(mean(difference^2)),
          max_absolute_difference = max(abs(difference))
        )
      },
      simplify = FALSE
    ))
  })
)

fitted_risk_pairwise_summary <- do.call(
  rbind,
  lapply(disease_names, function(disease) {
    agreement <- fitted_risk_pairwise_agreement[
      fitted_risk_pairwise_agreement$disease == disease,
      ,
      drop = FALSE
    ]
    finite_correlations <- is.finite(agreement$correlation)
    data.frame(
      disease = disease,
      min_pairwise_risk_map_correlation = if (any(finite_correlations)) {
        min(agreement$correlation[finite_correlations])
      } else {
        NA_real_
      },
      mean_pairwise_risk_map_correlation = if (any(finite_correlations)) {
        mean(agreement$correlation[finite_correlations])
      } else {
        NA_real_
      },
      max_pairwise_mean_absolute_difference = max(
        agreement$mean_absolute_difference
      ),
      max_county_absolute_difference = max(
        agreement$max_absolute_difference
      )
    )
  })
)

diagnostic_cat("\nAgreement of chain-specific fitted log-risk patterns\n")
diagnostic_print(fitted_risk_pairwise_summary, row.names = FALSE)

acceptance_summary <- do.call(
  rbind,
  lapply(seq_along(chain_fits), function(chain) {
    acceptance <- chain_fits[[chain]]$acceptance
    data.frame(
      chain = paste0("chain_", chain),
      beta = acceptance$beta_rate,
      theta = acceptance$theta_rate,
      r = acceptance$r_rate,
      V = acceptance$V_rate,
      rho = acceptance$rho_rate,
      eta = acceptance$eta_rate,
      A = acceptance$A_rate
    )
  })
)
diagnostic_cat("\nBaseline proposal acceptance rates\n")
diagnostic_print(acceptance_summary, row.names = FALSE)

deterministic_state_proposal_summary <- do.call(
  rbind,
  lapply(seq_along(chain_fits), function(chain) {
    acceptance <- chain_fits[[chain]]$acceptance
    data.frame(
      chain = paste0("chain_", chain),
      block = c("rho", "eta", "A"),
      acceptance_rate = c(
        acceptance$rho_rate,
        acceptance$eta_rate,
        acceptance$A_rate
      ),
      mean_proposed_allocation_changes = c(
        acceptance$rho_mean_proposed_allocation_changes,
        acceptance$eta_mean_proposed_allocation_changes,
        acceptance$A_mean_proposed_allocation_changes
      ),
      mean_accepted_allocation_changes = c(
        acceptance$rho_mean_accepted_allocation_changes,
        acceptance$eta_mean_accepted_allocation_changes,
        acceptance$A_mean_accepted_allocation_changes
      )
    )
  })
)
diagnostic_cat("\nAllocation changes induced by joint rho, eta, and A proposals\n")
diagnostic_print(deterministic_state_proposal_summary, row.names = FALSE)


# Boundary stability diagnostics ----------------------------------------------

boundary_probabilities <- array(
  NA_real_,
  dim = c(
    nrow(neighbor_list0),
    n_diseases,
    n_chains
  ),
  dimnames = list(
    boundary = neighbor_name,
    disease = disease_names,
    chain = paste0("chain_", seq_len(n_chains))
  )
)

for (chain in seq_len(n_chains)) {
  for (disease in seq_len(n_diseases)) {
    indices <- ((disease - 1L) * n_counties + 1L):
      (disease * n_counties)
    phi <- chain_fits[[chain]]$phi[, indices, drop = FALSE]
    phi <- phi[, order(final_perm), drop = FALSE]
    boundary_probabilities[, disease, chain] <- colMeans(
      phi[, neighbor_list0[, 1L], drop = FALSE] !=
        phi[, neighbor_list0[, 2L], drop = FALSE]
    )
  }
}

boundary_probability_range <- apply(
  boundary_probabilities,
  c(1L, 2L),
  function(x) diff(range(x))
)
boundary_probability_mean <- apply(
  boundary_probabilities,
  c(1L, 2L),
  mean
)
boundary_probability_mean_pairwise_absolute_difference <- apply(
  boundary_probabilities,
  c(1L, 2L),
  function(probabilities) {
    mean(abs(outer(probabilities, probabilities, "-"))[upper.tri(
      matrix(0, n_chains, n_chains)
    )])
  }
)

boundary_probability_range_summary <- do.call(
  rbind,
  lapply(seq_len(n_diseases), function(disease) {
    ranges <- boundary_probability_range[, disease]
    data.frame(
      disease = disease_names[[disease]],
      median_between_chain_range = median(ranges),
      mean_between_chain_range = mean(ranges),
      q90_between_chain_range = as.numeric(
        quantile(ranges, 0.90)
      ),
      max_between_chain_range = max(ranges),
      median_pairwise_absolute_difference = median(
        boundary_probability_mean_pairwise_absolute_difference[, disease]
      ),
      mean_pairwise_absolute_difference = mean(
        boundary_probability_mean_pairwise_absolute_difference[, disease]
      ),
      max_pairwise_absolute_difference = max(
        boundary_probability_mean_pairwise_absolute_difference[, disease]
      ),
      boundaries_range_ge_0.25 = sum(ranges >= 0.25),
      boundaries_range_ge_0.50 = sum(ranges >= 0.50),
      boundaries_range_ge_0.75 = sum(ranges >= 0.75)
    )
  })
)

boundary_agreement <- data.frame(
  boundary = rep(
    rownames(boundary_probability_range),
    times = n_diseases
  ),
  disease = rep(
    colnames(boundary_probability_range),
    each = nrow(boundary_probability_range)
  ),
  mean_probability = as.vector(boundary_probability_mean),
  between_chain_range = as.vector(boundary_probability_range),
  mean_pairwise_absolute_difference = as.vector(
    boundary_probability_mean_pairwise_absolute_difference
  )
)
boundary_agreement <- boundary_agreement[
  order(
    boundary_agreement$between_chain_range,
    decreasing = TRUE
  ),
]

safe_probability_correlation <- function(first, second) {
  if (stats::sd(first) == 0 || stats::sd(second) == 0) {
    return(if (isTRUE(all.equal(first, second))) 1 else NA_real_)
  }
  stats::cor(first, second)
}

boundary_probability_pairwise_comparison <- do.call(
  rbind,
  lapply(seq_len(n_diseases), function(disease) {
    do.call(rbind, combn(
      seq_len(n_chains),
      2L,
      function(pair) {
        first <- boundary_probabilities[, disease, pair[[1L]]]
        second <- boundary_probabilities[, disease, pair[[2L]]]
        absolute_difference <- abs(first - second)
        data.frame(
          disease = disease_names[[disease]],
          chain_1 = paste0("chain_", pair[[1L]]),
          chain_2 = paste0("chain_", pair[[2L]]),
          correlation = safe_probability_correlation(first, second),
          mean_absolute_difference = mean(absolute_difference),
          root_mean_squared_difference = sqrt(mean((first - second)^2)),
          max_absolute_difference = max(absolute_difference)
        )
      },
      simplify = FALSE
    ))
  })
)

boundary_probability_correlation_matrices <- setNames(
  lapply(seq_len(n_diseases), function(disease) {
    result <- matrix(
      NA_real_,
      n_chains,
      n_chains,
      dimnames = list(
        paste0("chain_", seq_len(n_chains)),
        paste0("chain_", seq_len(n_chains))
      )
    )
    diag(result) <- 1
    for (first in seq_len(n_chains - 1L)) {
      for (second in (first + 1L):n_chains) {
        result[first, second] <- result[second, first] <-
          safe_probability_correlation(
            boundary_probabilities[, disease, first],
            boundary_probabilities[, disease, second]
          )
      }
    }
    result
  }),
  disease_names
)

boundary_probability_pairwise_summary <- do.call(
  rbind,
  lapply(disease_names, function(disease) {
    comparison <- boundary_probability_pairwise_comparison[
      boundary_probability_pairwise_comparison$disease == disease,
      ,
      drop = FALSE
    ]
    data.frame(
      disease = disease,
      min_pairwise_correlation = if (any(is.finite(comparison$correlation))) {
        min(comparison$correlation[is.finite(comparison$correlation)])
      } else {
        NA_real_
      },
      mean_pairwise_correlation = if (any(is.finite(comparison$correlation))) {
        mean(comparison$correlation[is.finite(comparison$correlation)])
      } else {
        NA_real_
      },
      max_pairwise_mean_absolute_difference = max(
        comparison$mean_absolute_difference
      ),
      max_pairwise_edge_difference = max(comparison$max_absolute_difference)
    )
  })
)

boundary_fdr_selection <- function(probabilities, alpha = 0.05) {
  probabilities <- probabilities[!is.na(probabilities)]
  if (length(probabilities) == 0L) {
    return(list(
      selected_edges = 0L,
      threshold = NA_real_,
      estimated_fdr = NA_real_,
      cutoff_ties = 0L,
      mean_selected_probability = NA_real_,
      min_selected_probability = NA_real_,
      max_probability = NA_real_,
      selected_boundaries = character(0)
    ))
  }

  # Evaluate unique probability thresholds so all edges tied at the cutoff are
  # treated identically, matching the threshold-based manuscript FDR rule.
  thresholds <- sort(unique(probabilities), decreasing = TRUE)
  selected_counts <- vapply(
    thresholds,
    function(threshold) sum(probabilities >= threshold),
    integer(1L)
  )
  estimated_fdr <- vapply(
    thresholds,
    function(threshold) mean(1 - probabilities[probabilities >= threshold]),
    numeric(1L)
  )
  valid_thresholds <- which(estimated_fdr <= alpha)

  if (length(valid_thresholds) == 0L) {
    return(list(
      selected_edges = 0L,
      threshold = NA_real_,
      estimated_fdr = NA_real_,
      cutoff_ties = 0L,
      mean_selected_probability = NA_real_,
      min_selected_probability = NA_real_,
      max_probability = max(probabilities),
      selected_boundaries = character(0)
    ))
  }

  selected_threshold_index <- valid_thresholds[
    which.max(selected_counts[valid_thresholds])
  ]
  threshold <- thresholds[[selected_threshold_index]]
  selected_probabilities <- sort(
    probabilities[probabilities >= threshold],
    decreasing = TRUE
  )

  list(
    selected_edges = length(selected_probabilities),
    threshold = threshold,
    estimated_fdr = estimated_fdr[[selected_threshold_index]],
    cutoff_ties = sum(probabilities == threshold),
    mean_selected_probability = mean(selected_probabilities),
    min_selected_probability = min(selected_probabilities),
    max_probability = max(probabilities),
    selected_boundaries = names(selected_probabilities)
  )
}

selected_boundary_sets <- setNames(
  vector("list", n_diseases),
  disease_names
)
for (disease in disease_names) {
  selected_boundary_sets[[disease]] <- setNames(
    vector("list", n_chains),
    paste0("chain_", seq_len(n_chains))
  )
}

cardinality_summary_by_chain <- do.call(
  rbind,
  lapply(seq_len(n_chains), function(chain) {
    do.call(
      rbind,
      lapply(seq_len(n_diseases), function(disease) {
        adjacency_cardinality <-
          chain_fits[[chain]]$W_cardinality[, disease]
        cardinality_quantiles <- as.numeric(
          quantile(
            adjacency_cardinality,
            probs = c(0.05, 0.50, 0.95)
          )
        )

        data.frame(
          chain = paste0("chain_", chain),
          disease = disease_names[[disease]],
          mean_adjacency_edges = mean(adjacency_cardinality),
          sd_adjacency_edges = sd(adjacency_cardinality),
          q05_adjacency_edges = cardinality_quantiles[[1L]],
          median_adjacency_edges = cardinality_quantiles[[2L]],
          q95_adjacency_edges = cardinality_quantiles[[3L]],
          min_adjacency_edges = min(adjacency_cardinality),
          max_adjacency_edges = max(adjacency_cardinality)
        )
      })
    )
  })
)

boundary_count_summary <- do.call(
  rbind,
  lapply(seq_len(n_chains), function(chain) {
    do.call(
      rbind,
      lapply(seq_len(n_diseases), function(disease) {
        probabilities <- boundary_probabilities[
          ,
          disease,
          chain
        ]
        names(probabilities) <- dimnames(boundary_probabilities)$boundary
        selection <- boundary_fdr_selection(
          probabilities,
          alpha = 0.05
        )
        selected_boundary_sets[[disease_names[[disease]]]][[
          paste0("chain_", chain)
        ]] <<- selection$selected_boundaries

        adjacency_cardinality <-
          chain_fits[[chain]]$W_cardinality[, disease]

        data.frame(
          chain = paste0("chain_", chain),
          disease = disease_names[[disease]],
          selected_edges = selection$selected_edges,
          threshold = selection$threshold,
          estimated_fdr = selection$estimated_fdr,
          cutoff_ties = selection$cutoff_ties,
          mean_selected_probability =
            selection$mean_selected_probability,
          min_selected_probability =
            selection$min_selected_probability,
          max_probability = selection$max_probability,
          boundaries_prob_ge_0.95 = sum(probabilities >= 0.95),
          boundaries_prob_ge_0.90 = sum(probabilities >= 0.90),
          mean_adjacency_edges = mean(adjacency_cardinality),
          sd_adjacency_edges = sd(adjacency_cardinality)
        )
      })
    )
  })
)

pooled_selected_boundary_sets <- setNames(
  vector("list", n_diseases),
  disease_names
)
pooled_boundary_count_summary <- do.call(
  rbind,
  lapply(seq_len(n_diseases), function(disease) {
    probabilities <- boundary_probability_mean[, disease]
    names(probabilities) <- rownames(boundary_probability_mean)
    selection <- boundary_fdr_selection(
      probabilities,
      alpha = 0.05
    )
    pooled_selected_boundary_sets[[disease_names[[disease]]]] <<-
      selection$selected_boundaries

    data.frame(
      chain = "pooled_mean",
      disease = disease_names[[disease]],
      selected_edges = selection$selected_edges,
      threshold = selection$threshold,
      estimated_fdr = selection$estimated_fdr,
      cutoff_ties = selection$cutoff_ties,
      mean_selected_probability =
        selection$mean_selected_probability,
      min_selected_probability =
        selection$min_selected_probability,
      max_probability = selection$max_probability,
      boundaries_prob_ge_0.95 = sum(probabilities >= 0.95),
      boundaries_prob_ge_0.90 = sum(probabilities >= 0.90),
      mean_adjacency_edges = NA_real_,
      sd_adjacency_edges = NA_real_
    )
  })
)

boundary_count_ranges <- do.call(
  rbind,
  lapply(disease_names, function(disease) {
    chain_counts <- boundary_count_summary$selected_edges[
      boundary_count_summary$disease == disease
    ]
    pooled_count <- pooled_boundary_count_summary$selected_edges[
      pooled_boundary_count_summary$disease == disease
    ]

    data.frame(
      disease = disease,
      min_chain_selected_edges = min(chain_counts),
      median_chain_selected_edges = median(chain_counts),
      max_chain_selected_edges = max(chain_counts),
      pooled_selected_edges = pooled_count,
      max_minus_min_chain = max(chain_counts) - min(chain_counts),
      pooled_minus_median_chain =
        pooled_count - median(chain_counts)
    )
  })
)

boundary_prevalence_ranking <- rbind(
  boundary_count_summary[, c("chain", "disease", "selected_edges")],
  pooled_boundary_count_summary[, c("chain", "disease", "selected_edges")]
)
boundary_prevalence_ranking$boundary_prevalence_rank <- ave(
  -boundary_prevalence_ranking$selected_edges,
  boundary_prevalence_ranking$chain,
  FUN = function(values) rank(values, ties.method = "min")
)

boundary_prevalence_rank_agreement <- do.call(
  rbind,
  combn(
    unique(boundary_prevalence_ranking$chain),
    2L,
    function(pair) {
      first <- boundary_prevalence_ranking[
        boundary_prevalence_ranking$chain == pair[[1L]],
        c("disease", "boundary_prevalence_rank")
      ]
      second <- boundary_prevalence_ranking[
        boundary_prevalence_ranking$chain == pair[[2L]],
        c("disease", "boundary_prevalence_rank")
      ]
      second <- second[match(first$disease, second$disease), ]
      data.frame(
        analysis_1 = pair[[1L]],
        analysis_2 = pair[[2L]],
        identical_ranking = identical(
          first$boundary_prevalence_rank,
          second$boundary_prevalence_rank
        ),
        spearman_rank_correlation = suppressWarnings(stats::cor(
          first$boundary_prevalence_rank,
          second$boundary_prevalence_rank,
          method = "spearman"
        ))
      )
    },
    simplify = FALSE
  )
)

boundary_jaccard <- function(first, second) {
  union_size <- length(union(first, second))
  if (union_size == 0L) {
    return(NA_real_)
  }
  length(intersect(first, second)) / union_size
}

summarize_jaccard <- function(values, summary_function) {
  if (all(is.na(values))) {
    return(NA_real_)
  }
  summary_function(values, na.rm = TRUE)
}

boundary_pairwise_jaccard <- do.call(
  rbind,
  lapply(disease_names, function(disease) {
    chain_names <- names(selected_boundary_sets[[disease]])
    do.call(rbind, combn(
      chain_names,
      2L,
      function(pair) {
        first_set <- selected_boundary_sets[[disease]][[pair[[1L]]]]
        second_set <- selected_boundary_sets[[disease]][[pair[[2L]]]]
        data.frame(
          disease = disease,
          chain_1 = pair[[1L]],
          chain_2 = pair[[2L]],
          jaccard = boundary_jaccard(first_set, second_set),
          common_selected_edges =
            length(intersect(first_set, second_set)),
          union_selected_edges =
            length(union(first_set, second_set))
        )
      },
      simplify = FALSE
    ))
  })
)

boundary_pairwise_jaccard_summary <- do.call(
  rbind,
  lapply(disease_names, function(disease) {
    pairwise_jaccard <- boundary_pairwise_jaccard$jaccard[
      boundary_pairwise_jaccard$disease == disease
    ]

    data.frame(
      disease = disease,
      min_pairwise_jaccard =
        summarize_jaccard(pairwise_jaccard, min),
      median_pairwise_jaccard =
        summarize_jaccard(pairwise_jaccard, median),
      mean_pairwise_jaccard =
        summarize_jaccard(pairwise_jaccard, mean),
      max_pairwise_jaccard =
        summarize_jaccard(pairwise_jaccard, max)
    )
  })
)

boundary_jaccard_matrices <- setNames(
  lapply(disease_names, function(disease) {
    chain_names <- names(selected_boundary_sets[[disease]])
    result <- matrix(
      NA_real_,
      n_chains,
      n_chains,
      dimnames = list(chain_names, chain_names)
    )
    diag(result) <- 1
    disease_pairs <- boundary_pairwise_jaccard[
      boundary_pairwise_jaccard$disease == disease,
      ,
      drop = FALSE
    ]
    for (row in seq_len(nrow(disease_pairs))) {
      first <- disease_pairs$chain_1[[row]]
      second <- disease_pairs$chain_2[[row]]
      result[first, second] <- result[second, first] <-
        disease_pairs$jaccard[[row]]
    }
    result
  }),
  disease_names
)

boundary_selection_incidence <- do.call(
  rbind,
  lapply(disease_names, function(disease) {
    incidence <- vapply(
      paste0("chain_", seq_len(n_chains)),
      function(chain) {
        neighbor_name %in% selected_boundary_sets[[disease]][[chain]]
      },
      logical(length(neighbor_name))
    )
    colnames(incidence) <- paste0("chain_", seq_len(n_chains))
    selected_by_chains <- rowSums(incidence)
    disease_index <- match(disease, disease_names)

    data.frame(
      boundary = neighbor_name,
      disease = disease,
      incidence,
      selected_by_chains = selected_by_chains,
      selected_by_all_chains = selected_by_chains == n_chains,
      selected_by_majority = selected_by_chains >=
        (floor(n_chains / 2L) + 1L),
      selected_by_one_chain = selected_by_chains == 1L,
      pooled_probability = boundary_probability_mean[, disease_index],
      pooled_selected = neighbor_name %in%
        pooled_selected_boundary_sets[[disease]],
      row.names = NULL,
      check.names = FALSE
    )
  })
)

boundary_selection_frequency_summary <- do.call(
  rbind,
  lapply(disease_names, function(disease) {
    incidence <- boundary_selection_incidence[
      boundary_selection_incidence$disease == disease,
      ,
      drop = FALSE
    ]
    data.frame(
      disease = disease,
      boundaries_selected_by_all_chains = sum(
        incidence$selected_by_all_chains
      ),
      boundaries_selected_by_majority = sum(
        incidence$selected_by_majority
      ),
      boundaries_selected_by_one_chain = sum(
        incidence$selected_by_one_chain
      ),
      boundaries_selected_by_no_chain = sum(
        incidence$selected_by_chains == 0L
      ),
      pooled_boundaries_selected = sum(incidence$pooled_selected)
    )
  })
)

boundary_cutoff_sensitivity <- do.call(
  rbind,
  lapply(disease_names, function(disease) {
    pooled_row <- pooled_boundary_count_summary[
      pooled_boundary_count_summary$disease == disease,
      ,
      drop = FALSE
    ]
    incidence <- boundary_selection_incidence[
      boundary_selection_incidence$disease == disease,
      ,
      drop = FALSE
    ]
    threshold <- pooled_row$threshold[[1L]]
    incidence$pooled_fdr_threshold <- threshold
    incidence$absolute_distance_from_threshold <- if (is.finite(threshold)) {
      abs(incidence$pooled_probability - threshold)
    } else {
      NA_real_
    }
    if (is.finite(threshold)) {
      incidence <- incidence[
        order(incidence$absolute_distance_from_threshold),
        ,
        drop = FALSE
      ]
    } else {
      incidence <- incidence[
        order(incidence$pooled_probability, decreasing = TRUE),
        ,
        drop = FALSE
      ]
    }
    head(incidence, 20L)
  })
)

boundary_agreement_wide <- boundary_agreement
for (chain in seq_len(n_chains)) {
  boundary_agreement_wide[[paste0("chain_", chain)]] <- mapply(
    function(boundary, disease) {
      boundary_probabilities[boundary, disease, chain]
    },
    boundary_agreement_wide$boundary,
    boundary_agreement_wide$disease
  )
}
chain_probability_columns <- paste0("chain_", seq_len(n_chains))
boundary_agreement_wide$chains_above_0.95 <- rowSums(
  boundary_agreement_wide[, chain_probability_columns, drop = FALSE] >=
    0.95
)
boundary_agreement_wide$chains_above_0.90 <- rowSums(
  boundary_agreement_wide[, chain_probability_columns, drop = FALSE] >=
    0.90
)

diagnostic_cat("\nAdjacency-cardinality summaries by chain\n")
diagnostic_print(cardinality_summary_by_chain, row.names = FALSE)

diagnostic_cat("\nBoundary counts selected by each chain at estimated FDR <= 0.05\n")
diagnostic_print(boundary_count_summary, row.names = FALSE)

diagnostic_cat("\nPooled boundary counts selected at estimated FDR <= 0.05\n")
diagnostic_print(pooled_boundary_count_summary, row.names = FALSE)

diagnostic_cat("\nBetween-chain spread in selected boundary counts\n")
diagnostic_print(boundary_count_ranges, row.names = FALSE)

diagnostic_cat("\nCancer ranking by FDR-selected boundary prevalence\n")
diagnostic_print(boundary_prevalence_ranking, row.names = FALSE)

diagnostic_cat("\nBetween-chain spread in edge-level boundary probabilities\n")
diagnostic_print(boundary_probability_range_summary, row.names = FALSE)

diagnostic_cat("\nPairwise agreement of edge-level boundary probabilities\n")
diagnostic_print(boundary_probability_pairwise_summary, row.names = FALSE)

diagnostic_cat("\nPairwise Jaccard overlap among chain-specific selected boundaries\n")
diagnostic_print(boundary_pairwise_jaccard_summary, row.names = FALSE)

diagnostic_cat("\nFDR-selection frequency of boundaries across chains\n")
diagnostic_print(boundary_selection_frequency_summary, row.names = FALSE)

diagnostic_cat("\nEdges closest to each pooled FDR cutoff\n")
diagnostic_print(boundary_cutoff_sensitivity, row.names = FALSE)

diagnostic_cat("\nLargest boundary probability disagreements with chain-specific probabilities\n")
diagnostic_print(
  head(
    boundary_agreement_wide[
      ,
      c(
        "boundary", "disease", "mean_probability",
        "between_chain_range", chain_probability_columns,
        "chains_above_0.95", "chains_above_0.90"
      )
    ],
    20L
  ),
  row.names = FALSE
)

for (chain in seq_len(n_chains)) {
  chain_fits[[chain]]$posterior_summary <-
    posterior_summary_all_groups[
      posterior_summary_all_groups$chain == paste0("chain_", chain),
      ,
      drop = FALSE
    ]
  chain_fits[[chain]]$boundary_count_summary <-
    boundary_count_summary[boundary_count_summary$chain ==
                             paste0("chain_", chain), ]
  chain_fits[[chain]]$edge_boundary_probabilities <-
    boundary_probabilities[, , chain, drop = FALSE]
  chain_fits[[chain]]$selected_boundary_sets <- lapply(
    selected_boundary_sets,
    `[[`,
    paste0("chain_", chain)
  )
  saveRDS(chain_fits[[chain]], chain_files[[chain]])
}

pooled_boundary_map_file <- file.path(
  figures_dir,
  paste0("pooled_fdr_boundary_maps_", cvrts, ".pdf")
)
pooled_boundary_map_status <- list(
  requested = write_pooled_boundary_map,
  output_file = pooled_boundary_map_file,
  status = "not_requested",
  message = NA_character_
)

if (write_pooled_boundary_map) {
  pooled_boundary_map_status <- tryCatch(
    {
      grDevices::pdf(pooled_boundary_map_file, width = 16, height = 5)
      tryCatch(
        {
          graphics::par(
            mfrow = c(1L, n_diseases),
            mar = c(0.2, 0.2, 2.5, 0.2),
            oma = c(0, 0, 1.5, 0)
          )
          county_geometry <- sf::st_geometry(ca.poly)

          for (disease in seq_len(n_diseases)) {
            selected_names <- pooled_selected_boundary_sets[[
              disease_names[[disease]]
            ]]
            selected_indices <- match(selected_names, neighbor_name)
            selected_indices <- selected_indices[!is.na(selected_indices)]

            plot(
              county_geometry,
              col = "white",
              border = "grey65",
              lwd = 0.5,
              axes = FALSE,
              main = paste0(
                disease_names[[disease]],
                " (", length(selected_indices), ")"
              )
            )

            for (edge in selected_indices) {
              first_county <- neighbor_list0[edge, 1L]
              second_county <- neighbor_list0[edge, 2L]
              shared_boundary <- suppressWarnings(sf::st_intersection(
                sf::st_boundary(county_geometry[first_county]),
                sf::st_boundary(county_geometry[second_county])
              ))
              if (length(shared_boundary) > 0L &&
                  !all(sf::st_is_empty(shared_boundary))) {
                plot(
                  shared_boundary,
                  add = TRUE,
                  col = "#2166AC",
                  lwd = 3
                )
              }
            }
          }
          graphics::mtext(
            "Pooled posterior boundary selection (estimated FDR <= 0.05)",
            outer = TRUE,
            line = 0.2,
            cex = 1.1
          )
        },
        finally = grDevices::dev.off()
      )

      list(
        requested = TRUE,
        output_file = pooled_boundary_map_file,
        status = "written",
        message = NA_character_
      )
    },
    error = function(error) {
      warning(
        "Could not write the pooled boundary map: ",
        conditionMessage(error)
      )
      list(
        requested = TRUE,
        output_file = pooled_boundary_map_file,
        status = "error",
        message = conditionMessage(error)
      )
    }
  )
}

all_results <- list(
  settings = list(
    run_mode = run_mode,
    smoke_test = smoke_test,
    quick_test = quick_test,
    medium_pilot = medium_pilot,
    run_sampling = run_sampling,
    reuse_completed_chains = reuse_completed_chains,
    sampler_md5 = sampler_md5,
    write_pooled_boundary_map = write_pooled_boundary_map,
    n_chains = n_chains,
    chain_seeds = chain_seeds,
    initialization_regimes = initialization_regimes,
    warmup_iterations = warmup_iterations,
    sampling_iterations = sampling_iterations,
    retained_draws = retained_draws,
    sampling_thin = sampling_thin,
    transition_blocks = "non-tempered baseline blocks only",
    state_consistency_tolerance = state_consistency_tolerance,
    cvrts = cvrts,
    diagnostic_backend = diagnostic_backend
  ),
  chain_run_status = chain_run_status,
  chain_fits = chain_fits
)

diagnostics <- list(
  diagnostic_backend = diagnostic_backend,
  state_consistency_summary = state_consistency_summary,
  diagnostic_overview = diagnostic_overview,
  diagnostic_results = diagnostic_results,
  posterior_summary_all_groups = posterior_summary_all_groups,
  posterior_agreement_all_groups = posterior_agreement_all_groups,
  posterior_agreement_overview = posterior_agreement_overview,
  fitted_risk_pairwise_agreement = fitted_risk_pairwise_agreement,
  fitted_risk_pairwise_summary = fitted_risk_pairwise_summary,
  acceptance_summary = acceptance_summary,
  deterministic_state_proposal_summary =
    deterministic_state_proposal_summary,
  cardinality_summary_by_chain = cardinality_summary_by_chain,
  cardinality_posterior_summary = cardinality_posterior_summary,
  draw_boundary_count_posterior_summary =
    draw_boundary_count_posterior_summary
)

boundary_summaries <- list(
  boundary_probabilities = boundary_probabilities,
  boundary_probability_mean = boundary_probability_mean,
  boundary_probability_range = boundary_probability_range,
  boundary_probability_mean_pairwise_absolute_difference =
    boundary_probability_mean_pairwise_absolute_difference,
  boundary_probability_range_summary =
    boundary_probability_range_summary,
  boundary_probability_pairwise_comparison =
    boundary_probability_pairwise_comparison,
  boundary_probability_pairwise_summary =
    boundary_probability_pairwise_summary,
  boundary_probability_correlation_matrices =
    boundary_probability_correlation_matrices,
  boundary_agreement = boundary_agreement,
  boundary_agreement_wide = boundary_agreement_wide,
  boundary_count_summary = boundary_count_summary,
  pooled_boundary_count_summary = pooled_boundary_count_summary,
  boundary_count_ranges = boundary_count_ranges,
  boundary_prevalence_ranking = boundary_prevalence_ranking,
  boundary_prevalence_rank_agreement =
    boundary_prevalence_rank_agreement,
  selected_boundary_sets = selected_boundary_sets,
  pooled_selected_boundary_sets = pooled_selected_boundary_sets,
  boundary_pairwise_jaccard = boundary_pairwise_jaccard,
  boundary_pairwise_jaccard_summary =
    boundary_pairwise_jaccard_summary,
  boundary_jaccard_matrices = boundary_jaccard_matrices,
  boundary_selection_incidence = boundary_selection_incidence,
  boundary_selection_frequency_summary =
    boundary_selection_frequency_summary,
  boundary_cutoff_sensitivity = boundary_cutoff_sensitivity,
  pooled_boundary_map_status = pooled_boundary_map_status
)

# Save the complete binary objects before creating convenience exports, so the
# full diagnostics survive even if a later CSV or report write is interrupted.
saveRDS(
  all_results,
  file.path(output_dir, "multiple_chains_all_results.rds")
)
saveRDS(
  diagnostics,
  file.path(output_dir, "multiple_chains_diagnostics.rds")
)
saveRDS(
  boundary_summaries,
  file.path(output_dir, "multiple_chains_boundary_summaries.rds")
)

# Persist diagnostics in formats suitable for both programmatic analysis and
# direct reading. The RDS files preserve every object exactly; the CSV files
# expose the main tabular results without requiring console copy/paste.
diagnostic_tables <- list(
  chain_run_status = chain_run_status,
  state_consistency = state_consistency_summary,
  diagnostic_overview = diagnostic_overview,
  posterior_summary_all_groups = posterior_summary_all_groups,
  posterior_agreement_all_groups = posterior_agreement_all_groups,
  posterior_agreement_overview = posterior_agreement_overview,
  fitted_risk_pairwise_agreement = fitted_risk_pairwise_agreement,
  fitted_risk_pairwise_summary = fitted_risk_pairwise_summary,
  acceptance_rates = acceptance_summary,
  deterministic_state_proposals = deterministic_state_proposal_summary,
  adjacency_cardinality_by_chain = cardinality_summary_by_chain,
  adjacency_cardinality_posterior = cardinality_posterior_summary,
  boundary_count_posterior = draw_boundary_count_posterior_summary,
  boundary_probability_range_summary = boundary_probability_range_summary,
  boundary_probability_pairwise = boundary_probability_pairwise_comparison,
  boundary_probability_pairwise_summary =
    boundary_probability_pairwise_summary,
  boundary_agreement = boundary_agreement,
  boundary_agreement_by_chain = boundary_agreement_wide,
  boundary_counts_by_chain = boundary_count_summary,
  pooled_boundary_counts = pooled_boundary_count_summary,
  boundary_count_ranges = boundary_count_ranges,
  boundary_prevalence_ranking = boundary_prevalence_ranking,
  boundary_prevalence_rank_agreement =
    boundary_prevalence_rank_agreement,
  boundary_pairwise_jaccard = boundary_pairwise_jaccard,
  boundary_pairwise_jaccard_summary = boundary_pairwise_jaccard_summary,
  boundary_selection_incidence = boundary_selection_incidence,
  boundary_selection_frequency = boundary_selection_frequency_summary,
  boundary_cutoff_sensitivity = boundary_cutoff_sensitivity
)

for (group_name in names(diagnostic_results)) {
  diagnostic_tables[[paste0("rhat_ess_mcse_", group_name)]] <-
    diagnostic_results[[group_name]]$summary
  diagnostic_tables[[paste0("worst_rhat_", group_name)]] <-
    diagnostic_results[[group_name]]$worst
  diagnostic_tables[[paste0("worst_agreement_", group_name)]] <-
    diagnostic_results[[group_name]]$posterior_agreement_worst
}

boundary_probabilities_long <- as.data.frame.table(
  boundary_probabilities,
  responseName = "boundary_probability",
  stringsAsFactors = FALSE
)
diagnostic_tables$edge_boundary_probabilities_by_chain <-
  boundary_probabilities_long

selected_boundaries_to_data_frame <- function(boundary_sets, pooled = FALSE) {
  rows <- list()
  row_index <- 1L
  for (disease in names(boundary_sets)) {
    disease_sets <- boundary_sets[[disease]]
    if (pooled) {
      disease_sets <- list(pooled = disease_sets)
    }
    for (chain in names(disease_sets)) {
      boundaries <- disease_sets[[chain]]
      if (length(boundaries) > 0L) {
        rows[[row_index]] <- data.frame(
          disease = disease,
          chain = chain,
          boundary = as.character(boundaries),
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(
      disease = character(),
      chain = character(),
      boundary = character()
    ))
  }
  do.call(rbind, rows)
}

diagnostic_tables$selected_boundaries_by_chain <-
  selected_boundaries_to_data_frame(selected_boundary_sets)
diagnostic_tables$selected_boundaries_pooled <-
  selected_boundaries_to_data_frame(
    pooled_selected_boundary_sets,
    pooled = TRUE
  )

diagnostic_csv_manifest <- do.call(
  rbind,
  lapply(names(diagnostic_tables), function(table_name) {
    table_value <- diagnostic_tables[[table_name]]
    csv_file <- file.path(
      diagnostics_dir, paste0(table_name, "_", cvrts, ".csv")
    )
    utils::write.csv(table_value, csv_file, row.names = FALSE, na = "")
    data.frame(
      diagnostic = table_name,
      rows = nrow(table_value),
      columns = ncol(table_value),
      file = normalizePath(csv_file, winslash = "/", mustWork = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
utils::write.csv(
  diagnostic_csv_manifest,
  file.path(diagnostics_dir, paste0("diagnostics_manifest_", cvrts, ".csv")),
  row.names = FALSE
)

diagnostic_report_sections <- list(
  "Chain run status" = chain_run_status,
  "State consistency" = state_consistency_summary,
  "Rhat, ESS, and MCSE overview" = diagnostic_overview,
  "Posterior agreement overview" = posterior_agreement_overview,
  "Fitted-risk agreement" = fitted_risk_pairwise_summary,
  "Proposal acceptance" = acceptance_summary,
  "Allocation changes from joint proposals" =
    deterministic_state_proposal_summary,
  "Adjacency cardinality by chain" = cardinality_summary_by_chain,
  "Boundary counts by chain" = boundary_count_summary,
  "Pooled boundary counts" = pooled_boundary_count_summary,
  "Boundary-count ranges" = boundary_count_ranges,
  "Boundary probability disagreement" =
    boundary_probability_range_summary,
  "Boundary probability pairwise agreement" =
    boundary_probability_pairwise_summary,
  "Boundary-set Jaccard agreement" = boundary_pairwise_jaccard_summary,
  "Boundary selection frequency" = boundary_selection_frequency_summary,
  "Largest edge-level disagreements" = head(
    boundary_agreement_wide[
      order(
        boundary_agreement_wide$between_chain_range,
        decreasing = TRUE
      ),
    ],
    20L
  )
)

diagnostic_report_file <- file.path(
  diagnostics_dir,
  paste0("diagnostics_summary_", cvrts, ".txt")
)
report_connection <- file(diagnostic_report_file, open = "wt")
tryCatch(
  {
    writeLines(c(
      "Baseline multiple-chain SEER diagnostics",
      paste("Generated:", format(Sys.time())),
      paste("Diagnostic backend:", diagnostic_backend),
      paste("Run mode:", run_mode),
      paste("Chains:", n_chains),
      paste("Warm-up iterations per chain:", warmup_iterations),
      paste("Sampling iterations per chain:", sampling_iterations),
      paste("Retained draws per chain:", retained_draws),
      ""
    ), report_connection)
    for (section_name in names(diagnostic_report_sections)) {
      writeLines(c(
        paste0("=== ", section_name, " ==="),
        capture.output(print(
          diagnostic_report_sections[[section_name]],
          row.names = FALSE
        )),
        ""
      ), report_connection)
    }
  },
  finally = close(report_connection)
)

# Keep the machine-readable table inventory inside the saved diagnostics too.
# This makes every exported table discoverable without relying on directory
# listings and records its specification-specific file name.
diagnostics$diagnostic_csv_manifest <- diagnostic_csv_manifest

saveRDS(
  diagnostics,
  file.path(diagnostics_dir, "multiple_chains_diagnostics.rds")
)
saveRDS(
  boundary_summaries,
  file.path(
    diagnostics_dir,
    "multiple_chains_boundary_summaries.rds"
  )
)

cat("\nBaseline multiple-chain workflow complete.\n")
cat("Finish time:", format(Sys.time()), "\n")
cat("Saved outputs in:", output_dir, "\n")
cat("Saved diagnostics in:", diagnostics_dir, "\n")
