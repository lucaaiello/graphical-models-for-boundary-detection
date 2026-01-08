rm(list = ls())

library(maps)
ca.county <- map("county","california", fill=TRUE, plot=FALSE)
library(spdep)
library(maptools)
library(mapproj)
library(stringr)
library(classInt)
library(RColorBrewer)
library(dplyr)
library(caret)
library(Matrix)
library(mvtnorm)
library(LaplacesDemon)
library(dplyr)
library(caret)
library(matrixStats)
library(monomvn)

county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])
ca.poly <- map2SpatialPolygons(ca.county, IDs=county.ID)
ca.coords <- coordinates(ca.poly)

n_county <- length(county.ID)

###################################
### Assign true mean on the map ###
###################################

# discrete random effects for diseases

phi_true <- readRDS("Lee_Mitchell_comparison/RE_generation_DAGAR/phi_true.rds")

# spatial adjacency matrix for diseases

W_true <- readRDS("Lee_Mitchell_comparison/RE_generation_DAGAR/W_true.rds")

# covariates 

X <- readRDS("Lee_Mitchell_comparison/RE_generation_DAGAR/X.rds")

Z <- readRDS("Lee_Mitchell_comparison/RE_generation_DAGAR/Z.rds")

Y <- readRDS("Lee_Mitchell_comparison/RE_generation_DAGAR/Y_list.rds")

## Adjacency matrix
ca.neighbors <- poly2nb(ca.poly)
n <- length(ca.neighbors)

Adj <- sapply(ca.neighbors,function(x,n) {v=rep(0,n);v[x]=1;v},n)
colnames(Adj) <- county.ID

num_edge <- sum(Adj)/2

neighborvec0 <- NULL
neighbor_list0 <- NULL
for(i in 1:(n_county-1)){
  for(j in (i+1):n_county){
    if(Adj[i,j] == 1){
      neighborvec0 <- c(neighborvec0, paste(i, ",", j, sep=""))
      neighbor_list0 <- rbind(neighbor_list0, c(i,j))
    }
  }
}

## Reorder the map
ca.latrange <- round(quantile(ca.coords[,2],c(0.25,0.75)))
ca.albersproj <- mapproject(ca.coords[,1],ca.coords[,2],projection = "albers",param=ca.latrange)

perm <- order(ca.albersproj$x-ca.albersproj$y)
colnames(Adj)[perm]

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
index <- list()
for(i in 1:(n-2)){
  index[[i]] <- region[-(udnei[(cn[i+1] + 1):(cn[i+1] + ns[i+1])])]
}
index1 <- unlist(index)
mns <- max(dni) + 1

T_edge <- c(45, 50 ,55 ,60, 65, 70, 75, 80, 85, 90, 95, 100, 105)
# T_edge <- c(80, 90, 100, 110, 120)

spec_DAGAR <- matrix(0, 100, length(T_edge))
sens_DAGAR <- matrix(0, 100, length(T_edge))

spec_CARBayes <- matrix(0, 100, length(T_edge))
sens_CARBayes <- matrix(0, 100, length(T_edge))

T_edge_ad <- c(40, 50, 60, 70, 80, 90, 100, 110, 120)

specW <- rep(0, 100)
sensW <- rep(0, 100)

Z <- Z[order(final_perm), order(final_perm)]

X <- X[order(final_perm),]

seed <- 1

for(seed in 1:100){
  
  #######true difference boundary########
  
  true_diff <- NULL
  
  for(i in 1:nrow(neighbor_list0)){
    if(phi_true[[seed]][neighbor_list0[i,1]] != phi_true[[seed]][neighbor_list0[i,2]]){
      true_diff <- c(true_diff, 1)
    }else{
      true_diff <- c(true_diff, 0)
    }
  }
  
  print(seed)
  
  # DAGAR
  
  filename_DAGAR <- paste0("Lee_Mitchell_comparison/runs_DAGAR/mcmc_samples_", seed, ".rds")
  mcmc_samples_DAGAR <- readRDS(filename_DAGAR)
  
  names(mcmc_samples_DAGAR) <- c("beta", "taud", "phi", "theta", "u", "rho", "v", "r", 
                                 "F_r", "eta", "taus", "W")
  
  phis_DAGAR <-  mcmc_samples_DAGAR$phi
  phis_origin_DAGAR <- phis_DAGAR[,order(final_perm)]
  
  # CARBayes
  filename_CARBayes <- paste0("Lee_Mitchell_comparison/runs_CARBayes/mcmc_samples_", seed, ".rds")
  mcmc_samples_CARBayes <- readRDS(filename_CARBayes)
  
  W.prob.CARBayes <- mcmc_samples_CARBayes$localised.structure$W.border.prob[order(final_perm), order(final_perm)]
  
  # difference boundaries for each cancer individually ---------------------------
  
  bd <- function(x){
    diff_b <- NULL
    for(i in 1:nrow(neighbor_list0)){
      if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]]){
        diff_b <- c(diff_b, 1)
      }else{
        diff_b <- c(diff_b, 0)
      }
    }
    return(diff_b)
  }
  
  vij_samples <- t(apply(phis_DAGAR, 1, bd))
  
  # probabilities
  
  pvij_DAGAR <- apply(vij_samples, 2, mean)
  
  pvij_CARBayes <- vector(length = nrow(neighbor_list0))
  
  for(i in 1:nrow(neighbor_list0)){
    pvij_CARBayes[i] <- W.prob.CARBayes[neighbor_list0[i,1],neighbor_list0[i,2]]
  }
  
  # Estimated FDR curves
  
  threshold_DAGAR <- sort(pvij_DAGAR, decreasing = TRUE)[T_edge]
  threshold_CARBayes <- sort(pvij_CARBayes, decreasing = TRUE)[T_edge]
  
  # calculate sensitivity and specificity on boundary detection
  for(i in 1:length(T_edge)){
    est_diff_DAGAR <- factor(as.numeric(pvij_DAGAR >= threshold_DAGAR[i]), levels = c(0,1))
    conf_matrix_DAGAR <- table(est_diff_DAGAR,true_diff)
    
    est_diff_CARBayes <- factor(as.numeric(pvij_CARBayes >= threshold_CARBayes[i]), levels = c(0,1))
    conf_matrix_CARBayes <- table(est_diff_CARBayes,true_diff)
    
    spec_DAGAR[seed,i] <- sensitivity(conf_matrix_DAGAR)
    sens_DAGAR[seed,i] <- specificity(conf_matrix_DAGAR)
    
    spec_CARBayes[seed,i] <- sensitivity(conf_matrix_CARBayes)
    sens_CARBayes[seed,i] <- specificity(conf_matrix_CARBayes)
  }
  
}

spec_mean_DAGAR <- colMeans(spec_DAGAR)
sens_mean_DAGAR <- colMeans(sens_DAGAR)

spec_mean_CARBayes <- colMeans(spec_CARBayes)
sens_mean_CARBayes <- colMeans(sens_CARBayes)

table <- cbind(spec_mean_DAGAR, sens_mean_DAGAR, spec_mean_CARBayes, sens_mean_CARBayes)

colnames(table) <- c("spec", "sens", "spec", "sens")
rownames(table) <- as.character(c(45, 50 ,55 ,60, 65, 70, 75, 80, 85, 90, 95, 100, 105))

round(table,3)

################################################################################
seed <- 1 

specW_DAGAR <- rep(0, 100)
sensW_DAGAR <- rep(0, 100)
specW_CARBayes <- rep(0, 100)
sensW_CARBayes <- rep(0, 100)

for (seed in 1:100) {
  
  #######true adjacencies#######
  
  true_adj <- NULL
  
  for(i in 1:nrow(neighbor_list0)){
    if(W_true[neighbor_list0[i,1],neighbor_list0[i,2]] == 1){
      true_adj <- c(true_adj, 1)
    }else{
      true_adj <- c(true_adj, 0)
    }
  }
  
  print(seed)
  
  filename_DAGAR <- paste0("Lee_Mitchell_comparison/runs_DAGAR/mcmc_samples_", seed, ".rds")
  mcmc_samples_DAGAR <- readRDS(filename_DAGAR)
  
  names(mcmc_samples_DAGAR) <- c("beta", "taud", "phi", "theta", "u", "rho", "v", "r", 
                           "F_r", "eta", "taus", "W")
  
  W.post <- mcmc_samples_DAGAR$W
  
  W.mean <- apply(simplify2array(W.post), c(1, 2), mean)
  
  W.mean <- W.mean[order(final_perm), order(final_perm)]
  
  # calculate sensitivity and specificity on adjacency
  
  W.post.DAGAR <- Minc
  
  eta.mean <- apply(mcmc_samples_DAGAR$eta, 2, mean)
  
  Winc <- Minc[order(final_perm),order(final_perm)]
  Winc[upper.tri(Winc)] <- 0
  indices <- which(Winc == 1, arr.ind = TRUE)
  
  est_adj_DAGAR <- vector(length= nrow(indices))
  est_adj_CARBayes<- vector(length= nrow(indices))
  
  # CARBayes
  filename_CARBayes <- paste0("Lee_Mitchell_comparison/runs_CARBayes/mcmc_samples_", seed, ".rds")
  mcmc_samples_CARBayes <- readRDS(filename_CARBayes)
  
  W.post.CARBayes <- mcmc_samples_CARBayes$localised.structure$W.posterior[order(final_perm), order(final_perm)]
  
  for (i in 1:nrow(indices)) {
    row_idx <- indices[i, "row"]
    col_idx <- indices[i, "col"]
      
    W.post.DAGAR[row_idx,col_idx] <- as.numeric(exp(-Z[row_idx,col_idx]*eta.mean) >= 0.5)
      
    est_adj_DAGAR[i] <- W.post.DAGAR[row_idx,col_idx]
    
    est_adj_CARBayes[i] <- W.post.CARBayes[row_idx,col_idx]
  }
  
  est_adj_DAGAR <- factor(est_adj_DAGAR)
  conf_matrixW_DAGAR <- table(est_adj_DAGAR,true_adj)
  
  if (nrow(conf_matrixW_DAGAR) == 1) {
    conf_matrixW_DAGAR <- rbind(c(0,0),conf_matrixW_DAGAR)
    rownames(conf_matrixW_DAGAR)[1] <- "0"
    conf_matrixW_DAGAR <- as.table(conf_matrixW_DAGAR)
  }
  
  specW_DAGAR[seed] <- sensitivity(conf_matrixW_DAGAR)
  sensW_DAGAR[seed] <- specificity(conf_matrixW_DAGAR)
  
  est_adj_CARBayes <- factor(est_adj_CARBayes)
  conf_matrixW_CARBayes <- table(est_adj_CARBayes,true_adj)
  
  specW_CARBayes[seed] <- sensitivity(conf_matrixW_CARBayes)
  sensW_CARBayes[seed] <- specificity(conf_matrixW_CARBayes)
  
}

specW1_mean_DAGAR <- mean(specW_DAGAR)
sensW1_mean_DAGAR <- mean(sensW_DAGAR)

specW2_mean_CARBayes <- mean(specW_CARBayes)
sensW2_mean_CARBayes <- mean(sensW_CARBayes)

table_adj <- cbind(specW1_mean_DAGAR, sensW1_mean_DAGAR, 
                   specW2_mean_CARBayes, sensW2_mean_CARBayes)

colnames(table_adj) <- c("specW1", "sensW1", "specW2", "sensW2")
round(table_adj,3)

