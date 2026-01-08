rm(list = ls())
setwd("C:/Dati/Dottorato/Visiting UCLA/Spatial Disease Mapping/JASA_submission/Lee_Mitchell_comparison")

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
library(msos)
library(hesim)
library(LaplacesDemon)
library(blockmatrix)

### Functions to generate discrete random effects

# compute stick-breaking weights

makeprobs <- function(v){
  
  m <- length(v)
  p <- matrix(0,m,m)
  
  for(j in 2:m){
    p[1:(j-1),j] <- 1
  }
  
  probs <- exp(log(v)+log(1-v)%*%p)
  probs[m] <- 1-sum(probs[1:(m-1)])
  probs
}

# compute u weights from F_r and prob

makeu <- function(F_r,probs){
  
  m1 <- length(F_r)
  m2 <- length(probs)
  
  u <- rep(0,m1)
  
  for (k in 1:m1){
    for (l in 1:m2){
      if (sum(probs[1:(l-1)])<F_r[k] && F_r[k]<sum(probs[1:l])){
        u[k] <- l
      }
      if (u[k]==0){
        u[k] <- 1
      }
    }
  }
  
  return(u)
}

# Compute presicion matrices of DAGAR

Dinv_new <- function(rho, n, W){
  
  Winc <- W
  Winc[upper.tri(Winc)] <- 0
  
  ns <- rowSums(Winc)
  Tau <- diag((1 + (ns - 1) * rho^2) / (1 - rho^2))
  B <- matrix(0, n, n)
  
  indices <- which(Winc == 1, arr.ind = TRUE)
  
  # Iterate over the indices
  for (i in 1:nrow(indices)) {
    
    row_idx <- indices[i, "row"]
    col_idx <- indices[i, "col"]
    
    B[row_idx,col_idx] <- rho / (1 + (ns[row_idx] - 1) * rho^2)
    
  }
  
  invD <- t(diag(n) - B) %*% Tau %*% (diag(n) - B)
  
  return(invD)
  
}


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

county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])
ca.poly <- map2SpatialPolygons(ca.county, IDs=county.ID)
ca.coords <- coordinates(ca.poly)
n_county <- length(county.ID)

## Adjacency matrix
ca.neighbors <- poly2nb(ca.poly)
n <- length(ca.neighbors)

Adjs <- sapply(ca.neighbors,function(x,n) {v=rep(0,n);v[x]=1;v},n)

colnames(Adjs) <- county.ID

#rownames(Adjs) <- county.ID

num_edge <- sum(Adjs)/2

neighborvec0 <- NULL
neighbor_list0 <- NULL

for(i in 1:(n_county-1)){
  for(j in (i+1):n_county){
    if(Adjs[i,j] == 1){
      neighborvec0 <- c(neighborvec0, paste(i, ",", j, sep=""))
      neighbor_list0 <- rbind(neighbor_list0, c(i,j))
    }
  }
}

## Reorder the map
ca.latrange <- round(quantile(ca.coords[,2],c(0.25,0.75)))
ca.albersproj <- mapproject(ca.coords[,1],ca.coords[,2],projection = "albers",param=ca.latrange)
projmat <- cbind(ca.albersproj$x,ca.albersproj$y)

perm <- order(ca.albersproj$x-ca.albersproj$y)
colnames(Adjs)[perm]

Adj_new <- Adjs[perm,perm]

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

Minc <- Adjs[final_perm,final_perm]

### Specify Covariance matrix for spatial components

set.seed(160495)

X <- cbind(1, rnorm(n_county))

beta <- c(2,8)

taud <- 500

Z <- sd_diff_mat(X[,2],Minc)
# Z <- as.matrix(dist(X[,2], diag=TRUE, upper=TRUE))

M <- as.numeric(-log(0.5)/quantile(Z[which(Z!=0)],0.5))

eta <- 0.4

rho <- 0.2

W_true <- Minc

indices <- which(Minc == 1, arr.ind = TRUE)

for (i in 1:nrow(indices)) {
  row_idx <- indices[i, "row"]
  col_idx <- indices[i, "col"]
  
  W_true[row_idx,col_idx] <- as.numeric(exp(-Z[row_idx,col_idx]*eta) >= 0.5)
}

Q <- Dinv_new(rho, n, W_true)

Vr <- solve(Q)

alpha <- 5
taus <- 0.1
sqrt(1/taus)

K <- 15

W_true <- W_true[order(final_perm),order(final_perm)]
W <- Minc[order(final_perm),order(final_perm)]

Y_list <- list()

phi_true <- list()

# ADAGAR ------------------------------------------------------------------

library(Rcpp)
library(RcppArmadillo)

sourceCpp('sampler_sim_gaussian_DAGAR_fastest_univariate.cpp')

library(tictoc)

seed <- 1

for(seed in 1:100){
  print(seed)
  
  set.seed(seed)
  
  r <- rmvnorm(1, rep(0, n), Vr)
  v <- rbeta(K, 1, alpha)
  
  probvec <- makeprobs(v)
  F_r <- pnorm(r, 0, sqrt(diag(Vr)))
  u <- makeu(F_r, probvec)
  
  thetavec <- rnorm(K, 0, sqrt(1/taus))
  phivec <- thetavec[u]
  
  phi <- phivec[1:n][order(final_perm)]
  
  phi_true[[seed]] <- phi
  
  y <- X[order(final_perm),] %*% beta + (phi - mean(phi))  + sqrt(1/taud) * rnorm(n_county)
  
  yo <- y[final_perm]
  
  Y <- yo
  
  Y_list[[seed]] <- Y
  
  tic()
  
  mcmc_samples <- ADAGAR(y=Y, X=X, Z1=Z,
                          q=1, Minc=Minc,
                          alpha=1, n_atoms=15,
                          runs=10000, burn=10000, thin=1)
  
  toc()
  
  filename <- paste0("runs_DAGAR/mcmc_samples_", seed, ".rds")
  saveRDS(mcmc_samples, file = filename)
  
}

saveRDS(W_true, "RE_generation_DAGAR/W_true.rds")

saveRDS(X, "RE_generation_DAGAR/X.rds")

saveRDS(phi_true, "RE_generation_DAGAR/phi_true.rds")

saveRDS(Z, "RE_generation_DAGAR/Z.rds")

saveRDS(Y_list, "RE_generation_DAGAR/Y_list.rds")



# CARBayes ----------------------------------------------------------------

library(CARBayes)

Z[is.na(Z)] <- 0

for(seed in 1:100){
  print(seed)
  
  set.seed(seed)
  
  r <- rmvnorm(1, rep(0, n), Vr)
  v <- rbeta(K, 1, alpha)
  
  probvec <- makeprobs(v)
  F_r <- pnorm(r, 0, sqrt(diag(Vr)))
  u <- makeu(F_r, probvec)
  
  thetavec <- rnorm(K, 0, sqrt(1/taus))
  phivec <- thetavec[u]
  
  phi <- phivec[1:n][order(final_perm)]
  
  phi_true[[seed]] <- phi
  
  y <- X[order(final_perm),] %*% beta + (phi - mean(phi))  + sqrt(1/taud) * rnorm(n_county)
  
  yo <- y[final_perm]
  
  Y <- yo
  
  Y_list[[seed]] <- Y
  
  tic()
  
  simulated_data <- data.frame(y = Y, X = X[,2])
  
  formula <- y ~ X
  
  mcmc_samples <- S.CARdissimilarity(formula=formula, data=simulated_data,
                                     family="gaussian", W=W[final_perm, final_perm], Z=list(Z=Z),
                                     W.binary=TRUE, burnin=10000, n.sample=20000, thin=1)
  
  toc()
  
  filename <- paste0("runs_CARBayes/mcmc_samples_", seed, ".rds")
  saveRDS(mcmc_samples, file = filename)
  
}

saveRDS(W_true, "RE_generation_CARBayes/W_true.rds")

saveRDS(X, "RE_generation_CARBayes/X.rds")

saveRDS(phi_true, "RE_generation_CARBayes/phi_true.rds")

saveRDS(Z, "RE_generation_CARBayes/Z.rds")

saveRDS(Y_list, "RE_generation_CARBayes/Y_list.rds")

