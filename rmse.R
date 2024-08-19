rm(list = ls())
setwd("C:/Dati/Dottorato/Visiting UCLA/Spatial Disease Mapping/Review/Rcpp/Simulation/DAGAR")

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

Dinv_new <- function(Rho, n, q, W){
  Tau <- list()
  B <- list()
  invD <- list()
  
  for(d in 1:q){
    Winc <- W[[d]]
    Winc[upper.tri(Winc)] <- 0
    
    ns <- rowSums(Winc)
    Tau[[d]] <- diag((1 + (ns - 1) * Rho[d]^2) / (1 - Rho[d]^2))
    B[[d]] <- matrix(0, n, n)
    
    indices <- which(Winc == 1, arr.ind = TRUE)
    
    # Iterate over the indices
    for (i in 1:nrow(indices)) {
      
      row_idx <- indices[i, "row"]
      col_idx <- indices[i, "col"]
      
      B[[d]][row_idx,col_idx] <- Rho[d] / (1 + (ns[row_idx] - 1) * Rho[d]^2)
      
    }
    
    invD[[d]] <- t(diag(n) - B[[d]]) %*% Tau[[d]] %*% (diag(n) - B[[d]])
    
  }
  
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

X1 <- cbind(1, rnorm(n_county))
X2 <- X1
X3 <- X1
X4 <- X1

beta1 <- c(2,1)
beta2 <- c(1,2)
beta3 <- c(2,2)
beta4 <- c(1,2)

beta <- c(beta1,beta2,beta3,beta4)

taud1 <- 10
taud2 <- 10
taud3 <- 10
taud4 <- 10

taud <- c(taud1,taud2,taud3,taud4)

Z1 <- sd_diff_mat(X1[,2],Minc)

M <- as.numeric(-log(0.5)/quantile(Z1[which(Z1!=0)],0.5))

eta <- c(0.4, 0.2, 0.3, 0.5)

q <- 4

rho <- c(0.2, 0.8, 0.4, 0.6)

A <- as.matrix(cbind(c(1, 1, 1, 1), c(0, 1, 1, 1), c(0, 0, 1, 1), c(0, 0, 0, 1)))

kprod <- kronecker(A, diag(n))

W <- replicate(q, Minc, simplify = FALSE)

indices <- which(Minc == 1, arr.ind = TRUE)

for (d in 1:q) {
  
  for (i in 1:nrow(indices)) {
    row_idx <- indices[i, "row"]
    col_idx <- indices[i, "col"]
    
    W[[d]][row_idx,col_idx] <- as.numeric(exp(-Z1[row_idx,col_idx]*eta[d]) >= 0.5)
  }
  
}

Q <- Dinv_new(rho, n, q, W)
invQ1 <- solve(Q[[1]])
invQ2 <- solve(Q[[2]])
invQ3 <- solve(Q[[3]])
invQ4 <- solve(Q[[4]])

invQ <- as.matrix(bdiag(bdiag(invQ1, invQ2), bdiag(invQ3, invQ4)))

Vr <- as.matrix(kprod %*% invQ %*% t(kprod))

alpha <- 1
taus <- 1/4

K <- 15

nq <- n*q

mse_beta <- vector(length = 100) 
mse_phi1 <- vector(length = 100)
mse_phi2 <- vector(length = 100)
mse_phi3 <- vector(length = 100)
mse_phi4 <- vector(length = 100)
mse_r1 <- vector(length = 100)
mse_r2 <- vector(length = 100)
mse_r3 <- vector(length = 100)
mse_r4 <- vector(length = 100)
mse_theta <- vector(length = 100)
mse_taud <- vector(length = 100)
mse_v <- vector(length = 100)
mse_tau <- vector(length = 100)
mse_rho <- vector(length = 100)
mse_eta <- vector(length = 100)
mse_A <- vector(length = 100)

for(seed in 1:100){
  print(seed)
  set.seed(seed)
  
  r <- rmvnorm(1, rep(0, nq), Vr)
  v <- rbeta(K, 1, alpha)
  
  probvec <- makeprobs(v)
  F_r <- pnorm(r, 0, sqrt(diag(Vr)))
  u <- makeu(F_r, probvec)
  
  thetavec <- rnorm(K, 0, sqrt(1/taus))
  phivec <- thetavec[u]
  
  filename <- paste0("runs_DAGAR/mcmc_samples_", seed, ".rds")
  
  mcmc_samples <- readRDS(filename)
  names(mcmc_samples) <- c("beta", "taud", "phi", "theta", "u", "rho", "v", "r", 
                           "F_r", "eta", "taus", "W1", "W2", "W3", "W4", "A")
 
  post_mean_beta <- apply(mcmc_samples$beta, 2, mean)
  post_mean_phi1 <- apply(mcmc_samples$phi[,1:n], 2, mean)
  post_mean_phi2 <- apply(mcmc_samples$phi[,(n+1):(2*n)], 2, mean)
  post_mean_phi3 <- apply(mcmc_samples$phi[,(2*n+1):(3*n)], 2, mean)
  post_mean_phi4 <- apply(mcmc_samples$phi[,(3*n+1):(4*n)], 2, mean)
  post_mean_r1 <- apply(mcmc_samples$r[,1:n], 2, mean)
  post_mean_r2 <- apply(mcmc_samples$r[,(n+1):(2*n)], 2, mean)
  post_mean_r3 <- apply(mcmc_samples$r[,(2*n+1):(3*n)], 2, mean)
  post_mean_r4 <- apply(mcmc_samples$r[,(3*n+1):(4*n)], 2, mean)
  post_mean_theta <- apply(mcmc_samples$theta, 2, mean)
  post_mean_taud <- apply(mcmc_samples$taud, 2, mean)
  post_mean_v <- apply(mcmc_samples$v, 2, mean)
  post_mean_taus <- mean(mcmc_samples$taus)
  post_mean_rho <- apply(mcmc_samples$rho, 2, mean)
  post_mean_eta <- apply(mcmc_samples$eta, 2, mean)
  post_mean_A <- apply(mcmc_samples$A, 2, mean)
  
  mse_beta[seed] <- sum((post_mean_beta - beta)^2)
  mse_phi1[seed] <- sum((post_mean_phi1 - phivec[1:n])^2)
  mse_phi2[seed] <- sum((post_mean_phi2 - phivec[(n+1):(2*n)])^2)
  mse_phi3[seed] <- sum((post_mean_phi3 - phivec[(2*n+1):(3*n)])^2)
  mse_phi4[seed] <- sum((post_mean_phi4 - phivec[(3*n+1):(4*n)])^2)
  mse_r1[seed] <- sum((post_mean_r1 - r[1:n])^2)
  mse_r2[seed] <- sum((post_mean_r2 - r[(n+1):(2*n)])^2)
  mse_r3[seed] <- sum((post_mean_r3 - r[(2*n+1):(3*n)])^2)
  mse_r4[seed] <- sum((post_mean_r4 - r[(3*n+1):(4*n)])^2)
  mse_theta[seed] <- sum((post_mean_theta - thetavec)^2)
  mse_taud[seed] <- sum((post_mean_taud - taud)^2)
  mse_v[seed] <- sum((post_mean_v - v)^2)
  mse_tau[seed] <- (post_mean_taus - taus)^2
  mse_rho[seed] <- sum((post_mean_rho - rho)^2)
  mse_eta[seed] <- sum((post_mean_eta - eta)^2)
  mse_A[seed] <- sum((post_mean_A - as.vector(A)[-c(5,9,10,13,14,15)])^2)
  
  }

rmse_beta <- sqrt(median(mse_beta))
rmse_phi1 <- sqrt(median(mse_phi1))
rmse_phi2 <- sqrt(median(mse_phi2))
rmse_phi3 <- sqrt(median(mse_phi3))
rmse_phi4 <- sqrt(median(mse_phi4))
rmse_r1 <- sqrt(median(mse_r1))
rmse_r2 <- sqrt(median(mse_r2))
rmse_r3 <- sqrt(median(mse_r3))
rmse_r4 <- sqrt(median(mse_r4))
rmse_theta <- sqrt(median(mse_theta))
rmse_taud <- sqrt(median(mse_taud))
rmse_v <- sqrt(median(mse_v))
rmse_tau <- sqrt(median(mse_tau))
rmse_rho <- sqrt(median(mse_rho))
rmse_eta <- sqrt(median(mse_eta))
rmse_A <- sqrt(median(mse_A))

rmse <- c(rmse_beta,
          rmse_phi1,rmse_phi2,rmse_phi3,rmse_phi4,
          rmse_r1,rmse_r2,rmse_r3,rmse_r4,
          rmse_theta,
          rmse_taud,
          rmse_v,
          rmse_tau,
          rmse_rho,
          rmse_eta,
          rmse_A)

saveRDS(rmse, file = "Misspecified models/rmse/rmse_unstructured_unstructured.rds")
