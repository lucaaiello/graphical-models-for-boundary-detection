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

TauTau_new <- function(Rho, n, q, W){
  Tau <- list()
  
  for(d in 1:q){
    Winc <- W[[d]]
    Winc[upper.tri(Winc)] <- 0
    
    ns <- rowSums(Winc)
    Tau[[d]] <- diag((1 + (ns - 1) * Rho[d]^2) / (1 - Rho[d]^2))
    
  }
  
  return(Tau)
  
}

BB_new <- function(Rho, n, q, W){
  B <- list()
  
  for(d in 1:q){
    Winc <- W[[d]]
    Winc[upper.tri(Winc)] <- 0
    
    ns <- rowSums(Winc)
    B[[d]] <- matrix(0, n, n)
    
    indices <- which(Winc == 1, arr.ind = TRUE)
    
    # Iterate over the indices
    for (i in 1:nrow(indices)) {
      
      row_idx <- indices[i, "row"]
      col_idx <- indices[i, "col"]
      
      B[[d]][row_idx,col_idx] <- Rho[d] / (1 + (ns[row_idx] - 1) * Rho[d]^2)
      
    }
    
  }
  
  return(B)
  
}

# Jacobian matrix for updating A

Jacob_A <- function(A){
  
  prod <- 1
  
  for(i in 1:nrow(A)){
    prod <- prod*A[i,i]^(nrow(A)-i+1)
  }
  
  return(2^nrow(A)*prod)
}

# Functions for numerical issues
log1pexp <- function(x){
  
  q <- length(x)
  out <- vector(length = q)
  
  for (d in 1:q) {
    if(x[d] < -37){
      out[d] <- exp(x[d])
    }else if(x[d] > -37 & x[d] <= 18){
      out[d] <- log1p(exp(x[d]))
    }else if(x[d] > 18 & x[d] <= 33.3){
      out[d] <- x[d] + exp(-x[d])
    }else if(x[d] > 33.3){
      out[d] <- x[d]
    }
  }
  
  return(out)
}

inv_trans_par <-  function(x,lb,ub){
  
  if(x < -745){
    out <- (lb + ub*exp(-745))/(1+exp(-745))
  }else if(x > 16.81){
    out <- (lb + ub*exp(16.81))/(1+exp(16.81))
  }else if(x >= -745 & x <= 16.81){
    out <- (lb + ub*exp(x))/(1+exp(x))
  }
  
  return(out)
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

W_dis <- cbind(c(0,1,0,1),c(1,0,1,0),c(0,1,0,1),c(1,0,1,0))

D_dis <- diag(rowSums(W_dis))

DWD <- diag(1/sqrt(diag(D_dis)))%*%W_dis%*%diag(1/sqrt(diag(D_dis)))

eigen(DWD)$values

rho_dis <- 0.25

L_dis <- D_dis - rho_dis * W_dis

### Specify Covariance matrix for spatial components

set.seed(160495)

X1 <- cbind(1, rnorm(n_county))
X2 <- X1
X3 <- X1
X4 <- X1

beta1 <- c(2,8)
beta2 <- c(3,7)
beta3 <- c(4,6)
beta4 <- c(5,5)

taud1 <- 500
taud2 <- 500
taud3 <- 500
taud4 <- 500

Z1 <- sd_diff_mat(X1[,2],Minc)

M <- as.numeric(-log(0.5)/quantile(Z1[which(Z1!=0)],0.5))

eta <- c(0.4, 0.2, 0.3, 0.5)

q <- 4

rho <- c(0.2, 0.8, 0.4, 0.6)

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
B <- BB_new(rho, n, q, W)
Tau <- TauTau_new(rho, n, q, W)

U1 <- t(diag(n) - B[[1]])*sqrt(Tau[[1]])/sqrt(L_dis[1,1])
U2 <- t(diag(n) - B[[2]])*sqrt(Tau[[2]])/sqrt(L_dis[2,2])
U3 <- t(diag(n) - B[[3]])*sqrt(Tau[[3]])/sqrt(L_dis[3,3])
U4 <- t(diag(n) - B[[4]])*sqrt(Tau[[4]])/sqrt(L_dis[4,4])

Ublock <- as.matrix(bdiag(bdiag(U1,U2),
                          bdiag(U3,U4)))

Vr <- solve(Ublock%*%kronecker(L_dis,diag(n))%*%t(Ublock))

alpha <- 5
taus <- 0.1

K <- 15

nq <- n*q

W_true1 <- W[[1]][order(final_perm),order(final_perm)]
W_true2 <- W[[2]][order(final_perm),order(final_perm)]
W_true3 <- W[[3]][order(final_perm),order(final_perm)]
W_true4 <- W[[4]][order(final_perm),order(final_perm)]

Y_list <- list()

phi_true1 <- list()
phi_true2 <- list()
phi_true3 <- list()
phi_true4 <- list()

library(Rcpp)
library(RcppArmadillo)
sourceCpp('Misspecified models/sampler_sim_gaussian_DAGAR_fastest.cpp')

library(tictoc)

seed <- 1

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
  
  phi1 <- phivec[1:n][order(final_perm)]
  phi2 <- phivec[(n+1):(2*n)][order(final_perm)]
  phi3 <- phivec[(2*n+1):(3*n)][order(final_perm)]
  phi4 <- phivec[(3*n+1):(4*n)][order(final_perm)]
  
  phi_true1[[seed]] <- phi1
  phi_true2[[seed]] <- phi2
  phi_true3[[seed]] <- phi3
  phi_true4[[seed]] <- phi4
  
  y1 <- X1[order(final_perm),] %*% beta1 + (phi1 - mean(phi1))  + sqrt(1/taud1) * rnorm(n_county)
  y2 <- X2[order(final_perm),] %*% beta2 + (phi2 - mean(phi2))  + sqrt(1/taud2) * rnorm(n_county)
  y3 <- X3[order(final_perm),] %*% beta3 + (phi3 - mean(phi3))  + sqrt(1/taud3) * rnorm(n_county)
  y4 <- X4[order(final_perm),] %*% beta4 + (phi4 - mean(phi4))  + sqrt(1/taud4) * rnorm(n_county)
  
  X <- as.matrix(bdiag(bdiag(X1, X2),bdiag(X3, X4)))
  
  yo1 <- y1[final_perm]
  yo2 <- y2[final_perm]
  yo3 <- y3[final_perm]
  yo4 <- y4[final_perm]
  
  Y <- c(yo1, yo2, yo3, yo4)
  
  Y_list[[seed]] <- Y
  
  tic()

  mcmc_samples <- MADAGAR(y=Y, X=X, Z1=Z1,
                          q=4, Minc=Minc,
                          alpha=1, n_atoms=15,
                          runs=10000, burn=10000, thin=1)

  toc()

  filename <- paste0("Misspecified models/Undirected/Unstructured/runs_DAGAR/mcmc_samples_", seed, ".rds")
  saveRDS(mcmc_samples, file = filename)
  
}

saveRDS(W_true1, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/W_true1.rds")
saveRDS(W_true2, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/W_true2.rds")
saveRDS(W_true3, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/W_true3.rds")
saveRDS(W_true4, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/W_true4.rds")

saveRDS(W_dis, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/W_dis.rds")

saveRDS(Vr, file = "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/Vr_true.rds")

saveRDS(X1, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/X1.rds")
saveRDS(X2, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/X2.rds")
saveRDS(X3, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/X3.rds")
saveRDS(X4, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/X4.rds")

saveRDS(phi_true1, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/phi_true1.rds")
saveRDS(phi_true2, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/phi_true2.rds")
saveRDS(phi_true3, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/phi_true3.rds")
saveRDS(phi_true4, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/phi_true4.rds")

saveRDS(Z1, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/Z1.rds")

saveRDS(Y_list, "Misspecified models/Undirected/Unstructured/RE_generation_DAGAR/Y_list.rds")
