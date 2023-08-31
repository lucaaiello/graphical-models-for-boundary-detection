rm(list = ls())
setwd("C:/Dati/Dottorato/Visiting UCLA/Spatial Disease Mapping/Bibliografia/gao_2022/simulation/Disease graph/unstructured")

library(maps)
ca.county <- map("county","california", fill=TRUE, plot=FALSE)
library(spdep)
library(maptools)
library(mapproj)
library(stringr)
library(classInt)
library(RColorBrewer)
library(rjags)
library(R2jags)
library(dplyr)
library(caret)
library(Matrix)
library(mvtnorm)
library(LaplacesDemon)

county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])
ca.poly <- map2SpatialPolygons(ca.county, IDs=county.ID)
ca.coords <- coordinates(ca.poly)

n_county <- length(county.ID)

###################################
### Assign true mean on the map ###
###################################

# discrete random effects for diseases

phi_true1 <- readRDS("phi_true1_try.rds")
phi_true2 <- readRDS("phi_true2_try.rds")
phi_true3 <- readRDS("phi_true3_try.rds")
phi_true4 <- readRDS("phi_true4_try.rds")

# spatial adjacency matrix for diseases

W_true1 <- readRDS("W_true1_try.rds")
W_true2 <- readRDS("W_true2_try.rds")
W_true3 <- readRDS("W_true3_try.rds")
W_true4 <- readRDS("W_true4_try.rds")

# covariates 

Z <- readRDS("Z_try.rds")

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


#######true difference boundary########

true_diff1 <- NULL
true_diff2 <- NULL
true_diff3 <- NULL
true_diff4 <- NULL
for(i in 1:nrow(neighbor_list0)){
  if(phi_true1[neighbor_list0[i,1]] != phi_true1[neighbor_list0[i,2]]){
    true_diff1 <- c(true_diff1, 1)
  }else{
    true_diff1 <- c(true_diff1, 0)
  }
  if(phi_true2[neighbor_list0[i,1]] != phi_true2[neighbor_list0[i,2]]){
    true_diff2 <- c(true_diff2, 1)
  }else{
    true_diff2 <- c(true_diff2, 0)
  }
  if(phi_true1[neighbor_list0[i,1]] != phi_true2[neighbor_list0[i,2]]){
    true_diff3 <- c(true_diff3, 1)
  }else{
    true_diff3 <- c(true_diff3, 0)
  }
  if(phi_true2[neighbor_list0[i,1]] != phi_true1[neighbor_list0[i,2]]){
    true_diff4 <- c(true_diff4, 1)
  }else{
    true_diff4 <- c(true_diff4, 0)
  }
}

sum(true_diff1)
sum(true_diff2)
sum(true_diff3)
sum(true_diff4)

true_diffm <- rbind(true_diff1, true_diff2, true_diff3, true_diff4)

#######true adjacencies#######

true_adj1 <- NULL
true_adj2 <- NULL

for(i in 1:nrow(neighbor_list0)){
  if(W_true1[neighbor_list0[i,1],neighbor_list0[i,2]] == 1){
    true_adj1 <- c(true_adj1, 1)
  }else{
    true_adj1 <- c(true_adj1, 0)
  }
  if(W_true2[neighbor_list0[i,1],neighbor_list0[i,2]] == 1){
    true_adj2 <- c(true_adj2, 1)
  }else{
    true_adj2 <- c(true_adj2, 0)
  }
}


# Functions for MCMC updates ---------------------------------------------------

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

### MDAGAR with adjacency model with smoking covariate -------------------------

# Metropolis within Gibbs Sampler for MCMC updating 
ARDP_joint_diff <- function(y = Y, X = X, Z = Z, Minc = Minc, 
                            alpha = 1, q = 4, n.atoms = 15, 
                            runs = 1000, burn = 500){
  
  #y:       data
  #x:       covariates
  #n.atoms: number of atoms in the mixture dist.
  #theta:   the theta's (iid from baseline) in the model
  #alpha:   v~beta(1,alpha)
  #u:       the index indicator of spatial random effect
  #rho:     spatial autocorrelation parameter in DAGAR
  #Minc:    0-1 adjacency matrix
  #Winc:    0-1 lower triangular adjacency matrix
  #V_r:     covariance matrix of joint MDAGAR
  #Q:       precision matrix of DAGAR
  #r:       random effects following DAGAR
  #F_r:     Marginal CDF of r
  #taus:    precision for theta
  
  nq <- length(y)
  n <- nq/q
  
  p <- 4
  
  y1 <- y[1:n]
  y2 <- y[(n+1):(2*n)]
  y3 <- y[(2*n+1):(3*n)]
  y4 <- y[(3*n+1):(4*n)]
  
  sigmasq_beta <- 1
  keepbeta <- matrix(0,runs,p)
  keepphi <- matrix(0,runs,nq)
  keeptheta <- matrix(0,runs,n.atoms)
  keepu <- matrix(0,runs,nq)
  keeprho <- matrix(0,runs,q)
  keeptaus <- rep(0,runs)
  keepv <- matrix(0,runs,n.atoms)
  keepr <- matrix(0,runs,nq)
  keepFr <- matrix(0,runs,nq)
  keepeta <- matrix(0,runs,q)
  keepW1 <- array(0,dim = c(n,n,runs))
  keepW2 <- array(0,dim = c(n,n,runs))
  keepW3 <- array(0,dim = c(n,n,runs))
  keepW4 <- array(0,dim = c(n,n,runs))
  keepA <- matrix(0,runs,10)
  
  # initial values
  
  theta <- rep(0,n.atoms)
  
  beta <- rep(1,p) 
  
  taus <- 1
  
  c <- 2
  d <- 0.1
  d2 <- 1
  
  v <- rep(.1,n.atoms)
  vv <- log(v) - log(1-v)
  
  probs <- makeprobs(v)
  
  rho <- rep(0.4,q)
  rhorho <- log(rho) - log(1 - rho)
  
  M <- as.numeric(-log(0.5)/quantile(Z[which(Z!=0)],0.5))
  
  eta <- rep(0.25,q)
  etaeta <- log(eta) - log(M - eta)
  
  W <- replicate(q, Minc, simplify = FALSE)
  
  Q <- Dinv_new(rho, n, q, W)
  invQ1 <- solve(Q[[1]])
  invQ2 <- solve(Q[[2]])
  invQ3 <- solve(Q[[3]])
  invQ4 <- solve(Q[[4]])
  
  A <- matrix(0, q, q)
  for(i in 1:q){
    for(j in 1:i){
      if(j == i){
        A[i, j] <- exp(rnorm(1))
      }else{
        A[i, j] <- rnorm(1)
      }
    }
  }
  
  kprod <- kronecker(A, diag(n))
  invQ <- as.matrix(bdiag(bdiag(invQ1, invQ2),
                          bdiag(invQ3, invQ4)))
  
  Vr <- as.matrix(forceSymmetric(kprod %*% invQ %*% t(kprod)))
  
  r <- rmvnorm(1, rep(0, nq), Vr)
  
  F_r <- pnorm(r,0,sqrt(diag(Vr)))
  
  u <- makeu(F_r,probs)
  phi <- theta[u]
  
  nu <- 2
  R <- 0.1 * diag(q)
  
  s_theta <- 0.1
  s_beta <- 0.1
  s_r <- 0.1
  s_vv <- 0.1
  s_rhorho <- 0.1
  s_etaeta <- 0.1
  
  m_theta <- theta
  m_beta <- beta
  m_r <- r
  m_vv <- vv
  m_rhorho <- rhorho
  m_etaeta <- etaeta
  
  R_theta <- diag(0.1,n.atoms)
  R_beta <- diag(0.1,p)
  R_r <- diag(1,nq)
  R_vv <- diag(0.1,n.atoms)
  R_rhorho <- diag(0.1,q)
  R_etaeta <- diag(0.1,q)
  
  xi_theta <- exp(s_theta) * R_theta
  xi_beta <- exp(s_beta) * R_beta
  xi_r <- exp(s_r) * R_r
  xi_vv <- exp(s_vv) * R_vv
  xi_rhorho <- exp(s_rhorho) * R_rhorho
  xi_etaeta <- exp(s_etaeta) * R_etaeta
  
  accepttheta <- 0
  acceptbeta <- 0
  acceptr <- 0
  acceptv <- 0
  acceptrho <- 0
  accepteta <- 0
  acceptA <- 0
  
  count <- 0
  afterburn <- 0
  burn <- burn + 1
  
  for(iter in 1:runs){
    
    if(iter %% 100 == 0){
      print(iter)
      print(accepttheta/(iter-1))
      print(acceptbeta/(iter-1))
      print(acceptr/nq/(iter-1))
      print(acceptv/(iter-1))
      print(acceptrho/(iter-1))
      print(accepteta/(iter-1))
      print(acceptA/(iter-1))
    }
    
    ###################
    ### update beta ###
    ###################
    
    pro_beta <- as.numeric(rmvnorm(1,beta,xi_beta)) 
    
    MH_beta <- sum(- exp(X%*%pro_beta + theta[u]) + 
                     y *(X%*%pro_beta + theta[u])) - 
      sum(- exp(X%*%beta + theta[u]) + 
            y * (X%*%beta + theta[u])) 
    
    MH_beta <- min(0,MH_beta)
    
    if(runif(1,0,1) < exp(MH_beta)){
      beta <- pro_beta
      acceptbeta <- acceptbeta + 1
    }
    
    s_beta <- s_beta + (iter+1)^(-0.7) * (exp(MH_beta) - 0.234)
    
    dbeta_2 <- as.vector(beta - m_beta)
    
    # Covariance
    R_beta <- (1-(iter+1)^(-0.7))*R_beta + 
      (iter+1)^(-0.7)*(dbeta_2%*%t(dbeta_2))
    
    # Mean update
    m_beta <- m_beta + (iter+1)^(-0.7)*dbeta_2
    
    xi_beta <- exp(s_beta) * R_beta
    
    diag(xi_beta) <- diag(xi_beta) + 1e-6 * max(diag(xi_beta))
    
    ####################
    ### update theta ###
    ####################
    
    pro_theta <- rmvnorm(1,theta,xi_theta) 
    
    MH_theta <- sum(- exp(X%*%beta + pro_theta[u]) + y * pro_theta[u]) - 
      sum(taus/2*pro_theta^2) - 
      sum(- exp(X%*%beta + theta[u]) + y * theta[u]) + 
      sum(taus/2*theta^2)
    
    MH_theta <- min(0,MH_theta)
    
    if(runif(1,0,1) < exp(MH_theta)){
      theta <- pro_theta
      accepttheta <- accepttheta + 1
    }
    
    s_theta <- s_theta + (iter+1)^(-0.7) * (exp(MH_theta) - 0.234)
    
    dtheta_2 <- as.vector(theta - m_theta)
    
    # Covariance
    R_theta <- (1-(iter+1)^(-0.7))*R_theta + 
      (iter+1)^(-0.7)*(dtheta_2%*%t(dtheta_2))
    
    # Mean update
    m_theta <- m_theta + (iter+1)^(-0.7)*dtheta_2
    
    xi_theta <- exp(s_theta) * R_theta
    
    diag(xi_theta) <- diag(xi_theta) + 1e-6 * max(diag(xi_theta))
    
    ################
    ### update r ###
    ################
    
    for (k in 1:nq){
      
      pro_r <- r
      pro_Fr <- F_r
      pro_u <- u
      
      pro_r[k] <- rnorm(1,r[k],0.6)
      pro_Fr[k] <- pnorm(pro_r[k],0,sqrt(Vr[k,k]))
      pro_u[k] <- makeu(pro_Fr[k],probs)
      
      MH_r <- dmvnorm(pro_r, mean = rep(0, nq), sigma = Vr, log = T) +
        dpois(y[k], exp(X[k,]%*%beta + theta[pro_u[k]]), log = T) -
        dmvnorm(r, mean = rep(0, nq), sigma = Vr, log = T) -
        dpois(y[k], exp(X[k,]%*%beta + theta[u[k]]), log = T)
      
      MH_r <- min(0,MH_r)
      
      if(runif(1,0,1)<exp(MH_r)){
        r[k] <- pro_r[k]
        F_r[k] <- pro_Fr[k]
        u[k] <- pro_u[k]
        acceptr <- acceptr+1
      }
      
    }
    
    ################
    ### update v ###
    ################
    
    pro_vv <- rmvnorm(1,vv,xi_vv)
    
    pro_v <- apply(pro_vv, 2, inv_trans_par, 0, 1)
    
    pro_probs <- makeprobs(pro_v)
    
    pro_u <- makeu(F_r,pro_probs)
    
    MH_vv <- sum(dbeta(pro_v, 1, alpha, log = T)) +
      sum(dpois(y, exp(X%*%beta + theta[pro_u]), log = T)) +
      sum(log(pro_v) + log(1-pro_v)) -
      sum(dbeta(v, 1, alpha, log = T)) -
      sum(dpois(y, exp(X%*%beta + theta[u]), log = T)) -
      sum(log(v) + log(1-v))
    
    MH_vv <- min(0,MH_vv)
    
    if(runif(1,0,1) < exp(MH_vv)){
      vv <- pro_vv
      v <- pro_v
      probs <- pro_probs
      u <- pro_u
      acceptv <- acceptv + 1
    }
    
    s_vv <- s_vv + (iter+1)^(-0.7) * (exp(MH_vv) - 0.234)
    
    dvv_2 <- as.vector(vv - m_vv)
    
    # Covariance
    R_vv <- (1-(iter+1)^(-0.7)) * R_vv + (iter+1)^(-0.7) * (dvv_2%*%t(dvv_2))
    
    # Mean update
    m_vv <- m_vv + (iter+1)^(-0.7) * dvv_2
    
    xi_vv <- exp(s_vv) * R_vv
    
    diag(xi_vv) <- diag(xi_vv) + 1e-6 * max(diag(xi_vv))
    
    ###################
    ### update taus ###
    ###################
    
    taus <- rgamma(1, shape = n.atoms/2 + c, rate = sum(theta^2)/2 + d2)
    
    ##################
    ### update rho ###
    ##################
    
    pro_rhorho <- rmvnorm(1, rhorho, xi_rhorho)
    
    pro_rho <- apply(pro_rhorho, 2, inv_trans_par, 0, 1)
    
    pro_Q <- Dinv_new(pro_rho, n, q, W)
    pro_invQ1 <- solve(pro_Q[[1]])
    pro_invQ2 <- solve(pro_Q[[2]])
    pro_invQ3 <- solve(pro_Q[[3]])
    pro_invQ4 <- solve(pro_Q[[4]])
    
    kprod <- kronecker(A, diag(n))
    pro_invQ <- as.matrix(bdiag(bdiag(pro_invQ1, pro_invQ2),
                                bdiag(pro_invQ3, pro_invQ4)))
    
    pro_Vr <- as.matrix(forceSymmetric(kprod %*% pro_invQ %*% t(kprod)))
    
    MH_rhorho <- dmvnorm(r, mean = rep(0, nq), sigma = pro_Vr, log = T) +
      sum(log(pro_rho) + log(1-pro_rho)) - 
      dmvnorm(r, mean = rep(0, nq), sigma = Vr, log = T) - 
      sum(log(rho) + log(1-rho)) 
    
    MH_rhorho <- min(0, MH_rhorho)
    
    if(runif(1,0,1) < exp(MH_rhorho)){
      rhorho <- pro_rhorho
      rho <- pro_rho
      Vr <- pro_Vr
      acceptrho <- acceptrho + 1
    }
    
    s_rhorho <- s_rhorho + (iter+1)^(-0.7) * (exp(MH_rhorho) - 0.234)
    
    drhorho_2 <- as.vector(rhorho - m_rhorho)
    
    # Covariance
    R_rhorho <- (1-(iter+1)^(-0.7)) * R_rhorho + (iter+1)^(-0.7) * (drhorho_2%*%t(drhorho_2))
    
    # Mean update
    m_rhorho <- m_rhorho + (iter+1)^(-0.7) * drhorho_2
    
    xi_rhorho <- exp(s_rhorho) * R_rhorho
    
    diag(xi_rhorho) <- diag(xi_rhorho) + 1e-6 * max(diag(xi_rhorho))
    
    ##################
    ### update eta ###
    ##################
    
    pro_etaeta <- rmvnorm(1, etaeta, xi_etaeta)
    
    pro_eta <- apply(pro_etaeta, 2, inv_trans_par, 0, M)
    
    pro_W <- replicate(q, matrix(0, nrow = n, ncol = n), simplify = FALSE)
    
    indices <- which(Minc == 1, arr.ind = TRUE)
    
    for (d in 1:q) {
      
      for (i in 1:nrow(indices)) {
        row_idx <- indices[i, "row"]
        col_idx <- indices[i, "col"]
        
        pro_W[[d]][row_idx,col_idx] <- as.numeric(exp(-Z[row_idx,col_idx]*pro_eta[d]) >= 0.5)
      }
      
    }
    
    pro_Q <- Dinv_new(rho, n, q, pro_W)
    pro_invQ1 <- solve(pro_Q[[1]])
    pro_invQ2 <- solve(pro_Q[[2]])
    pro_invQ3 <- solve(pro_Q[[3]])
    pro_invQ4 <- solve(pro_Q[[4]])
    
    kprod <- kronecker(A, diag(n))
    pro_invQ <- as.matrix(bdiag(bdiag(pro_invQ1, pro_invQ2),
                                bdiag(pro_invQ3, pro_invQ4)))
    
    pro_Vr <- as.matrix(forceSymmetric(kprod %*% pro_invQ %*% t(kprod)))
    
    MH_etaeta <- dmvnorm(r, mean = rep(0, nq), sigma = pro_Vr, log = T) +
      sum(pro_etaeta - 2*log1pexp(pro_etaeta)) -
      dmvnorm(r, mean = rep(0, nq), sigma = Vr, log = T) -
      sum(etaeta - 2*log1pexp(etaeta))
    
    MH_etaeta <- min(0, MH_etaeta)
    
    if(runif(1,0,1) < exp(MH_etaeta)){
      eta <- pro_eta
      etaeta <- pro_etaeta
      W <- pro_W
      Vr <- pro_Vr
      accepteta <- accepteta + 1
    }
    
    s_etaeta <- s_etaeta + (iter+1)^(-0.7) * (exp(MH_etaeta) - 0.234)
    
    detaeta_2 <- as.vector(etaeta - m_etaeta)
    
    # Covariance
    R_etaeta <- (1-(iter+1)^(-0.7)) * R_etaeta + (iter+1)^(-0.7)*(detaeta_2%*%t(detaeta_2))
    
    # Mean update
    m_etaeta <- m_etaeta + (iter+1)^(-0.7) * detaeta_2
    
    xi_etaeta <- exp(s_etaeta) * R_etaeta
    
    diag(xi_etaeta) <- diag(xi_etaeta) + 1e-6 * max(diag(xi_etaeta))
    
    ################
    ### update A ###
    ################
    
    pro_A <- matrix(0, q, q)
    
    for(i in 1:q){
      for(j in 1:i){
        if(j == i){
          pro_A[i, j] <- exp(rnorm(1, log(A[i, j]), 0.015))
          
        }else{
          pro_A[i, j] <- rnorm(1, A[i, j], 0.015)
        }
      }
    }
    
    Q <- Dinv_new(rho, n, q, W)
    invQ1 <- solve(Q[[1]])
    invQ2 <- solve(Q[[2]])
    invQ3 <- solve(Q[[3]])
    invQ4 <- solve(Q[[4]])
    
    Sigma <-  A %*% t(A)
    
    pro_kprod <- kronecker(pro_A, diag(n))
    invQ <- as.matrix(bdiag(bdiag(invQ1, invQ2),
                            bdiag(invQ3, invQ4)))
    
    pro_Vr <- as.matrix(forceSymmetric(pro_kprod %*% invQ %*% t(pro_kprod)))
    pro_Sigma <-  pro_A %*% t(pro_A)
    
    lpA <- -(nu+q+1) / 2 * logdet(Sigma) - 1/2 * sum(diag(nu * R%*%solve(Sigma))) +
      log(Jacob_A(A)) + sum(log(diag(A)))
    pro_lpA <- -(nu+q+1) / 2 * logdet(pro_Sigma) - 1/2 * sum(diag(nu * R%*%solve(pro_Sigma))) +
      log(Jacob_A(pro_A)) + sum(log(diag(pro_A)))
    
    MH_A <- dmvnorm(r, mean=rep(0, nq), sigma=pro_Vr, log=T) + pro_lpA -
      dmvnorm(r, mean=rep(0, nq), sigma=Vr, log=T) - lpA
    
    MH_A <- min(0, MH_A)
    
    if(runif(1,0,1)<exp(MH_A)){
      A <- pro_A
      Vr <- pro_Vr
      acceptA <- acceptA+1
    }
    
    ######################
    ### record samples ###
    ######################
    
    keeptheta[iter,] <- theta
    keepphi[iter,] <- theta[u]
    keepbeta[iter,] <- beta
    keeptaus[iter] <- taus
    keeprho[iter,] <- rho
    keepu[iter,] <- u
    keepv[iter,] <- v
    keepr[iter,] <- r
    keepFr[iter,] <- F_r
    keepeta[iter,] <- eta
    keepW1[,,iter] <- W[[1]]
    keepW2[,,iter] <- W[[2]]
    keepW3[,,iter] <- W[[3]]
    keepW4[,,iter] <- W[[4]]
    keepA[iter,] <- as.vector(A)[-c(5,9,10,13:15)]
    
  }
  
  list(beta = keepbeta[burn:runs,], phi=keepphi[burn:runs,], 
       theta = keeptheta[burn:runs,], u = keepu[burn:runs,], 
       v = keepv[burn:runs,], r = keepr[burn:runs,], 
       Fr = keepFr[burn:runs,], taus=keeptaus[burn:runs], 
       rho = keeprho[burn:runs,], eta = keepeta[burn:runs,], 
       W1 = keepW1[,,burn:runs], W2 = keepW2[,,burn:runs], 
       W3 = keepW3[,,burn:runs], W4 = keepW4[,,burn:runs],
       A = keepA[burn:runs,],
       acc.beta = acceptbeta/runs, acc.theta = accepttheta/runs, 
       acc.r = acceptr/nq/runs, acc.v = acceptv/runs, 
       acc.rho = acceptrho/runs, acc.eta = accepteta/runs, 
       acc.A = acceptA/runs)
}

Z <- Z[final_perm,final_perm]

# Simulate 50 datasets
for(seed in 1:5){
  print(seed)
  set.seed(seed)
  
  beta1 <- -2
  beta2 <- 2
  beta3 <- 1
  beta4 <- -1
  
  X <- as.matrix(bdiag(bdiag(rep(1,n), rep(1,n)),bdiag(rep(1,n), rep(1,n))))
  
  y1 <- rpois(n_county, exp(beta1 + phi_true1))
  y2 <- rpois(n_county, exp(beta2 + phi_true2))
  y3 <- rpois(n_county, exp(beta3 + phi_true3))
  y4 <- rpois(n_county, exp(beta4 + phi_true4))
  
  yo1 <- y1[final_perm]
  yo2 <- y2[final_perm]
  yo3 <- y3[final_perm]
  yo4 <- y4[final_perm]
  
  Y <- c(yo1, yo2, yo3, yo4)
  
  mcmc_samples <- ARDP_joint_diff(y = Y, X = X, Z = Z, Minc = Minc,
                                  alpha = 1, q = 4, n.atoms = 15, 
                                  runs = 5000, burn = 2500)
  
  filename <- paste0("runs/mcmc_samples_", seed, ".rds")
  saveRDS(mcmc_samples, file = filename)
  
}
