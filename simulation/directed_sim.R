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

# disease graph adjacency matrices

W_dis <- readRDS("W_dis_try.rds")

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
ARDP_joint_diff <- function(y = Y, X = X, Z = Z, Minc = Minc, W_dis = W_dis,
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
  
  sigmasq_beta <- 1000
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
  keepVr <- array(0,dim = c(nq,nq,runs))
  keepalphas <- matrix(0,runs,8)
  
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
  
  rho <- rep(0.1,q)
  rhorho <- log(rho) - log(1 - rho)
  
  M <- as.numeric(-log(0.5)/quantile(Z[which(Z!=0)],0.5))
  
  eta <- rep(0.1,q)
  etaeta <- log(eta) - log(M - eta)
  
  W <- replicate(q, Minc, simplify = FALSE)
  
  Q <- Dinv_new(rho, n, q, W)
  invQ1 <- solve(Q[[1]])
  invQ2 <- solve(Q[[2]])
  invQ3 <- solve(Q[[3]])
  invQ4 <- solve(Q[[4]])
  
  A <- matrix(0,nq,nq)
  
  alpha0 <- cbind(c(0,0.25,0,0.25),c(0,0,0.25,0),c(0,0,0,0.25),c(0,0,0,0))
  alpha1 <- cbind(c(0,0.25,0,0.25),c(0,0,0.25,0),c(0,0,0,0.25),c(0,0,0,0))
  
  indices <- which(W_dis == 1, arr.ind = TRUE)
  
  # Iterate over the indices
  for (i in 1:nrow(indices)) {
    
    row_idx <- indices[i, "row"]
    col_idx <- indices[i, "col"]
    
    A[((row_idx-1)*n + 1):(row_idx*n),((col_idx-1)*n + 1):(col_idx*n)] <- 
      alpha0[row_idx,col_idx]*diag(n) + alpha1[row_idx,col_idx]*Minc
    
  }
  
  invQ <- as.matrix(bdiag(bdiag(invQ1, invQ2), bdiag(invQ3, invQ4)))
  
  Vr <- as.matrix(forceSymmetric(solve(diag(nq)-A) %*% invQ %*% solve(diag(nq)-t(A))))
  
  r <- rmvnorm(1, rep(0, nq), Vr)
  
  F_r <- pnorm(r,0,sqrt(diag(Vr)))
  
  u <- makeu(F_r,probs)
  phi <- theta[u]
  
  nu <- 2
  R <- 0.1 * diag(q)
  
  s_theta <- 1
  s_beta <- 1
  s_r <- 1
  s_vv <- 1
  s_rhorho <- 1e-8
  s_etaeta <- 1e-8
  
  m_theta <- theta
  m_beta <- beta
  m_r <- r
  m_vv <- vv
  m_rhorho <- rhorho
  m_etaeta <- etaeta
  
  R_theta <- diag(1,n.atoms)
  R_beta <- diag(1,p)
  R_r <- diag(1,nq)
  R_vv <- diag(1,n.atoms)
  R_rhorho <- diag(1e-8,q)
  R_etaeta <- diag(1e-8,q)
  
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
      print("#########")
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
    R_theta <- (1-(iter+1)^(-0.7))*R_theta + (iter+1)^(-0.7)*(dtheta_2%*%t(dtheta_2))
    
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
      
      pro_r[k] <- rnorm(1,r[k],2.5)
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
    
    pro_invQ <- as.matrix(bdiag(bdiag(pro_invQ1, pro_invQ2), 
                                bdiag(pro_invQ3, pro_invQ4)))
    
    pro_Vr <- as.matrix(forceSymmetric(solve(diag(nq)-A) %*% pro_invQ %*% solve(diag(nq)-t(A))))
    
    MH_rhorho <- dmvnorm(r, mean = rep(0, nq), sigma = pro_Vr, log = T) +
      sum(log(pro_rho) + log(1-pro_rho)) - 
      dmvnorm(r, mean = rep(0, nq), sigma = Vr, log = T) - 
      sum(log(rho) + log(1-rho)) 
    
    MH_rhorho <- min(0, MH_rhorho)
    
    if(runif(1,0,1) < exp(MH_rhorho)){
      rhorho <- pro_rhorho
      rho <- pro_rho
      Q <- pro_Q
      invQ <- pro_invQ
      Vr <- pro_Vr
      acceptrho <- acceptrho + 1
    }
    
    s_rhorho <- s_rhorho + (iter+1)^(-0.7) * (exp(MH_rhorho) - 0.234)
    
    drhorho_2 <- as.vector(rhorho - m_rhorho)
    
    # Covariance
    R_rhorho <- (1-(iter+1)^(-0.7)) * R_rhorho + 
      (iter+1)^(-0.7) * (drhorho_2%*%t(drhorho_2))
    
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
    
    pro_invQ <- as.matrix(bdiag(bdiag(pro_invQ1, pro_invQ2), 
                                bdiag(pro_invQ3, pro_invQ4)))
    
    pro_Vr <- as.matrix(forceSymmetric(solve(diag(nq)-A) %*% pro_invQ %*% solve(diag(nq)-t(A))))
    
    MH_eta <- dmvnorm(r,mean = rep(0, nq), sigma = pro_Vr, log = T) +
      sum(pro_etaeta - 2*log1pexp(pro_etaeta)) -
      dmvnorm(r,mean = rep(0, nq), sigma = Vr, log = T) -
      sum(etaeta - 2*log1pexp(etaeta))
    
    MH_eta <- min(0, MH_eta)
    
    if(runif(1,0,1) < exp(MH_eta)){
      eta <- pro_eta
      etaeta <- pro_etaeta
      W <- pro_W
      Q <- pro_Q
      invQ <- pro_invQ
      Vr <- pro_Vr
      accepteta <- accepteta + 1
    }
    
    s_etaeta <- s_etaeta + (iter+1)^(-0.7) * (exp(MH_eta) - 0.234)
    
    detaeta_2 <- as.vector(etaeta - m_etaeta)
    
    # Covariance
    R_etaeta <- (1-(iter+1)^(-0.7)) * R_etaeta + 
      (iter+1)^(-0.7)*(detaeta_2%*%t(detaeta_2))
    
    # Mean update
    m_etaeta <- m_etaeta + (iter+1)^(-0.7) * detaeta_2
    
    xi_etaeta <- exp(s_etaeta) * R_etaeta
    
    diag(xi_etaeta) <- diag(xi_etaeta) + 1e-6 * max(diag(xi_etaeta))
    
    #####################
    ### update alphas ###
    #####################
    
    # proposal parameters for alpha21
    
    F1 <- Minc%*%r[1:n]
    
    F1_mat <- as.matrix(cbind(r[1:n], F1))
    
    H2 <- solve(t(F1_mat)%*%Q[[2]]%*%F1_mat + diag(rep(1/100,2)))
    
    h2 <- t(F1_mat)%*%Q[[2]]%*%r[(n+1):(2*n)]
    
    alpha21 <- rmvnorm(1, mean = H2%*%h2, sigma = H2)
    
    # proposal parameters for alpha32
    
    F2 <- Minc%*%r[(n+1):(2*n)]
    
    F2_mat <- as.matrix(cbind(r[(n+1):(2*n)], F2))
    
    H3 <- solve(t(F2_mat)%*%Q[[3]]%*%F2_mat + diag(rep(1/100,2)))
    
    h3 <- t(F2_mat)%*%Q[[3]]%*%r[(2*n+1):(3*n)]
    
    alpha32 <- rmvnorm(1, mean = H3%*%h3, sigma = H3)
    
    # proposal parameters for c(alpha41,alpha43)
    
    F3 <- Minc%*%r[(2*n+1):(3*n)]
    
    F13_mat <- as.matrix(cbind(F1_mat, r[(2*n+1):(3*n)], F3))
    
    H4 <- solve(t(F13_mat)%*%Q[[4]]%*%F13_mat + diag(rep(1/100,4)))
    
    h4 <- t(F13_mat)%*%Q[[4]]%*%r[(3*n+1):(4*n)]
    
    alpha4143 <- rmvnorm(1, mean = H4%*%h4, sigma = H4)
    
    # matrix construction
    
    alphas <- c(alpha21,alpha32,alpha4143)
    
    alpha0 <- cbind(c(0,alpha21[1],0,alpha4143[1]),c(0,0,alpha32[1],0),
                    c(0,0,0,alpha4143[3]),c(0,0,0,0))
    alpha1 <- cbind(c(0,alpha21[2],0,alpha4143[2]),c(0,0,alpha32[2],0),
                    c(0,0,0,alpha4143[4]),c(0,0,0,0))
    
    indices <- which(W_dis == 1, arr.ind = TRUE)
    
    A <- matrix(0, nq, nq)
    
    # Iterate over the indices
    for (i in 1:nrow(indices)) {
      
      row_idx <- indices[i, "row"]
      col_idx <- indices[i, "col"]
      
      A[((row_idx-1)*n + 1):(row_idx*n),((col_idx-1)*n + 1):(col_idx*n)] <-
        alpha0[row_idx,col_idx]*diag(n) + alpha1[row_idx,col_idx]*Minc
      
    }
    
    Vr <- as.matrix(forceSymmetric(solve(diag(nq)-A) %*% invQ %*% solve(diag(nq)-t(A))))
    
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
    keepVr[,,iter] <- Vr
    keepalphas[iter,] <- alphas
    
  }
  
  list(beta = keepbeta[burn:runs,], phi=keepphi[burn:runs,], 
       theta = keeptheta[burn:runs,], u = keepu[burn:runs,], 
       v = keepv[burn:runs,], r = keepr[burn:runs,], 
       Fr = keepFr[burn:runs,], taus=keeptaus[burn:runs],
       rho = keeprho[burn:runs,], eta = keepeta[burn:runs,], 
       alphas = keepalphas[burn:runs,],
       W1 = keepW1[,,burn:runs], W2 = keepW2[,,burn:runs],
       W3 = keepW3[,,burn:runs], W4 = keepW4[,,burn:runs],
       Vr = keepVr[,,burn:runs],
       acc.beta = acceptbeta/runs, acc.theta = accepttheta/runs, 
       acc.r = acceptr/nq/runs, acc.v = acceptv/runs, acc.rho = acceptrho/runs, 
       acc.eta = accepteta/runs)
}

T_edge <- c(60, 65, 70, 75, 80, 85, 90, 95, 100, 105) 
spec1 <- matrix(0, 50, length(T_edge))
sens1 <- matrix(0, 50, length(T_edge))
spec2 <- matrix(0, 50, length(T_edge))
sens2 <- matrix(0, 50, length(T_edge))
spec3 <- matrix(0, 50, length(T_edge))
sens3 <- matrix(0, 50, length(T_edge))
spec4 <- matrix(0, 50, length(T_edge))
sens4 <- matrix(0, 50, length(T_edge))

specW1 <- rep(0, 50)
sensW1 <- rep(0, 50)
specW2 <- rep(0, 50)
sensW2 <- rep(0, 50)
specW3 <- rep(0, 50)
sensW3 <- rep(0, 50)
specW4 <- rep(0, 50)
sensW4 <- rep(0, 50)

mcmc_samples <- list()

Z <- Z[final_perm,final_perm]

# Simulate 50 datasets
for(seed in 1:50){
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
  
  mcmc_samples[[seed]] <- ARDP_joint_diff(y = Y, X = X, Z = Z, Minc = Minc, W_dis = W_dis,
                                          alpha = 1, q = 4, n.atoms = 15, 
                                          runs = 10000, burn = 5000)
  
  phis <-  mcmc_samples[[seed]]$phi
  phis1 <- phis[,1:58][,order(final_perm)]
  phis2 <- phis[,59:116][,order(final_perm)]
  phis3 <- phis[,117:174][,order(final_perm)]
  phis4 <- phis[,175:232][,order(final_perm)]
  phis_origin <- cbind(phis1, phis2, phis3, phis4)
  
  W1 <- mcmc_samples[[seed]]$W1[order(final_perm),order(final_perm),]
  W2 <- mcmc_samples[[seed]]$W2[order(final_perm),order(final_perm),]
  W3 <- mcmc_samples[[seed]]$W3[order(final_perm),order(final_perm),]
  W3 <- mcmc_samples[[seed]]$W4[order(final_perm),order(final_perm),]
  
  W1.mean <- apply(W1, c(1,2), mean)
  W2.mean <- apply(W2, c(1,2), mean)
  W3.mean <- apply(W3, c(1,2), mean)
  W4.mean <- apply(W4, c(1,2), mean)
  
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
  
  vij_samples1 <- t(apply(phis1, 1, bd))
  vij_samples2 <- t(apply(phis2, 1, bd))
  vij_samples3 <- t(apply(phis3, 1, bd))
  vij_samples4 <- t(apply(phis4, 1, bd))
  
  # probabilities
  
  pvij1 <- apply(vij_samples1, 2, mean)
  pvij2 <- apply(vij_samples2, 2, mean)
  pvij3 <- apply(vij_samples3, 2, mean)
  pvij4 <- apply(vij_samples4, 2, mean)
  
  names(pvij1) <- neighbor_name
  names(pvij2) <- neighbor_name
  names(pvij3) <- neighbor_name
  names(pvij4) <- neighbor_name
  
  # Estimated FDR curves
  T_edge1 <- seq(0, 135, 1)[-1]
  threshold1 <- sort(pvij1, decreasing = TRUE)[T_edge1]
  threshold2 <- sort(pvij2, decreasing = TRUE)[T_edge1]
  threshold3 <- sort(pvij3, decreasing = TRUE)[T_edge1]
  threshold4 <- sort(pvij4, decreasing = TRUE)[T_edge1]
  
  pvijm <- rbind(pvij1, pvij2, pvij3, pvij4)
  
  # calculate sensitivity and specificity on boundary detection
  for(i in 1:length(T_edge)){
    est_diff1 <- factor(as.numeric(pvijm[1,] >= threshold1[i]), levels = c(0,1))
    conf_matrix1 <- table(est_diff1,true_diff1)
    spec1[seed,i] <- sensitivity(conf_matrix1)
    sens1[seed,i] <- specificity(conf_matrix1)
    
    est_diff2 <- factor(as.numeric(pvijm[2,] >= threshold2[i]), levels = c(0,1))
    conf_matrix2 <- table(est_diff2,true_diff2)
    spec2[seed,i] <- sensitivity(conf_matrix2)
    sens2[seed,i] <- specificity(conf_matrix2)
    
    est_diff3 <- factor(as.numeric(pvijm[3,] >= threshold3[i]), levels = c(0,1))
    conf_matrix3 <- table(est_diff3,true_diff3)
    spec3[seed,i] <- sensitivity(conf_matrix3)
    sens3[seed,i] <- specificity(conf_matrix3)
    
    est_diff4 <- factor(as.numeric(pvijm[4,] >= threshold4[i]), levels = c(0,1))
    conf_matrix4 <- table(est_diff4,true_diff4)
    spec4[seed,i] <- sensitivity(conf_matrix4)
    sens4[seed,i] <- specificity(conf_matrix4)
  }
  
  # calculate sensitivity and specificity on adjacency
  
  W <- replicate(4, Minc, simplify = FALSE)
  
  eta.mean <- apply(mcmc_samples[[seed]]$eta, 2, mean)
  
  Winc <- Minc[order(final_perm),order(final_perm)]
  Winc[upper.tri(Winc)] <- 0
  
  indices <- which(Winc == 1, arr.ind = TRUE)
  
  est_adj <- matrix(0,nrow = nrow(indices), ncol = 4)
  
  for (d in 1:4) {
    
    for (i in 1:nrow(indices)) {
      row_idx <- indices[i, "row"]
      col_idx <- indices[i, "col"]
      
      W[[d]][row_idx,col_idx] <- as.numeric(exp(-Z[order(final_perm),order(final_perm)][row_idx,col_idx]*eta.mean[d]) >= 0.5)
      
      est_adj[i,d] <- W[[d]][row_idx,col_idx]
    }
    
  }
  
  est_adj1 <- factor(est_adj[,1])
  conf_matrixW1 <- table(est_adj1,true_adj1)
  specW1[seed] <- sensitivity(conf_matrixW1)
  sensW1[seed] <- specificity(conf_matrixW1)
  
  est_adj2 <- factor(est_adj[,2])
  conf_matrixW2 <- table(est_adj2,true_adj2)
  specW2[seed] <- sensitivity(conf_matrixW2)
  sensW2[seed] <- specificity(conf_matrixW2)
  
  est_adj3 <- factor(est_adj[,3])
  conf_matrixW3 <- table(est_adj3,true_adj3)
  specW3[seed] <- sensitivity(conf_matrixW3)
  sensW3[seed] <- specificity(conf_matrixW3)
  
  est_adj4 <- factor(est_adj[,4])
  conf_matrixW4 <- table(est_adj4,true_adj4)
  specW4[seed] <- sensitivity(conf_matrixW4)
  sensW4[seed] <- specificity(conf_matrixW4)
  
}
