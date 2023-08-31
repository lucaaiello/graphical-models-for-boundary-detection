rm(list = ls())
setwd("C:/Dati/Dottorato/Visiting UCLA/Spatial Disease Mapping/Bibliografia/gao_2022")

# Libraries --------------------------------------------------------------------

library(maps)
#Import California map
ca.county <- map("county","california", fill=TRUE, plot=FALSE)
library(readr)
library(spdep)
library(maptools)
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

rate_5y<- read_csv("data/SIR_adjusted.csv")

covariates <- read_csv("data/covariates.csv")
race <- read_csv("data/race.csv")
sex <- read_csv("data/sex.csv")
insurance <- read_csv("data/insurance.csv")
smoking <- read.csv("data/smoking.csv")
smoking$smoking <- as.numeric(substr(smoking$Cigarette.Smoking.Rate, 1,4))

# Import incidence data for 4 cancers in California ----------------------------

rate_lung <- rate_5y %>% filter(Site.code == "Lung and Bronchus")
rate_esophagus <- rate_5y %>% filter(Site.code == "Esophagus")
rate_larynx <- rate_5y %>% filter(Site.code == "Larynx")
rate_colrect <- rate_5y %>% filter(Site.code == "Colon and Rectum")

county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])

ca.poly <- map2SpatialPolygons(ca.county, IDs=county.ID)
ca.poly$rate_lung <- rate_lung$standard_ratio
ca.poly$rate_esophagus <- rate_esophagus$standard_ratio
ca.poly$rate_larynx <- rate_larynx$standard_ratio
ca.poly$rate_colrect <- rate_colrect$standard_ratio
ca.poly$smoking <- smoking$smoking

ca.coords <- coordinates(ca.poly)

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

ca.neighbors <- poly2nb(ca.poly)
n <- length(ca.neighbors)

Adj <- sapply(ca.neighbors, function(x,n) { v <- rep(0,n); v[x] <- 1; v}, n)
colnames(Adj) <- county.ID

ca.coord <- coordinates(ca.poly)
ca.latrange <- round(quantile(ca.coord[,2],c(0.25,0.75)))
ca.albersproj <- mapproject(ca.coord[,1],ca.coord[,2],
                            projection = "albers",param=ca.latrange)

# Calculate Moran's I for each cancer using albers projection ------------------

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

ggplot(df, aes(r, moran_value)) + geom_point(size = 3) +
  ylab("Moran's I") + facet_wrap(~cancer) + theme_bw() +
  scale_x_continuous(breaks = 1:11) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 20), axis.text = element_text(size = 15),
        strip.text = element_text(size = 15))

# Neighboring counties information ---------------------------------------------

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
X <- as.matrix(bdiag(bdiag(X1[,c(1,2,4,6)], X2[,c(1,2,4,6)]), 
                     bdiag(X3[,c(1,2,4,6)], X4[,c(1,2,4,6)])))

### MADAGAR model with smoking, old and poverty covariates and differences -----

# Metropolis within Gibbs Sampler for MCMC updating 

MADAGAR <- function(y=Y, X = X, Z1 = Z1, Z2 = Z2, Z3 = Z3, E = E, q = 4,
                    Winc = Winc, Minc = Minc,
                    alpha = 1, n.atoms = 15, 
                    runs = 100000, burn = 75000){
  
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
  
  p <- ncol(X)
  
  y1 <- y[1:n]
  y2 <- y[(n+1):(2*n)]
  y3 <- y[(2*n+1):(3*n)]
  y4 <- y[(3*n+1):(4*n)]
  
  E1 <- E[1:n]
  E2 <- E[(n+1):(2*n)]
  E3 <- E[(2*n+1):(3*n)]
  E4 <- E[(3*n+1):(4*n)]
  
  sigmasq_beta <- 10000
  keepbeta <- matrix(0,runs,p)
  keepphi <- matrix(0,runs,nq)
  keeptheta <- matrix(0,runs,n.atoms)
  keepu <- matrix(0,runs,nq)
  keeprho <- matrix(0,runs,q)
  keeptaus <- rep(0,runs)
  keepv <- matrix(0,runs,n.atoms)
  keepr <- matrix(0,runs,nq)
  keepFr <- matrix(0,runs,nq)
  keepeta <- matrix(0,runs,3*q)
  keepW1 <- array(0,dim = c(n,n,runs))
  keepW2 <- array(0,dim = c(n,n,runs))
  keepW3 <- array(0,dim = c(n,n,runs))
  keepW4 <- array(0,dim = c(n,n,runs))
  keepA <- matrix(0,runs,10)
  
  # initial values
  
  theta <- rep(0,n.atoms)
  
  beta <- rep(0,ncol(X))
  
  taus <- 1
  
  c <- 2
  d <- 0.1
  d2 <- 1
  
  v <- rep(.1,n.atoms)
  vv <- log(v) - log(1-v)
  
  probs <- makeprobs(v)
  
  rho <- rep(0.1,q)
  rhorho <- log(rho) - log(1 - rho)
  
  M1 <- as.numeric(-log(0.5)/quantile(Z1[which(Z1!=0)],0.5))
  M2 <- as.numeric(-log(0.5)/quantile(Z2[which(Z2!=0)],0.5))
  M3 <- as.numeric(-log(0.5)/quantile(Z3[which(Z3!=0)],0.5))
  
  M <- c(rep(M1,q),rep(M2,q),rep(M3,q))
  
  eta <- rep(0.1,3*q)
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
  invQ <- as.matrix(bdiag(bdiag(invQ1, invQ2), bdiag(invQ3, invQ4)))
  
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
  R_etaeta <- diag(0.1,3*q)
  
  xi_theta <- exp(s_theta) * R_theta
  xi_beta <- exp(s_beta) * R_beta
  xi_r <- exp(s_r) * R_r
  xi_vv <- exp(s_vv) * R_vv
  xi_rhorho <- exp(s_rhorho) * R_rhorho
  xi_etaeta <- exp(s_etaeta) * R_etaeta
  
  acceptr <- 0
  acceptv <- 0
  acceptrho <- 0
  accepttheta <- 0
  acceptbeta <- 0
  accepteta <- 0
  acceptA <- 0
  
  count <- 0
  afterburn <- 0
  burn <- burn + 1
  
  for(iter in 1:runs){
    
    if(iter %% 100 == 0){
      print(iter)
      print(acceptbeta/(iter-1))
      print(accepttheta/(iter-1))
      print(acceptr/nq/(iter-1))
      print(acceptv/(iter-1))
      print(acceptrho/(iter-1))
      print(accepteta/(iter-1))
      print(acceptA/(iter-1))
    }
    
    ## update beta -------------------------------------------------------------
    
    pro_beta <- as.numeric(rmvnorm(1,beta,xi_beta)) 
    
    MH_beta <- sum(- E * exp(X%*%pro_beta + theta[u]) + 
                     y *(X%*%pro_beta + theta[u])) - 
      sum(- E * exp(X%*%beta + theta[u]) + 
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
    
    ## update theta ------------------------------------------------------------
    
    pro_theta <- rmvnorm(1,theta,xi_theta) 
    
    MH_theta <- sum(- E * exp(X%*%beta + pro_theta[u]) + 
                      y *(X%*%beta + pro_theta[u])) - 
      sum(taus/2*pro_theta^2) - 
      sum(- E * exp(X%*%beta + theta[u]) + 
            y * (X%*%beta + theta[u])) + 
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
    
    ## update r ----------------------------------------------------------------
    
    for (k in 1:nq){
      
      pro_r <- r
      pro_Fr <- F_r
      pro_u <- u
      
      pro_r[k] <- rnorm(1,r[k],0.6)
      pro_Fr[k] <- pnorm(pro_r[k],0,sqrt(Vr[k,k]))
      pro_u[k] <- makeu(pro_Fr[k],probs)
      
      MH_r <- dmvnorm(pro_r, mean = rep(0, nq), sigma = Vr, log = T) +
        dpois(y[k], E[k] * exp(X[k,]%*%beta + theta[pro_u[k]]), log = T) -
        dmvnorm(r, mean = rep(0, nq), sigma = Vr, log = T) -
        dpois(y[k] ,E[k] * exp(X[k,]%*%beta + theta[u[k]]), log = T)
      
      MH_r <- min(0,MH_r)
      
      if(runif(1,0,1) < exp(MH_r)){
        r[k] <- pro_r[k]
        F_r[k] <- pro_Fr[k]
        u[k] <- pro_u[k]
        acceptr <- acceptr+1
      }
      
    }
    
    ## update v ----------------------------------------------------------------
    
    pro_vv <- rmvnorm(1,vv,xi_vv)
    
    pro_v <- apply(pro_vv, 2, inv_trans_par, 0, 1)
    
    pro_probs <- makeprobs(pro_v)
    
    pro_u <- makeu(F_r,pro_probs)
    
    MH_vv <- sum(dbeta(pro_v, 1, alpha, log = T)) +
      sum(dpois(y, E * exp(X%*%beta + theta[pro_u]), log = T)) +
      sum(log(pro_v) + log(1-pro_v)) -
      sum(dbeta(v, 1, alpha, log = T)) -
      sum(dpois(y, E * exp(X%*%beta + theta[u]), log = T)) -
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
    R_vv <- (1-(iter+1)^(-0.7)) * R_vv + 
      (iter+1)^(-0.7) * (dvv_2%*%t(dvv_2))
    
    # Mean update
    m_vv <- m_vv + (iter+1)^(-0.7) * dvv_2
    
    xi_vv <- exp(s_vv) * R_vv
    
    diag(xi_vv) <- diag(xi_vv) + 1e-6 * max(diag(xi_vv))
    
    ## update taus -------------------------------------------------------------
    
    taus <- rgamma(1, shape = 1/2 * n.atoms + c, rate = sum(theta^2)/2 + d2)
    
    ## update rho --------------------------------------------------------------
    
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
    R_rhorho <- (1-(iter+1)^(-0.7)) * R_rhorho + 
      (iter+1)^(-0.7) * (drhorho_2%*%t(drhorho_2))
    
    # Mean update
    m_rhorho <- m_rhorho + (iter+1)^(-0.7) * drhorho_2
    
    xi_rhorho <- exp(s_rhorho) * R_rhorho
    
    diag(xi_rhorho) <- diag(xi_rhorho) + 1e-6 * max(diag(xi_rhorho))
    
    ## update eta --------------------------------------------------------------
    
    pro_etaeta <- rmvnorm(1, etaeta, xi_etaeta)
    
    pro_eta1 <- apply(matrix(pro_etaeta[,1:4],ncol=4), 2, inv_trans_par, 0, M1)
    pro_eta2 <- apply(matrix(pro_etaeta[,5:8],ncol=4), 2, inv_trans_par, 0, M2)
    pro_eta3 <- apply(matrix(pro_etaeta[,9:12],ncol=4), 2, inv_trans_par, 0, M3)
    
    pro_eta <- c(pro_eta1,pro_eta2,pro_eta3)
    
    pro_W <- replicate(q, matrix(0, nrow = n, ncol = n), simplify = FALSE)
    
    indices <- which(Minc == 1, arr.ind = TRUE)
    
    for (d in 1:q) {
      
      for (i in 1:nrow(indices)) {
        row_idx <- indices[i, "row"]
        col_idx <- indices[i, "col"]
        
        pro_W[[d]][row_idx,col_idx] <- as.numeric(exp(-Z1[row_idx,col_idx]*pro_eta1[d]-
                                                       Z2[row_idx,col_idx]*pro_eta2[d]-
                                                       Z3[row_idx,col_idx]*pro_eta3[d]) >= 0.5)
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
    
    MH_eta <- dmvnorm(r,mean = rep(0, nq), sigma = pro_Vr, log = T) +
      sum(pro_etaeta - 2*log1pexp(pro_etaeta)) -
      dmvnorm(r,mean = rep(0, nq), sigma = Vr, log = T) -
      sum(etaeta - 2*log1pexp(etaeta))
    
    MH_eta <- min(0, MH_eta)
    
    if(runif(1,0,1) < exp(MH_eta)){
      eta <- pro_eta
      etaeta <- pro_etaeta
      W <- pro_W
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
    
    ## update A ----------------------------------------------------------------
    
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
    invQ <- as.matrix(bdiag(bdiag(invQ1, invQ2), bdiag(invQ3, invQ4)))
    pro_Vr <- as.matrix(forceSymmetric(pro_kprod %*% invQ %*% t(pro_kprod)))
    pro_Sigma <-  pro_A %*% t(pro_A)
    
    lpA <- -(nu+q+1) / 2 * logdet(Sigma) - 
      1/2 * sum(diag(nu * R%*%solve(Sigma))) +
      log(Jacob_A(A)) + sum(log(diag(A)))
    pro_lpA <- -(nu+q+1) / 2 * logdet(pro_Sigma) - 
      1/2 * sum(diag(nu * R%*%solve(pro_Sigma))) +
      log(Jacob_A(pro_A)) + sum(log(diag(pro_A)))
    
    MH_A <- dmvnorm(r, mean=rep(0, nq), sigma=pro_Vr, log=T) + pro_lpA -
      dmvnorm(r, mean=rep(0, nq), sigma=Vr, log=T) - lpA
    
    MH_A <- min(0, MH_A)
    
    if(runif(1,0,1)<exp(MH_A)){
      A <- pro_A
      Vr <- pro_Vr
      acceptA <- acceptA+1
    }
    
    ## record samples ----------------------------------------------------------
    
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

set.seed(123)

mcmc_samples <- MADAGAR(y=Y, X = X, Z1 = Z1, Z2 = Z2, Z3 = Z3, E = E, 
                        q = 4, Winc = Winc, Minc = Minc,
                        alpha = 1, n.atoms = 15, 
                        runs = 30000, burn = 20000)


samples.mcmc <- mcmc.list(mcmc(data.frame(beta = mcmc_samples$beta,
                                          theta = mcmc_samples$theta,
                                          gamma = mcmc_samples$r,
                                          v = mcmc_samples$v,
                                          taus = mcmc_samples$taus,
                                          rho = mcmc_samples$rho,
                                          etaa = mcmc_samples$eta,
                                          A = mcmc_samples$A)))

# Chain diagnostics ------------------------------------------------------------

samples.ggs <- ggs(samples.mcmc, keep_original_order = TRUE)

ggs_traceplot(samples.ggs, family = "beta", hpd = TRUE) + 
  facet_wrap(~Parameter, ncol = 4, scales = "free")
ggs_autocorrelation(samples.ggs, family = "beta") + 
  facet_wrap(~Parameter, ncol = 4, scales = "free")
ggs_running(samples.ggs, family = "beta") + 
  facet_wrap(~Parameter, ncol = 4, scales = "free")

ggs_traceplot(samples.ggs, family = "theta") + 
  facet_wrap(~Parameter, ncol = 3, scales = "free")
ggs_autocorrelation(samples.ggs, family = "theta") + 
  facet_wrap(~Parameter, ncol = 3, scales = "free")
ggs_running(samples.ggs, family = "theta") + 
  facet_wrap(~Parameter, ncol = 3, scales = "free")

ggs_traceplot(samples.ggs, family = "v") + 
  facet_wrap(~Parameter, ncol = 3, scales = "free")
ggs_autocorrelation(samples.ggs, family = "v") + 
  facet_wrap(~Parameter, ncol = 3, scales = "free")
ggs_running(samples.ggs, family = "v") + 
  facet_wrap(~Parameter, ncol = 3, scales = "free")

ggs_traceplot(samples.ggs, family = "gamma") + 
  facet_wrap(~Parameter, ncol = 5, scales = "free")
ggs_autocorrelation(samples.ggs, family = "gamma") + 
  facet_wrap(~Parameter, ncol = 5, scales = "free")

ggs_traceplot(samples.ggs, family = "taus")
ggs_autocorrelation(samples.ggs, family = "taus")
ggs_running(samples.ggs, family = "taus")

ggs_traceplot(samples.ggs, family = "rho")
ggs_autocorrelation(samples.ggs, family = "rho")
ggs_running(samples.ggs, family = "rho")

ggs_traceplot(samples.ggs, family = "etaa") + 
  facet_wrap(~Parameter, ncol = 3, scales = "free")
ggs_autocorrelation(samples.ggs, family = "etaa") + 
  facet_wrap(~Parameter, ncol = 3, scales = "free")
ggs_running(samples.ggs, family = "etaa") +
  facet_wrap(~Parameter, ncol = 3, scales = "free")

ggs_traceplot(samples.ggs, family = "A") + 
  facet_wrap(~Parameter, ncol = 2, scales = "free")
ggs_autocorrelation(samples.ggs, family = "A") + 
  facet_wrap(~Parameter, ncol = 2, scales = "free")
ggs_running(samples.ggs, family = "A") + 
  facet_wrap(~Parameter, ncol = 2, scales = "free")


# Monte Carlo Standard Error estimation ----------------------------------------

library(mcmcse)

# as.data.frame(samples.mcmc[[1]])

# Batch means estimator
mcerror_bm <- mcse.multi(x = as.data.frame(samples.mcmc[[1]]), 
                         method = "bm", r = 1,
                         size = NULL, g = NULL, adjust = TRUE,
                         blather = TRUE)

mcse.mat(x = as.data.frame(samples.mcmc[[1]]), method = "bm", g = NULL)

ess(as.data.frame(samples.mcmc[[1]]))

plot(confRegion(mcerror_bm, which = c(2,3), level = .95), type = 'l', asp = 1)
lines(confRegion(mcerror_bm, which = c(2,3), level = .90), col = "red")

multiESS(as.data.frame(samples.mcmc[[1]]), covmat = mcerror_bm$cov)

minESS(p = 305, alpha = .05, eps = .05) 
minESS(p = 305, alpha = .05, ess = 1538.565)

# extracting spatial random effects --------------------------------------------

phis <- mcmc_samples$phi
phis1 <- phis[,1:58][,order(final_perm)]
phis2 <- phis[,59:116][,order(final_perm)]
phis3 <- phis[,117:174][,order(final_perm)]
phis4 <- phis[,175:232][,order(final_perm)]
phis_origin <- cbind(phis1, phis2, phis3, phis4)

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
FDR_est1 <- rep(0, length(T_edge1))
FDR_est2 <- rep(0, length(T_edge1))
FDR_est3 <- rep(0, length(T_edge1))
FDR_est4 <- rep(0, length(T_edge1))


for(i in 1:length(threshold1)){
  th1 <- threshold1[i]
  est_diff1 <- as.numeric(pvij1 >= th1)
  FDR_est1[i] <- sum((1-pvij1) * est_diff1)  / sum(est_diff1)
  
  th2 <- threshold2[i]
  est_diff2 <- as.numeric(pvij2 >= th2)
  FDR_est2[i] <- sum((1-pvij2) * est_diff2)  / sum(est_diff2)
  
  th3 <- threshold3[i]
  est_diff3 <- as.numeric(pvij3 >= th3)
  FDR_est3[i] <- sum((1-pvij3) * est_diff3)  / sum(est_diff3)
  
  th4 <- threshold4[i]
  est_diff4 <- as.numeric(pvij4 >= th4)
  FDR_est4[i] <- sum((1-pvij4) * est_diff4)  / sum(est_diff4)
}

plot(FDR_est1,type="l", lty = 1, xlab="Number of edges selected", 
     ylab = "Estimated FDR", ylim=c(0,0.5))
lines(FDR_est2,type="l", lty = 2)
lines(FDR_est3,type="l", lty = 3)
lines(FDR_est4,type="l", lty = 4)
abline(h=0.05, col = "red")
legend("topleft", legend=c("Lung", "Esophageal", "Layrnx", "Colorectal"),
       lty=1:4, cex=0.7)

# ggplot of the previous plot --------------------------------------------------

library(reshape2)

df <- data.frame(
  x = 1:length(FDR_est1),
  FDR_est1 = FDR_est1,
  FDR_est2 = FDR_est2,
  FDR_est3 = FDR_est3,
  FDR_est4 = FDR_est4
)

# Reshape the data to long format
df_long <- melt(df, id.vars = "x", variable.name = "Variable", value.name = "Value")

# Plot the data
ggplot(df_long, aes(x = x, y = Value, linetype = Variable)) +
  geom_line(size = 2) +
  labs(x = "Number of edges selected", y = "Estimated FDR") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash"),
                        labels = c("Lung", "Esophageal", "Larynx", "Colorectal"),
                        name = "Cancer") +
  geom_hline(yintercept = 0.05, color = "red", size = 2) +
  theme_bw() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        plot.title = element_text(size=40, face="bold",hjust = 0.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1,"cm"))


# thresholds -------------------------------------------------------------------

alpha_n <- 0.05
T1 <- sum(FDR_est1<=alpha_n)
T2 <- sum(FDR_est2<=alpha_n)
T3 <- sum(FDR_est3<=alpha_n)
T4 <- sum(FDR_est4<=alpha_n)

est_diff1 <- as.numeric(pvij1 >= threshold1[T1])
name_diff1 <- names(pvij1[order(pvij1,decreasing = T)][1:T1])

est_diff2 <- as.numeric(pvij2 >= threshold2[T2])
name_diff2 <- names(pvij2[order(pvij2,decreasing = T)][1:T1])

est_diff3 <- as.numeric(pvij3 >= threshold3[T3])
name_diff3 <- names(pvij3[order(pvij3,decreasing = T)][1:T1])

est_diff4 <- as.numeric(pvij4 >= threshold4[T4])
name_diff4 <- names(pvij4[order(pvij4,decreasing = T)][1:T1])

# Table for top T1 pairs of neighbors
name_diff = cbind(name_diff1, name_diff2, name_diff3, name_diff4)
colnames(name_diff) = c("Lung", "Esophageal", "Layrnx", "Coloretal")

# Plot boundaries
if (!require(gpclib)) install.packages("gpclib", type="source")
#gpclibPermit()
p <- list()
p_est <- list()

library(plyr)
ca.poly$county <- county.ID
ca.poly$centroid_x <- ca.coords[,1]
ca.poly$centroid_y <- ca.coords[,2]
ca.poly$old <- county_attribute1$V_Persons_age_65_ACS_2012_2016
ca.poly$poverty <- county_attribute1$VFamiliesbelowpovertyACS201220

library(rgeos)

ca.poly@data$id <- rownames(ca.poly@data)
ca.poly.f <- fortify(ca.poly, region = "id")

# ca.poly.f <- coords
ca.poly.df <- join(ca.poly.f, ca.poly@data, by = "id")

# Correct boundary data for plots
path <- list()
for(i in 1: nrow(neighbor_list0)){
  r1 <- ca.poly.df[ca.poly.df$id %in% neighbor_list0[i,1],1:2]
  r2 <- ca.poly.df[ca.poly.df$id %in% neighbor_list0[i,2],1:2]
  edges <- generics::intersect(r1, r2)
  path[[i]] <- edges
}
path[[2]][nrow(path[[2]])+1,] <- path[[2]][1,]
path[[2]] <- path[[2]][-1,]
path[[11]][nrow(path[[11]])+1,] <- path[[11]][1,]
path[[11]] <- path[[11]][-1,]
path[[19]][nrow(path[[19]])+1,] <- path[[19]][1,]
path[[19]] <- path[[19]][-1,]
path[[25]][nrow(path[[25]])+1,] <- path[[25]][1,]
path[[25]] <- path[[25]][-1,]
path[[27]][nrow(path[[27]])+1,] <- path[[27]][1,]
path[[27]] <- path[[27]][-1,]
path[[31]][nrow(path[[31]])+1,] <- path[[31]][1,]
path[[31]] <- path[[31]][-1,]
path[[39]][nrow(path[[39]])+1,] <- path[[39]][1,]
path[[39]] <- path[[39]][-1,]
path[[43]][nrow(path[[43]])+1,] <- path[[43]][1,]
path[[43]] <- path[[43]][-1,]
path[[46]][nrow(path[[46]])+1,] <- path[[46]][1,]
path[[46]] <- path[[46]][-1,]
path[[64]][nrow(path[[64]])+1,] <- path[[64]][1,]
path[[64]] <- path[[64]][-1,]
path[[75]][nrow(path[[75]])+1,] <- path[[75]][1,]
path[[75]] <- path[[75]][-1,]
path[[76]][nrow(path[[76]])+1,] <- path[[76]][1,]
path[[76]] <- path[[76]][-1,]
path[[81]][nrow(path[[81]])+1,] <- path[[81]][1,]
path[[81]] <- path[[81]][-1,]
path[[100]][nrow(path[[100]])+1,] <- path[[100]][1,]
path[[100]] <- path[[100]][-1,]
path[[125]][nrow(path[[125]])+1,] <- path[[125]][1,]
path[[125]] <- path[[125]][-1,]
path[[130]][nrow(path[[130]])+1,] <- path[[130]][1,]
path[[130]] <- path[[130]][-1,]
path[[133]][nrow(path[[133]])+1,] <- path[[133]][1,]
path[[133]] <- path[[133]][-1,]


# California map with county names ---------------------------------------------

ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black", fill = "white", size = 0.5,
               alpha = 0.6) +
  geom_text(data = ca.poly.df, aes(label = county, x = centroid_x, y = centroid_y),
            color = "darkblue", size = 4, family = "sans", hjust = 0.5, vjust = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Plot SIR in areal map --------------------------------------------------------

color <- rev(brewer.pal(5,"RdBu"))

breaks1 <- quantile(ca.poly$rate_lung, c(0.2, 0.4, 0.6, 0.8))
breaks2 <- quantile(ca.poly$rate_esophagus, c(0.2, 0.4, 0.6, 0.8))
breaks3 <- quantile(ca.poly$rate_larynx, c(0.2, 0.4, 0.6, 0.8))
breaks4 <- quantile(ca.poly$rate_colrect, c(0.2, 0.4, 0.6, 0.8))

lung_SIR <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group, fill = rate_lung),
               color = "black",
               alpha = 0.6) +
  scale_fill_stepsn(colors = color,
                    breaks = breaks1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank()) +
  ggtitle("Lung cancer")

esophageal_SIR <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group, fill = rate_esophagus),
               color = "black",
               alpha = 0.6) +
  scale_fill_stepsn(colors = color,
                    breaks = breaks2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank()) +
  ggtitle("Esophageal cancer")

larynx_SIR <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group, fill = rate_larynx),
               color = "black",
               alpha = 0.6) +
  scale_fill_stepsn(colors = color,
                    breaks = breaks3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank()) +
  ggtitle("Larynx cancer")

colorectal_SIR <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group, fill = rate_colrect),
               color = "black",
               alpha = 0.6) +
  scale_fill_stepsn(colors = color,
                    breaks = breaks4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank()) +
  ggtitle("Colorectal cancer")

ggarrange(lung_SIR, esophageal_SIR, larynx_SIR, colorectal_SIR, nrow = 2, ncol = 2)


# Covariates maps --------------------------------------------------------------

brks_fit_smoking <- quantile(ca.poly$smoking, c(0.2, 0.4, 0.6, 0.8))

color <- brewer.pal(5,"YlOrRd")

smoking_plot <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group, fill = smoking),
               color = "black",
               alpha = 0.6) +
  scale_fill_stepsn(colors = color,
                    breaks = brks_fit_smoking) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank()) +
  ggtitle("Smoking rate")


brks_fit_old <- quantile(ca.poly$old, c(0.2, 0.4, 0.6, 0.8))

old_plot <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group, fill = old),
               color = "black",
               alpha = 0.6) +
  scale_fill_stepsn(colors = color,
                    breaks = brks_fit_old) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank()) +
  ggtitle("Over 65 rate")

brks_fit_poverty <- quantile(ca.poly$poverty, c(0.2, 0.4, 0.6, 0.8))

poverty_plot <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group, fill = poverty),
               color = "black",
               alpha = 0.6) +
  scale_fill_stepsn(colors = color,
                    breaks = brks_fit_poverty) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank()) +
  ggtitle("Below poverty threshold rate")

ggarrange(smoking_plot, old_plot, poverty_plot, nrow = 1, ncol = 3)

#saveRDS(path, "path.rds")

# Boundary detection -----------------------------------------------------------

color <- rev(brewer.pal(5,"RdBu"))

ch_edge1 <- pvij1[which(est_diff1 == 1)]
edge_plot1 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group, fill = rate_lung),
               #color = ifelse(ca.poly.df$zero_pos == 1, "white", "black"),
               #alpha = ifelse(ca.poly.df$zero_pos == 1, 0.5, 1),
               color = "black",
               alpha = 0.6
               #fill = "white"
               #linetype = ifelse(ca.poly.df$zero_pos == 1, "dotted", "solid")
  ) +
  scale_fill_stepsn(colors = color,
                    breaks = breaks1)

for(i in which(est_diff1 == 1)){
  edge_plot1 = edge_plot1 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                      size = ((pvij1[i]- min(ch_edge1))/(max(ch_edge1) - min(ch_edge1)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5),
          legend.text = element_text(size = 15),
          legend.key.size = unit(0.75,"cm"),
          legend.title = element_blank()) +
    ggtitle(paste("Lung (T = ", T1, ")", sep=""))
}

edge_plot1

ch_edge2 <- pvij2[which(est_diff2 == 1)]
edge_plot2 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group, fill = rate_esophagus),
               #color = ifelse(ca.poly.df$zero_pos == 1, "white", "black"),
               #alpha = ifelse(ca.poly.df$zero_pos == 1, 0.5, 1),
               color = "black",
               alpha = 0.6
               #fill = "white"
               #linetype = ifelse(ca.poly.df$zero_pos == 1, "dotted", "solid")
  ) + 
  scale_fill_stepsn(colors = color,
                    breaks = breaks2)

for(i in which(est_diff2 == 1)){
  edge_plot2 = edge_plot2 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                      size = ((pvij2[i]- min(ch_edge2))/(max(ch_edge2) - min(ch_edge2)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5),
          legend.text = element_text(size = 15),
          legend.key.size = unit(0.75,"cm"),
          legend.title = element_blank()) +
    ggtitle(paste("Esophageal (T = ", T2, ")", sep=""))
}

edge_plot2

ch_edge3 <- pvij3[which(est_diff3 == 1)]
edge_plot3 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group, fill = rate_larynx),
               #color = ifelse(ca.poly.df$zero_pos == 1, "white", "black"),
               #alpha = ifelse(ca.poly.df$zero_pos == 1, 0.5, 1),
               color = "black",
               alpha = 0.6
               #fill = "white"
               #linetype = ifelse(ca.poly.df$zero_pos == 1, "dotted", "solid")
  ) +
  scale_fill_stepsn(colors = color,
                    breaks = breaks3)

for(i in which(est_diff3 == 1)){
  edge_plot3 = edge_plot3 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                      size = ((pvij3[i]- min(ch_edge3))/(max(ch_edge3) - min(ch_edge3)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5),
          legend.text = element_text(size = 15),
          legend.key.size = unit(0.75,"cm"),
          legend.title = element_blank()) +
    ggtitle(paste("Larynx (T = ", T3, ")", sep=""))
}

edge_plot3

ch_edge4 <- pvij4[which(est_diff4 == 1)]
edge_plot4 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group, fill = rate_colrect),
               #color = ifelse(ca.poly.df$zero_pos == 1, "white", "black"),
               #alpha = ifelse(ca.poly.df$zero_pos == 1, 0.5, 1),
               color = "black",
               alpha = 0.6
               #fill = "white"
               #linetype = ifelse(ca.poly.df$zero_pos == 1, "dotted", "solid")
  ) +
  scale_fill_stepsn(colors = color,
                    breaks = breaks4)

for(i in which(est_diff4 == 1)){
  edge_plot4 = edge_plot4 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                      size = ((pvij4[i]- min(ch_edge4))/(max(ch_edge4) - min(ch_edge4)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5),
          legend.text = element_text(size = 15),
          legend.key.size = unit(0.75,"cm"),
          legend.title = element_blank()) +
    ggtitle(paste("Colorectal (T = ", T4, ")", sep=""))
}

edge_plot4

#pdf("est_diff_cancer_dagar_pois.pdf", height = 10, width = 12)
ggarrange(edge_plot1, edge_plot2, edge_plot3, edge_plot4, nrow = 2, ncol = 2)

# Shared difference boundary for each pair of cancers --------------------------

bc <- function(x){
  diff_bc <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]] & x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]+58]){
      diff_bc <- c(diff_bc, 1)
    }else{
      diff_bc <- c(diff_bc, 0)
    }
  }
  return(diff_bc)
}

vij_samplesc12 <- t(apply(cbind(phis1,phis2), 1, bc))
vij_samplesc13 <- t(apply(cbind(phis1,phis3), 1, bc))
vij_samplesc14 <- t(apply(cbind(phis1,phis4), 1, bc))
vij_samplesc23 <- t(apply(cbind(phis2,phis3), 1, bc))
vij_samplesc24 <- t(apply(cbind(phis2,phis4), 1, bc))
vij_samplesc34 <- t(apply(cbind(phis3,phis4), 1, bc))

# probabilities
pvijc12 <- apply(vij_samplesc12, 2, mean)
pvijc13 <- apply(vij_samplesc13, 2, mean)
pvijc14 <- apply(vij_samplesc14, 2, mean)
pvijc23 <- apply(vij_samplesc23, 2, mean)
pvijc24 <- apply(vij_samplesc24, 2, mean)
pvijc34 <- apply(vij_samplesc34, 2, mean)


# Estimated FDR curves
thresholdc12 <- sort(pvijc12, decreasing = TRUE)[T_edge1]
thresholdc13 <- sort(pvijc13, decreasing = TRUE)[T_edge1]
thresholdc14 <- sort(pvijc14, decreasing = TRUE)[T_edge1]
thresholdc23 <- sort(pvijc23, decreasing = TRUE)[T_edge1]
thresholdc24 <- sort(pvijc24, decreasing = TRUE)[T_edge1]
thresholdc34 <- sort(pvijc34, decreasing = TRUE)[T_edge1]

FDR_estc12 <- rep(0, length(T_edge1))
FDR_estc13 <- rep(0, length(T_edge1))
FDR_estc14 <- rep(0, length(T_edge1))
FDR_estc23 <- rep(0, length(T_edge1))
FDR_estc24 <- rep(0, length(T_edge1))
FDR_estc34 <- rep(0, length(T_edge1))


for(i in 1:length(threshold1)){
  thc12 <- thresholdc12[i]
  est_diffc12 <- as.numeric(pvijc12 >= thc12)
  FDR_estc12[i] <- sum((1-pvijc12) * est_diffc12)  / sum(est_diffc12)
  
  thc13 <- thresholdc13[i]
  est_diffc13 <- as.numeric(pvijc13 >= thc13)
  FDR_estc13[i] <- sum((1-pvijc13) * est_diffc13)  / sum(est_diffc13)
  
  thc14 <- thresholdc14[i]
  est_diffc14 <- as.numeric(pvijc14 >= thc14)
  FDR_estc14[i] <- sum((1-pvijc14) * est_diffc14)  / sum(est_diffc14)
  
  thc23 <- thresholdc23[i]
  est_diffc23 <- as.numeric(pvijc23 >= thc23)
  FDR_estc23[i] <- sum((1-pvijc23) * est_diffc23)  / sum(est_diffc23)
  
  thc24 <- thresholdc24[i]
  est_diffc24 <- as.numeric(pvijc24 >= thc24)
  FDR_estc24[i] <- sum((1-pvijc24) * est_diffc24)  / sum(est_diffc24)
  
  thc34 <- thresholdc34[i]
  est_diffc34 <- as.numeric(pvijc34 >= thc34)
  FDR_estc34[i] <- sum((1-pvijc34) * est_diffc34)  / sum(est_diffc34)
  
}

# Use the same threshold as above and identify shared boundaries

alpha_nc <- 0.05
T12c <- sum(FDR_estc12<=alpha_nc)
T13c <- sum(FDR_estc13<=alpha_nc)
T14c <- sum(FDR_estc14<=alpha_nc)
T23c <- sum(FDR_estc23<=alpha_nc)
T24c <- sum(FDR_estc24<=alpha_nc)
T34c <- sum(FDR_estc34<=alpha_nc)
#Tallc = sum(FDR_estcall<=alpha_nc)


est_diffc12 <- as.numeric(pvijc12 >= thresholdc12[T12c])
neighbor_list_diffc12 <- neighbor_list0[est_diffc12 == 1, ]

est_diffc13 <- as.numeric(pvijc13 >= thresholdc13[T13c])
neighbor_list_diffc13 <- neighbor_list0[est_diffc13 == 1, ]

est_diffc14 <- as.numeric(pvijc14 >= thresholdc14[T14c])
neighbor_list_diffc14 <- neighbor_list0[est_diffc14 == 1, ]

est_diffc23 <- as.numeric(pvijc23 >= thresholdc23[T23c])
neighbor_list_diffc23 <- neighbor_list0[est_diffc23 == 1, ]

est_diffc24 <- as.numeric(pvijc24 >= thresholdc24[T24c])
neighbor_list_diffc24 <- neighbor_list0[est_diffc24 == 1, ]

est_diffc34 <- as.numeric(pvijc34 >= thresholdc34[T34c])
neighbor_list_diffc34 <- neighbor_list0[est_diffc34 == 1, ]

ch_edge12 <- pvijc12[est_diffc12==1]
edge_plotc12 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )

for(i in which(est_diffc12==1)){
  edge_plotc12 <- edge_plotc12 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                           size = ((pvijc12[i]- min(ch_edge12))/(max(ch_edge12) - min(ch_edge12)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Lung, Esophageal (T = ", T12c, ")", sep=""))
}

edge_plotc12

ch_edge13 <- pvijc13[est_diffc13==1]
edge_plotc13 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )

for(i in which(est_diffc13==1)){
  edge_plotc13 <- edge_plotc13 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                           size = ((pvijc13[i]- min(ch_edge13))/(max(ch_edge13) - min(ch_edge13)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Lung, Layrnx (T = ", T13c, ")", sep=""))
}

edge_plotc13

ch_edge14 <- pvijc14[est_diffc14==1]
edge_plotc14 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )

for(i in which(est_diffc14==1)){
  edge_plotc14 <- edge_plotc14 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                           size = ((pvijc14[i]- min(ch_edge14))/(max(ch_edge14) - min(ch_edge14)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Lung, Colorectal (T = ", T14c, ")", sep=""))
}
edge_plotc14

ch_edge23 <- pvijc23[est_diffc23==1]
edge_plotc23 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )

for(i in which(est_diffc23==1)){
  edge_plotc23 <- edge_plotc23 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                           size = ((pvijc23[i]- min(ch_edge23))/(max(ch_edge23) - min(ch_edge23)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Esophageal, Layrnx (T = ", T23c, ")", sep=""))
}
edge_plotc23

ch_edge24 <- pvijc24[est_diffc24==1]
edge_plotc24 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )

for(i in which(est_diffc24==1)){
  edge_plotc24 <- edge_plotc24 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                           size = ((pvijc24[i]- min(ch_edge24))/(max(ch_edge24) - min(ch_edge24)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Esophageal, Colorectal (T = ", T24c, ")", sep=""))
}
edge_plotc24

ch_edge34 <- pvijc34[est_diffc34==1]
edge_plotc34 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )

for(i in which(est_diffc34==1)){
  edge_plotc34 <- edge_plotc34 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                           size = ((pvijc34[i]- min(ch_edge34))/(max(ch_edge34) - min(ch_edge34)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Larynx, Colorectal (T = ", T34c, ")", sep=""))
}
edge_plotc34


#pdf("est_share_dagar_pois.pdf", height = 10, width = 15)
ggarrange(edge_plotc12, edge_plotc13, edge_plotc14, 
          edge_plotc23, edge_plotc24, edge_plotc34, nrow = 2, ncol = 3)
#dev.off()

# Mutual cross diseases difference boundary ------------------------------------

bc1 <- function(x){
  diff_bc1 <- NULL
  for(i in 1:nrow(neighbor_list0)){
    if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]+58] &
       x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]]){
      diff_bc1 <- c(diff_bc1, 1)
    }else{
      diff_bc1 <- c(diff_bc1, 0)
    }
  }
  return(diff_bc1)
}

vij_samples12 <- t(apply(cbind(phis1,phis2), 1, bc1))
vij_samples13 <- t(apply(cbind(phis1,phis3), 1, bc1))
vij_samples14 <- t(apply(cbind(phis1,phis4), 1, bc1))
vij_samples23 <- t(apply(cbind(phis2,phis3), 1, bc1))
vij_samples24 <- t(apply(cbind(phis2,phis4), 1, bc1))
vij_samples34 <- t(apply(cbind(phis3,phis4), 1, bc1))

#probability
pvij12 <- apply(vij_samples12, 2, mean)
pvij13 <- apply(vij_samples13, 2, mean)
pvij14 <- apply(vij_samples14, 2, mean)
pvij23 <- apply(vij_samples23, 2, mean)
pvij24 <- apply(vij_samples24, 2, mean)
pvij34 <- apply(vij_samples34, 2, mean)

# Estimated FDR curves
threshold12 <- sort(pvij12, decreasing = TRUE)[T_edge1]
threshold13 <- sort(pvij13, decreasing = TRUE)[T_edge1]
threshold14 <- sort(pvij14, decreasing = TRUE)[T_edge1]
threshold23 <- sort(pvij23, decreasing = TRUE)[T_edge1]
threshold24 <- sort(pvij24, decreasing = TRUE)[T_edge1]
threshold34 <- sort(pvij34, decreasing = TRUE)[T_edge1]

FDR_est12 <- rep(0, length(T_edge1))
FDR_est13 <- rep(0, length(T_edge1))
FDR_est14 <- rep(0, length(T_edge1))
FDR_est23 <- rep(0, length(T_edge1))
FDR_est24 <- rep(0, length(T_edge1))
FDR_est34 <- rep(0, length(T_edge1))

for(i in 1:length(threshold1)){
  th12 <- threshold12[i]
  est_diff12<- as.numeric(pvij12 >= th12)
  FDR_est12[i] <- sum((1-pvij12) * est_diff12)  / sum(est_diff12)
  
  th13 <- threshold13[i]
  est_diff13 <- as.numeric(pvij13 >= th13)
  FDR_est13[i] <- sum((1-pvij13) * est_diff13)  / sum(est_diff13)
  
  th14 <- threshold14[i]
  est_diff14 <- as.numeric(pvij14 >= th14)
  FDR_est14[i] <- sum((1-pvij14) * est_diff14)  / sum(est_diff14)
  
  th23 <- threshold23[i]
  est_diff23 <- as.numeric(pvij23 >= th23)
  FDR_est23[i] <- sum((1-pvij23) * est_diff23)  / sum(est_diff23)
  
  th24 <- threshold24[i]
  est_diff24 <- as.numeric(pvij24 >= th24)
  FDR_est24[i] <- sum((1-pvij24) * est_diff24)  / sum(est_diff24)
  
  th34 <- threshold34[i]
  est_diff34 <- as.numeric(pvij34 >= th34)
  FDR_est34[i] <- sum((1-pvij34) * est_diff34)  / sum(est_diff34)
  
}


# Use the same threshold as above and identify cross-disease boundaries

alpha_n1 <- 0.05
T12 <- sum(FDR_est12<=alpha_n1)
T13 <- sum(FDR_est13<=alpha_n1)
T14 <- sum(FDR_est14<=alpha_n1)
T23 <- sum(FDR_est23<=alpha_n1)
T24 <- sum(FDR_est24<=alpha_n1)
T34 <- sum(FDR_est34<=alpha_n1)

est_diff12 <- as.numeric(pvij12 >= threshold12[T12])
est_diff13 <- as.numeric(pvij13 >= threshold13[T13])
est_diff14 <- as.numeric(pvij14 >= threshold23[T14])
est_diff23 <- as.numeric(pvij23 >= threshold23[T23])
est_diff24 <- as.numeric(pvij24 >= threshold24[T24])
est_diff34 <- as.numeric(pvij34 >= threshold34[T34])

ch_edge12 <- pvij12[which(est_diff12 == 1)]
edge_plot12 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diff12==1)){
  edge_plot12 <- edge_plot12 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                         size = ((pvij12[i]- min(ch_edge12))/(max(ch_edge12) - min(ch_edge12)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Lung, Esophageal (T = ", T12, ")", sep=""))
}
edge_plot12

ch_edge13 <- pvij13[which(est_diff13 == 1)]
edge_plot13 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diff13==1)){
  edge_plot13 <- edge_plot13 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                         size = ((pvij13[i]- min(ch_edge13))/(max(ch_edge13) - min(ch_edge13)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Lung, Layrnx (T = ", T13, ")", sep=""))
}
edge_plot13

ch_edge14 <- pvij14[which(est_diff14 == 1)]
edge_plot14 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diff14==1)){
  edge_plot14 <- edge_plot14 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                         size = ((pvij14[i]- min(ch_edge14))/(max(ch_edge14) - min(ch_edge14)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Lung, Colorectal (T = ", T14, ")", sep=""))
}
edge_plot14

ch_edge23 <- pvij23[which(est_diff23 == 1)]
edge_plot23 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diff23==1)){
  edge_plot23 <- edge_plot23 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                         size = ((pvij23[i]- min(ch_edge23))/(max(ch_edge23) - min(ch_edge23)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Esophageal, Layrnx (T = ", T23, ")", sep=""))
}
edge_plot23

ch_edge24 <- pvij24[which(est_diff24 == 1)]
edge_plot24 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diff24==1)){
  edge_plot24 <- edge_plot24 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                         size = ((pvij24[i]- min(ch_edge24))/(max(ch_edge24) - min(ch_edge24)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Esophageal, Colorectal (T = ", T24, ")", sep=""))
}
edge_plot24

ch_edge34 <- pvij34[which(est_diff34 == 1)]
edge_plot34 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white"
  )
for(i in which(est_diff34==1)){
  edge_plot34 <- edge_plot34 + geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "red",
                                         size = ((pvij34[i]- min(ch_edge34))/(max(ch_edge34) - min(ch_edge34)))*0.5 + 0.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
    ggtitle(paste("Larynx, Colorectal (T = ", T34, ")", sep=""))
}
edge_plot34

#pdf("est_crossd_car_pois.pdf", height = 10, width = 15)
ggarrange(edge_plot12, edge_plot13, edge_plot14, 
          edge_plot23, edge_plot24, edge_plot34, nrow = 2, ncol = 3)
#dev.off()

# Adjacency modeling -----------------------------------------------------------

W1 <- mcmc_samples$W1
W2 <- mcmc_samples$W2
W3 <- mcmc_samples$W3
W4 <- mcmc_samples$W4

W1.mean <- apply(W1, c(1,2), mean)

W2.mean <- apply(W2, c(1,2), mean)

W3.mean <- apply(W3, c(1,2), mean)

W4.mean <- apply(W4, c(1,2), mean)

W1.mean.aux <- W1.mean[order(final_perm),order(final_perm)]
W2.mean.aux <- W2.mean[order(final_perm),order(final_perm)]
W3.mean.aux <- W3.mean[order(final_perm),order(final_perm)]
W4.mean.aux <- W4.mean[order(final_perm),order(final_perm)]

adj_det1 <- NULL
adj_det2 <- NULL
adj_det3 <- NULL
adj_det4 <- NULL

for(i in 1:nrow(neighbor_list0)){
  adj_det1 <- c(adj_det1,W1.mean.aux[neighbor_list0[i,1],neighbor_list0[i,2]])
  adj_det2 <- c(adj_det2,W2.mean.aux[neighbor_list0[i,1],neighbor_list0[i,2]])
  adj_det3 <- c(adj_det3,W3.mean.aux[neighbor_list0[i,1],neighbor_list0[i,2]])
  adj_det4 <- c(adj_det4,W4.mean.aux[neighbor_list0[i,1],neighbor_list0[i,2]])
}

# breaks_smoking <- quantile(ca.poly$smoking, c(0.2,0.4,0.6,0.8))
# 
# color <- brewer.pal(5,"YlOrRd")

adj_plot1 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white",
               alpha = 0.6)

for(i in 1:nrow(neighbor_list0)){
  if(adj_det1[i]!=1){
    adj_plot1 <- adj_plot1 + 
      geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "blue",
                size = (1-adj_det1[i])*1.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
      ggtitle("Lung non-adjacencies")  
  }
}

adj_plot1

adj_plot2 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white",
               alpha = 0.6)

for(i in 1:nrow(neighbor_list0)){
  if(adj_det2[i]!=1){
    adj_plot2 <- adj_plot2 + 
      geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "blue",
                size = (1-adj_det2[i])*1.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
      ggtitle("Esophageal non-adjacencies")
  }
}

adj_plot2

adj_plot3 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white",
               alpha = 0.6)

for(i in 1:nrow(neighbor_list0)){
  if(adj_det3[i]!=1){
    adj_plot3 <- adj_plot3 + 
      geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "blue",
                size = (1-adj_det3[i])*1.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
      ggtitle("Larynx non-adjacencies")  
  }
}

adj_plot3

adj_plot4 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white",
               alpha = 0.6)

for(i in 1:nrow(neighbor_list0)){
  if(adj_det4[i]!=1){
    adj_plot4 <- adj_plot4 + 
      geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "blue",
                size = (1-adj_det4[i])*1.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
      ggtitle("Colorectal non-adjacencies")  
  }
}

adj_plot4

ggarrange(adj_plot1,adj_plot2,adj_plot3,adj_plot4, nrow = 2, ncol = 2)
