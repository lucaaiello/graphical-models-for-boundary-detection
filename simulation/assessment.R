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

phi_true1 <- readRDS("RE_generation_DAGAR/phi_true1.rds")
phi_true2 <- readRDS("RE_generation_DAGAR/phi_true2.rds")
phi_true3 <- readRDS("RE_generation_DAGAR/phi_true3.rds")
phi_true4 <- readRDS("RE_generation_DAGAR/phi_true4.rds")

# spatial adjacency matrix for diseases

W_true1 <- readRDS("RE_generation_DAGAR/W_true1.rds")
W_true2 <- readRDS("RE_generation_DAGAR/W_true2.rds")
W_true3 <- readRDS("RE_generation_DAGAR/W_true3.rds")
W_true4 <- readRDS("RE_generation_DAGAR/W_true4.rds")

# covariates 

X1 <- readRDS("RE_generation_DAGAR/X1.rds")
X2 <- readRDS("RE_generation_DAGAR/X2.rds")
X3 <- readRDS("RE_generation_DAGAR/X3.rds")
X4 <- readRDS("RE_generation_DAGAR/X4.rds")

Z1 <- readRDS("RE_generation_DAGAR/Z1.rds")

Y <- readRDS("RE_generation_DAGAR/Y_list.rds")

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

spec1 <- matrix(0, 100, length(T_edge))
sens1 <- matrix(0, 100, length(T_edge))
spec2 <- matrix(0, 100, length(T_edge))
sens2 <- matrix(0, 100, length(T_edge))
spec3 <- matrix(0, 100, length(T_edge))
sens3 <- matrix(0, 100, length(T_edge))
spec4 <- matrix(0, 100, length(T_edge))
sens4 <- matrix(0, 100, length(T_edge))

T_edge_ad <- c(40, 50, 60, 70, 80, 90, 100, 110, 120)

spec12 <- matrix(0, 100, length(T_edge_ad))
sens12 <- matrix(0, 100, length(T_edge_ad))
spec21 <- matrix(0, 100, length(T_edge_ad))
sens21 <- matrix(0, 100, length(T_edge_ad))
spec13 <- matrix(0, 100, length(T_edge_ad))
sens13 <- matrix(0, 100, length(T_edge_ad))
spec31 <- matrix(0, 100, length(T_edge_ad))
sens31 <- matrix(0, 100, length(T_edge_ad))
spec14 <- matrix(0, 100, length(T_edge_ad))
sens14 <- matrix(0, 100, length(T_edge_ad))
spec41 <- matrix(0, 100, length(T_edge_ad))
sens41 <- matrix(0, 100, length(T_edge_ad))
spec23 <- matrix(0, 100, length(T_edge_ad))
sens23 <- matrix(0, 100, length(T_edge_ad))
spec32 <- matrix(0, 100, length(T_edge_ad))
sens32 <- matrix(0, 100, length(T_edge_ad))
spec24 <- matrix(0, 100, length(T_edge_ad))
sens24 <- matrix(0, 100, length(T_edge_ad))
spec42 <- matrix(0, 100, length(T_edge_ad))
sens42 <- matrix(0, 100, length(T_edge_ad))
spec34 <- matrix(0, 100, length(T_edge_ad))
sens34 <- matrix(0, 100, length(T_edge_ad))
spec43 <- matrix(0, 100, length(T_edge_ad))
sens43 <- matrix(0, 100, length(T_edge_ad))

specW1 <- rep(0, 100)
sensW1 <- rep(0, 100)
specW2 <- rep(0, 100)
sensW2 <- rep(0, 100)
specW3 <- rep(0, 100)
sensW3 <- rep(0, 100)
specW4 <- rep(0, 100)
sensW4 <- rep(0, 100)

Z1 <- Z1[order(final_perm), order(final_perm)]

X1 <- X1[order(final_perm),]
X2 <- X2[order(final_perm),]
X3 <- X3[order(final_perm),]
X4 <- X4[order(final_perm),]

seed <- 1

summaryArray <- array(dim = c(28,4,100))

for (seed in 1:100) {
  print(seed)
  
  filename <- paste0("runs_DAGAR/mcmc_samples_", seed, ".rds")
  mcmc_samples <- readRDS(filename)
  
  names(mcmc_samples) <- c("beta", "taud", "phi", "theta", "u", "rho", "v", "r", 
                           "F_r", "eta", "taus", "W1", "W2", "W3", "W4", "A")
  
  summaryArray[,c(1,2),seed] <- cbind(c(colMeans(mcmc_samples$beta),
                                        colMeans(mcmc_samples$taud),
                                        colMeans(mcmc_samples$theta),
                                        mean(mcmc_samples$taus)), 
                                      c(sqrt(colVars(mcmc_samples$beta)),
                                        sqrt(colVars(mcmc_samples$taud)),
                                        sqrt(colVars(mcmc_samples$theta)),
                                        sd(mcmc_samples$taus)))
  
}

idx <- 1:100

idx <- idx[-c(11,54,56,57,66,87,90,97)]

for(seed in 1:100){
  
  #######true difference boundary########

  true_diff1 <- NULL
  true_diff2 <- NULL
  true_diff3 <- NULL
  true_diff4 <- NULL

  for(i in 1:nrow(neighbor_list0)){
    if(phi_true1[[seed]][neighbor_list0[i,1]] != phi_true1[[seed]][neighbor_list0[i,2]]){
      true_diff1 <- c(true_diff1, 1)
    }else{
      true_diff1 <- c(true_diff1, 0)
    }
    if(phi_true2[[seed]][neighbor_list0[i,1]] != phi_true2[[seed]][neighbor_list0[i,2]]){
      true_diff2 <- c(true_diff2, 1)
    }else{
      true_diff2 <- c(true_diff2, 0)
    }
    if(phi_true3[[seed]][neighbor_list0[i,1]] != phi_true3[[seed]][neighbor_list0[i,2]]){
      true_diff3 <- c(true_diff3, 1)
    }else{
      true_diff3 <- c(true_diff3, 0)
    }
    if(phi_true4[[seed]][neighbor_list0[i,1]] != phi_true4[[seed]][neighbor_list0[i,2]]){
      true_diff4 <- c(true_diff4, 1)
    }else{
      true_diff4 <- c(true_diff4, 0)
    }
  }

  # #######true difference boundaries across diseases#######
  # 
  # true_diff12 <- NULL
  # true_diff21 <- NULL
  # true_diff13 <- NULL
  # true_diff31 <- NULL
  # true_diff14 <- NULL
  # true_diff41 <- NULL
  # true_diff23 <- NULL
  # true_diff32 <- NULL
  # true_diff24 <- NULL
  # true_diff42 <- NULL
  # true_diff34 <- NULL
  # true_diff43 <- NULL
  # 
  # for(i in 1:nrow(neighbor_list0)){
  #   if(phi_true1[[seed]][neighbor_list0[i,1]] != phi_true2[[seed]][neighbor_list0[i,2]]){
  #     true_diff12 <- c(true_diff12, 1)
  #   }else{
  #     true_diff12 <- c(true_diff12, 0)
  #   }
  #   if(phi_true2[[seed]][neighbor_list0[i,1]] != phi_true1[[seed]][neighbor_list0[i,2]]){
  #     true_diff21 <- c(true_diff21, 1)
  #   }else{
  #     true_diff21 <- c(true_diff21, 0)
  #   }
  #   if(phi_true1[[seed]][neighbor_list0[i,1]] != phi_true3[[seed]][neighbor_list0[i,2]]){
  #     true_diff13 <- c(true_diff13, 1)
  #   }else{
  #     true_diff13 <- c(true_diff13, 0)
  #   }
  #   if(phi_true3[[seed]][neighbor_list0[i,1]] != phi_true1[[seed]][neighbor_list0[i,2]]){
  #     true_diff31 <- c(true_diff31, 1)
  #   }else{
  #     true_diff31 <- c(true_diff31, 0)
  #   }
  #   if(phi_true1[[seed]][neighbor_list0[i,1]] != phi_true4[[seed]][neighbor_list0[i,2]]){
  #     true_diff14 <- c(true_diff14, 1)
  #   }else{
  #     true_diff14 <- c(true_diff14, 0)
  #   }
  #   if(phi_true4[[seed]][neighbor_list0[i,1]] != phi_true1[[seed]][neighbor_list0[i,2]]){
  #     true_diff41 <- c(true_diff41, 1)
  #   }else{
  #     true_diff41 <- c(true_diff41, 0)
  #   }
  #   if(phi_true2[[seed]][neighbor_list0[i,1]] != phi_true3[[seed]][neighbor_list0[i,2]]){
  #     true_diff23 <- c(true_diff23, 1)
  #   }else{
  #     true_diff23 <- c(true_diff23, 0)
  #   }
  #   if(phi_true3[[seed]][neighbor_list0[i,1]] != phi_true2[[seed]][neighbor_list0[i,2]]){
  #     true_diff32 <- c(true_diff32, 1)
  #   }else{
  #     true_diff32 <- c(true_diff32, 0)
  #   }
  #   if(phi_true2[[seed]][neighbor_list0[i,1]] != phi_true4[[seed]][neighbor_list0[i,2]]){
  #     true_diff24 <- c(true_diff24, 1)
  #   }else{
  #     true_diff24 <- c(true_diff24, 0)
  #   }
  #   if(phi_true4[[seed]][neighbor_list0[i,1]] != phi_true2[[seed]][neighbor_list0[i,2]]){
  #     true_diff42 <- c(true_diff42, 1)
  #   }else{
  #     true_diff42 <- c(true_diff42, 0)
  #   }
  #   if(phi_true3[[seed]][neighbor_list0[i,1]] != phi_true4[[seed]][neighbor_list0[i,2]]){
  #     true_diff34 <- c(true_diff34, 1)
  #   }else{
  #     true_diff34 <- c(true_diff34, 0)
  #   }
  #   if(phi_true4[[seed]][neighbor_list0[i,1]] != phi_true3[[seed]][neighbor_list0[i,2]]){
  #     true_diff43 <- c(true_diff43, 1)
  #   }else{
  #     true_diff43 <- c(true_diff43, 0)
  #   }
  # }

  print(seed)
  
  filename <- paste0("runs_DAGAR/mcmc_samples_", seed, ".rds")
  mcmc_samples <- readRDS(filename)
  
  names(mcmc_samples) <- c("beta", "taud", "phi", "theta", "u", "rho", "v", "r", 
                           "F_r", "eta", "taus", "W1", "W2", "W3", "W4", "A")
  
  phis <-  mcmc_samples$phi
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
  
  # Estimated FDR curves
  
  threshold1 <- sort(pvij1, decreasing = TRUE)[T_edge]
  threshold2 <- sort(pvij2, decreasing = TRUE)[T_edge]
  threshold3 <- sort(pvij3, decreasing = TRUE)[T_edge]
  threshold4 <- sort(pvij4, decreasing = TRUE)[T_edge]
  
  pvijm <- rbind(pvij1, pvij2, pvij3, pvij4)
  
  # calculate sensitivity and specificity on boundary detection
  for(i in 1:length(T_edge)){
    est_diff1 <- factor(as.numeric(pvijm[1,] >= threshold1[i]), levels = c(0,1))
    conf_matrix1 <- table(est_diff1,true_diff1)
    if (ncol(conf_matrix1) == 1) {
      conf_matrix1 <- cbind(conf_matrix1,c(0,0))
      colnames(conf_matrix1)[2] <- "1"
      conf_matrix1 <- as.table(conf_matrix1)
    }
    spec1[seed,i] <- sensitivity(conf_matrix1)
    sens1[seed,i] <- specificity(conf_matrix1)
    
    est_diff2 <- factor(as.numeric(pvijm[2,] >= threshold2[i]), levels = c(0,1))
    conf_matrix2 <- table(est_diff2,true_diff2)
    if (ncol(conf_matrix2) == 1) {
      conf_matrix2 <- cbind(conf_matrix2,c(0,0))
      colnames(conf_matrix2)[2] <- "1"
      conf_matrix2 <- as.table(conf_matrix2)
    }
    spec2[seed,i] <- sensitivity(conf_matrix2)
    sens2[seed,i] <- specificity(conf_matrix2)
    
    est_diff3 <- factor(as.numeric(pvijm[3,] >= threshold3[i]), levels = c(0,1))
    conf_matrix3 <- table(est_diff3,true_diff3)
    if (ncol(conf_matrix3) == 1) {
      conf_matrix3 <- cbind(conf_matrix3,c(0,0))
      colnames(conf_matrix3)[2] <- "1"
      conf_matrix3 <- as.table(conf_matrix3)
    }
    spec3[seed,i] <- sensitivity(conf_matrix3)
    sens3[seed,i] <- specificity(conf_matrix3)
    
    est_diff4 <- factor(as.numeric(pvijm[4,] >= threshold4[i]), levels = c(0,1))
    conf_matrix4 <- table(est_diff4,true_diff4)
    if (ncol(conf_matrix4) == 1) {
      conf_matrix4 <- cbind(conf_matrix4,c(0,0))
      colnames(conf_matrix4)[2] <- "1"
      conf_matrix4 <- as.table(conf_matrix4)
    }
    spec4[seed,i] <- sensitivity(conf_matrix4)
    sens4[seed,i] <- specificity(conf_matrix4)
  }
  
  # # difference boundaries across cancers
  # 
  # vij_samples12 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]+58]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # vij_samples21 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # vij_samples13 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]+116]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # vij_samples31 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]+116] != x[neighbor_list0[i,2]]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # vij_samples14 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]] != x[neighbor_list0[i,2]+174]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # vij_samples41 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]+174] != x[neighbor_list0[i,2]]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # vij_samples23 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]+116]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # vij_samples32 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]+116] != x[neighbor_list0[i,2]+58]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # vij_samples24 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]+58] != x[neighbor_list0[i,2]+174]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # vij_samples42 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]+174] != x[neighbor_list0[i,2]+58]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # vij_samples34 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]+116] != x[neighbor_list0[i,2]+174]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # vij_samples43 <- t(apply(phis_origin, 1, function(x){
  #   #x <- phis[1,]
  #   diff_b1 <- NULL
  #   for(i in 1:nrow(neighbor_list0)){
  #     if(x[neighbor_list0[i,1]+174] != x[neighbor_list0[i,2]+116]){
  #       diff_b1 <- c(diff_b1, 1)
  #     }else{
  #       diff_b1 <- c(diff_b1, 0)
  #     }
  #   }
  #   return(diff_b1)
  # }))
  # 
  # pvij12 <- apply(vij_samples12, 2, mean)
  # pvij21 <- apply(vij_samples21, 2, mean)
  # pvij13 <- apply(vij_samples13, 2, mean)
  # pvij31 <- apply(vij_samples31, 2, mean)
  # pvij14 <- apply(vij_samples14, 2, mean)
  # pvij41 <- apply(vij_samples41, 2, mean)
  # pvij23 <- apply(vij_samples23, 2, mean)
  # pvij32 <- apply(vij_samples32, 2, mean)
  # pvij24 <- apply(vij_samples12, 2, mean)
  # pvij42 <- apply(vij_samples21, 2, mean)
  # pvij34 <- apply(vij_samples12, 2, mean)
  # pvij43 <- apply(vij_samples21, 2, mean)
  # pvijm <- rbind(pvij12, pvij21, pvij13, pvij31, pvij14, pvij41,
  #                pvij23, pvij32, pvij24, pvij42,
  #                pvij34, pvij43)
  # 
  # threshold12 <- sort(pvij12, decreasing = TRUE)[T_edge_ad]
  # threshold21 <- sort(pvij21, decreasing = TRUE)[T_edge_ad]
  # threshold13 <- sort(pvij13, decreasing = TRUE)[T_edge_ad]
  # threshold31 <- sort(pvij31, decreasing = TRUE)[T_edge_ad]
  # threshold14 <- sort(pvij14, decreasing = TRUE)[T_edge_ad]
  # threshold41 <- sort(pvij41, decreasing = TRUE)[T_edge_ad]
  # threshold23 <- sort(pvij23, decreasing = TRUE)[T_edge_ad]
  # threshold32 <- sort(pvij32, decreasing = TRUE)[T_edge_ad]
  # threshold24 <- sort(pvij24, decreasing = TRUE)[T_edge_ad]
  # threshold42 <- sort(pvij42, decreasing = TRUE)[T_edge_ad]
  # threshold34 <- sort(pvij34, decreasing = TRUE)[T_edge_ad]
  # threshold43 <- sort(pvij43, decreasing = TRUE)[T_edge_ad]
  # 
  # for(i in 1:length(T_edge_ad)){
  #   est_diff12 <- factor(as.numeric(pvijm[1,] >= threshold12[i]), levels = c(0,1))
  #   conf_matrix12 <- table(est_diff12,true_diff12)
  #   spec12[seed,i] <- sensitivity(conf_matrix12)
  #   sens12[seed,i] <- specificity(conf_matrix12)
  #   
  #   est_diff21 <- factor(as.numeric(pvijm[2,] >= threshold21[i]), levels = c(0,1))
  #   conf_matrix21 <- table(est_diff2,true_diff21)
  #   spec21[seed,i] <- sensitivity(conf_matrix21)
  #   sens21[seed,i] <- specificity(conf_matrix21)
  #   
  #   est_diff13 <- factor(as.numeric(pvijm[3,] >= threshold13[i]), levels = c(0,1))
  #   conf_matrix13 <- table(est_diff13,true_diff13)
  #   spec13[seed,i] <- sensitivity(conf_matrix13)
  #   sens13[seed,i] <- specificity(conf_matrix13)
  #   
  #   est_diff31 <- factor(as.numeric(pvijm[4,] >= threshold31[i]), levels = c(0,1))
  #   conf_matrix31 <- table(est_diff31,true_diff31)
  #   spec31[seed,i] <- sensitivity(conf_matrix31)
  #   sens31[seed,i] <- specificity(conf_matrix31)
  #   
  #   est_diff14 <- factor(as.numeric(pvijm[5,] >= threshold14[i]), levels = c(0,1))
  #   conf_matrix14 <- table(est_diff14,true_diff14)
  #   spec14[seed,i] <- sensitivity(conf_matrix14)
  #   sens14[seed,i] <- specificity(conf_matrix14)
  #   
  #   est_diff41 <- factor(as.numeric(pvijm[6,] >= threshold41[i]), levels = c(0,1))
  #   conf_matrix41 <- table(est_diff41,true_diff41)
  #   spec41[seed,i] <- sensitivity(conf_matrix41)
  #   sens41[seed,i] <- specificity(conf_matrix41)
  #   
  #   est_diff23 <- factor(as.numeric(pvijm[7,] >= threshold23[i]), levels = c(0,1))
  #   conf_matrix23 <- table(est_diff23,true_diff23)
  #   spec23[seed,i] <- sensitivity(conf_matrix23)
  #   sens23[seed,i] <- specificity(conf_matrix23)
  #   
  #   est_diff32 <- factor(as.numeric(pvijm[8,] >= threshold32[i]), levels = c(0,1))
  #   conf_matrix32 <- table(est_diff32,true_diff32)
  #   spec32[seed,i] <- sensitivity(conf_matrix32)
  #   sens32[seed,i] <- specificity(conf_matrix32)
  #   
  #   est_diff24 <- factor(as.numeric(pvijm[9,] >= threshold24[i]), levels = c(0,1))
  #   conf_matrix24 <- table(est_diff24,true_diff24)
  #   spec24[seed,i] <- sensitivity(conf_matrix24)
  #   sens24[seed,i] <- specificity(conf_matrix24)
  #   
  #   est_diff42 <- factor(as.numeric(pvijm[10,] >= threshold42[i]), levels = c(0,1))
  #   conf_matrix42 <- table(est_diff42,true_diff42)
  #   spec42[seed,i] <- sensitivity(conf_matrix42)
  #   sens42[seed,i] <- specificity(conf_matrix42)
  #   
  #   est_diff34 <- factor(as.numeric(pvijm[11,] >= threshold34[i]), levels = c(0,1))
  #   conf_matrix34 <- table(est_diff34,true_diff34)
  #   spec34[seed,i] <- sensitivity(conf_matrix34)
  #   sens34[seed,i] <- specificity(conf_matrix34)
  #   
  #   est_diff43 <- factor(as.numeric(pvijm[12,] >= threshold43[i]), levels = c(0,1))
  #   conf_matrix43 <- table(est_diff43,true_diff43)
  #   spec43[seed,i] <- sensitivity(conf_matrix43)
  #   sens43[seed,i] <- specificity(conf_matrix43)
  # }
  # 
  
}

# table for within disease boundary detection

spec1 <- spec1[-c(11,54,56,57,66,87,90,97),]
sens1 <- sens1[-c(11,54,56,57,66,87,90,97),]

spec2 <- spec2[-c(11,54,56,57,66,87,90,97),]
sens2 <- sens2[-c(11,54,56,57,66,87,90,97),]

spec3 <- spec3[-c(11,54,56,57,66,87,90,97),]
sens3 <- sens3[-c(11,54,56,57,66,87,90,97),]

spec4 <- spec4[-c(11,54,56,57,66,87,90,97),]
sens4 <- sens4[-c(11,54,56,57,66,87,90,97),]

spec1_mean_dagar <- colMeans(spec1)
sens1_mean_dagar <- colMeans(sens1)

spec2_mean_dagar <- colMeans(spec2)
sens2_mean_dagar <- colMeans(sens2)

spec3_mean_dagar <- colMeans(spec3)
sens3_mean_dagar <- colMeans(sens3)

spec4_mean_dagar <- colMeans(spec4)
sens4_mean_dagar <- colMeans(sens4)

table <- cbind(spec1_mean_dagar, sens1_mean_dagar, 
               spec2_mean_dagar, sens2_mean_dagar,
               spec3_mean_dagar, sens3_mean_dagar, 
               spec4_mean_dagar, sens4_mean_dagar)

rownames(table[9:13,]) <- c("85","90","95","100","105")

round(table[9:13,],3)

################################################################################

for (seed in 1:100) {
  
  #######true adjacencies#######

  true_adj1 <- NULL
  true_adj2 <- NULL
  true_adj3 <- NULL
  true_adj4 <- NULL

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
    if(W_true3[neighbor_list0[i,1],neighbor_list0[i,2]] == 1){
      true_adj3 <- c(true_adj3, 1)
    }else{
      true_adj3 <- c(true_adj3, 0)
    }
    if(W_true4[neighbor_list0[i,1],neighbor_list0[i,2]] == 1){
      true_adj4 <- c(true_adj4, 1)
    }else{
      true_adj4 <- c(true_adj4, 0)
    }
  }
  
  print(seed)
  
  filename <- paste0("runs_DAGAR/mcmc_samples_", seed, ".rds")
  mcmc_samples <- readRDS(filename)
  
  names(mcmc_samples) <- c("beta", "taud", "phi", "theta", "u", "rho", "v", "r", 
                           "F_r", "eta", "taus", "W1", "W2", "W3", "W4", "A")
  
  W1 <- mcmc_samples$W1#[order(final_perm),order(final_perm)]
  W2 <- mcmc_samples$W2#[order(final_perm),order(final_perm)]
  W3 <- mcmc_samples$W3#[order(final_perm),order(final_perm)]
  W4 <- mcmc_samples$W4#[order(final_perm),order(final_perm)]
  
  W1.mean <- apply(simplify2array(W1), c(1, 2), mean)
  W2.mean <- apply(simplify2array(W2), c(1, 2), mean)
  W3.mean <- apply(simplify2array(W3), c(1, 2), mean)
  W4.mean <- apply(simplify2array(W4), c(1, 2), mean)
  
  W1.mean <- W1.mean[order(final_perm), order(final_perm)]
  W2.mean <- W2.mean[order(final_perm), order(final_perm)]
  W3.mean <- W3.mean[order(final_perm), order(final_perm)]
  W4.mean <- W4.mean[order(final_perm), order(final_perm)]
  
  # calculate sensitivity and specificity on adjacency

  W <- replicate(4, Minc, simplify = FALSE)

  eta.mean <- apply(mcmc_samples$eta, 2, mean)

  Winc <- Minc[order(final_perm),order(final_perm)]
  Winc[upper.tri(Winc)] <- 0

  indices <- which(Winc == 1, arr.ind = TRUE)

  est_adj <- matrix(0,nrow = nrow(indices), ncol = 4)

  for (d in 1:4) {

    for (i in 1:nrow(indices)) {
      row_idx <- indices[i, "row"]
      col_idx <- indices[i, "col"]

      W[[d]][row_idx,col_idx] <- as.numeric(exp(-Z1[row_idx,col_idx]*eta.mean[d]) >= 0.5)

      est_adj[i,d] <- W[[d]][row_idx,col_idx]
    }

  }

  est_adj1 <- factor(est_adj[,1])
  conf_matrixW1 <- table(est_adj1,true_adj1)
  
  if (nrow(conf_matrixW1) == 1) {
    conf_matrixW1 <- rbind(c(0,0),conf_matrixW1)
    rownames(conf_matrixW1)[1] <- "0"
    conf_matrixW1 <- as.table(conf_matrixW1)
  } 

  specW1[seed] <- sensitivity(conf_matrixW1)
  sensW1[seed] <- specificity(conf_matrixW1)

  est_adj2 <- factor(est_adj[,2])
  conf_matrixW2 <- table(est_adj2,true_adj2)

  if (nrow(conf_matrixW2) == 1) {
    conf_matrixW2 <- rbind(c(0,0),conf_matrixW2)
    rownames(conf_matrixW2)[1] <- "0"
    conf_matrixW2 <- as.table(conf_matrixW2)
  } 

  specW2[seed] <- sensitivity(conf_matrixW2)
  sensW2[seed] <- specificity(conf_matrixW2)

  est_adj3 <- factor(est_adj[,3])
  conf_matrixW3 <- table(est_adj3,true_adj3)

  if (nrow(conf_matrixW3) == 1) {
    conf_matrixW3 <- rbind(c(0,0),conf_matrixW3)
    rownames(conf_matrixW3)[1] <- "0"
    conf_matrixW3 <- as.table(conf_matrixW3)
  } 

  specW3[seed] <- sensitivity(conf_matrixW3)
  sensW3[seed] <- specificity(conf_matrixW3)

  est_adj4 <- factor(est_adj[,4])
  conf_matrixW4 <- table(est_adj4,true_adj4)
  
  if (nrow(conf_matrixW4) == 1) {
    conf_matrixW4 <- rbind(c(0,0),conf_matrixW4)
    rownames(conf_matrixW4)[1] <- "0"
    conf_matrixW4 <- as.table(conf_matrixW4)
  }

  specW4[seed] <- sensitivity(conf_matrixW4)
  sensW4[seed] <- specificity(conf_matrixW4)
  
}




#table for between diseases boundary detection

spec12_mean_dagar <- colMeans(spec12)
sens12_mean_dagar <- colMeans(sens12)

spec21_mean_dagar <- colMeans(spec21)
sens21_mean_dagar <- colMeans(sens21)

spec13_mean_dagar <- colMeans(spec13)
sens13_mean_dagar <- colMeans(sens13)

spec31_mean_dagar <- colMeans(spec31)
sens31_mean_dagar <- colMeans(sens31)

spec14_mean_dagar <- colMeans(spec14)
sens14_mean_dagar <- colMeans(sens14)

spec41_mean_dagar <- colMeans(spec41)
sens41_mean_dagar <- colMeans(sens41)

spec23_mean_dagar <- colMeans(spec23)
sens23_mean_dagar <- colMeans(sens23)

spec32_mean_dagar <- colMeans(spec32)
sens32_mean_dagar <- colMeans(sens32)

spec24_mean_dagar <- colMeans(spec24)
sens24_mean_dagar <- colMeans(sens24)

spec42_mean_dagar <- colMeans(spec42)
sens42_mean_dagar <- colMeans(sens42)

spec34_mean_dagar <- colMeans(spec34)
sens34_mean_dagar <- colMeans(sens34)

spec43_mean_dagar <- colMeans(spec43)
sens43_mean_dagar <- colMeans(sens43)

table_ad <- cbind(spec12_mean_dagar, sens12_mean_dagar, 
                  spec21_mean_dagar, sens21_mean_dagar,
                  spec13_mean_dagar, sens13_mean_dagar, 
                  spec31_mean_dagar, sens31_mean_dagar,
                  spec14_mean_dagar, sens14_mean_dagar, 
                  spec41_mean_dagar, sens41_mean_dagar,
                  spec23_mean_dagar, sens23_mean_dagar, 
                  spec32_mean_dagar, sens32_mean_dagar,
                  spec24_mean_dagar, sens24_mean_dagar, 
                  spec42_mean_dagar, sens42_mean_dagar,
                  spec34_mean_dagar, sens34_mean_dagar, 
                  spec43_mean_dagar, sens43_mean_dagar)

table_ad

# adjacency detection tables

# specW1 <- specW1[-c(11,54,56,57,66,87,90,97)]
# sensW1 <- sensW1[-c(11,54,56,57,66,87,90,97)]
# 
# specW2 <- specW2[-c(11,54,56,57,66,87,90,97)]
# sensW2 <- sensW2[-c(11,54,56,57,66,87,90,97)]
# 
# specW3 <- specW3[-c(11,54,56,57,66,87,90,97)]
# sensW3 <- sensW3[-c(11,54,56,57,66,87,90,97)]
# 
# specW4 <- specW4[-c(11,54,56,57,66,87,90,97)]
# sensW4 <- sensW4[-c(11,54,56,57,66,87,90,97)]

specW1_mean_dagar <- mean(specW1)
sensW1_mean_dagar <- mean(sensW1)

specW2_mean_dagar <- mean(specW2)
sensW2_mean_dagar <- mean(sensW2)

specW3_mean_dagar <- mean(specW3)
sensW3_mean_dagar <- mean(sensW3)

specW4_mean_dagar <- mean(specW4)
sensW4_mean_dagar <- mean(sensW4)

table_adj <- cbind(specW1_mean_dagar, sensW1_mean_dagar, 
                   specW2_mean_dagar, sensW2_mean_dagar,
                   specW3_mean_dagar, sensW3_mean_dagar, 
                   specW4_mean_dagar, sensW4_mean_dagar)

round(table_adj,3)

#WAIC

WAIC <- rep(0,100)

for(seed in 1:100){
  
  filename <- paste0("runs_DAGAR/mcmc_samples_", seed, ".rds")
  mcmc_samples <- readRDS(filename)
  
  names(mcmc_samples) <- c("beta", "taud", "phi", "theta", "u", "rho", "v", "r", 
                           "F_r", "eta", "taus", "W1", "W2", "W3", "W4", "A")
  
  print(seed)
  set.seed(seed)
  
  beta1 <- c(2,1)
  beta2 <- c(1,2)
  beta3 <- c(2,2)
  beta4 <- c(1,2)
  
  taud1 <- 10
  taud2 <- 10
  taud3 <- 10
  taud4 <- 10
  
  y1 <- X1 %*% beta1 + (phi_true1[[seed]] - mean(phi_true1[[seed]]))  + sqrt(1/taud1) * rnorm(n_county)
  y2 <- X2 %*% beta2 + (phi_true2[[seed]] - mean(phi_true2[[seed]]))  + sqrt(1/taud2) * rnorm(n_county)
  y3 <- X3 %*% beta3 + (phi_true3[[seed]] - mean(phi_true3[[seed]]))  + sqrt(1/taud3) * rnorm(n_county)
  y4 <- X4 %*% beta4 + (phi_true4[[seed]] - mean(phi_true4[[seed]]))  + sqrt(1/taud4) * rnorm(n_county)
  
  X <- as.matrix(bdiag(bdiag(X1[final_perm,], X2[final_perm,]),
                       bdiag(X3[final_perm,], X4[final_perm,])))
  
  yo1 <- y1[final_perm]
  yo2 <- y2[final_perm]
  yo3 <- y3[final_perm]
  yo4 <- y4[final_perm]

  Y <- c(yo1, yo2, yo3, yo4)
  
  LL <- matrix(0, nrow = 4*n, ncol = 10000)
  
  for(i in 1:ncol(LL)){
    
    beta.post <- mcmc_samples$beta
    phi.post <- mcmc_samples$phi
    
    LL[1:n,i] <- dnorm(Y[1:n], X[1:n,1:2]%*%beta.post[i,1:2] + phi.post[i,1:n], sqrt(1/taud1), log = T)
    
    LL[(n+1):(2*n),i] <- dnorm(Y[(n+1):(2*n)], X[(n+1):(2*n),3:4]%*%beta.post[i,3:4] + phi.post[i,(n+1):(2*n)], sqrt(1/taud2), log = T)
    
    LL[(2*n+1):(3*n),i] <- dnorm(Y[(2*n+1):(3*n)], X[(2*n+1):(3*n),5:6]%*%beta.post[i,5:6] + phi.post[i,(2*n+1):(3*n)], sqrt(1/taud3), log = T)
    
    LL[(3*n+1):(4*n),i] <- dnorm(Y[(3*n+1):(4*n)], X[(3*n+1):(4*n),7:8]%*%beta.post[i,7:8] + phi.post[i,(3*n+1):(4*n)], sqrt(1/taud4), log = T)
    
  }
  
  WAIC[seed] <- WAIC(LL)$WAIC
  
}

saveRDS(WAIC, file = "WAIC_unstructured.rds")

#WAIC loo package

library(loo)

WAIC <- rep(0,100)

for(seed in 1:100){
  
  filename <- paste0("runs_DAGAR/mcmc_samples_", seed, ".rds")
  mcmc_samples <- readRDS(filename)
  
  names(mcmc_samples) <- c("beta", "taud", "phi", "theta", "u", "rho", "v", "r", 
                           "F_r", "eta", "taus", "W1", "W2", "W3", "W4", "A")
  
  print(seed)
  set.seed(seed)
  
  beta1 <- c(2,1)
  beta2 <- c(1,2)
  beta3 <- c(2,2)
  beta4 <- c(1,2)
  
  taud1 <- 10
  taud2 <- 10
  taud3 <- 10
  taud4 <- 10
  
  y1 <- X1 %*% beta1 + (phi_true1[[seed]] - mean(phi_true1[[seed]]))  + sqrt(1/taud1) * rnorm(n_county)
  y2 <- X2 %*% beta2 + (phi_true2[[seed]] - mean(phi_true2[[seed]]))  + sqrt(1/taud2) * rnorm(n_county)
  y3 <- X3 %*% beta3 + (phi_true3[[seed]] - mean(phi_true3[[seed]]))  + sqrt(1/taud3) * rnorm(n_county)
  y4 <- X4 %*% beta4 + (phi_true4[[seed]] - mean(phi_true4[[seed]]))  + sqrt(1/taud4) * rnorm(n_county)
  
  X <- as.matrix(bdiag(bdiag(X1[final_perm,], X2[final_perm,]),
                       bdiag(X3[final_perm,], X4[final_perm,])))
  
  yo1 <- y1[final_perm]
  yo2 <- y2[final_perm]
  yo3 <- y3[final_perm]
  yo4 <- y4[final_perm]
  
  Y <- c(yo1, yo2, yo3, yo4)
  
  LL <- matrix(0, nrow = 4*n, ncol = 10000)
  
  for(i in 1:ncol(LL)){
    
    beta.post <- mcmc_samples$beta
    phi.post <- mcmc_samples$phi
    
    LL[1:n,i] <- dnorm(Y[1:n], X[1:n,1:2]%*%beta.post[i,1:2] + phi.post[i,1:n], sqrt(1/taud1), log = T)
    
    LL[(n+1):(2*n),i] <- dnorm(Y[(n+1):(2*n)], X[(n+1):(2*n),3:4]%*%beta.post[i,3:4] + phi.post[i,(n+1):(2*n)], sqrt(1/taud2), log = T)
    
    LL[(2*n+1):(3*n),i] <- dnorm(Y[(2*n+1):(3*n)], X[(2*n+1):(3*n),5:6]%*%beta.post[i,5:6] + phi.post[i,(2*n+1):(3*n)], sqrt(1/taud3), log = T)
    
    LL[(3*n+1):(4*n),i] <- dnorm(Y[(3*n+1):(4*n)], X[(3*n+1):(4*n),7:8]%*%beta.post[i,7:8] + phi.post[i,(3*n+1):(4*n)], sqrt(1/taud4), log = T)
    
  }
  
  WAIC[seed] <- waic(t(LL))[[1]][3,1]
  
}

saveRDS(WAIC, file = "WAIC_unstructured.rds")

################################################################################

