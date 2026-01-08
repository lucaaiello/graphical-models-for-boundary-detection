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

phi_true1 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/phi_true1.rds")
phi_true2 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/phi_true2.rds")
phi_true3 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/phi_true3.rds")
phi_true4 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/phi_true4.rds")

# spatial adjacency matrix for diseases

W_true1 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/W_true1.rds")
W_true2 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/W_true2.rds")
W_true3 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/W_true3.rds")
W_true4 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/W_true4.rds")

# covariates 

X1 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/X1.rds")
X2 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/X2.rds")
X3 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/X3.rds")
X4 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/X4.rds")

Z1 <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/Z1.rds")

Y <- readRDS("Multivariate_Li_comparison/Our_method/RE_generation_DAGAR/Y_list.rds")

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

  print(seed)
  
  filename <- paste0("Multivariate_Li_comparison/Our_method/runs_DAGAR/mcmc_samples_", seed, ".rds")
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
  
}

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

colnames(table) <- c("spec1", "sens1", "spec2", "sens2", "spec3", "sens3", "spec4", "sens4")
rownames(table) <- as.character(c(45, 50 ,55 ,60, 65, 70, 75, 80, 85, 90, 95, 100, 105))

round(table,3)

################################################################################

A_mean_replicates <- array(dim = c(4 , 4, 100))

for(seed in 1:100){
  
  print(seed)
  
  filename <- paste0("Multivariate_Li_comparison/Our_method/runs_DAGAR/mcmc_samples_", seed, ".rds")
  mcmc_samples <- readRDS(filename)
  
  names(mcmc_samples) <- c("beta", "taud", "phi", "theta", "u", "rho", "v", "r", 
                           "F_r", "eta", "taus", "W1", "W2", "W3", "W4", "A")
  
  A_mean_replicates[,,seed] <- Reduce("+", mcmc_samples$A) / length(mcmc_samples$A)
  
}

# Compute mean across replicates
A_mean <- apply(A_mean_replicates, c(1,2), mean)

# Compute ratio: |A_ij| / median(diag)
diag_median <- median(diag(A_mean))
A_ratio <- abs(A_mean) / diag_median

# Optional: set diagonal to NA for labeling but keep numeric value for scale
diag(A_ratio) <- NA  # this will not affect scale if we set max properly

# Convert to long format
A_long <- melt(A_ratio)
colnames(A_long) <- c("Row", "Col", "Ratio")

# If all off-diagonals are NA, set a small max to avoid -Inf
max_ratio <- max(A_long$Ratio, na.rm = TRUE)
if(is.infinite(max_ratio)) max_ratio <- 1e-3  # small number for color scale

# Plot
ggplot(A_long, aes(x = factor(Col), y = factor(Row), fill = Ratio)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(Ratio), "", round(Ratio,3)))) +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, max_ratio)) +
  labs(x = "Column", y = "Row", fill = "Ratio (|offdiag| / median(diag))",
       title = "Relative magnitude of off-diagonal elements of A") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

