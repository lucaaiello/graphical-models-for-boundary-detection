rm(list = ls())
setwd("C:/Dati/Dottorato/Visiting UCLA/Spatial Disease Mapping/Bibliografia/gao_2022/simulation/Disease graph/directed")

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

W_dis <- cbind(c(0,1,0,1),c(0,0,1,0),c(0,0,0,1),c(0,0,0,0))

saveRDS(W_dis, "W_dis_try.rds")

alpha0 <- cbind(c(0,0.3,0,0.5),c(0,0,0.4,0),c(0,0,0,0.8),c(0,0,0,0))
  
alpha1 <- cbind(c(0,0.5,0,0.4),c(0,0,0.4,0),c(0,0,0,0.1),c(0,0,0,0))
  
### Specify Covariance matrix for spatial components

q <- 4

rho <- c(0.2, 0.8, 0.4, 0.6)

A <- matrix(0,n*q,n*q)

indices <- which(W_dis == 1, arr.ind = TRUE)

# Iterate over the indices
for (i in 1:nrow(indices)) {
  
  row_idx <- indices[i, "row"]
  col_idx <- indices[i, "col"]
  
  A[((row_idx-1)*n + 1):(row_idx*n),((col_idx-1)*n + 1):(col_idx*n)] <- 
    alpha0[row_idx,col_idx]*diag(n) + alpha1[row_idx,col_idx]*Minc
  
}

alpha <- 1
sigmas_sq <- 4

#

K <- 15
phi_true1 <- list()
phi_true2 <- list()
nq <- n*q

### Generate discrete random effects
seed <- 16
set.seed(seed)

z <- rnorm(n,mean = 15, sd = 5)

Z <- sd_diff_mat(z,Minc)

# Z <- matrix(NA, nrow = n, ncol = n)
# 
# for (i in 2:n) {
#   for (j in 1:(i-1)) {
# 
#     Z[i, j] <- abs(rnorm(1,3,1))
#     Z[j, i] <- Z[i,j]
#      
#   }
# }

M <- as.numeric(-log(0.5)/quantile(Z[which(Z!=0)],0.5))

om <- c(0.5,0.25,0.33,0.6)

W <- replicate(q, Minc, simplify = FALSE)

indices <- which(Minc == 1, arr.ind = TRUE)

for (d in 1:q) {
  
  for (i in 1:nrow(indices)) {
    row_idx <- indices[i, "row"]
    col_idx <- indices[i, "col"]
    
    W[[d]][row_idx,col_idx] <- as.numeric(exp(-Z[row_idx,col_idx]*om[d]) >= 0.5)
  }
  
}

Q <- Dinv_new(rho, n, q, W)
invQ1 <- solve(Q[[1]])
invQ2 <- solve(Q[[2]])
invQ3 <- solve(Q[[3]])
invQ4 <- solve(Q[[4]])

invQ <- as.matrix(bdiag(bdiag(invQ1, invQ2), bdiag(invQ3, invQ4)))
Vr <- as.matrix(forceSymmetric(solve(diag(nq)-A) %*% invQ %*% solve(diag(nq)-t(A))))

saveRDS(Vr, file = "Vr_true_try.rds")

r <- rmvnorm(1, rep(0, nq), Vr)
v <- rbeta(K,1,alpha)

probvec <- makeprobs(v)
F_r <- pnorm(r,0,sqrt(diag(Vr)))
u <- makeu(F_r,probvec)

thetavec <- rnorm(K, 0, sqrt(sigmas_sq))
phivec <- thetavec[u]

W_true1 <- W[[1]][order(final_perm),order(final_perm)]   
W_true2 <- W[[2]][order(final_perm),order(final_perm)]
W_true3 <- W[[3]][order(final_perm),order(final_perm)]   
W_true4 <- W[[4]][order(final_perm),order(final_perm)]

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

true_adj1

true_adj2

true_adj3

true_adj4

saveRDS(W_true1, "W_true1_try.rds")
saveRDS(W_true2, "W_true2_try.rds")
saveRDS(W_true3, "W_true3_try.rds")
saveRDS(W_true4, "W_true4_try.rds")

phi_true1 <- phivec[1:n][order(final_perm)]
phi_true2 <- phivec[(n+1):(2*n)][order(final_perm)]
phi_true3 <- phivec[(2*n+1):(3*n)][order(final_perm)]
phi_true4 <- phivec[(3*n+1):(4*n)][order(final_perm)]

saveRDS(phi_true1, "phi_true1_try.rds")
saveRDS(phi_true2, "phi_true2_try.rds")
saveRDS(phi_true3, "phi_true3_try.rds")
saveRDS(phi_true4, "phi_true4_try.rds")

Z <- Z[order(final_perm),order(final_perm)]

saveRDS(Z, "Z_try.rds")

ca.poly$phi_true1 <- phi_true1
ca.poly$phi_true2 <- phi_true2
ca.poly$phi_true3 <- phi_true3
ca.poly$phi_true4 <- phi_true4

values <- sort(unique(c(phi_true1,phi_true2,phi_true3,phi_true4)))
color.pallete <- brewer.pal(length(values),"Blues")

col1 <- rep(0, length(phi_true1))
for(i in 1:length(values)){
  col1[phi_true1 == values[i]] <- color.pallete[i]
}

col2 <- rep(0, length(phi_true2))
for(i in 1:length(values)){
  col2[phi_true2 == values[i]] <- color.pallete[i]
}

col3 <- rep(0, length(phi_true3))
for(i in 1:length(values)){
  col3[phi_true3 == values[i]] <- color.pallete[i]
}

col4 <- rep(0, length(phi_true4))
for(i in 1:length(values)){
  col4[phi_true4 == values[i]] <- color.pallete[i]
}

################################################################################

# Plot boundaries
if (!require(gpclib)) install.packages("gpclib", type="source")
#gpclibPermit()
# p <- list()
# p_est <- list()

library(plyr)
library(rgeos)

rownames(ca.poly@data) <- colnames(Minc[order(final_perm),order(final_perm)])
ca.poly@data$id <- rownames(ca.poly@data)
ca.poly.f <- fortify(ca.poly, region = "id")
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

ca.poly.df_aux <- rbind(ca.poly.df[,1:7],ca.poly.df[,1:7],
                        ca.poly.df[,1:7],ca.poly.df[,1:7])
ca.poly.df_aux <- cbind(ca.poly.df_aux,
                        c(ca.poly.df$phi_true1,ca.poly.df$phi_true2,
                          ca.poly.df$phi_true2,ca.poly.df$phi_true4)) 
ca.poly.df_aux <- cbind(ca.poly.df_aux,c(rep("Cancer 1",2977),rep("Cancer 2",2977),
                                         rep("Cancer 3",2977),rep("Cancer 4",2977)))

colnames(ca.poly.df_aux)[c(8,9)] <- c("phi","cancer")

ggplot() +
  geom_polygon(data = ca.poly.df_aux,  aes(long, lat, group = group, fill = factor(phi)),
               color = "black") +
  facet_wrap(~ cancer) + 
  scale_fill_manual(values = color.pallete,
                    labels = c("-5.51", "-1.50", "-0.77", "-0.60", "0.87", "2.03")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=16, face="bold",hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.key.size = unit(0.75,"cm"),
        legend.position = "top")

count_1 <- 0
count_2 <- 0
count_3 <- 0
count_4 <- 0

for (i in 1:dim(neighbor_list0)[1]) {
  if (phi_true1[neighbor_list0[i,1]]!=phi_true1[neighbor_list0[i,2]]){
    count_1 <- count_1 + 1
  }
  if (phi_true2[neighbor_list0[i,1]]!=phi_true2[neighbor_list0[i,2]]){
    count_2 <- count_2 + 1
  }
  if (phi_true3[neighbor_list0[i,1]]!=phi_true3[neighbor_list0[i,2]]){
    count_3 <- count_3 + 1
  }
  if (phi_true4[neighbor_list0[i,1]]!=phi_true4[neighbor_list0[i,2]]){
    count_4 <- count_4 + 1
  }
}

Winc_aux <- Minc[order(final_perm),order(final_perm)]

Winc_aux[upper.tri(Winc_aux)] <- 0

indices <- which(Winc_aux == 1, arr.ind = TRUE)

count_12 <- 0
count_21 <- 0

for (i in 1:nrow(indices)) {
  
  row_idx <- indices[i, "row"]
  col_idx <- indices[i, "col"]
  
  if (phi_true1[col_idx]!=phi_true2[row_idx]){
    count_12 <- count_12 + 1
  }
  if (phi_true2[col_idx]!=phi_true1[row_idx]){
    count_21 <- count_21 + 1
  }
  
}

################################################################################

