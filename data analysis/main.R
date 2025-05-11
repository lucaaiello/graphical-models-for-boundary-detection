rm(list = ls())

setwd("")

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

# Neighboring counties information 

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

cvrts = "adj" # "adj" for covariates only in the adjacency
              # "mean" for covariates only in the mean structure
              # "meanadj" for covariates both in the mean and in the adjacency model                   
            
if (cvrts == "mean" | cvrts = "meanadj"){
  X <- as.matrix(bdiag(bdiag(X1[,c(1,2,4,6)], X2[,c(1,2,4,6)]),
                       bdiag(X3[,c(1,2,4,6)], X4[,c(1,2,4,6)])))  
} else if (cvrts == "adj") {
  X <- as.matrix(bdiag(bdiag(X1[,1], X2[,1]),
                       bdiag(X3[,1], X4[,1])))  
}

library(Rcpp)
library(RcppArmadillo)

sourceCpp('sampler.cpp')

set.seed(12345)

library(tictoc)

tic()

mcmc_samples <- MADAGAR(y=Y, X=X, Z1=Z1, Z2=Z2, Z3=Z3, E=E, cvrts = "adj",
                        q=4, Winc=Winc, Minc=Minc,
                        alpha=1, n_atoms=15,
                        runs=10000, burn=50000, thin=5)

toc()

names(mcmc_samples) <- c("beta", "phi", "theta", "u", "rho", "V", "r", 
                         "F_r", "eta", "tau", "W1", "W2", "W3", "W4", "A")

samples.mcmc <- mcmc.list(mcmc(data.frame(beta = mcmc_samples$beta,
                                          theta = mcmc_samples$theta,
                                          gamma = mcmc_samples$r,
                                          V = mcmc_samples$V,
                                          tau = mcmc_samples$tau,
                                          rho = mcmc_samples$rho,
                                          eta = mcmc_samples$eta,
                                          A = mcmc_samples$A)))

################################################################################

cbind(round(apply(as.matrix(samples.mcmc[[1]][,c(1:19,267)]),2,mean),3),
      round(apply(as.matrix(samples.mcmc[[1]][,c(1:19,267)]),2,sd),3))

# cbind(round(apply(as.matrix(samples.mcmc[[1]][,c(1:31,279)]),2,mean),3),
#       round(apply(as.matrix(samples.mcmc[[1]][,c(1:31,279)]),2,sd),3))

library(mcmcse)

# Batch means estimator
mcerror_bm <- mcse.multi(x = as.data.frame(samples.mcmc[[1]][,-which(colVars(samples.mcmc[[1]])==0)]), 
                         method = "bm", r = 1,
                         size = NULL, g = NULL, adjust = TRUE,
                         blather = TRUE)

which(colVars(samples.mcmc[[1]])==0)

round(mcse.mat(x = as.data.frame(samples.mcmc[[1]][,c(1:19,267)]), method = "bm", g = NULL),3)
# round(mcse.mat(x = as.data.frame(samples.mcmc[[1]][,c(1:31,279)]), method = "bm", g = NULL),3)

ess(as.data.frame(samples.mcmc[[1]]))

plot(confRegion(mcerror_bm, which = c(2,3), level = .95), type = 'l', asp = 1)
lines(confRegion(mcerror_bm, which = c(2,3), level = .90), col = "red")

multiESS(as.data.frame(samples.mcmc[[1]][,-which(colVars(samples.mcmc[[1]])==0)]), covmat = mcerror_bm$cov)
multiESS(as.data.frame(samples.mcmc[[1]]), covmat = mcerror_bm$cov)

minESS(p = 293, alpha = .05, eps = .05) 
minESS(p = 293, alpha = .05, ess = 1847.068)


# Chain diagnostics ------------------------------------------------------------

samples.ggs <- ggs(samples.mcmc, keep_original_order = TRUE)

greek_label <- function(variable, base_greek) {
  sapply(variable, function(v) {
    num <- sub(paste0(base_greek, "\\."), "", v)
    bquote(.(as.name(base_greek))[.(num)])
  })
}

ggs_caterpillar(samples.ggs, family = "beta") + 
  scale_y_discrete(labels = function(labels) greek_label(labels, "beta")) + 
  theme_bw() +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) 

ggs_caterpillar(samples.ggs, family = "theta") +
  scale_y_discrete(labels = function(labels) greek_label(labels, "theta")) + 
  theme_bw() +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) 

samples.gamma1 <- samples.ggs[samples.ggs$Parameter %in% c(paste0("gamma.", 1:58)), ]  
samples.gamma2 <- samples.ggs[samples.ggs$Parameter %in% c(paste0("gamma.", 59:116)), ]
samples.gamma3 <- samples.ggs[samples.ggs$Parameter %in% c(paste0("gamma.", 117:174)), ]
samples.gamma4 <- samples.ggs[samples.ggs$Parameter %in% c(paste0("gamma.", 175:232)), ]

ggs_caterpillar(samples.gamma1) +
  scale_y_discrete(labels = function(labels) greek_label(labels, "gamma")) +
  theme_bw() +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) 

ggs_caterpillar(samples.gamma2) +
  scale_y_discrete(labels = function(labels) greek_label(labels, "gamma")) +
  theme_bw() +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) 

ggs_caterpillar(samples.gamma3) +
  scale_y_discrete(labels = function(labels) greek_label(labels, "gamma")) +
  theme_bw() +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) 

ggs_caterpillar(samples.gamma4) +
  scale_y_discrete(labels = function(labels) greek_label(labels, "gamma")) + 
  theme_bw() +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) 

ggs_caterpillar(samples.ggs, family = "V") + 
  scale_y_discrete(labels = function(labels) greek_label(labels, "V")) +
  theme_bw() +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) 

ggs_caterpillar(samples.ggs, family = "rho") + 
  scale_y_discrete(labels = function(labels) greek_label(labels, "rho")) +
  theme_bw() +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) 

samples.eta <- samples.ggs[samples.ggs$Parameter %in% c(paste0("eta.", 1:12)), ]  

ggs_caterpillar(samples.eta) + 
  scale_y_discrete(labels = function(labels) greek_label(labels, "eta")) + 
  theme_bw() +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) 

ggs_caterpillar(samples.ggs, family = "A") + 
  scale_y_discrete(labels = function(labels) greek_label(labels, "A")) +
  theme_bw() +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) 

# A parameters

AL <- array(0,dim = c(4,4,10000))

AAT <- array(0,dim = c(4,4,10000))

A_aux <- matrix(0,nrow = 10000, ncol = 10)

for (g in 1:dim(AL)[3]) {
  print(g)
  AL[,,g][lower.tri(AL[,,g], diag = TRUE)] <- samples.mcmc[[1]][g,284:293]
  # AL[,,g][lower.tri(AL[,,g], diag = TRUE)] <- samples.mcmc[[1]][g,296:305]
  AAT[,,g] <- AL[,,g]%*%t(AL[,,g])
  A_aux[g,] <- AAT[,,g][lower.tri(AL[,,g], diag = TRUE)]                        
}

A_aux.mcmc <- mcmc.list(mcmc(data.frame(AAT = A_aux)))

A_aux.ggs <- ggs(A_aux.mcmc, keep_original_order = TRUE)

ggs_caterpillar(A_aux.ggs) + 
  scale_y_discrete(labels = c(expression(AA[41]^T), expression(AA[31]^T), expression(AA[21]^T),
                              expression(AA[32]^T), expression(AA[42]^T), expression(AA[22]^T),
                              expression(AA[43]^T), expression(AA[33]^T), expression(AA[44]^T),
                              expression(AA[11]^T))) +
  # scale_y_discrete(labels = c(expression(AA[21]^T), expression(AA[31]^T), expression(AA[41]^T),
  #                             expression(AA[11]^T), expression(AA[42]^T), expression(AA[32]^T),
  #                             expression(AA[22]^T), expression(AA[43]^T), expression(AA[33]^T),
  #                             expression(AA[44]^T))) +
  # scale_y_discrete(labels = c(expression(AA[21]^T), expression(AA[42]^T), expression(AA[43]^T),
  #                             expression(AA[22]^T), expression(AA[32]^T), expression(AA[41]^T),
  #                             expression(AA[31]^T), expression(AA[44]^T), expression(AA[33]^T),
  #                             expression(AA[11]^T))) +
  theme_bw() +
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15)) 

# extracting spatial random effects --------------------------------------------

phis <- mcmc_samples$phi
phis1 <- phis[,1:58][,order(final_perm)]
phis2 <- phis[,59:116][,order(final_perm)]
phis3 <- phis[,117:174][,order(final_perm)]
phis4 <- phis[,175:232][,order(final_perm)]
phis_origin <- cbind(phis1, phis2, phis3, phis4)

# Set options to use tigris package
options(tigris_use_cache = TRUE)

# Download the counties shapefile for California
ca_counties <- counties(state = "CA", cb = TRUE)

# Convert the counties data to an sf object
ca_counties_sf <- st_as_sf(ca_counties)
ca_counties_sf <- ca_counties_sf[order(ca_counties_sf$NAME),]

color <- brewer.pal(9,"YlOrRd")

lung_SIR <- ggplot(data = ca_counties_sf) +
  geom_sf(aes(fill = rate_lung1$standard_ratio), color = "black", alpha = 0.6) +
  scale_fill_gradientn(colours = color) +
  ggtitle("Lung") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, face="bold",hjust = 0.5),
        legend.text = element_text(size = 30),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank())

esophageal_SIR <- ggplot(data = ca_counties_sf) +
  geom_sf(aes(fill = rate_esophagus1$standard_ratio), color = "black", alpha = 0.6) +
  scale_fill_gradientn(colours = color) +
  ggtitle("Esophageal") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, face="bold",hjust = 0.5),
        legend.text = element_text(size = 30),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank())

larynx_SIR <- ggplot(data = ca_counties_sf) +
  geom_sf(aes(fill = rate_larynx1$standard_ratio), color = "black", alpha = 0.6) +
  scale_fill_gradientn(colours = color) +
  ggtitle("Larynx") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, face="bold",hjust = 0.5),
        legend.text = element_text(size = 30),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank())

colorectal_SIR <- ggplot(data = ca_counties_sf) +
  geom_sf(aes(fill = rate_colrect1$standard_ratio), color = "black", alpha = 0.6) +
  scale_fill_gradientn(colours = color) +
  ggtitle("Colorectal") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, face="bold",hjust = 0.5),
        legend.text = element_text(size = 30),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank())

ggarrange(lung_SIR, esophageal_SIR, larynx_SIR, colorectal_SIR, nrow = 2, ncol = 2)

################################################################################

smoking_plot <- ggplot(ca_counties_sf) +
  geom_sf(aes(fill = rate_lung1$smoking), color = "black", alpha = 0.6) +
  scale_fill_gradientn(colours = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, face="bold",hjust = 0.5),
        legend.text = element_text(size = 30),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank()) +
  ggtitle("Smoking rate")

old_plot <- ggplot(ca_counties_sf) +
  geom_sf(aes(fill = rate_lung1$old), color = "black", alpha = 0.6) +
  scale_fill_gradientn(colours = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, face="bold",hjust = 0.5),
        legend.text = element_text(size = 30),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank()) +
  ggtitle("Over 65 rate")

poverty_plot <- ggplot(ca_counties_sf) +
  geom_sf(aes(fill = rate_lung1$poverty), color = "black", alpha = 0.6) +
  scale_fill_gradientn(colours = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, face="bold",hjust = 0.5),
        legend.text = element_text(size = 30),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank()) +
  ggtitle("Below poverty rate")

ggarrange(smoking_plot, old_plot, poverty_plot, nrow = 1, ncol = 3)

################################################################################

library(tidyverse)
library(cowplot)

library(rmapshaper)
bord <- ca_counties_sf %>% ms_innerlines()

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

for(i in 1:length(threshold1)){
  th1 <- threshold1[i]
  est_diff1 <- as.numeric(pvij1 >= th1)
  
  th2 <- threshold2[i]
  est_diff2 <- as.numeric(pvij2 >= th2)
  
  th3 <- threshold3[i]
  est_diff3 <- as.numeric(pvij3 >= th3)
  
  th4 <- threshold4[i]
  est_diff4 <- as.numeric(pvij4 >= th4)
}

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

################################################################################

ca_counties_sf$NAME <- tolower(ca_counties_sf$NAME)

ch_edge1 <- pvij1[which(est_diff1 == 1)]

borders_ca1 <- names(pvij1) %>% 
  str_split(", ") %>% 
  map_dfr(~data.frame(col1 = .x[1], col2 = .x[2]))

edge_plot1 <- ggplot(data = ca_counties_sf) +
  geom_sf(aes(fill = apply(phis1,2,mean)), color = "black", alpha = 0.6) + 
  scale_fill_gradientn(colours = color) +
  ggtitle(paste("Lung (", T1, ")", sep=""))

for(i in which(est_diff1 == 1)){
  shared_border1 <- st_intersection(ca_counties_sf[ca_counties_sf$NAME==borders_ca1[i,1],],
                                    ca_counties_sf[ca_counties_sf$NAME==borders_ca1[i,2],],
                                    model = "closed") 
  edge_plot1 <- edge_plot1 + geom_sf(data = shared_border1, color = "blue", 
                                     linewidth = ((pvij1[i]- min(ch_edge1))/(max(ch_edge1) - min(ch_edge1))) + 1)
}

edge_plot1 <- edge_plot1 +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, face="bold",hjust = 0.5),
        legend.text = element_text(size = 30),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank())

edge_plot1

## 

ch_edge2 <- pvij2[which(est_diff2 == 1)]

borders_ca2 <- names(pvij2) %>% 
  str_split(", ") %>% 
  map_dfr(~data.frame(col1 = .x[1], col2 = .x[2]))

edge_plot2 <- ggplot(data = ca_counties_sf) +
  geom_sf(aes(fill = apply(phis2,2,mean)), color = "black", alpha = 0.6) + 
  scale_fill_gradientn(colours = color) +
  ggtitle(paste("Esophageal (", T2, ")", sep=""))

for(i in which(est_diff2 == 1)){
  shared_border2 <- st_intersection(ca_counties_sf[ca_counties_sf$NAME==borders_ca2[i,1],],
                                    ca_counties_sf[ca_counties_sf$NAME==borders_ca2[i,2],],
                                    model = "closed") 
  edge_plot2 <- edge_plot2 + geom_sf(data = shared_border2, color = "blue", 
                                     linewidth = ((pvij2[i]- min(ch_edge2))/(max(ch_edge2) - min(ch_edge2))) + 1)
}

edge_plot2 <- edge_plot2 +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, face="bold",hjust = 0.5),
        legend.text = element_text(size = 30),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank())

edge_plot2

##

ch_edge3 <- pvij3[which(est_diff3 == 1)]

borders_ca3 <- names(pvij3) %>% 
  str_split(", ") %>% 
  map_dfr(~data.frame(col1 = .x[1], col2 = .x[2]))

edge_plot3 <- ggplot(data = ca_counties_sf) +
  geom_sf(aes(fill = apply(phis3,2,mean)), color = "black", alpha = 0.6) + 
  scale_fill_gradientn(colours = color) +
  ggtitle(paste("Larynx (", T3, ")", sep=""))

for(i in which(est_diff3 == 1)){
  shared_border3 <- st_intersection(ca_counties_sf[ca_counties_sf$NAME==borders_ca3[i,1],],
                                    ca_counties_sf[ca_counties_sf$NAME==borders_ca3[i,2],],
                                    model = "closed") 
  edge_plot3 <- edge_plot3 + geom_sf(data = shared_border3, color = "blue", 
                                     linewidth = ((pvij3[i]- min(ch_edge3))/(max(ch_edge3) - min(ch_edge3))) + 1)
}

edge_plot3 <- edge_plot3 +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, face="bold",hjust = 0.5),
        legend.text = element_text(size = 30),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank())

edge_plot3

##

ch_edge4 <- pvij4[which(est_diff4 == 1)]

borders_ca4 <- names(pvij4) %>% 
  str_split(", ") %>% 
  map_dfr(~data.frame(col1 = .x[1], col2 = .x[2]))

edge_plot4 <- ggplot(data = ca_counties_sf) +
  geom_sf(aes(fill = apply(phis4,2,mean)), color = "black", alpha = 0.6) + 
  scale_fill_gradientn(colours = color) +
  ggtitle(paste("Colorectal (", T4, ")", sep=""))

for(i in which(est_diff4 == 1)){
  shared_border4 <- st_intersection(ca_counties_sf[ca_counties_sf$NAME==borders_ca4[i,1],],
                                    ca_counties_sf[ca_counties_sf$NAME==borders_ca4[i,2],],
                                    model = "closed") 
  edge_plot4 <- edge_plot4 + geom_sf(data = shared_border4, color = "blue", 
                                     linewidth = ((pvij4[i]- min(ch_edge4))/(max(ch_edge4) - min(ch_edge4))) + 1)
}

edge_plot4 <- edge_plot4 +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, face="bold",hjust = 0.5),
        legend.text = element_text(size = 30),
        legend.key.size = unit(0.75,"cm"),
        legend.title = element_blank())

edge_plot4


#pdf("est_diff_cancer_dagar_pois.pdf", height = 10, width = 12)
ggarrange(edge_plot1, edge_plot2, edge_plot3, edge_plot4, nrow = 2, ncol = 2)

################################################################################

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
library(sf)

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

edge_plotc13 <- edge_plotc13 +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
  ggtitle(paste("Lung, Layrnx (T = ", T13c, ")", sep=""))

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

edge_plotc23 <- edge_plotc23 +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
  ggtitle(paste("Lung, Layrnx (T = ", T23c, ")", sep=""))

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

edge_plotc34 <- edge_plotc34 +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
  ggtitle(paste("Lung, Layrnx (T = ", T34c, ")", sep=""))

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
est_diff14 <- as.numeric(pvij14 >= threshold14[T14])
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

edge_plot23 <- edge_plot23 + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
  ggtitle(paste("Esophageal, Layrnx (T = ", T23, ")", sep=""))

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

edge_plot34 <- edge_plot34 + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20, face="bold",hjust = 0.5)) +
  ggtitle(paste("Esophageal, Layrnx (T = ", T34, ")", sep=""))

edge_plot34

#pdf("est_crossd_car_pois.pdf", height = 10, width = 15)
ggarrange(edge_plot12, edge_plot13, edge_plot14, 
          edge_plot23, edge_plot24, edge_plot34, nrow = 2, ncol = 3)

# Adjacency modeling -----------------------------------------------------------

W1 <- mcmc_samples$W1
W2 <- mcmc_samples$W2
W3 <- mcmc_samples$W3
W4 <- mcmc_samples$W4

W1.mean <- Reduce("+", W1)/length(W1)

W2.mean <- Reduce("+", W2)/length(W2)

W3.mean <- Reduce("+", W3)/length(W3)

W4.mean <- Reduce("+", W4)/length(W4)

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

adj_plot1 <- ggplot() +
  geom_polygon(data = ca.poly.df,  aes(long, lat, group = group),
               color = "black",
               fill = "white",
               alpha = 0.6)

for(i in 1:nrow(neighbor_list0)){
  if(adj_det1[i]!=1){
    adj_plot1 <- adj_plot1 + 
      geom_path(aes_string(x = path[[i]][,1], y = path[[i]][,2]), color = "blue",
                size = (1-adj_det1[i])*2.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size=35, face="bold",hjust = 0.5)) +
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
                size = (1-adj_det2[i])*2.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size=35, face="bold",hjust = 0.5)) +
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
                size = (1-adj_det3[i])*2.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size=35, face="bold",hjust = 0.5)) +
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
                size = (1-adj_det4[i])*2.5) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size=35, face="bold",hjust = 0.5)) +
      ggtitle("Colorectal non-adjacencies")  
  }
}

adj_plot4

ggarrange(adj_plot1,adj_plot2,adj_plot3,adj_plot4, nrow = 1, ncol = 4)

