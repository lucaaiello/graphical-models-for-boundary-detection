rm(list = ls())

library(dplyr)
library(caret)
library(Matrix)
library(LaplacesDemon)
library(matrixStats)
library(monomvn)

setwd("C:/Dati/Dottorato/Visiting UCLA/Spatial Disease Mapping/Bibliografia/gao_2022/simulation/Disease graph/directed")

###### Extracting objects from the simulation branches #####

# load .rds files with the results

# discrete random effects for diseases

phi_true1 <- readRDS("phi_true1_try.rds")
phi_true2 <- readRDS("phi_true2_try.rds")
phi_true3 <- readRDS("phi_true3_try.rds")
phi_true4 <- readRDS("phi_true4_try.rds")

# "disease difference" in the same region

disease_diff12 <- as.numeric(phi_true1 != phi_true2)
disease_diff13 <- as.numeric(phi_true1 != phi_true3)
disease_diff14 <- as.numeric(phi_true1 != phi_true4)
disease_diff23 <- as.numeric(phi_true2 != phi_true3)
disease_diff24 <- as.numeric(phi_true2 != phi_true4)
disease_diff34 <- as.numeric(phi_true3 != phi_true4)

sum(disease_diff12)
sum(disease_diff13)
sum(disease_diff14)
sum(disease_diff23)
sum(disease_diff24)
sum(disease_diff34)

# specificity and specificity for difference boundaries

spec1_aux <- matrix(0,50,10)
spec2_aux <- matrix(0,50,10)
spec3_aux <- matrix(0,50,10)
spec4_aux <- matrix(0,50,10)

sens1_aux <- matrix(0,50,10)
sens2_aux <- matrix(0,50,10)
sens3_aux <- matrix(0,50,10)
sens4_aux <- matrix(0,50,10)

for(seed in 1:50){
  print(seed)
  
  if(seed >=1 & seed<=10){
    attach(my_sim_1_10_results)
  } else if(seed >=11 & seed<=20){
    attach(my_sim_11_20_results)
  } else if(seed >=21 & seed<=30){
    attach(my_sim_21_30_results)
  } else if(seed >=31 & seed<=40){
    attach(my_sim_31_40_results)
  } else if(seed >=41 & seed<=50){
    attach(my_sim_41_50_results)
  }
  
  spec1_aux[seed,] <- spec1[seed,]
  spec2_aux[seed,] <- spec2[seed,]
  spec3_aux[seed,] <- spec3[seed,]
  spec4_aux[seed,] <- spec4[seed,]
  
  sens1_aux[seed,] <- sens1[seed,]
  sens2_aux[seed,] <- sens2[seed,]
  sens3_aux[seed,] <- sens3[seed,]
  sens4_aux[seed,] <- sens4[seed,]
  
  if(seed >=1 & seed<=10){
    detach(my_sim_1_10_results)
  } else if(seed >=11 & seed<=20){
    detach(my_sim_11_20_results)
  } else if(seed >=21 & seed<=30){
    detach(my_sim_21_30_results)
  } else if(seed >=31 & seed<=40){
    detach(my_sim_31_40_results)
  } else if(seed >=41 & seed<=50){
    detach(my_sim_41_50_results)
  }
  
}

spec1_mean_dagar <- colMeans(spec1_aux)
sens1_mean_dagar <- colMeans(sens1_aux)

spec2_mean_dagar <- colMeans(spec2_aux)
sens2_mean_dagar <- colMeans(sens2_aux)

spec3_mean_dagar <- colMeans(spec3_aux)
sens3_mean_dagar <- colMeans(sens3_aux)

spec4_mean_dagar <- colMeans(spec4_aux)
sens4_mean_dagar <- colMeans(sens4_aux)

table <- cbind(spec1_mean_dagar, sens1_mean_dagar, 
               spec2_mean_dagar, sens2_mean_dagar,
               spec3_mean_dagar, sens3_mean_dagar, 
               spec4_mean_dagar, sens4_mean_dagar)

# sensitivity and specificity for adjacencies

specW1_aux <- vector(length = 50)
specW2_aux <- vector(length = 50)
specW3_aux <- vector(length = 50)
specW4_aux <- vector(length = 50)

sensW1_aux <- vector(length = 50)
sensW2_aux <- vector(length = 50)
sensW3_aux <- vector(length = 50)
sensW4_aux <- vector(length = 50)



for(seed in 1:50){
  print(seed)
  
  if(seed >=1 & seed<=10){
    attach(my_sim_1_10_results)
  } else if(seed >=11 & seed<=20){
    attach(my_sim_11_20_results)
  } else if(seed >=21 & seed<=30){
    attach(my_sim_21_30_results)
  } else if(seed >=31 & seed<=40){
    attach(my_sim_31_40_results)
  } else if(seed >=41 & seed<=50){
    attach(my_sim_41_50_results)
  }
  
  specW1_aux[seed] <- specW1[seed]
  specW2_aux[seed] <- specW2[seed]
  specW3_aux[seed] <- specW3[seed]
  specW4_aux[seed] <- specW4[seed]
  
  sensW1_aux[seed] <- sensW1[seed]
  sensW2_aux[seed] <- sensW2[seed]
  sensW3_aux[seed] <- sensW3[seed]
  sensW4_aux[seed] <- sensW4[seed]
  
  if(seed >=1 & seed<=10){
    detach(my_sim_1_10_results)
  } else if(seed >=11 & seed<=20){
    detach(my_sim_11_20_results)
  } else if(seed >=21 & seed<=30){
    detach(my_sim_21_30_results)
  } else if(seed >=31 & seed<=40){
    detach(my_sim_31_40_results)
  } else if(seed >=41 & seed<=50){
    detach(my_sim_41_50_results)
  }
  
}

specW1_mean_dagar <- mean(specW1_aux)
sensW1_mean_dagar <- mean(sensW1_aux)

specW2_mean_dagar <- mean(specW2_aux)
sensW2_mean_dagar <- mean(sensW2_aux)

specW3_mean_dagar <- mean(specW3_aux)
sensW3_mean_dagar <- mean(sensW3_aux)

specW4_mean_dagar <- mean(specW4_aux)
sensW4_mean_dagar <- mean(sensW4_aux)

