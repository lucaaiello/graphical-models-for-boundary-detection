rm(list = ls())

WAIC_unstructured_unstructured <- readRDS("WAIC/WAIC_unstructured_unstructured.rds")
WAIC_unstructured_directed <- readRDS("WAIC/WAIC_unstructured_directed.rds")
WAIC_unstructured_undirected <- readRDS("WAIC/WAIC_unstructured_undirected.rds")

WAIC_directed_unstructured <- readRDS("WAIC/WAIC_directed_unstructured.rds")
WAIC_directed_directed <- readRDS("WAIC/WAIC_directed_directed.rds")
WAIC_directed_undirected <- readRDS("WAIC/WAIC_directed_undirected.rds")

WAIC_undirected_unstructured <- readRDS("WAIC/WAIC_undirected_unstructured.rds")
WAIC_undirected_directed <- readRDS("WAIC/WAIC_undirected_directed.rds")
WAIC_undirected_undirected <- readRDS("WAIC/WAIC_undirected_undirected.rds")



# Unstructured ------------------------------------------------------------

count_unstr_dir <- 0
count_unstr_undir <- 0

for (i in 1:100) {
  if (WAIC_unstructured_directed[i] < WAIC_unstructured_unstructured[i]){
    count_unstr_dir <- count_unstr_dir + 1
  }
  
  if (WAIC_unstructured_undirected[i] < WAIC_unstructured_unstructured[i]){
    count_unstr_undir <- count_unstr_undir + 1
  }
  
}

count_dir_unstr <- 0
count_dir_undir <- 0

for (i in 1:100) {
  if (WAIC_directed_unstructured[i] < WAIC_directed_directed[i]){
    count_dir_unstr <- count_dir_unstr + 1
  }
  
  if (WAIC_directed_undirected[i] < WAIC_directed_directed[i]){
    count_dir_undir <- count_dir_undir + 1
  }
  
}

count_undir_unstr <- 0
count_undir_dir <- 0

for (i in 1:100) {
  if (WAIC_undirected_unstructured[i] < WAIC_undirected_undirected[i]){
    count_undir_unstr <- count_undir_unstr + 1
  }
  
  if (WAIC_undirected_directed[i] < WAIC_undirected_undirected[i]){
    count_undir_dir <- count_undir_dir + 1
  }
  
}

################################################################################

(mean(WAIC_unstructured_directed) - mean(WAIC_unstructured_unstructured))/mean(WAIC_unstructured_unstructured)
(mean(WAIC_unstructured_undirected) - mean(WAIC_unstructured_unstructured))/mean(WAIC_unstructured_unstructured)

(mean(WAIC_directed_unstructured) - mean(WAIC_directed_directed))/mean(WAIC_directed_directed)
(mean(WAIC_directed_undirected) - mean(WAIC_directed_directed))/mean(WAIC_directed_directed)

(mean(WAIC_undirected_unstructured) - mean(WAIC_undirected_undirected))/mean(WAIC_undirected_undirected)
(mean(WAIC_undirected_directed) - mean(WAIC_undirected_undirected))/mean(WAIC_undirected_undirected)

################################################################################

# Plots for WAIC
WAIC_value <- c(WAIC_unstructured, WAIC_directed, WAIC_undirected)
Type <- c(rep("Unstructured",100),rep("Directed",100),rep("Undirected",100))
df <- data.frame(Type)
df$WAIC <- WAIC_value

library(plyr)
WAIC_mu <- ddply(df, "Type", summarise, WAIC.mean=mean(WAIC))
WAIC_median <- ddply(df, "Type", summarise, WAIC.median=median(WAIC))
df1 <- merge(df, WAIC_median, by = "Type")

library(ggplot2)
ggplot(df1, aes(x = WAIC, color = Type, fill = Type)) +
  geom_density(alpha = 0.4, adjust = 1.5, size = 2) + 
  geom_vline(data = df1, aes(xintercept = WAIC.median, color = Type),
             linetype="dashed", size = 2) +
  theme_bw() +
  xlim(-6000,2500) +
  xlab("WAIC") + ylab("Density") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1,"cm"))
