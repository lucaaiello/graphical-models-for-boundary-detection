rm(list = ls())

WAIC_unstructured <- readRDS("CAR/WAIC_unstructured.rds")
WAIC_directed <- readRDS("CAR/WAIC_directed.rds")
WAIC_undirected <- readRDS("CAR/WAIC_undirected.rds")

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
  xlim(-50000,2.5e5) +
  xlab("WAIC") + ylab("Density") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1,"cm"))

