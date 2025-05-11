rm(list = ls())

rmse_unstructured_unstructured <- readRDS("rmse/rmse_unstructured_unstructured.rds")
rmse_unstructured_directed <- readRDS("rmse/rmse_unstructured_directed.rds")
rmse_unstructured_undirected <- readRDS("rmse/rmse_unstructured_undirected.rds")

rmse_directed_unstructured <- readRDS("rmse/rmse_directed_unstructured.rds")
rmse_directed_directed <- readRDS("rmse/rmse_directed_directed.rds")
rmse_directed_undirected <- readRDS("rmse/rmse_directed_undirected.rds")

rmse_undirected_unstructured <- readRDS("rmse/rmse_undirected_unstructured.rds")
rmse_undirected_directed <- readRDS("rmse/rmse_undirected_directed.rds")
rmse_undirected_undirected <- readRDS("rmse/rmse_undirected_undirected.rds")

rmse <- cbind(rmse_unstructured_unstructured, rmse_unstructured_directed, rmse_unstructured_undirected,
              rmse_directed_unstructured, rmse_directed_directed, rmse_directed_undirected,
              rmse_undirected_unstructured, rmse_undirected_directed, rmse_undirected_undirected)

rownames(rmse) <- c("beta",
                    "phi1", "phi2", "phi3", "phi4",
                    "gamma1", "gamma2", "gamma3", "gamma4",
                    "theta", 
                    "taud",
                    "V",
                    "tau",
                    "rho",
                    "eta",
                    "rhodis")

colnames(rmse) <- rep(c("Unstructured", "Directed", "Undirected"),3)

round(rmse,3)
