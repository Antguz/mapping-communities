#'------------------------------------------------------------------------------
#' @title Codes to visualize the optimal sigma
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(ggplot2)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Read information


root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme Pro/Projects/DESIS_mapping"

#Read sigma models
models <- readRDS(paste0(root_path, "/06_model/optimal_sigma_models.rds"))


optimal <- models[, .(Rsq = mean(Rsquared), Rsq_SD = sd(Rsquared), 
                      RMSE = mean(RMSE), RMSE_SD = sd(RMSE)), by = c("MDS", "sigma")]
sigma_optimal <- optimal[, which.min(RMSE), by = "MDS"]$V1
sigma_optimal <- seq(0.001, 0.015, by = 0.0005)[sigma_optimal]

ggplot(optimal) +
  geom_ribbon(aes(x = sigma, y = RMSE, ymax = RMSE + RMSE_SD, ymin = RMSE - RMSE_SD, fill = MDS), alpha=0.2) +
  geom_line(aes(x = sigma, y = RMSE, colour = MDS)) +
  geom_point(aes(x = sigma, y = RMSE, colour = MDS, shape = MDS), size = 2) + 
  geom_vline(xintercept = sigma_optimal[1], colour = "red", linetype = "dotted") +
  geom_vline(xintercept = sigma_optimal[2], colour = "green", linetype = "dotted") +
  geom_vline(xintercept = sigma_optimal[3], colour = "blue", linetype = "dotted") +
  coord_cartesian(xlim = c(0.001, 0.015)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab("Sigma") + ylab("RMSE (Spatial Cross-Validation)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

