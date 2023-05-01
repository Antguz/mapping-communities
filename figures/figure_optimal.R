#'------------------------------------------------------------------------------
#' @title Codes to visualize optimal number of components of the models
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Read information

#' Folder root
root_path <- "E:/Projects/DESIS_mapping/06_model/SW_beta_PLS"
#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"

# Get RMSEP
rmsep <- readRDS(paste0(root_path, "/optimal_models.rds"))
rmsep[MDS == "MDS1", MDS := "Axis 1"]
rmsep[MDS == "MDS2", MDS := "Axis 2"]
rmsep[MDS == "MDS3", MDS := "Axis 3"]

# Summary by MDS
summary_rmsep <- rmsep[, .(mean = mean(RMSE), sd = sd(RMSE)), 
                       by = c("MDS", "ncomp")]

# Get summary of coefficients

#-------------------------------------------------------------------------------
# Get plot

col <- c("#fc8d62", "#66c2a5", "#8da0cb")

th <- theme_bw() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill= "transparent"),
            panel.background = element_rect(fill = "transparent"),
            rect = element_rect(fill = "transparent"))

plot_optimal <- ggplot(summary_rmsep) +
  geom_ribbon(aes(x = ncomp, ymin = mean - sd, ymax = mean + sd, fill = MDS), alpha = 0.25) +
  geom_line(aes(x = ncomp, y = mean, colour = MDS), linewidth = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 20)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.15, 0.23)) +
  geom_vline(xintercept = c(5,4,15), colour = col, linetype = "dotted", linewidth = 1.0) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name= "", values = col) +
  xlab("Number of components") +
  ylab("RMSEP") + 
  th + theme(legend.position = c(0.93, 0.90), 
             legend.background = element_rect(fill= "transparent"),
             panel.background = element_rect(fill = "transparent"),
             rect = element_rect(fill = "transparent"))


png("figures/pls_optimal.png", width = 180, height = 90, units = "mm", res = 300, bg = "transparent")
plot_optimal
dev.off()
