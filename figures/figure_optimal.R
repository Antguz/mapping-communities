#'------------------------------------------------------------------------------
#' @title Codes to visualize optimal number of components of the models
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)
library(bit64)
library(ggdark)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
sensor <- "DESIS"
out_path <- paste0(root_path, "/06_model/SW_beta_PLS/", sensor)

#'------------------------------------------------------------------------------
#' Read information

# Get RMSEP
rmsep <- fread(paste0(out_path, "/", sensor, "_optimal_models_test.csv"))
rmsep[MDS == "MDS1", MDS := "Axis 1"]
rmsep[MDS == "MDS2", MDS := "Axis 2"]
rmsep[MDS == "MDS3", MDS := "Axis 3"]

# Summary by MDS
summary_rmsep <- rmsep[, .(mean = mean(RMSE), sd = sd(RMSE)), 
                       by = c("MDS", "ncomp")]

# Mode best tune
# Mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

bestTune <- rmsep[, .(best = getmode(bestTune)), by = c("MDS")]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Plot details

col <- c("#fc8d62", "#66c2a5", "#8da0cb")

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.93, 0.90),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"))

th_trans <- dark_theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.08, 0.90), 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth = 0))

#-------------------------------------------------------------------------------
# Get plot

DESIS <- c(12, 14, 17)
EMIT <- c(20, 12, 9)

plot_optimal <- ggplot(summary_rmsep) +
  geom_ribbon(aes(x = ncomp, ymin = mean - sd, ymax = mean + sd, fill = MDS), alpha = 0.25) +
  geom_line(aes(x = ncomp, y = mean, colour = MDS), linewidth = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 30)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.11, 0.23), n.breaks = 4) +
  geom_vline(xintercept = DESIS, colour = col, linetype = "dotted", linewidth = 1.0) +
  #geom_vline(xintercept = EMIT, colour = col, linetype = "dotted", linewidth = 1.0) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name= "", values = col) +
  xlab("Number of components") +
  ylab("RMSEP") +  
  #th_trans + 
  theme(legend.position = c(0.93, 0.90))

png(paste0(root_path, "/figures/", sensor, "_optimal_transparent.png"), 
    width = 150, 
    height = 90, 
    pointsize = 12, 
    units = "mm", 
    res = 300, 
    bg = "transparent")

plot_optimal

dev.off()
