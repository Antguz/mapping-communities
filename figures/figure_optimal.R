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
#' Read information

#' Folder root
#root_path <- "E:/Projects/DESIS_mapping/06_model/SW_beta_PLS"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping/06_model/SW_beta_PLS"

# Get RMSEP
rmsep <- fread(paste0(root_path, "/optimal_models.csv"))
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
  scale_x_continuous(expand = c(0, 0), limits = c(1, 60)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0.18, 0.218)) +
  geom_vline(xintercept = c(55,22,21), colour = col, linetype = "dotted", linewidth = 1.0) +
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




# Black background
plot_optimal <- ggplot(summary_rmsep) +
  geom_ribbon(aes(x = ncomp, ymin = mean - sd, ymax = mean + sd, fill = MDS), alpha = 0.25) +
  geom_line(aes(x = ncomp, y = mean, colour = MDS), linewidth = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 60)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.31, 0.49)) +
  geom_vline(xintercept = c(55,22,21), colour = col, linetype = "dotted", linewidth = 1.0) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name= "", values = col) +
  xlab("Number of components") +
  ylab("RMSEP") + 
  dark_theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.90, 0.90), 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text=element_text(color="white"), 
        plot.background = element_rect(colour = "transparent", 
                                       fill = "transparent", 
                                       size=0)) +
  labs(subtitle = "Optimal number of components")

png("figures/pls_optimal_black.png", width = 140, height = 90, units = "mm", res = 600, bg = "transparent", pointsize = 12)
plot_optimal
dev.off()

ggplot(rmsep[MDS == "Axis 1"]) +
  geom_line(aes(x = ncomp, y = Rsquared, group = interaction(MDS, iteration)), colour = "black") +
  geom_vline(xintercept = c(38,22,21), colour = col, linetype = "dotted", linewidth = 1.0) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 60))

ggplot(rmsep[MDS == "Axis 1"]) +
  geom_line(aes(x = ncomp, y = RMSE, group = interaction(MDS, iteration)), colour = "black") +
  geom_vline(xintercept = c(38,22,14), colour = col, linetype = "dotted", linewidth = 1.0) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 60))

ggplot(rmsep[MDS == "Axis 1"]) +
  geom_line(aes(x = ncomp, y = MAE, group = interaction(MDS, iteration)), colour = "black") +
  geom_vline(xintercept = c(38,22,14), colour = col, linetype = "dotted", linewidth = 1.0) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 60))
