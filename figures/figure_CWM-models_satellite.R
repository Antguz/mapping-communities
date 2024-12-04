#'------------------------------------------------------------------------------
#' @title Visualize performance of trait models of CWM
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(bit64)
library(data.table)
library(ggplot2)
library(ggdark)
library(ggh4x)

#'------------------------------------------------------------------------------
#' @Read-performance files
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"

# Read files 
data <- fread(paste0(root_path, "/Models/04_CWM_performance/performance_summary_CWM.csv"))

# Rename parameters
data$perRMSE <- data$perRMSE*100
data <- data[order(plants),]
data[plants == "HTq", plants := "Maximun tree height"]
data[plants == "SLq", plants := "Slenderness"]
data[plants == "WD", plants := "Wood density"]
data[plants == "shade", plants := "Shade tolerance index"]
data[plants == "drought", plants := "Drought tolerance index"]
data[plants == "waterlogging", plants := "Waterlogging tolerance index"]
data[plants == "leaf_area", plants := "Leaf area"]
data[plants == "LMA", plants := "Leaf mass area"]
data[plants == "N", plants := "Nitrogen"]

data$plants <- as.factor(data$plants)
data$plants <- factor(data$plants, levels = rev(c("Maximun tree height",  "Slenderness",
                                                  "Shade tolerance index", "Drought tolerance index", "Waterlogging tolerance index",
                                                  "Leaf area", "Leaf mass area", "Nitrogen", "Wood density")))

#-------------------------------------------------------------------------------
# Plot details

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.5, 2.0, 0.5, 0.5), "pt"),
        legend.position= c("top"), 
        strip.background = element_rect(color="black", 
                                        fill="black", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "white"))

th_trans <- dark_theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position= c("top"),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth=0),
        plot.margin = unit(c(0.5, 2.0, 0.5, 0.5), "pt"),
        strip.background = element_rect(color="white", 
                                        fill="white", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "black"))

th <- th_trans
th <- th_black
line <- "black"
text <- "black"

#-------------------------------------------------------------------------------
#' @Plot 

# Melt dataset
summary_performance <- melt(data, 
                            id.vars = c("sensor", "dataset", "plants"),
                            measure.vars = c("R2", "BIAS", "RMSE", "perRMSE"),
                            variable.name = "metric",
                            value.name = "performance")

summary_performance$dataset <- as.factor(summary_performance$dataset)
summary_performance$dataset <- factor(summary_performance$dataset, 
                                      levels = c("training", "testing"),
                                      labels = c("Training", "Testing"))

summary_performance$metric <- as.factor(summary_performance$metric)
summary_performance$metric <- factor(summary_performance$metric, 
                                     levels = c("R2", "BIAS", "RMSE", "perRMSE"))

# Summary regardless of trait
summary_all <- summary_performance[, .(mean = mean(performance)),
                                   by = c("sensor", "dataset", "metric")]

# Summary of results
summary_performance <- summary_performance[, .(mean = mean(performance)),
                                           by = c("sensor", "dataset", "metric", "plants")]


my_labeller <- as_labeller(c(R2= "italic(R)^2", BIAS ="BIAS", RMSE = "RMSE", perRMSE = "\"%\"*RMSE"),
                           default = label_parsed)

p <- ggplot(summary_performance) +
  geom_point(aes(x = mean, 
                 y = plants,
                 colour = interaction(sensor, dataset),
                 fill = interaction(sensor, dataset),
                 shape = interaction(sensor, dataset)),
             size = 3,
             #colour = "grey",
             position = position_jitterdodge(dodge.width=0.4, jitter.width=0),
             alpha = 0.75) +
  geom_vline(data = summary_all,
             aes(xintercept = mean, 
                 colour = interaction(sensor, dataset),
                 linetype = interaction(sensor, dataset))) +
  scale_shape_manual("Dataset", values = c(24, 21, 24, 21)) +
  scale_colour_manual("Dataset", values = c("#e66101", "#5e3c99", "#e66101", "#5e3c99")) +
  scale_fill_manual("Dataset", values = c("#e66101", "#5e3c99", "transparent", "transparent")) +
  scale_linetype_manual("Dataset", values = c("dotted", "dotted", "dashed", "dashed")) +
  xlab("Values") + ylab("Community-weigted mean") +
  facet_grid(. ~ metric, scales = "free", labeller = my_labeller) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(-0.07, 0.07), breaks = c(-0.05, 0.0, 0.05)),
                               scale_x_continuous(limits = c(0, 2), breaks = c(0, 1.0, 2)),
                               scale_x_continuous(limits = c(15, 25), breaks = c(15, 20, 25)))) +
  th

p 

jpeg(paste0(root_path, "/figures/Manuscript_black/satellite_CWM-models.jpg"), 
     width = 210, 
     height = 150, 
     units = "mm", 
     res = 600, 
     bg = "white",
     #bg = "transparent", 
     pointsize = 14)
p
dev.off()


png(paste0(root_path, "/figures/Manuscript_white/satellite_CWM-models.png"), 
    width = 210, 
    height = 200, 
    units = "mm", 
    res = 600, 
    bg = "transparent",
    pointsize = 14)
p
dev.off()

