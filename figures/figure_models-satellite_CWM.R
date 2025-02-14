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
root_path <- "E:/MMDB"

# Read files 
data <- fread(paste0(root_path, "/Models/04_CWM_performance/performance_summary_CWM.csv"))

# Rename parameters
data$perRMSE <- data$perRMSE*100
data[plant == "LMA",]$RMSE <- data[plant == "LMA",]$RMSE/10
data[plant == "LMA",]$BIAS <- data[plant == "LMA",]$BIAS/10

data <- data[order(plant),]
data[plant == "HTq", plant := "Tree height (max)"]
data[plant == "SLq", plant := "Slenderness"]
data[plant == "WD", plant := "Wood density"]
data[plant == "shade", plant := "Shade tolerance index"]
data[plant == "drought", plant := "Drought tolerance index"]
data[plant == "waterlogging", plant := "Waterlogging tolerance index"]
data[plant == "LMA", plant := "Leaf Mass Area"]
data[plant == "N", plant := "Nitrogen"]

data$plant <- as.factor(data$plant)
data$plant <- factor(data$plant, levels = rev(c("Tree height (max)", 
                                                "Slenderness", 
                                                "Shade tolerance index", 
                                                "Drought tolerance index",
                                                "Waterlogging tolerance index",
                                                "Leaf Mass Area",
                                                "Nitrogen",
                                                "Wood density")))

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
                            id.vars = c("sensor", "dataset", "plant"),
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

# CI
get_lower <- function(x, confidence = 0.95) {
  n <- length(x)  # Sample size
  mean_x <- mean(x)  # Sample mean
  std_error <- sd(x) / sqrt(n)  # Standard error
  t_value <- qt((1 - confidence) / 2, df = n - 1, lower.tail = TRUE)  # Critical t-value
  
  lower_ci <- mean_x + t_value * std_error  # Lower confidence bound
  return(lower_ci)
}

# Summary of results
summary_performance <- summary_performance[, .(mean = mean(performance), 
                                               sd = sd(performance),
                                               lower = get_lower(performance, confidence = 0.05),
                                               upper = get_lower(performance, confidence = 0.95)),
                                           by = c("sensor", "dataset", "metric", "plant")]


my_labeller <- as_labeller(c(R2= "italic(R)^2", BIAS ="BIAS", RMSE = "RMSE", perRMSE = "\"%\"*RMSE"),
                           default = label_parsed)

p <- ggplot(summary_performance) +
  geom_linerange(aes(xmin = mean - sd,
                     xmax = mean + sd,
                     y = plant,
                     colour = interaction(sensor, dataset)),
                 position = position_dodge(width = 0.4),
                 alpha = 0.75) +
  geom_point(aes(x = mean, 
                 y = plant,
                 colour = interaction(sensor, dataset),
                 fill = interaction(sensor, dataset),
                 shape = interaction(sensor, dataset)),
             size = 2,
             position = position_dodge(width = 0.4),
             alpha = 0.75) +
  geom_vline(data = summary_all,
             aes(xintercept = mean, 
                 colour = interaction(sensor, dataset),
                 linetype = interaction(sensor, dataset))) +
  scale_shape_manual("Dataset", values = c(24, 21, 24, 21)) +
  scale_colour_manual("Dataset", values = c("#8c510a", "#01665e", "#8c510a", "#01665e")) +
  scale_fill_manual("Dataset", values = c("#8c510a", "#01665e", "transparent", "transparent")) +
  scale_linetype_manual("Dataset", values = c("dotted", "dotted", "dashed", "dashed")) +
  xlab("Values") + ylab("Community-weigted mean") +
  facet_grid(. ~ metric, scales = "free", labeller = my_labeller) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(-0.07, 0.07), breaks = c(-0.05, 0.0, 0.05)),
                               scale_x_continuous(limits = c(0, 2), breaks = c(0, 1.0, 2)),
                               scale_x_continuous(limits = c(10, 25), breaks = c(10, 17, 25)))) +
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

