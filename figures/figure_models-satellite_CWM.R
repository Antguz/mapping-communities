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
data[plant == "LMA",]$RMSE <- data[plant == "LMA",]$RMSE/10
data[plant == "LMA",]$BIAS <- data[plant == "LMA",]$BIAS/10

data <- data[order(plant),]
data[plant == "HTq", plant := "Tree height (max)"]
data[plant == "SLq", plant := "Slenderness"]
data[plant == "WD", plant := "Wood density"]
data[plant == "shade", plant := "Shade tolerance index"]
data[plant == "drought", plant := "Drought tolerance index"]
data[plant == "waterlogging", plant := "Waterlogging tolerance index"]
data[plant == "LMA", plant := "Leaf Mass Area*"]
data[plant == "N", plant := "Nitrogen"]

data$plant <- as.factor(data$plant)
data$plant <- factor(data$plant, levels = rev(c("Tree height (max)", 
                                                "Slenderness", 
                                                "Shade tolerance index", 
                                                "Drought tolerance index",
                                                "Waterlogging tolerance index",
                                                "Leaf Mass Area*",
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
        plot.margin = unit(c(0.5, 2.0, 0.5, 0.5), "mm"),
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
        plot.margin = unit(c(0.5, 2.0, 0.5, 0.5), "mm"),
        strip.background = element_rect(color="white", 
                                        fill="white", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "black"))

th <- th_trans
line <- "white"
text <- "white"

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
get_ci <- function(x, type, confidence = 0.95) {
  n <- length(x)
  mean_x <- mean(x)
  std_error <- sd(x) / sqrt(n)
  z_value <- qnorm((1 - 0.95) / 2)  # Z-score for normal distribution
  
  if(type == "lower") {
    ci <- mean_x - z_value * std_error
  } else if(type == "upper") {
    ci <- mean_x + z_value * std_error
  }
  
  return(ci)
}

# Summary of results
summary_performance <- summary_performance[, .(mean = mean(performance), 
                                               sd = sd(performance),
                                               lower = get_ci(performance, "lower"),
                                               upper = get_ci(performance, "upper")),
                                           by = c("sensor", "dataset", "metric", "plant")]

# Plot labels
my_labeller <- as_labeller(c(R2= "italic(R)^2", BIAS ="BIAS", RMSE = "RMSE", perRMSE = "\"%\"*RMSE"),
                           default = label_parsed)

# Create rects
rects <- data.frame(
  ymin = seq(0.5, length(unique(summary_performance$plant)), by = 2),
  ymax = seq(1.5, length(unique(summary_performance$plant)) + 1, by = 2),
  xmin = -Inf,
  xmax = Inf
)

p <- ggplot(summary_performance) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +  
  geom_linerange(aes(xmin = lower,
                     xmax = upper,
                     y = plant,
                     colour = interaction(sensor, dataset)),
                 position = position_dodge(width = 0.4),
                 alpha = 0.75) +
  geom_point(aes(x = mean, 
                 y = plant,
                 colour = interaction(sensor, dataset),
                 fill = interaction(sensor, dataset),
                 shape = interaction(sensor, dataset)),
             size = 1.5,
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
  xlab("Values") + ylab("Community-weighted mean plant traits") +
  facet_grid(. ~ metric, scales = "free", labeller = my_labeller) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(-0.04, 0.04), breaks = c(-0.03, 0, 0.03)),
                               scale_x_continuous(limits = c(0, 2.5), breaks = c(0, 1, 2)),
                               scale_x_continuous(limits = c(10, 20), breaks = c(10, 15, 20)))) +
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

