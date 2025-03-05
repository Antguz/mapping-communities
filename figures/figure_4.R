#'------------------------------------------------------------------------------
#' @title Figure of summary 
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

# Forest type
data <- fread(paste0(root_path, "/Models/04_fortyp_performance/performance_summary_fortyp.csv"))
data <- data[order(plant),]
data$plant <- as.factor(data$plant)
data$plant <- factor(data$plant, levels = rev(as.character(unique(data$plant))))
data$perRMSE <- data$perRMSE*100
data$BIAS <- data$BIAS/100000
data$RMSE <- data$RMSE/100000
data$intercept <- data$intercept/100000
fortyp <- data
fortyp$dimension <- "Forest types"

# Lineage occurrance
data <- fread(paste0(root_path, "/Models/04_occurence_performance/performance_summary_lineages.csv"))
data <- data[order(plant),]
data$perRMSE <- data$perRMSE*100
data$BIAS <- data$BIAS/100000
data$RMSE <- data$RMSE/100000
data$intercept <- data$intercept/100000
lineage <- data
lineage$dimension <- "Plant lineages"

# CWM of plant traits
# Read files 
data <- fread(paste0(root_path, "/Models/04_CWM_performance/performance_summary_CWM.csv"))
data$perRMSE <- data$perRMSE*100
data[plant == "LMA",]$RMSE <- data[plant == "LMA",]$RMSE/100
data[plant == "LMA",]$BIAS <- data[plant == "LMA",]$BIAS/100
data[plant == "N",]$RMSE <- data[plant == "N",]$RMSE/10
data[plant == "N",]$BIAS <- data[plant == "N",]$BIAS/10
data <- data[order(plant),]
data[plant == "HTq", plant := "Tree height (max)"]
data[plant == "SLq", plant := "Slenderness"]
data[plant == "WD", plant := "Wood density"]
data[plant == "shade", plant := "Shade tolerance index"]
data[plant == "drought", plant := "Drought tolerance index"]
data[plant == "waterlogging", plant := "Waterlogging tolerance index"]
data[plant == "LMA", plant := "Leaf Mass Area*"]
data[plant == "N", plant := "Nitrogen"]
cwm <- data
cwm$dimension <- "CWM of plant traits"

data <- rbind(fortyp, lineage, cwm)

#-------------------------------------------------------------------------------
# Plot details

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.5, 2, 0.5, 0.5), "mm"),
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
        plot.margin = unit(c(0.5, 2, 0.5, 0.5), "mm"),
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
                            id.vars = c("sensor", "dataset", "dimension", "plant"),
                            measure.vars = c("R2", "BIAS", "RMSE", "perRMSE"),
                            variable.name = "metric",
                            value.name = "performance")

summary_performance$dataset <- as.factor(summary_performance$dataset)
summary_performance$dataset <- factor(summary_performance$dataset, 
                                      levels = c("training", "testing"),
                                      labels = c("Training", "Testing"))

# Summary all
summary_all <- summary_performance[, .(mean = mean(performance)),
                                   by = c("sensor", "dataset", "metric")]

summary_performance$metric <- as.factor(summary_performance$metric)
summary_performance$metric <- factor(summary_performance$metric, 
                                     levels = c("R2", "BIAS", "RMSE", "perRMSE"))

summary_performance <- summary_performance[plant == "Black spruce" |
                                             plant == "Red maple / upland"|
                                             plant == "Sugar maple / beech / yellow birch"|
                                             plant == "Pinus" |
                                             plant == "Acer" |
                                             plant == "Fagus" |
                                             plant == "Nitrogen" | 
                                             plant == "Leaf Mass Area*" |
                                             plant == "Wood density", ]

summary_performance$plant <- factor(summary_performance$plant, levels = rev(c("Black spruce",
                                                                        "Red maple / upland",
                                                                        "Sugar maple / beech / yellow birch",
                                                                        "Pinus", "Acer", "Fagus",
                                                                        "Nitrogen", 
                                                                        "Leaf Mass Area*",
                                                                        "Wood density")))

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
                                           by = c("dimension", "sensor", "dataset", "metric", "plant")]
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
                 position = position_dodge(width = 0.6),
                 alpha = 1) +
  geom_point(aes(x = mean, 
                 y = plant,
                 colour = interaction(sensor, dataset),
                 fill = interaction(sensor, dataset),
                 shape = interaction(sensor, dataset)),
             size = 1.5,
             position = position_dodge(width = 0.6),
             alpha = 0.75) +
  geom_vline(data = summary_all,
             aes(xintercept = mean,
                 colour = interaction(sensor, dataset),
                 linetype = interaction(sensor, dataset))) +
  scale_shape_manual("Dataset", values = c(24, 21, 24, 21)) +
  scale_colour_manual("Dataset", values = c("#8c510a", "#01665e", "#8c510a", "#01665e")) +
  scale_fill_manual("Dataset", values = c("#8c510a", "#01665e", "transparent", "transparent")) +
  scale_linetype_manual("Dataset", values = c("dotted", "dotted", "dashed", "dashed")) +
  xlab("Performance metric") + ylab("Forest attributes") +
  facet_grid(. ~ metric, scales = "free", labeller = my_labeller) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(-0.12, 0.12), breaks = c(-0.10, 0.0, 0.10)),
                               scale_x_continuous(limits = c(0, 0.3), breaks = c(0, 0.15, 0.3)),
                               scale_x_continuous(limits = c(5, 35), breaks = c(5, 20, 35)))) +
  th

p 

jpeg(paste0(root_path, "/figures/Manuscript_black/forest-attributes.jpeg"), 
     width = 210, 
     height = 120, 
     units = "mm", 
     res = 600, 
     bg = "white",
     #bg = "transparent", 
     pointsize = 14)

p

dev.off()


png(paste0(root_path, "/figures/Manuscript_white/satellite_FORTYP-models.png"), 
    width = 220, 
    height = 250, 
    units = "mm", 
    res = 600, 
    bg = "transparent", 
    pointsize = 14)

p

dev.off()
