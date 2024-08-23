#'------------------------------------------------------------------------------
#' @title Visualize performance of occurrence of plant groups
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(bit64)
library(data.table)
library(rcartocolor)
library(viridis)
library(ggplot2)
library(ggdark)

#'------------------------------------------------------------------------------
#' @Read-performance files
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/DESIS_mapping"

# Categorical 
data <- fread(paste0(root_path, "/08_CWM_performance/performance_summary_CWM.csv"))
data$perRMSE <- data$perRMSE*100
data$name <- data$plants
data[plants == "HTq", name := "Tree height (max)"]
data[plants == "SLq", name := "Slenderness"]
data[plants == "CCLCD_mode", name := "Crown class (mode)"]
data[plants == "WOOD_SPGR_GREENVOL_DRYWT", name := "Wood density"]
data[plants == "MC_PCT_GREEN_WOOD", name := "Moisture content"]
data[plants == "CARBON_RATIO_LIVE", name := "Carbon fraction"]
data[plants == "shade", name := "Shade tolerance index"]
data[plants == "drought", name := "Drought tolerance index"]
data[plants == "waterlogging", name := "Waterlogging tolerance index"]
data[plants == "leaf_area", name := "Leaf area"]
data[plants == "LMA", name := "Leaf mass area"]
data[plants == "N", name := "Nitrogen"]

data$name <- as.factor(data$name)
data$name <- factor(data$name, levels = rev(c("Tree height (max)",  "Crown class (mode)", "Slenderness",
                                                  "Carbon fraction", "Moisture content", 
                                                  "Shade tolerance index", "Drought tolerance index", "Waterlogging tolerance index",
                                                  "Leaf area", "Leaf mass area", "Nitrogen", "Wood density")))

# Downscale leaf area
data <- data[plants != "leaf_area"]

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

#th <- th_trans
th <- th_black
line <- "black"
text <- "black"

#-------------------------------------------------------------------------------
#' @Plot 

# Melt dataset
summary_performance <- melt(data, 
                            id.vars = c("dataset", "sensor", "name"),
                            measure.vars = c("R2", "BIAS", "RMSE", "perRMSE"),
                            variable.name = "metric",
                            value.name = "performance")

summary_performance$dataset <- as.factor(summary_performance$dataset)
summary_performance$dataset <- factor(summary_performance$dataset, 
                                      levels = c("training", "testing"),
                                      labels = c("Training", "Testing"))

summary_performance$metric <- as.factor(summary_performance$metric)
summary_performance$metric <- factor(summary_performance$metric, 
                                     levels = c("R2", "BIAS", "RMSE", "perRMSE"),
                                     labels = c("R2", "BIAS", "RMSE", "perRMSE"))

# Summary regardless of trait
summary_all <- summary_performance[, .(mean = mean(performance)),
                                   by = c("dataset", "sensor", "metric")]

# Summary of results
summary_performance <- summary_performance[, .(mean = mean(performance),
                                               sd = sd(performance)),
                            by = c("dataset", "sensor", "metric", "name")]

my_labeller <- as_labeller(c(R2= "italic(R)^2", BIAS ="BIAS", RMSE = "RMSE", perRMSE = "\"%\"*RMSE"),
                           default = label_parsed)


p <- ggplot(summary_performance) +
  geom_vline(data = summary_all,
             aes(xintercept = mean, 
                 colour = sensor,
                 linetype = dataset)) +
  geom_point(aes(x = mean, 
                 y = name,
                 colour = sensor,
                 fill = dataset,
                 shape = dataset),
             size = 2,
             alpha = 0.75) +
  scale_shape_manual("Dataset", values = c(24, 25)) +
  scale_colour_manual("Sensor", values = c("#a6611a", "#018571")) +
  scale_fill_manual("Dataset", values = c("black", "grey")) +
  scale_linetype_manual("Dataset", values = c("dotted", "dashed")) +
  xlab("Values") + ylab("Community-weigted mean") +
  facet_grid(. ~ metric, scales = "free", labeller = my_labeller) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(-0.7, 0.7), breaks = c(-0.5, 0.0, 0.5)),
                               scale_x_continuous(limits = c(0, 25), breaks = c(0, 10, 20)),
                               scale_x_continuous(limits = c(10, 30), breaks = c(10, 20, 30)))) +
  th

p

jpeg(paste0(root_path, "/figures/sensors-model_CWM.jpg"), 
    width = 210, 
    height = 150, 
    units = "mm", 
    res = 600, 
    #bg = "transparent",
    bg = "white",
    pointsize = 14)
p
dev.off()

png(paste0(root_path, "/figures/sensors-model_CWM.png"), 
    width = 210, 
    height = 200, 
    units = "mm", 
    res = 600, 
    #bg = "transparent",
    bg = "white",
    pointsize = 14)
p
dev.off()

