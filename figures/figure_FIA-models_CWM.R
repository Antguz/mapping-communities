#'------------------------------------------------------------------------------
#' @title Visualize performance of trait models of CWM
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
library(ggh4x)

#'------------------------------------------------------------------------------
#' @Read-performance files
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/DESIS_mapping"

# Categorical 
data <- fread(paste0(root_path, "/07_model_CWM/func_model-CWM.csv"))
data$perRMSE <- data$perRMSE*100
data <- data[order(plant),]
data[plant == "HTq", plant := "Tree height (max)"]
data[plant == "SLq", plant := "Slenderness"]
data[plant == "CCLCD_mode", plant := "Crown class (mode)"]
data[plant == "WOOD_SPGR_GREENVOL_DRYWT", plant := "Wood density"]
data[plant == "MC_PCT_GREEN_WOOD", plant := "Moisture content"]
data[plant == "CARBON_RATIO_LIVE", plant := "Carbon fraction"]
data[plant == "shade", plant := "Shade tolerance index"]
data[plant == "drought", plant := "Drought tolerance index"]
data[plant == "waterlogging", plant := "Waterlogging tolerance index"]
data[plant == "leaf_area", plant := "Leaf area"]
data[plant == "LMA", plant := "Leaf mass area"]
data[plant == "N", plant := "Nitrogen"]

data$plant <- as.factor(data$plant)
data$plant <- factor(data$plant, levels = rev(c("Tree height (max)",  "Crown class (mode)", "Slenderness",
                                              "Carbon fraction", "Moisture content", 
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

#th <- th_trans
th <- th_black
line <- "black"
text <- "black"

#-------------------------------------------------------------------------------
#' @Plot 

# Melt dataset
summary_performance <- melt(data[plant != "Leaf area"], 
                            id.vars = c("dataset", "plant"),
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
                                   by = c("dataset", "metric")]

# Summary of results
summary_performance <- summary_performance[, .(mean = mean(performance)),
                                           by = c("dataset", "metric", "plant")]


my_labeller <- as_labeller(c(R2= "italic(R)^2", BIAS ="BIAS", RMSE = "RMSE", perRMSE = "\"%\"*RMSE"),
                           default = label_parsed)

p <- ggplot(summary_performance) +
  geom_point(aes(x = mean, 
                 y = plant,
                 #colour = dataset,
                 fill = dataset,
                 shape = dataset),
             size = 2,
             colour = "grey",
             alpha = 0.75) +
  geom_vline(data = summary_all,
             aes(xintercept = mean, 
                 colour = dataset,
                 linetype = dataset)) +
  scale_shape_manual("Dataset", values = c(24, 25)) +
  scale_colour_manual("Dataset", values = c("#e66101", "#5e3c99")) +
  scale_fill_manual("Dataset", values = c("#e66101", "#5e3c99")) +
  scale_linetype_manual("Dataset", values = c("dotted", "dashed")) +
  xlab("Values") + ylab("Community-weigted mean") +
  facet_grid(. ~ metric, scales = "free", labeller = my_labeller) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(-0.25, 0.25), breaks = c(-0.2, 0.0, 0.2)),
                               scale_x_continuous(limits = c(0, 12), breaks = c(0, 5, 10)),
                               scale_x_continuous(limits = c(5, 25), breaks = c(5, 15, 25)))) +
  th

p 

jpeg(paste0(root_path, "/figures/FIA-models_CWM.jpg"), 
     width = 210, 
     height = 150, 
     units = "mm", 
     res = 600, 
     bg = "white",
     #bg = "transparent", 
     pointsize = 14)
p
dev.off()


png(paste0(root_path, "/figures/FIA-models_CWM.png"), 
    width = 210, 
    height = 200, 
    units = "mm", 
    res = 600, 
    #bg = "transparent",
    bg = "white",
    pointsize = 14)
p
dev.off()

