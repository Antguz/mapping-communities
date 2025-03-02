#'------------------------------------------------------------------------------
#' @title Visualize performance of community-weighted mean plant traits
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

root_path <- "/media/antonio/Work/MMDB"
root_path <- "E:/MMDB"

# Model performance results 
data <- fread(paste0(root_path, "/Models/03_model_CWM/func_model-CWM.csv"))

# Change of names
data <- data[order(plant),]
data[plant == "HTq", name := "Tree height (max)"]
data[plant == "SLq", name := "Slenderness"]
data[plant == "WD", name := "Wood density"]
data[plant == "shade", name := "Shade tolerance index"]
data[plant == "drought", name := "Drought tolerance index"]
data[plant == "waterlogging", name := "Waterlogging tolerance index"]
data[plant == "LMA", name := "Leaf Mass Area"]
data[plant == "N", name := "Nitrogen"]

data$name <- as.factor(data$name)
data$name <- factor(data$name, 
                      levels = rev(c("Tree height (max)", 
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

#th <- th_trans
th <- th_black
line <- "black"
text <- "black"

#-------------------------------------------------------------------------------
#' @Plot 

# Melt dataset
summary_performance <- melt(data, 
                            id.vars = c("dataset", "name"),
                            measure.vars = c("R2", "BIAS", "RMSE", "perRMSE"),
                            variable.name = "metric",
                            value.name = "performance")

summary_performance$dataset <- as.factor(summary_performance$dataset)
summary_performance$dataset <- factor(summary_performance$dataset, 
                                      levels = c("training", "testing"),
                                      labels = c("Training", "Testing"))

#summary_performance <- summary_performance[metric == "R2" | metric == "perRMSE", ]

summary_performance$metric <- as.factor(summary_performance$metric)
summary_performance$metric <- factor(summary_performance$metric, 
                                     levels = c("R2", "BIAS", "RMSE", "perRMSE"))

# Summary regardless of trait
summary_all <- summary_performance[, .(mean = mean(performance)),
                                   by = c("dataset", "metric")]

# Summary of results
summary_performance <- summary_performance[, .(mean = mean(performance)),
                                           #lower = t.test(performance)$conf.int[1],
                                           #upper = t.test(performance)$conf.int[2]),
                                           by = c("dataset", "metric", "name")]

my_labeller <- as_labeller(c(R2= "italic(R)^2", BIAS ="BIAS", RMSE = "RMSE", perRMSE = "\"%\"*RMSE"),
                           default = label_parsed)

p <- ggplot(summary_performance) +
  geom_point(aes(x = mean, 
                 y = name,
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
  xlab("Values") + ylab("Community-weighted mean plant traits") +
  facet_grid(. ~ metric, scales = "free", labeller = my_labeller) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(0, 1.0), breaks = c(0.0, 0.5, 1.0)))) +
  th

p

jpeg(paste0(root_path, "/figures/FIA-models_lineages.jpg"), 
     width = 210, 
     height = 200, 
     units = "mm", 
     res = 600, 
     bg = "white",
     #bg = "transparent", 
     pointsize = 14)
p
dev.off()


png(paste0(root_path, "/figures/FIA-models_lineages.png"), 
    width = 210, 
    height = 200, 
    units = "mm", 
    res = 600, 
    bg = "white",
    #bg = "transparent", 
    pointsize = 14)
p
dev.off()