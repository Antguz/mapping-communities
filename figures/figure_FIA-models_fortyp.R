#'------------------------------------------------------------------------------
#' @title Visualize performance of occurrence models of forest types
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
data <- fread(paste0(root_path, "/07_model_fortyp/taxo_model-FORTYP.csv"))
data <- data[order(plant),]
data$plants <- as.factor(data$plant)
data$plants <- factor(data$plants, levels = rev(as.character(unique(data$plants))))

#-------------------------------------------------------------------------------
# Plot details

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.5, 2, 0.5, 0.5), "pt"),
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
        plot.margin = unit(c(0.5, 2, 0.5, 0.5), "pt"),
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
                            id.vars = c("dataset", "plants"),
                            measure.vars = c("Accuracy", "Sens", "Spec", "ROC"),
                            variable.name = "metric",
                            value.name = "performance")

summary_performance$dataset <- as.factor(summary_performance$dataset)
summary_performance$dataset <- factor(summary_performance$dataset, 
                                      levels = c("training", "testing"),
                                      labels = c("Training", "Testing"))

#summary_performance <- summary_performance[metric == "R2" | metric == "perRMSE", ]

summary_performance$metric <- as.factor(summary_performance$metric)
summary_performance$metric <- factor(summary_performance$metric, 
                                     levels = c("Accuracy", "Sens", "Spec", "ROC"),
                                     labels = c("Accuracy", "Sensitivity", "Specificity", "AUC"))

# Summary regardless of trait
summary_all <- summary_performance[, .(mean = mean(performance)),
                                   by = c("dataset", "metric")]

# Summary of results
summary_performance <- summary_performance[, .(mean = mean(performance)),
                                           by = c("dataset", "plants", "metric")]

p <- ggplot(summary_performance) +
  geom_point(aes(x = mean, 
                 y = plants,
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
  coord_cartesian(xlim = c(0.4, 1.0)) +
  scale_x_continuous(breaks = c(0.4, 0.7, 1.0)) +
  scale_shape_manual("Dataset", values = c(24, 25)) +
  scale_colour_manual("Dataset", values = c("#e66101", "#5e3c99")) +
  scale_fill_manual("Dataset", values = c("#e66101", "#5e3c99")) +
  scale_linetype_manual("Dataset", values = c("dotted", "dashed")) +
  facet_grid(. ~ metric, scales = "free") +
  xlab("Values") + ylab("Forest type") +
  th

jpeg(paste0(root_path, "/figures/FIA-models_fortype.jpeg"), 
     width = 210, 
     height = 240, 
     units = "mm", 
     res = 600, 
     bg = "white",
     #bg = "transparent", 
     pointsize = 14)
p
dev.off()


png(paste0(root_path, "/figures/FIA-models-fortype.png"), 
    width = 220, 
    height = 240, 
    units = "mm", 
    res = 600, 
    bg = "white",
    #bg = "transparent", 
    pointsize = 14)
p
dev.off()

