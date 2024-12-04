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

root_path <- "/media/antonio/Work/MMDB"

# Categorical 
cwm <- fread(paste0(root_path, "/Models/04_CWM_performance/performance_summary_CWM.csv"))

# Melt
data <- melt(cwm, 
             id.vars = c("sensor", "plants", "dataset"),
             measure.vars = c("R2", "BIAS", "RMSE", "perRMSE"),
             value.name = "performance",
             variable.name = "metric")

data$plants <- as.factor(data$plants)
data$plants <- factor(data$plants, 
                      levels = rev(c("HTq", "SLq", "WD", "shade", "drought", "waterlogging",
                                     "LMA", "N")),
                      labels = rev(c("Tree height (max)", "Slenderness",
                                                  "Wood density",
                                                  "Shade tolerance index", 
                                                  "Drought tolerance index", 
                                                  "Waterlogging tolerance index",
                                                  "Leaf mass area", "Nitrogen")))

data$dataset <- as.factor(data$dataset)
data$dataset <- factor(data$dataset, 
                       levels = c("training", "testing"),
                       labels = c("Training", "Testing"))

unique(data$plants)

#-------------------------------------------------------------------------------
# Plot details

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"),
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
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"),
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


p <- ggplot(data) +
  geom_point(aes(x = performance, 
                 y = plants, 
                 colour = sensor,
                 fill = interaction(sensor, dataset),
                 shape = interaction(sensor, dataset)),
             size = 2) +
  #coord_cartesian(xlim = c(0.04, 0.5), expand = TRUE) +
  scale_shape_manual("Dataset", values = c(16, 17, 21, 22)) +
  scale_fill_manual("Dataset", values = c(NA, NA, "#e08214", "#8073ac")) +
  scale_colour_manual("Sensor", values = c("#e08214", "#8073ac")) +
  facet_grid(. ~ metric, scales="free") +
  xlab("%RMSE") + ylab("Community-weigted mean") +
  th

png(paste0(root_path, "/figures/CWM_performance.png"), 
    width = 210, 
    height = 200, 
    units = "mm", 
    res = 600, 
    bg = "transparent", 
    pointsize = 14)
p
dev.off()

