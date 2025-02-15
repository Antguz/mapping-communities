#'------------------------------------------------------------------------------
#' @title Visualize performance of occurrence models of forest types
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

# Categorical 
data <- fread(paste0(root_path, "/Models/04_occurence_performance/performance_summary_lineages.csv"))
data <- data[order(plant),]
data$perRMSE <- data$perRMSE*100

# Define levels
data[plant == "AM", plant := "Arbuscular mycorrhizal"]
data[plant == "ECM", plant := "Ectomycorrhizal"]
data[plant == "gymnosperm", plant := "Gymnosperm"]
data[plant == "angiosperm", plant := "Angiosperm"]
data <- data[plant != "deciduous",]
data <- data[plant != "evergreen",]
data <- data[plant != "needle",]
data <- data[plant != "simple",]
data <- data[plant != "compound",]
data <- data[plant != "broadleaf",]

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

# Define factors
summary_performance$plant <- as.factor(summary_performance$plant)
summary_performance$plant <- factor(summary_performance$plant, levels = rev(c("Angiosperm", "Gymnosperm",
                                                                               "Altingiaceae", "Betulaceae", "Cupressaceae", 
                                                                               "Fagaceae", "Juglandaceae", "Magnoliaceae", 
                                                                               "Nyssaceae", "Oleaceae", "Pinaceae", "Rosaceae",
                                                                               "Salicaceae", "Sapindaceae", "Ulmaceae", 
                                                                               "Abies", "Acer", "Betula", "Carya", "Fagus", 
                                                                               "Fraxinus", "Liquidambar", "Liriodendron", 
                                                                               "Nyssa", "Picea", "Pinus", "Populus",
                                                                               "Prunus", "Quercus", "Ulmus",
                                                                               "Arbuscular mycorrhizal", "Ectomycorrhizal")))
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
             size = 1,
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
  xlab("Values") + ylab("Plant group") +
  facet_grid(. ~ metric, scales = "free", labeller = my_labeller) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(-0.15, 0.15), breaks = c(-0.13, 0.0, 0.13)),
                               scale_x_continuous(limits = c(0, 0.4), breaks = c(0, 0.20, 0.40)),
                               scale_x_continuous(limits = c(5, 40), breaks = c(5, 20, 35)))) +
  th

p 

jpeg(paste0(root_path, "/figures/Manuscript_black/satellite_occurrence-models.jpeg"), 
     width = 210, 
     height = 240, 
     units = "mm", 
     res = 600, 
     bg = "white",
     #bg = "transparent", 
     pointsize = 14)
p
dev.off()


png(paste0(root_path, "/figures/Manuscript_white/satellite_occurrence-models.png"), 
    width = 220, 
    height = 240, 
    units = "mm", 
    res = 600, 
    bg = "transparent", 
    pointsize = 14)
p
dev.off()

