#'------------------------------------------------------------------------------
#' @title Visualize performance of community-weighted mean plant traits (Fig S18)
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

root_path <- "/media/antonio/Work/MMDB"

# Model performance results 
climate_topography <- fread(paste0(root_path, "/Models/03_model_CWM/func_climate-topography_model-CWM.csv"))
axes <- fread(paste0(root_path, "/Models/03_model_CWM/func_axes_model-CWM.csv"))
altogether <- fread(paste0(root_path, "/Models/03_model_CWM/func_altogether_model-CWM.csv"))

data <- rbindlist(list(climate_topography = climate_topography,
                       axes = axes,
                       altogether = altogether),
                  fill = TRUE,
                  use.names = TRUE,
                  idcol = "Model")
data <- data[, .SD, .SDcols = c("Model", colnames(altogether))]

# Change of model
data$Model <- as.factor(data$Model)
data$Model <- factor(data$Model,
                     levels = c("climate_topography", 
                                "axes", 
                                "altogether"),
                     labels = c("Climate and topography", 
                                "MDS axes", 
                                "MDS axes, climate, and topography"))

# Change of names
data <- data[order(plant),]
data[plant == "HTq", name := "Tree height"]
data[plant == "SLq", name := "Slenderness"]
data[plant == "WD", name := "Wood density"]
data[plant == "shade", name := "Shade tolerance index"]
data[plant == "drought", name := "Drought tolerance index"]
data[plant == "waterlogging", name := "Waterlogging tolerance index"]
data[plant == "LMA", name := "Leaf mass per area*"]
data[plant == "N", name := "Nitrogen"]

data$name <- as.factor(data$name)
data$name <- factor(data$name, 
                      levels = rev(c("Tree height", 
                                     "Slenderness", 
                                     "Shade tolerance index", 
                                     "Drought tolerance index",
                                     "Waterlogging tolerance index",
                                     "Leaf mass per area*",
                                     "Nitrogen",
                                     "Wood density")))

# Scale LMA
data[plant == "LMA", "RMSE"] <- data[plant == "LMA", "RMSE"]/10
data[plant == "LMA", "BIAS"] <- data[plant == "LMA", "BIAS"]/10

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
                            id.vars = c("Model", "dataset", "name"),
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
summary_all <- summary_performance[, .(mean = mean(performance), sd = sd(performance)),
                                   by = c("Model", "dataset", "metric")]

# Confidence
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
                                               lower = get_ci(performance, "lower"),
                                               upper = get_ci(performance, "upper")),
                                           by = c("Model", "dataset", "metric", "name")]
# Plot labels
my_labeller <- as_labeller(c(R2= "italic(R)^2", BIAS ="BIAS", RMSE = "RMSE", perRMSE = "\"%\"*RMSE"),
                           default = label_parsed)

# Create rects
rects <- data.frame(
  ymin = seq(0.5, length(unique(summary_performance$name)), by = 2),
  ymax = seq(1.5, length(unique(summary_performance$name)) + 1, by = 2),
  xmin = -Inf,
  xmax = Inf
)

selection <- "Testing"
pd <- position_dodge2(width = 0.6, preserve = "single")
p <- ggplot(summary_performance[dataset == selection,]) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +  
  geom_linerange(aes(xmin = lower,
                     xmax = upper,
                     y = name,
                     colour = Model,
                     group = Model),
                 position = pd,
                 alpha = 0.75) +
  geom_point(aes(x = mean,
                 y = name,
                 shape = Model,
                 fill = Model),
             size = 2,
             # colour = "grey",
             #shape = 21,
             position = pd,
             alpha = 0.75) +
  geom_vline(data = summary_all[dataset == selection,],
             aes(xintercept = mean, 
                 colour = Model,
                 linetype = Model)) +
  scale_shape_manual("Model", values = c(25, 21, 24)) +
  scale_colour_manual("Model", values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_fill_manual("Model", values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_linetype_manual("Model", values = c("dotted", "dotted", "dotted")) +
  xlab("Values") + ylab("Community-weighted mean of plant traits") +
  facet_grid(. ~ metric, scales = "free", labeller = my_labeller) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(-0.05, 0.05), breaks = c(-0.04, 0, 0.04)),
                               scale_x_continuous(limits = c(0, 10), breaks = c(0, 5, 10)),
                               scale_x_continuous(limits = c(10, 21), breaks = c(10, 15, 20)))) +
  th

p

jpeg(paste0(root_path, "/figures/Manuscript_black/FIA-models_CWM.jpg"), 
     width = 210, 
     height = 150, 
     units = "mm", 
     res = 600, 
     bg = "white",
     #bg = "transparent", 
     pointsize = 14)
p
dev.off()

#-----------------------------------------------------
# Testing

selection <- "Testing"
pd <- position_dodge2(width = 0.6, preserve = "single")
p <- ggplot(summary_performance[dataset == selection,]) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +  
  geom_linerange(aes(xmin = lower,
                     xmax = upper,
                     y = name,
                     colour = Model,
                     group = Model),
                 position = pd,
                 alpha = 0.75) +
  geom_point(aes(x = mean,
                 y = name,
                 shape = Model,
                 fill = Model),
             size = 2,
             # colour = "grey",
             #shape = 21,
             position = pd,
             alpha = 0.75) +
  geom_vline(data = summary_all[dataset == selection,],
             aes(xintercept = mean, 
                 colour = Model,
                 linetype = Model)) +
  scale_shape_manual("Dataset", values = c(25, 21, 24)) +
  scale_colour_manual("Dataset", values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_fill_manual("Dataset", values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_linetype_manual("Dataset", values = c("dotted", "dotted", "dotted")) +
  xlab("Values") + ylab("Community-weighted mean of plant traits") +
  facet_grid(. ~ metric, scales = "free", labeller = my_labeller) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(-0.045, 0.045), breaks = c(-0.04, 0, 0.04)),
                               scale_x_continuous(limits = c(0, 10), breaks = c(0, 5, 10)),
                               scale_x_continuous(limits = c(10, 21), breaks = c(10, 15, 20)))) +
  th

p

jpeg(paste0(root_path, "/figures/Manuscript_black/FIA-models_CWM.jpg"), 
     width = 210, 
     height = 150, 
     units = "mm", 
     res = 600, 
     bg = "white",
     #bg = "transparent", 
     pointsize = 14)
p
dev.off()

# For transparencies
png(paste0(root_path, "/figures/Manuscript_white/FIA-models_CWM.png"), 
    width = 210, 
    height = 150, 
    units = "mm", 
    res = 600, 
    bg = "transparent", 
    pointsize = 14)
p
dev.off()
