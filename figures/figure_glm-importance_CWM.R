#'------------------------------------------------------------------------------
#' @title Visualize the importance of variables on GLM models of occurrence
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

# Categorical 
data <- fread(paste0(root_path, "/Models/03_model_CWM/func_model-CWM.csv"))
data <- data[, .SD, .SDcols = c(1:3, 23:33)]
colnames(data)[4:14] <- gsub("imp_","", colnames(data)[4:14])
colnames(data)[4:6] <- c("Axis 1", "Axis 2", "Axis 3")

data[plant == "HTq", name := "Tree height (max)"]
data[plant == "SLq", name := "Slenderness"]
data[plant == "WD", name := "Wood density"]
data[plant == "shade", name := "Shade tolerance index"]
data[plant == "drought", name := "Drought tolerance index"]
data[plant == "waterlogging", name := "Waterlogging tolerance index"]
data[plant == "LMA", name := "Leaf Mass Area"]
data[plant == "N", name := "Nitrogen"]

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
summary_importance <- melt(data, 
                           id.vars = c("dataset", "name"),
                           measure.vars = c("Axis 1", "Axis 2", "Axis 3", "MAT", 
                                            "TS", "TAR", "MAP", "PS", "DEM", "SL","CTI"),
                           variable.name = "metric",
                           value.name = "importance")

summary_importance <- summary_importance[dataset == "training"]

summary_importance$metric <- as.factor(summary_importance$metric)
summary_importance$metric <- factor(summary_importance$metric, 
                                    levels = c("Axis 1", "Axis 2", "Axis 3", "MAT", 
                                               "TS", "TAR", "MAP", "PS", "DEM", "SL","CTI"))

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

# Summary of all
summary_all <- summary_importance[, .(mean = mean(importance),
                                      lower = get_ci(importance, "lower"),
                                      upper = get_ci(importance, "upper")),
                                  by = c("metric")]
summary_all <- cbind(name = "Average", summary_all)
summary_all <- summary_all[, c(2,1,3:5)]

# Summary of results
summary_importance <- summary_importance[, .(mean = mean(importance),
                                             lower = get_ci(importance, "lower"),
                                             upper = get_ci(importance, "upper")),
                                         by = c("metric", "name")]

# Order forest types
summary_importance <- summary_importance[order(name)]
names <- c(rev(c("Tree height (max)", 
                 "Slenderness", 
                 "Shade tolerance index", 
                 "Drought tolerance index",
                 "Waterlogging tolerance index",
                 "Leaf Mass Area*",
                 "Nitrogen",
                 "Wood density")), "Average")
summary_importance <- rbind(summary_importance, summary_all)
summary_importance$name <- factor(summary_importance$name, 
                                   levels = names)

# Create rects
rects <- data.frame(
  ymin = seq(0.5, length(unique(summary_importance$name)), by = 2),
  ymax = seq(1.5, length(unique(summary_importance$name)) + 1, by = 2),
  xmin = -Inf,
  xmax = Inf
)

p <- ggplot(summary_importance) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_point(aes(x = metric, 
                 y = name,
                 size = mean,
                 fill = mean),
             shape = 21,
             colour = "grey",
             alpha = 0.75) +
  scale_size_continuous(" ",
                        limits = c(0, 180),
                        breaks = c(0, 60, 120, 180)) +
  scale_fill_viridis(name = "Variable of Importance", 
                     option = "turbo",
                     limits = c(0, 180),
                     breaks = c(0, 60, 120, 180),
                     oob = scales::squish) +
  xlab("Variables") + ylab("Community-weighted mean of plant traits") +
  th

p


jpeg(paste0(root_path, "/figures/Manuscript_black/glm-models_CWM.jpg"), 
     width = 210, 
     height = 150, 
     units = "mm", 
     res = 600, 
     bg = "white",
     #bg = "transparent", 
     pointsize = 14)
p
dev.off()


png(paste0(root_path, "/figures/Manuscript_white/glm-models_CWM.png"), 
    width = 210, 
    height = 150, 
    units = "mm", 
    res = 600, 
    bg = "transparent", 
    pointsize = 14)
p
dev.off()
