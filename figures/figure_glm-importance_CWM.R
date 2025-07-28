#'------------------------------------------------------------------------------
#' @title Visualize the importance of variables on GLM models of plant traits (Fig. S21)
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(bit64)
library(data.table)
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
colnames(data)[4:14] <- c("Axis 1", "Axis 2", "Axis 3", "MAT", "TS", "TAR", "MAP", "PS", "ELE", "SLO", "CTI")

# Estimate relative importance
total <- rowSums(data[,.SD, .SDcols = 4:ncol(data)])
data[, 4:ncol(data)] <- data[,.SD, .SDcols = 4:ncol(data)]/total

# Change names
data[plant == "HTq", name := "Tree height"]
data[plant == "SLq", name := "Slenderness"]
data[plant == "WD", name := "Wood density"]
data[plant == "shade", name := "Shade tolerance index"]
data[plant == "drought", name := "Drought tolerance index"]
data[plant == "waterlogging", name := "Waterlogging tolerance index"]
data[plant == "LMA", name := "Leaf mass per area"]
data[plant == "N", name := "Nitrogen"]

#-------------------------------------------------------------------------------
# Plot details

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.5, 4, 0.5, 0.5), "mm"),
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
        plot.margin = unit(c(0.5, 4, 0.5, 0.5), "mm"),
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
                                            "TS", "TAR", "MAP", "PS", "ELE", "SLO", "CTI"),
                           variable.name = "metric",
                           value.name = "importance")

summary_importance <- summary_importance[dataset == "training"]

summary_importance$metric <- as.factor(summary_importance$metric)
summary_importance$metric <- factor(summary_importance$metric, 
                                    levels = rev(c("Axis 1", "Axis 2", "Axis 3", "MAT", 
                                                   "TS", "TAR", "MAP", "PS", "ELE", "SLO", "CTI")))

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
names <- c(rev(c("Tree height", 
                 "Slenderness", 
                 "Shade tolerance index", 
                 "Drought tolerance index",
                 "Waterlogging tolerance index",
                 "Leaf mass per area",
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


p <- ggplot(summary_importance, aes(x = mean, 
                                    y = name,
                                    fill = metric)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_vline(xintercept = 0.5, colour = "grey", linetype = "dotted") +
  scale_fill_viridis(name = "Relative importance", 
                     option = "mako",
                     direction = -1,
                     discrete = TRUE,
                     breaks = c("Axis 1", "Axis 2", "Axis 3", "MAT", 
                                "TS", "TAR", "MAP", "PS", "ELE", "SLO", "CTI")) +
  theme(legend.position="top", legend.box = "horizontal") +
  guides(fill = guide_legend(nrow = 1)) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  xlab("Relative importance") + ylab("Community-weighted mean of plant traits") +
  th
p

jpeg(paste0(root_path, "/figures/Manuscript_black/glm-models_CWM.jpg"), 
     width = 270, 
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

