#'------------------------------------------------------------------------------
#' @title Visualize performance of occurrence models of forest types (Fig S16)
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
data <- fread(paste0(root_path, "/Models/03_models_fortyp/taxo_model-FORTYP.csv"))
data$MDS <- "non-rotated"
data <- data[order(forty),]
data$forty <- as.factor(data$forty)
data$forty <- factor(data$forty, levels = rev(as.character(unique(data$forty))))
data$TSS <- data$Sens + data$Spec - 1


data1 <- fread(paste0(root_path, "/Models/OLD/03_models_fortyp/taxo_model-FORTYP.csv"))
data1$MDS <- "rotated"
data1 <- data1[order(forty),]
data1$forty <- as.factor(data1$forty)
data1$forty <- factor(data1$forty, levels = rev(as.character(unique(data1$forty))))
data1$TSS <- data1$Sens + data1$Spec - 1

data <- rbind(data, data1)

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
                            id.vars = c("MDS", "dataset", "forty"),
                            measure.vars = c("TSS", "Sens", "Spec", "ROC"),
                            variable.name = "metric",
                            value.name = "performance")

summary_performance$dataset <- as.factor(summary_performance$dataset)
summary_performance$dataset <- factor(summary_performance$dataset, 
                                      levels = c("training", "testing"),
                                      labels = c("Training", "Testing"))

summary_performance$metric <- as.factor(summary_performance$metric)
summary_performance$metric <- factor(summary_performance$metric, 
                                     levels = c("TSS", "Sens", "Spec", "ROC"),
                                     labels = c("TSS", "Sensitivity", "Specificity", "ROC"))

# Summary regardless of trait
summary_all <- summary_performance[, .(mean = mean(performance), sd = sd(performance)),
                                   by = c("MDS", "dataset", "metric")]

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
                                           by = c("MDS", "dataset", "metric", "forty")]
# Plot labels
my_labeller <- as_labeller(c(TSS = "TSS", 
                             Sensitivity = "Sensitivity", 
                             Specificity = "Specificity", 
                             ROC = "AUC[ROC]"),
                           default = label_parsed)

# Create rects
rects <- data.frame(
  ymin = seq(0.5, length(unique(summary_performance$forty)), by = 2),
  ymax = seq(1.5, length(unique(summary_performance$forty)) + 1, by = 2),
  xmin = -Inf,
  xmax = Inf
)

p <- ggplot(summary_performance) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +  
  # geom_linerange(aes(xmin = lower,
  #                    xmax = upper,
  #                    y = forty,
  #                    colour = dataset),
  #                position = position_dodge(width = -0.4),
  #                alpha = 0.75) +
  geom_point(aes(x = mean, 
                 y = forty,
                 colour = interaction(dataset, MDS),
                 fill = interaction(dataset, MDS),
                 shape = interaction(dataset, MDS)),
             size = 1.5,
             #colour = "grey",
             alpha = 0.75,
             position = position_dodge(width = -0.4)) +
  # geom_vline(data = summary_all[dataset == "Testing",],
  #             aes(xintercept = mean,
  #                 colour = interaction(dataset, MDS),
  #                 linetype = interaction(dataset, MDS))) +
  scale_shape_manual("Dataset", values = c(24, 25, 24, 25)) +
  scale_colour_manual("Dataset", values = c("#e66101", "#5e3c99", "#e66101", "#5e3c99")) +
  scale_fill_manual("Dataset", values = c("white", "white", "#e66101", "#5e3c99")) +
  #scale_linetype_manual("Dataset", values = c("dotted", "dotted", "dashed", "dashed")) +
  xlab("Values") + ylab("Forest type") +
  facet_grid(. ~ metric, scales = "free", labeller = my_labeller) +
  facetted_pos_scales(x = list(scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1.0)),
                               scale_x_continuous(limits = c(0.6, 1.0), breaks = c(0.6, 0.8, 1.0)),
                               scale_x_continuous(limits = c(0.6, 1.0), breaks = c(0.6, 0.8, 1.0)),
                               scale_x_continuous(limits = c(0.6, 1.0), breaks = c(0.6, 0.80, 1.0)))) +
  th

p

jpeg(paste0(root_path, "/figures/Manuscript_black/FIA-models_fortype.jpeg"), 
     width = 210, 
     height = 240, 
     units = "mm", 
     res = 600, 
     bg = "white",
     #bg = "transparent", 
     pointsize = 14)

p

dev.off()


png(paste0(root_path, "/figures/Manuscript_white/FIA-models-fortype.png"), 
    width = 220, 
    height = 240, 
    units = "mm", 
    res = 600, 
    bg = "transparent", 
    pointsize = 14)

p

dev.off()

