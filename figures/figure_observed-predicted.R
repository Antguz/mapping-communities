#'------------------------------------------------------------------------------
#' @title Codes to visualize the NDMS_predicted_mean values
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(ggplot2)
library(ggpubr)
library(bit64)
library(ggdark)
library(ggpointdensity)
library(rcartocolor)
library(viridis)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Working directory

#root_path <- "E:/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
sensor <- "EMIT"
folder <- paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor)

#'------------------------------------------------------------------------------
#' Read information

taxo_performance <- fread(paste0(folder, "/", sensor, "_", "taxo", "_performance.csv"))
taxo_performance$dimension <- "Taxonomic"
phylo_performance <- fread(paste0(folder, "/", sensor, "_", "phylo-log", "_performance.csv"))
phylo_performance$dimension <- "Phylogenetic"
func_performance <- fread(paste0(folder, "/", sensor, "_", "func", "_performance.csv"))
func_performance$dimension <- "Functional"

dataset <- "training"

if(dataset == "training") {
  
  performance <- rbind(taxo_performance, phylo_performance, func_performance)
  performance$dimension <- as.factor(performance$dimension)
  performance$dimension <- factor(performance$dimension, levels = c("Taxonomic", 
                                                                    "Phylogenetic", 
                                                                    "Functional"))
  performance <- performance[type == "Training", ]
  
  taxo <- fread(paste0(folder, "/", sensor, "_", "taxo", "_predicted_training.csv"))
  taxo$dimension <- "Taxonomic"
  phylo <- fread(paste0(folder, "/", sensor, "_", "phylo-log", "_predicted_training.csv"))
  phylo$dimension <- "Phylogenetic"
  func <- fread(paste0(folder, "/", sensor, "_", "func", "_predicted_training.csv"))
  func$dimension <- "Functional"
  
  predicted <- rbind(taxo, phylo, func)
  predicted$dimension <- as.factor(predicted$dimension)
  predicted$dimension <- factor(predicted$dimension, levels = c("Taxonomic", 
                                                                "Phylogenetic", 
                                                                "Functional"))
  
} else if(dataset == "testing") {
  
  performance <- rbind(taxo_performance, phylo_performance, func_performance)
  performance$dimension <- as.factor(performance$dimension)
  performance$dimension <- factor(performance$dimension, levels = c("Taxonomic", 
                                                                    "Phylogenetic", 
                                                                    "Functional"))
  performance <- performance[type == "Testing", ]
  
  taxo <- fread(paste0(folder, "/", sensor, "_", "taxo", "_predicted_testing.csv"))
  taxo$dimension <- "Taxonomic"
  phylo <- fread(paste0(folder, "/", sensor, "_", "phylo-log", "_predicted_testing.csv"))
  phylo$dimension <- "Phylogenetic"
  func <- fread(paste0(folder, "/", sensor, "_", "func", "_predicted_testing.csv"))
  func$dimension <- "Functional"
  
  predicted <- rbind(taxo, phylo, func)
  predicted$dimension <- as.factor(predicted$dimension)
  predicted$dimension <- factor(predicted$dimension, levels = c("Taxonomic", 
                                                                "Phylogenetic", 
                                                                "Functional"))
  
}

#-------------------------------------------------------------------------------
# Plot details

col <- c("#fc8d62", "#66c2a5", "#8da0cb")

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = c(0.93, 0.90),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"))

th_trans <- dark_theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position.inside = c(0.90, 0.90), 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth=0))

#th <- th_trans
th <- th_black
#line_col <- "white"
line_col <- "black"

point_size <- 2
point_shape <- 20
point_alpha <- 0.9

#'------------------------------------------------------------------------------
#' Plot function

#Test
#performance <- performance[type == "Testing" & dimension == "Phylogenetic" & Axes == "Axis1",]
#predicted <- predicted[type == "Testing" & dimension == "Phylogenetic" & Axes == "Axis1",]


plot_function <- function(performance, 
                          predicted, 
                          limits_range = c(-1.5, 1.5),
                          breaks = c(-0.5, 0, 0.5),
                          axis = 1) {
  
  print(paste0("Minimun: ", round(min(range(predicted$observed), range(predicted$predicted)), 2),
               " Maximun: ", round(max(range(predicted$observed), range(predicted$predicted)), 2)))
  
  #xval <- limits_range[2] * 0.95
  #yval <- limits_range[1] * 0.95
  
  
  p <- ggplot() +
    geom_abline(aes(intercept = 0, slope = 1),
                colour = line_col, linetype = "dotted") +
    geom_errorbarh(data = predicted,
                   aes(xmin = predicted - predicted_sd, 
                       xmax = predicted + predicted_sd,
                       y = observed),
                   colour = "grey") +
    geom_pointdensity(data = predicted,
                      aes(x = predicted, y = observed),
                      size = point_size,
                      shape = point_shape,
                      alpha = point_alpha) +
    geom_abline(data = performance,
                aes(intercept = intercept, slope = slope),
                colour = "grey", alpha = 0.5) +
    geom_abline(data = performance,
                aes(intercept = mean(intercept), slope = mean(slope)),
                colour = line_col) +
    ylab(" ") + 
    xlab(" ") +
    #ylab("Observed") + 
    #xlab("Predicted") +
    geom_text(data = performance,
              aes(label = paste("R2 = ", round(mean(performance$R2), 2), "\n",
                                "MAE = ", round(mean(performance$MAE), 2), "\n",
                                "RMSE = ", round(mean(performance$RMSE), 2), "\n",
                                "%RMSE = ", round(mean(performance$perRMSE), 2)),
                  x = Inf, 
                  y = -Inf, 
                  hjust = 1.1, 
                  vjust = -0.2),
              colour = "grey25",
              size = 2.5) +
    #labs(subtitle = paste0("PCoA", axis)) + 
    th
  
  if(is.null(limits_range) == FALSE) {
    p <- p + coord_cartesian(xlim = limits_range, ylim = limits_range, expand = FALSE) +
      scale_x_continuous(breaks = breaks, labels = as.character(breaks)) +
      scale_y_continuous(breaks = breaks, labels = as.character(breaks))
  }
  
  if(axis == 1) {
    p <- p + scale_colour_carto_c(type = "sequential", 
                                  palette = "Peach", 
                                  direction = 1,
                                  guide = "none")
    
  } else if(axis == 2) {
    p <- p + scale_colour_carto_c(type = "sequential", 
                                  palette = "Mint", 
                                  direction = 1,
                                  guide = "none")
    
  } else if(axis == 3) {
    p <- p + scale_colour_carto_c(type = "sequential", 
                                  palette = "Teal", 
                                  direction = 1,
                                  guide = "none")
    
  }
  
  return(p)
  
}


taxo1 <- plot_function(performance[dimension == "Taxonomic" & Axes == "Axis1",], 
                       predicted[dimension == "Taxonomic" & Axes == "Axis1",], 
                       limits_range = c(-0.55, 0.55),
                       breaks = c(-0.5, 0, 0.5),
                       axis = 1)

phylo1 <- plot_function(performance[dimension == "Phylogenetic" & Axes == "Axis1",], 
                        predicted[dimension == "Phylogenetic" & Axes == "Axis1",], 
                        limits_range = c(-0.55, 0.55),
                        breaks = c(-0.5, 0, 0.5),
                        axis = 1)

func1 <- plot_function(performance[dimension == "Functional" & Axes == "Axis1",], 
                       predicted[dimension == "Functional" & Axes == "Axis1",], 
                       limits_range = c(-0.55, 0.55),
                       breaks = c(-0.5, 0, 0.5),
                       axis = 1)

taxo2 <- plot_function(performance[dimension == "Taxonomic" & Axes == "Axis2",], 
                       predicted[dimension == "Taxonomic" & Axes == "Axis2",], 
                       limits_range = c(-0.55, 0.55),
                       breaks = c(-0.5, 0, 0.5),
                       axis = 2)

phylo2 <- plot_function(performance[dimension == "Phylogenetic" & Axes == "Axis2",], 
                        predicted[dimension == "Phylogenetic" & Axes == "Axis2",], 
                        limits_range = c(-0.55, 0.55),
                        breaks = c(-0.5, 0, 0.5),
                        axis = 2)

func2 <- plot_function(performance[dimension == "Functional" & Axes == "Axis2",], 
                       predicted[dimension == "Functional" & Axes == "Axis2",], 
                       limits_range = c(-0.55, 0.55),
                       breaks = c(-0.5, 0, 0.5),
                       axis = 2)

taxo3 <- plot_function(performance[dimension == "Taxonomic" & Axes == "Axis3",], 
                       predicted[dimension == "Taxonomic" & Axes == "Axis3",], 
                       limits_range = c(-0.55, 0.55),
                       breaks = c(-0.5, 0, 0.5),
                       axis = 3)

phylo3 <- plot_function(performance[dimension == "Phylogenetic" & Axes == "Axis3",], 
                        predicted[dimension == "Phylogenetic" & Axes == "Axis3",], 
                        limits_range = c(-0.55, 0.55),
                        breaks = c(-0.5, 0, 0.5),
                        axis = 3)

func3 <- plot_function(performance[dimension == "Functional" & Axes == "Axis3",], 
                       predicted[dimension == "Functional" & Axes == "Axis3",], 
                       limits_range = c(-0.55, 0.55),
                       breaks = c(-0.5, 0, 0.5),
                       axis = 3)

# Merge plots
axis1 <- ggarrange(taxo1, phylo1, func1, 
                   ncol=3, nrow=1,
                   labels = c("A", "B", "C"),
                   font.label =  list(size = 9, face = "plain", color ="grey25"),
                   label.x = 0.08,
                   label.y = 1.05)

axis2 <- ggarrange(taxo2, phylo2, func2, ncol=3, nrow=1,
                   labels = c("D", "E", "F"),
                   font.label =  list(size = 9, face = "plain", color ="grey25"),
                   label.x = 0.08,
                   label.y = 1.05)

axis3 <- ggarrange(taxo3, phylo3, func3, ncol=3, nrow=1,
                   labels = c("G", "H", "I"),
                   font.label =  list(size = 9, face = "plain", color ="grey25"),
                   label.x = 0.08,
                   label.y = 1.05)

# Export
jpeg(paste0(root_path, "/figures/", sensor, "_predicted_", dataset, ".jpeg"), 
     width = 210, 
     height = 190, 
     units = "mm", 
     res = 300, 
     bg = "white")

annotate_figure(ggarrange(axis1, axis2, axis3, ncol= 1, nrow=3),
                top = text_grob("Taxonomic                        Phylogenetic                        Functional", face = "bold", color = "black"),
                left = "Observed",
                right = text_grob("Axis 1                        Axis 2                        Axis 3", face = "bold", color = "black", rot = -90),
                bottom = "Predicted")

dev.off()

