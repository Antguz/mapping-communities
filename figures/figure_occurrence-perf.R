#'------------------------------------------------------------------------------
#' @title Codes to visualize the performance of probability of occurrence
#'------------------------------------------------------------------------------

#' @description Rerun .jpeg of observed and predicted probabilities

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(parallel)
library(pbmcapply)
library(ggplot2)
library(ggpubr)
library(bit64)
library(ggdark)
library(ggpointdensity)
library(viridis)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

#root_path <- "E:/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping/08_occurence_performance"
out_path <- paste0(root_path, "/figures")

#'------------------------------------------------------------------------------
#' @Plot aesthetics
#-------------------------------------------------------------------------------

th_black <- theme_bw(base_size = 11) + 
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

point_size <- 0.5
point_shape <- 16
point_alpha <- 0.9

#'------------------------------------------------------------------------------
#' @Plot-function
#-------------------------------------------------------------------------------

performance_per_plant <- function(X,
                                  file,
                                  sensor,
                                  dataset,
                                  type,
                                  limits_range = c(-0.01, 0.51),
                                  breaks = c(0, 0.25, 0.5),
                                  out_path) {
  
  # Plot features
  if(type == "taxo") {
    title <- "Taxonomic"
  } else if(type == "phylo-log") {
    title <- "Phylogenetic"
  } else if(type == "func") {
    title <- "Functional"
  }
  
  if(dataset == "training") {
    dat <- "Training"
  } else if(dataset == "testing") {
    dat <- "Testing"
  }
  
  # Make frame
  frame <- data.table(obs = file$obs[, .SD, .SDcols = X],
                      pred = file$pred[, .SD, .SDcols = X])
  colnames(frame) <- c("observed", "predicted")
  perf <- file$perf[X, ]
  
  p <- ggplot(frame) +
    geom_pointdensity(data = frame,
                      aes(x = predicted, y = observed),
                      size = point_size,
                      shape = point_shape,
                      alpha = point_alpha,
                      adjust = .1,
                      show.legend = FALSE) +
    geom_abline(aes(intercept = 0, slope = 1),
                colour = "grey50", linetype = "dotted") +
    geom_abline(data = perf,
                aes(intercept = intercept, slope = slope),
                colour = line_col) +
    ylab("Observed probability") + 
    xlab("Predicted probability") +
    geom_text(data = perf,
              aes(label = paste("R2 = ", round(mean(perf$R2), 2), "\n",
                                "MAE = ", round(mean(perf$MAE), 2), "\n",
                                "RMSE = ", round(mean(perf$RMSE), 2), "\n",
                                "%RMSE = ", round(mean(perf$perRMSE), 2)),
                  x = Inf, 
                  y = -Inf, 
                  hjust = 1.1, 
                  vjust = -0.2),
              colour = "grey25",
              size = 2.5) +
    labs(title = paste0(sensor, "  -  ", title),
         subtitle = perf$plants[1],
         caption = paste0(dat, " dataset")) + 
    scale_colour_viridis() +
    th
  
  p
  
  if(is.null(limits_range) == FALSE) {
    p <- p + coord_cartesian(xlim = limits_range, ylim = limits_range, expand = FALSE) +
      scale_x_continuous(breaks = breaks, labels = as.character(breaks)) +
      scale_y_continuous(breaks = breaks, labels = as.character(breaks))
  }
  
  # Export
  export_name <- paste0(out_path, "/", 
                        sensor, "/",
                        type, "/",
                        dataset, "/",
                        sensor, "_",
                        type, "_",
                        dataset, "_",
                        perf$plants[1], "_prob.jpeg")
  
  ggsave(export_name, p, 
         width = 100,
         height = 100,
         units = c("mm"),
         dpi = 300)
  
}

#'------------------------------------------------------------------------------
#' @Read information and apply plot function
#-------------------------------------------------------------------------------

sensor <- "EMIT"
dataset <- "training"
type <- "phylo-log"
file <- readRDS(paste0(root_path, "/", sensor, "_", type, "_", dataset, "_probabilities.rds"))

pbmclapply(X = 1:nrow(file$perf),
           FUN = performance_per_plant,
           file = file,
           sensor = sensor,
           dataset = dataset,
           type = type,
           limits_range = c(-0.01, 1.01),
           breaks = c(0, 0.5, 1.0),
           out_path = out_path,
           mc.cores = 25)
