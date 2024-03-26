#'------------------------------------------------------------------------------
#' @title Codes to visualize the stressplot from NMDS
#'------------------------------------------------------------------------------

#' @description Export a stressplot from NMDS. Our NMDS are huge, so this task 
#' most be completed in MSI.

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)
library(bit64)
library(ggdark)
library(viridis)
library(scales)
library(ggpubr)
library(vegan)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Working directory

#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"

#'------------------------------------------------------------------------------
#' Read information

#type <- "taxo"
#NMDS_path <- paste0(root_path, "/05_NMDS/NMDS_", type, ".rds")
#text <- "Taxonomic"

#-------------------------------------------------------------------------------
# Plot details

th_black <- theme_bw(base_size = 11, base_family = "Helvetica") + 
  theme(legend.position = "top",
        legend.title.align = 0.5,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"))

th_trans <- dark_theme_bw(base_size = 11) +
  theme(legend.position = "top",
        legend.title.align = 0.5,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        #legend.position = c(0.90, 0.90), 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth = 0))

th <- th_black

#'------------------------------------------------------------------------------
#' Function

stress_figure <- function(object, text, out_path) {
  
  th_black <- theme_bw(base_size = 11, base_family = "Helvetica") + 
    theme(legend.position = "top",
          legend.title.align = 0.5,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill= "transparent"),
          panel.background = element_rect(fill= "transparent"),
          rect = element_rect(fill = "transparent"))
  
  th_trans <- dark_theme_bw(base_size = 11) +
    theme(legend.position = "top",
          legend.title.align = 0.5,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          #legend.position = c(0.90, 0.90), 
          legend.background = element_rect(fill= "transparent"),
          panel.background = element_rect(fill = "transparent"),
          rect = element_rect(fill = "transparent"),
          legend.title = element_text(colour="white"),
          legend.text = element_text(color = "white"),
          axis.text = element_text(color="white"),
          plot.background = element_rect(colour = "transparent",
                                         fill = "transparent",
                                         linewidth = 0))
  
  th <- th_black
  
  # Read NMDS
  #object <- readRDS(NMDS_path)
  istart <- object$istart
  ndis <- object$ndis
  iregn <- object$iregn
  ngrp <- object$ngrp
  stress <- object$stress
  dist <- object$dist
  diss <- object$diss
  dhat <- object$dhat
  #rm(object)
  gc()
  
  # Collect values for 'linear fit'
  ralscal <- 0
  k <- NA
  
  # Monotone line
  if (iregn != 2) {
    ist <- c(istart, ndis + 1)
    
    for(j in 1:ngrp) {
      k <- seq(ist[j], ist[j+1]-1)
      ralscal <- ralscal + cor(dist[k], dhat[k])^2
    }
  }
  
  # Stress as R2
  rstress <- 1 - stress^2
  Rst <- format(rstress, digits = 2)
  Ral <- format(ralscal, digits = 2)
  
  a <- sample(1:ndis, 100000)
  b <- seq(from = 1, to = ndis, length.out = 10000000)
  
  gc()
  
  # Resample line
  
  p <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey") +
    geom_hex(aes(x = diss,#[a], 
                 y = dist),#[a]),
             bins = 100) +
    geom_line(aes(x = diss[k][b], 
                  y = dhat[k][b]), 
              colour = "grey75") +
    geom_abline(slope = 1, 
                intercept = 0, 
                linetype = "dotted", 
                linewidth = 0.5,
                colour = "grey") +
    xlab("Observed dissimilarity") + ylab("Ordination distance") +
    coord_cartesian(xlim = c(0, 1), 
                    ylim = c(0, 1), 
                    expand = FALSE) +
    scale_x_continuous(n.breaks = 3) +
    scale_y_continuous(n.breaks = 3) +
    scale_fill_viridis_c(name = "Count", 
                         option = "mako",
                         trans = "log10", 
                         #limits = c(1, 10000000),  
                         #breaks = c(1, 3162.278, 10000000),
                         labels = trans_format("log10", math_format(10^.x))) +
    guides(fill = guide_colourbar(barwidth = 15, 
                                  barheight = 0.7, 
                                  title.position = "bottom")) +
    geom_text(aes(label = text,
                  x = -Inf, 
                  y = Inf, 
                  hjust = -0.1, 
                  vjust = 1.50),
              colour = "black",
              size = 3,
              family = "Helvetica") +
    geom_text(aes(label = paste0("Non-metric fit, R^2 = ", Rst, "\n",
                                 "Linear fit, R^2 = ", Ral),
                  x = -Inf, 
                  y = Inf, 
                  hjust = -0.1, 
                  vjust = 2.50),
              colour = "grey25",
              size = 2,
              family = "Helvetica") +
    th +
    theme(panel.grid = element_blank())

  ggsave(paste0("figure_stressplot_", type, ".jpeg"),
         plot = p,
         path = out_path,
         width = 90, 
         height = 90, 
         units = "mm", 
         dpi = 600, 
         bg = "white")
  
}
