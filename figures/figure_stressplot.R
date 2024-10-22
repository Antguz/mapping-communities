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
library(rTLS)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Working directory

#root_path <- "/media/antonio/Work/MMDB"

#'------------------------------------------------------------------------------
#' Read information

#type <- "taxo"
#object <- readRDS(paste0(root_path, "/04_NMDS/NMDS_", type, "_total.rds"))
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
  
th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(1, 10, 1, 1), "pt"),
        legend.position= c("top"), 
        strip.background = element_rect(color="black", 
                                        fill="black", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "white"))
  
  th_trans <- dark_theme_bw(base_size = 11) +
    theme(legend.position = "top",
          legend.title.align = 0.5,
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"),
          #legend.position = "top", 
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
  rm(object)
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
  
  # Line
  b <- seq(from = 1, to = ndis, length.out = 10000)
  x_line <- diss[k][b]
  y_line <- dhat[k][b]
  
  rm(list = c("k", "b", "dhat"))
  gc()
  
  # Stress as R2
  rstress <- 1 - stress^2
  Rst <- format(rstress, digits = 2)
  Ral <- format(ralscal, digits = 2)
  
  # Vox
  frame <- data.table(x = diss, y = dist, z = 1)
  rm(list = c("diss", "dist"))
  gc()
  vox <- voxels(frame, edge_length = c(0.01, 0.01, 0.01), threads = 16, obj.voxels = FALSE)
  plot(vox$X-0.00025, vox$Y-0.00025)
  rm(list = c("frame"))
  gc()
  
  p <- ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey") +
    geom_raster(data = vox, aes(x = X, y = Y, fill = N)) +
    #geom_hex(aes(x = diss[a],#[a], 
    #             y = dist[a]),#[a]),
    #         bins = 100) +
    geom_line(aes(x = x_line, 
                  y = y_line), 
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
                         limits = c(1, 1995262),  
                         breaks = c(1, 1000, 1000000),
                         labels = trans_format("log10", math_format(10^.x))) +
    guides(fill = guide_colourbar(barwidth = 15, 
                                  barheight = 0.7, 
                                  title.position = "bottom")) +
    geom_text(aes(label = paste0("Non-metric fit, R2 = ", Rst, "\n",
                                 "Linear fit, R2 = ", Ral),
                  x = -Inf, 
                  y = Inf, 
                  hjust = -0.1, 
                  vjust = 1.50),
              colour = "grey25",
              size = 2,
              family = "Helvetica") +
    th +
    theme(panel.grid = element_blank())

  jpeg(paste0(out_path, "/figure_stressplot_", type, ".jpeg"), 
      width = 90, 
      height = 90, 
      pointsize = 12, 
      units = "mm", 
      res = 600, 
      bg = "white")
  
  p
  
  dev.off()
  
}
