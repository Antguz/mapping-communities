#'------------------------------------------------------------------------------
#' @title Visualize the NDMS as a video
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(gg3D)
library(stringr)
library(bit64)
library(plotly)
library(plot3D)
library(data.table)
library(ggpubr)
library(Rcpp)
library(RcppArmadillo)
library(pbmcapply)
library(viridis)
library(scales)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

sourceCpp("src/probability_rcpp.cpp")
sourceCpp("src/cwm_rcpp.cpp")

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"

#'------------------------------------------------------------------------------
#' Read NMDS axes

# Communities
taxo_comm <- fread(paste0(root_path, "/Diversity/02_LMDS/MDS_scores_comm_taxo_total.csv"))
phylo_comm <- fread(paste0(root_path, "/Diversity/02_LMDS/MDS_scores_comm_phylo_total.csv"))
func_comm <- fread(paste0(root_path, "/Diversity/02_LMDS/MDS_scores_comm_func_total.csv"))

#-------------------------------------------------------------------------------
# Occurrence models
# 
# # Read occurrence models
# taxo_occurance <- fread(paste0(root_path, "/07_model_occurrence/taxo_model-occurrence.csv"))
# taxo_occurance <- taxo_occurance[dataset == "training", ]
# phylo_occurance <- fread(paste0(root_path, "/07_model_occurrence/phylo-log_model-occurrence.csv"))
# phylo_occurance <- phylo_occurance[dataset == "training", ]
# func_occurance <- fread(paste0(root_path, "/07_model_occurrence/func_model-occurrence.csv"))
# func_occurance <- func_occurance[dataset == "training", ]
# 
# # Get probabilities
# taxo_prob <- multi_probability(as.matrix(taxo_comm[,-1]),
#                                as.matrix(taxo_occurance[, c(14:17)]))
# taxo_prob <- as.data.table(taxo_prob)/10000
# colnames(taxo_prob) <- taxo_occurance$plant
# 
# phylo_prob <- multi_probability(as.matrix(phylo_comm[,-1]),
#                                 as.matrix(phylo_occurance[, c(14:17)]))
# phylo_prob <- as.data.table(phylo_prob)/10000
# colnames(phylo_prob) <- phylo_occurance$plant
# 
# func_prob <- multi_probability(as.matrix(func_comm[,-1]),
#                                as.matrix(func_occurance[, c(14:17)]))
# func_prob <- as.data.table(func_prob)/10000
# colnames(func_prob) <- func_occurance$plant
# 
# #-------------------------------------------------------------------------------
# # CWM models
# 
# # Read CWM models
# taxo_CWM <- fread(paste0(root_path, "/07_model_CWM/taxo_model-CWM.csv"))
# taxo_CWM <- taxo_CWM[dataset == "training", ]
# phylo_CWM <- fread(paste0(root_path, "/07_model_CWM/phylo-log_model-CWM.csv"))
# phylo_CWM <- phylo_CWM[dataset == "training", ]
# func_CWM <- fread(paste0(root_path, "/07_model_CWM/func_model-CWM.csv"))
# func_CWM <- func_CWM[dataset == "training", ]
# 
# # Get CWM estimates
# taxo_cwm <- multi_cwm(as.matrix(taxo_comm[,-1]),
#                       as.matrix(taxo_CWM[, c(10:13)]))
# taxo_cwm <- as.data.table(taxo_cwm)
# colnames(taxo_cwm) <- taxo_CWM$plant
# 
# phylo_cwm <- multi_cwm(as.matrix(phylo_comm[,-1]),
#                        as.matrix(phylo_CWM[, c(10:13)]))
# phylo_cwm <- as.data.table(phylo_cwm)
# colnames(phylo_cwm) <- phylo_CWM$plant
# 
# func_cwm <- multi_cwm(as.matrix(func_comm[,-1]),
#                       as.matrix(func_CWM[, c(10:13)]))
# func_cwm <- as.data.table(func_cwm)
# colnames(func_cwm) <- func_CWM$plant

#'------------------------------------------------------------------------------
#' @Plot-functions
#-------------------------------------------------------------------------------

# RGB color
RGB <- function(X, Y, Z, min_range = -0.50, max_range = 0.50) {
  
  if(is.null(min_range) == TRUE & is.null(max_range) == TRUE) {
    
    x_new <-  ((X - min(X)) / ((max(X)) - (min(X)))) * (255 - 0) + 0
    y_new <-  ((Y - min(Y)) / ((max(Y)) - (min(Y)))) * (255 - 0) + 0
    z_new <-  ((Z - min(Z)) / ((max(Z)) - (min(Z)))) * (255 - 0) + 0
    
  } else {
    
    x_new <-  ((X - min_range) / (max_range - min_range)) * (255 - 0) + 0
    y_new <-  ((Y - min_range) / (max_range - min_range)) * (255 - 0) + 0
    z_new <-  ((Z - min_range) / (max_range - min_range)) * (255 - 0) + 0
    
  }
  
  # Convert colour channels to hex
  hex_value <- paste0(
    "#",
    # Convert red
    str_pad(
      string = as.hexmode(round(x_new, 0)),
      width = 2, pad = "0", side = "left"
    ),
    # Convert green
    str_pad(
      string = as.hexmode(round(y_new, 0)),
      width = 2, pad = "0", side = "left"
    ),
    # Convert blue
    str_pad(
      string = as.hexmode(round(z_new, 0)),
      width = 2, pad = "0", side = "left"
    )
  )
  
  return(hex_value)
  
}

taxo_comm$color <- RGB(taxo_comm$Axis1, 
                       taxo_comm$Axis2, 
                       taxo_comm$Axis3, 
                       min_range = -0.50, 
                       max_range = 0.50)

phylo_comm$color <- RGB(phylo_comm$Axis1, 
                        phylo_comm$Axis2, 
                        phylo_comm$Axis3, 
                        min_range = -0.50, 
                        max_range = 0.50)

func_comm$color <- RGB(func_comm$Axis1, 
                       func_comm$Axis2, 
                       func_comm$Axis3, 
                       min_range = -0.50, 
                       max_range = 0.50)

#'------------------------------------------------------------------------------
#' @Plot-communities
#-------------------------------------------------------------------------------

# Plot collector
plot_list <- list()

# Viewing angles
theta <- c(seq(0, 180, length.out = 180), seq(-180, 0, length.out = 180))  # 36 positions in the full circle
phi <- c(seq(0, 0, length.out = length(theta)))
h_adjust_axis1 <- c(seq(1, 0, length.out = 180), seq(0, 1, length.out = 180))
h_adjust_axis2 <- c(seq(1, 1, length.out = 180), #0-180
                    seq(1, 0, length.out = 45), #180-225
                    seq(0, 0, length.out = 45), #225-270
                    seq(0, 0, length.out = 60),
                    seq(0, 1, length.out = 30)) #270-360

# Define plot features
th <- theme_void(base_size = 15) +
  theme(panel.background = element_rect(fill = "black",
                                        colour = "black",
                                        linewidth = 0.5, 
                                        linetype = "solid"),
        plot.background = element_rect(fill = "black"),
        legend.title = element_text(colour="white", size = 12),
        legend.text = element_text(color = "white", size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(5,5,5,5), 'pt'))


# Function to apply in parallel
application <- function(X,
                        taxo_comm, 
                        phylo_comm, 
                        func_comm, 
                        theta, 
                        phi, 
                        h_adjust_axis1,
                        h_adjust_axis2,
                        th,
                        export_folder) {
  
  i <- X
  
  taxo_p <- ggplot(taxo_comm, aes(x = Axis1, 
                                  y = Axis2, 
                                  z = Axis3,
                                  color = color)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            linewidth = 0.5,
            colour = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    scale_colour_identity() +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    th
    
  phylo_p <- ggplot(phylo_comm, aes(x = Axis1, 
                                    y = Axis2, 
                                    z = Axis3,
                                    colour = color)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            linewidth = 0.5,
            colour = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    scale_colour_identity() +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    th
  
  func_p <- ggplot(func_comm, aes(x = Axis1, 
                                  y = Axis2, 
                                  z = Axis3,
                                  color = color)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            linewidth = 0.5,
            colour = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    scale_colour_identity() +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    th
  
  name <- str_pad(i, 4, pad = "0")
  file_name <- paste0(export_folder, "/plot_communities_", name, ".png")
  #combined <- 
  png(file_name,
      width = 300, 
      height = 100, 
      units = "mm", 
      pointsize = 32,
      bg = NA,  
      res = 300)
  print(ggarrange(taxo_p, phylo_p, func_p, ncol = 3))
  dev.off()
  
  return(NULL)
  
}

# Export in parallel
complete <- pbmclapply(X = 300:length(theta), 
                       FUN = application, 
                       taxo_comm = taxo_comm, 
                       phylo_comm = phylo_comm, 
                       func_comm = func_comm, 
                       theta = theta, 
                       phi = phi, 
                       h_adjust_axis1 = h_adjust_axis1,
                       h_adjust_axis2 = h_adjust_axis2,
                       th = th,
                       export_folder = "/media/antonio/Work/MMDB/figures/video/LMDS",
                       mc.preschedule = TRUE, 
                       mc.set.seed = FALSE,
                       mc.silent = FALSE, 
                       mc.cores = 30,
                       mc.cleanup = TRUE, 
                       mc.allow.recursive = TRUE, 
                       affinity.list = NULL)

#'------------------------------------------------------------------------------
#' @Plot-probabilities
#-------------------------------------------------------------------------------

theta <- c(seq(0, 180, length.out = 180), seq(-180, 0, length.out = 180))  # 36 positions in the full circle
phi <- c(seq(0, 0, length.out = length(theta)))
h_adjust_axis1 <- c(seq(1, 0, length.out = 180), seq(0, 1, length.out = 180))
h_adjust_axis2 <- c(seq(1, 1, length.out = 180), #0-180
                    seq(1, 0, length.out = 45), #180-225
                    seq(0, 0, length.out = 45), #225-270
                    seq(0, 0, length.out = 60),
                    seq(0, 1, length.out = 30)) #270-360

th <- theme_void(base_size = 15) +
  theme(panel.background = element_rect(fill = "black",
                                        colour = "black",
                                        linewidth = 0.5, 
                                        linetype = "solid"),
        plot.background = element_rect(fill = "black"),
        legend.title = element_text(colour="white", size = 20),
        legend.text = element_text(color = "white", size = 15),
        legend.position = "bottom",
        plot.margin = unit(c(0,0,0,0), 'pt'),
        legend.title.align=0.5)


application <- function(X,
                        taxo_comm, 
                        phylo_comm, 
                        func_comm, 
                        taxo_prob,
                        phylo_prob,
                        func_prob,
                        legend, 
                        pallete,
                        lim,
                        export, 
                        class,
                        theta, 
                        phi, 
                        h_adjust_axis1,
                        h_adjust_axis2,
                        th) {
  
  i <- X
  
  taxo_p <- ggplot(taxo_comm, aes(x = Axis1, 
                                  y = Axis2, 
                                  z = Axis3,
                                  color = taxo_prob)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            linewidth = 0.5,
            colour = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    scale_colour_gradientn(legend, colours = pallete, limits = lim, oob = scales::squish) +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    guides(colour = guide_colourbar(barwidth = 20, barheight = 0.7, title.position = "top")) +
    th 
  
  phylo_p <- ggplot(phylo_comm, aes(x = Axis1, 
                                    y = Axis2, 
                                    z = Axis3,
                                    color = phylo_prob)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            linewidth = 0.5,
            colour = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    scale_colour_gradientn(legend, colours = pallete, limits = lim, oob = scales::squish) +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    guides(colour = guide_colourbar(barwidth = 20, barheight = 0.7, title.position = "bottom")) +
    th 
  
  func_p <- ggplot(func_comm, aes(x = Axis1, 
                                  y = Axis2, 
                                  z = Axis3,
                                  color = func_prob)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            linewidth = 0.5,
            colour = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    scale_colour_gradientn(legend, colours = pallete, limits = lim, oob = scales::squish) +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    guides(colour = guide_colourbar(barwidth = 20, barheight = 0.7, title.position = "bottom")) +
    th 
  
  name <- str_pad(i, 4, pad = "0")
  file_name = paste0(export, "/", class, "_", name, ".png", sep="")
  png(file_name,
      width = 600, 
      height = 210, 
      units = "mm", 
      pointsize = 32,
      bg = NA,  
      res = 300)
  print(ggarrange(taxo_p, phylo_p, func_p, ncol = 3, common.legend = TRUE,
                  legend = "bottom"))
  dev.off()
  
  return(NULL)
  
}

# Quercus
pbmclapply(X = 1:length(theta),  #length(theta)
           FUN = application, 
           taxo_comm, 
           phylo_comm, 
           func_comm, 
           taxo_prob = taxo_prob$Quercus,
           phylo_prob = phylo_prob$Quercus,
           func_prob = func_prob$Quercus,
           legend = "Occurrence probability of Quercus", 
           pallete = viridis(1000),
           lim = c(0, 1),
           export = "/home/antonio/Documents/3d/Quercus", 
           class = "Quercus",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 30,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

# Pinus
pbmclapply(X = 1:length(theta),  #
           FUN = application, 
           taxo_comm, 
           phylo_comm, 
           func_comm, 
           taxo_prob = taxo_prob$Pinus,
           phylo_prob = phylo_prob$Pinus,
           func_prob = func_prob$Pinus,
           legend = "Occurrence probability of Pinus", 
           pallete = viridis(1000),
           lim = c(0, 1),
           export = "/home/antonio/Documents/3d/Pinus", 
           class = "Pinus",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 30,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

# Acer
pbmclapply(X = 1:length(theta),  #
           FUN = application, 
           taxo_comm, 
           phylo_comm, 
           func_comm, 
           taxo_prob = taxo_prob$Acer,
           phylo_prob = phylo_prob$Acer,
           func_prob = func_prob$Acer,
           legend = "Occurrence probability of Acer", 
           pallete = viridis(1000),
           lim = c(0, 1),
           export = "/home/antonio/Documents/3d/Acer", 
           class = "Acer",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 30,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

#LMA
pbmclapply(X = 1:length(theta),  #
           FUN = application, 
           taxo_comm, 
           phylo_comm, 
           func_comm, 
           taxo_prob = taxo_cwm$LMA,
           phylo_prob = phylo_cwm$LMA,
           func_prob = func_cwm$LMA,
           legend = "CWM of Leaf Mass Area", 
           pallete = viridis::inferno(1000),
           lim = c(25, 175),
           export = "/home/antonio/Documents/3d/LMA", 
           class = "LMA",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 30,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

#WD
pbmclapply(X = 1:length(theta),  #
           FUN = application, 
           taxo_comm, 
           phylo_comm, 
           func_comm, 
           taxo_prob = taxo_cwm$WOOD_SPGR_GREENVOL_DRYWT,
           phylo_prob = phylo_cwm$WOOD_SPGR_GREENVOL_DRYWT,
           func_prob = func_cwm$WOOD_SPGR_GREENVOL_DRYWT,
           legend = "CWM of Wood Density", 
           pallete = viridis::inferno(1000),
           lim = c(0.3, 0.6),
           export = "/home/antonio/Documents/3d/WD", 
           class = "WD",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 30,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)
#N
pbmclapply(X = 1:length(theta),  #
           FUN = application, 
           taxo_comm, 
           phylo_comm, 
           func_comm, 
           taxo_prob = taxo_cwm$N,
           phylo_prob = phylo_cwm$N,
           func_prob = func_cwm$N,
           legend = "CWM of Nitrogen", 
           pallete = viridis::inferno(1000),
           lim = c(10, 20),
           export = "/home/antonio/Documents/3d/nitrogen", 
           class = "Nitrogen",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 30,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)
