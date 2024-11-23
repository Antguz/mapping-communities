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

#'------------------------------------------------------------------------------
#' Read information

# Folder
root_path <- "/media/antonio/Work/MMDB"


taxo_comm <- fread(paste0(root_path, "/04_NMDS/NMDS_scores_comm_taxo_total.csv"))
phylo_comm <- fread(paste0(root_path, "/04_NMDS/NMDS_scores_comm_phylo_total.csv"))
func_comm <- fread(paste0(root_path, "/04_NMDS/NMDS_scores_comm_func.csv"))

#'------------------------------------------------------------------------------
#' @Plot-functions
#-------------------------------------------------------------------------------

# RGB color
RGB <- function(X, Y, Z, min_range = -0.55, max_range = 0.55) {
  
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

# Compute RGB colors

taxo_comm$color <- RGB(taxo_comm$Axis1, 
                       taxo_comm$Axis2, 
                       taxo_comm$Axis3, 
                       min_range = -0.55, 
                       max_range = 0.55)

phylo_comm$color <- RGB(phylo_comm$Axis1, 
                        phylo_comm$Axis2, 
                        phylo_comm$Axis3, 
                        min_range = -0.55, 
                        max_range = 0.55)

func_comm$color <- RGB(func_comm$Axis1, 
                       func_comm$Axis2, 
                       func_comm$Axis3, 
                       min_range = -0.55, 
                       max_range = 0.55)

#'------------------------------------------------------------------------------
#' @Plot-communities
#-------------------------------------------------------------------------------

th <- theme_void(base_size = 15) +
  theme(panel.background = element_rect(fill = "black",
                                        colour = "black",
                                        linewidth = 0.5, 
                                        linetype = "solid"),
        plot.background = element_rect(fill = "black"),
        legend.title = element_text(colour="white", size = 12),
        legend.text = element_text(color = "white", size = 12),
        legend.position = "bottom",
        plot.margin = unit(c(0,0,0,0), 'pt'))


# Collector of plots to export
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

# Apply loop on viewing angles
for (i in 1:length(theta)) {
  
  taxo_p <- ggplot(taxo_comm, aes(x = Axis1, 
                         y = Axis2, 
                         z = -Axis3,
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
                                  z = -Axis3,
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
  
  func_p <- ggplot(func_comm, aes(x = Axis1, 
                                  y = Axis2, 
                                  z = -Axis3,
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
  
  p <- ggarrange(taxo_p, phylo_p, func_p, ncol = 3)# +
  
  plot_list[[i]] = p
  
}

# Define folder to export plots
export <- "/home/antonio/Documents/3d/communities"

# Export plots in a loop
for (i in 1:length(theta)) { #length(theta)
  
  name <- str_pad(i, 4, pad = "0")
  file_name = paste0(export, "/plot_communities_", name, ".png", sep="")
  png(file_name,
      width = 600, 
      height = 200, 
      units = "mm", 
      pointsize = 32,
      bg = NA,  
      res = 300)
  p
  dev.off()
  
}

################################################################################
## Just communities and species
################################################################################

plot_list <- list()

theta <- c(seq(0, 180, length.out = 180), seq(-180, 0, length.out = 180))  # 36 positions in the full circle
phi <- c(seq(0, 0, length.out = length(theta)))
h_adjust_axis1 <- c(seq(1, 0, length.out = 180), seq(0, 1, length.out = 180))

for (i in 1:length(theta)) {
  
  p <- ggplot(frame, aes(x = MDS1, 
                         y = MDS2, 
                         z = MDS3, 
                         colour = value, 
                         size = value)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            size = 0.5,
            color = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            alpha = 0.6) +
    scale_colour_manual(" ", values = c("grey", "darkred")) +
    scale_size_manual(" ", values = c(0.6, 2)) +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    theme_void() +
    theme(panel.background = element_rect(fill = "black",
                                          colour = "black",
                                          linewidth = 0.5, 
                                          linetype = "solid"),
          plot.background = element_rect(fill = "black"),
          legend.title = element_text(colour="white", size = 12),
          legend.text = element_text(color = "white", size = 12),
          legend.position = "bottom") 

  plot_list[[i]] = p
  
}

export <- "/home/antonio/Documents/3d/communities_species"

for (i in 1:length(theta)) { #length(theta)
  
  name <- str_pad(i, 4, pad = "0")
  file_name = paste0(export, "/plot_communities_species_", name, ".png", sep="")
  png(file_name,
      width = 210, 
      height = 210, 
      units = "mm", 
      pointsize = 12,
      bg = NA,  
      res = 600)
  print(plot_list[[i]])
  dev.off()
  
}



################################################################################
## Just angiosperms
################################################################################

plot_list <- list()

theta <- c(seq(0, 180, length.out = 180), seq(-180, 0, length.out = 180))  # 36 positions in the full circle
phi <- c(seq(0, 0, length.out = length(theta)))
h_adjust_axis1 <- c(seq(1, 0, length.out = 180), seq(0, 1, length.out = 180))

for (i in 1:length(theta)) {
  
  p <- ggplot(plots_fia, aes(x = MDS1, 
                             y = MDS2, 
                             z = MDS3, 
                             colour = A_to_T)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            size = 0.5,
            color = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    theme_void() +
    scale_colour_carto_c("Proportion of Angiosperms", 
                         type = "diverging",
                         palette = "Fall") +
    theme(panel.background = element_rect(fill = "black",
                                          colour = "black",
                                          size = 0.5, 
                                          linetype = "solid"),
          plot.background = element_rect(fill = "black"),
          legend.title = element_text(colour="white", size = 12),
          legend.text = element_text(color = "white", size = 12),
          legend.position = "bottom") +
    guides(color = guide_colorbar(title.position = "top",
                                  title.hjust = 0.5,
                                  barwidth = 20,
                                  barheight = 1))
  
  plot_list[[i]] = p
  
}

export <- "/home/antonio/Documents/3d/angiosperms"

for (i in 1:length(theta)) { #length(theta)
  
  name <- str_pad(i, 4, pad = "0")
  file_name = paste0(export, "/plot_angiosperms_", name, ".png", sep="")
  png(file_name,
      width = 210, 
      height = 210, 
      units = "mm", 
      pointsize = 12,
      bg = NA,  
      res = 600)
  print(plot_list[[i]])
  dev.off()
  
}



################################################################################
## Just mycorrhizal association
################################################################################

theta <- c(seq(0, 180, length.out = 180), seq(-180, 0, length.out = 180))  # 36 positions in the full circle
phi <- c(seq(0, 0, length.out = length(theta)))
h_adjust_axis1 <- c(seq(1, 0, length.out = 180), seq(0, 1, length.out = 180))

plot_list <- list()

for (i in 1:length(theta)) {
  
  p <- ggplot(plots_fia, aes(x = MDS1, 
                             y = MDS2, 
                             z = MDS3, 
                             colour = AM_to_T)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            size = 0.5,
            color = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    theme_void() +
    scale_colour_carto_c("Proportion of species with arbuscular mycorrhizal", 
                         type = "diverging",
                         palette = "TealRose") +
    theme(panel.background = element_rect(fill = "black",
                                          colour = "black",
                                          size = 0.5, 
                                          linetype = "solid"),
          plot.background = element_rect(fill = "black"),
          legend.title = element_text(colour="white", size = 12),
          legend.text = element_text(color = "white", size = 12),
          legend.position = "bottom") +
    guides(color = guide_colorbar(title.position = "top",
                                  title.hjust = 0.5,
                                  barwidth = 20,
                                  barheight = 1))
  
  plot_list[[i]] = p
  
}

export <- "/home/antonio/Documents/3d/mycorrizal"

for (i in 1:length(theta)) { #length(theta)
  
  name <- str_pad(i, 4, pad = "0")
  file_name = paste0(export, "/plot_mycorrizal_", name, ".png", sep="")
  png(file_name,
      width = 210, 
      height = 210, 
      units = "mm", 
      pointsize = 12,
      bg = NA,  
      res = 600)
  print(plot_list[[i]])
  dev.off()
  
}

#4000
#3240



root_path <- "/media/antonio/Work/Projects/DESIS_mapping"

taxo_comm <- fread(paste0(root_path, "/05_MDS/PCoA_scores_comm_taxo.csv"))
phylo_comm <- fread(paste0(root_path, "/05_MDS/PCoA_scores_comm_phylo-log.csv"))
func_comm <- fread(paste0(root_path, "/05_MDS/PCoA_scores_comm_func.csv"))

species_scores <- readRDS(paste0(root_path, "/05_ordination/species_scores_BA.rds"))
dimensions_plots <- readRDS(paste0(root_path, "/05_ordination/communities_scores_BA.rds"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"))
proportions <- fread(paste0(root_path, "/00_FIA/proportions_BA.csv"))


fig <- plot_ly(taxo_comm, x = ~PCoA1, y = ~PCoA2, z = ~PCoA3)

scatter3D(x = taxo_comm$PCoA1, y = taxo_comm$PCoA2, z = taxo_comm$PCoA3, 
          data = taxo_comm, 
          colvar = taxo_comm$PCoA1,
          phi = 20, bty = "g", 
          pch = 20, cex = 0.5, ticktype = "detailed",
          xlab = "PCoA1",
          ylab = "PCoA2",
          zlab = "PCoA3",
          xlim = c(-0.60, 0.60),
          ylim = c(-0.60, 0.60),
          zlim = c(-0.60, 0.60),
          nticks = 3)
