library(gg3D)
library(stringr)
library(bit64)
library(plotly)
library(plot3D)


root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme Pro/Projects/DESIS_mapping"
#root_path <- "F:/Projects/DESIS_mapping"
species_scores <- readRDS(paste0(root_path, "/05_ordination/species_scores_BA.rds"))
dimensions_plots <- readRDS(paste0(root_path, "/05_ordination/communities_scores_BA.rds"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"))
proportions <- fread(paste0(root_path, "/00_FIA/proportions_BA.csv"))

#'------------------------------------------------------------------------------
#' Arrange information

species <- merge(species_scores, SPCD, by = "Taxon2", all.x = TRUE, all.y = FALSE)
plots_fia <- merge(proportions, dimensions_plots, by = "PLT_CN", all.x = FALSE, all.y = TRUE)

normalization <- dimensions_plots[, c("MDS1", "MDS2", "MDS3")]
normalization$MDS1 <- (normalization$MDS1-min(normalization$MDS1))/(max(normalization$MDS1)-min(normalization$MDS1))
normalization$MDS2 <- (normalization$MDS2-min(normalization$MDS2))/(max(normalization$MDS2)-min(normalization$MDS2))
normalization$MDS3 <- (normalization$MDS3-min(normalization$MDS3))/(max(normalization$MDS3)-min(normalization$MDS3))

RGB <- rgb(normalization$MDS1, normalization$MDS2, normalization$MDS3)

################################################################################
## Just communities
################################################################################

species$value <- "Species"
plots_fia$value <- "Communities"
plots_fia$RGB <- RGB

frame <- rbind(species[, c("MDS1", "MDS2", "MDS3", "value")],
               plots_fia[, c("MDS1", "MDS2", "MDS3", "value")])

plot_list <- list()

theta <- c(seq(0, 180, length.out = 180), seq(-180, 0, length.out = 180))  # 36 positions in the full circle
phi <- c(seq(0, 0, length.out = length(theta)))
h_adjust_axis1 <- c(seq(1, 0, length.out = 180), seq(0, 1, length.out = 180))
h_adjust_axis2 <- c(seq(1, 1, length.out = 180), #0-180
                    seq(1, 0, length.out = 45), #180-225
                    seq(0, 0, length.out = 45), #225-270
                    seq(0, 0, length.out = 60),
                    seq(0, 1, length.out = 30)) #270-360

for (i in 1:length(theta)) {
  
  p <- ggplot(plots_fia, aes(x = MDS1, 
                         y = MDS2, 
                         z = MDS3, 
                         color = value)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            size = 0.5,
            colour = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    scale_color_manual(" ", values = "black") +
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

export <- "/home/antonio/Documents/3d/communities"

for (i in 1:length(theta)) { #length(theta)
  
  name <- str_pad(i, 4, pad = "0")
  file_name = paste0(export, "/plot_communities_", name, ".png", sep="")
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
