#'------------------------------------------------------------------------------
#' @title Codes to visualize the NDMS
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(bigmds)
library(ggplot2)
library(bit64)
library(ggrepel)
options(scipen = 100)
library(viridis)
library(ggpointdensity)
library(ggpubr)
library(plot3D)
library(plot3Drgl)
library(rgl)
library(magick)
library(RColorBrewer)
library(scales)
library(ggpubr)

#'------------------------------------------------------------------------------
#' Read information

root_path <- "/media/antonio/Extreme Pro/Projects/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"
species_scores <- readRDS(paste0(root_path, "/05_ordination/species_scores_phylogenetic_DRYBIO.rds"))
dimensions_plots <- readRDS(paste0(root_path, "/05_ordination/dimensions_plots_phylogenetic_DRYBIO.rds"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"))
proportions <- fread(paste0(root_path, "/00_FIA/proportions.csv"))

#'------------------------------------------------------------------------------
#' Arrange information

species <- merge(species_scores, SPCD, by = "Taxon2", all.x = TRUE, all.y = FALSE)
plots_fia <- merge(proportions, dimensions_plots, by = "PLT_CN", all.x = FALSE, all.y = TRUE)
plots_fia <- na.exclude(plots_fia)

th <- theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
tam <- 12
col <- "black"
sp_fill <- "darkred"
sp_color <- "white"
sp_size <- 2
plot_range <- c(0.5, 2)

### Create plots ---------------------------------------------------------------
### Grey

rNMDS1 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS2, size = goodness), shape = 21, fill = "grey75", colour = "grey95", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS2), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  geom_text_repel(data = species, aes(x = MDS1, y = MDS2),
                  label= species$`PLANTS-Code`, 
                  direction = "both", 
                  colour = sp_fill,
                  box.padding = 0.5,
                  max.overlaps = 20,
                  size = sp_size) + 
  xlab("NMDS1") + ylab("NMDS2") +
  scale_fill_distiller("Proportion of Angiosperms", palette = "BrBG")

rNMDS2 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS3, size = goodness), shape = 21, fill = "grey75", colour = "grey95", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  geom_text_repel(data = species, aes(x = MDS1, y = MDS3),
                  label= species$`PLANTS-Code`, 
                  direction = "both", 
                  colour = sp_fill,
                  box.padding = 0.5,
                  max.overlaps = 20,
                  size = sp_size) + 
  xlab("NMDS1") + ylab("NMDS3") +
  scale_fill_distiller("Proportion of Angiosperms", palette = "BrBG")

rNMDS3 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS2, y = MDS3, size = goodness), shape = 21, fill = "grey75", colour = "grey95", alpha = 0.5) +
  geom_point(data = species, aes(MDS2, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  geom_text_repel(data = species, aes(x = MDS2, y = MDS3),
                  label= species$`PLANTS-Code`, 
                  direction = "both", 
                  colour = sp_fill,
                  box.padding = 0.5,
                  max.overlaps = 20,
                  size = sp_size) + 
  xlab("NMDS2") + ylab("NMDS3") +
  scale_fill_distiller("Proportion of Angiosperms", palette = "BrBG")

### Angiosperms ----------------------------------------------------------------

aNMDS1 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS2, fill = A_to_T, size = goodness), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS2), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("NMDS1") + ylab("NMDS2") +
  scale_fill_distiller("Proportion of Angiosperms", palette = "BrBG")

aNMDS2 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS3, fill = A_to_T, size = goodness), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("NMDS1") + ylab("NMDS3") +
  scale_fill_distiller("Proportion of Angiosperms", palette = "BrBG")

aNMDS3 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS2, y = MDS3, fill = A_to_T, size = goodness), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS2, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("NMDS2") + ylab("NMDS3") +
  scale_fill_distiller("Proportion of Angiosperms", palette = "BrBG")

### Micorrizal -----------------------------------------------------------------

mNMDS1 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS2, fill = AM_to_T, size = goodness), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS2), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("NMDS1") + ylab("NMDS2") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn")

mNMDS2 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS3, fill = AM_to_T, size = goodness), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("NMDS1") + ylab("NMDS3") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn")

mNMDS3 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS2, y = MDS3, fill = AM_to_T, size = goodness), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS2, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("NMDS2") + ylab("NMDS3") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn")

#Multiple panel -----------------------------------------------------------------------

p <- ggarrange(rNMDS1, rNMDS2, rNMDS3, 
          aNMDS1, aNMDS2, aNMDS3,
          mNMDS1, mNMDS2, mNMDS3, 
          ncol=3, nrow=3, common.legend = TRUE, legend="bottom")


jpeg("figures/NMDS_taxonomic.jpg", width = 210, height = 210, units = "mm", res = 300)
p
dev.off()



#3D plot -----------------------------------------------------------------------
colors <- function(x, 
                   gradient = rev(brewer.pal(11, "BrBG")), 
                   colsteps= 100) {
  
  col <- colorRampPalette(gradient)(colsteps)[findInterval(x, seq(min(x), max(x), length.out = colsteps))]
  return(col)
}


# Static chart
plot3d(species$MDS1, species$MDS2, species$MDS3, 
       xlab = "NMDS1", ylab ="NMDS2", zlab = "NMDS3",
       col = rep("darkred", nrow(species)),
       type = "s", 
       aspect = TRUE,
       decorate = FALSE,
       radius = 0.05,
       fov = 1)
text3d(species$MDS1, species$MDS2, (species$MDS3+0.10), text = species$`PLANTS-Code`)

plot3d(plots_fia$MDS1, plots_fia$MDS2, plots_fia$MDS3, 
       xlab = "NMDS1", ylab ="NMDS2", zlab = "NMDS3",
       col = rep("grey75", nrow(species)),
       type = "s", 
       aspect = TRUE,
       decorate = FALSE,
       add = TRUE,
       radius = plots_fia$goodness*3000)


plot3d(plots_fia$MDS1, plots_fia$MDS2, plots_fia$MDS3, 
       xlab = "NMDS1", ylab ="NMDS2", zlab = "NMDS3",
       col = colors(plots_fia$A_to_T, gradient = rev(brewer.pal(11, "BrBG"))),
       type = "s", 
       aspect = TRUE,
       decorate = FALSE,
       add = TRUE,
       radius = plots_fia$goodness*3000)

plot3d(plots_fia$MDS1, plots_fia$MDS2, plots_fia$MDS3, 
       xlab = "NMDS1", ylab ="NMDS2", zlab = "NMDS3",
       col = colors(plots_fia$AM_to_T, gradient = rev(brewer.pal(11, "PRGn"))),
       type = "s", 
       aspect = TRUE,
       decorate = FALSE,
       add = TRUE,
       radius = plots_fia$goodness*3000)

# We can indicate the axis and the rotation velocity
Sys.sleep(5)
play3d( spin3d( axis = c(0, 0, 1), rpm = 2), duration = 30)

# Save like gif
movie3d(
  movie = "NDMS", 
  spin3d( axis = c(0, 0, 1), rpm = 2),
  duration = 30, 
  dir = root_path,
  type = "gif", 
  webshot = TRUE,
)
