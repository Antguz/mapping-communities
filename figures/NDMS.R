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

#'------------------------------------------------------------------------------
#' Read information

root_path <- "E:/Projects/DESIS_mapping"
species_scores <- readRDS(paste0(root_path, "/05_ordination/species_scores_taxonomic_DRYBIO.rds"))
dimensions_plots <- readRDS(paste0(root_path, "/05_ordination/dimensions_plots_taxonomic_DRYBIO.rds"))
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


### Create plots ---------------------------------------------------------------

### Angiosperms

NMDS1 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS2, fill = A_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS2), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  geom_text_repel(data = species, aes(x = MDS1, y = MDS2),
                  label= species$`PLANTS-Code`, 
                  direction = "both", 
                  colour = col,
                  box.padding = 0.5,
                  max.overlaps = 20,
                  size = 3) + 
  xlab("NMDS1") + ylab("NMDS2") +
  scale_fill_distiller("Proportion of Angiosperms", palette = "BrBG")

NMDS2 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS3, fill = A_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS3), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  geom_text_repel(data = species, aes(x = MDS1, y = MDS3),
                  label= species$`PLANTS-Code`, 
                  direction = "both", 
                  colour = col,
                  box.padding = 0.5,
                  max.overlaps = 20,
                  size = 3) + 
  xlab("NMDS1") + ylab("NMDS3") +
  scale_fill_distiller("Proportion of Angiosperms", palette = "BrBG")

NMDS3 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS2, y = MDS3, fill = A_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS2, MDS3), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  geom_text_repel(data = species, aes(x = MDS2, y = MDS3),
                  label= species$`PLANTS-Code`, 
                  direction = "both", 
                  colour = col,
                  box.padding = 0.5,
                  max.overlaps = 20,
                  size = 3) + 
  xlab("NMDS2") + ylab("NMDS3") +
  scale_fill_distiller("Proportion of Angiosperms", palette = "BrBG")


jpeg("figures/NMDS_angiosperm_taxonomic.jpg", width = 280, height = 210, units = "mm", res = 300)
ggarrange(NMDS1, NMDS2, NMDS3, 
          labels = c("a", "b", "c"),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "right")
dev.off()

### Micorrizal

NMDS1 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS2, fill = AM_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS2), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  geom_text_repel(data = species, aes(x = MDS1, y = MDS2),
                  label= species$`PLANTS-Code`, 
                  direction = "both", 
                  colour = col,
                  box.padding = 0.5,
                  max.overlaps = 20,
                  size = 3) + 
  xlab("NMDS1") + ylab("NMDS2") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn")

NMDS2 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS3, fill = AM_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS3), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  geom_text_repel(data = species, aes(x = MDS1, y = MDS3),
                  label= species$`PLANTS-Code`, 
                  direction = "both", 
                  colour = col,
                  box.padding = 0.5,
                  max.overlaps = 20,
                  size = 3) + 
  xlab("NMDS1") + ylab("NMDS3") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn")

NMDS3 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS2, y = MDS3, fill = AM_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS2, MDS3), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  geom_text_repel(data = species, aes(x = MDS2, y = MDS3),
                  label= species$`PLANTS-Code`, 
                  direction = "both", 
                  colour = col,
                  box.padding = 0.5,
                  max.overlaps = 20,
                  size = 3) + 
  xlab("NMDS2") + ylab("NMDS3") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn")

jpeg("figures/NMDS_mycorrhizal.jpg", width = 210, height = 210, units = "mm", res = 300)
ggarrange(NMDS1, NMDS2, NMDS3, 
          labels = c("a", "b", "c"),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "bottom")
dev.off()


scatter3D(plots_fia$MDS1, plots_fia$MDS2, plots_fia$MDS3, clab = c("Proportion of", "Angiosperms"),
          col.var = plots_fia$A_to_T)

with(plots_fia, scatter3D(MDS1, MDS2, MDS3, 
                          labels = c("NMDS1", "NMDS2", "NMDS3"),
                          colvar = A_to_T,
                          col = brewer.pal(10, "BrBG"),
                          clab = c("Proportion of", "Angiosperms"),
                          pch = 20,
                          cex = 1,
                          xlab = "NMDS1", ylab ="NMDS2", zlab = "NMDS3"),
     points3D(species$MDS1, species$MDS2, species$MDS3, add = TRUE, col = "red"))

with(plots_fia, scatter3D(MDS1, MDS2, MDS3, 
                          labels = c("NMDS1", "NMDS2", "NMDS3"),
                          colvar = AM_to_T,
                          col = brewer.pal(10, "Spectral") ,
                          clab = c("Proportion of", "Angiosperms"),
                          pch = 20,
                          cex = 1,
                          xlab = "NMDS1", ylab ="NMDS2", zlab = "NMDS3"),
     points3D(species$MDS1, species$MDS2, species$MDS3, add = TRUE, col = "red"))

plotrgl()
rgl.snapshot(filename="figures/3dplot.png")

with(plots_fia, scatter3D(MDS1, MDS2, MDS3, 
                          labels = c("NMDS1", "NMDS2", "NMDS3"),
                          colvar = AM_to_T,
                          col = viridis(100, option = "magma", direction = 1),
                          clab = c("Proportion of", "trees with AM"),
                          pch = 20))




scatter3D(plots_fia$MDS1, plots_fia$MDS2, plots_fia$MDS3, 
          xlab = "NMDS1", ylab ="NMDS2", zlab = "NMDS3",
          colvar = plots_fia$A_to_T,
          col = viridis(100, option = "mako", direction = -1),
          clab = c("Proportion of", "Angiosperms"),
          pch = 20,
          cex = 0.75)

scatter3D(plots_fia$MDS1, plots_fia$MDS2, plots_fia$MDS3, 
          xlab = "NMDS1", ylab ="NMDS2", zlab = "NMDS3",
          colvar = plots_fia$AM_to_T,
          col = viridis(100, option = "cividis", direction = 1),
          clab = c("Proportion of", "trees with AM"),
          pch = 20,
          cex = 0.75)

text3D(species$MDS1, species$MDS2, species$MDS3, labels = species$`PLANTS-Code`, add = TRUE, col = "#990000")
points3D(species$MDS1, species$MDS2, species$MDS3, add = TRUE, col = "red")
plotrgl()






colors <- function(x, 
                   gradient = brewer.pal(11, "BrBG"), 
                   colsteps= 100) {
  
  col <- colorRampPalette(gradient)(colsteps)[findInterval(x, seq(min(x), max(x), length.out = colsteps))]
  
  return(col)
  
}

colors(plots_fia$A_to_T)





color <- data.frame( COLOR_VALUE = seq(-1000,1000,1), color.name=my.colors(2001))
reg1.with.color<-merge(reg1, color.df, by="COLOR_VALUE")


# Static chart
plot3d(plots_fia$MDS1, plots_fia$MDS2, plots_fia$MDS3, 
       xlab = "NMDS1", ylab ="NMDS2", zlab = "NMDS3",
       col = colors(plots_fia$A_to_T),
       type = "s", 
       radius = plots_fia$goodness*2000)

# We can indicate the axis and the rotation velocity
play3d( spin3d( axis = c(0, 0, 1), rpm = 20), duration = 10 )

# Save like gif
movie3d(
  movie = "NDMS", 
  spin3d( axis = c(0, 0, 1), rpm = 7),
  duration = 10, 
  dir = root_path,
  type = "gif", 
  webshot = TRUE,
)
