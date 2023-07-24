#'------------------------------------------------------------------------------
#' @title Codes to visualize the NDMS
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(ggplot2)
library(bit64)
options(scipen = 100)
library(ggrepel)
library(rcartocolor)


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

th <- theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill= "transparent"),
            panel.background = element_rect(fill = "transparent"),
            rect = element_rect(fill = "transparent"),
            legend.position="top")

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
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS2), shape = 21, fill = "grey75", colour = "grey95", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS2), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-0.55, 0.40), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-0.47, 0.40), expand = c(0, 0)) +
  scale_size(range = plot_range) +
  geom_text_repel(data = species, aes(x = MDS1, y = MDS2),
                  label= species$`PLANTS-Code-new`, 
                  direction = "both", 
                  colour = sp_fill,
                  box.padding = 0.2,
                  max.overlaps = 5,
                  size = sp_size) + 
  xlab("Axis 1") + ylab("Axis 2")

rNMDS2 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS3), shape = 21, fill = "grey75", colour = "grey95", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-0.55, 0.40), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-0.53, 0.35), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  geom_text_repel(data = species, aes(x = MDS1, y = MDS3),
                  label= species$`PLANTS-Code-new`, 
                  direction = "both", 
                  colour = sp_fill,
                  box.padding = 0.2,
                  max.overlaps = 5,
                  size = sp_size) + 
  xlab("Axis 1") + ylab("Axis 3")

rNMDS3 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS2, y = MDS3), shape = 21, fill = "grey75", colour = "grey95", alpha = 0.5) +
  geom_point(data = species, aes(MDS2, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-0.47, 0.40), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-0.53, 0.35), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  geom_text_repel(data = species, aes(x = MDS2, y = MDS3),
                  label= species$`PLANTS-Code-new`, 
                  direction = "both", 
                  colour = sp_fill,
                  box.padding = 0.2,
                  max.overlaps = 5,
                  size = sp_size) + 
  xlab("Axis 2") + ylab("Axis 3")

### Angiosperms ----------------------------------------------------------------

aNMDS1 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS2, fill = A_to_T), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS2), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-0.55, 0.40), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-0.47, 0.40), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("Axis 1") + ylab("Axis 2") +
  scale_fill_carto_c("Proportion of Angiosperms", 
                     type = "diverging", palette = "Fall", direction = 1) +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7))

aNMDS2 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS3, fill = A_to_T), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-0.55, 0.40), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-0.53, 0.35), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("Axis 1") + ylab("Axis 3") +
  scale_fill_carto_c("Proportion of Angiosperms", 
                     type = "diverging", palette = "Fall", direction = 1) +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7))

aNMDS3 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS2, y = MDS3, fill = A_to_T), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS2, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-0.47, 0.40), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-0.53, 0.35), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("Axis 2") + ylab("Axis 3") +
  scale_fill_carto_c("Proportion of Angiosperms", 
                     type = "diverging", palette = "Fall", direction = 1) +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7))

### Micorrizal -----------------------------------------------------------------

mNMDS1 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS2, fill = AM_to_T), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS2), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-0.65, 0.65), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-0.65, 0.65), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("Axis 1") + ylab("Axis 2") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn") +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7))

mNMDS2 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS3, fill = AM_to_T), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-0.65, 0.65), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-0.65, 0.65), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("Axis 1") + ylab("Axis 3") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn") +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7))

mNMDS3 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1, colour = "grey") + 
  geom_vline(xintercept= 0, linewidth = 0.1, colour = "grey") +
  geom_point(data = plots_fia, aes(x = MDS2, y = MDS3, fill = AM_to_T), shape = 21, colour = "grey75", alpha = 0.5) +
  geom_point(data = species, aes(MDS2, MDS3), colour = sp_color, fill = sp_fill, shape = 21, size = sp_size) +
  theme_bw(base_size = tam) + th +
  #scale_x_continuous(limits = c(-0.65, 0.65), expand = c(0, 0)) +
  #scale_y_continuous(limits = c(-0.65, 0.65), expand = c(0, 0)) +
  scale_size(range = plot_range) + 
  xlab("Axis 2") + ylab("Axis 3") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn") +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7))

#Multiple panel -----------------------------------------------------------------------

p <- ggarrange(rNMDS1, rNMDS2, rNMDS3, 
               aNMDS1, aNMDS2, aNMDS3,
               mNMDS1, mNMDS2, mNMDS3, 
               ncol=3, nrow=3, common.legend = TRUE, legend="bottom")

a <- ggarrange(rNMDS1, rNMDS2, rNMDS3, ncol=3, nrow=1, common.legend = TRUE)
b <- ggarrange(aNMDS1, aNMDS2, aNMDS3, ncol=3, nrow=1, common.legend = TRUE)
c <- ggarrange(mNMDS1, mNMDS2, mNMDS3, ncol=3, nrow=1, common.legend = TRUE)

p <- ggarrange(a,b,c, ncol=1, nrow=3)


jpeg("figures/NMDS_betaB.jpg", width = 210, height = 210, units = "mm", res = 300)
p
dev.off()

png("figures/NMDS_betaA.png", width = 210, height = 210, units = "mm", res = 300, bg = "transparent")
p
dev.off()

png("figures/NMDS_betaB.png", width = 210, height = 210, units = "mm", res = 300, bg = "transparent")
p
dev.off()

