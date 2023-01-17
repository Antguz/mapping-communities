#'------------------------------------------------------------------------------
#' @title Non-metric Multidimensional Scaling
#'------------------------------------------------------------------------------

#' @description It creates a new ordination of species abundance of plots
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
library(doSNOW)
library(foreach)
library(ecodist)
library(bigmds)
library(ape)
library(ggplot2)
library(bit64)
library(ggrepel)
#library(BiodiversityR)
options(scipen = 100)
library(viridis)
library(ggpointdensity)
library(ggpubr)

#'------------------------------------------------------------------------------
#' Arguments
#' @param stand
#' @param out_path

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
ic_MDP <- readRDS(paste0(root_path, "/03_diversity/MDP_distance.rds"))
stand <- readRDS(paste0(root_path, "/03_diversity/master_community.rds"))
out_path <- paste0(root_path, "/05_ordination")

#'------------------------------------------------------------------------------
#' @analysis  

#-------------------------------------------------------------------------------
# Run the NMDS on taxonomic
dist <- vegdist(stand, method="bray")
NMDS_taxonomic <- metaMDS(dist, model = "linear", autotransform = FALSE, k = 3, try = 10, 
                          trymax = 30, trace = 1, plot = FALSE, engine = "monoMDS", 
                          maxit = 99, pc = TRUE, parallel = 10)

saveRDS(NMDS_taxonomic, paste0(out_path, "/NMDS_taxonomy.rds"))

# Run the NMDS on phylogenetic
NMDS_phylogenetic <- metaMDS(ic_MDP, model = "linear", autotransform = FALSE, k = 3, try = 10, 
                             trymax = 30, trace = 1, plot = FALSE, engine = "monoMDS", 
                             maxit = 50, pc = TRUE, parallel = 5)

saveRDS(NMDS_phylogenetic, paste0(out_path, "/NMDS_phylogenetic.rds"))

#-------------------------------------------------------------------------------
# Read and inspect
NMDS_taxonomic <- readRDS(paste0(out_path, "/NMDS_taxonomy.rds"))
NMDS_phylogenetic <- readRDS(paste0(out_path, "/NMDS_phylogenetic.rds"))

# Plot for inspection
plot(NMDS_taxonomic, choices = c(1, 2))
plot(NMDS_taxonomic, choices = c(1, 3))
plot(NMDS_taxonomic, choices = c(2, 3))

# Plot for inspection
plot(NMDS_phylogenetic, choices = c(1, 2))
plot(NMDS_phylogenetic, choices = c(1, 3))
plot(NMDS_phylogenetic, choices = c(2, 3))

#-------------------------------------------------------------------------------
# Get goodness
god_phylogenetic <- goodness(NMDS_phylogenetic)
saveRDS(god_phylogenetic, paste0(out_path, "/goodness_phylogenetic.rds"))

god_taxonomic <- goodness(NMDS_taxonomic)
saveRDS(god_taxonomic, paste0(out_path, "/goodness_taxonomic.rds"))

#-------------------------------------------------------------------------------
# Get stress plot
stressplot(NMDS_phylogenetic)
stressplot(NMDS_taxonomic)

#-------------------------------------------------------------------------------
# Get weighted averages scores of species
species_scores_taxonomic <- wascores(NMDS_taxonomic$points, stand, expand = FALSE)
species_names_taxonomic <- rownames(species_scores_taxonomic)
species_scores_taxonomic <- cbind(Taxon2 = species_names_taxonomic, as.data.table(species_scores_taxonomic))
saveRDS(species_scores_taxonomic, paste0(out_path, "/species_scores_taxonomic.rds"))

species_scores_phylogenetic <- wascores(NMDS_phylogenetic$points, stand, expand = FALSE)
species_names_phylogenetic <- rownames(species_scores_phylogenetic)
species_scores_phylogenetic <- cbind(Taxon2 = species_names_phylogenetic, as.data.table(species_scores_phylogenetic))
saveRDS(species_scores_phylogenetic, paste0(out_path, "/species_scores_phylogenetic.rds"))

#-------------------------------------------------------------------------------
# Get dimensions
plot_names_taxonomic <- as.integer64(rownames(NMDS_taxonomic$points))
dimensions_plots_taxonomic <- as.data.table(NMDS_taxonomic$points)
dimensions_plots_taxonomic$PLT_CN = plot_names_taxonomic
saveRDS(dimensions_plots_taxonomic, paste0(out_path, "/dimensions_plots_taxonomic.rds"))

plot_names_phylogenetic <- as.integer64(rownames(NMDS_phylogenetic$points))
dimensions_plots_phylogenetic <- as.data.table(NMDS_phylogenetic$points)
dimensions_plots_phylogenetic$PLT_CN = plot_names_phylogenetic
saveRDS(dimensions_plots_phylogenetic, paste0(out_path, "/dimensions_plots_phylogenetic.rds"))

### Create plots ---------------------------------------------------------------

species_scores <- species_scores_taxonomic
dimensions_plots <- dimensions_plots_taxonomic

SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"))
proportions <- fread(paste0(root_path, "/00_FIA/proportions.csv"))

species <- merge(species_scores, SPCD, by = "Taxon2", all.x = TRUE, all.y = FALSE)
plots_fia <- merge(proportions, dimensions_plots, by = "PLT_CN", all.x = FALSE, all.y = TRUE)

th <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
tam <- 12
col <- "black"

### Angiosperms

NMDS1 <- ggplot() +
  geom_hline(yintercept= 0, linewidth = 0.1) + geom_vline(xintercept= 0, linewidth = 0.1) +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS2, fill = A_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS2), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
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
  geom_hline(yintercept= 0, linewidth = 0.1) + geom_vline(xintercept= 0, linewidth = 0.1) +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS3, fill = A_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS3), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
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
  geom_hline(yintercept= 0, linewidth = 0.1) + geom_vline(xintercept= 0, linewidth = 0.1) +
  geom_point(data = plots_fia, aes(x = MDS2, y = MDS3, fill = A_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS2, MDS3), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
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
  geom_hline(yintercept= 0, linewidth = 0.1) + geom_vline(xintercept= 0, linewidth = 0.1) +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS2, fill = AM_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS2), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
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
  geom_hline(yintercept= 0, linewidth = 0.1) + geom_vline(xintercept= 0, linewidth = 0.1) +
  geom_point(data = plots_fia, aes(x = MDS1, y = MDS3, fill = AM_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS1, MDS3), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
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
  geom_hline(yintercept= 0, linewidth = 0.1) + geom_vline(xintercept= 0, linewidth = 0.1) +
  geom_point(data = plots_fia, aes(x = MDS2, y = MDS3, fill = AM_to_T), shape = 21, colour = "grey75", size = 2, alpha = 0.5) +
  geom_point(data = species, aes(MDS2, MDS3), colour = "white", fill = col, shape = 21, size = 2) +
  theme_bw(base_size = tam) + th +
  scale_x_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
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





install.packages("plot3D")
library(plot3D)
library(plot3Drgl)

scatter3D(plots_fia$MDS1, plots_fia$MDS2, plots_fia$MDS3, clab = c("Proportion of", "Angiosperms"),
          col.var = plots_fia$A_to_T)

with(plots_fia, scatter3D(MDS1, MDS2, MDS3, 
                          labels = c("NMDS1", "NMDS2", "NMDS3"),
                          colvar = A_to_T,
                          col = brewer.pal(10, "BrBG") ,
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
          labels = c("NMDS1", "NMDS2", "NMDS3"),
          colvar = plots_fia$A_to_T,
          col = viridis(100, option = "mako", direction = 1),
          clab = c("Proportion of", "Angiosperms"),
          pch = 20)
text3D(species$MDS1, species$MDS2, species$MDS3, labels = species$`PLANTS-Code`, add = TRUE, col = "#990000")
points3D(species$MDS1, species$MDS2, species$MDS3, add = TRUE, col = "red")
plotrgl()
