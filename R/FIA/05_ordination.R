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
library(ggrepel)
#library(BiodiversityR)
options(scipen = 100)
library(viridis)
library(ggpointdensity)

#'------------------------------------------------------------------------------
#' Arguments
#' @param stand
#' @param out_path

root_path <- "F:/DESIS"
root_path <- "/media/antonio/Extreme SSD/DESIS"
stand <- fread(paste0(root_path, "/04_match_spectra-stand/species.csv"))
out_path <- paste0(root_path, "/05_NMDS")

#'------------------------------------------------------------------------------
#' @analysis  

#Metadata
meta <- stand[,1:5]
data <- stand[, .SD, .SDcols=c(6:ncol(stand))]

#Remove plots without species and less than 5 species
nsp <- apply(data, 1, function(x) sum(x > 0))
meta <- meta[nsp != 0]
data <- data[nsp != 0]

#Remove rare species (presence <10 in plots)
presence <- colSums( data != 0)
presence <- names(presence[presence < 5])
remove <- NA

for(i in 1:nrow(data)) {
  
  value <- sum(data[i, ..presence])
  
  if(value > 0) {
    remove[i] <- TRUE
  } else {
    remove[i] <- FALSE
  }
  
}

data <- data[isFALSE(remove) == FALSE]
meta <- meta[isFALSE(remove) == FALSE]
data <- data[, !..presence]

#Transform
data_transfrom <- sqrt(data)
data_dist <- vegdist(data_transfrom, "bray")
data_dist <- readRDS("MPS.rds")
labels <- as.matrix(data_dist)
PCOA <- pcoa(data_dist, rn = rownames(labels))
saveRDS(PCOA, "PCoA_PD.rds")
PCOA <- readRDS("PCoA_PD.rds")

biplot(PCOA, plot.axes = c(1,2))

frame <- as.data.frame(matched$comm)
frame$PLOT_CN <- rownames(frame)



plots_fia <- as.data.table(PCOA$vectors)
plots_fia$PLOT_CN <- rownames(PCOA$vectors)

species <- wascores(x = PCOA$vectors, w = data_transfrom)
SPCD <- rownames(species)
species <- as.data.table(species)

plots_fia <- cbind(meta, plots_fia)
species <- cbind(data.table(species = SPCD), species)

#Export
fwrite(species, paste0(out_path, "/PCoA_species.csv"))
fwrite(plots_fia, paste0(out_path, "/PCoA_plots.csv"))

species <- fread(paste0(out_path, "/PCoA_species_PD.csv"), header = TRUE)
SPCD <- fread("/media/antonio/Extreme SSD/DESIS/00_FIA/SPCD.csv")
species <- merge(species, SPCD, by.x = "species", by.y = "Taxon2", all.x = TRUE)

plots_fia <- fread(paste0(out_path, "/PCoA_plots_PD.csv"), header = TRUE)
plots <- fread("/media/antonio/Extreme SSD/DESIS/00_FIA/plots_WI_MN_20221102.csv")
plots_fia <- merge(plots_fia, plots, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
proportions <- fread("/media/antonio/antonio4/Github/mapping-species/proportions.csv")
plots_fia <- merge(plots_fia, proportions, by = "PLT_CN", all.x = TRUE, all.y = FALSE)

plot(plots_fia$Axis.1, plots_fia$Axis.2)
points(species$Axis.20, species$Axis.10, col = "red")

jpeg("MPS_PCoA_2-3.jpg", width = 210, height = 150, units = "mm", res = 300)
ggplot() +
  geom_point(data = plots_fia, aes(x = Axis.2, y = Axis.3, fill = AM_to_T), shape = 21, colour = "grey", size = 2) +
  #geom_pointdensity(data = plots_fia, aes(x = Axis.1, y = Axis.2), shape = 19, alpha = 0.9, adjust = 40) +
  #stat_density_2d(data = plots_fia, aes(x = Axis.1, y = Axis.2, fill = ..level..), alpha = 0.5, geom = "polygon", contour = TRUE) +
  #stat_density_2d(data = plots_fia, aes(x = Axis.1, y = Axis.2, fill = after_stat(level)), alpha = 0.2, geom = "polygon", contour = TRUE) +
  #scale_fill_continuous(low="lavenderblush", high="red") +
  #geom_jitter(alpha=0.5, size = 1.1) +
  geom_point(data = species, aes(Axis.2, Axis.3), colour = "white", fill = "black", shape = 21, size = 3) +
  theme_bw() +
  #coord_fixed(xlim = c(-420, 420), ylim = c(-150, 150)) +
  scale_x_continuous(limits = c(-420, 420)) +
  scale_y_continuous(limits = c(-150, 150)) +
  geom_text_repel(data = species, aes(x = Axis.2, y = Axis.3),
    label= species$`PLANTS-Code`, 
    direction = "both", 
    box.padding = 0.5,
    max.overlaps = 40) + xlab("PCoA2") + ylab("PCoA3") +
  #theme(legend.position="none") +
  #scale_color_viridis(option= "mako", direction = -1)
  scale_fill_distiller("Proportion of Angiosperms", palette = "BrBG")
  #scale_fill_viridis("Proportion of Angiosperms", option = "mako")
dev.off()


jpeg("MPS_PCoA_1-2_AM.jpg", width = 210, height = 150, units = "mm", res = 300)
ggplot() +
  geom_point(data = plots_fia, aes(x = Axis.1, y = Axis.2, fill = AM_to_T), shape = 21, colour = "grey", size = 2) +
  #geom_pointdensity(data = plots_fia, aes(x = Axis.1, y = Axis.2), shape = 19, alpha = 0.9, adjust = 40) +
  #stat_density_2d(data = plots_fia, aes(x = Axis.1, y = Axis.2, fill = ..level..), alpha = 0.5, geom = "polygon", contour = TRUE) +
  #stat_density_2d(data = plots_fia, aes(x = Axis.1, y = Axis.2, fill = after_stat(level)), alpha = 0.2, geom = "polygon", contour = TRUE) +
  #scale_fill_continuous(low="lavenderblush", high="red") +
  #geom_jitter(alpha=0.5, size = 1.1) +
  geom_point(data = species, aes(Axis.1, Axis.2), colour = "white", fill = "black", shape = 21, size = 3) +
  theme_bw() +
  #coord_fixed(xlim = c(-420, 420), ylim = c(-150, 150)) +
  scale_x_continuous(limits = c(-420, 420)) +
  scale_y_continuous(limits = c(-150, 150)) +
  geom_text_repel(data = species, aes(x = Axis.1, y = Axis.2),
                  label= species$`PLANTS-Code`, 
                  direction = "both", 
                  box.padding = 0.5,
                  max.overlaps = 40) + xlab("PCoA1") + ylab("PCoA2") +
  #theme(legend.position="none") +
  #scale_color_viridis(option= "mako", direction = -1)
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn")
#scale_fill_viridis("Proportion of Angiosperms", option = "mako")
dev.off()




