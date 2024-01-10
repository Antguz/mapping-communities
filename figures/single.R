#'------------------------------------------------------------------------------
#' @title Figure to plot spectra of inventories with a species
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)
library(bit64)
library(ggdark)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
sensor <- "EMIT"

#'------------------------------------------------------------------------------
#' Read information

spectra_SG <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor, "_spectra_SG.csv"))
spectra_VN <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor, "_spectra_VN.csv"))
spectra_SW <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor, "_spectra_SW.csv"))
communities <- fread(paste0(root_path, "/03_diversity/", "proportions.csv"))

#'------------------------------------------------------------------------------
#' Data reshaping

# Select communities with a species
communities <- communities[SR == 1, ]

spectra_SG <- merge(communities, spectra_SG, by = "PLT_CN", all.x = FALSE, all.y = FALSE)
spectra_SW <- merge(communities, spectra_SW, by = "PLT_CN", all.x = FALSE, all.y = FALSE)

# Collector
results_SG <- data.table()
results_SW <- data.table()

# Unique species
unique_species <- unique(spectra_SG$species)

for(i in 1:length(unique_species)) {
  
  species_selection <- unique_species[i]
  species_SG <- subset(spectra_SG, species == species_selection)
  species_SW <- subset(spectra_SW, species == species_selection)
  
  if(sensor == "DESIS") {
    
    species_SG <- species_SG[, .SD, .SDcols = c(47:254)] #
    species_SW <- species_SW[, .SD, .SDcols = c(45:222)] #440 - 971
    
    species_SG <- colMeans(species_SG)
    species_SW <- colMeans(species_SW)
    
    sub_results_SG <- as.data.table(matrix(c(species_selection, species_SG), nrow = 1))
    colnames(sub_results_SG) <- c("species", names(species_SG))
    
    sub_results_SW <- as.data.table(matrix(c(species_selection, species_SW), nrow = 1))
    colnames(sub_results_SW) <- c("species", names(species_SW))
    
    results_SG <- rbind(results_SG, sub_results_SG)
    results_SW <- rbind(results_SW, sub_results_SW)
    
  } else if(sensor == "EMIT") {
    
    # Bands from 447.74 to 1265.17, from 1488.07 to 1710.97, from 2015.60 to 2424.25
    species_SG <- species_SG[, .SD, .SDcols = c(38:148,
                                                177:208,
                                                248:304)] #
    species_SW <- species_SW[, .SD, .SDcols = c(38:148,
                                                178:208,
                                                249:304)] #440 - 971
    
    species_SG <- colMeans(species_SG)
    species_SW <- colMeans(species_SW)
    
    sub_results_SG <- as.data.table(matrix(c(species_selection, species_SG), nrow = 1))
    colnames(sub_results_SG) <- c("species", names(species_SG))
    
    sub_results_SW <- as.data.table(matrix(c(species_selection, species_SW), nrow = 1))
    colnames(sub_results_SW) <- c("species", names(species_SW))
    
    results_SG <- rbind(results_SG, sub_results_SG)
    results_SW <- rbind(results_SW, sub_results_SW)
    
  }
}

fwrite(results_SG, paste0(root_path, "/03_reshape/", sensor, "_species_SG.csv"))
fwrite(results_SW, paste0(root_path, "/03_reshape/", sensor, "_species_SW.csv"))

# Melt for plot

SG_frame <- melt(results_SG, id.vars = c("species"), 
                 variable.name = "wavelength",
                 value.name = "value")
SG_frame$wavelength <- as.numeric(as.character(SG_frame$wavelength))
SG_frame$value <- as.numeric(as.character(SG_frame$value))
SW_frame <- melt(results_SW, id.vars = c("species"), 
                 variable.name = "wavelength",
                 value.name = "value")
SW_frame$wavelength <- as.numeric(as.character(SW_frame$wavelength))
SW_frame$value <- as.numeric(as.character(SW_frame$value))

if(sensor == "EMIT") {
  
  add_na <- data.table(species = unique_species,
                       wavelength = c(rep(1300, length(unique_species)),
                                      rep(1500, length(unique_species)),
                                      rep(1800, length(unique_species))),
                       value = rep(NA, length(unique_species)*3))
  
  SG_frame <- rbind(SG_frame, add_na)
  SW_frame <- rbind(SW_frame, add_na)
  
  SG_frame <- SG_frame[order(species, wavelength)]
  SW_frame <- SW_frame[order(species, wavelength)]
  
  band <- 805.36
  
} else {

  band <- 801.16

}

order_sp <- subset(SG_frame, wavelength == band) 
order_sp <- order_sp[order(value)]

SG_frame$species <- as.factor(SG_frame$species)
SG_frame$species <- factor(SG_frame$species, levels = order_sp$species)

SW_frame$species <- as.factor(SW_frame$species)
SW_frame$species <- factor(SW_frame$species, levels = order_sp$species)

wavelet <- ggplot(SW_frame, aes(x = wavelength, 
                                y = value, 
                                colour = species)) +
  geom_line() + 
  #scale_color_carto_d(type = "sequential", palette = "Peach", direction = -1) +
  scale_color_carto_d(type = "sequential", palette = "Magenta", direction = -1) +
  theme_bw(base_size = 12) + ylab("Wavelet transformed") + xlab("Wavelength (nm)") +
  #coord_cartesian(xlim = c(440, 900), ylim = c(-0.5, 0.5), expand = FALSE) +
  coord_cartesian(xlim = c(447, 2425), ylim = c(-0.5, 0.5), expand = FALSE) +
  scale_y_continuous(n.breaks = 3) +
  #scale_y_continuous(limits = c(-1.0, 1.0), expand = c(0, 0), n.breaks = 3) +
  dark_theme_bw(18) +
  theme(legend.position="none")

reflectance <- ggplot(SG_frame, aes(x = wavelength, 
                                    y = value, 
                                    colour = species)) +
  geom_line() + 
  #scale_color_carto_d(type = "sequential", palette = "Peach", direction = -1) +
  scale_color_carto_d(type = "sequential", palette = "Magenta", direction = -1) +
  theme_bw(base_size = 12) + ylab("Reflectance") + xlab(" ") +
  #coord_cartesian(xlim = c(440, 900), ylim = c(-0.01, 0.65), expand = FALSE) + 
  coord_cartesian(xlim = c(447, 2425), ylim = c(-0.01, 0.65), expand = FALSE) +
  scale_y_continuous(n.breaks = 3) +
  dark_theme_bw(18) +
  theme(legend.position="none")

spectra_plot <- ggarrange(reflectance, wavelet, 
                          ncol = 1, nrow = 2,
                          common.legend = FALSE)


png(paste0(root_path, "/figures/", sensor, "_spectra.png"), width = 210, height = 210, units = "mm", res = 600, pointsize = 22)
spectra_plot
dev.off()
