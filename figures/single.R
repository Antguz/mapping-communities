spectra <- fread(paste0(root_path, "/04_match_spectra-stand/spectra.csv"))
spectra_SG <- fread(paste0(root_path, "/04_match_spectra-stand/spectra_SG.csv"))
species <- fread(paste0(root_path, "/04_match_spectra-stand/species.csv"))

spectra <- spectra[, .SD, .SDcols = c(1, 7:ncol(spectra))]
spectra_SG <- spectra_SG[, .SD, .SDcols = c(1, 7:ncol(spectra_SG))]

inventory <- species[, .SD, .SDcols = c(7:ncol(species))]

nspecies <- NA

for(i in 1:nrow(inventory)) {
  
  sub <- as.numeric(inventory[i, ])
  sub <- sub[sub > 0]
  nspecies[i] <- length(sub)

}

species$n <- nspecies 

sub_species <- subset(species, n == 1)
inventory <- sub_species[, .SD, .SDcols = c(7:(ncol(sub_species)-1))]

sp_names <- colnames(inventory)
sp <- NA

for(i in 1:nrow(inventory)) {
  
  sub <- as.numeric(inventory[i, ])
  sub <- sp_names[sub > 0]
  sp[i] <- sub
  
}

sub_species <- sub_species[, 1:6]
sub_species$sp <- sp

frame <- merge(sub_species, spectra, by = c("PLT_CN"), all.x = TRUE, all.y = FALSE)
frame  <- frame[, .SD, .SDcols = c(7:ncol(frame))]

melt_frame <- melt(frame, id.vars = c("sp"), 
                 variable.name = "wavelength",
                 value.name = "value")
frame_summary <- summarySE(data = melt_frame, measurevar = "value", groupvars = c("sp", "wavelength"),
                         na.rm = TRUE, conf.interval = 0.95, .drop = TRUE)
frame_summary$wavelength <- as.numeric(as.character(frame_summary$wavelength))
frame_summary <- as.data.table(frame_summary)

for_order <- subset(frame_summary, wavelength == 746.00)
for_order <- for_order[order(-value)]
for_order <- for_order$sp

frame_summary$sp <- as.factor(frame_summary$sp)
frame_summary$sp <- factor(frame_summary$sp, levels = for_order)

wavelet <- ggplot(frame_summary, aes(x = wavelength, y = value/10000, colour = sp)) +
  geom_line() + 
  scale_color_viridis(discrete = TRUE, option = "A", direction = -1) +
  theme_bw(base_size = 12) + ylab("Wavelet") + xlab("Wavelength (nm)") +
  scale_x_continuous(limits = c(460, 900), expand = c(0, 0)) 



frame <- merge(sub_species, spectra_SG, by = c("PLT_CN"), all.x = TRUE, all.y = FALSE)
frame  <- frame[, .SD, .SDcols = c(7:ncol(frame))]

melt_frame <- melt(frame, id.vars = c("sp"), 
                   variable.name = "wavelength",
                   value.name = "value")
frame_summary <- summarySE(data = melt_frame, measurevar = "value", groupvars = c("sp", "wavelength"),
                           na.rm = TRUE, conf.interval = 0.95, .drop = TRUE)
frame_summary$wavelength <- as.numeric(as.character(frame_summary$wavelength))
frame_summary <- as.data.table(frame_summary)

frame_summary$sp <- as.factor(frame_summary$sp)
frame_summary$sp <- factor(frame_summary$sp, levels = for_order)

reflectance <- ggplot(frame_summary, aes(x = wavelength, y = value/10000, colour = sp)) +
  geom_line() + 
  scale_color_viridis("Monocultures", discrete = TRUE, option = "A", direction = -1) +
  theme_bw(base_size = 12) + ylab("Reflectance") + xlab("Wavelength (nm)") +
  scale_x_continuous(limits = c(460, 900), expand = c(0, 0)) 
  #scale_x_continuous(limits = c(460, 900), expand = c(0, 0)) 

jpeg("figures/spectra.jpg", width = 210, height = 210, units = "mm", res = 300)
ggarrange(reflectance, wavelet, 
          labels = c("a", "b"),
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "right")
dev.off()
