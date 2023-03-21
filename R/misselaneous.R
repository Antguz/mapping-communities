#Plot DESIS spectra

spectral <- fread("/media/antonio/Extreme SSD/Hyperspectral/AGU/AGU/AGU.csv")
meta <- spectral[, c(1:7)]
spectral <- spectral[, c(13:242)]

smooth <- as.data.frame(spectral)

for(i in 1:nrow(smooth)) {
  
  sample <- as.numeric(smooth[i, ])
  filter <- sgolayfilt(sample, p = 1, n = 5, m = 0, ts = 1)
  smooth[i, ] <- filter
  
}

smooth <- as.data.table(smooth)
spectra <- smooth[, .SD, .SDcols = c(3:(ncol(smooth)-2))]
smooth <- as.data.frame(spectra)


frame <- cbind(meta[, 2], smooth)

frame_mean <- frame[, lapply(.SD, mean, na.rm = TRUE), by= SPCD]
frame_sd <- frame[, lapply(.SD, sd, na.rm = TRUE), by= SPCD]

frame_mean <- frame_mean[order(frame_mean$`811.6`)]
species <- frame_mean$SPCD

frame_mean <- melt(frame_mean, id.vars = c("SPCD"), variable.name = "Wavelength", value.name = "Reflectance")
frame_mean$Wavelength <- as.numeric(as.character(frame_mean$Wavelength))
frame_mean$SPCD <- as.factor(frame_mean$SPCD)
frame_mean$SPCD <- factor(frame_mean$SPCD, levels = species)

jpeg("Monocultures_DESIS.jpg", width = 210, height = 150, units = "mm", res = 300)
ggplot(frame_mean) +
  geom_line(aes(x = Wavelength, y = Reflectance/10000, colour = SPCD)) +
  theme_bw(base_size = 12) +
  scale_colour_viridis("Monocultures", option = "magma", discrete = TRUE) + ylab("Reflectance") + xlab("Wavelength (nm)")
dev.off()

#PLot model prediction


folder <- "/media/antonio/Extreme SSD/DESIS/00_imagery/North_2022/DESIS-HSI-L2A-DT0752621296 032-20220804T192231-V0215"

plsr1 <- readRDS("pls1.rds")
plsr2 <- readRDS("pls2.rds")
plsr3 <- readRDS("pls3.rds")

image <- rast(paste0(folder, "/", "DESIS-HSI-L2A-DT0752621296_032-20220804T192231-V0215-SPECTRAL_IMAGE.tif"))
bands
saveRDS(bands, "bands.rds")
bands <- readRDS("bands.rds")
names(image) <- bands
image <- subset(image, 10:233)

image_a2 <- image^2
image_sum <- sum(image_a2)
image_sqrt <- sqrt(image_sum)
image_normalize <- image/image_sqrt


axis1 <- predict(image_normalize, plsr1, cores= 1, cpkgs="pls", ncomp = 9)
axis2 <- predict(image_normalize, plsr2, cores= 1, cpkgs="pls", ncomp = 8)
axis3 <- predict(image_normalize, plsr3, cores= 1, cpkgs="pls", ncomp = 8)

predicted <- c(axis1, axis2, axis3)
names(predicted) <- c("PCoA1", "PCoA2", "PCoA3")
writeRaster(predicted, "model.tif", overwrite=TRUE)

#Get proportions

root_path <- "/media/antonio/Extreme SSD/DESIS"
fia <- fread(paste0(root_path, "/00_FIA/trees_WI_MN_20221102.csv"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"), na.strings= "NA")
out_path <- paste0(root_path, "/03_reshape")
SPCD <- na.exclude(SPCD)

fia <- subset(fia, STATECD == 27 | STATECD == 55)
fia <- fia[order(PLT_CN, SUBP)]

#Merge species code
fia <- merge(fia, SPCD[, c("SPCD", "PLANTS-Code", "A_G", "AM_EM")], by = "SPCD", all.x = TRUE, all.y = TRUE)

#Get Basal Area
fia$BA <- pi*(fia$DIA/2)^2

#Remove dead trees
fia <- fia[is.na(STANDING_DEAD_CD)]
fia$date <- as.IDate(paste0(fia$MEASYEAR, "-", fia$MEASMON, "-", fia$MEASDAY))

#Summarize
fia_summary <- fia[, .(ntrees = .N,
                       BA = sum(BA, na.rm = TRUE), 
                       HT = mean(HT, na.rm = TRUE), 
                       TPA_UNADJ = sum(TPA_UNADJ, na.rm = TRUE), 
                       DRYBIO_BG = sum(DRYBIO_BG, na.rm = TRUE),
                       CARBON_AG = sum(CARBON_AG, na.rm = TRUE), 
                       DRYBIO_AG = sum(DRYBIO_AG, na.rm = TRUE), 
                       CARBON_BG = sum(CARBON_BG, na.rm = TRUE)), 
                   by = c("PLOT", "PLT_CN", "A_G", "AM_EM")]
fia_summary <- fia_summary[order(PLOT, PLT_CN)]

unique_PLT_CN <- unique(fia_summary$PLT_CN)

results <- data.table()

for(i in 1:length(unique_PLT_CN)) {
  
  plot <- subset(fia_summary, PLT_CN == unique_PLT_CN[i])
  
  A <- as.numeric(sum(plot[A_G == "A", "BA"]))
  G <- as.numeric(sum(plot[A_G == "G", "BA"]))
  
  if(is.na(A)) {
    proportion_A <- 0
  } else if(is.na(G)) {
    proportion_A <- 1
  } else {
    Total <- A + G
    proportion_A <- round(A/Total, 2)
  }
  
  
  AM <- as.numeric(sum(plot[AM_EM == "AM", "BA"]))
  EM <- as.numeric(sum(plot[AM_EM == "EM", "BA"]))
  
  if(is.na(AM)) {
    proportion_AM <- 0
  } else if(is.na(EM)) {
    proportion_AM <- 1
  } else {
    Total_M <- AM + EM
    proportion_AM <- round(AM/Total_M, 2)
  }
  
  
  to_results <- data.table(PLT_CN = unique_PLT_CN[i],
                           A_to_T = proportion_A,
                           AM_to_T = proportion_AM)
  
  results <- rbind(results, to_results)
  
}

fwrite(results, "proportions.csv")


#

