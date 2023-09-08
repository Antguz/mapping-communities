#'------------------------------------------------------------------------------
#' @title Quality assessment of spectra for DESIS
#'------------------------------------------------------------------------------

#' @description It helps to do a quality assessment of the spectra extracted
#' 
#' @return A .txt file (master file) with the spectral observations that can be 
#' used.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(ggplot2)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments

# Select root path
root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping/02_compile-files"

# Select sensor
sensor <- "DESIS"

# Read smoothed spectra (DESIS or EMIT)
spectra <- fread(paste0(root_path, "/02_", sensor, "_spectra-SG.txt"))

# Filter by flags
spectra <- spectra[Shadow != 1]
spectra <- spectra[`Haze (land)` != 1]
spectra <- spectra[`Cloud (land)` != 1]
spectra <- spectra[Water != 1]

# Filter by NDVI
NDVI <- (spectra$`801.16` - spectra$`629.66`) / (spectra$`801.16` + spectra$`629.66`)
hist(NDVI)
spectra <- spectra[NDVI >= 0.6 & NDVI <= 1.0,]
NDVI <- (spectra$`801.16` - spectra$`629.66`) / (spectra$`801.16` + spectra$`629.66`)
hist(NDVI)

# Filter by NIR
NIR <- spectra$`801.16`
hist(NIR)
spectra <- spectra[NIR >= 1500 & NIR <= 8000,]
NIR <- spectra$`801.16`
hist(NIR)

# Filter by CCI
CCI <- (spectra$`532.33` - spectra$`629.66`) / (spectra$`532.33` + spectra$`629.66`)
hist(CCI)
spectra <- spectra[CCI >= 0.05 & CCI <= 1.0]
CCI <- (spectra$`532.33` - spectra$`629.66`) / (spectra$`532.33` + spectra$`629.66`)
hist(CCI)

# Filter by AOT_550
AOT_550 <- spectra$`AOT (550 nm) * 100`
hist(AOT_550)
spectra <- spectra[AOT_550 <= 100,]
AOT_550 <- spectra$`AOT (550 nm) * 100`
hist(AOT_550)

# Filter by water contend
WV <- spectra$`WV [cm] * 42`
hist(WV)

# Get metrics to export
NDVI <- (spectra$`801.16` - spectra$`629.66`) / (spectra$`801.16` + spectra$`629.66`)
NDVI <- round(NDVI, 2)
CCI <- (spectra$`532.33` - spectra$`629.66`) / (spectra$`532.33` + spectra$`629.66`)
CCI <- round(CCI, 2)
NIR <- round(spectra$`801.16`, 2)

# Final file
master <- spectra[, c("folder", "scene", "PLT_CN")]
master$NDVI <- NDVI
master$NIR <- NIR
master$CCI <- CCI

# Export
fwrite(master, paste0(root_path, "/DESIS_master_spectra.csv"))
