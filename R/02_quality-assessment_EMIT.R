#'------------------------------------------------------------------------------
#' @title Quality assessment of spectra
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
root_path <- "/media/antonio/Extreme_Pro/Projects/MMDB/02_compile-files"
root_path <- "G:/Projects/MMBD/02_compile-files"

# Select sensor
sensor <- "EMIT"

# Read smoothed spectra (DESIS or EMIT)
spectra <- fread(paste0(root_path, "/02_", sensor, "_spectra-SG.txt"))

# Filter by flags
spectra <- spectra[`Aggregate Flag` == 1]
spectra_mean <- rowMeans(spectra[, 11:295])
spectra <- spectra[spectra_mean != -9999.00, ]
  
# Filter by NDVI
NDVI <- (spectra$`797.9` - spectra$`626.37`) / (spectra$`797.9` + spectra$`626.37`)
hist(NDVI)
spectra <- spectra[NDVI >= 0.4 & NDVI <= 1.0,]
NDVI <- (spectra$`797.9` - spectra$`626.37`) / (spectra$`797.9` + spectra$`626.37`)
hist(NDVI)

# Filter by NIR
NIR <- spectra$`797.9`
hist(NIR)
spectra <- spectra[NIR >= 0.1 & NIR <= 0.8,]
NIR <- spectra$`797.9`
hist(NIR)

# Filter by CCI
CCI <- (spectra$`529.53` - spectra$`626.37`) / (spectra$`529.53` + spectra$`626.37`)
hist(CCI)
spectra <- spectra[CCI >= 0.00 & CCI <= 0.6]
CCI <- (spectra$`529.53` - spectra$`626.37`) / (spectra$`529.53` + spectra$`626.37`)
hist(CCI)

# Filter by AOT_550
AOD550 <- spectra$AOD550
hist(AOD550)
spectra <- spectra[AOD550 <= 1,]
AOD550 <- spectra$AOD550
hist(AOD550)

# Filter by water contend
WV <- spectra$`H2O (g cm-2)`
hist(WV)

# Get metrics to export
NDVI <- (spectra$`797.9` - spectra$`626.37`) / (spectra$`797.9` + spectra$`626.37`)
NDVI <- round(NDVI, 2)
CCI <- (spectra$`529.53` - spectra$`626.37`) / (spectra$`529.53` + spectra$`626.37`)
CCI <- round(CCI, 2)
NIR <- round(spectra$`797.9`, 2)

# Final file
master <- spectra[, c("scene", "PLT_CN")]
master$NDVI <- NDVI
master$NIR <- NIR
master$CCI <- CCI

# Export
fwrite(master, paste0(root_path, "/EMIT_master_spectra.csv"))
