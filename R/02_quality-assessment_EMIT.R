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
library(Rcpp)
library(RcppArmadillo)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Source code
sourceCpp("src/sam.cpp")

#'------------------------------------------------------------------------------
#' Arguments

# Select root path
root_path <- "/media/antonio/Work/MMDB/02_compile-files"
#root_path <- "G:/Projects/MMBD/02_compile-files"

# Select sensor
sensor <- "EMIT"

# Read smoothed spectra (DESIS or EMIT)
spectra <- fread(paste0(root_path, "/02_", sensor, "_spectra-SG.txt"))

# Filter by flags
spectra <- spectra[`Cloud flag` != 1,]
spectra <- spectra[`Cirrus flag` != 1,]
spectra <- spectra[`Water flag` != 1,]
spectra <- spectra[`Spacecraft Flag` != 1,]
#spectra <- spectra[`Aggregate Flag` == 1]
#spectra <- spectra[`Aggregate Flag` == 1]

# Remove date
scene <- substr(spectra$scene, 1, 44)
date <- substr(scene, 22, 23)
spectra <- spectra[date != "05"]

# Remove NA pixels
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
#CCI <- (spectra$`529.53` - spectra$`626.37`) / (spectra$`529.53` + spectra$`626.37`)
#hist(CCI)
#spectra <- spectra[CCI >= 0.00 & CCI <= 0.7]
#CCI <- (spectra$`529.53` - spectra$`626.37`) / (spectra$`529.53` + spectra$`626.37`)
#hist(CCI)

# Filter by AOT_550
AOD550 <- spectra$AOD550
hist(AOD550)
spectra <- spectra[AOD550 <= 1,]
AOD550 <- spectra$AOD550
hist(AOD550)

# Filter by water contend
WV <- spectra$`H2O (g cm-2)`
hist(WV)

# Filter using SMA
ref_sam <- sam_ref(as.matrix(spectra[, .SD, .SDcols = c(20:135,157:195,226:295)]),
                   colMeans(as.matrix(spectra[, .SD, .SDcols = c(20:135,157:195,226:295)])),
                   threads = 12)[,1]
hist(ref_sam)
plot(as.matrix(spectra[which.max(ref_sam), .SD, .SDcols = c(20:135,156:195,226:295)])[1,])

spectra <- spectra[ref_sam < quantile(ref_sam, 0.95), ]
ref_sam <- ref_sam[ref_sam < quantile(ref_sam, 0.95)]

# Filter using SMA in the visible
ref_sam_vis <- sam_ref(as.matrix(spectra[, .SD, .SDcols = c(20:50)]),
                       colMeans(as.matrix(spectra[, .SD, .SDcols = c(20:50)])),
                       threads = 12)[,1]
hist(ref_sam_vis)
plot(as.matrix(spectra[which.max(ref_sam_vis), .SD, .SDcols = c(20:50)])[1,])

spectra <- spectra[ref_sam_vis < quantile(ref_sam_vis, 0.999), ]
ref_sam_vis <- ref_sam_vis[ref_sam_vis < quantile(ref_sam_vis, 0.999)]
ref_sam <- ref_sam[ref_sam_vis < quantile(ref_sam_vis, 0.999)]

# Get metrics to export
master <- spectra[, c("scene", "PLT_CN", "AOD550", "H2O (g cm-2)")]
master$NDVI <- (spectra$`797.9` - spectra$`626.37`) / (spectra$`797.9` + spectra$`626.37`)
master$CCI <- (spectra$`529.53` - spectra$`626.37`) / (spectra$`529.53` + spectra$`626.37`)
master$NIR <- spectra$`797.9`
master$CRI700 <- (1/spectra$`507.21`) - (1/spectra$`700.93`) #(1/R510) - (1/R700)
master$VOG2 <- (spectra$`738.22` - spectra$`745.68`)/(spectra$`715.84` + spectra$`723.3`)  #(R734 - R747)/(R715 + R726)
master$NPCI <- (spectra$`678.55` - spectra$`432.89`)/(spectra$`678.55` + spectra$`432.89`) #(R680 - R430)/(R680 + R430)
master$WBI <- spectra$`969.53`/spectra$`902.37`# R970/R900
master$PRI570 <- (spectra$`574.21` - spectra$`529.53`) / (spectra$`574.21` + spectra$`529.53`) #(R570 - R531)/(R570 + R531)
master$CMS <- (spectra$`693.47`/spectra$`760.6`) #R694/R760
master$RDVI <- (spectra$`797.9` - spectra$`671.1`)/ sqrt(spectra$`797.9` + spectra$`671.1`) #(R800 - R670)/sqrt(R800 + R670)
  
# Remove points based on mahalanobis distance
mahalanobis_dist <- mahalanobis(master[, 5:14], colMeans(master[, 5:14]), cov(master[, 5:14]))
hist(mahalanobis_dist)
threshold <- qchisq(0.95, df = ncol(master[, 5:14]))
hist(mahalanobis_dist[mahalanobis_dist <= threshold])
outliers <- mahalanobis_dist <= threshold
master <- master[outliers == TRUE, ]

# Export
fwrite(master, paste0(root_path, "/EMIT_master_spectra.csv"))
