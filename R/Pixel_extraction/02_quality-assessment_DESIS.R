#'------------------------------------------------------------------------------
#' @title Quality assessment of spectra for DESIS
#'------------------------------------------------------------------------------

#' @description It helps to do a quality assessment of the spectra extracted
#' 
#' @return A .txt file (master file) with the spectral observations that can be 
#' used.

#'------------------------------------------------------------------------------
#' @Library

library(data.table)
library(ggplot2)
library(outliers)
library(Rcpp)
library(RcppArmadillo)

#'------------------------------------------------------------------------------
#' @Source-code
sourceCpp("src/sam.cpp")

#'------------------------------------------------------------------------------
#' @Cleaning-of-potential-outliers

# Select root path
root_path <- "/media/antonio/Work/MMDB/Pixel_extraction/02_compile-files"

# Select sensor
sensor <- "DESIS"

# Read smoothed spectra (DESIS or EMIT)
spectra <- fread(paste0(root_path, "/02_", sensor, "_spectra-SG.txt"))

# Filter by flags
spectra <- spectra[Shadow != 1]
spectra <- spectra[Snow != 1]
spectra <- spectra[`Haze (water)` != 1]
spectra <- spectra[`Cloud (land)` != 1]
spectra <- spectra[`Cloud (water)` != 1]
spectra <- spectra[Water != 1]

# Filter by NDVI
NDVI <- (spectra$`801.16` - spectra$`629.66`) / (spectra$`801.16` + spectra$`629.66`)
hist(NDVI)
spectra <- spectra[NDVI >= 0.4 & NDVI <= 1.0,]
NDVI <- (spectra$`801.16` - spectra$`629.66`) / (spectra$`801.16` + spectra$`629.66`)
hist(NDVI)

# Filter by NIR
NIR <- spectra$`801.16`
hist(NIR)
spectra <- spectra[NIR >= 1000 & NIR <= 8000,]
NIR <- spectra$`801.16`
hist(NIR)

# Filter by EVI
EVI <- 2.5*((spectra$`801.16` - spectra$`629.66`)/(spectra$`801.16` + 6 * spectra$`629.66` - 7.5 * spectra$`450.34` + 1))
hist(EVI)
spectra <- spectra[EVI >= quantile(EVI, 0.05) & EVI <= quantile(EVI, 0.95),]
EVI <- 2.5*((spectra$`801.16` - spectra$`629.66`)/(spectra$`801.16` + 6 * spectra$`629.66` - 7.5 * spectra$`450.34` + 1))
hist(EVI)

# Filter by AOT_550
AOT_550 <- spectra$`AOT (550 nm) * 100`
hist(AOT_550)
spectra <- spectra[AOT_550 <= 100,]
AOT_550 <- spectra$`AOT (550 nm) * 100`
hist(AOT_550)

# Filter by water contend
WV <- spectra$`WV [cm] * 42`
hist(WV)

# Filter using SMA
ref_sam <- sam_ref(as.matrix(spectra[, 33:209]), #450 - 900
                   colMeans(as.matrix(spectra[, 33:209])),
                   threads = 12)[,1]

hist(ref_sam)
plot(as.matrix(spectra[which.max(ref_sam), .SD, .SDcols = c(33:209)])[1,])

spectra <- spectra[ref_sam < quantile(ref_sam, 0.99), ]
ref_sam <- ref_sam[ref_sam < quantile(ref_sam, 0.99)]

# Filter using SMA in the visible
ref_sam_vis <- sam_ref(as.matrix(spectra[, .SD, .SDcols = c(33:120)]),
                       colMeans(as.matrix(spectra[, .SD, .SDcols = c(33:120)])),
                       threads = 12)[,1]
hist(ref_sam_vis)
plot(names(colMeans(as.matrix(spectra[, .SD, .SDcols = c(33:120)]))),
     as.matrix(spectra[which.max(ref_sam_vis), .SD, .SDcols = c(33:120)])[1,])

spectra <- spectra[ref_sam_vis < quantile(ref_sam_vis, 0.95), ]
ref_sam_vis <- ref_sam_vis[ref_sam_vis < quantile(ref_sam_vis, 0.95)]
ref_sam <- ref_sam[ref_sam_vis < quantile(ref_sam_vis, 0.95)]

# Get metrics to export
master <- spectra[, c("folder", "scene", "PLT_CN", "AOT (550 nm) * 100", "WV [cm] * 42")]
master$NDVI <- (spectra$`801.16` - spectra$`629.66`) / (spectra$`801.16` + spectra$`629.66`)
master$NIR <- spectra$`801.16`
master$EVI <- EVI <- 2.5*((spectra$`801.16` - spectra$`629.66`)/(spectra$`801.16` + 6 * spectra$`629.66` - 7.5 * spectra$`450.34` + 1))
master$CCI <- (spectra$`532.33` - spectra$`629.66`) / (spectra$`532.33` + spectra$`629.66`)
master$CRI700 <- (1/spectra$`509.54`) - (1/spectra$`701.41`) #(1/R510) - (1/R700)
master$VOG2 <- (spectra$`734.37` - spectra$`747.13`)/(spectra$`714.02` + spectra$`726.85`)  #(R734 - R747)/(R715 + R726)
master$NPCI <- (spectra$`680.91` - spectra$`429.86`)/(spectra$`680.91` + spectra$`429.86`) #(R680 - R430)/(R680 + R430)
master$WBI <- spectra$`970.43`/spectra$`901.2` # R970/R900
master$PRI570 <- (spectra$`571` - spectra$`532.33`) / (spectra$`571` + spectra$`532.33`) #(R570 - R531)/(R570 + R531)
master$CMS <- (spectra$`693.47`/spectra$`760.33`) #R694/R760
master$RDVI <- (spectra$`801.16` - spectra$`670.48`)/ sqrt(spectra$`801.16` + spectra$`670.48`) #(R800 - R670)/sqrt(R800 + R670)

master <- master[!is.infinite(NPCI)]

# Remove points based on mahalanobis distance
mahalanobis_dist <- mahalanobis(master[, 6:16], colMeans(master[, 6:16]), cov(master[, 6:16]))
hist(mahalanobis_dist)
threshold <- qchisq(0.95, df = ncol(master[, 6:16]))
hist(mahalanobis_dist[mahalanobis_dist <= threshold])
outliers <- mahalanobis_dist <= threshold
master <- master[outliers == TRUE, ]

# Export
fwrite(master, paste0(root_path, "/DESIS_master_spectra.csv"))
