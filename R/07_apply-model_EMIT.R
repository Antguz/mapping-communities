#'------------------------------------------------------------------------------
#' @title Apply PLS model to EMIT scenes
#'------------------------------------------------------------------------------

#' @description It applies the model to EMIT scenes.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(terra)
library(Rcpp)
library(RcppArmadillo)
library(parallel)

#'------------------------------------------------------------------------------
#' Source code

source("R/EMIT-transform.R")
sourceCpp("src/model-application.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param scene_main_path
#' @param coefficients
#' @param output_raster
#' @param threads

#'------------------------------------------------------------------------------
#' Function

apply_model <- function(scene_path,
                        scene_frame, 
                        coefficients, 
                        output_path,
                        threads = 28) {
  
  #---------------------------------------------------------
  # Select and read products
  
  # Select products
  spectra <- paste0(scene_path, "/", 
                    scene_frame[type == "RFL"]$scene)
  
  mask <- paste0(scene_path, "/", 
                 scene_frame[type == "MASK" & band_mask == "mask"]$scene)
  
  # Read products
  rast_spectra <- rast(spectra)
  names(rast_spectra) <- EMIT_bands
  rast_mask <- rast(mask)
  
  #---------------------------------------------------------
  # Mask reflectance
  
  # Mask for quality
  rast_spectra <- mask(rast_spectra, 
                       rast_mask$EMIT_L2A_MASK_001_20230527T194414_2314713_001_mask_8,
                       maskvalues = 1)
  
  # Mask for NDVI
  NDVI <- (rast_spectra$`805.361755371093` - rast_spectra$`686.010314941406`) / 
          (rast_spectra$`805.361755371093` + rast_spectra$`686.010314941406`)
  NDVI[NDVI == Inf] <- NA
  NDVI[NDVI < 0.50] <- 0
  NDVI[NDVI >= 0.50] <- 1
  rast_spectra <- mask(rast_spectra, NDVI, maskvalues = 1, inverse = TRUE)
  
  # Mask for NIR
  NIR <- rast_spectra$`805.361755371093`
  NIR[NIR == Inf] <- NA
  NIR[NIR < 0.3] <- 0
  NIR[NIR >= 0.3 & NIR <= 0.8] <- 1
  NIR[NIR > 0.8] <- 1
  rast_spectra <- mask(rast_spectra, NIR, maskvalues = 1, inverse = TRUE)
  
  # As frame
  data_spectra <- as.data.table(as.data.frame(rast_spectra, cells = TRUE, na.rm = TRUE))
  colnames(data_spectra)[2:ncol(data_spectra)] <- as.character(EMIT_bands)
  
  # Release memory
  rm(list = c('rast_spectra', 'rast_mask','NDVI', 'NIR'))
  
  #---------------------------------------------------------
  # Transform raster files
  
  cell <- data_spectra$cell
  spectra <- data_spectra[, .SD, .SDcols = 2:ncol(data_spectra)]
  
  sw <- EMIT_transform(spectra = data_spectra[, .SD, .SDcols = 2:ncol(data_spectra)], 
                       bands = EMIT_bands,
                       cell = cell,
                       new_bands = seq(min(EMIT_bands), max(EMIT_bands), by = 7.43),
                       threads = threads)
  
  cell <- sw$cell
  coef_bands <- colnames(coefficients)[4:ncol(coefficients)]
  sw <- sw[, ..coef_bands]
  
  #--------------------------------------------------------
  # Predict raster files
  
  # As matrices
  spectra <- as.matrix(sw)
  MDS1 <- as.matrix(coefficients[MDS == "1", .SD, .SDcols = 3:ncol(coefficients)])
  MDS2 <- as.matrix(coefficients[MDS == "2", .SD, .SDcols = 3:ncol(coefficients)])
  MDS3 <- as.matrix(coefficients[MDS == "3", .SD, .SDcols = 3:ncol(coefficients)])
  
  # Predict
  predicted_ordination <- predict_pls(spectra,
                                      MDS1,
                                      MDS2,
                                      MDS3)
  
  predicted_ordination <- as.data.table(predicted_ordination)
  colnames(predicted_ordination) <- c("MDS1", "MDS2", "MDS3",
                                      "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                                      "amp-total")
  
  #--------------------------------------------------------
  #Get a scene frame
  predicted_scene <- subset(rast_spectra, 1:7)
  predicted_scene[cell] <- as.matrix(predicted_ordination)
  names(predicted_scene) <- c("MDS1", "MDS2", "MDS3",
                              "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                              "amp-total")
  
  #-----------------------------------------------------------
  # Scene export
  
  #Export
  writeRaster(x = predicted_scene,
              filename = output_raster,
              names = c("MDS1", "MDS2", "MDS3",
                        "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                        "amp-total"),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
}

#'------------------------------------------------------------------------------
#' @example 

# Scenes path
#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "E:/Projects/DESIS_mapping"
scene_path <- paste0(root_path, "/00_imagery/EMIT", "")

# List scenes
files <- list.files(scene_path)

# Arrange path in a frame
frame <- data.table(scene = files)

# Get product information
frame[, type := strsplit(scene, "_")[[1]][3], by = seq_along(1:nrow(frame))]
frame[, time := strsplit(scene, "_")[[1]][5], by = seq_along(1:nrow(frame))]
frame[, band_mask := strsplit(scene, "_")[[1]][8], by = seq_along(1:nrow(frame))]
frame[, keep := unlist(gregexpr('.hdr', scene))[1], by = seq_along(1:nrow(frame))]
frame <- frame[keep == -1, 1:4]

list_scenes <- data.table(time = unique(frame$time))
list_scenes$ID <- 1:nrow(list_scenes)
frame <- merge(frame, list_scenes, by = "time", all.x = TRUE, all.y = TRUE)

output_raster <- paste0(scene_main_path, "_predicted.tif")
threads <- 28

# Apply in a single scene
apply_model(scene_path,
            scene_frame = frame, 
            coefficients, 
            output_path,
            threads = 28)

# Apply in a batch of scenes
root_path <- "/media/antonio/Work/Projects/EMIT_mapping"
scenes <- list.dirs(path = paste0(root_path, "/00_imagery"),
                    recursive = FALSE)
coefficients <- fread(paste0(root_path, "/06_model/SW_beta_PLS/EMIT/EMIT_", "coefficients.csv"))
output_scenes <- paste0(scenes, "_predicted.tif")

for(i in 1:length(scenes)) {
  
  apply_model(scenes[i], 
              coefficients, 
              output_scenes[i],
              threads = 28)
  
}
