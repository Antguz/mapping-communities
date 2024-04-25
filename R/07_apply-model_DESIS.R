#'------------------------------------------------------------------------------
#' @title Apply PLS model to DESIS scenes
#'------------------------------------------------------------------------------

#' @description It applies the model to DESIS scenes.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(terra)
library(Rcpp)
library(RcppArmadillo)
library(parallel)

#'------------------------------------------------------------------------------
#' Source code

source("R/DESIS-transform.R")
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

apply_model <- function(scene_main_path, 
                        type, 
                        coefficients, 
                        output_path,
                        threads = 28) {

  #---------------------------------------------------------
  # Select and read products
  files <- list.files(path = scene_main_path, 
                      pattern = ".tif", 
                      recursive = TRUE,
                      include.dirs = TRUE)
  
  # Arrange path in a frame
  frame <- data.table(matrix(unlist(strsplit(files, "/")), 
                             nrow= length(files), 
                             byrow=TRUE), 
                      stringsAsFactors=FALSE)
  colnames(frame) <- c("scene")
  
  # Extract Product type
  frame[, type := unlist(strsplit(scene, "V"))[2], 
        by = seq_along(1:nrow(frame))]
  frame[, type := substr(type, 6, 30),
        by = seq_along(1:nrow(frame))]
  frame[, type := unlist(strsplit(type, ".tif"))[1], 
        by = seq_along(1:nrow(frame))]
  
  # Select products
  spectra <- paste0(scene_main_path, "/", frame[type == "SPECTRAL_IMAGE"]$scene[1])
  quality_cube <- paste0(scene_main_path, "/", frame[type == "QL_QUALITY"]$scene[1])
  quality_flags <- paste0(scene_main_path, "/", frame[type == "QL_QUALITY-2"]$scene[1])
  
  # Read products
  rast_spectra <- rast(spectra)
  names(rast_spectra) <- DESIS_bands
  rast_quality_cube <- rast(quality_cube)
  rast_quality_flags <- rast(quality_flags)
  
  #---------------------------------------------------------
  # Define mask
  
  # Mask for NDVI
  NDVI <- (rast_spectra$`801.16` - rast_spectra$`680.91`) / 
          (rast_spectra$`801.16` + rast_spectra$`680.91`)
  NDVI[NDVI == Inf] <- 0
  NDVI[NDVI >= 0.50] <- 2
  NDVI[NDVI < 0.50] <- 0
  NDVI[NDVI == 0] <- 1
  NDVI[NDVI == 2] <- 0
  
  # NIR mask
  NIR <- rast_spectra$`801.16`
  NIR[NIR < 3000] <- 1
  NIR[NIR >= 3000 & NIR <= 8000] <- 0
  NIR[NIR > 8000] <- 1
  
  # Quality flags
  mask_quality <- subset(rast_quality_flags, 1)
  names(mask_quality) <- "mask"
  mask_quality[] <- 0
  
  # Quality masks
  mask_quality <- mask(mask_quality, rast_quality_flags$Shadow, maskvalues = 1)
  #mask_quality <- mask(mask_quality, rast_quality_flags$Land, maskvalues = 0)
  mask_quality <- mask(mask_quality, rast_quality_flags$Snow, maskvalues = 1)
  #mask_quality <- mask(mask_quality, rast_quality_flags$`Haze (land)`, maskvalues = 1)
  mask_quality <- mask(mask_quality, rast_quality_flags$`Cloud (land)`, maskvalues = 1)
  mask_quality <- mask(mask_quality, rast_quality_flags$Water, maskvalues = 1)
  mask_quality <- mask(mask_quality, NDVI, maskvalues = 1)
  mask_quality <- mask(mask_quality, NIR, maskvalues = 1)
  
  # Mask spectra
  rast_spectra <- mask(rast_spectra, mask_quality, maskvalues = 0, inverse = TRUE)
  
  # As frame
  data_spectra <- as.data.table(as.data.frame(rast_spectra, cells = TRUE, na.rm = TRUE))
  colnames(data_spectra)[2:ncol(data_spectra)] <- as.character(DESIS_bands)
  
  # To predict
  predicted_scene <- subset(rast_spectra, 1:7)
  predicted_scene[] <- NA
  
  # Release memory
  rm(list = c('rast_spectra', 'mask_quality','NDVI', 'NIR', 'rast_quality_flags',
              'rast_quality_cube'))
  
  #---------------------------------------------------------
  # Transform raster files
  
  cell <- data_spectra$cell
  
  sw <- DESIS_transform(spectra = data_spectra[, .SD, .SDcols = 2:ncol(data_spectra)], 
                        DESIS_bands = DESIS_bands,
                        cell = cell,
                        new_bands = seq(min(DESIS_bands), max(DESIS_bands), by = 3.0),
                        threads = threads)
  
  # Subset bands
  cell <- sw$cell
  coef_bands <- colnames(coefficients)[4:ncol(coefficients)]
  sw <- sw[, ..coef_bands]
  
  # Release memory
  rm(list = c('data_spectra'))
  
  #--------------------------------------------------------
  # Predict raster files
  
  # As matrices
  Axis1 <- as.matrix(coefficients[Axes == "Axis1", .SD, .SDcols = 3:ncol(coefficients)])
  Axis2 <- as.matrix(coefficients[Axes == "Axis2", .SD, .SDcols = 3:ncol(coefficients)])
  Axis3 <- as.matrix(coefficients[Axes == "Axis3", .SD, .SDcols = 3:ncol(coefficients)])
  
  # Predict
  predicted_ordination <- predict_pls(as.matrix(sw),
                                      Axis1,
                                      Axis2,
                                      Axis3)
  
  predicted_ordination <- as.data.table(predicted_ordination)
  colnames(predicted_ordination) <- c("Axis1", "Axis2", "Axis3",
                                      "amp-Axis1", "amp-Axis2", "amp-Axis3", 
                                      "amp-total")
  
  #--------------------------------------------------------
  #Get a scene frame
  predicted_scene[cell] <- as.matrix(predicted_ordination)
  names(predicted_scene) <- c("Axis1", "Axis2", "Axis3",
                              "amp-Axis1", "amp-Axis2", "amp-Axis3", 
                              "amp-total")
  
  #-----------------------------------------------------------
  # Scene export
  
  scene_name <- strsplit(frame$scene[1], "DESIS-HSI-L2A-")[[1]][2]
  scene_name <- strsplit(scene_name, "-V")[[1]][1]
  scene_name <- paste0("DESIS_L3_", scene_name, "_", type, "_predicted.tif")
  
  output_raster <- paste0(output_path, "/", scene_name)
  
  #Export
  writeRaster(x = predicted_scene,
              filename = output_raster,
              names = c("Axis1", "Axis2", "Axis3",
                        "amp-Axis1", "amp-Axis2", "amp-Axis3", 
                        "amp-total"),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")

}

#'------------------------------------------------------------------------------
#' @example 

# Apply in a batch of scenes
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
scenes <- list.dirs(path = "/media/antonio/Work/Projects/DESIS_mapping/00_imagery/DESIS",
                    recursive = FALSE)
type <- "taxo"
coefficients <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/DESIS/DESIS_", type, "_coefficients.csv"))
output_path <- paste0(root_path, "/00_imagery/predicted/DESIS/", type)

pb <- txtProgressBar(min = 0, max = length(scenes), style = 3)

for(i in 1:length(scenes)) {
  
  setTxtProgressBar(pb, i)
  
  apply_model(scenes[i], 
              type,
              coefficients, 
              output_path,
              threads = 28)
  
}

# Apply in a batch of scenes
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
scenes <- list.dirs(path = "/media/antonio/Work/Projects/DESIS_mapping/00_imagery/DESIS",
                    recursive = FALSE)
type <- "phylo-log"
coefficients <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/DESIS/DESIS_", type, "_coefficients.csv"))
output_path <- paste0(root_path, "/00_imagery/predicted/DESIS/", type)

pb <- txtProgressBar(min = 0, max = length(scenes), style = 3)

for(i in 1:length(scenes)) {
  
  setTxtProgressBar(pb, i)
  
  apply_model(scenes[i], 
              type,
              coefficients, 
              output_path,
              threads = 28)
  
}

# Apply in a batch of scenes
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
scenes <- list.dirs(path = "/media/antonio/Work/Projects/DESIS_mapping/00_imagery/DESIS",
                    recursive = FALSE)
type <- "func"
coefficients <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/DESIS/DESIS_", type, "_coefficients.csv"))
output_path <- paste0(root_path, "/00_imagery/predicted/DESIS/", type)

pb <- txtProgressBar(min = 0, max = length(scenes), style = 3)

for(i in 1:length(scenes)) {
  
  setTxtProgressBar(pb, i)
  
  apply_model(scenes[i], 
              type,
              coefficients, 
              output_path,
              threads = 28)
  
}