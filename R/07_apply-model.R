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

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"
scene_main_path <- paste0(root_path, "/00_imagery/", "DESIS-HSI-L2A-DT0465807416_009-20200617T193003-V0217")
coefficients <- fread(paste0(root_path, "/06_model/SW_beta_PLS/", "coefficients.csv"))
output_raster <- paste0(scene_main_path, "_predicted.tif")
threads <- 28

apply_model(scene_main_path, 
            coefficients, 
            output_raster,
            threads = 28)

#'------------------------------------------------------------------------------
#' Function

apply_model <- function(scene_main_path, 
                        coefficients, 
                        output_raster,
                        threads = 28) {
  
  #Colnames
  bands <- colnames(coefficients)
  bands <- bands[4:length(bands)]
  
  #---------------------------------------------------------
  # Raster files processing
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
  frame[, type := unlist(strsplit(scene, "V0217-"))[2], 
        by = seq_along(1:nrow(frame))]
  #frame[, type := unlist(strsplit(type, "-"))[2], 
  #      by = seq_along(1:nrow(frame))]
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
  
  # As frame
  data_spectra <- as.data.table(as.data.frame(rast_spectra, cells = TRUE, na.rm = TRUE))
  colnames(data_spectra)[2:ncol(data_spectra)] <- as.character(DESIS_bands)
  #data_quality_cube <- as.data.table(as.data.frame(rast_quality_cube, cells = TRUE))
  #colnames(data_quality_cube)[-1] <- as.character(DESIS_bands, cells = TRUE)
  data_quality_flags <- as.data.table(as.data.frame(rast_quality_flags, cells = TRUE, na.rm = TRUE))
  
  #---------------------------------------------------------
  # Transform raster files
  
  cell <- data_spectra$cell
  spectra <- data_spectra[, .SD, .SDcols = 2:ncol(data_spectra)]
  
  sw <- desis_transform(spectra = spectra, 
                        DESIS_bands = DESIS_bands,
                        cell = cell,
                        threads = threads)
  
  cell <- sw$cell
  sw <- sw[, .SD, .SDcols = 25:195]
  colnames(sw) <- bands
  
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
  # Select products
  mask_quality <- subset(rast_quality_flags, 1)
  names(mask_quality) <- "mask"
  mask_quality[] <- 1
  
  #Masks
  mask_quality <- mask(mask_quality, rast_quality_flags$Shadow, maskvalues = 1)
  mask_quality <- mask(mask_quality, rast_quality_flags$Land, maskvalues = 0)
  mask_quality <- mask(mask_quality, rast_quality_flags$Snow, maskvalues = 1)
  mask_quality <- mask(mask_quality, rast_quality_flags$`Haze (land)`, maskvalues = 1)
  mask_quality <- mask(mask_quality, rast_quality_flags$`Haze (water)`, maskvalues = 1)
  mask_quality <- mask(mask_quality, rast_quality_flags$`Cloud (land)`, maskvalues = 1)
  mask_quality <- mask(mask_quality, rast_quality_flags$`Cloud (water)`, maskvalues = 1)
  mask_quality <- mask(mask_quality, rast_quality_flags$Water, maskvalues = 1)
  mask_quality <- mask(mask_quality, rast_quality_flags$Water, maskvalues = 1)
  NDVI <- (rast_spectra$`801` - rast_spectra$`680.8`) / (rast_spectra$`801` + rast_spectra$`680.8`)
  NDVI[NDVI == Inf] <- NA
  NDVI[NDVI < 0.60] <- 0
  NDVI[NDVI >= 0.60] <- 1
  
  
  #Final mask
  predicted_scene <- mask(predicted_scene, mask_quality, maskvalues = 1, inverse = TRUE)
  predicted_scene <- mask(predicted_scene, NDVI, maskvalues = 1, inverse = TRUE)
  
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
