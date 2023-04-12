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

#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "E:/Projects/DESIS_mapping"
scene_main_path <- paste0(root_path, "/00_imagery/", "DESIS-HSI-L2A-DT0599346912 015-20210615T161539-V0215")
coefficients <- fread(paste0(root_path, "/06_model/SW_beta_PLS/", "coefficients.csv"))
output_raster <- paste0(scene_main_path, "_predicted.tif")

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
  frame[, type := unlist(strsplit(scene, "V0215-"))[2], 
        by = seq_along(1:nrow(frame))]
  frame[, type := unlist(strsplit(type, ".tif"))[1], 
        by = seq_along(1:nrow(frame))]
  
  # Select products
  spectra <- paste0(scene_main_path, "/", frame[type == "SPECTRAL_IMAGE"]$scene[1])
  quality_cube <- paste0(scene_main_path, "/", frame[type == "QL_QUALITY"]$scene[1])
  quality_flags <- paste0(scene_main_path, "/", frame[type == "QL_QUALITY-2"]$scene[1])
  
  # Read products
  rast_spectra <- rast(spectra)
  rast_quality_cube <- rast(quality_cube)
  rast_quality_flags <- rast(quality_flags)
  
  # As frame
  data_spectra <- as.data.table(as.data.frame(rast_spectra, cells = TRUE))
  colnames(data_spectra)[-1] <- as.character(DESIS_bands)
  data_quality_cube <- as.data.table(as.data.frame(rast_quality_cube, cells = TRUE))
  colnames(data_quality_cube)[-1] <- as.character(DESIS_bands, cells = TRUE)
  data_quality_flags <- as.data.table(as.data.frame(rast_quality_flags, cells = TRUE))
  
  #---------------------------------------------------------
  # Transform raster files
  
  sw <- desis_transform(scene_main_path, 
                        DESIS_bands,
                        threads = threads)
  
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
  #-----------------------------------------------------------
  # Select products
  quality_flags <- paste0(scene_main_path, "/", frame[type == "QL_QUALITY-2"]$scene[1])
  scene_flag <- rast(quality_flags)
  scene_frame <- as.data.table(as.data.frame(scene_flag, cells = TRUE))
  
  #Quality flags
  scene_frame <- scene_frame[Shadow != 1, kepp := 1]
  scene_frame <- scene_frame[`Haze (land)` != 1, kepp := 1]
  scene_frame <- scene_frame[`Cloud (land)` != 1, kepp := 1]
  scene_frame <- scene_frame[Water != 1, kepp := 1]
  scene_frame <- scene_frame[kepp != 1, kepp := 0]
  
  #-----------------------------------------------------------
  # Merge
  predicted_ordination$cell <- scene_frame[, cell]
  
  #Get an scene of predictions
  scene_predictions <- c(scene_flag$`AOT (550 nm) * 100`,
                         scene_flag$`AOT (550 nm) * 100`,
                         scene_flag$`AOT (550 nm) * 100`,
                         scene_flag$`AOT (550 nm) * 100`,
                         scene_flag$`AOT (550 nm) * 100`,
                         scene_flag$`AOT (550 nm) * 100`,
                         scene_flag$`AOT (550 nm) * 100`)
  names(scene_predictions) <- c("MDS1", "MDS2", "MDS3",
                                "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                                "amp-total")
  scene_predictions <- 
    scene_predictions[] <- as.matrix(scene_frame)
    
}