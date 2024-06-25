#'------------------------------------------------------------------------------
#' @title Apply PLS model to DESIS scenes (MSI)
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
                        coefficients_taxonomic,
                        coefficients_phylogenetic,
                        coefficients_functional,
                        output_path,
                        temp,
                        threads = 28) {
  
  #---------------------------------------------------------
  # Unzip files
  name_product <- strsplit(scene_main_path, "/")[[1]][length(strsplit(scene_main_path, "/")[[1]])]
  name_product <- strsplit(name_product, ".zip")[[1]]
  scene_path <- paste0(temp, "/", name_product)
  unzip(scene_main_path, exdir = paste0(temp, "/", name_product), overwrite = TRUE)
  
  #---------------------------------------------------------
  # Select and read products
  files <- list.files(path = scene_path, 
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
  spectra <- paste0(scene_path, "/", frame[type == "SPECTRAL_IMAGE"]$scene[1])
  quality_cube <- paste0(scene_path, "/", frame[type == "QL_QUALITY"]$scene[1])
  quality_flags <- paste0(scene_path, "/", frame[type == "QL_QUALITY-2"]$scene[1])
  
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
  mask_quality <- mask(mask_quality, rast_quality_flags$Land, maskvalues = 0)
  mask_quality <- mask(mask_quality, rast_quality_flags$Snow, maskvalues = 1)
  mask_quality <- mask(mask_quality, rast_quality_flags$`Haze (land)`, maskvalues = 1)
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
  
  # Final names
  scene_name <- strsplit(name_product, "DESIS-HSI-L2A-")[[1]][2]
  tscene_name <- paste0("DESIS_L3_", scene_name, "_taxonomic-composition.tif")
  pscene_name <- paste0("DESIS_L3_", scene_name, "_phylogenetic-composition.tif")
  fscene_name <- paste0("DESIS_L3_", scene_name, "_functional-composition.tif")
  
  ### TAXONOMIC ###
  # As matrices
  tMDS1 <- as.matrix(coefficients_taxonomic[MDS == "1", .SD, .SDcols = 3:ncol(coefficients_taxonomic)])
  tMDS2 <- as.matrix(coefficients_taxonomic[MDS == "2", .SD, .SDcols = 3:ncol(coefficients_taxonomic)])
  tMDS3 <- as.matrix(coefficients_taxonomic[MDS == "3", .SD, .SDcols = 3:ncol(coefficients_taxonomic)])
  
  # Predict
  tpredicted_ordination <- predict_pls(as.matrix(sw),
                                      tMDS1,
                                      tMDS2,
                                      tMDS3)
  
  tpredicted_ordination <- as.data.table(tpredicted_ordination)
  colnames(tpredicted_ordination) <- c("MDS1", "MDS2", "MDS3",
                                      "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                                      "amp-total")
  
  #--------------------------------------------------------
  #Get a scene frame
  tpredicted_scene <- predicted_scene
  tpredicted_scene[cell] <- as.matrix(tpredicted_ordination)
  names(tpredicted_scene) <- c("MDS1", "MDS2", "MDS3",
                              "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                              "amp-total")
  
  #-----------------------------------------------------------
  # Scene export
  
  toutput_raster <- paste0(output_path, "/", tscene_name)
  
  #Export
  writeRaster(x = tpredicted_scene,
              filename = toutput_raster,
              names = c("MDS1", "MDS2", "MDS3",
                        "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                        "amp-total"),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
  ### PHYLOGENETIC ###
  # As matrices
  pMDS1 <- as.matrix(coefficients_phylogenetic[MDS == "1", .SD, .SDcols = 3:ncol(coefficients_phylogenetic)])
  pMDS2 <- as.matrix(coefficients_phylogenetic[MDS == "2", .SD, .SDcols = 3:ncol(coefficients_phylogenetic)])
  pMDS3 <- as.matrix(coefficients_phylogenetic[MDS == "3", .SD, .SDcols = 3:ncol(coefficients_phylogenetic)])
  
  # Predict
  ppredicted_ordination <- predict_pls(as.matrix(sw),
                                       pMDS1,
                                       pMDS2,
                                       pMDS3)
  
  ppredicted_ordination <- as.data.table(ppredicted_ordination)
  colnames(ppredicted_ordination) <- c("MDS1", "MDS2", "MDS3",
                                      "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                                      "amp-total")
  
  #--------------------------------------------------------
  #Get a scene frame
  ppredicted_scene <- predicted_scene
  ppredicted_scene[cell] <- as.matrix(ppredicted_ordination)
  names(ppredicted_scene) <- c("MDS1", "MDS2", "MDS3",
                               "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                               "amp-total")
  
  #-----------------------------------------------------------
  # Scene export
  
  poutput_raster <- paste0(output_path, "/", pscene_name)
  
  #Export
  writeRaster(x = ppredicted_scene,
              filename = poutput_raster,
              names = c("MDS1", "MDS2", "MDS3",
                        "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                        "amp-total"),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
  ### FUNCTIONAL ###
  # As matrices
  fMDS1 <- as.matrix(coefficients_functional[MDS == "1", .SD, .SDcols = 3:ncol(coefficients_functional)])
  fMDS2 <- as.matrix(coefficients_functional[MDS == "2", .SD, .SDcols = 3:ncol(coefficients_functional)])
  fMDS3 <- as.matrix(coefficients_functional[MDS == "3", .SD, .SDcols = 3:ncol(coefficients_functional)])
  
  # Predict
  fpredicted_ordination <- predict_pls(as.matrix(sw),
                                       fMDS1,
                                       fMDS2,
                                       fMDS3)
  
  fpredicted_ordination <- as.data.table(fpredicted_ordination)
  colnames(fpredicted_ordination) <- c("MDS1", "MDS2", "MDS3",
                                       "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                                       "amp-total")
  
  #--------------------------------------------------------
  #Get a scene frame
  fpredicted_scene <- predicted_scene
  fpredicted_scene[cell] <- as.matrix(fpredicted_ordination)
  names(fpredicted_scene) <- c("MDS1", "MDS2", "MDS3",
                               "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                               "amp-total")
  
  #-----------------------------------------------------------
  # Scene export
  
  foutput_raster <- paste0(output_path, "/", fscene_name)
  
  #Export
  writeRaster(x = fpredicted_scene,
              filename = foutput_raster,
              names = c("MDS1", "MDS2", "MDS3",
                        "amp-MDS1", "amp-MDS2", "amp-MDS3", 
                        "amp-total"),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
}


#'------------------------------------------------------------------------------
#' @example 

# Apply in a batch of scenes
root_path <- "C:/Users/guzman/Documents/test"
scenes <- list.files(path = root_path,
                     recursive = TRUE,
                     full.names = TRUE)
temp <- "C:/Users/guzman/Documents/temp"
coefficients <- fread(paste0(root_path, "/06_model/SW_beta_PLS/DESIS/DESIS_", "coefficients.csv"))
output_path <- paste0(root_path, "/00_imagery/predicted/DESIS")

for(i in 1:length(scenes)) {
  
  apply_model(scenes[i], 
              scene_main_path, 
              coefficients_taxonomic,
              coefficients_phylogenetic,
              coefficients_functional,
              output_path,
              temp,
              threads = 64)
  
}
