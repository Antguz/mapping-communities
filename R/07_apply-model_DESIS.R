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
library(pbmcapply)

#'------------------------------------------------------------------------------
#' Source code

#source("R/DESIS-transform.R")
#sourceCpp("src/beta_predict_pls.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param scene_main_path
#' @param coefficients
#' @param output_raster
#' @param threads

#'------------------------------------------------------------------------------
#' Function

apply_model <- function(X,
                        scene_main_path, 
                        coefficients_taxonomic,
                        coefficients_phylogenetic,
                        coefficients_functional,
                        output_path,
                        temp,
                        threads = 28) {
  
  #---------------------------------------------------------
  # Unzip files
  name_product <- strsplit(scene_main_path[X], "/")[[1]][length(strsplit(scene_main_path[X], "/")[[1]])]
  name_product <- strsplit(name_product, ".zip")[[1]]
  scene_path <- paste0(temp, "/", name_product)
  unzip(scene_main_path[X], exdir = paste0(temp, "/", name_product), overwrite = TRUE)
  
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
  coef_bands <- colnames(coefficients_taxonomic)[4:ncol(coefficients_taxonomic)]
  sw <- sw[, ..coef_bands]
  
  # Release memory
  rm(list = c('data_spectra'))
  
  #--------------------------------------------------------
  # Predict raster files
  
  # Final names
  scene_name <- strsplit(name_product, "DESIS-HSI-L2A-")[[1]][2]
  tscene_name <- paste0("DESIS-L3A_", scene_name, "_NDMS-taxonomic-composition.tif")
  pscene_name <- paste0("DESIS-L3A_", scene_name, "_NDMS-phylogenetic-composition.tif")
  fscene_name <- paste0("DESIS-L3A_", scene_name, "_NDMS-functional-composition.tif")
  
  ######## TAXONOMIC ########
  # As matrices
  tAxis1 <- as.matrix(coefficients_taxonomic[Axes == "Axis1", .SD, .SDcols = 3:ncol(coefficients_taxonomic)])
  tAxis2 <- as.matrix(coefficients_taxonomic[Axes == "Axis2", .SD, .SDcols = 3:ncol(coefficients_taxonomic)])
  tAxis3 <- as.matrix(coefficients_taxonomic[Axes == "Axis3", .SD, .SDcols = 3:ncol(coefficients_taxonomic)])
  
  # Predict
  tpredicted_ordination <- beta_predict_pls(as.matrix(sw),
                                            tAxis1,
                                            tAxis2,
                                            tAxis3)
  
  #Get a scene frame
  tpredicted_scene <- predicted_scene
  tpredicted_scene[cell] <- tpredicted_ordination
  names(tpredicted_scene) <- c("Axis1", "Axis2", "Axis3",
                               "amp-Axis1", "amp-Axis2", "amp-Axis3", 
                               "amp-total")
  
  # Scene export
  toutput_raster <- paste0(output_path, "/", tscene_name)
  writeRaster(x = tpredicted_scene,
              filename = toutput_raster,
              names = c("Axis1", "Axis2", "Axis3",
                        "amp-Axis1", "amp-Axis2", "amp-Axis3", 
                        "amp-total"),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
  # Release memory
  rm(list = c('scene_name', 'tscene_name', 'tAxis1', 'tAxis2', 'tAxis3',
              'tpredicted_ordination', 'tpredicted_scene', 'toutput_raster'))
  
  ######## PHYLOGENETIC ########
  # As matrices
  pAxis1 <- as.matrix(coefficients_phylogenetic[Axes == "Axis1", .SD, .SDcols = 3:ncol(coefficients_phylogenetic)])
  pAxis2 <- as.matrix(coefficients_phylogenetic[Axes == "Axis2", .SD, .SDcols = 3:ncol(coefficients_phylogenetic)])
  pAxis3 <- as.matrix(coefficients_phylogenetic[Axes == "Axis3", .SD, .SDcols = 3:ncol(coefficients_phylogenetic)])
  
  # Predict
  ppredicted_ordination <- beta_predict_pls(as.matrix(sw),
                                            pAxis1,
                                            pAxis2,
                                            pAxis3)
  
  #Get a scene frame
  ppredicted_scene <- predicted_scene
  ppredicted_scene[cell] <- ppredicted_ordination
  names(ppredicted_scene) <- c("Axis1", "Axis2", "Axis3",
                               "amp-Axis1", "amp-Axis2", "amp-Axis3", 
                               "amp-total")
  
  # Scene export
  poutput_raster <- paste0(output_path, "/", pscene_name)
  writeRaster(x = ppredicted_scene,
              filename = poutput_raster,
              names = c("Axis1", "Axis2", "Axis3",
                        "amp-Axis1", "amp-Axis2", "amp-Axis3", 
                        "amp-total"),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
  # Release memory
  rm(list = c('pscene_name', 'pAxis1', 'pAxis2', 'pAxis3',
              'ppredicted_ordination', 'ppredicted_scene', 'poutput_raster'))
  
  ######## FUNCTIONAL ########
  # As matrices
  fAxis1 <- as.matrix(coefficients_functional[Axes == "Axis1", .SD, .SDcols = 3:ncol(coefficients_functional)])
  fAxis2 <- as.matrix(coefficients_functional[Axes == "Axis2", .SD, .SDcols = 3:ncol(coefficients_functional)])
  fAxis3 <- as.matrix(coefficients_functional[Axes == "Axis3", .SD, .SDcols = 3:ncol(coefficients_functional)])
  
  # Predict
  fpredicted_ordination <- beta_predict_pls(as.matrix(sw),
                                            fAxis1,
                                            fAxis2,
                                            fAxis3)
  
  # Get a scene frame
  fpredicted_scene <- predicted_scene
  fpredicted_scene[cell] <- fpredicted_ordination
  names(fpredicted_scene) <- c("Axis1", "Axis2", "Axis3",
                               "amp-Axis1", "amp-Axis2", "amp-Axis3", 
                               "amp-total")
  
  # Scene export
  foutput_raster <- paste0(output_path, "/", fscene_name)
  writeRaster(x = fpredicted_scene,
              filename = foutput_raster,
              names = c("Axis1", "Axis2", "Axis3",
                        "amp-Axis1", "amp-Axis2", "amp-Axis3", 
                        "amp-total"),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
  # Release memory
  rm(list = c('fscene_name', 'fAxis1', 'fAxis2', 'fAxis3',
              'fpredicted_ordination', 'fpredicted_scene', 'foutput_raster'))
  
  #--------------------------------------------------------
  # Remove temporal folder
  unlink(paste0(temp, "/", name_product), recursive = TRUE)
  
  # Return
  gc()
  return(paste0(name_product, ": success"))
  
}

#'------------------------------------------------------------------------------
#' @example 

################################################################################
#' @Local predictions from .zip files
################################################################################

# # Apply in a batch of scenes
# root_path <- "/media/antonio/Work/DESIS_mapping"
# scenes <- list.files(path = paste0(root_path, "/00_imagery/DESIS/test"),
#                      full.names = TRUE,
#                      recursive = FALSE)
# type <- "taxo"
# coefficients_taxonomic <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/DESIS/DESIS_", type, "_coefficients.csv"))
# type <- "phylo-log"
# coefficients_phylogenetic <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/DESIS/DESIS_", type, "_coefficients.csv"))
# type <- "func"
# coefficients_functional <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/DESIS/DESIS_", type, "_coefficients.csv"))
# 
# temp <- paste0(root_path, "/00_imagery/predicted/DESIS/temp")
# output_path <- paste0(root_path, "/00_imagery/Predicted/EMIT/L3")
# 
# # Parallel processing
# run <- pbmclapply(X = 1:length(scenes), 
#                   FUN = apply_model, 
#                   scene_main_path = scenes, 
#                   coefficients_taxonomic = coefficients_taxonomic,
#                   coefficients_phylogenetic = coefficients_phylogenetic,
#                   coefficients_functional = coefficients_functional,
#                   output_path = output_path,
#                   temp= temp,
#                   threads = 1,
#                   mc.preschedule = TRUE, 
#                   mc.set.seed = FALSE,
#                   mc.silent = FALSE, 
#                   mc.cores = 12,
#                   mc.cleanup = TRUE, 
#                   mc.allow.recursive = TRUE, 
#                   affinity.list = NULL)
# 
# complete <- do.call(rbind, run)
# fwrite(complete, paste0(output_path, "/files ", Sys.time(), ".txt"))

################################################################################
#' @MSI predictions from .zip files
################################################################################

# Apply in a batch of scenes
# Path
root_path <- "/panfs/jay/groups/17/cavender/guzman/Imagery"

# Source code
source(paste0(root_path, "/Source/DESIS-transform.R"))
sourceCpp(paste0(root_path, "/Source/beta_predict_pls.cpp"))

# Scenes path
scenes <- list.files(path = paste0(root_path, "/DESIS"),
                     full.names = TRUE,
                     recursive = TRUE)

# Coefficients
type <- "taxo"
coefficients_taxonomic <- fread(paste0(root_path, 
                                       "/Model/DESIS_", 
                                       type, 
                                       "_coefficients.csv"))
type <- "phylo-log"
coefficients_phylogenetic <- fread(paste0(root_path, 
                                          "/Model/DESIS_", 
                                          type, 
                                          "_coefficients.csv"))
type <- "func"
coefficients_functional <- fread(paste0(root_path, 
                                        "/Model/DESIS_", 
                                        type, 
                                        "_coefficients.csv"))

# Out paths
temp <- paste0(root_path, "/temp")
output_path <- paste0(root_path, "/Predicted/DESIS/L3")

# Parallel processing
run <- pbmclapply(X = 1:length(scenes), 
                  FUN = apply_model, 
                  scene_main_path = scenes, 
                  coefficients_taxonomic = coefficients_taxonomic,
                  coefficients_phylogenetic = coefficients_phylogenetic,
                  coefficients_functional = coefficients_functional,
                  output_path = output_path,
                  temp= temp,
                  threads = 1,
                  mc.preschedule = TRUE, 
                  mc.set.seed = FALSE,
                  mc.silent = FALSE, 
                  mc.cores = 32,
                  mc.cleanup = TRUE, 
                  mc.allow.recursive = TRUE, 
                  affinity.list = NULL)

# Export report
complete <- do.call(rbind, run)
fwrite(complete, paste0(output_path, "/files_DESIS_", Sys.time(), ".txt"))
