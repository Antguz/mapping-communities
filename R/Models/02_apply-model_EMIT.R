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
library(pbmcapply)

#'------------------------------------------------------------------------------
#' Source code

#source("R/EMIT-transform.R")
#sourceCpp("src/model-application.cpp")

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
                        frame_scenes, 
                        coefficients_taxonomic,
                        coefficients_phylogenetic,
                        coefficients_functional,
                        output_path,
                        threads = 28) {
  
  #---------------------------------------------------------
  # Select and read products
  scene <- frame_scenes[ID == X,]
  
  spectra <- scene[type == "RFL"]$scenes_path
  mask <- scene[type == "MASK" & band_mask == "mask"]$scenes_path
  
  # Read products
  rast_spectra <- rast(spectra)
  names(rast_spectra) <- EMIT_bands
  rast_mask <- rast(mask)
  names(rast_mask)<- c("cloud", "cirrus", "water", "spacecraft", "dilated",
                       "AOD550", "H20", "mask")
  
  #---------------------------------------------------------
  # Define mask
  
  # Mask for NIR
  NIR <- rast_spectra$`805.361755371093`
  NIR[NIR < 0.1] <- 1
  NIR[NIR >= 0.1 & NIR <= 0.8] <- 0
  NIR[NIR > 0.8] <- 1
  
  # Mask for NDVI
  NDVI <- (rast_spectra$`805.361755371093` - rast_spectra$`686.010314941406`) / 
          (rast_spectra$`805.361755371093` + rast_spectra$`686.010314941406`)
  NDVI[NDVI == Inf] <- 0
  NDVI[NDVI >= 0.40] <- 2
  NDVI[NDVI < 0.40] <- 0
  NDVI[NDVI == 0] <- 1
  NDVI[NDVI == 2] <- 0
  
  # Quality mask
  rast_mask <- rast_mask$mask
  rast_mask[rast_mask != 1] <- 0
  rast_mask[rast_mask == 1] <- 1
  
  # Final mask 
  rast_mask <- sum(rast_mask, NDVI, NIR, na.rm = TRUE)
  rast_mask[rast_mask >= 1] <- 1
  rast_mask[rast_mask < 1] <- 0
  data_mask <- as.data.table(as.data.frame(rast_mask, cells = TRUE, na.rm = TRUE))
  data_mask <- data_mask[mask == 0,]
  
  # As frame
  data_spectra <- as.data.table(as.data.frame(rast_spectra, cells = TRUE, na.rm = TRUE))
  colnames(data_spectra)[2:ncol(data_spectra)] <- as.character(EMIT_bands)
  
  # Mask spectra
  cells_match <- data_spectra$cell %in% data_mask$cell
  data_spectra <- data_spectra[cells_match == TRUE,]
  
  # To predict
  predicted_scene <- subset(rast_spectra, 1:7)
  predicted_scene[] <- NA
  
  # Release memory
  rm(list = c('rast_spectra', 'rast_mask', 'data_mask','NDVI', 'NIR', 'cells_match'))
  gc()
  
  #---------------------------------------------------------
  # Transform raster files
  
  cell <- data_spectra$cell
  
  sw <- EMIT_transform(spectra = data_spectra[, .SD, .SDcols = 2:ncol(data_spectra)], 
                       EMIT_bands = EMIT_bands,
                       cell = cell,
                       new_bands = seq(min(EMIT_bands), max(EMIT_bands), by = 7.5),
                       threads = 1)
  
  # Subset bands
  cell <- sw$cell
  coef_bands <- colnames(coefficients_taxonomic)[4:ncol(coefficients_taxonomic)]
  sw <- sw[, ..coef_bands]
  
  # Release memory
  rm(list = c('data_spectra'))
  
  #--------------------------------------------------------
  # Predict raster files
  
  # Final names
  scene_name <- strsplit(scene[type == "RFL"]$scenes, "EMIT_L2A_RFL_")[[1]][2]
  scene_name <- strsplit(scene_name, "_reflectance")[[1]][1]
  tscene_name <- paste0("L3A/EMIT-L3A_", scene_name, "_NDMS-taxonomic-composition.tif")
  pscene_name <- paste0("L3B/EMIT-L3B_", scene_name, "_NDMS-phylogenetic-composition.tif")
  fscene_name <- paste0("L3C/EMIT-L3C_", scene_name, "_NDMS-functional-composition.tif")
  
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
  
  # Return
  gc()
  return(paste0(scene[type == "RFL"]$scenes, ": success"))
  
  #tmpFiles(current=TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
  
}

#'------------------------------------------------------------------------------
#' @example 

################################################################################
#' @Local predictions from .hdr files
################################################################################

# Apply in a batch of scenes
# Path
root_path <- "/media/antonio/Work/MMDB"

# Source code
source("R/auxiliary/EMIT-transform.R")
sourceCpp("src/beta_predict_pls.cpp")

# Scenes frame
frame_scenes <- list.files(path = paste0(root_path, "/00_imagery/EMIT"),
                           full.names = TRUE,
                           recursive = TRUE)

# Arrange path in a frame
frame_scenes <- data.table(scenes_path = frame_scenes)

# Get product information
n <- length(strsplit(frame_scenes$scenes_path[1], "/")[[1]])
frame_scenes[, scenes := strsplit(scenes_path, "/")[[1]][n], by = seq_along(1:nrow(frame_scenes))]
frame_scenes[, type := strsplit(scenes, "_")[[1]][3], by = seq_along(1:nrow(frame_scenes))]
frame_scenes[, time := strsplit(scenes, "_")[[1]][5], by = seq_along(1:nrow(frame_scenes))]
frame_scenes[, band_mask := strsplit(scenes, "_")[[1]][8], by = seq_along(1:nrow(frame_scenes))]
frame_scenes[, keep := unlist(gregexpr('.hdr', scenes))[1], by = seq_along(1:nrow(frame_scenes))]
frame_scenes <- frame_scenes[keep == -1, 1:5]

# Unique scenes and merge
list_scenes <- data.table(time = unique(frame_scenes$time))
list_scenes$ID <- 1:nrow(list_scenes)
frame_scenes <- merge(frame_scenes,
                      list_scenes,
                      by = "time",
                      all.x = TRUE,
                      all.y = TRUE)

# Coefficients
type <- "taxo"
coefficients_taxonomic <- fread(paste0(root_path,
                                       "/Models/01_model_beta/EMIT/EMIT_",
                                       type,
                                       "_coefficients.csv"))
type <- "phylo"
coefficients_phylogenetic <- fread(paste0(root_path,
                                          "/Models/01_model_beta/EMIT/EMIT_",
                                          type,
                                          "_coefficients.csv"))
type <- "func"
coefficients_functional <- fread(paste0(root_path,
                                        "/Models/01_model_beta/EMIT/EMIT_",
                                        type,
                                        "_coefficients.csv"))
# Out paths
output_path <- paste0(root_path, "/00_imagery/Predicted/EMIT/L3")

# Parallel processing
run <- pbmclapply(X = 1:length(unique(frame_scenes$ID)),
                  FUN = apply_model,
                  frame_scenes = frame_scenes,
                  coefficients_taxonomic = coefficients_taxonomic,
                  coefficients_phylogenetic = coefficients_phylogenetic,
                  coefficients_functional = coefficients_functional,
                  output_path = output_path,
                  threads = 1,
                  mc.preschedule = TRUE,
                  mc.set.seed = FALSE,
                  mc.silent = FALSE,
                  mc.cores = 2,
                  mc.cleanup = TRUE,
                  mc.allow.recursive = TRUE,
                  affinity.list = NULL)

# Export report
complete <- do.call(rbind, run)
fwrite(complete, paste0(output_path, "/files ", Sys.time(), ".txt"))

################################################################################
#' @MSI predictions from .hdr files
################################################################################

# Apply in a batch of scenes
# Path
root_path <- "/panfs/jay/groups/17/cavender/guzman/Imagery"

# Source code
source(paste0(root_path, "/Source/EMIT-transform.R"))
sourceCpp(paste0(root_path, "/Source/beta_predict_pls.cpp"))

# Scenes frame
frame_scenes <- list.files(path = paste0(root_path, "/EMIT/2023"),
                           full.names = TRUE,
                           recursive = TRUE)

# Arrange path in a frame
frame_scenes <- data.table(scenes_path = frame_scenes)

# Get product information
n <- length(strsplit(frame_scenes$scenes_path[1], "/")[[1]])
frame_scenes[, scenes := strsplit(scenes_path, "/")[[1]][n], by = seq_along(1:nrow(frame_scenes))]
frame_scenes[, type := strsplit(scenes, "_")[[1]][3], by = seq_along(1:nrow(frame_scenes))]
frame_scenes[, time := strsplit(scenes, "_")[[1]][5], by = seq_along(1:nrow(frame_scenes))]
frame_scenes[, band_mask := strsplit(scenes, "_")[[1]][8], by = seq_along(1:nrow(frame_scenes))]
frame_scenes[, keep := unlist(gregexpr('.hdr', scenes))[1], by = seq_along(1:nrow(frame_scenes))]
frame_scenes <- frame_scenes[keep == -1, 1:5]

# Unique scenes and merge
list_scenes <- data.table(time = unique(frame_scenes$time))
list_scenes$ID <- 1:nrow(list_scenes)
frame_scenes <- merge(frame_scenes, 
                      list_scenes, 
                      by = "time", 
                      all.x = TRUE, 
                      all.y = TRUE)

# Coefficients
type <- "taxo"
coefficients_taxonomic <- fread(paste0(root_path, 
                                       "/Model/EMIT_", 
                                       type, 
                                       "_coefficients.csv"))
type <- "phylo-log"
coefficients_phylogenetic <- fread(paste0(root_path, 
                                          "/Model/EMIT_", 
                                          type, 
                                          "_coefficients.csv"))
type <- "func"
coefficients_functional <- fread(paste0(root_path, 
                                        "/Model/EMIT_", 
                                        type, 
                                        "_coefficients.csv"))
# Out paths
output_path <- paste0(root_path, "/Predicted/EMIT/L3")

# Parallel processing
run <- pbmclapply(X = 1:length(unique(frame_scenes$ID)), 
                  FUN = apply_model, 
                  frame_scenes = frame_scenes, 
                  coefficients_taxonomic = coefficients_taxonomic,
                  coefficients_phylogenetic = coefficients_phylogenetic,
                  coefficients_functional = coefficients_functional,
                  output_path = output_path,
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
fwrite(complete, paste0(output_path, "/files_EMIT_", Sys.time(), ".txt"))
