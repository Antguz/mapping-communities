#'------------------------------------------------------------------------------
#' @title Apply occurrence of plant groups to scenes
#'------------------------------------------------------------------------------

#' @description It applies the occurrence models (i.e., probability) of plant
#' lineages (03_model-lineages.R) to the predicted MDS axes of phylogenetic
#' composition within scenes.

#' @return GeoTIFF files per scene (L4B) with the mean occurrence probability
#' per lineage (_plant-occurrence_mean.tif) and the width of its 95% confidence
#' interval among models (_plant-occurrence_error.tif), and a .txt file with the
#' list of processed scenes.

#' @Note This script was not implemented in the manuscript.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(terra)
library(Rcpp)
library(RcppArmadillo)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' Source code

# sourceCpp("src/probability_rcpp.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param predicted_scene_path
#' @param coefficients
#' @param output_raster
#' @param threads

#'------------------------------------------------------------------------------
#' Function

apply_model_occurence <- function(X,
                                  sensor,
                                  phylogenetic_composition, 
                                  coefficients, 
                                  climatic,
                                  topographic,
                                  scaling,
                                  output_path) {
  
  # ----------------------------------------------------------------------------
  # Read predicted scene
  predicted_scene <- rast(phylogenetic_composition[X])
  
  #Ask data frame
  frame_scene <- as.data.table(as.data.frame(predicted_scene, 
                                             cells = TRUE, 
                                             na.rm = TRUE,
                                             xy = TRUE))
  ordination <- frame_scene[, 4:6]/10000
  
  # Create a vector for extraction of climatic and topographic data
  vector <- vect(x = frame_scene[, 1:3], 
                 geom=c("x", "y"), 
                 crs = crs(predicted_scene), 
                 keepgeom = FALSE)
  frame_scene <- frame_scene[,1]
  
  # ----------------------------------------------------------------------------
  # Get climatic data
  
  # Test projection
  if(crs(vector) != crs(climatic)) {
    points <- project(vector, climatic)
    climatic_extract <- extract(climatic, points, cells=FALSE, xy=FALSE)
    
  } else {
    climatic_extract <- extract(climatic, vector, cells=FALSE, xy=FALSE)
  }
  
  ordination <- cbind(ordination, climatic_extract[,-1])
  
  # ----------------------------------------------------------------------------
  # Get topographic data
  
  # Test projection
  if(crs(vector) != crs(topographic)) {
    points <- project(vector, topographic)
    topographic_extract <- extract(topographic, points, cells=FALSE, xy=FALSE)
    
  } else {
    topographic_extract <- extract(topographic, vector, cells=FALSE, xy=FALSE)
  }
  
  ordination <- cbind(ordination, topographic_extract[,-1])
  
  # Remove residuals
  rm(list = c('vector', 'points', 'climatic_extract', 'topographic_extract'))
  
  # ----------------------------------------------------------------------------
  # Check for potential NA
  
  # Replace potential NaN with NA
  ordination[, (names(ordination)) := lapply(.SD, function(x) replace(x, is.nan(x), NA))]
  
  # Check for potential NA
  NA_row <- apply(ordination, 1, anyNA)
  
  # Replace NA for zeros
  ordination[, (names(ordination)) := lapply(.SD, function(x) fifelse(is.na(x), 0, x))]
  
  # ----------------------------------------------------------------------------
  # Scaling of values
  
  # Scale values
  ordination <- (ordination - as.list(as.numeric(scaling[1]))) / as.list(as.numeric(scaling[2]))
  
  # ----------------------------------------------------------------------------
  # Reshape coefficients
  models <- unique(coefficients$plant)
  coef_array <- array(dim = c(100, 11, length(models)))
  intercept_mat <- matrix(nrow = 100, ncol = length(models))

  for(i in 1:length(models)) {

    sub <- coefficients[plant == models[i]]
    coef_array[,,i] <- as.matrix(sub[, 3:13])
    intercept_mat[,i] <- sub$`(Intercept)`

  }

  #---------------------------------------------------------
  # Estimate probabilities based on coefficients
  plant_probabilities <- summary_prob(as.matrix(ordination),
                                      coef_array,
                                      intercept_mat)
  rm(list = c('ordination'))
  
  #---------------------------------------------------------
  # Replace NA
  
  if(length(NA_row[NA_row == TRUE]) != 0) {
    
    plant_probabilities[[1]][NA_row == TRUE,] <- NA
    plant_probabilities[[2]][NA_row == TRUE,] <- NA
    
  }
  
  #--------------------------------------------------------
  #Get a scene frame
  species_scene_mean <- subset(predicted_scene, 1)
  species_scene_mean[] <- NA
  species_scene_mean <- rep(species_scene_mean, length(models))
  species_scene_error <- species_scene_mean
  species_scene_mean[frame_scene$cell] <- plant_probabilities[[1]]
  species_scene_error[frame_scene$cell] <- plant_probabilities[[2]]
  names(species_scene_mean) <- models
  names(species_scene_error) <- models
  rm(list = c('plant_probabilities'))
  
  #-----------------------------------------------------------
  # Scene export
  
  scene_name <- strsplit(phylogenetic_composition[X], "/")[[1]]
  scene_name <- scene_name[length(scene_name)]
  
  if(sensor == "DESIS") {
    scene_name <- strsplit(scene_name, "DESIS-L3B")[[1]][2]
    scene_name <- paste0("DESIS-L4B", scene_name)
    
  } else if(sensor == "EMIT") {
    scene_name <- strsplit(scene_name, "EMIT-L3B")[[1]][2]
    scene_name <- paste0("EMIT-L4B", scene_name)
  }
  
  scene_name <- strsplit(scene_name, "MDS-phylogenetic-composition.tif")[[1]][1]
  scene_name_mean <- paste0(scene_name, "plant-occurrence_mean.tif")
  scene_name_error <- paste0(scene_name, "plant-occurrence_error.tif")
  
  output_raster_mean <- paste0(output_path, "/", scene_name_mean)
  output_raster_error <- paste0(output_path, "/", scene_name_error)

  #Export
  writeRaster(x = species_scene_mean,
              filename = output_raster_mean,
              names = models,
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
  writeRaster(x = species_scene_error,
              filename = output_raster_error,
              names = models,
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
}

#'------------------------------------------------------------------------------
#' @example

################################################################################
#' @Local Predictions from phylogenetic-composition.tif files
################################################################################

# sourceCpp("src/probability_rcpp.cpp")
# 
# # Define root path, sensor, and coefficients
# root_path <- "/media/antonio/Work/MMDB"
# sensor <- "EMIT"
# type <- "phylo"
# 
# # Read coefficients
# coefficients <- fread(paste0(root_path, "/Models/03_model_occurrence/", 
#                              type, "_model-occurrence.csv"))
# coefficients <- coefficients[dataset == "training", 
#                              .SD,
#                              .SDcols = c(1, 15:26)]
# 
# # Read climatic and topographic data
# MAT <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_1.tif"))
# TS <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_4.tif"))
# TAR <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_7.tif"))
# MAP <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_12.tif"))
# PS <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_15.tif"))
# climatic <- c(MAT, TS, TAR, MAP, PS)
# names(climatic) <- c("MAT", "TS", "TAR", "MAP", "PS")
# rm(list = c("MAT", "TS", "TAR", "MAP", "PS"))
# 
# # Read DEM information
# DEM <- rast(paste0(root_path, "/env/gt30h1kna/na_dem.bil"))
# DEM[DEM == -9999] <- NA
# names(DEM) <- "DEM"
# 
# # Read HYDRO1k information
# slope <- rast(paste0(root_path, "/env/gt30h1kna/na_slope.bil"))
# slope[slope == -9999] <- NA
# CTI <- rast(paste0(root_path, "/env/gt30h1kna/na_cti.bil"))
# CTI[CTI == -9999] <- NA
# HYDRO1k <- c(slope, CTI)
# names(HYDRO1k) <- c("SL", "CTI")
# crs(HYDRO1k) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"
# 
# # Merge
# topographic <- c(DEM, HYDRO1k)
# rm(list = c("DEM", "CTI", "HYDRO1k"))
# 
# # Scaling coefficients
# scaling <- fread(paste0(root_path, "/Models/03_model_occurrence/", 
#                         type, "_scaling-occurrence.csv"))
# 
# # Define inputs and outputs folder
# folder_path <- paste0(root_path, "/00_imagery/Predicted/", sensor, "/L3/L3B")
# output_path <- paste0(root_path, "/00_imagery/Predicted/", sensor, "/L4/Occurrence")
# 
# # Select files to process
# phylogenetic_composition <- list.files(path = folder_path, 
#                                        pattern = "phylogenetic-composition.tif", 
#                                        full.names = TRUE,
#                                        recursive = TRUE)
# 
# # Parallel processing
# run <- pbmclapply(X = 1:length(phylogenetic_composition), 
#                   FUN = apply_model_occurence, 
#                   sensor = sensor,
#                   phylogenetic_composition = phylogenetic_composition, 
#                   coefficients = coefficients, 
#                   climatic = climatic,
#                   topographic = topographic,
#                   scaling = scaling,
#                   output_path = output_path,
#                   mc.preschedule = TRUE, 
#                   mc.set.seed = FALSE,
#                   mc.silent = FALSE, 
#                   mc.cores = 1,
#                   mc.cleanup = TRUE, 
#                   mc.allow.recursive = TRUE, 
#                   affinity.list = NULL)
# 
# # Export report
# complete <- do.call(rbind, run)
# fwrite(complete, paste0(output_path, "/files_", sensor, "_Occurrence_", Sys.time(), ".txt"))

################################################################################
#' @MSI predictions
################################################################################

# Define root path, sensor, and coefficients
setwd("/panfs/jay/groups/17/cavender/guzman/MMDB")
root_path <- "/panfs/jay/groups/17/cavender/guzman"
type <- "phylo"

# Source code
sourceCpp("src/probability_rcpp.cpp")

# Read coefficients
coefficients <- fread(paste0(root_path, "/Models/03_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- coefficients[dataset == "training", 
                             .SD,
                             .SDcols = c(1, 15:26)]

# Read climatic and topographic data
MAT <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_1.tif"))
TS <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_4.tif"))
TAR <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_7.tif"))
MAP <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_12.tif"))
PS <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_15.tif"))
climatic <- c(MAT, TS, TAR, MAP, PS)
names(climatic) <- c("MAT", "TS", "TAR", "MAP", "PS")
rm(list = c("MAT", "TS", "TAR", "MAP", "PS"))

# Read DEM information
DEM <- rast(paste0(root_path, "/env/gt30h1kna/na_dem.bil"))
DEM[DEM == -9999] <- NA
names(DEM) <- "DEM"

# Read HYDRO1k information
slope <- rast(paste0(root_path, "/env/gt30h1kna/na_slope.bil"))
slope[slope == -9999] <- NA
CTI <- rast(paste0(root_path, "/env/gt30h1kna/na_cti.bil"))
CTI[CTI == -9999] <- NA
HYDRO1k <- c(slope, CTI)
names(HYDRO1k) <- c("SL", "CTI")
crs(HYDRO1k) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"

# Merge
topographic <- c(DEM, HYDRO1k)
rm(list = c("DEM", "CTI", "HYDRO1k"))

# Scaling coefficients
scaling <- fread(paste0(root_path, "/Models/03_model_occurrence/", 
                        type, "_scaling-occurrence.csv"))

#-------------------------------------------------------------------------------
#' @EMIT

sensor <- "EMIT"

# Define inputs and outputs folder
folder_path <- paste0(root_path, "/Imagery/Predicted/", sensor, "/L3/L3B")
output_path <- paste0(root_path, "/Imagery/Predicted/", sensor, "/L4/Occurrence")

# Select files to process
phylogenetic_composition <- list.files(path = folder_path, 
                                       pattern = "phylogenetic-composition.tif", 
                                       full.names = TRUE,
                                       recursive = TRUE)

# Parallel processing
run <- pbmclapply(X = 1:length(phylogenetic_composition), 
                  FUN = apply_model_occurence, 
                  sensor = sensor,
                  phylogenetic_composition = phylogenetic_composition, 
                  coefficients = coefficients, 
                  climatic = climatic,
                  topographic = topographic,
                  scaling = scaling,
                  output_path = output_path,
                  mc.preschedule = TRUE, 
                  mc.set.seed = FALSE,
                  mc.silent = FALSE, 
                  mc.cores = 10,
                  mc.cleanup = TRUE, 
                  mc.allow.recursive = TRUE, 
                  affinity.list = NULL)

# Export report
complete <- do.call(rbind, run)
fwrite(complete, paste0(output_path, "/files_", sensor, "_Occurrence_", Sys.time(), ".txt"))

#-------------------------------------------------------------------------------
#' @DESIS

sensor <- "DESIS"

# Define inputs and outputs folder
folder_path <- paste0(root_path, "/Imagery/Predicted/", sensor, "/L3/L3B")
output_path <- paste0(root_path, "/Imagery/Predicted/", sensor, "/L4/Occurrence")

# Select files to process
phylogenetic_composition <- list.files(path = folder_path, 
                                       pattern = "phylogenetic-composition.tif", 
                                       full.names = TRUE,
                                       recursive = TRUE)

# Parallel processing
run <- pbmclapply(X = 1:length(phylogenetic_composition), 
                  FUN = apply_model_occurence, 
                  sensor = sensor,
                  phylogenetic_composition = phylogenetic_composition, 
                  coefficients = coefficients, 
                  climatic = climatic,
                  topographic = topographic,
                  scaling = scaling,
                  output_path = output_path,
                  mc.preschedule = TRUE, 
                  mc.set.seed = FALSE,
                  mc.silent = FALSE, 
                  mc.cores = 10,
                  mc.cleanup = TRUE, 
                  mc.allow.recursive = TRUE, 
                  affinity.list = NULL)

# Export report
complete <- do.call(rbind, run)
fwrite(complete, paste0(output_path, "/files_", sensor, "_Occurrence_", Sys.time(), ".txt"))

