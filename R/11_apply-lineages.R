#'------------------------------------------------------------------------------
#' @title Apply occurrence of plant groups to scenes
#'------------------------------------------------------------------------------

#' @description It apply the occurrence (i.e., probability) of plant groups to the
#' predicted axes within scenes

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(terra)
library(Rcpp)
library(RcppArmadillo)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/probability_rcpp.cpp")

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
                                  phylogenetic_composition, 
                                  coefficients, 
                                  output_path) {
  
  # Read predicted scene
  predicted_scene <- rast(phylogenetic_composition[X])
  
  #Ask data frame
  frame_scene <- as.data.table(as.data.frame(predicted_scene, cells = TRUE, na.rm = TRUE))
  ordination <- as.matrix(frame_scene[, 2:4]/1000)
  
  #---------------------------------------------------------
  # Estimate probabilities based on coefficients
  plant_probabilities <- multi_probability(ordination, 
                                           as.matrix(coefficients[, 2:5]),
                                           as.numeric(as.factor(coefficients$plant)))
  
  #--------------------------------------------------------
  #Get a scene frame
  species_scene_mean <- subset(predicted_scene, 1)
  species_scene_mean[] <- NA
  species_scene_mean <- rep(species_scene_mean, length(unique(coefficients$plant)))
  species_scene_error <- species_scene_mean
  species_scene_mean[frame_scene$cell] <- plant_probabilities[[1]]
  species_scene_error[frame_scene$cell] <- plant_probabilities[[2]]
  names(species_scene_mean) <- unique(coefficients$plant)
  names(species_scene_error) <- unique(coefficients$plant)
  
  #-----------------------------------------------------------
  # Scene export
  
  scene_name <- strsplit(phylogenetic_composition[X], "/")[[1]]
  scene_name <- scene_name[length(scene_name)]
  scene_name <- strsplit(scene_name, "DESIS-L3A")[[1]][2]
  scene_name <- paste0("DESIS-L4B", scene_name)
  scene_name <- strsplit(scene_name, "NDMS-phylogenetic-composition.tif")[[1]][1]
  scene_name_mean <- paste0(scene_name, "plant-groups_mean.tif")
  scene_name_error <- paste0(scene_name, "plant-groups_error.tif")
  output_raster_mean <- paste0(output_path, "/", scene_name_mean)
  output_raster_error <- paste0(output_path, "/", scene_name_error)
  
  #Export
  writeRaster(x = species_scene_mean,
              filename = output_raster_mean,
              names = unique(coefficients$plant),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
  writeRaster(x = species_scene_error,
              filename = output_raster_error,
              names = unique(coefficients$plant),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
}

#'------------------------------------------------------------------------------
#' @example

################################################################################
#' @Local predictions from phylogenetic-composition.tif files
################################################################################

# Define root path, sensor, and coefficients
root_path <- "/media/antonio/Work/MMDB"
sensor <- "DESIS"
type <- "phylo-log"
coefficients <- fread(paste0(root_path, "/07_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- coefficients[dataset == "training", c("plant", "(Intercept)", "Axis1", "Axis2", "Axis3")]

# Define inputs and outputs folder
folder_path <- paste0(root_path, "/00_imagery/predicted/", sensor)
output_path <- paste0(root_path, "/00_imagery/lineages/", sensor)

# Select files to process
phylogenetic_composition <- list.files(path = folder_path, 
                                       pattern = "phylogenetic-composition.tif", 
                                       full.names = TRUE,
                                       recursive = TRUE)

# Parallel processing
run <- pbmclapply(X = 1:length(phylogenetic_composition), 
                  FUN = apply_model_occurence, 
                  phylogenetic_composition = phylogenetic_composition, 
                  coefficients = coefficients, 
                  output_path = output_path,
                  mc.preschedule = TRUE, 
                  mc.set.seed = FALSE,
                  mc.silent = FALSE, 
                  mc.cores = 1,
                  mc.cleanup = TRUE, 
                  mc.allow.recursive = TRUE, 
                  affinity.list = NULL)

# Export report
complete <- do.call(rbind, run)
fwrite(complete, paste0(output_path, "/files_", type, " ", Sys.time(), ".txt"))