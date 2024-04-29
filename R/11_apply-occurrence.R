#'------------------------------------------------------------------------------
#' @title Apply occurance of plant groups to scenes
#'------------------------------------------------------------------------------

#' @description It apply the occurance (i.e., probability) of plant groups to the
#' predicted axes within scenes

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(terra)
library(Rcpp)
library(RcppArmadillo)

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

apply_model_occurence <- function(predicted_scene_path, 
                                  coefficients, 
                                  output_path) {
  
  # Read predicted scene
  predicted_scene <- rast(predicted_scene_path)
  
  #Ask data frame
  frame_scene <- as.data.table(as.data.frame(predicted_scene, cells = TRUE, na.rm = TRUE))
  ordination <- as.matrix(frame_scene[, 2:4]/1000)
  
  #---------------------------------------------------------
  # Estimate probabilities based on coefficients
  plant_probabilities <- multi_probability(ordination, 
                                           as.matrix(coefficients[, 2:5]))
  
  #--------------------------------------------------------
  #Get a scene frame
  species_scene <- subset(predicted_scene, 1)
  species_scene[] <- NA
  species_scene <- rep(species_scene, nrow(coefficients))
  species_scene[frame_scene$cell] <- plant_probabilities
  names(species_scene) <- coefficients$plant
  
  #-----------------------------------------------------------
  # Scene export
  
  scene_name <- strsplit(predicted_scene_path, "/")[[1]]
  scene_name <- scene_name[length(scene_name)]
  scene_name <- strsplit(scene_name, "_predicted.tif")[[1]][1]
  scene_name <- paste0(scene_name, "_plant-groups.tif")
  output_raster <- paste0(output_path, "/", scene_name)
  
  #Export
  writeRaster(x = species_scene,
              filename = output_raster,
              names = coefficients$plant,
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
}

#'------------------------------------------------------------------------------
#' @example

# Run in a batch
root_path <- "E:/Projects/DESIS_mapping"
sensor <- "DESIS"
type <- "func"
coefficients <- fread(paste0(root_path, "/07_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- coefficients[dataset == "training", c("plant", "(Intercept)", "Axis1", "Axis2", "Axis3")]
folder_path <- paste0(root_path, "/00_imagery/predicted/", sensor, "/", type)
output_path <- paste0(root_path, "/00_imagery/occurrence/", sensor, "/", type)

files <- list.files(path = folder_path, 
                    pattern = ".tif", 
                    recursive = TRUE,
                    include.dirs = TRUE)

pb <- txtProgressBar(min = 0, max = length(files), style = 3)

for(i in 1:length(files)) {
  
  setTxtProgressBar(pb, i)
  
  scene_name <- paste0(folder_path, "/", files[i])  
  
  apply_model_occurence(predicted_scene_path = scene_name, 
                        coefficients = coefficients,
                        output_path)
  
}
