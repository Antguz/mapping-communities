#'------------------------------------------------------------------------------
#' @title Compute species approximation of the predicted values
#'------------------------------------------------------------------------------

#' @description It computes the approximation of the predicted values to species.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(terra)
library(Rcpp)
library(RcppArmadillo)

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/approximation_sp.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param predicted_scene
#' @param wascores
#' @param output_raster
#' @param threads

#'------------------------------------------------------------------------------
#' Function

apply_model <- function(predicted_scene, 
                        wascores, 
                        output_raster,
                        threads = 28) {
  
  #Ask data frame
  frame_scene <- as.data.table(as.data.frame(predicted_scene, cells = TRUE, na.rm = TRUE))
  ordination <- as.matrix(frame_scene[, 2:4]/1000)
  
  #---------------------------------------------------------
  # Estimate distance
  distance <- approximation_sp(ordination, as.matrix(wascores[, 2:4]))
  
  #Get the inverse
  distance <- 1/distance
  total <- rowSums(distance, na.rm=FALSE)
  
  #Get fraction
  fraction <- distance/total
  fraction <- round(fraction*10000, 0)
  
  max <- NA
  
  for(i in 1:nrow(fraction)) {
    
    max[i] <- which.max(fraction[i,])
    
  }
  
  fraction <- cbind(fraction, max)
  
  #--------------------------------------------------------
  #Get a scene frame
  species_scene <- subset(predicted_scene, 1)
  species_scene <- rep(species_scene, (nrow(wascores)+1))
  species_scene[frame_scene$cell] <- fraction
  names(species_scene) <- c(wascores$Taxon2, "sp-max")
  
  #-----------------------------------------------------------
  # Scene export
  
  #Export
  writeRaster(x = species_scene,
              filename = output_raster,
              names = c(wascores$Taxon2, "sp-max"),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
}

#'------------------------------------------------------------------------------
#' @example

# Run in a single scene

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"
predicted_scene <- rast(paste0(root_path, "/00_imagery/", "DESIS-HSI-L2A-DT0621221713_022-20220621T171230-V0217_predicted.tif"))
wascores <- readRDS(paste0(root_path, "/05_ordination/species_scores_BA.rds"))
output_raster <- paste0(root_path, "/00_imagery/", "DESIS-HSI-L2A-DT0621221713_022-20220621T171230-V0217_species.tif")
threads <- 28

apply_model(predicted_scene, 
            wascores, 
            output_raster,
            threads = 28)

# Run in a batch

folder_path <- paste0(root_path, "/00_imagery/predicted")
output_path <- paste0(root_path, "/00_imagery/species")

files <- list.files(path = folder_path, 
                    pattern = ".tif", 
                    recursive = TRUE,
                    include.dirs = TRUE)

for(i in 1:length(files)) {
  
  species_name <- substr(files[i], 0, 52)
  scene <- rast(paste0(folder_path, "/", files[i]))
  
  apply_model(predicted_scene = scene, 
              wascores = wascores, 
              output_raster = paste0(output_path, "/", species_name, "_species.tif"),
              threads = 28)
  
}
  
  



