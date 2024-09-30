#'------------------------------------------------------------------------------
#' @title Apply predictions of CWM to scenes
#'------------------------------------------------------------------------------

#' @description It apply CWM predictions of plant traits to the predicted axes 
#' within scenes

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(terra)
library(Rcpp)
library(RcppArmadillo)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/cwm_rcpp.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param predicted_scene_path
#' @param coefficients
#' @param output_raster

#'------------------------------------------------------------------------------
#' Function

apply_model_CWM <- function(X,
                            sensor,
                            functional_composition, 
                            coefficients,
                            scaling,
                            output_path) {
  
  # Read predicted scene
  predicted_scene <- rast(functional_composition[X])
  cwn_scene <- subset(predicted_scene, 1)
  cwn_scene[] <- NA
  
  #Ask data frame
  frame_scene <- as.data.table(as.data.frame(predicted_scene, 
                                             cells = TRUE, 
                                             na.rm = TRUE))
  ordination <- as.matrix(frame_scene[, 2:4]/1000)
  
  # Order coefficients for naming
  coefficients <- coefficients[order(plant)]
  cwm_traits <- unique(coefficients$plant)
  coefficients$plant <- factor(coefficients$plant, levels = cwm_traits)
  
  #---------------------------------------------------------
  # Estimate CWM based on coefficients
  plant_cwm <- multi_cwm(ordination,
                         as.matrix(coefficients[, 2:5]),
                         as.numeric(coefficients$plant))
  
  #---------------------------------------------------------
  #Scaling results
  
  for(i in 1:length(cwm_traits)) {
    scale <- scaling[plant == cwm_traits[i]]
    plant_cwm[[1]][,i] <- plant_cwm[[1]][,i]*scale$scale[1]
    plant_cwm[[2]][,i] <- plant_cwm[[2]][,i]*scale$scale[1]*1000
  }
  
  #--------------------------------------------------------
  #Get a scene frame
  
  cwm_scene_mean <- subset(predicted_scene, 1)
  cwm_scene_mean[] <- NA
  cwm_scene_mean <- rep(cwm_scene_mean, length(cwm_traits))
  cwm_scene_error <- cwm_scene_mean
  cwm_scene_mean[frame_scene$cell] <- plant_cwm[[1]]
  cwm_scene_error[frame_scene$cell] <- plant_cwm[[2]]
  names(cwm_scene_mean) <- cwm_traits
  names(cwm_scene_error) <- cwm_traits
  
  #-----------------------------------------------------------
  # Scene export
  
  scene_name <- strsplit(functional_composition[X], "/")[[1]]
  scene_name <- scene_name[length(scene_name)]
  
  if(sensor == "DESIS") {
    scene_name <- strsplit(scene_name, "DESIS-L3A")[[1]][2]
    scene_name <- paste0("DESIS-L4C", scene_name)
    
  } else if(sensor == "EMIT") {
    scene_name <- strsplit(scene_name, "EMIT-L3C")[[1]][2]
    scene_name <- paste0("DESIS-L4C", scene_name)
  }
  
  scene_name <- strsplit(scene_name, "NDMS-functional-composition.tif")[[1]][1]
  scene_name_mean <- paste0(scene_name, "CWM_mean.tif")
  scene_name_error <- paste0(scene_name, "CWM_error.tif")
  
  output_raster_mean <- paste0(output_path, "/", scene_name_mean)
  output_raster_error <- paste0(output_path, "/", scene_name_error)
  
  #Export
  writeRaster(x = cwm_scene_mean,
              filename = output_raster_mean,
              names = cwm_traits,
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
  writeRaster(x = cwm_scene_error,
              filename = output_raster_error,
              names = cwm_traits,
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
}

#'------------------------------------------------------------------------------
#' @example

################################################################################
#' @Local predictions from taxonomic-composition.tif files
################################################################################

# Define root path, sensor, and coefficients
root_path <- "/media/antonio/Work/MMDB"
sensor <- "DESIS"
type <- "func"

# Read coefficients
coefficients <- fread(paste0(root_path, "/07_model_CWM/", 
                             type, "_model-CWM.csv"))
coefficients <- coefficients[dataset == "training", 
                             c("plant", "(Intercept)", "Axis1", "Axis2", "Axis3")]

# Read scaling factor
scaling <- fread(paste0(root_path, "/07_model_CWM/CWM_scale.csv"))

# Define inputs and outputs folder
folder_path <- paste0(root_path, "/00_imagery/Predicted/", sensor, "/L3/Functional")
output_path <- paste0(root_path, "/00_imagery/Predicted/", sensor, "/L4/CWM")

# Select files to process
functional_composition <- list.files(path = folder_path, 
                                     pattern = "functional-composition.tif", 
                                     full.names = TRUE,
                                     recursive = TRUE)

# Parallel processing
run <- pbmclapply(X = 1:length(functional_composition), 
                  FUN = apply_model_cwm, 
                  sensor = sensor,
                  functional_composition = functional_composition, 
                  coefficients = coefficients, 
                  scaling = scaling,
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