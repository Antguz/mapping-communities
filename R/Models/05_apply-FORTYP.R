#'------------------------------------------------------------------------------
#' @title Apply the occurrence probability of forest type to scenes
#'------------------------------------------------------------------------------

#' @description It apply the occurrence (i.e., probability) of forest types to the
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
#' @param X
#' @param sensor
#' @param taxonomic_composition
#' @param coefficients
#' @param output_path

#'------------------------------------------------------------------------------
#' Function

apply_model_occurence <- function(X,
                                  sensor,
                                  taxonomic_composition, 
                                  coefficients, 
                                  output_path) {
  
  # Read predicted scene
  predicted_scene <- rast(taxonomic_composition[X])
  
  #Ask data frame
  frame_scene <- as.data.table(as.data.frame(predicted_scene, cells = TRUE, na.rm = TRUE))
  ordination <- as.matrix(frame_scene[, 2:4]/1000)
  
  # Reshape coefficients
  models <- unique(coefficients$forty)
  coef_array <- array(dim = c(100, 3, length(models)))
  intercept_mat <- matrix(nrow = 100, ncol = length(models))
  
  for(i in 1:length(models)) {
    
    sub <- coefficients[forty == models[i]]
    coef_array[,,i] <- as.matrix(sub[, 3:5])
    intercept_mat[,i] <- sub$`(Intercept)`
    
  }
  
  # Order coefficients for naming
  # coefficients <- coefficients[order(plant)]
  # FORTYP <- unique(coefficients$plant)
  # coefficients$plant <- factor(coefficients$plant, levels = FORTYP)
  
  #---------------------------------------------------------
  # Estimate probabilities based on coefficients
  plant_probabilities <- summary_prob(ordination,
                                      coef_array,
                                      intercept_mat)
  
  # plant_probabilities <- multi_probability(ordination, 
  #                                          as.matrix(coefficients[, 2:5]),
  #                                          as.numeric(FORTYP))
  
  #--------------------------------------------------------
  #Get a scene frame
  fortyp_scene_mean <- subset(predicted_scene, 1)
  fortyp_scene_mean[] <- NA
  fortyp_scene_mean <- rep(fortyp_scene_mean, length(models))
  fortyp_scene_error <- fortyp_scene_mean
  fortyp_scene_mean[frame_scene$cell] <- plant_probabilities[[1]]
  fortyp_scene_error[frame_scene$cell] <- plant_probabilities[[2]]
  names(fortyp_scene_mean) <- models
  names(fortyp_scene_error) <- models
  
  #-----------------------------------------------------------
  # Scene export
  
  scene_name <- strsplit(taxonomic_composition[X], "/")[[1]]
  scene_name <- scene_name[length(scene_name)]
  
  if(sensor == "DESIS") {
    scene_name <- strsplit(scene_name, "DESIS-L3A")[[1]][2]
    scene_name <- paste0("DESIS-L4A", scene_name)

  } else if(sensor == "EMIT") {
    scene_name <- strsplit(scene_name, "EMIT-L3A")[[1]][2]
    scene_name <- paste0("EMIT-L4A", scene_name)
  }
  
  scene_name <- strsplit(scene_name, "NDMS-taxonomic-composition.tif")[[1]][1]
  scene_name_mean <- paste0(scene_name, "FORTYP_mean.tif")
  scene_name_error <- paste0(scene_name, "FORTYP_error.tif")
  
  output_raster_mean <- paste0(output_path, "/", scene_name_mean)
  output_raster_error <- paste0(output_path, "/", scene_name_error)
  
  #Export
  writeRaster(x = fortyp_scene_mean,
              filename = output_raster_mean,
              names = models,
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
  writeRaster(x = fortyp_scene_error,
              filename = output_raster_error,
              names = models,
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
sensor <- "EMIT"
type <- "taxo"

# Read coefficients
coefficients <- fread(paste0(root_path, "/Models/03_models_fortyp/", 
                             type, "_model-FORTYP.csv"))

coefficients <- coefficients[dataset == "training", 
                             c("forty", "(Intercept)", "Axis1", "Axis2", "Axis3")]
select <- c("Sugar maple / beech / yellow birch", "Red spruce / balsam fir", "Red spruce", "Red maple / upland", "Paper birch", "Northern red oak", "Eastern white pine / eastern hemlock", "Eastern hemlock", "Black spruce", "Balsam fir")
coefficients <- coefficients[forty %in% select, ]

# Define inputs and outputs folder
folder_path <- paste0(root_path, "/00_imagery/Predicted/", sensor, "/L3/L3A")
output_path <- paste0(root_path, "/00_imagery/Predicted/", sensor, "/L4/FORTYP")

# Select files to process
taxonomic_composition <- list.files(path = folder_path, 
                                       pattern = "taxonomic-composition.tif", 
                                       full.names = TRUE,
                                       recursive = TRUE)

# Parallel processing
run <- pbmclapply(X = 1:length(taxonomic_composition), 
                  FUN = apply_model_occurence, 
                  sensor = sensor,
                  taxonomic_composition = taxonomic_composition, 
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
fwrite(complete, paste0(output_path, "/files_", sensor, "_FORTYR_", Sys.time(), ".txt"))