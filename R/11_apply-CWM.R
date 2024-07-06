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

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/cwm_rcpp.cpp")
sourceCpp("src/summary_rcpp.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param predicted_scene_path
#' @param coefficients
#' @param output_raster
#' @param threads

#'------------------------------------------------------------------------------
#' Function

apply_model_CWM <- function(predicted_scene_path, 
                            coefficients, 
                            output_path) {
  
  # Read predicted scene
  predicted_scene <- rast(predicted_scene_path)
  cwn_scene <- subset(predicted_scene, 1)
  cwn_scene[] <- NA
  
  #Ask data frame
  frame_scene <- as.data.table(as.data.frame(predicted_scene, cells = TRUE, na.rm = TRUE))
  cell <- frame_scene$cell
  ordination <- as.matrix(frame_scene[, 2:4]/1000)
  
  #---------------------------------------------------------
  # Estimate probabilities based on coefficients
  plant_cwm <- multi_cwm(ordination,
                         as.matrix(coefficients[, 2:5]))
  
  coefficients$col <- 1:nrow(coefficients)
  unique_coeff <- unique(coefficients$plant)
  result <- data.table()
  
  for(i in 1:length(unique_coeff)) {
    
    trait <- coefficients[plant == unique_coeff[i]]
    summary_trait <- as.data.table(summary_rcpp(plant_cwm[, trait$col]))
    colnames(summary_trait) <- c(unique_coeff[i], 
                                 paste0(unique_coeff[i], "_amplitud"))
    result <- cbind(result, summary_trait)
    
  }
  
  result <- as.data.frame(result)
  
  # Remove elements
  rm(list = c("predicted_scene", "frame_scene", "ordination", "plant_cwm",
              "trait", "summary_trait"))
  
  #---------------------------------------------------------
  #Scaling factor for exporting (ugly but works)
  result[, 1] <- round(result[, 1] * 100)
  result[, 2] <- round(result[, 2] * 100000000)
  result[, 3] <- round(result[, 3] * 100)
  result[, 4] <- round(result[, 4] * 10000000)
  result[, 5] <- round(result[, 5] * 1000)
  result[, 6] <- round(result[, 6] * 1000000000)
  result[, 7] <- round(result[, 7] * 10000)
  result[, 8] <- round(result[, 8] * 10000000000)
  result[, 9] <- round(result[, 9] * 100)
  result[, 10] <- round(result[, 10] * 10000000)
  result[, 11] <- round(result[, 11] * 10000)
  result[, 12] <- round(result[, 12] * 10000000000)
  result[, 13] <- round(result[, 13] * 1000)
  result[, 14] <- round(result[, 14] * 1000000000)
  result[, 15] <- round(result[, 15] * 1000)
  result[, 16] <- round(result[, 16] * 1000000000)
  result[, 17] <- round(result[, 17] * 1000)
  result[, 18] <- round(result[, 18] * 1000000000)
  result[, 19] <- round(result[, 19] * 0.1)
  result[, 20] <- round(result[, 20] * 100000)
  result[, 21] <- round(result[, 21] * 10)
  result[, 22] <- round(result[, 22] * 10000000)
  result[, 23] <- round(result[, 23] * 100)
  result[, 24] <- round(result[, 24] * 100000000)
  
  result <- as.data.table(result)
  
  #a <- apply(result, 2, max)
  #number <- 24
  #(a[number] * 10) < 15000
  #(a[number] * 100) < 15000
  #(a[number] * 1000) < 15000
  #(a[number] * 10000) < 15000
  #(a[number] * 100000) < 15000
  #(a[number] * 1000000) < 15000
  #(a[number] * 10000000) < 15000
  #(a[number] * 100000000) < 15000
  #(a[number] * 1000000000) < 15000
  #(a[number] * 10000000000) < 15000
  
  #--------------------------------------------------------
  #Get a scene frame
  cwn_scene <- rep(cwn_scene, length(unique_coeff)*2)
  cwn_scene[cell] <- result
  names(cwn_scene) <- colnames(result)
  
  # Remove elements
  rm(list = c("result", "cell", "unique_coeff"))
  
  #-----------------------------------------------------------
  # Scene export
  
  scene_name <- strsplit(predicted_scene_path, "/")[[1]]
  scene_name <- scene_name[length(scene_name)]
  scene_name <- strsplit(scene_name, "-functional-composition.tif")[[1]][1]
  scene_name <- paste0(scene_name, "-CWM.tif")
  scene_name <- strsplit(scene_name, "DESIS-L3A_")[[1]][2]
  scene_name <- paste0("DESIS-L4C_", scene_name)
  output_raster <- paste0(output_path, "/", scene_name)
  
  #Export
  writeRaster(x = cwn_scene,
              filename = output_raster,
              names = names(cwn_scene),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
}

#'------------------------------------------------------------------------------
#' @example

# Run in a batch
#root_path <- "E:/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/DESIS_mapping"
sensor <- "DESIS"
type <- "func"
coefficients <- fread(paste0(root_path, "/07_model_CWM/", 
                             type, "_model-CWM.csv"))
coefficients <- coefficients[dataset == "training", c("plant", "(Intercept)", "Axis1", "Axis2", "Axis3")]
files <- list.files(path = paste0(root_path, "/00_imagery/predicted/", sensor),
                          pattern = "functional-composition.tif",
                          full.names = TRUE)
output_path <- paste0(root_path, "/00_imagery/CWM/", sensor, "/", type)

pb <- txtProgressBar(min = 0, max = length(files), style = 3)

for(i in 1:length(files)) {
  
  setTxtProgressBar(pb, i)
  
  apply_model_CWM(predicted_scene_path = files, 
                  coefficients = coefficients,
                  output_path)
  
}
