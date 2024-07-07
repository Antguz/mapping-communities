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
library(pbmcapply)

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/probability_rcpp.cpp")
sourceCpp("src/summary_rcpp.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param predicted_scene_path
#' @param coefficients
#' @param output_raster

#'------------------------------------------------------------------------------
#' Function

apply_model_probability <- function(X,
                                    predicted_scene_path,
                                    coefficients,
                                    output_path) {
  
  # Read predicted scene
  predicted_scene <- rast(predicted_scene_path[X])
  occur_scene <- subset(predicted_scene, 1)
  occur_scene[] <- NA
  
  #Ask data frame
  frame_scene <- as.data.table(as.data.frame(predicted_scene, cells = TRUE, na.rm = TRUE))
  cell <- frame_scene$cell
  ordination <- as.matrix(frame_scene[, 2:4]/1000)
  
  #---------------------------------------------------------
  # Estimate probabilities based on coefficients
  plant_probabilities <- multi_probability(ordination, 
                                           as.matrix(coefficients[, 2:5]))
  coefficients$col <- 1:nrow(coefficients)
  unique_coeff <- unique(coefficients$plant)
  result <- data.table()
  
  for(i in 1:length(unique_coeff)) {
    
    occur <- coefficients[plant == unique_coeff[i]]
    summary_occur <- as.data.table(summary_rcpp(plant_probabilities[, occur$col]))
    colnames(summary_occur) <- c(unique_coeff[i], 
                                 paste0(unique_coeff[i], "_amplitud"))
    result <- cbind(result, summary_occur)
    
  }
  
  result <- as.data.frame(result)
  
  # Remove elements
  rm(list = c("predicted_scene", "frame_scene", "ordination", "plant_probabilities",
              "occur", "summary_occur"))
  
  #--------------------------------------------------------
  #Get a scene frame
  occur_scene <- rep(occur_scene, length(unique_coeff)*2)
  occur_scene[cell] <- result
  names(occur_scene) <- colnames(result)
  
  # Remove elements
  rm(list = c("result", "cell", "unique_coeff"))
  
  #-----------------------------------------------------------
  # Scene export
  
  scene_name <- strsplit(predicted_scene_path[X], "/")[[1]]
  scene_name <- scene_name[length(scene_name)]
  fl <- scene_name
  type <- strsplit(scene_name, "_NDMS-")[[1]][2]
  
  if(type == "phylogenetic-composition.tif") {
    
    scene_name <- strsplit(scene_name, "-phylogenetic-composition.tif")[[1]][1]
    scene_name <- paste0(scene_name, "-lineages.tif")
    scene_name <- strsplit(scene_name, "DESIS-L3A_")[[1]][2]
    scene_name <- paste0("DESIS-L4B_", scene_name)
    output_raster <- paste0(output_path, "/", scene_name)
    
  } else if(type == "taxonomic-composition.tif") {
    
    scene_name <- strsplit(scene_name, "-phylogenetic-composition.tif")[[1]][1]
    scene_name <- paste0(scene_name, "-fortyp.tif")
    scene_name <- strsplit(scene_name, "DESIS-L3A_")[[1]][2]
    scene_name <- paste0("DESIS-L4A_", scene_name)
    output_raster <- paste0(output_path, "/", scene_name)
    
  }
  
  #Export
  writeRaster(x = occur_scene,
              filename = output_raster,
              names = names(occur_scene),
              NAflag = -9999,
              overwrite= TRUE,
              datatype = "INT4S")
  
  # Return
  gc()
  return(paste0(fl, ": success"))
  
}

#'------------------------------------------------------------------------------
#' @example

# Run in a batch
root_path <- "/media/antonio/Work/DESIS_mapping"
sensor <- "DESIS"
type <- "phylo-log"
coefficients <- fread(paste0(root_path, "/07_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- coefficients[dataset == "training", c("plant", "(Intercept)", "Axis1", "Axis2", "Axis3")]
files <- list.files(path = paste0(root_path, "/00_imagery/predicted/", sensor),
                    pattern = "phylogenetic-composition.tif",
                    full.names = TRUE)
output_path <- paste0(root_path, "/00_imagery/occurrence/", sensor)

# Parallel processing
run <- pbmclapply(X = 1:length(files), 
                  FUN = apply_model_probability, 
                  predicted_scene_path = files, 
                  coefficients = coefficients,
                  output_path = output_path,
                  mc.preschedule = TRUE, 
                  mc.set.seed = FALSE,
                  mc.silent = FALSE, 
                  mc.cores = 4,
                  mc.cleanup = TRUE, 
                  mc.allow.recursive = TRUE, 
                  affinity.list = NULL)

complete <- do.call(rbind, run)
fwrite(complete, paste0(output_path, "/files ", Sys.time(), ".txt"))
