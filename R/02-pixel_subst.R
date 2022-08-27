#'------------------------------------------------------------------------------
#' @title Subtraction of pixels based on low variability 
#'------------------------------------------------------------------------------

#' @description Subtract pixels from Geotiffs based on low variability
#' @return data.table with names according the Geotiff file, pixel coordinates, 
#' and spectral bands

#'------------------------------------------------------------------------------
#' Library

library(terra)
library(data.table)
library(RcppHNSW)

#'------------------------------------------------------------------------------
#' Arguments
#' @param root_path Root path of the cropped rasters
#' @param out_path Path to export the outputs
#' @param n_pixels Number of pixels to extract
#' @param NDVI_threshold Threshold of NDVI to mask

#'------------------------------------------------------------------------------
#' @example 

root_path <- "/media/antonio/antonio_ssd/mapping_species/level3_crop"
out_path <- "/media/antonio/antonio_ssd/mapping_species/level4_pixel-selection"

#'------------------------------------------------------------------------------
#' Function

pixel_subst <- function(root_path, out_path, n_pixels = 10, NDVI_threshold = 0.5) {
  
  #Search for paths
  files <- list.files(path = paste0(root_path), 
                      pattern = ".tif", 
                      all.files = TRUE,
                      full.names = FALSE, 
                      recursive = TRUE)
  
  #Arrange path in frame
  frame <- data.table(matrix(unlist(strsplit(files, "_")), 
                             nrow= length(files), 
                             byrow=TRUE), stringsAsFactors=FALSE)
  
  #Colnames
  colnames(frame) <- c("x_title", "y_title", "stand_key", "mn_sp")
  
  #Species
  mn_sp <- frame[, strsplit(mn_sp, ".tif")[[1]][1], 
                 by = seq_along(1:nrow(frame))]
  
  frame$mn_sp <-  mn_sp$V1
  
  #For data export
  results <- data.table()
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = nrow(frame), style = 3)
  
  #Function on each raster
  for(i in 1:nrow(frame)) {
    
    #Progress
    setTxtProgressBar(pb, i)
    
    #Read file
    raster <- rast(paste0(root_path, "/", files[i]))
    
    #Transform to frame
    rast_frame <- as.data.frame(raster, 
                                xy = TRUE,
                                cells = TRUE,
                                na.rm = TRUE)
    
    #Get NDVI and apply threshold
    NDVI <- (rast_frame$NIR-rast_frame$RED) / (rast_frame$NIR + rast_frame$RED)
    
    #Subsample frame
    rast_frame <- rast_frame[NDVI >= NDVI_threshold, ]
    
    #Minimun number of pixels
    if(nrow(rast_frame) <= n_pixels) {
      next
    }
    
    #Frame for PCA
    frame_pca <- rast_frame[, 4:13]
    
    #Apply PCA
    apply_pca <- prcomp(frame_pca, center = TRUE, scale. = FALSE)
    apply_pca <- apply_pca$x
    
    #Get center points of density
    density_points <- NA
    
    #Find the pick of density
    for(ii in 1:ncol(apply_pca)) {
      
      den <- stats::density(apply_pca[, ii], adjust = 1)
      max_y <- which.max(den$y)
      max_x <- den$x[max_y]
      
      density_points[ii] <- max_x
      
    }
    
    #Transform to matrix
    density_points <- matrix(density_points, ncol = ncol(apply_pca))
    
    #Nearest neighbor index
    to_search <- hnsw_build(apply_pca)
    pixels_nn <- hnsw_search(density_points, to_search, n_pixels)
    pixels_nn <- pixels_nn$idx[1,]
    
    #Get pixels
    pixels <- as.data.table(rast_frame[pixels_nn,])
    
    #Merge information
    pixels <- cbind(frame[i, ], pixels)
    
    #Export
    name_export <- paste0(out_path, "/",
                          frame$x_title[i], "_",
                          frame$y_title[i], "_",
                          frame$stand_key[i], "_",
                          frame$mn_sp[i], ".txt")
    fwrite(pixels, name_export, sep = "\t")
    
    #Merge results
    results <- rbind(results, pixels)
    
  }
  
  #Close progress
  close(pb)

  #Return
  return(results)
  
}

results <- pixel_subst(root_path, out_path)
fwrite(results, paste0(out_path, "/results.txt"), sep = "\t")


