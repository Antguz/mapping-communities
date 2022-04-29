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
library(Rcpp)

#'------------------------------------------------------------------------------
#' Load source code

sourceCpp("src/euclidean_distance.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' @param root_path 
#' @param out_path
#' 

#'------------------------------------------------------------------------------
#' @example 

root_path <- "F:/mapping_species/level3_crop"
base_name <- "20210715_LEVEL3_SEN2L_BAP.tif"

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
  
  #Function on each raster
  for(i in 1:nrow(frame)) {
    
    #Read file
    raster <- rast(paste0(root_path, "/", files[i]))
    
    #Transform to frame
    frame <- as.data.frame(raster, 
                           xy = TRUE,
                           cells = TRUE,
                           na.rm = TRUE)
    
    #Get NDVI and apply threshold
    NDVI <- (frame$NIR-frame$RED) / (frame$NIR + frame$RED)
    
    #Subsample frame
    frame <- frame[NDVI >= NDVI_threshold, ]
    
    #Frame for PCA
    frame_pca <- frame[, 4:13]
    
    #Apply PCA
    apply_pca <- prcomp(frame_pca, center = TRUE, scale. = FALSE)
    apply_pca <- apply_pca$x
    
    density_points <- NA
    
    for(ii in 1:ncol(apply_pca)) {
      
      den <- density(apply_pca[, ii])
      max_y <- which.max(den$y)
      max_x <- den$x[max_y]
      
      density_points[ii] <- max_x
      
    }
    
    apply_pca <- as.data.table(apply_pca)
    
    #Estimate distances
    distances <- euclidean_distance(density_points, as.matrix(apply_pca))
    distances <- order(distances)
  
    
    #Get pixels
    pixels <- frame[(distances <= n_pixels) == TRUE,]
    
  }
  
  
}

plot(apply_pca$PC1, apply_pca$PC2)
points(density_points[1], density_points[2], col = "red")
pca <- apply_pca[(distances <= n_pixels) == TRUE,]
points(pca$PC1, pca$PC2, col = "blue", cex = 2)
