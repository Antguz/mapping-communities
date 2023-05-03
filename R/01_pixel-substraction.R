#'------------------------------------------------------------------------------
#' @tile Pixel Subtraction
#'------------------------------------------------------------------------------

#' @description Subtract pixels from Geotiffs of DESIS L2. During the extraction
#' use a file of coordinates and plot labels (PLT_CN, "LON", "LAT") and transform
#' the extracted bands using three methods.
#' 
#' @return several .txt files with names according the DESIS file of pixels 
#' transformed.

#'------------------------------------------------------------------------------
#' Library

library(terra)
library(data.table)
library(utils)
library(parallel)
library(bit64)
options(scipen = 100000)

#'------------------------------------------------------------------------------
#' Load source code

source("R/spectral_transform.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param scene_path Path to the root folder where the scenes are located (compressed).
#' @param FIA_path Path where the FIA plot file with labels is located.
#' @param temp Path to a folder to storage temporal information during extraction.
#' @param out_path Path to export the extracted and transformed pixels.
#' @param DESIS_bands A vector with the center of bands from DESIS

#'------------------------------------------------------------------------------
#' @example 

# List of folders: "Appalachians", "Great_Lakes", "Mid_Atlantic", "Northeast".
# "Ozarks", "Southeast".

root_path <- "/media/antonio/Extreme SSD/Imagery"
folder <- "Appalachians"
scene_path <- paste0("/media/antonio/Extreme SSD/Imagery", "/", folder)
temp <- paste0("/media/antonio/Extreme SSD/Imagery", "/temp")
FIA_path <- paste0(root_path, "/plots_eastern_20230427.csv")
out_path <- paste0("/media/antonio/Extreme SSD/Imagery", "/01_pixel-subtraction/", folder)

rast_subst(scene_path, FIA_path, out_path, DESIS_bands)

#'------------------------------------------------------------------------------
#' Function
rast_subst <- function(scene_path, 
                       FIA_path, 
                       out_path, 
                       DESIS_bands = DESIS_bands) {
  
  ### Vector file processing 
  # Read the FI
  FIA <- fread(FIA_path)
  FIA <- FIA[, c("PLT_CN", "LAT", "LON")]
  
  # Create a virtual vector file
  FIA_vector <- vect(x = FIA, 
                     geom=c("LON", "LAT"), 
                     crs = "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs", 
                     keepgeom = FALSE)
  
  ### Raster files processing
  #List files
  files <- list.files(path = scene_path, 
                      pattern = ".zip", 
                      recursive = TRUE,
                      include.dirs = TRUE)
  
  # Arrange path in a frame
  frame <- data.table(matrix(unlist(strsplit(files, ".zip")), 
                             nrow= length(files), 
                             byrow=TRUE), 
                      stringsAsFactors=FALSE)
  colnames(frame) <- c("scene")
  
  # Unique scenes
  n_scenes <- nrow(frame)
  
  #Function for pixel extraction
  extraction <- function(X, 
                         frame = frame, 
                         scene_path = scene_path, 
                         out_path = out_path, 
                         temp = temp, 
                         FIA_vector = FIA_vector, 
                         DESIS_bands = DESIS_bands) {
    
    scene <- frame$scene[X]
    
    #Unzip files
    zip_scene <- paste0(scene_path, "/", scene, ".zip")
    out_zip <- paste0(temp, "/", scene)
    unzip(zip_scene, exdir = out_zip)

    # Select products
    spectra <- paste0(out_zip, "/", 
                      scene, "-", 
                      "SPECTRAL_IMAGE.tif")
    
    quality_cube <- paste0(out_zip, "/", 
                           scene, "-", 
                           "QL_QUALITY.tif")
    
    quality_flags <- paste0(out_zip, "/", 
                            scene, "-", 
                            "QL_QUALITY-2.tif")
    
    # Read products
    rast_spectra <- rast(spectra)
    rast_quality_cube <- rast(quality_cube)
    rast_quality_flags <- rast(quality_flags)
    
    # Test for projections
    if(crs(FIA_vector) != crs(rast_quality_flags)) {
      points <- project(FIA_vector, rast_quality_flags)
    } else {
      points <- FIA_vector
    }
    
    #Crop
    e <- ext(rast_quality_flags)
    points <- crop(points, e)
    
    # If there is overlap
    if(length(points) != 0) {
      
      # Pixel extraction
      point_spectra <- as.data.frame(extract(rast_spectra, points, cells=FALSE, xy=FALSE))
      point_quality_cube <- as.data.frame(extract(rast_quality_cube, points, cells=FALSE, xy=FALSE))
      point_quality_flags <-as.data.frame(extract(rast_quality_flags, points, cells=FALSE, xy=FALSE))
      
      # Colnames
      colnames(point_spectra) <- c("PLT_CN", DESIS_bands)
      colnames(point_quality_cube) <- c("PLT_CN", DESIS_bands)
      
      #Change ID
      point_spectra$PLT_CN <- points$PLT_CN
      point_quality_cube$PLT_CN <- points$PLT_CN
      colnames(point_quality_flags)[1] <- "PLT_CN"
      point_quality_flags$PLT_CN <- points$PLT_CN
      
      # To remove
      remove <- NA
      
      # Test values for bands
      for(j in 1:nrow(point_spectra)) {
        for(k in 2:ncol(point_spectra)) {
          if(is.na(point_quality_cube[j, k])) {
            remove[j] <- TRUE
          } else if(point_quality_cube[j, k] == 32) {
            point_spectra[j, k] <- NA
          }
        }
      }
      
      # NA observations
      remove[is.na(remove)] <- FALSE
      point_spectra <- point_spectra[remove != TRUE,]
      point_quality_flags <- point_quality_flags[remove != TRUE,]
      
      # If there are points at all
      if(nrow(point_spectra) > 0) {
        
        ### Perform transformations
        VN <- cbind(point_quality_flags, vector_normalization(point_spectra[, 13:207]))
        SG <- cbind(point_quality_flags, smooth(point_spectra[, 2:ncol(point_spectra)]))
        SW <- cbind(point_quality_flags, CWT(point_spectra[, 2:ncol(point_spectra)], 
                                             sampling = 3.0, 
                                             scales = 5))
        #ME <- cbind(point_quality_flags, mean_extraction(spectra))
        #SN <- cbind(point_quality_flags, spectra_normalization(spectra))
        
        ### Export 
        # Base file name
        export_name <- paste0(out_path, "/",
                              scene)
        
        # Export files
        fwrite(VN, paste0(export_name, "_VN.txt"), sep = "\t")
        fwrite(SG, paste0(export_name, "_SG.txt"), sep = "\t")
        fwrite(SW, paste0(export_name, "_SW.txt"), sep = "\t")
        #fwrite(ME, paste0(export_name, "_ME.txt"), sep = "\t")
        #fwrite(SN, paste0(export_name, "_SN.txt"), sep = "\t")
        
      }
    }
    
    #Remove temporal directory
    unlink(out_zip, recursive = TRUE)
    return(X)
    
  }
  
  #Add progress bar
  pb <- txtProgressBar(min = 1, max = n_scenes, style = 3)
  
  for(i in 1:n_scenes) {
    
    setTxtProgressBar(pb, i)
    
    extraction(X = i,
               frame = frame, 
               scene_path = scene_path, 
               out_path = out_path, 
               temp = temp, 
               FIA_vector = FIA_vector, 
               DESIS_bands = DESIS_bands)
  }
}



