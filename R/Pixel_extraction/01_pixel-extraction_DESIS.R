#'------------------------------------------------------------------------------
#' @title Pixel Subtraction
#'------------------------------------------------------------------------------

#' @description Subtract pixels from Geotiffs of DESIS L2. During the extraction
#' use a file of coordinates and plot labels (PLT_CN, "LON", "LAT") and transform
#' the extracted bands using three methods.
#' 
#' @return several .txt files with names according to the DESIS file of pixels 
#' transformed.

#'------------------------------------------------------------------------------
#' @Library

library(terra)
library(data.table)
library(utils)
library(parallel)
library(bit64)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' @Load-source-code

source("R/auxiliary/spectral_transform.R")

#'------------------------------------------------------------------------------
#' @Arguments

#' @param scene_path Path to the root folder where the scenes are located (compressed).
#' @param FIA_path Path where the FIA plot file with labels is located.
#' @param out_path Path to export the extracted and transformed pixels.
#' @param DESIS_bands A vector with the center of bands from DESIS

#'------------------------------------------------------------------------------
#' @example 

# List of folders: "Appalachians", 
#                  "Great_Lakes", 
#                  "Mid_Atlantic", 
#                  "Northeast"
#                  "Ozarks", 
#                  "Southeast"

root_path <- "E:/DESIS"
folder <- "Appalachians"
scene_path <- paste0(root_path, "/", folder)
FIA_path <- "F:/Projects/DESIS_mapping/00_FIA/plots_eastern_20230427.csv"
out_path <- paste0("E:/pixel_extraction/DESIS/", folder)

rast_subst_DESIS(scene_path, 
                 FIA_path, 
                 out_path, 
                 DESIS_bands,
                 DESIS_mask,
                 threads = 1) # No threads in windows

#'------------------------------------------------------------------------------
#' @Function

rast_subst_DESIS <- function(scene_path, 
                       FIA_path, 
                       out_path, 
                       DESIS_bands = DESIS_bands,
                       DESIS_mask = DESIS_mask,
                       threads = 6) {
  
  ### Vector file processing 
  # Read the FI
  FIA <- fread(FIA_path)
  FIA <- FIA[, c("PLT_CN", "LAT", "LON")]
  
  # Create a virtual vector file
  FIA_vector <- vect(x = FIA, 
                     geom=c("LON", "LAT"), 
                     crs = "EPSG:4269", 
                     keepgeom = FALSE)
  
  ### Raster files processing
  #List files
  files <- dir(scene_path)
  
  # Arrange path in a frame
  frame <- data.table(scene = files)
  
  # Unique scenes
  n_scenes <- nrow(frame)
  
  #Function for pixel extraction
  extraction <- function(X, 
                         frame = frame, 
                         scene_path = scene_path, 
                         out_path = out_path, 
                         FIA_vector = FIA_vector, 
                         DESIS_bands = DESIS_bands,
                         DESIS_mask = DESIS_mask) {
    
    # Select scene of interest
    scene <- frame$scene[X]

    # Select products
    spectra <- paste0(scene_path, "/", 
                      scene, "/",
                      scene, "-", 
                      "SPECTRAL_IMAGE.tif")
    
    quality_cube <- paste0(scene_path, "/", 
                           scene, "/",
                           scene, "-", 
                           "QL_QUALITY.tif")
    
    quality_flags <- paste0(scene_path, "/", 
                            scene, "/",
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
      
      # Remove pixels outside of the scene
      remove <- rowMeans(point_quality_cube[, -1], na.rm = TRUE)
      point_spectra <- point_spectra[!is.nan(remove),]
      point_quality_cube <- point_quality_cube[!is.nan(remove),]
      point_quality_flags <- point_quality_flags[!is.nan(remove),]
      
      # If there are points at all
      if(nrow(point_spectra) > 0) {
        
        # Test values for bands
        for(j in 1:nrow(point_spectra)) {
          for(k in 2:ncol(point_spectra)) {
            if(point_quality_cube[j, k] == 32) {
              point_spectra[j, k] <- NA
            }
          }
        }
        
        # Interpolation of potential bad pixels
        spectra_interpolated <- cbind(point_spectra[, 1], 
                                      interpolation(point_spectra[, 2:ncol(point_spectra)]))
        
        ### Perform transformations
        VN <- cbind(point_quality_flags, vector_normalization(spectra_interpolated[, 13:207]))
        SG <- cbind(point_quality_flags, smooth(spectra_interpolated[, 2:ncol(spectra_interpolated)]))
        SW <- cbind(point_quality_flags, CWT(spectra_interpolated[, 2:ncol(spectra_interpolated)], 
                                             sampling = 3.0, 
                                             scales = c(3, 4, 5),
                                             sensor_bands = DESIS_bands))

        ### Export 
        # Base file name
        export_name <- paste0(out_path, "/",
                              scene)
        
        # Export files
        fwrite(VN, paste0(export_name, "_VN.txt"), sep = "\t")
        fwrite(SG, paste0(export_name, "_SG.txt"), sep = "\t")
        fwrite(SW, paste0(export_name, "_SW.txt"), sep = "\t")
        
      }
    }
  }
  
  # On Windows, set cores to be 1
  if (.Platform$OS.type == "windows") {
    cores = 1
  } else {
    cores = threads
  }
  
  # Batch application
  pbmclapply(X = 1:n_scenes, 
             FUN = extraction, 
             frame = frame, 
             scene_path = scene_path, 
             out_path = out_path, 
             FIA_vector = FIA_vector, 
             DESIS_bands = DESIS_bands,
             DESIS_mask = DESIS_mask,
             mc.cores = cores,
             mc.preschedule = FALSE,
             mc.cleanup = FALSE)

}



