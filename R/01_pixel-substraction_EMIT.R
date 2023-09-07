#'------------------------------------------------------------------------------
#' @title Pixel Subtraction
#'------------------------------------------------------------------------------

#' @description Subtract pixels from ENVI files of EMIT L2. During the extraction
#' use a file of coordinates and plot labels (PLT_CN, "LON", "LAT") and transform
#' the extracted bands using three methods.
#' 
#' @return several .txt files with names according to the EMIT file of pixels 
#' transformed.

#'------------------------------------------------------------------------------
#' Library

library(terra)
library(data.table)
library(parallel)
library(bit64)
library(pbmcapply)
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

root_path <- "E:"
folder <- "EMIT"
scene_path <- paste0(root_path, "/", folder)
FIA_path <- "F:/Projects/DESIS_mapping/00_FIA/plots_eastern_20230427.csv"
out_path <- "E:/pixel_extraction/EMIT"

rast_subst_EMIT(scene_path,
                FIA_path,
                out_path,
                EMIT_bands = EMIT_bands,
                EMIT_mask = EMIT_mask,
                EMIT_good_bands = EMIT_good_bands,
                threads = 1) # No threads in windows

#'------------------------------------------------------------------------------
#' Function
rast_subst_EMIT <- function(scene_path,
                            FIA_path,
                            out_path,
                            EMIT_bands = EMIT_bands,
                            EMIT_mask = EMIT_mask,
                            EMIT_good_bands = EMIT_good_bands,
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
  files <- list.files(scene_path)
  
  # Arrange path in a frame
  frame <- data.table(scene = files)
  
  # Get product information
  frame[, type := strsplit(scene, "_")[[1]][3], by = seq_along(1:nrow(frame))]
  frame[, time := strsplit(scene, "_")[[1]][5], by = seq_along(1:nrow(frame))]
  frame[, band_mask := strsplit(scene, "_")[[1]][8], by = seq_along(1:nrow(frame))]
  frame[, keep := unlist(gregexpr('.hdr', scene))[1], by = seq_along(1:nrow(frame))]
  frame <- frame[keep == -1, 1:4]
  
  list_scenes <- data.table(time = unique(frame$time))
  list_scenes$ID <- 1:nrow(list_scenes)
  frame <- merge(frame, list_scenes, by = "time", all.x = TRUE, all.y = TRUE)
  
  # Unique scenes
  n_scenes <- length(unique(frame$time))
  
  #Function for pixel extraction
  extraction <- function(X, 
                         frame = frame, 
                         scene_path = scene_path, 
                         out_path = out_path, 
                         FIA_vector = FIA_vector, 
                         EMIT_bands = EMIT_bands,
                         EMIT_mask = EMIT_mask,
                         EMIT_good_bands = EMIT_good_bands) {
    
    # Select scene of interest
    scene <- subset(frame, ID == X)
    
    # Select products
    spectra <- paste0(scene_path, "/", 
                      scene[type == "RFL"]$scene)
    
    mask <- paste0(scene_path, "/", 
                   scene[type == "MASK" & band_mask == "mask"]$scene)
    
    # Read products
    rast_spectra <- rast(spectra)
    rast_mask <- rast(mask)
    
    # Test for projections
    if(crs(FIA_vector) != crs(rast_mask)) {
      points <- project(FIA_vector, rast_mask)
    } else {
      points <- FIA_vector
    }
    
    # Crop points
    e <- ext(rast_mask)
    points <- crop(points, e)
    
    # If there is overlap
    if(length(points) != 0) {
      
      # Pixel extraction
      point_spectra <- as.data.frame(extract(rast_spectra, points, cells=FALSE, xy=FALSE))
      point_mask <- as.data.frame(extract(rast_mask, points, cells=FALSE, xy=FALSE))
      
      # Colnames
      colnames(point_spectra) <- c("PLT_CN", round(EMIT_bands, 2))
      colnames(point_mask) <- c("PLT_CN", EMIT_mask)
      
      #Change ID
      point_spectra$PLT_CN <- points$PLT_CN
      point_mask$PLT_CN <- points$PLT_CN
      
      # Check if that is inside
      outside <- rowMeans(point_spectra[, -1])
      point_spectra <- point_spectra[outside != 0,]
      point_mask <- point_mask[outside != 0,]
      
      # If there are points at all
      if(nrow(point_spectra) > 0) {
        
        # Remove bad bands
        good <- as.character(round(EMIT_bands[EMIT_good_bands != 0], 2))
        
        # Perform transformations
        VN <- cbind(point_mask, vector_normalization(point_spectra[, good]))
        SG <- cbind(point_mask, smooth(point_spectra[, 2:ncol(point_spectra)]))
        SW <- cbind(point_mask, CWT(point_spectra[, 2:ncol(point_spectra)], 
                                    sampling = 7.43,
                                    scales = c(2, 3),
                                    sensor_bands = EMIT_bands))
        ### Export 
        # Base file name
        export_name <- paste0(out_path, "/",
                              scene[type == "RFL", scene])
        
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
             EMIT_bands = EMIT_bands,
             EMIT_mask = EMIT_mask,
             EMIT_good_bands = EMIT_good_bands,
             mc.cores = cores,
             mc.preschedule = FALSE,
             mc.cleanup = FALSE)
  
}