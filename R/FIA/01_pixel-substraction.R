#'------------------------------------------------------------------------------
#' @tile Pixel Subtraction
#'------------------------------------------------------------------------------

#' @description Subtract pixels from Geotiffs of DESIS L2 from vector file
#' @return several .txt files with names according the DESIS file

#'------------------------------------------------------------------------------
#' Library

library(terra)
library(data.table)

#'------------------------------------------------------------------------------
#' Arguments
#' @param in_path
#' @param fsi_plot
#' @param out_path

#'------------------------------------------------------------------------------
#' @example 

root_path <- "E:/DESIS"
in_path <- paste0(root_path, "/00_imagery")
frame_plot <- fread(paste0(root_path, "/00_FIA/plots_WI_MN_20221102.csv"))
frame_plot <- frame_plot[, c("PLT_CN", "LAT", "LON")]
fsi_plot <- vect(frame_plot,  geom=c("LON", "LAT"), crs = "EPSG:4326")
out_path <- paste0(root_path, "/01_pixel-subtraction")
rast_subst(in_path, fsi_plot, out_path)

#'------------------------------------------------------------------------------
#' Function
rast_subst <- function(in_path, fsi_plot, out_path) {
  
  #For dinamicy projection
  fsi_vect <- fsi_plot
  
  #List files
  files <- list.files(path = in_path, 
                      pattern = ".tif", 
                      recursive = TRUE,
                      include.dirs = TRUE)
  
  #Arrange path in frame
  frame <- data.table(matrix(unlist(strsplit(files, "/")), 
                             nrow= length(files), 
                             byrow=TRUE), stringsAsFactors=FALSE)
  colnames(frame) <- c("folder", "scene", "product")
  
  #Product type
  frame[, type := unlist(strsplit(product, "V0215-"))[2], 
        by = seq_along(1:nrow(frame))]
  frame[, type := unlist(strsplit(type, ".tif"))[1], 
        by = seq_along(1:nrow(frame))]
  
  #Unique scenes
  unique_scenes <- frame[, .N, by= c("folder", "scene")]
  n_scenes <- nrow(unique_scenes)
  unique_scenes$N <- 1:n_scenes
  
  #Merge with frame
  frame <- merge(unique_scenes, frame, by = c("folder", "scene"))
  
  # Point extraction -----------------------------------------------------------
  for(i in 1:n_scenes) {
    
    #Subframe
    sub_frame <- frame[N == i]
    
    #Create name
    spectra <- paste0(in_path, "/", 
                       sub_frame$folder[1], "/", 
                       sub_frame$scene[1], "/",
                       sub_frame[type == "SPECTRAL_IMAGE"]$product[1])
    
    quality <- paste0(in_path, "/", 
                      sub_frame$folder[1], "/", 
                      sub_frame$scene[1], "/",
                      sub_frame[type == "QL_QUALITY"]$product[1])
    
    quality2 <- paste0(in_path, "/", 
                       sub_frame$folder[1], "/", 
                       sub_frame$scene[1], "/",
                       sub_frame[type == "QL_QUALITY-2"]$product[1])
    
    #Read rast
    rast_spectra <- rast(spectra)
    rast_quality <- rast(quality)
    rast_quality2 <- rast(quality2)
    
    #Test projections
    if(crs(fsi_vect) != crs(rast_quality2)) {
      fsi_vect <- project(fsi_vect, rast_quality2)
    }
    
    #Crop points
    points <- fsi_vect 
    e <- ext(rast_quality2)
    points <- crop(points, e)
    
    if(length(points) != 0) {
      
      #Extraction
      point_spectra <- as.data.frame(extract(rast_spectra, points, cells=FALSE, xy=FALSE))
      point_quality <- as.data.frame(extract(rast_quality, points, cells=FALSE, xy=FALSE))
      point_quality2 <-as.data.frame(extract(rast_quality2, points, cells=FALSE, xy=FALSE))
      
      #Colnames
      colnames(point_spectra) <- c("ID", 1:235)
      colnames(point_quality) <- c("ID", 1:235)
      
      #Change ID
      point_spectra$ID <- points$PLT_CN
      point_quality$ID <- points$PLT_CN
      point_quality2$ID <- points$PLT_CN
      
      #To remove
      remove <- NA
      
      for(j in 1:nrow(point_spectra)) {
        for(k in 2:ncol(point_spectra)) {
          
          if(is.na(point_quality[j, k])) {
            remove[j] <- TRUE
          } else if(point_quality[j, k] == 32) {
            point_spectra[j, k] <- NA
          }
        }
      }
      
      #NA observations
      remove[is.na(remove)] <- FALSE
      point_spectra <- point_spectra[remove != TRUE,]
      point_quality2 <- point_quality2[remove != TRUE,]
      
      export_spectra <- merge(point_quality2, point_spectra, by = c("ID"))
      
      #Export
      export_name <- paste0(out_path, "/",
                            sub_frame$folder[1], "/",
                            unlist(strsplit(sub_frame[type == "QL_QUALITY-2"]$product[1], "-V0215-"))[1], 
                            ".txt")
      
      fwrite(export_spectra, export_name, sep = "\t")
      
    }
  }
}
