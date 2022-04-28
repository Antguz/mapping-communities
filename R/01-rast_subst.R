#'------------------------------------------------------------------------------
#' @title Raster subtraction
#'------------------------------------------------------------------------------

#' @description Subtract pixels from Geotiffs in titles index and a vector file
#' @return Geotif with names according the vector file

#'------------------------------------------------------------------------------
#' Library

library(terra)
library(data.table)

#'------------------------------------------------------------------------------
#' Arguments
#' @param title
#' @param root_path
#' @param base_name
#' @param fsi
#' @param out_path
#' 
 
#'------------------------------------------------------------------------------
#' @example 

title <- vect("data/title_clipped.gpkg")
root_path <- "/media/antonio/antonio_ssd/mapping_species/level3_BAP"
base_name <- "20210715_LEVEL3_SEN2L_BAP.tif"
fsi <- vect("data/FSI.gpkg")
out_path <- "/media/antonio/antonio_ssd/mapping_species/level3_crop"

#'------------------------------------------------------------------------------
#' Function

rast_subst <- function(title, root_path, base_name, fsi, out_path) {
  
  #Check projection
  fsi <- project(fsi, title)
  
  #See interception of layers
  title_fsi <- intersect(fsi, title)
  
  #As frame
  fsi_frame <- as.data.table(title_fsi)
  
  #See in how many titles the vector falls
  n_titles <- fsi_frame[, .(count = .N), by = STAND_KEY]
  
  #Merge to frame
  fsi_frame <- merge(fsi_frame, n_titles, by = "STAND_KEY", all.x = TRUE)
  
  #Set progress bar
  pb <- txtProgressBar(min = 0, max = length(unique(fsi$STAND_KEY)), style = 3)
  n <- 1
  
  #Do vectors that fall in a unique title first---------------------------------
  
  #One intercept
  one_intercep <- subset(fsi_frame, count == 1)
  
  #N titles to read
  n_titles <- unique(one_intercep$Tile_ID)
  
  #Loop over the number of titles
  for(i in 1:length(n_titles)) {
    
    # Name of rast to read
    rast_name <- paste0(root_path, "/", n_titles[i], "/", base_name)
    
    #Read raster
    raster <- rast(rast_name)
    
    #Subset vector
    vector <- subset(title_fsi, title_fsi$Tile_ID == n_titles[i])
    
    #Unique values
    stand <- unique(vector$STAND_KEY)
    
    #Loop over vectors
    for(ii in 1:length(stand)) {
      
      # Update progress bar
      setTxtProgressBar(pb, n)
      n <- n + 1
      
      #Select plot
      plot <- subset(fsi, fsi$STAND_KEY == stand[ii])
      
      #Crop raster
      raster_crop <- crop(raster, plot, snap = "in", mask = TRUE)
      
      #Name to export
      name_export <- paste0(out_path, "/", 
                            n_titles[i], "_", 
                            stand[ii], "_",
                            plot$MN_SPP, ".tif")
      
      #Export
      writeRaster(raster_crop, 
                  name_export, 
                  overwrite = TRUE,
                  names = names(raster_crop))
      
      rm(list = c("plot", "raster_crop", "name_export"))
      
    }
    
    rm(list = c("rast_name", "raster", "vector", "stand"))
    
  }
  
  #Do vectors that fall in several tiles ---------------------------------------
  
  #Two or more intercepts
  tom_intercept <- subset(fsi_frame, count > 1)
  
  #N stands to read
  stands <- unique(tom_intercept$STAND_KEY)
  
  for(i in 1:length(stands)) {
    
    #Subset stands
    sub_stand <- subset(tom_intercept, STAND_KEY == stands[i])
    
    #N of tiles
    titles <- unique(sub_stand$Tile_ID)
    
    #Empty raster
    raster <- NA
    
    for(ii in 1:length(titles)) {
      
      if(ii == 1) {
        
        # Name of rast to read
        sub_raster_name <- paste0(root_path, "/", titles[ii], "/", base_name)
        
        #Read raster
        raster <- rast(sub_raster_name)
        
      } else {
      
      # Name of rast to read
      sub_raster_name <- paste0(root_path, "/", titles[ii], "/", base_name)
      
      #Read raster
      sub_raster <- rast(sub_raster_name)
      
      #Mosaic
      raster <- terra::mosaic(raster, sub_raster)
      
      }
    }
    
    # Update progress bar
    setTxtProgressBar(pb, n)
    n <- n + 1
    
    #Select plot
    plot <- subset(fsi, fsi$STAND_KEY == stands[i])
    
    #Crop raster
    raster_crop <- crop(raster, plot, snap = "in", mask = TRUE)
    
    #Name to export
    name_export <- paste0(out_path, "/", 
                          titles[1], "_", 
                          stand[i], "_",
                          plot$MN_SPP, ".tif")
    
    #Export
    writeRaster(raster_crop, 
                name_export, 
                overwrite = TRUE,
                names = names(raster_crop))
    
    rm(list = c("plot", "raster_crop", "name_export", 
                "sub_stand", "titles", "raster"))
    
  }
}
