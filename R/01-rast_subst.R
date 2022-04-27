#'------------------------------------------------------------------------------
#' @title Raster substraction
#'------------------------------------------------------------------------------

#' @description Substract pixels from rasters in tiles and a vector file
#' @return Croping rasters with names according the vector file

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
root_path <- "F:/move_here/BAP/home/cavender/shared/oak-wilt/level3_BAP"
base_name <- "20210715_LEVEL3_SEN2L_BAP.tif"
fsi <- vect("data/FSI.gpkg")
out_path <- "C:/Users\guzman\Documents\Goal_2\GIS\MN_forest_inventory\Modify"

#'------------------------------------------------------------------------------
#' Function

rast_subst <- function(title, root_path, base_name, fsi, out_path) {
  
  #Check projection
  if(all.equal(crs(title), crs(fsi))) {
    fsi <- project(fsi, title)
  }
  
  #See interception of layers
  title_fsi <- intersect(fsi, title)
  
  #As frame
  fsi_frame <- as.data.table(title_fsi)
  
  #See in how many titles the vector falls
  n_titles <- fsi_frame[, .(count = .N), by = STAND_KEY]
  
  #Merge to frame
  fsi_frame <- merge(fsi_frame, n_titles, by = "STAND_KEY", all.x = TRUE)
  
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
      
      #Select plot
      plot <- subset(fsi, fsi$STAND_KEY == stand[ii])
      
      #Crop raster
      raster_crop <- crop(raster, plot, snap = "in", mask = TRUE)
      
      #Name to export
      name_export <- paste0(out_path, "/", 
                            n_titles[i], "_", 
                            stand[ii], "_",
                            plot$MN_SPP, "tif")
      
    }
    
  }
  
  
  #Do vectors that fall in several tiles ---------------------------------------
  
  
}
