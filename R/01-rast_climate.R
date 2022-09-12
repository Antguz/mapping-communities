#'------------------------------------------------------------------------------
#' @title Climate metrics
#'------------------------------------------------------------------------------

#' @description It extract climatic metrics
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(terra)

#'------------------------------------------------------------------------------
#' Arguments
#' @param master
#' @param climate_path
#' @param out_path

#'------------------------------------------------------------------------------
#' @example 

fia <- fread("/media/antonio/antonio_ssd/mapping_species/results/master_FIA-stand.txt")
climate_path <- "/media/antonio/antonio_ssd/mapping_species/climate"
out_path <- "/media/antonio/antonio_ssd/mapping_species/results"

#In case
options("scipen"=100, "digits"=6)

#'------------------------------------------------------------------------------
#' Function
extraction_climate <- function(fia, climate_path, out_path) {
  
  #Make a vect
  fia_vect <- vect(fia, geom=c("LON", "LAT"), crs= "+proj=longlat +datum=WGS84 +no_defs +type=crs")
  
  #Search for climatic rasters ------------------------------------------------- 
  climate <- list.files(path = climate_path, 
                        pattern = ".tif", 
                        all.files = TRUE,
                        full.names = FALSE, 
                        recursive = TRUE)
  
  results <- data.table()
  
  for(i in 1:length(climate)) {
    
    raster <- rast(paste0(climate_path, "/", climate[i]))
    values <- as.data.table(extract(raster, fia_vect))
    name <- names(values)[2]
    results <- cbind(results, values[, ..name])
    
  }
  
  fia_climate <- cbind(fia[, c(1:4)], results)
  
  #Export
  fwrite(fia_climate, paste0(out_path, "/master_FIA-climate.txt"), sep = "\t")
  
}