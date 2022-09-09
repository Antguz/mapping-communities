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

fsi <- fread("data/FIA_SPCD.csv")
climate_path <- "E:/mapping_species/climate"
out_path <- "E:/mapping_species"

#In case
options("scipen"=100, "digits"=6)

#'------------------------------------------------------------------------------
#' Function
extraction_climate <- function(fsi, climate_path, out_path) {
  
  #Clean fsi
  #Clean FIA
  fsi <- subset(fsi, STATECD == 27 | STATECD == 55)
  fsi <- fsi[order(PLOTID, SUBP)]
  
  #FIA cordinates
  fia_coordinates <- fsi[, .(LAT = mean(LAT, na.rm = TRUE), LON = mean(LON, na.rm = TRUE)), by = c("PLOTID", "SUBP")]
  fia_coordinates <- fia_coordinates[SUBP == "1"]
  
  #Make a vect
  fsi_vect <- vect(fia_coordinates, geom=c("LON", "LAT"), crs= "+proj=longlat +datum=WGS84 +no_defs +type=crs")
  
  #Search for climatic rasters ------------------------------------------------- 
  climate <- list.files(path = climate_path, 
                        pattern = ".tif", 
                        all.files = TRUE,
                        full.names = FALSE, 
                        recursive = TRUE)
  
  results <- data.table()
  
  for(i in 1:length(climate)) {
    
    raster <- rast(paste0(climate_path, "/", climate[i]))
    values <- as.data.table(extract(raster, fsi_vect))
    name <- names(values)[2]
    results <- cbind(results, values[, ..name])
    
  }
  
  fia_climate <- cbind(fia_coordinates, results)
  
  #Export
  fwrite(fia_climate, paste0(out_path, "/master_FIA-climate.txt"), sep = "\t")
  
}