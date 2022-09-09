#'------------------------------------------------------------------------------
#' @tile Raster subtraction
#'------------------------------------------------------------------------------

#' @description Subtract pixels from Geotiffs in tiles index and a vector file
#' @return .txt files with names according the vector file

#'------------------------------------------------------------------------------
#' Library

library(terra)
library(data.table)
library(bit64)

#'------------------------------------------------------------------------------
#' Arguments
#' @param fsi
#' @param mask
#' @param root_path
#' @param base_name
#' @param out_path

#fsi <- fread("data/FIA_SPCD.csv")
#mask_plotID <- rast("data/mask_plotID.tif")
#root_path <- "E:/mapping_species/level3"
#base_name <- "2021-2021_001-365_HL_TSA_SEN2L_"
#out_path <- "E:/mapping_species/level4"

fsi <- fread("/panfs/jay/groups/17/cavender/shared/oak-wilt/misc/FIA/FIA_SPCD.csv")
mask_plotID <- rast("/panfs/jay/groups/17/cavender/shared/oak-wilt/misc/FIA/mask_plotID.tif")
root_path <- "/panfs/jay/groups/17/cavender/shared/oak-wilt/level3"
base_name <- "2021-2021_001-365_HL_TSA_SEN2L_"
out_path <- "/panfs/jay/groups/17/cavender/shared/oak-wilt/level4_rast"
options("scipen"=100, "digits"=6)

#'------------------------------------------------------------------------------
#' Function
rast_mask <- function(fsi, mask_plotID, root_path, base_name, out_path) {

  #Search for VI paths
  files <- list.files(path = root_path, 
                      pattern = ".tif", 
                      all.files = TRUE,
                      full.names = FALSE, 
                      recursive = TRUE)
  
  #Arrange path in frame
  frame <- data.table(matrix(unlist(strsplit(files, "/")), 
                             nrow= length(files), 
                             byrow=TRUE), stringsAsFactors=FALSE)
  colnames(frame) <- c("tile", "VI")
  tiles <- unique(frame$tile)
  tiles <- tiles[order(tiles)]
  
  #Clean fsi
  fsi <- subset(fsi, STATECD == 27 |
                  STATECD == 55)
  
  fsi <- subset(fsi, SUBP == 1)
  fsi$ID <- 1:nrow(fsi)
  
  #Make a vect
  fsi <- vect(fsi, geom=c("LON", "LAT"), crs= "+proj=longlat +datum=WGS84")
  fsi <- project(fsi, mask_plotID)
  fsi <- buffer(fsi, width = 100)
  
  #plot
  plot <- data.table(ID = fsi$ID, PLOTID = fsi$PLOTID)

  #Extraction function
  extraction <- function(i, tiles, base_name, mask_plotID, plot, root_path, out_path) {
    
    #Tile
    tile_name <- tiles[i]
    
    # Name of rast to read
    rast_name <- paste0(root_path, "/", tile_name, "/", base_name)
    
    #Read raster
    CCI <- terra::rast(paste0(rast_name, "CCI_TSI.tif"))
    CRE <- terra::rast(paste0(rast_name, "CRE_TSI.tif"))
    NDW <- terra::rast(paste0(rast_name, "NDW_TSI.tif"))
    KNV <- terra::rast(paste0(rast_name, "KNV_TSI.tif"))
    
    #Crop raster
    mask_crop <- crop(mask_plotID, ext(CCI))
    
    #Crop raster
    CCI$ID <- mask_crop
    CRE$ID <- mask_crop
    NDW$ID <- mask_crop
    KNV$ID <- mask_crop
    
    #Values
    CCI_values <- as.data.frame(CCI, 
                                xy = TRUE,
                                na.rm = FALSE)
    
    CRE_values <- as.data.frame(CRE, 
                                xy = TRUE,
                                na.rm = FALSE)
    
    NDW_values <- as.data.frame(NDW, 
                                xy = TRUE,
                                na.rm = FALSE)
    
    KNV_values <- as.data.frame(KNV, 
                                xy = TRUE,
                                na.rm = FALSE)
    
    #Remove NAN
    crop_nan <- is.nan(CCI_values$ID)
    
    if(any(crop_nan == FALSE) == TRUE) {
      
      CCI_values <- CCI_values[crop_nan != TRUE,]
      CRE_values <- CRE_values[crop_nan != TRUE,]
      NDW_values <- NDW_values[crop_nan != TRUE,]
      KNV_values <- KNV_values[crop_nan != TRUE,]
      
      #Melt frames
      CCI_melted <- melt(as.data.table(CCI_values), 
                         id.vars = c("ID", "x", "y"), 
                         variable.name = "date",
                         value.name = "CCI")
      
      CRE_melted <- melt(as.data.table(CRE_values), 
                         id.vars = c("ID", "x", "y"), 
                         variable.name = "date",
                         value.name = "CRE")
      
      NDW_melted <- melt(as.data.table(NDW_values), 
                         id.vars = c("ID", "x", "y"), 
                         variable.name = "date",
                         value.name = "NDW")
      
      KNV_melted <- melt(as.data.table(KNV_values), 
                         id.vars = c("ID", "x", "y"), 
                         variable.name = "date",
                         value.name = "KNV")
      
      #Merge
      values <- merge(CCI_melted, CRE_melted, by = c("ID", "x", "y", "date"))
      values <- merge(values, NDW_melted, by = c("ID", "x", "y", "date"))
      values <- merge(values, KNV_melted, by = c("ID", "x", "y", "date"))
      values <- merge(plot, values, by = "ID", all.x = FALSE, all.y = TRUE)
      
      #Modify date
      values$date <- as.Date(paste0(substr(values$date, 1, 4), "-", 
                                    substr(values$date, 5, 6), "-", 
                                    substr(values$date, 7, 8)))
      
      #Unique stand_key
      name_export <- paste0(out_path, "/", 
                            tile_name, ".txt")
      
      #Export
      fwrite(values, 
             name_export, 
             sep = "\t")
      
      #Remove residuals
      rm(list = c("tile_name", "rast_name", 
                  "CCI", "CRE", "NDW", "KNV", "mask_crop", "crop_nan",
                  "CCI_values", "CRE_values", "NDW_values", "KNV_values", 
                  "CCI_melted", "CRE_melted", "NDW_melted", "KNV_melted",
                  "values", "name_export"))
      
    }
     
  }
  
  # Loop for chucks in parallel  -----------------------------------------------
  for(k in 1:length(tiles)) {
    
    extraction(k, tiles = tiles, base_name = base_name, mask_plotID = mask_plotID, 
            plot = plot, root_path = root_path, out_path = out_path)
    
  }
}

rast_mask(fsi, mask_plotID, root_path, base_name, out_path)

