#'------------------------------------------------------------------------------
#' @tile Raster subtraction
#'------------------------------------------------------------------------------

#' @description Subtract pixels from Geotiffs in tiles index and a vector file
#' @return .txt files with names according the vector file

#'------------------------------------------------------------------------------
#' Library

library(terra)
library(data.table)
library(foreach)
library(doParallel)

#'------------------------------------------------------------------------------
#' Arguments
#' @param tile
#' @param fsi
#' @param root_path
#' @param base_name
#' @param out_path

#'------------------------------------------------------------------------------
#' @example 

tile <- vect("/panfs/jay/groups/17/cavender/shared/oak-wilt/misc/GIS/title_clipped.shp")
fsi <- fread("/panfs/jay/groups/17/cavender/shared/oak-wilt/misc/FIA/FIA_SPCD.csv")
root_path <- "/panfs/jay/groups/17/cavender/shared/oak-wilt/level3"
base_name <- "2021-2021_001-365_HL_TSA_SEN2L_"
out_path <- "/panfs/jay/groups/17/cavender/shared/oak-wilt/level4"

#tile <- vect("data/title_clipped.shp")
#fsi <- vect("data/fsi_clean.shp")
#root_path <- "/media/antonio/antonio_ssd/mapping_species/level3"
#base_name <- "2021-2021_001-365_HL_TSA_SEN2L_"
#out_path <- "/media/antonio/antonio_ssd/mapping_species/level4"
#threads <- 24

#tile <- vect("data/title_clipped.gpkg")
#fsi <- fread("data/FIA_SPCD.csv")
#root_path <- "E:/mapping_species/level3"
#base_name <- "2021-2021_001-365_HL_TSA_SEN2L_"
#out_path <- "E:/mapping_species/level4"
#threads <- 24

#'------------------------------------------------------------------------------
#' Function
rast_subst <- function(tile, fsi, root_path, base_name, out_path, threads = 24) {
  
  #Clean fsi
  fsi <- subset(fsi, STATECD == 27 |
                     STATECD == 55)
  
  fsi <- subset(fsi, SUBP == 1)
  
  #Make a vect
  fsi <- vect(fsi, geom=c("LON", "LAT"), crs= "+proj=longlat +datum=WGS84")
  fsi <- project(fsi, tile)
  fsi <- buffer(fsi, width = 100)
  tile_fsi <- intersect(fsi, tile)
  
  #Order tile
  tile <- unique(tile_fsi$Tile_ID)
  tile <- tile[order(tile)]
  ids <- data.table(ID = 1:dim(tile_fsi)[1],
                    PLOTID = tile_fsi$PLOTID)
  
  # Set up cluster
  #cl <- makeCluster(threads, type = "FORK")
  #registerDoParallel(cl)
  
  # Loop over tiles ------------------------------------------------------------
  extraction <- foreach(i = 1:length(tile),
                        .packages = c("terra", "data.table"),
                        .inorder = F) %do% {
                          
                          #Tile
                          tile_name <- tile[i]
                          
                          # Name of rast to read
                          rast_name <- paste0(root_path, "/", tile_name, "/", base_name)
                          
                          #Read raster
                          CCI <- terra::rast(paste0(rast_name, "CCI_TSI.tif"))
                          CRE <- terra::rast(paste0(rast_name, "CRE_TSI.tif"))
                          NDW <- terra::rast(paste0(rast_name, "NDW_TSI.tif"))
                          KNV <- terra::rast(paste0(rast_name, "KNV_TSI.tif"))
                          
                          #Extract values
                          CCI_values <- terra::extract(CCI, 
                                                       tile_fsi, 
                                                       method = "simple",
                                                       factors = FALSE,
                                                       cells = TRUE,
                                                       xy = TRUE,
                                                       ID = TRUE,
                                                       exact = FALSE,
                                                       weights = TRUE, 
                                                       touches = TRUE)
                          
                          CRE_values <- terra::extract(CRE, 
                                                       tile_fsi, 
                                                       method = "simple",
                                                       factors = FALSE,
                                                       cells = TRUE,
                                                       xy = TRUE,
                                                       ID = TRUE,
                                                       exact = FALSE,
                                                       weights = TRUE, 
                                                       touches = TRUE)
                          
                          NDW_values <- terra::extract(NDW, 
                                                       tile_fsi, 
                                                       method = "simple",
                                                       factors = FALSE,
                                                       cells = TRUE,
                                                       xy = TRUE,
                                                       ID = TRUE,
                                                       exact = FALSE,
                                                       weights = TRUE, 
                                                       touches = TRUE)
                          
                          KNV_values <- terra::extract(KNV, 
                                                       tile_fsi, 
                                                       method = "simple",
                                                       factors = FALSE,
                                                       cells = TRUE,
                                                       xy = TRUE,
                                                       ID = TRUE,
                                                       exact = FALSE,
                                                       weights = TRUE, 
                                                       touches = TRUE)
                          
                          #Melt frames
                          CCI_melted <- melt(as.data.table(CCI_values), 
                                             id.vars = c("ID", "cell", "weight", "x", "y"), 
                                             variable.name = "date",
                                             value.name = "value")
                          
                          CRE_melted <- melt(as.data.table(CRE_values), 
                                             id.vars = c("ID", "cell", "weight", "x", "y"), 
                                             variable.name = "date",
                                             value.name = "value")
                          
                          NDW_melted <- melt(as.data.table(NDW_values), 
                                             id.vars = c("ID", "cell", "weight", "x", "y"), 
                                             variable.name = "date",
                                             value.name = "value")
                          
                          KNV_melted <- melt(as.data.table(KNV_values), 
                                             id.vars = c("ID", "cell", "weight", "x", "y"), 
                                             variable.name = "date",
                                             value.name = "value")
                          
                          #Merge
                          values <- merge(CCI_melted, CRE_melted, by = c("ID", "cell", "weight", "x", "y", "date"))
                          colnames(values) <- c("ID", "cell", "weight", "x", "y", "date", "CCI", "CRE")
                          values <- merge(values, NDW_melted, by = c("ID", "cell", "weight", "x", "y", "date"))
                          colnames(values) <- c("ID", "cell", "weight", "x", "y", "date", "CCI", "CRE", "NDW")
                          values <- merge(values, KNV_melted, by = c("ID", "cell", "weight", "x", "y", "date"))
                          colnames(values) <- c("ID", "cell", "x", "y", "weight", "date", "CCI", "CRE", "NDW", "KNV")
                          
                          #Merge ID stand
                          values <- merge(ids, values, by = "ID", all.x = FALSE, all.y = TRUE)
                          
                          #Weight subset
                          values <- subset(values, weight >= 0.5)
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
                                      "CCI", "CRE", "NDW", "KNV",
                                      "CCI_values", "CRE_values", "NDW_values", "KNV_values", 
                                      "CCI_melted", "CRE_melted", "NDW_melted", "KNV_melted",
                                      "values", "name_export"))
                          
                        }
  
  #Stop cluster
  #stopCluster(cl)
  #gc()
  
}

rast_subst(tile, fsi, root_path, base_name, out_path, threads = 100)
