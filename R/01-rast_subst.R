#'------------------------------------------------------------------------------
#' @title Raster subtraction
#'------------------------------------------------------------------------------

#' @description Subtract pixels from Geotiffs in titles index and a vector file
#' @return .txt files with names according the vector file

#'------------------------------------------------------------------------------
#' Library

library(terra)
library(data.table)
library(foreach)
#library(doParallel)

#'------------------------------------------------------------------------------
#' Arguments
#' @param title
#' @param fsi
#' @param root_path
#' @param base_name
#' @param out_path

#'------------------------------------------------------------------------------
#' @example 

title <- vect("data/title_clipped.gpkg")
fsi <- vect("data/fsi_clean.gpkg")
root_path <- "/media/antonio/antonio_ssd/mapping_species/level3_new"
base_name <- "2021-2021_001-365_HL_TSA_SEN2L_"
out_path <- "/media/antonio/antonio_ssd/mapping_species/level3_crop"

#'------------------------------------------------------------------------------
#' Function
rast_subst <- function(title, fsi, root_path, base_name, out_path, threads = 22) {
  
  #Check projection
  fsi <- project(fsi, title)
  fsi <- buffer(fsi, width=0)
  
  #See interception of layers
  title_fsi <- intersect(fsi, title)
  title_fsi <- subset(title_fsi, title_fsi$Tile_ID == "X0015_Y0019")
  
  #As frame
  fsi_frame <- as.data.table(title_fsi)
  
  #Unique stand_key
  tile <- unique(fsi_frame$Tile_ID)
  
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
                          
                          #Unique forest stand
                          tile_vector <- terra::subset(title_fsi, fsi_frame$Tile_ID == tile_name)
                          stands <- unique(tile_vector$STAND_KEY) 
                          ids <- data.table(ID = 1:length(stands),
                                            STAND_KEY = stands)
                          
                          #Extract values
                          CCI_values <- extract(CCI, 
                                                tile_vector, 
                                                method = "simple",
                                                factors = TRUE,
                                                cells = FALSE,
                                                xy = TRUE,
                                                exact = TRUE,
                                                weights = TRUE, 
                                                touches = TRUE)
                          
                          CRE_values <- extract(CRE, 
                                                tile_vector, 
                                                method = "simple",
                                                factors = TRUE,
                                                cells = FALSE,
                                                xy = TRUE,
                                                exact = TRUE,
                                                weights = TRUE, 
                                                touches = TRUE)
                          
                          NDW_values <- extract(NDW, 
                                                tile_vector, 
                                                method = "simple",
                                                factors = TRUE,
                                                cells = FALSE,
                                                xy = TRUE,
                                                exact = TRUE,
                                                weights = TRUE, 
                                                touches = TRUE)
                          
                          KNV_values <- extract(KNV, 
                                                tile_vector, 
                                                method = "simple",
                                                factors = TRUE,
                                                cells = FALSE,
                                                xy = TRUE,
                                                exact = TRUE,
                                                weights = TRUE, 
                                                touches = TRUE)
                          
                          #Melt frames
                          CCI_melted <- melt(as.data.table(CCI_values), 
                                             id.vars = c("ID", "weight", "x", "y"), 
                                             variable.name = "date",
                                             value.name = "value")
                          
                          CRE_melted <- melt(as.data.table(CRE_values), 
                                             id.vars = c("ID", "weight", "x", "y"), 
                                             variable.name = "date",
                                             value.name = "value")
                          
                          NDW_melted <- melt(as.data.table(NDW_values), 
                                             id.vars = c("ID", "weight", "x", "y"), 
                                             variable.name = "date",
                                             value.name = "value")
                          
                          KNV_melted <- melt(as.data.table(KNV_values), 
                                             id.vars = c("ID", "weight", "x", "y"), 
                                             variable.name = "date",
                                             value.name = "value")
                          
                          #Merge
                          values <- merge(CCI_melted, CRE_melted, by = c("ID", "weight", "x", "y", "date"))
                          colnames(values) <- c("ID", "weight", "x", "y", "date", "CCI", "CRE")
                          values <- merge(values, NDW_melted, by = c("ID", "weight", "x", "y", "date"))
                          colnames(values) <- c("ID", "weight", "x", "y", "date", "CCI", "CRE", "NDW")
                          values <- merge(values, KNV_melted, by = c("ID", "weight", "x", "y", "date"))
                          colnames(values) <- c("ID", "x", "y", "weight", "date", "CCI", "CRE", "NDW", "KNV")
                          
                          #Merge ID stand
                          values <- merge(ids, values, by = "ID", all.x = TRUE, all.y = TRUE)
                          
                          #Weight subset
                          values <- subset(values, weight >= 1.0)
                          values <- values[, c(2:4, 6:10)]
                          values$date <- as.Date(paste0(substr(values$date, 1, 4), "-", 
                                                        substr(values$date, 5, 6), "-", 
                                                        substr(values$date, 7, 8)))
                          
                          for(ii in 1:length(stands)) {
                            
                            sub_value <- subset(values, STAND_KEY == stands[ii])
                            
                            #Name to export
                            name_export <- paste0(out_path, "/", 
                                                  stands[ii], "_",
                                                  tile_name,".txt")
                            
                            #Export
                            fwrite(sub_value, 
                                   name_export, 
                                   sep = "\t")
                            
                          }
                        }
  
  #Stop cluster
  #stopCluster(cl)
  #gc()
  
}

rast_subst(title, fsi, root_path, base_name, out_path, threads = 22)
