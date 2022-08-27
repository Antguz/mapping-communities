#'------------------------------------------------------------------------------
#' @title DNR Forest Stand Inventory
#'------------------------------------------------------------------------------

#' @description Select the polygons to use from the DNR Forest Stand Inventory
#' @return A .gpkg with the selected polygons

#'------------------------------------------------------------------------------
#' Library

library(terra)
library(data.table)

#'------------------------------------------------------------------------------
#' Arguments
#' @param fsi
#' @param out_path

#'------------------------------------------------------------------------------
#' @example 

fsi <- vect("data/biota_dnr_forest_stand_inventory.gpkg")
out_path <- "data/fsi_clean.gpkg"

#'------------------------------------------------------------------------------
#' Function

fsi_clean <- function(fsi, out_path) {
  
  #Subset by year
  fsi_clean <- subset(fsi, fsi$SURVEY_YR >= 2017 & fsi$SURVEY_YR <= 2021)
  
  #Subset by stand age
  fsi_clean <- subset(fsi_clean, fsi_clean$STAND_AGE >= 5)
  
  #Export
  writeVector(fsi_clean, out_path, filetype = "GPKG")
  
}
