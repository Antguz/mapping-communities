#'------------------------------------------------------------------------------
#' @title Deep clean and merging of the GIS files
#'------------------------------------------------------------------------------

#' @description It merge the functional diversity stack with the GIS metadata
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(terra)
library(data.table)

#'------------------------------------------------------------------------------
#' Arguments
#' @param gis
#' @param fd

#'------------------------------------------------------------------------------
#' @example 

fsi <- vect("data/fsi_clean.gpkg")
fd <- fread("E:/mapping_species/FD.txt")

#'------------------------------------------------------------------------------
#' Function

fd_meta <- function(fsi, fd) {
  
  #Deep cleanning of the FSI file ----------------------------------------------
  #Remove by Main Cover Type
  fsi <- subset(fsi, fsi$MN_CTYPE != 0 & #None
                     fsi$MN_CTYPE != 82 & #Cutover Area
                     fsi$MN_CTYPE != 83 & #Lowland Grass
                     fsi$MN_CTYPE != 84 & #Upland Grass
                     fsi$MN_CTYPE != 85 & #Lowland Brush
                     fsi$MN_CTYPE != 86 & #Duffs
                     fsi$MN_CTYPE != 88 & #Moss
                     fsi$MN_CTYPE != 90 & #unknown
                     fsi$MN_CTYPE != 91 & #Agriculture
                     fsi$MN_CTYPE != 92 & #Industrial Development
                     fsi$MN_CTYPE != 93 & #Recreation Develop
                     fsi$MN_CTYPE != 94 & #Roads
                     fsi$MN_CTYPE != 95 & #Rock Outcrop
                     fsi$MN_CTYPE != 96 & #Permanent Water
                     fsi$MN_CTYPE != 97 & #Nonpermanent Water
                     fsi$MN_CTYPE != 98 & #Marsh
                     fsi$MN_CTYPE != 99) #Muskeg
  
  #Remove by Main Cover Type Size Class
  fsi <- subset(fsi, fsi$MN_SIZE != 0) #N/A
                  
  #Remove by Timber Status
  fsi <- subset(fsi, fsi$TMBR_STAT != 0) #N/A
  
  #Remove by Level of Reconnaissance
  fsi <- subset(fsi, fsi$RECON != 0) #None
  
  #Remove by Main Species
  fsi <- subset(fsi, fsi$MN_SPP != 0 & #None
                     fsi$MN_SPP != 99) #Miscellaneous
  
  #Remove by Unit of Measure for the volume entered for Main Species
  fsi <- subset(fsi, fsi$MN_SPP_UOM != "N/A")
  
  #Remove by Percentage of Main Species Affected by Damage
  fsi <- subset(fsi, fsi$PCT_AFFECT <= 2) #<= 11-25%
  
  #Remove by Percentage of Main Species with Mortality
  fsi <- subset(fsi, fsi$PCT_MORTAL <= 2) #<= 11-25%
  
  #Remove by Age of Stand
  fsi <- subset(fsi, fsi$STAND_AGE >= 5 &
                     fsi$STAND_AGE <= 250) #<= 11-25%
  
  #Remove if more than 10 species are present
  fsi <- subset(fsi, fsi$MISC_CORDS == 0 &
                  fsi$MISC_MBF == 0) #<= 11-25%
  
  #Remove if more than 10 species are present
  fsi <- subset(fsi, fsi$MISC_CORDS == 0 |
                  fsi$MISC_MBF == 0) #<= 11-25%
  
  #Remove if Miscellaneous species are present
  fsi <- subset(fsi, fsi$MN_SPP != 0 & #0 =	No trees
                     fsi$MN_SPP != 99 & #99 =	Miscellaneous
                     fsi$SPP2 != 99 &
                     fsi$SPP3 != 99 &
                     fsi$SPP4 != 99 &
                     fsi$SPP5 != 99 &
                     fsi$SPP6 != 99 &
                     fsi$SPP7 != 99 &
                     fsi$SPP8 != 99 &
                     fsi$SPP9 != 99 &
                     fsi$SPP10 != 99)
  
  #Remove by area
  fsi$SHAPE_Area <- expanse(fsi)
  fsi <- subset(fsi, fsi$SHAPE_Area <= 500000) #50 hec
  
  #Export
  writeVector(fsi, "fsi_final.gpkg", filetype = "GPKG")
  
  #Reframe for data analysis
  frame_fsi <- as.data.table(fsi)
  
  sp1 <- frame_fsi[, c("STAND_KEY", "MN_SPP", "MN_SPP_DBH", "MN_SPP_VOL", "MN_SPP_UOM", "MN_SPP_DAM")]
  sp1$SPP1_DST <- NA
  sp1 <- sp1[, c("STAND_KEY", "MN_SPP", "MN_SPP_DBH", "MN_SPP_VOL", "MN_SPP_UOM", "SPP1_DST", "MN_SPP_DAM")]
  sp2 <- frame_fsi[, c("STAND_KEY", "SPP2", "SPP2_DBH", "SPP2_VOL", "SPP2_UOM", "SPP2_DST", "SPP2_DAM")]
  sp3 <- frame_fsi[, c("STAND_KEY", "SPP3", "SPP3_DBH", "SPP3_VOL", "SPP3_UOM", "SPP3_DST", "SPP3_DAM")]
  sp4 <- frame_fsi[, c("STAND_KEY", "SPP4", "SPP4_DBH", "SPP4_VOL", "SPP4_UOM", "SPP4_DST", "SPP4_DAM")]
  sp5 <- frame_fsi[, c("STAND_KEY", "SPP5", "SPP5_DBH", "SPP5_VOL", "SPP5_UOM", "SPP5_DST", "SPP5_DAM")]
  sp6 <- frame_fsi[, c("STAND_KEY", "SPP6", "SPP6_DBH", "SPP6_VOL", "SPP6_UOM", "SPP6_DST", "SPP6_DAM")]
  sp7 <- frame_fsi[, c("STAND_KEY", "SPP7", "SPP7_DBH", "SPP7_VOL", "SPP7_UOM", "SPP7_DST", "SPP7_DAM")]
  sp8 <- frame_fsi[, c("STAND_KEY", "SPP8", "SPP8_DBH", "SPP8_VOL", "SPP8_UOM", "SPP8_DST", "SPP8_DAM")]
  sp9 <- frame_fsi[, c("STAND_KEY", "SPP9", "SPP9_DBH", "SPP9_VOL", "SPP9_UOM", "SPP9_DST", "SPP9_DAM")]
  sp10 <- frame_fsi[, c("STAND_KEY", "SPP10", "SPP10_DBH", "SPP10_VOL", "SPP10_UOM", "SPP10_DST", "SPP10_DAM")]
  
  colnames(sp1) <- c("STAND_KEY", "SP", "DBH", "VOL", "UOM", "DST", "DAM")
  colnames(sp2) <- c("STAND_KEY", "SP", "DBH", "VOL", "UOM", "DST", "DAM")
  colnames(sp3) <- c("STAND_KEY", "SP", "DBH", "VOL", "UOM", "DST", "DAM")
  colnames(sp4) <- c("STAND_KEY", "SP", "DBH", "VOL", "UOM", "DST", "DAM")
  colnames(sp5) <- c("STAND_KEY", "SP", "DBH", "VOL", "UOM", "DST", "DAM")
  colnames(sp6) <- c("STAND_KEY", "SP", "DBH", "VOL", "UOM", "DST", "DAM")
  colnames(sp7) <- c("STAND_KEY", "SP", "DBH", "VOL", "UOM", "DST", "DAM")
  colnames(sp8) <- c("STAND_KEY", "SP", "DBH", "VOL", "UOM", "DST", "DAM")
  colnames(sp9) <- c("STAND_KEY", "SP", "DBH", "VOL", "UOM", "DST", "DAM")
  colnames(sp10) <- c("STAND_KEY", "SP", "DBH", "VOL", "UOM", "DST", "DAM")
  
  sp1$SPP_RANK <- 1
  sp2$SPP_RANK <- 2
  sp3$SPP_RANK <- 3
  sp4$SPP_RANK <- 4
  sp5$SPP_RANK <- 5
  sp6$SPP_RANK <- 6
  sp7$SPP_RANK <- 7
  sp8$SPP_RANK <- 8
  sp9$SPP_RANK <- 9
  sp10$SPP_RANK <- 10
  
  SPP <- rbind(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9, sp10)
  SPP <- subset(SPP, SP != 0)
  SPP <- SPP[, c("STAND_KEY", "SPP_RANK", "SP", "DBH", "VOL", "UOM", "DST", "DAM")]
  
  to_merge <- frame_fsi[, c("STAND_KEY", "ADMIN1", "MN_CTYPE", "MN_SIZE", "MN_DENS", 
                            "UTYPE", "USIZE", "UDENSE", "SURVEY_DT", "SURVEY_YR", 
                            "STAND_AGE", "TOTAL_ACRE", "TMBR_STAT", "DEV_STAT", "TOPO",
                            "PHYS", "COND", "ORIGIN", "RECON", "BASAL_AREA", "TOTAL_CORD",
                            "TOTAL_MBF", "SITE_INDEX", "SI_SPECIES", "MGMT1", "OG1_ID", 
                            "SHAPE_Area", "CUR_AGE")]
  final_frame <- merge(to_merge, SPP, by = "STAND_KEY", all.x = TRUE, all.y = TRUE)
  final_frame[UOM == "1", UOM := "Cords"]
  final_frame[UOM == "2", UOM := "MBF"]
  final_frame[UOM == "3", UOM := "Stems"]
  
  final_frame <- final_frame[order(STAND_KEY, SPP_RANK)]
  
  #N species
  final_frame[, N_SPP := max(SPP_RANK), by = "STAND_KEY"]
  
  #Functional diversity by volumen
  final_frame[, fd_fric(VOL), by = "STAND_KEY"]
  
  fwrite(final_frame, "final_frame.txt", sep = "\t")
  
  #Merge all
  frame_master <- merge(final_frame, fd, by = "STAND_KEY", all.x = TRUE, allow.cartesian=TRUE)
  na_date <- is.na(frame_master$date)
  frame_master <- frame_master[na_date == FALSE]
  
  fwrite(frame_master, "frame_master.txt", sep = "\t")
  
}
