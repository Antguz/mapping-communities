#'------------------------------------------------------------------------------
#' @title Match spectra and stands
#'------------------------------------------------------------------------------

#' @description Its match the spectral information with the reshape stands
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param stand
#' @param spectra
#' @param out_path

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
fia <- fread(paste0(root_path, "/00_FIA/plots_WI_MN_20221102.csv"))
stand <- fread(paste0(root_path, "/03_reshape/BA_TPA-reshape.csv"))
spectra_assesment <- fread(paste0(root_path, "/02_compile-files/master_spectra.csv"), header = TRUE)
spectra <- fread(paste0(root_path, "/02_compile-files/02-spectra_SW.txt"), header = TRUE)
out_path <- paste0(root_path, "/04_match_spectra-stand")

#'------------------------------------------------------------------------------
#' @example 

match_datasets(fia, spectra, stand, out_path)

#'------------------------------------------------------------------------------
#' Script

match_datasets <- function(fia, spectra, stand, out_path) {
  
  #------------------------------------
  #Clean stand
  
  #Remove stands with unknown trees
  stand[REMOVE1 > 0, REMOVE1 := NA]
  stand[REMOVE2 > 0, REMOVE1 := NA]
  stand[REMOVE3 > 0, REMOVE1 := NA]
  stand <- na.exclude(stand)

  #------------------------------------
  #Clean spectra
  
  #Use samples that were checked
  spectra_assesment$scene <- substr(spectra_assesment$scene, 1, 46)
  spectra$scene <- substr(spectra$scene, 1, 46)
  
  spectra <- merge(spectra_assesment, spectra, 
                   by = c("folder", "scene", "PLT_CN"),
                   all.x = TRUE, all.y = FALSE)
  
  #Remove flaged pixels
  spectra <- spectra[Shadow != 1]
  spectra <- spectra[`Haze (land)` != 1]
  spectra <- spectra[`Cloud (land)` != 1]
  spectra <- spectra[Water != 1]
  spectra$sample <- 1:nrow(spectra)
  
  #------------------------------------
  #Master file for merge
  unique_plots <- stand[, .N, by = c("PLT_CN", "date")]
  colnames(unique_plots)[3] <- "ID"
  unique_plots$ID <- 1:nrow(unique_plots)
  unique_plots$spectra <- 0
  unique_plots$dtime <- 0
  
  #Search for spectra
  for(i in 1:nrow(unique_plots)) {
    
    sub_spectra <- subset(spectra, PLT_CN == unique_plots$PLT_CN[i])
    
    #If no observations
    if(nrow(sub_spectra) == 0) {
      
      unique_plots$spectra[i] <- NA
      unique_plots$dtime[i] <- NA
    
    #If there is one observation
    } else {
      
      date <- as.IDate(paste0(substr(sub_spectra$scene, 32, 35), "-",
                              substr(sub_spectra$scene, 36, 37), "-",
                              substr(sub_spectra$scene, 38, 39)))
      
      dtime <- date - unique_plots$date[i]
      
      if(nrow(sub_spectra) == 1) {
        
        unique_plots$spectra[i] <- sub_spectra$sample[1]
        unique_plots$dtime[i] <- dtime[1]
      
      } else {
        
        unique_plots$spectra[i] <- sub_spectra$sample[which.min(dtime)]
        unique_plots$dtime[i] <- dtime[which.min(dtime)]
        
      }
    }
  }
  
  #Plots with usefull spectra
  unique_plots <- na.exclude(unique_plots)
  
  #------------------------------------
  # Matching files
  
  #Stands
  #sub_stands <- stand[unique_plots$ID, ]
  sub_stands <- stand[, .SD, .SDcols=c(2, 4:ncol(stand)),]
  species_number <- colSums(sub_stands[, .SD, .SDcols=c(2:ncol(sub_stands)),])
  species_number <- species_number[species_number >= 1]
  species_number <- c("PLT_CN", names(species_number))
  sub_stands <- sub_stands[, ..species_number]
  
  #Spectra
  sub_spectra <- spectra[unique_plots$spectra, ]
  sub_spectra <- sub_spectra[, .SD, .SDcols=c(c(3, 14:(ncol(sub_spectra)-1)))]
  
  #Unique plots
  unique_plots <- unique_plots[, c(1, 2, 5)]
  unique_plots <- merge(unique_plots, fia[, c(1, 13, 14, 15)], by = "PLT_CN", all.x = TRUE, all.y = FALSE)
  
  #-------------------------------------
  # Second step mergin (just in case)
  sub_stands <- merge(unique_plots, sub_stands, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
  sub_spectra <- merge(unique_plots, sub_spectra, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
  
  #Export
  fwrite(unique_plots, paste0(out_path, "/plots_SW.csv"))
  fwrite(sub_stands, paste0(out_path, "/species_SW.csv"))
  fwrite(sub_spectra, paste0(out_path, "/spectra_SW.csv"))

}
