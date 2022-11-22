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

root_path <- "E:/DESIS"
#root_path <- "/media/antonio/Extreme SSD/DESIS"
fia <- fread(paste0(root_path, "/00_FIA/plots_WI_MN_20221102_reproject.csv"))
stand <- fread(paste0(root_path, "/03_reshape/DRYBIO_AG_reshape.csv"))
spectra <- fread(paste0(root_path, "/02_compile-files/02-spectra.csv"), header = TRUE)
out_path <- paste0(root_path, "/04_match_spectra-stand")

#'------------------------------------------------------------------------------
#' @example 

#'------------------------------------------------------------------------------
#' Script

match_datasets <- function(stand, spectra, out_path) {
  
  #------------------------------------
  #Clean stand
  
  #Remove stands with unknown trees
  stand[REMOVE1 > 0, REMOVE1 := NA]
  stand[REMOVE2 > 0, REMOVE1 := NA]
  stand[REMOVE3 > 0, REMOVE1 := NA]
  stand <- na.exclude(stand)

  #------------------------------------
  #Clean spectra
  
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
    
    sub_spectra <- subset(spectra, ID == unique_plots$PLT_CN[i])
    
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
  sub_stands <- stand[unique_plots$ID, ]
  sub_stands <- sub_stands[, c(2, 4:81)]
  
  species_number <- colSums(sub_stands[, 2:79])
  species_number <- species_number[species_number >= 1]
  species_number <- c("PLT_CN", names(species_number))
  sub_stands <- sub_stands[, ..species_number]
  
  #Spectra
  sub_spectra <- spectra[unique_plots$spectra, ]
  sub_spectra <- sub_spectra[, c(3, 14:248)]
  colnames(sub_spectra)[1] <- "PLT_CN"
  
  #Unique plots
  unique_plots <- unique_plots[, c(1, 2, 5)]
  unique_plots <- merge(unique_plots, fia[, c(2, 19, 20)], by = "PLT_CN", all.x = TRUE, all.y = FALSE)
  
  #-------------------------------------
  # Second step mergin (just in case)
  sub_stands <- merge(unique_plots, sub_stands, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
  sub_spectra <- merge(unique_plots, sub_spectra, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
  
  #Export
  fwrite(unique_plots, paste0(out_path, "/plots.csv"))
  fwrite(sub_stands[, c(1, 6:56)], paste0(out_path, "/species.csv"))
  fwrite(sub_spectra[, c(1, 4:5, 6:240)], paste0(out_path, "/spectra.csv"))

}