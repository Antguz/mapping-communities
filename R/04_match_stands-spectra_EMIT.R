#'------------------------------------------------------------------------------
#' @title Match spectra from EMIT and stands
#'------------------------------------------------------------------------------

#' @description Its match the spectral information from EMIT with the reshape 
#' stands
#' 
#' @return Frame files with best batch between datasets

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(bit64)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param stand
#' @param spectra_sw
#' @param master_spectra
#' @param sensor
#' @param out_path

#'------------------------------------------------------------------------------
#' @example 

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
fia <- fread(paste0(root_path, "/00_FIA/plots_eastern_20230427.csv"))
stand <- fread(paste0(root_path, "/03_reshape/AGB_TPA-reshape.csv"))
master_spectra <- fread(paste0(root_path, "/02_compile-files/EMIT_master_spectra.csv"), header = TRUE)
spectra_sw <- fread(paste0(root_path, "/02_compile-files/02_EMIT_spectra-SW.txt"), header = TRUE)
out_path <- paste0(root_path, "/04_match_spectra-stand")

match_datasets_EMIT(fia, stand, spectra_sw, master_spectra, out_path)

#'------------------------------------------------------------------------------
#' Script
match_datasets_EMIT <- function(fia, 
                                stand,
                                spectra_sw,
                                master_spectra,
                                out_path) {
  
  #------------------------------------
  # Clean spectra
  # Use samples that were checked
  master_spectra$scene <- substr(master_spectra$scene, 1, 44)
  master_spectra$date <- as.IDate(paste0(substr(master_spectra$scene, 18, 21), "-",
                                         substr(master_spectra$scene, 22, 23), "-",
                                         substr(master_spectra$scene, 24, 25)))
  spectra_sw$scene <- substr(spectra_sw$scene, 1, 44)
  spectra <- merge(master_spectra[, c("scene", "PLT_CN", "date")], spectra_sw, 
                   by = c("scene", "PLT_CN"),
                   all.x = TRUE, all.y = FALSE)
  spectra$sample <- 1:nrow(spectra)
  
  # Error of SW
  remove <- spectra$`678.205584716796`
  spectra <- spectra[is.na(remove) == FALSE, ]
  spectra$sample <- 1:nrow(spectra)
  
  #------------------------------------
  # Master file for merge
  unique_plots <- stand[, .N, by = c("Domain", "PLT_CN", "date")]
  unique_plots$sample <- 0
  unique_plots$plot <- 0
  unique_plots$dtime <- 0
  unique_plots$nobservations <- 0
  
  # Progress bar
  pb <- txtProgressBar(min = 1, max = nrow(unique_plots), style = 3)
  
  # Search for spectra #3970
  for(i in 1:nrow(unique_plots)) {
    
    # update progress bar
    setTxtProgressBar(pb, i)
    
    # Subset plot
    sub_spectra <- subset(spectra, PLT_CN == unique_plots$PLT_CN[i])
    
    # If no observations
    if(nrow(sub_spectra) == 0) {
      
      unique_plots$sample[i] <- NA
      unique_plots$plot[i] <- NA
      unique_plots$dtime[i] <- NA
      
      # If there is one or mode observation
    } else {
      
      # If there is one observation
      if(nrow(sub_spectra) == 1) {
        
        unique_plots$sample[i] <- sub_spectra$sample[1]
        unique_plots$plot[i] <- i
        unique_plots$dtime[i] <- sub_spectra$date[1] - unique_plots$date[i]
        unique_plots$nobservations[i] <- 1
        
        # If there is two or more observation
      } else {
        
        unique_plots$nobservations[i] <- nrow(sub_spectra)
        
        # Select observations closest to the inventory year
        delta_year <- abs(year(sub_spectra$date) - year(unique_plots$date[i]))
        sub_spectra <- sub_spectra[delta_year == delta_year[which.min(delta_year)]]
        
        # Select observations closest to the peak of season
        peak <- as.IDate(paste0(unique(year(sub_spectra$date)), "-07-15"))
        delta_time <- sub_spectra$date-peak
        
        unique_plots$sample[i] <- sub_spectra$sample[which.min(abs(delta_time))]
        unique_plots$plot[i] <- i
        unique_plots$dtime[i] <- delta_time[which.min(delta_time)]
        
      }
    }
  }
  
  #Plots with useful spectra
  unique_plots <- na.exclude(unique_plots)
  
  #------------------------------------
  # Matching files
  
  #Stands
  sub_stands <- stand[unique_plots$plot, ]
  
  nspecies <- apply(sub_stands[, .SD, .SDcols=c(4:ncol(sub_stands)),], 1, function(value) {
    return(length(value[value > 0]))
  })
  
  biomass <- apply(sub_stands[, .SD, .SDcols=c(4:ncol(sub_stands)),], 1, sum)
  
  sub_stands <- sub_stands[, .SD, .SDcols=c(2, 1, 3)]
  sub_stands$nspecies <- nspecies
  sub_stands$biomass <- biomass
  sub_stands$dtime <- unique_plots$dtime
  sub_stands$nobservations <- unique_plots$nobservations
  
  #Spectra
  sub_spectra <- spectra[unique_plots$sample, ]
  sub_spectra <- sub_spectra[, .SD, .SDcols=c(c(2, 16:(ncol(sub_spectra)-1)))]
  
  all.equal(sub_spectra$PLT_CN, sub_stands$PLT_CN)
  
  #Export
  fwrite(sub_stands, paste0(out_path, "/EMIT_stands_SW.csv"))
  fwrite(sub_spectra, paste0(out_path, "/EMIT_spectra_SW.csv"))
  
}
