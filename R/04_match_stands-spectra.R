#'------------------------------------------------------------------------------
#' @title Match spectra and stands
#'------------------------------------------------------------------------------

#' @description Its match the spectral information with the reshape stands
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(bit64)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param stand
#' @param spectra
#' @param out_path

#'------------------------------------------------------------------------------
#' Script

# 10123
# 13089

match_datasets <- function(fia, 
                           stand,
                           master_spectra,
                           spectra_sw, 
                           out_path) {
  
  #------------------------------------
  # Clean spectra
  
  # Use samples that were checked
  master_spectra$scene <- substr(master_spectra$scene, 1, 46)
  master_spectra$date <- as.IDate(paste0(substr(master_spectra$scene, 32, 35), "-",
                                         substr(master_spectra$scene, 36, 37), "-",
                                         substr(master_spectra$scene, 38, 39)))
  spectra_sw$scene <- substr(spectra_sw$scene, 1, 46)
  
  spectra <- merge(master_spectra[, c("folder", "scene", "PLT_CN", "date")], spectra_sw, 
                   by = c("folder", "scene", "PLT_CN"),
                   all.x = TRUE, all.y = FALSE)
  spectra$sample <- 1:nrow(spectra)
  
  # Error of SW
  #remove <- spectra$`680.91`
  remove <- spectra$`692.53`
  spectra <- spectra[is.na(remove) == FALSE, ]
  spectra$sample <- 1:nrow(spectra)
  
  #------------------------------------
  # Master file for merge
  unique_plots <- stand[, .N, by = c("Domain", "PLT_CN", "date")]
  unique_plots$sample <- 0
  unique_plots$plot <- 0
  unique_plots$dtime <- 0
  
  # Progress bar
  pb <- txtProgressBar(min = 1, max = nrow(unique_plots), style = 3)
  
  # Search for spectra
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
    
    # If there is one observation
    } else {
      
      dtime <- sub_spectra$date - unique_plots$date[i]
      
      if(nrow(sub_spectra) == 1) {
        
        unique_plots$sample[i] <- sub_spectra$sample[1]
        unique_plots$plot[i] <- i
        unique_plots$dtime[i] <- dtime[1]
      
      } else {
        
        unique_plots$sample[i] <- sub_spectra$sample[which.min(dtime)]
        unique_plots$plot[i] <- i
        unique_plots$dtime[i] <- dtime[which.min(dtime)]
        
      }
    }
  }
  
  #Plots with usefull spectra
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
  
  #Spectra
  sub_spectra <- spectra[unique_plots$sample, ]
  sub_spectra <- sub_spectra[, .SD, .SDcols=c(c(3, 15:(ncol(sub_spectra)-1)))]
  
  all.equal(sub_spectra$PLT_CN, sub_stands$PLT_CN)
  
  #Export
  fwrite(sub_stands, paste0(out_path, "/stands_SW.csv"))
  fwrite(sub_spectra, paste0(out_path, "/spectra_SW.csv"))

}

#'------------------------------------------------------------------------------
#' @example 

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
fia <- fread(paste0(root_path, "/00_FIA/plots_eastern_20230427.csv"))
stand <- fread(paste0(root_path, "/03_reshape/AGB_TPA-reshape.csv"))
master_spectra <- fread(paste0(root_path, "/02_compile-files/master_spectra.csv"), header = TRUE)
spectra_sw <- fread(paste0(root_path, "/02_compile-files/02_spectra-SW.txt"), header = TRUE)
out_path <- paste0(root_path, "/04_match_spectra-stand")

match_datasets(fia, spectra, stand, out_path)
