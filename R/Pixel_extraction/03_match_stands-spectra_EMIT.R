#'------------------------------------------------------------------------------
#' @title Match spectra from EMIT and stands
#'------------------------------------------------------------------------------

#' @description A function to match the spectral information from EMIT with the reshape 
#' stands
#' 
#' @return Frame files with best batch between datasets

#'------------------------------------------------------------------------------
#' @Library

library(data.table)
library(bit64)

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param fia A data.table with metadata about the inventories
#' @param master_comm A data.table with the curated inventories from FAI
#' @param master_spectra A data.table with the curated spectra extracted from pixels
#' @param spectra_sg A data.table with the savitzky-golay spectra
#' @param spectra_vn A data.table with the vector normalized spectra
#' @param spectra_sw A data.table with the summed wavelet spectra
#' @param out_path A path to export the files

#'------------------------------------------------------------------------------
#' @example 

root_path <- "/media/antonio/Work/MMDB"
fia <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))
master_comm <- readRDS(paste0(root_path, "/FIA/01_FIA-reshape/fia_reshape.rds"))
master_spectra <- fread(paste0(root_path, "/Pixel_extraction/02_compile-files/EMIT_master_spectra.csv"), header = TRUE)
spectra_sg <- fread(paste0(root_path, "/Pixel_extraction/02_compile-files/02_EMIT_spectra-SG.txt"), header = TRUE)
spectra_vn <- fread(paste0(root_path, "/Pixel_extraction/02_compile-files/02_EMIT_spectra-VN.txt"), header = TRUE)
spectra_sw <- fread(paste0(root_path, "/Pixel_extraction/02_compile-files/02_EMIT_spectra-SW.txt"), header = TRUE)
out_path <- paste0(root_path, "/Pixel_extraction/03_match-spectra")

match_datasets_EMIT(fia, 
                    master_comm,
                    master_spectra,
                    spectra_sg,
                    spectra_vn,
                    spectra_sw,
                    out_path)

#'------------------------------------------------------------------------------
#' @Function

match_datasets_EMIT <- function(fia, 
                                master_comm,
                                master_spectra,
                                spectra_sg,
                                spectra_vn,
                                spectra_sw,
                                out_path) {
  
  #------------------------------------
  # Use samples that were checked
  master_spectra$scene <- substr(master_spectra$scene, 1, 44)
  master_spectra$date <- as.IDate(paste0(substr(master_spectra$scene, 18, 21), "-",
                                         substr(master_spectra$scene, 22, 23), "-",
                                         substr(master_spectra$scene, 24, 25)))
  master_spectra$sample <- 1:nrow(master_spectra)
  
  #------------------------------------
  # Master file for merge
  # Master community file to merge
  master_comm <- master_comm[, .(n_sp = length(unique(SPCD))),
                             by = "PLT_CN"]
  
  fia <- fia[MEASYEAR >= 2018]
  fia$date <- as.IDate(paste0(fia$MEASYEAR, "-",
                              fia$MEASMON, "-",
                              fia$MEASDAY))
  fia <- fia[, .N, by = c("PLT_CN", "date", "LAT", "LON")]
  fia <- fia[, 1:4]
  
  unique_plots <- merge(master_comm, fia, 
                        by = "PLT_CN", 
                        all.x = TRUE, 
                        all.y = FALSE)
  
  #------------------------------------
  #Create unique plots to fill
  unique_plots$sample <- 0
  unique_plots$plot <- 0
  unique_plots$dtime <- 0
  unique_plots$nobservations <- 0
  unique_plots$observation_date <- as.IDate(NULL)
  
  # Progress bar
  pb <- txtProgressBar(min = 1, max = nrow(unique_plots), style = 3)
  
  # Search for spectra #3970
  for(i in 1:nrow(unique_plots)) {
    
    # update progress bar
    setTxtProgressBar(pb, i)
    
    # Subset plot
    sub_spectra <- subset(master_spectra, PLT_CN == unique_plots$PLT_CN[i])
    
    # If no observations
    if(nrow(sub_spectra) == 0) {
      
      unique_plots$sample[i] <- NA
      unique_plots$plot[i] <- NA
      unique_plots$dtime[i] <- NA
      unique_plots$observation_date[i] <- NA
      
      # If there is one or mode observation
    } else {
      
      # If there is one observation
      if(nrow(sub_spectra) == 1) {
        
        unique_plots$sample[i] <- sub_spectra$sample[1]
        unique_plots$plot[i] <- i
        unique_plots$dtime[i] <- sub_spectra$date[1] - unique_plots$date[i]
        unique_plots$nobservations[i] <- 1
        unique_plots$observation_date[i] <- sub_spectra$date[1]
        
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
        unique_plots$dtime[i] <- sub_spectra$date[which.min(abs(delta_time))] - unique_plots$date[i]
        unique_plots$observation_date[i] <- sub_spectra$date[which.min(abs(delta_time))]
        
      }
    }
  }
  
  # Export match spectra
  fwrite(unique_plots, paste0(out_path, "/EMIT/EMIT_plot_match_all.csv"))
  # unique_plots <- fread(paste0(out_path, "/EMIT_plot_match_all.csv"))
  
  #Plots with useful spectra
  unique_plots <- na.exclude(unique_plots)
  master_spectra <- master_spectra[unique_plots$sample, ]
  master_spectra <- master_spectra[, c("scene", "AOD550", "H2O (g cm-2)",
                                       "NDVI", "NIR", "EVI", "CCI", "CRI700", "VOG2", "NPCI",
                                       "WBI", "PRI570", "CMS", "RDVI")]
  unique_plots <- cbind(unique_plots, master_spectra)
  
  #------------------------------------
  # Matching files
  
  # Remove labels
  spectra_sg$scene <- substr(spectra_sg$scene, 1, 44)
  spectra_vn$scene <- substr(spectra_vn$scene, 1, 44)
  spectra_sw$scene <- substr(spectra_sw$scene, 1, 44)
  
  # Merge
  spectra_sg <- merge(unique_plots, spectra_sg, by = c("scene", "PLT_CN"),
                      all.x = TRUE, all.y = FALSE)
  spectra_vn <- merge(unique_plots, spectra_vn, by = c("scene", "PLT_CN"),
                      all.x = TRUE, all.y = FALSE)
  spectra_sw <- merge(unique_plots, spectra_sw, by = c("scene", "PLT_CN"),
                      all.x = TRUE, all.y = FALSE)
  
  #Export
  fwrite(unique_plots, paste0(out_path, "/EMIT/EMIT_match_usefull.csv"))
  fwrite(spectra_sg, paste0(out_path, "/EMIT/EMIT_spectra_SG.csv"))
  fwrite(spectra_vn, paste0(out_path, "/EMIT/EMIT_spectra_VN.csv"))
  fwrite(spectra_sw, paste0(out_path, "/EMIT/EMIT_spectra_SW.csv"))
  
}