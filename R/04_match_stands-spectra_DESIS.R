#'------------------------------------------------------------------------------
#' @title Match spectra from DESIS and stands
#'------------------------------------------------------------------------------

#' @description Its match the spectral information from DESIS with the reshape 
#' stands
#' 
#' @return Frame files with best batch between datasets

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(bit64)
#options(scipen = 10000000000)

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
master_comm <- readRDS(paste0(root_path, "/03_reshape/fia_reshape.rds"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"))
master_spectra <- fread(paste0(root_path, "/02_compile-files/DESIS_master_spectra.csv"), header = TRUE)
spectra_sg <- fread(paste0(root_path, "/02_compile-files/02_DESIS_spectra-SG.txt"), header = TRUE)
spectra_vn <- fread(paste0(root_path, "/02_compile-files/02_DESIS_spectra-VN.txt"), header = TRUE)
spectra_sw <- fread(paste0(root_path, "/02_compile-files/02_DESIS_spectra-SW.txt"), header = TRUE)
out_path <- paste0(root_path, "/04_match_spectra-stand")

match_datasets_DESIS(fia, stand, spectra_sw, master_spectra, out_path)

#'------------------------------------------------------------------------------
#' Script
match_datasets_DESIS <- function(fia, 
                                 stand,
                                 spectra_sg,
                                 spectra_vn,
                                 spectra_sw,
                                 master_spectra,
                                 out_path) {
  
  #------------------------------------
  # Use samples that were checked
  master_spectra$scene <- substr(master_spectra$scene, 1, 46)
  master_spectra$date <- as.IDate(paste0(substr(master_spectra$scene, 32, 35), "-",
                                         substr(master_spectra$scene, 36, 37), "-",
                                         substr(master_spectra$scene, 38, 39)))
  master_spectra$sample <- 1:nrow(master_spectra)

  #------------------------------------
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
  unique_plots$observation_date <- as.IDate("2000-01-01")
  unique_plots$folder <- ""
  
  # Progress bar
  pb <- txtProgressBar(min = 1, max = nrow(unique_plots), style = 3)
  
  # Search for spectra #5307
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
      unique_plots$folder[i] <- NA
      
      # If there is one or mode observation
    } else {
      
      # If there is one observation
      if(nrow(sub_spectra) == 1) {
        
        unique_plots$sample[i] <- sub_spectra$sample[1]
        unique_plots$plot[i] <- i
        unique_plots$dtime[i] <- sub_spectra$date[1] - unique_plots$date[i]
        unique_plots$nobservations[i] <- 1
        unique_plots$observation_date[i] <- sub_spectra$date[1]
        unique_plots$folder[i] <- sub_spectra$folder[1]
        
        # If there is two or more observation
      } else if(nrow(sub_spectra) > 1) {
        
        unique_plots$nobservations[i] <- nrow(sub_spectra)
        
        # For DESIS avoid observations on 2023
        if(length(unique(sub_spectra[year(date) != 2023, year(date)])) > 0) {
          sub_spectra <- subset(sub_spectra, year(date) != 2023)
        }
        
        # Select observations closest to the inventory year
        delta_year <- abs(year(sub_spectra$date) - year(unique_plots$date[i]))
        sub_spectra <- sub_spectra[delta_year == delta_year[which.min(delta_year)]]
        
        # Select observations closest to the peak of season
        peak <- as.IDate(paste0(unique(year(sub_spectra$date)), "-07-15"))
        delta_time <- sub_spectra$date-peak
        
        unique_plots$sample[i] <- sub_spectra$sample[which.min(abs(delta_time))]
        unique_plots$plot[i] <- i
        unique_plots$dtime[i] <- delta_time[which.min(delta_time)]
        unique_plots$observation_date[i] <- sub_spectra$date[which.min(abs(delta_time))]
        unique_plots$folder[i] <- sub_spectra$folder[which.min(abs(delta_time))]
        
      }
    }
  }
  
  # Export match spectra
  fwrite(unique_plots, paste0(out_path, "/DESIS_plot_match.csv"))
  
  #Plots with useful spectra
  unique_plots <- na.exclude(unique_plots)
  master_spectra <- master_spectra[unique_plots$sample, ]
  master_spectra <- master_spectra[, c(2, 4:6)]
  unique_plots <- cbind(unique_plots, master_spectra)
  
  #------------------------------------
  # Matching files
  
  # Remove labels
  spectra_sg$scene <- substr(spectra_sg$scene, 1, 46)
  spectra_vn$scene <- substr(spectra_vn$scene, 1, 46)
  spectra_sw$scene <- substr(spectra_sw$scene, 1, 46)

  # Merge
  spectra_sg <- merge(unique_plots, spectra_sg, by = c("folder", "scene", "PLT_CN"),
                      all.x = TRUE, all.y = FALSE)
  spectra_vn <- merge(unique_plots, spectra_vn, by = c("folder", "scene", "PLT_CN"),
                      all.x = TRUE, all.y = FALSE)
  spectra_sw <- merge(unique_plots, spectra_sw, by = c("folder", "scene", "PLT_CN"),
                      all.x = TRUE, all.y = FALSE)
  
  spectra_sg <- spectra_sg[, -1]
  spectra_vn <- spectra_vn[, -1]
  spectra_sw <- spectra_sw[, -1]
  
  #Export
  fwrite(unique_plots, paste0(out_path, "/DESIS_match_usefull.csv"))
  fwrite(spectra_sg, paste0(out_path, "/DESIS_spectra_SG.csv"))
  fwrite(spectra_vn, paste0(out_path, "/DESIS_spectra_VN.csv"))
  fwrite(spectra_sw, paste0(out_path, "/DESIS_spectra_SW.csv"))
  
}
