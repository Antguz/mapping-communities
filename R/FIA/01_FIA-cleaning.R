#'------------------------------------------------------------------------------
#' @title Cleaning and reshaping FIA data
#'------------------------------------------------------------------------------

#' @description It clean and reshape FIA to estimate multiple dimensions of beta
#' diversity
#' 
#' @return Return a frame as a .rds file

#'------------------------------------------------------------------------------
#' @Library

library(data.table)
library(bit64)

#'------------------------------------------------------------------------------
#' @Arguments
#' @param fia A data.table of the FIA dataset
#' @param SPCD A data.table with information of SPCD of FIA and species names
#' @param domain A data.table describing the NEON domain of each plot (optional)
#' @param outpath A path to export

#'------------------------------------------------------------------------------
#' @Function
fia_reshaping <- function(fia,
                          SPCD,
                          domain,
                          outpath) {
  
  ### Prepare plots (Cleaning) -------------------------------------------------
  
  # Read fia
  fia <- fia[order(PLT_CN, SUBP)]
  
  # Inventories from 2018 until now
  fia <- fia[MEASYEAR >= 2018]
  
  # Remove dead trees
  fia <- fia[STATUSCD == 1]
  
  # Get forest metrics by TPA
  fia$BA <- pi*(fia$DIA/2)^2
  fia$BA_TPA <- fia$BA*fia$TPA_UNADJ
  fia$DRYBIO_AG_TPA <- fia$DRYBIO_AG*fia$TPA_UNADJ
  fia$VOL <- pi*(round(((fia$DIA*0.0833333333)/2), 6)^2)*fia$HT
  fia[VOL >= 600, VOL := NA,] #Potential errors
  fia$VOL_TPA <- fia$VOL*fia$TPA_UNADJ
  
  # Select plots with complete volume information
  complete_plots <- fia[, sum(is.na(VOL_TPA)), by = "PLT_CN"]
  complete_plots <- complete_plots[V1 == 0]
  
  fia <- merge(complete_plots[, 1], 
               fia, 
               by = "PLT_CN", 
               all.x = TRUE, 
               all.y = FALSE)
  
  # Select plots with complete tree identification
  fia <- merge(fia, 
               SPCD[, c("SPCD", "Taxon2")], 
               by = "SPCD", 
               all.x = TRUE, 
               all.y = TRUE)
  
  fia[is.na(Taxon2), Remove := 1]
  fia[is.na(Remove), Remove := 0]
  
  complete_plots <- fia[, sum(Remove), by = "PLT_CN"]
  complete_plots <- complete_plots[V1 == 0]
  
  fia <- merge(complete_plots[, 1], 
               fia, 
               by = "PLT_CN", 
               all.x = TRUE, 
               all.y = FALSE)
  
  # Summarize
  fia_summary <- fia[, .(ntrees = .N,
                         BA_TPA = sum(BA, na.rm = TRUE), 
                         DRYBIO_AG_TPA = sum(DRYBIO_AG_TPA, na.rm = TRUE),
                         VOL_TPA = sum(VOL_TPA, na.rm = TRUE)),
                     by = c("PLT_CN", "SPCD")]
  fia_summary <- fia_summary[order(PLT_CN)]
  
  # Select plots with proper boundaries within domains
  colnames(domain)[20] <- "Domain"
  
  fia_summary <- merge(domain[, c("PLT_CN")], 
                       fia_summary, 
                       by = "PLT_CN", 
                       all.x = FALSE, 
                       all.y = FALSE)
  
  # Export 
  saveRDS(fia_summary, paste0(out_path, "/fia_reshape.rds"))
  
}

#'------------------------------------------------------------------------------
#' @Example 

root_path <- "/media/antonio/Work/MMDB/FIA"
fia <- fread(paste0(root_path, "/00_FIA/trees_eastern_20230427.csv"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"), na.strings="NA")
domain <- fread(paste0(root_path, "/00_FIA/plot_domain.csv"))
out_path <- paste0(root_path, "/01_FIA-reshape")

fia_reshaping(fia,
              SPCD,
              domain,
              outpath)