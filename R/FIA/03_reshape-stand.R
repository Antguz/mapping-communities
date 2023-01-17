#'------------------------------------------------------------------------------
#' @title Reshape tree stands
#'------------------------------------------------------------------------------

#' @description Its reshape tree stands to match with the spectral information
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param out_path
#' @param threads

#root_path <- "E:/DESIS"
root_path <-  "/media/antonio/Work/Projects/DESIS_mapping"
fia <- fread(paste0(root_path, "/00_FIA/trees_WI_MN_20221102.csv"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"))
out_path <- paste0(root_path, "/03_reshape")

#'------------------------------------------------------------------------------
#' @example 

#'------------------------------------------------------------------------------
#' Script

reshape_stands <- function(fia, SPCD, out_path) { 
  
  #Clean FIA
  fia <- subset(fia, STATECD == 27 | STATECD == 55)
  fia <- fia[order(PLT_CN, SUBP)]
  
  #Merge species code
  fia <- merge(fia, SPCD[, c("SPCD", "PLANTS-Code")], by = "SPCD", all.x = TRUE, all.y = TRUE)
  
  #Get Basal Area
  fia$BA <- pi*(fia$DIA/2)^2
  
  #Basal area multiplied by TPA
  fia$BA_TPA <- fia$BA*fia$TPA_UNADJ
  
  #Remove dead trees
  fia <- fia[STATUSCD == 1]
  
  #Add date
  fia$date <- as.IDate(paste0(fia$MEASYEAR, "-", fia$MEASMON, "-", fia$MEASDAY))
  
  #Get plots with complete information
  complete_plots <- fia[, sum(is.na(BA_TPA)), by = "PLT_CN"]
  complete_plots <- complete_plots[V1 == 0]
  
  #Select plots with complete information
  fia <- merge(complete_plots[, 1], 
               fia, 
               by = "PLT_CN", 
               all.x = TRUE, 
               all.y = FALSE)
  
  #Summarize
  fia_summary <- fia[, .(ntrees = .N,
                         BA_TPA = sum(BA_TPA, na.rm = TRUE),
                         BA = sum(BA, na.rm = TRUE), 
                         HT = mean(HT, na.rm = TRUE), 
                         TPA_UNADJ = sum(TPA_UNADJ, na.rm = TRUE), 
                         DRYBIO_BG = sum(DRYBIO_BG, na.rm = TRUE),
                         CARBON_AG = sum(CARBON_AG, na.rm = TRUE), 
                         DRYBIO_AG = sum(DRYBIO_AG, na.rm = TRUE), 
                         CARBON_BG = sum(CARBON_BG, na.rm = TRUE)), 
                     by = c("PLOT", "PLT_CN", "PLANTS-Code", "date")]
  fia_summary <- fia_summary[order(PLOT, PLT_CN)]
  
  #Reshapes
  #ntress
  #ntrees <- fia_summary[, c("PLOT", "PLT_CN", "date", "PLANTS-Code", "ntrees")]
  #ntrees <- dcast(ntrees, PLOT + PLT_CN + date ~ `PLANTS-Code`)
  #ntrees <- setnafill(ntrees, fill = 0)
  
  #BA_TPA
  BA_TPA <- fia_summary[, c("PLOT", "PLT_CN", "date", "PLANTS-Code", "BA_TPA")]
  BA_TPA <- dcast(BA_TPA, PLOT + PLT_CN + date ~ `PLANTS-Code`)
  BA_TPA <- setnafill(BA_TPA, fill = 0)
  
  #BA
  #BA <- fia_summary[, c("PLOT", "PLT_CN", "date", "PLANTS-Code", "BA")]
  #BA <- dcast(BA, PLOT + PLT_CN + date ~ `PLANTS-Code`)
  #BA <- setnafill(BA, fill = 0)
  
  #DryBio
  #DRYBIO_AG <- fia_summary[, c("PLOT", "PLT_CN", "date", "PLANTS-Code", "DRYBIO_AG")]
  #DRYBIO_AG <- dcast(DRYBIO_AG, PLOT + PLT_CN + date ~ `PLANTS-Code`)
  #DRYBIO_AG <- setnafill(DRYBIO_AG, fill = 0)
  
  #Relative abundance
  #ntrees_sum <- rowSums(ntrees[, 4:81])
  #BA_sum <- rowSums(BA[, 4:81])
  #AG_sum <- rowSums(DRYBIO_AG[, 4:81])
  
  #rntrees <- ntrees[, 4:81]/ntrees_sum
  #rBA <- BA[, 4:81]/BA_sum
  #rDRYBIO_AG <- DRYBIO_AG[, 4:81]/AG_sum
  
  #Final to export
  #rntrees <- cbind(ntrees[, 1:3], rntrees)
  #rBA <- cbind(BA[, 1:3], rBA)
  #rDRYBIO_AG <- cbind(DRYBIO_AG[, 1:3], rDRYBIO_AG)
  
  #Export
  fwrite(BA_TPA, paste0(out_path, "/BA_TPA-reshape.csv"))
  #fwrite(BA, paste0(out_path, "/BA_reshape.csv"))
  #fwrite(DRYBIO_AG, paste0(out_path, "/DRYBIO_AG_reshape.csv"))
  
}
