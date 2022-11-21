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

root_path <- "E:/DESIS"
root_path <- "/media/antonio/Extreme SSD/DESIS"
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
  
  #Remove dead trees
  fia <- fia[is.na(STANDING_DEAD_CD)]
  fia$date <- as.IDate(paste0(fia$MEASYEAR, "-", fia$MEASMON, "-", fia$MEASDAY))
  
  #Summarize
  fia_summary <- fia[, .(ntrees = .N,
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
  ntrees <- fia_summary[, c("PLOT", "PLT_CN", "date", "PLANTS-Code", "ntrees")]
  ntrees <- dcast(ntrees, PLOT + PLT_CN + date ~ `PLANTS-Code`)
  ntrees <- setnafill(ntrees, fill = 0)
  
  #BA
  BA <- fia_summary[, c("PLOT", "PLT_CN", "date", "PLANTS-Code", "BA")]
  BA <- dcast(BA, PLOT + PLT_CN + date ~ `PLANTS-Code`)
  BA <- setnafill(BA, fill = 0)
  
  #ntress
  DRYBIO_AG <- fia_summary[, c("PLOT", "PLT_CN", "date", "PLANTS-Code", "DRYBIO_AG")]
  DRYBIO_AG <- dcast(DRYBIO_AG, PLOT + PLT_CN + date ~ `PLANTS-Code`)
  DRYBIO_AG <- setnafill(DRYBIO_AG, fill = 0)
  
  #Export
  fwrite(ntrees, paste0(out_path, "/ntrees_reshape.csv"))
  fwrite(BA, paste0(out_path, "/BA_reshape.csv"))
  fwrite(DRYBIO_AG, paste0(out_path, "/DRYBIO_AG_reshape.csv"))
  
}
