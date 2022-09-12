#'------------------------------------------------------------------------------
#' @title Metrics structural diversity
#'------------------------------------------------------------------------------

#' @description It estimates metrics of the structural diversity
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(fundiversity)
library(foreach)
library(doParallel)

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param out_path
#' @param threads

fia <- fread("data/FIA_SPCD.csv")
out_path <- "/media/antonio/antonio_ssd/mapping_species/results"
threads <- 30

#'------------------------------------------------------------------------------
#' @example 

#'------------------------------------------------------------------------------
#' Function

fia_diversity <- function(fia, out_path, threads) {
  
  #Clean FIA
  fia <- subset(fia, STATECD == 27 | STATECD == 55)
  fia <- fia[order(PLOTID, SUBP)]
  
  #FIA cordinates
  fia_coordinates <- fia[, .(LAT = mean(LAT, na.rm = TRUE), LON = mean(LON, na.rm = TRUE)), by = c("PLOTID", "SUBP")]
  fia_coordinates <- fia_coordinates[SUBP == "1"]
  fia_coordinates <- fia_coordinates[order(PLOTID)]
  fia_coordinates$ID <- 1:nrow(fia_coordinates)
  
  #Summarize
  fia_summary <- fia[, .(BA = sum(BA, na.rm = TRUE), 
                         TPA = sum(TPA, na.rm = TRUE), 
                         SAPLING = sum(SAPLING, na.rm = TRUE), 
                         BIGTREE = sum(BIGTREE, na.rm = TRUE),
                         TREETALLY = sum(TREETALLY, na.rm = TRUE), 
                         TOTALBIOMASS = sum(TOTALBIOMASS, na.rm = TRUE)), 
                     by = c("PLOTID", "SPCD")]
  fia_summary <- fia_summary[order(PLOTID)]
  
  #Unique stands
  plots <- unique(fia_coordinates$PLOTID)
  
  #No memory
  options(fundiversity.memoise = TRUE)
  
  #shannon
  shannon <- function(n_points) {
    p.i <- n_points/sum(n_points)
    H <- (-1) * sum(p.i * log(p.i))
    return(H)
  }
  
  # Set up cluster
  cl <- makeCluster(threads, type = "FORK")
  registerDoParallel(cl)
  
  # Loop over STAND_KEY---------------------------------------------------------
  stand <- foreach(i = 1:length(plots),
                   .combine = rbind,
                   .export = c("shannon"),
                   .packages = c("data.table", "fundiversity"),
                   .inorder = F) %dopar% {
                     
                     sub_plot <- subset(fia_summary, PLOTID == fia_coordinates$PLOTID[i])
                     
                     final <- data.table(PLOTID = fia_coordinates$PLOTID[i],
                                         ID = fia_coordinates$ID[i],
                                         LAT = fia_coordinates$LAT[i],
                                         LON = fia_coordinates$LON[i],
                                         nsp = nrow(sub_plot),
                                         TBA = sum(sub_plot$BA),
                                         TPA = sum(sub_plot$TPA),
                                         SAPLING = sum(sub_plot$SAPLING),
                                         BIGTREE = sum(sub_plot$BIGTREE),
                                         TREETALLY = sum(sub_plot$TREETALLY),
                                         TB = sum(sub_plot$TOTALBIOMASS),
                                         H_BA = shannon(sub_plot$BA),
                                         H_TPA = shannon(sub_plot$BA),
                                         H_TOTALBIOMASS = shannon(sub_plot$BA),
                                         H_BIGTREE = shannon(sub_plot$BIGTREE),
                                         FRic_TB = fd_fric(sub_plot$TOTALBIOMASS)[1, 2],
                                         FDiv_TB = fd_fdiv(sub_plot$TOTALBIOMASS)[1, 2],
                                         FEve_TB = fd_feve(sub_plot$TOTALBIOMASS)[1, 2],
                                         FDis_TB = fd_fdis(sub_plot$TOTALBIOMASS)[1, 2],
                                         Rao_TB  = fd_raoq(sub_plot$TOTALBIOMASS)[1, 2],
                                         FRic_BA = fd_fric(sub_plot$BA)[1, 2],
                                         FDiv_BA = fd_fdiv(sub_plot$BA)[1, 2],
                                         FEve_BA = fd_feve(sub_plot$BA)[1, 2],
                                         FDis_BA = fd_fdis(sub_plot$BA)[1, 2],
                                         Rao_BA  = fd_raoq(sub_plot$BA)[1, 2],
                                         flag_sp = any(sub_plot$SPCD == 299) |  any(sub_plot$SPCD == 998) | any(sub_plot$SPCD == 999),
                                         flag_ba = any(is.na(sub_plot$BA) == TRUE),
                                         flag_tpa = any(is.na(sub_plot$tpa) == TRUE),
                                         flag_tb = any(is.na(sub_plot$tb) == TRUE))
                     
                     return(final)
                     
                   }
  
  #Stop cluster
  stopCluster(cl)
  gc()

  #Export
  fwrite(stand, paste0(out_path, "/master_FIA-diversity.txt"), sep = "\t")
  
}

#Cords 12484
#MBF    20
#Both  3093
