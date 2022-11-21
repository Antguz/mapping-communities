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

root_path <- "E:/DESIS"
root_path <- "/media/antonio/Extreme SSD/DESIS"
fia <- fread(paste0(root_path, "/00_FIA/trees_WI_MN_20221102.csv"))
out_path <- paste0(root_path, "/03_diversity")
threads <- 8

#'------------------------------------------------------------------------------
#' @example 

#'------------------------------------------------------------------------------
#' Function

fia_diversity <- function(fia, out_path, threads) {
  
  #Clean FIA
  fia <- subset(fia, STATECD == 27 | STATECD == 55)
  fia <- fia[order(PLT_CN, SUBP)]
  
  #Get Basal Area
  fia$BA <- pi*(fia$DIA/2)^2
  
  #Remove dead trees
  fia <- fia[is.na(STANDING_DEAD_CD)]
  
  #Summarize
  fia_summary <- fia[, .(BA = sum(BA, na.rm = TRUE), 
                         HT = mean(HT, na.rm = TRUE), 
                         TPA_UNADJ = sum(TPA_UNADJ, na.rm = TRUE), 
                         DRYBIO_BG = sum(DRYBIO_BG, na.rm = TRUE),
                         CARBON_AG = sum(CARBON_AG, na.rm = TRUE), 
                         DRYBIO_AG = sum(DRYBIO_AG, na.rm = TRUE), 
                         CARBON_BG = sum(CARBON_BG, na.rm = TRUE)), 
                     by = c("PLT_CN", "SPCD")]
  fia_summary <- fia_summary[order(PLT_CN)]
  fwrite(fia_summary, paste0(out_path, "/master_FIA-summary.txt"), sep = "\t")
  
  #Unique stands
  plots <- unique(fia$PLT_CN)
  
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
                     
                     sub_plot <- subset(fia, PLT_CN == plots[i])
                     
                     final <- data.table(PLT_CN = plots[i],
                                         n_trees = nrow(sub_plot),
                                         n_sp = length(unique(sub_plot$SPCD)),
                                         T_BA = sum(sub_plot$BA),
                                         T_DRYBIO_BG = sum(sub_plot$DRYBIO_BG),
                                         T_CARBON_AG = sum(sub_plot$CARBON_AG),
                                         T_DRYBIO_AG = sum(sub_plot$DRYBIO_AG),
                                         T_CARBON_BG = sum(sub_plot$CARBON_BG),
                                         H_BA = shannon(sub_plot$BA),
                                         H_DRYBIO_AG = shannon(sub_plot$DRYBIO_AG),
                                         FRic_AG = fd_fric(sub_plot$DRYBIO_AG)[1, 2],
                                         FDiv_AG = fd_fdiv(sub_plot$DRYBIO_AG)[1, 2],
                                         FEve_AG = fd_feve(sub_plot$DRYBIO_AG)[1, 2],
                                         FDis_AG = fd_fdis(sub_plot$DRYBIO_AG)[1, 2],
                                         Rao_AG  = fd_raoq(sub_plot$DRYBIO_AG)[1, 2],
                                         FRic_BA = fd_fric(sub_plot$BA)[1, 2],
                                         FDiv_BA = fd_fdiv(sub_plot$BA)[1, 2],
                                         FEve_BA = fd_feve(sub_plot$BA)[1, 2],
                                         FDis_BA = fd_fdis(sub_plot$BA)[1, 2],
                                         Rao_BA  = fd_raoq(sub_plot$BA)[1, 2])
                     
                     return(final)
                     
                   }
  
  #Stop cluster
  stopCluster(cl)
  gc()

  #Export
  fwrite(stand, paste0(out_path, "/master_FIA-diversity.txt"), sep = "\t")
  
}

