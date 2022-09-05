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
#' @param master
#' @param out_path
#' @param threads

master <- fread("/media/antonio/antonio_ssd/mapping_species/results/master_clean.txt")
out_path <- "/media/antonio/antonio_ssd/mapping_species/results"
threads <- 30

#'------------------------------------------------------------------------------
#' @example 

#'------------------------------------------------------------------------------
#' Function

stand_diversity <- function(master, threads) {
  
  #Subset Stem and flags
  master <- subset(master, Unique_UOM != "Stems")
  master <- subset(master, DBH_flag != 0)
  master <- subset(master, VOL_flag != 0)
  
  #Subset by Volume error
  master <- subset(master, TOTAL_VOL <= 70)
  
  #Unique stands
  plots <- unique(master$STAND_KEY)
  
  #No memory
  options(fundiversity.memoise = TRUE)
  
  # Set up cluster
  cl <- makeCluster(threads, type = "FORK")
  registerDoParallel(cl)
  
  # Loop over STAND_KEY---------------------------------------------------------
  stand <- foreach(i = 1:length(plots),
                   .combine = rbind,
                   .packages = c("data.table", "fundiversity"),
                   .inorder = F) %dopar% {
                     
                     sub_plot <- subset(master, STAND_KEY == plots[i])
                     final <- sub_plot[1, c(1:29, 34, 35, 36, 37)]
                     
                     #VOL
                     final$FRic_VOL <- fd_fric(sub_plot$VOL)[1, 2]
                     final$FDiv_VOL <- fd_fdiv(sub_plot$VOL)[1, 2]
                     final$FEve_VOL <- fd_feve(sub_plot$VOL)[1, 2]
                     final$FDis_VOL <- fd_fdis(sub_plot$VOL)[1, 2]
                     final$Rao_VOL  <- fd_raoq(sub_plot$VOL)[1, 2]
                     
                     #DBH
                     final$FRic_DBH <- fd_fric(sub_plot$DBH)[1, 2]
                     final$FDiv_DBH <- fd_fdiv(sub_plot$DBH)[1, 2]
                     final$FEve_DBH <- fd_feve(sub_plot$DBH)[1, 2]
                     final$FDis_DBH <- fd_fdis(sub_plot$DBH)[1, 2]
                     final$Rao_DBH  <- fd_raoq(sub_plot$DBH)[1, 2]
                     
                     return(final)
                     
                   }
  
  #Stop cluster
  stopCluster(cl)
  gc()

  fwrite(stand, paste0(out_path, "/master_stand-diversity.txt"), sep = "\t")
  
}

stand[, .(count = .N), by = Unique_UOM]

#Cords 12484
#MBF    20
#Both  3093
