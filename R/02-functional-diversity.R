#'------------------------------------------------------------------------------
#' @title Functional Diversity
#'------------------------------------------------------------------------------

#' @description Estimate the functional diversity of a set of plots
#' @return .txt files with names according the vector file

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(fundiversity)
library(foreach)
library(doParallel)

#'------------------------------------------------------------------------------
#' Arguments
#' @param input_path
#' @param out_path
#' @param threads

#'------------------------------------------------------------------------------
#' @example 

input_path <- "/media/antonio/antonio_ssd/mapping_species/level4"
out_path <- "/media/antonio/antonio_ssd/mapping_species/level5"

input_path <- "E:/mapping_species/level4"
out_path <- "E:/mapping_species/level5"
threads = 6

#'------------------------------------------------------------------------------
#' Function
functional_diversity <- function(input_path, out_path, threads = 24) {
  
  #Search for paths
  files <- list.files(path = input_path, 
                      pattern = ".txt", 
                      all.files = TRUE,
                      full.names = FALSE, 
                      recursive = TRUE)
  files <- files[order(files)]
  
  #Arrange files frame
  frame <- data.table(matrix(unlist(strsplit(files, "-")), 
                             nrow= length(files), 
                             byrow=TRUE), stringsAsFactors=FALSE)
  colnames(frame) <- c("tile", "stand")
  frame[, stand := strsplit(stand, ".txt")[[1]][1], 
        by = seq_along(1:nrow(frame))]
  
  #Unique values
  unique_stands <- unique(frame$stand)
  
  #No memory
  options(fundiversity.memoise = TRUE)
  
  # Set up cluster
  cl <- makeCluster(threads) #, type = "FORK"
  registerDoParallel(cl)
  
  # Loop over tiles ------------------------------------------------------------
  extraction <- foreach(i = 1:length(unique_stands),
                        .combine = rbind,
                        .packages = c("data.table", "fundiversity"),
                        .inorder = F) %dopar% {
                          
                          sub_stand <- subset(frame, stand == unique_stands[i])
                          
                          plot <- NA
                          
                          #Read stand
                          if(nrow(sub_stand) == 1) {
                            
                            plot <- fread(paste0(input_path, "/",
                                                 sub_stand$tile[1], "-",
                                                 sub_stand$stand[1], ".txt"))
                            
                          } else {
                            
                            for(ii in 1:nrow(sub_stand)) {
                              
                              sub_plot <- fread(paste0(input_path, "/",
                                                       sub_stand$tile[ii], "-",
                                                       sub_stand$stand[ii], ".txt"))
                              
                              plot <- rbind(plot, sub_plot)
                              
                              rm("sub_plot")
                              
                            }
                          }
                          
                          #Scale back
                          plot$CCI <- plot$CCI/10000
                          plot$CRE <- plot$CRE/1000
                          plot$NDW <- plot$NDW/10000
                          plot$KNV <- plot$KNV/10000
                          
                          #Subset by NDVI
                          mature <- subset(plot, date == as.IDate("2021-07-12"))
                          mature$KNV <- mature$KNV >= tanh(0.5^2)
                          mature <- mature[KNV == TRUE]
                          mature <- mature[, c(1, 2)]
                          
                          #Remove pixels with low NDVI
                          plot <- merge(mature, plot, by = c("STAND_KEY", "cell"), all.x = TRUE, all.y = FALSE)
                          
                          #Unique pixels and dates
                          unique_date <- unique(plot$date)
                          
                          #Collect results
                          results <- data.table()
                          
                          #Loop of metrics over dates
                          for(ii in 1:length(unique_date)) {
                            
                            #Subset date
                            sub_date <- subset(plot, date == unique_date[ii])
                            sub_date <- na.exclude(sub_date)
                            
                            #Basic information
                            date <- unique_date[ii]
                            npixels <- nrow(sub_date)
                            
                            #Summary mean
                            mean_KND <- mean(sub_date$KNV)
                            mean_CCI <- mean(sub_date$CCI)
                            mean_CRE <- mean(sub_date$CRE)
                            mean_NDW <- mean(sub_date$NDW)
                            
                            #Summary sd
                            sd_KND <- sd(sub_date$KNV)
                            sd_CCI <- sd(sub_date$CCI)
                            sd_CRE <- sd(sub_date$CRE)
                            sd_NDW <- sd(sub_date$NDW)
                            
                            #Functional Richness
                            FRic_CCI <- fd_fric(sub_date[, c(4)])[1, 2]
                            FRic_CRE <- fd_fric(sub_date[, c(5)])[1, 2]
                            FRic_NDW <- fd_fric(sub_date[, c(6)])[1, 2]
                            FRic_all <- fd_fric(sub_date[, c(4:6)])[1, 2]
                            
                            #Functional Divergence
                            FDiv_CCI <- fd_fdiv(sub_date[, c(4)])[1, 2]
                            FDiv_CRE <- fd_fdiv(sub_date[, c(5)])[1, 2]
                            FDiv_NDW <- fd_fdiv(sub_date[, c(6)])[1, 2]
                            FDiv_all <- fd_fdiv(sub_date[, c(4:6)])[1, 2]
                            
                            #Functional Evenness
                            FEve_CCI <- fd_feve(sub_date[, c(4)])[1, 2]
                            FEve_CRE <- fd_feve(sub_date[, c(5)])[1, 2]
                            FEve_NDW <- fd_feve(sub_date[, c(6)])[1, 2]
                            FEve_all <- fd_feve(sub_date[, c(4:6)])[1, 2]
                            
                            #Functional Dispersion
                            FDis_CCI <- fd_fdis(sub_date[, c(4)])[1, 2]
                            FDis_CRE <- fd_fdis(sub_date[, c(5)])[1, 2]
                            FDis_NDW <- fd_fdis(sub_date[, c(6)])[1, 2]
                            FDis_all <- fd_fdis(sub_date[, c(4:6)])[1, 2]
                            
                            #Pairwise mean distance
                            Pmd_CCI <- mean(dist(sub_date[, c(4)], method = "euclidean", diag = FALSE, upper = FALSE))
                            Pmd_CRE <- mean(dist(sub_date[, c(5)], method = "euclidean", diag = FALSE, upper = FALSE))
                            Pmd_NDW <- mean(dist(sub_date[, c(6)], method = "euclidean", diag = FALSE, upper = FALSE))
                            Pmd_all <- mean(dist(sub_date[, c(4:6)], method = "euclidean", diag = FALSE, upper = FALSE))
                            
                            #Rao’s Quadratic Entropy
                            Rao_all <- fd_raoq(sub_date[, c(4:6)])[1, 2]
                            
                            #All together
                            sub_results <- data.table(STAND_KEY = sub_date$STAND_KEY[1],
                                                      date = date,
                                                      npixels = npixels,
                                                      mean_KND = mean_KND,
                                                      mean_CCI = mean_CCI,
                                                      mean_CRE = mean_CRE,
                                                      mean_NDW = mean_NDW,
                                                      sd_KND = sd_KND,
                                                      sd_CCI = sd_CCI,
                                                      sd_CRE = sd_CRE,
                                                      sd_NDW = sd_NDW,
                                                      FRic_CCI = FRic_CCI,
                                                      FRic_CRE = FRic_CRE,
                                                      FRic_NDW = FRic_NDW,
                                                      FRic_all = FRic_all,
                                                      FDiv_CCI = FDiv_CCI,
                                                      FDiv_CRE = FDiv_CRE,
                                                      FDiv_NDW = FDiv_NDW,
                                                      FDiv_all = FDiv_all,
                                                      FEve_CCI = FEve_CCI,
                                                      FEve_CRE = FEve_CRE,
                                                      FEve_NDW = FEve_NDW,
                                                      FEve_all = FEve_all,
                                                      FDis_CCI = FDis_CCI,
                                                      FDis_CRE = FDis_CRE,
                                                      FDis_NDW = FDis_NDW,
                                                      FDis_all = FDis_all,
                                                      Pmd_CCI = Pmd_CCI,
                                                      Pmd_CRE = Pmd_CRE,
                                                      Pmd_NDW = Pmd_NDW, 
                                                      Pmd_all = Pmd_all,
                                                      Rao_all = Rao_all)
                            
                            results <- rbind(results, sub_results)
                            
                            rm(list = c("sub_date", "date", "npixels",
                                        "mean_KND", "mean_CCI", "mean_CRE", "mean_NDW",
                                        "sd_KND", "sd_CCI", "sd_CRE", "sd_NDW",
                                        "FRic_CCI", "FRic_CRE", "FRic_NDW", "FRic_all",
                                        "FDiv_CCI", "FDiv_CRE", "FDiv_NDW", "FDiv_all",
                                        "FEve_CCI", "FEve_CRE", "FEve_NDW", "FEve_all",
                                        "FDis_CCI", "FDis_CRE", "FDis_NDW", "FDis_all",
                                        "Pmd_CCI", "Pmd_CRE", "Pmd_NDW", "Pmd_all",
                                        "Rao_all", "sub_results"))
                            
                          }
                          
                          #Name to export
                          name_export <- paste0(out_path, "/", 
                                                "stand_", 
                                                unique_stands[i],".txt")
                          
                          #Export
                          fwrite(results, 
                                 name_export, 
                                 sep = "\t")
                          
                          
                          #Remove residuals
                          rm(list = c("sub_stand", "plot", "mature", "unique_date", 
                                      "name_export", "results"))
                          gc()
                          
                        }
  
  #Stop cluster
  stopCluster(cl)
  gc()
  
}

functional_diversity(input_path, out_path, threads = 7)



