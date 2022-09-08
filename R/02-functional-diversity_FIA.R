#'------------------------------------------------------------------------------
#' @title Functional Diversity
#'------------------------------------------------------------------------------

#' @description Estimate the functional diversity of a set of plots
#' @return .txt files with names according the vector file

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(fundiversity)
library(doParallel)
library(terra)

#'------------------------------------------------------------------------------
#' Arguments
#' @param input_path
#' @param out_path
#' @param threads
#' @param chuck

#'------------------------------------------------------------------------------
#' @example 

input_path <- "E:/mapping_species/level4"
out_path <- "E:/mapping_species/level5"
threads <- 100

#input_path <- "/panfs/jay/groups/17/cavender/shared/oak-wilt/level4"
#out_path <- "/panfs/jay/groups/17/cavender/shared/oak-wilt/level5"
#chuck <- 1000
#threads <- 110

tile <- vect("data/title_clipped.gpkg")
fsi <- fread("data/FIA_SPCD.csv")
input_path <- "/media/antonio/antonio_ssd/mapping_species/level4"
out_path <- "/media/antonio/antonio_ssd/mapping_species/level5"
chuck <- 100
threads <- 30

#'------------------------------------------------------------------------------
#' Function
functional_diversity <- function(tile, input_path, out_path, chuck = 40, threads = 30) {
  
  #Clean fsi
  fsi <- subset(fsi, STATECD == 27 | STATECD == 55)
  fsi <- subset(fsi, SUBP == 1)
  
  #Make a vect
  fsi <- vect(fsi, geom=c("LON", "LAT"), crs= "+proj=longlat +datum=WGS84")
  fsi <- project(fsi, tile)
  fsi <- buffer(fsi, width = 100)
  tile_fsi <- intersect(fsi, tile)
  ids <- data.table(ID = 1:dim(tile_fsi)[1],
                    PLOTID = tile_fsi$PLOTID)
  frame <- as.data.table(as.data.frame(tile_fsi))
  frame <- frame[, c(1, 29:31)]
  frame <- frame[order(Tile_ID)]
  
  #Unique tiles
  tile <- unique(frame$Tile_ID)
  
  #dates
  unique_dates <- as.IDate(c("2021-01-01", "2021-01-17", "2021-02-02", "2021-02-18", "2021-03-06",
                      "2021-03-22", "2021-04-07", "2021-04-23", "2021-05-09", "2021-05-25",
                      "2021-06-10", "2021-06-26", "2021-07-12", "2021-07-28", "2021-08-13",
                      "2021-08-29", "2021-09-14", "2021-09-30", "2021-10-16", "2021-11-01",
                      "2021-11-17", "2021-12-03"))

  
  #No memory
  options(fundiversity.memoise = TRUE)
  
  # Loop over tiles ------------------------------------------------------------
  cl <- makeCluster(threads, type = "FORK")
  registerDoParallel(cl)
  
  # i <- 15
  # Loop over tiles ------------------------------------------------------------
  extraction <- foreach(i = 1:length(tile),
                        .packages = ("data.table"),
                        .inorder = F) %dopar% {
                          
                          #Get tile
                          X <- as.numeric(substr(tile[i], 2, 5))
                          Y <- as.numeric(substr(tile[i], 8, 11))
                          
                          #Get plots
                          sub_frame <- subset(frame, Tile_X <= (X+1) &
                                                     Tile_X >= (X-1) &
                                                     Tile_Y <= (Y+1) &
                                                     Tile_Y >= (Y-1))
                          
                          unique_section <- unique(sub_frame$Tile_ID)
                          
                          #Loop to compile tiles
                          pixels <- data.table()
                          
                          for(ii in 1:length(unique_section)) {
                            
                            extraction <- fread(paste0(input_path, "/", unique_section[ii], ".txt"))
                            pixels <- rbind(pixels, extraction)
                            
                          }
                          
                          pixels <- merge(pixels[, c(1, 3:11)], ids, all.x = TRUE, all.y = FALSE)
                          pixels <- pixels[, c(1, 11, 2:10)]
                          colnames(pixels) <- c("ID", "PLOTID", "cell", "weight", "x", "y", "date", "CCI", "CRE","NDW", "KNV")
                          
                          #Clean pixels
                          pixels[, c(8, 10:11)] <- pixels[, c(8, 10:11)]/10000
                          pixels[, c(9)] <- pixels[, c(9)]/1000
                          
                          #Unique PLOTID for tile of interest
                          plotid <- subset(frame, Tile_ID == tile[i])
                          plotid <- unique(plotid$PLOTID)
                          
                          #Loop over plots
                          for(k in 1:length(plotid)) {
                            
                            mask <- subset(pixels, PLOTID == plotid[k])
                            
                            #Compile results
                            results <- data.table()
                            
                            #Loop over dates
                            for(j in 1:length(unique_dates)) {
                              
                              mask_date <- subset(mask, date == unique_dates[j])
                              
                              #N_pixels
                              if(nrow(mask_date) >= 3) {
                                
                                #Basic information
                                date <- unique_dates[j]
                                npixels <- nrow(mask_date)
                                
                                #Summary mean
                                mean_KND <- mean(mask_date$KNV, na.rm = TRUE)
                                mean_CCI <- mean(mask_date$CCI, na.rm = TRUE)
                                mean_CRE <- mean(mask_date$CRE, na.rm = TRUE)
                                mean_NDW <- mean(mask_date$NDW, na.rm = TRUE)
                                
                                #Summary sd
                                sd_KND <- sd(mask_date$KNV, na.rm = TRUE)
                                sd_CCI <- sd(mask_date$CCI, na.rm = TRUE)
                                sd_CRE <- sd(mask_date$CRE, na.rm = TRUE)
                                sd_NDW <- sd(mask_date$NDW, na.rm = TRUE)
                                
                                #Functional Richness
                                FRic_CCI <- fd_fric(mask_date[, c(8)])[1, 2]
                                FRic_CRE <- fd_fric(mask_date[, c(9)])[1, 2]
                                FRic_NDW <- fd_fric(mask_date[, c(10)])[1, 2]
                                FRic_all <- fd_fric(mask_date[, c(8:10)])[1, 2]
                                
                                #Functional Divergence
                                FDiv_CCI <- fd_fdiv(mask_date[, c(8)])[1, 2]
                                FDiv_CRE <- fd_fdiv(mask_date[, c(9)])[1, 2]
                                FDiv_NDW <- fd_fdiv(mask_date[, c(10)])[1, 2]
                                FDiv_all <- fd_fdiv(mask_date[, c(8:10)])[1, 2]
                                
                                #Functional Evenness
                                FEve_CCI <- fd_feve(mask_date[, c(8)])[1, 2]
                                FEve_CRE <- fd_feve(mask_date[, c(9)])[1, 2]
                                FEve_NDW <- fd_feve(mask_date[, c(10)])[1, 2]
                                FEve_all <- fd_feve(mask_date[, c(8:10)])[1, 2]
                                
                                #Functional Dispersion
                                FDis_CCI <- fd_fdis(mask_date[, c(8)])[1, 2]
                                FDis_CRE <- fd_fdis(mask_date[, c(9)])[1, 2]
                                FDis_NDW <- fd_fdis(mask_date[, c(10)])[1, 2]
                                FDis_all <- fd_fdis(mask_date[, c(8:10)])[1, 2]
                                
                                #Pairwise mean distance
                                Pmd_CCI <- mean(dist(mask_date[, c(8)], method = "euclidean", diag = FALSE, upper = FALSE))
                                Pmd_CRE <- mean(dist(mask_date[, c(9)], method = "euclidean", diag = FALSE, upper = FALSE))
                                Pmd_NDW <- mean(dist(mask_date[, c(10)], method = "euclidean", diag = FALSE, upper = FALSE))
                                Pmd_all <- mean(dist(mask_date[, c(8:10)], method = "euclidean", diag = FALSE, upper = FALSE))
                                
                                #Rao’s Quadratic Entropy
                                Rao_all <- fd_raoq(mask_date[, c(8:10)])[1, 2]
                                
                                #All together
                                sub_results <- data.table(PLOTID = mask_date$PLOTID[1],
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
                                
                              }
                              
                              rm(list = c("npixels", "mean_KND", "mean_CCI", "mean_CRE", "mean_NDW",
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
                                                  "PLOTID_", 
                                                  plotid[k],".txt")
                            
                            #Export
                            fwrite(results, 
                                   name_export, 
                                   sep = "\t")
                            
                            rm(list = c("j"))
                          
                          }
                          
                          rm(list = c("X", "Y", "sub_frame", "unique_section", 
                                      "pixels", "ii", "extraction", "plotid",
                                      "k"))

    }
  
}
