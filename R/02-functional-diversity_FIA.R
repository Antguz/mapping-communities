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
library(bit64)
library(foreach)

#'------------------------------------------------------------------------------
#' Arguments
#' @param input_path
#' @param out_path
#' @param threads
#' @param chuck

#'------------------------------------------------------------------------------
#' @example 

#fia <- fread("data/FIA_SPCD.csv")
#input_path <- "E:/mapping_species/level4"
#out_path <- "E:/mapping_species/level5"
#chuck <- 4
#threads <- 4

#fia <- fread("/panfs/jay/groups/17/cavender/shared/oak-wilt/misc/FIA/FIA_SPCD.csv")
#input_path <- "/panfs/jay/groups/17/cavender/shared/oak-wilt/level4_rast"
#out_path <- "/panfs/jay/groups/17/cavender/shared/oak-wilt/level5"
#chuck <- 100
#threads <- 100

fia <- fread("data/FIA_SPCD.csv")
input_path <- "/media/antonio/antonio_ssd/mapping_species/level4"
out_path <- "/media/antonio/antonio_ssd/mapping_species/level5"
chuck <- 500
threads <- 30

#In case
options("scipen"=100, "digits"=6)

#'------------------------------------------------------------------------------
#' Function
functional_diversity <- function(fia, input_path, out_path, chuck = 100, threads = 100) {
  
  #Search for the tiles -------------------------------------------------------- 
  tiles <- list.files(path = input_path, 
                      pattern = ".txt", 
                      all.files = TRUE,
                      full.names = FALSE, 
                      recursive = TRUE)
  tiles <- tiles[order(tiles)]
  
  #dates
  unique_dates <- as.IDate(c("2021-01-01", "2021-01-17", "2021-02-02", "2021-02-18", "2021-03-06",
                      "2021-03-22", "2021-04-07", "2021-04-23", "2021-05-09", "2021-05-25",
                      "2021-06-10", "2021-06-26", "2021-07-12", "2021-07-28", "2021-08-13",
                      "2021-08-29", "2021-09-14", "2021-09-30", "2021-10-16", "2021-11-01",
                      "2021-11-17", "2021-12-03"))
  
  #No memory
  options(fundiversity.memoise = TRUE)

  # Functional diversity function ----------------------------------------------
  extraction <- function(i, tiles, unique_dates, input_path, out_path) {
    
    print(tiles[i])
    
    #Read tile
    pixels <- fread(paste0(input_path, "/", tiles[i]))
    
    #Clean pixels
    pixels[, c(6, 8:9)] <- pixels[, c(6, 8:9)]/10000
    pixels[, c(7)] <- pixels[, c(7)]/1000
    
    #Unique PLOTID for tile of interest
    plotid <- unique(pixels$PLOTID)
    
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
          FRic_CCI <- fd_fric(mask_date[, c(6)])[1, 2]
          FRic_CRE <- fd_fric(mask_date[, c(7)])[1, 2]
          FRic_NDW <- fd_fric(mask_date[, c(8)])[1, 2]
          FRic_all <- fd_fric(mask_date[, c(6:8)])[1, 2]
          
          #Functional Divergence
          FDiv_CCI <- fd_fdiv(mask_date[, c(6)])[1, 2]
          FDiv_CRE <- fd_fdiv(mask_date[, c(7)])[1, 2]
          FDiv_NDW <- fd_fdiv(mask_date[, c(8)])[1, 2]
          FDiv_all <- fd_fdiv(mask_date[, c(6:8)])[1, 2]
          
          #Functional Evenness
          FEve_CCI <- fd_feve(mask_date[, c(6)])[1, 2]
          FEve_CRE <- fd_feve(mask_date[, c(7)])[1, 2]
          FEve_NDW <- fd_feve(mask_date[, c(8)])[1, 2]
          FEve_all <- fd_feve(mask_date[, c(6:8)])[1, 2]
          
          #Functional Dispersion
          FDis_CCI <- fd_fdis(mask_date[, c(6)])[1, 2]
          FDis_CRE <- fd_fdis(mask_date[, c(7)])[1, 2]
          FDis_NDW <- fd_fdis(mask_date[, c(8)])[1, 2]
          FDis_all <- fd_fdis(mask_date[, c(6:8)])[1, 2]
          
          #Pairwise mean distance
          Pmd_CCI <- mean(dist(mask_date[, c(6)], method = "euclidean", diag = FALSE, upper = FALSE))
          Pmd_CRE <- mean(dist(mask_date[, c(7)], method = "euclidean", diag = FALSE, upper = FALSE))
          Pmd_NDW <- mean(dist(mask_date[, c(8)], method = "euclidean", diag = FALSE, upper = FALSE))
          Pmd_all <- mean(dist(mask_date[, c(6:8)], method = "euclidean", diag = FALSE, upper = FALSE))
          
          #Rao’s Quadratic Entropy
          Rao_all <- fd_raoq(mask_date[, c(6:8)])[1, 2]
          
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
    
    #Remove residuals
    rm(list = c("pixels", "plotid", "k"))

  }
  
  # Loop for chucks ------------------------------------------------------------
  x <- 1
  x_max <- length(tiles)
  n_times <- ceiling(x_max/chuck)
  
  #n_files <- c(14, 15, 16, 28)
  for(j in 5:n_times) {
    
    print(j)
    
    #n_files <- c(14, 15, 16, 28)
    if((j*chuck) <= x_max) {
      n_files <- x:(j*chuck) 
    } else {
      n_files <- ((j-1)*chuck):x_max
    }
    
    # Set up cluster
    cl <- makeCluster(threads, type="FORK") #
    registerDoParallel(cl)
    clusterEvalQ(cl , c(library(data.table), library(fundiversity)))
    
    parLapply(cl, 
              n_files, 
              fun = extraction, 
              tile = tiles, 
              unique_dates = unique_dates,
              input_path = input_path, 
              out_path = out_path)  
    
    stopCluster(cl)
    gc()

    Sys.sleep(10)
    
    x <- (j*chuck)
    
  }
  
  #Compile results  ------------------------------------------------------------
  files <- list.files(path = out_path, 
                      pattern = ".txt", 
                      all.files = TRUE,
                      full.names = TRUE, 
                      recursive = TRUE)

  # Loop over files ------------------------------------------------------------
  results <- data.table()
  
  for(i in 1:length(files)) {
    
    stand <- fread(files[i])
    results <- rbind(results, stand)
    
  }
  
  results <- results[order(PLOTID, date)]
  
  fwrite(results, paste0(out_path, "/master_FIA-functional.txt"), sep = "\t")
  
}

functional_diversity(fia, input_path, out_path, chuck = 90, threads = 28)
