#'------------------------------------------------------------------------------
#' @title Results stacking of functional diversity
#'------------------------------------------------------------------------------

#' @description It stack the multiple files of results from functional diversity
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(foreach)
library(doParallel)

#'------------------------------------------------------------------------------
#' Arguments
#' @param input_path
#' @param path_name

#'------------------------------------------------------------------------------
#' @example 

input_path <- "/media/antonio/antonio_ssd/mapping_species/level4"
out_path <- "/media/antonio/antonio_ssd/mapping_species/level5"

input_path <- "E:/mapping_species/level5"
path_name <- "E:/mapping_species/FD.txt"

#'------------------------------------------------------------------------------
#' Function
file_stacking <- function(input_path, path_name, threads = 14) {
  
  #Search for paths
  files <- list.files(path = input_path, 
                      pattern = ".txt", 
                      all.files = TRUE,
                      full.names = TRUE, 
                      recursive = TRUE)
  
  # Set up cluster
  cl <- makeCluster(threads)
  registerDoParallel(cl)
  
  # Loop over files ------------------------------------------------------------
  results <- foreach(i = 1:length(files),
                     .combine = rbind,
                     .packages = c("data.table"),
                     .inorder = F) %dopar% {
                       
                       stand <- fread(files[i])
                       
                       return(stand)
    
                     }
  
  #Stop cluster
  stopCluster(cl)
  gc()
  
  results <- results[order(STAND_KEY)]

  fwrite(results, path_name, sep = "\t")

}  