#'------------------------------------------------------------------------------
#' @title Results stacking of functional diversity
#'------------------------------------------------------------------------------

#' @description It stack the multiple files of results from functional diversity
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)

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
file_stacking <- function(input_path, path_name) {
  
  #Search for paths
  files <- list.files(path = input_path, 
                      pattern = ".txt", 
                      all.files = TRUE,
                      full.names = TRUE, 
                      recursive = TRUE)
  
  result <- data.table()
  
  for(i in 1:length(files)) {
    
    stand <- fread(files[i])
    
    result <- rbind(result, stand)
    
  }
  
  fwrite(result, path_name, sep = "\t")

}  