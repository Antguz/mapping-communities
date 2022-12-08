#'------------------------------------------------------------------------------
#' @tile Compile pixel extraction
#'------------------------------------------------------------------------------

#' @description Compile the pixel extraction from previews function
#' @return Return a .csv.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param in_path
#' @param out_path

#'------------------------------------------------------------------------------
#' @example 

root_path <- "F:/DESIS"
root_path <- "/media/antonio/Extreme SSD/DESIS"

in_path <- paste0(root_path, "/01_pixel-subtraction")
bands <- fread(paste0(root_path, "/00_imagery/DESIS_bands.csv"))
out_name_path <- paste0(root_path, "/02_compile-files/02-spectra.csv")

compile_files(in_path, out_path)

#'------------------------------------------------------------------------------
#' Function
compile_files <- function(in_path, bands, out_path) {
  
  #List files
  files <- list.files(path = in_path, 
                      pattern = ".txt", 
                      recursive = TRUE,
                      include.dirs = TRUE)
  
  #Arrange path in frame
  frame <- data.table(matrix(unlist(strsplit(files, "/")), 
                             nrow= length(files), 
                             byrow=TRUE), stringsAsFactors=FALSE)
  colnames(frame) <- c("folder", "scene")
  
  #Compile return
  compile <- data.table()
  
  # Loop file read  -----------------------
  for(i in 1:nrow(frame)) {
    
    #read
    file <- fread(paste0(in_path, "/", files[i]), header = TRUE)
    file <- cbind(data.table(folder = frame$folder[i],
                             scene = strsplit(frame$scene[i], ".txt")[[1]][1]), file)
    compile <- rbind(compile, file)
      
  }
  
  compile <- compile[!is.na(ID)]
  
  #Add band names
  colnames(compile)[14:248] <- bands$Wavelength
  
  fwrite(compile, out_name_path)
  
}
