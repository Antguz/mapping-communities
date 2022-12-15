#'------------------------------------------------------------------------------
#' @tile Compile pixel extraction
#'------------------------------------------------------------------------------

#' @description Compile the files of pixel extraction by scene and spectral 
#' transformation from previews function.

#' @return Return a .csv.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param in_path Path of the main folder where all the extracted files are located
#' @param out_path A path and name of a file to export.

#'------------------------------------------------------------------------------
#' @example 

root_path <- "/media/antonio/Extreme SSD/DESIS"
in_path <- paste0(root_path, "/01_pixel-subtraction")
out_path <- paste0(root_path, "/02_compile-files/02-spectra")

compile_files(in_path, out_path)

#'------------------------------------------------------------------------------
#' Function
compile_files <- function(in_path, out_path) {
  
  # List files
  files <- list.files(path = in_path, 
                      pattern = ".txt", 
                      recursive = TRUE,
                      include.dirs = TRUE)
  
  # Arrange path in a frame
  frame <- data.table(matrix(unlist(strsplit(files, "/")), 
                             nrow= length(files), 
                             byrow=TRUE), stringsAsFactors=FALSE)
  colnames(frame) <- c("folder", "scene")
  
  # Look for types of transformations
  frame[, type := strsplit(scene, "_")[[1]][3], by = seq_along(1:nrow(frame))]
  unique_type <- unique(frame$type)
  
  # Loop file read per type  -----------------------
  for(j in 1:length(unique_type)) {
    
    # Select type
    type <- subset(frame, type == unique_type[j])
    
    # Compile return
    compile <- data.table()
    
    for(i in 1:nrow(type)) {
      
      # Read
      file <- fread(paste0(in_path, "/", type$folder[i], "/", type$scene[i]), header = TRUE)
      file <- cbind(data.table(folder = type$folder[i],
                               scene = strsplit(type$scene[i], ".txt")[[1]][1]), file)
      compile <- rbind(compile, file)
      
    }
    
    # Name to export and export
    name_export <- paste0(out_path, "_", unique_type[j])
    fwrite(compile, name_export, sep = "\t")
    
  }
}
