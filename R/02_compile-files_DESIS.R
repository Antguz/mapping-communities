#'------------------------------------------------------------------------------
#' @tile Compile pixel extraction files from DESIS
#'------------------------------------------------------------------------------

#' @description Compile the files of pixel extraction by scene and spectral 
#' transformation from previews function.

#' @return Export a frame

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(bit64)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param in_path Path of the main folder where all the extracted files are located
#' @param out_path A path and name of a file to export.

#'------------------------------------------------------------------------------
#' @example 

root_path <- "/media/antonio/Extreme Pro/Imagery" #Main folder (you just need to change this)
folder <- "01_pixel-subtraction"
in_path <- paste0(root_path, "/", folder)
out_path <- paste0(root_path, "/", folder, "/02_spectra") #This will return 3 files

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
  frame[, type := strsplit(type, ".txt")[[1]][1], by = seq_along(1:nrow(frame))]
  unique_type <- unique(frame$type)
  
  # Loop file read per type  -----------------------
  for(j in 1:length(unique_type)) {
    
    # Select type
    type <- subset(frame, type == unique_type[j])
    
    # Compile return
    compile <- data.table()
    
    # Add progress bar
    cat(paste0("Getting ", unique_type[j], " files"))
    pb <- txtProgressBar(min = 1, max = nrow(type), style = 3)
    
    for(i in 1:nrow(type)) {
      
      # Progress bar
      setTxtProgressBar(pb, i)
      
      # Read
      file <- fread(paste0(in_path, "/", type$folder[i], "/", type$scene[i]), header = TRUE)
      file <- cbind(data.table(folder = type$folder[i],
                               scene = strsplit(type$scene[i], ".txt")[[1]][1]), file)
      compile <- rbind(compile, file)
      
    }
    
    # Name to export and export
    name_export <- paste0(out_path, "-", unique_type[j], ".txt")
    fwrite(compile, name_export, sep = "\t")
    
  }
}
