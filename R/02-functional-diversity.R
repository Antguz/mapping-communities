#'------------------------------------------------------------------------------
#' @title Functional Diversity
#'------------------------------------------------------------------------------

#' @description Estimate the functional diversity of a set of plots
#' @return .txt files with names according the vector file

#'------------------------------------------------------------------------------
#' Library

library(data.table)
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

#'------------------------------------------------------------------------------
#' Function
functional_diversity <- function(input_path, out_path, threads = 24) {
  
  #Search for paths
  files <- list.files(path = input_path, 
                      pattern = ".txt", 
                      all.files = TRUE,
                      full.names = FALSE, 
                      recursive = TRUE)
  
  #Arrange files frame
  frame <- data.table(matrix(unlist(strsplit(files, "-")), 
                             nrow= length(files), 
                             byrow=TRUE), stringsAsFactors=FALSE)
  colnames(frame) <- c("X", "Y", "stand")
  frame[, stand := strsplit(stand, ".txt")[[1]][1], 
        by = seq_along(1:nrow(frame))]
  
  #Unique values
  unique_stands <- unique(frame$stand)
  
  
  # Set up cluster
  cl <- makeCluster(threads, type = "FORK")
  registerDoParallel(cl)
  
  # Loop over tiles ------------------------------------------------------------
  extraction <- foreach(i = 1:length(unique_stands),
                        .packages = c("terra", "data.table"),
                        .inorder = F) %dopar% {
                          
                          sub_stand <- subset(frame, stand == unique_stands[i])
                          
                          plot <- NA
                          
                          #Read stand
                          if(nrow(sub_stand) == 1) {
                            
                            plot <- fread(paste0(input_path, "/",
                                                 sub_stand$X[1], "_",
                                                 sub_stand$Y[1], "_",
                                                 sub_stand$stand[1], ".txt"))
                            
                          } else {
                            
                            for(ii in 1:nrow(sub_stand)) {
                              
                              sub_plot <- fread(paste0(input_path, "/",
                                                       sub_stand$X[1], "_",
                                                       sub_stand$Y[1], "_",
                                                       sub_stand$stand[1], ".txt"))
                              
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
                          mature <- mature[, c(1, 2, 3)]
                          
                          #Remove pixels with low NDVI
                          plot <- merge(mature, plot, by = c("STAND_KEY", "x", "y"), all.x = TRUE, all.y = FALSE)
                          
                          #Unique pixels and dates
                          unique_date <- unique(plot$date)
                          
                          #Loop of metrics over dates
                          for(ii in 1:length(unique_date)) {
                            
                            sub_date <- subset(plot, date == unique_date[ii])
                            
                            
                          }
                        }
}




