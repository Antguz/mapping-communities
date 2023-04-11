#'------------------------------------------------------------------------------
#' @title Optimal k for Ordination
#'------------------------------------------------------------------------------

#' @description It estimate the optimal k dimension for the NDMS
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)

#'------------------------------------------------------------------------------
#' Arguments

root_path <- "E:/Projects/DESIS_mapping/05_ordination"
folder <- "phylogenetic"
files_path <- paste0(root_path, "/", folder)
files <- dir(files_path)
#Get stress

stress <- NA

for(i in 1:length(files)) {
  
  name <- paste0(root_path, "/",
                 folder, "/",
                 "NMDS_", folder, "_DRYBIO_k", i, ".rds")
  
  NDS <- readRDS(name)
  
  stress[i] <- NDS$stress
  
}


frame <- data.table(type = "taxonomic", k = 1:25, stress = stress)
framep <- data.table(type = "phylogenetic", k = 1:27, stress = stress)

frame <- rbind(frame, framep)
fwrite(frame, paste0(root_path, "/strees_optimal.csv"))
