#'------------------------------------------------------------------------------
#' @title Non-metric Multidimensional Scaling
#'------------------------------------------------------------------------------

#' @description It creates a NMDS ordination for multiple dimensions of beta 
#' diversity

#' @return .rds file with scores per communities and species

################################################################################
#' @warming If you have to many communities this need to be run using HPC
################################################################################

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
library(bit64)

#'------------------------------------------------------------------------------
#' Arguments
#' @param beta
#' @param stand
#' @param type
#' @param diversity
#' @param out_path
#' @param threads

#'------------------------------------------------------------------------------
#' Function
ordination <- function(beta, stand, type, diversity, out_path, threads = 1) {
  
  # Read distances
  distance <- readRDS(beta)
  
  # Select total distances
  if(diversity == "total") {
    Beta <- distance$Btotal/10000
  } else if(diversity == "replacement") {
    Beta <- distance$Brepl/10000
  } else if(diversity == "richness") {
    Beta <- distance$Brich/10000
  }
  
  rm("distance")
  gc
  
  # Run NMDS
  NMDS <- metaMDS(comm = as.dist(Beta),
                  k = 3,
                  try = 4, #15
                  trymax = 8, #30
                  maxit = 30, #50
                  engine = "monoMDS",
                  autotransform = FALSE,
                  plot = FALSE,
                  parallel = threads,
                  wascores = FALSE,
                  trace = 3,
                  model = "linear",
                  center = FALSE,
                  pc = FALSE,
                  stress = 1,
                  halfchange = FALSE,
                  maxdist = 1,
                  weakties = FALSE,
                  scaling = FALSE)
  
  # Ordination export
  saveRDS(NMDS, paste0(out_path, "/NMDS_", type, "_", diversity, ".rds"))
  
}

#'------------------------------------------------------------------------------
#' @example 

#root_path <- "E:/Projects/MMDB"
#root_path <- "/media/antonio/Work/MMDB"
root_path <- "/panfs/jay/groups/17/cavender/guzman/MMDB"

# Functional
beta <- paste0(root_path, "/03_diversity/beta_functional_VOL.rds")
out_path <- paste0(root_path, "/04_NMDS")
type <- "func"
threads <- 1

diversity <- "total"
if(!file.exists(paste0(out_path, "/NMDS_", type, "_", diversity, ".rds"))) {
  ordination(beta, stand, type, diversity, out_path, threads)
} 

diversity <- "replacement"
if(!file.exists(paste0(out_path, "/NMDS_", type, "_", diversity, ".rds"))) {
  ordination(beta, stand, type, diversity, out_path, threads)
} 

diversity <- "richness"
if(!file.exists(paste0(out_path, "/NMDS_", type, "_", diversity, ".rds"))) {
  ordination(beta, stand, type, diversity, out_path, threads)
}
