#'------------------------------------------------------------------------------
#' @title Non-metric Multidimensional Scaling
#'------------------------------------------------------------------------------

#' @description It creates a NMDS ordination for multiple dimensions of beta 
#' diversity

#' @return .rds file with scores per communities

################################################################################
#' @warming If you have to many communities this need to be run using HPC
################################################################################

#'------------------------------------------------------------------------------
#' @Library

library(data.table)
library(vegan)
library(bit64)

#'------------------------------------------------------------------------------
#' @Arguments

#' @param beta A path to an .rds file with the estimations of beta diversity.
#' @param type A character of the dimension of diversity to export.
#' @param diversity A character describing the metric of beta diversity to use.
#' @param out_path A path to a folder to export the NMDS result.
#' @param threads A number of threads to use.

#'------------------------------------------------------------------------------
#' @Function

ordination <- function(beta, type, diversity, out_path, threads = 1) {
  
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
  gc()
  
  # Run NMDS
  NMDS <- metaMDS(comm = as.dist(Beta),
                  k = 3,
                  try = 0, #15
                  trymax = 0, #30
                  maxit = 0, #50
                  engine = "monoMDS",
                  autotransform = FALSE,
                  plot = FALSE,
                  parallel = threads,
                  wascores = FALSE,
                  trace = 3,
                  model = "global",
                  center = FALSE,
                  pc = FALSE,
                  weakties = FALSE,
                  halfchange = TRUE,
                  threshold = 1.0,
                  scaling = FALSE)
  
  # Ordination export
  saveRDS(NMDS, paste0(out_path, "/NMDS_", type, "_", diversity, ".rds"))
  
}

#'------------------------------------------------------------------------------
#' @example 

#root_path <- "/media/antonio/Work/MMDB"
root_path <- "/panfs/jay/groups/17/cavender/guzman/MMDB"

# Phylogenetic
beta <- paste0(root_path, "/Diversity/beta/beta_phylo_VOL.rds")
out_path <- paste0(root_path, "/Diversity/NMDS")
type <- "phylo"
threads <- 1

diversity <- "total"
if(!file.exists(paste0(out_path, "/NMDS_", type, "_", diversity, ".rds"))) {
  ordination(beta, type, diversity, out_path, threads)
} 

# diversity <- "replacement"
# if(!file.exists(paste0(out_path, "/NMDS_", type, "_", diversity, ".rds"))) {
#   ordination(beta, type, diversity, out_path, threads)
# } 
# 
# diversity <- "richness"
# if(!file.exists(paste0(out_path, "/NMDS_", type, "_", diversity, ".rds"))) {
#   ordination(beta, type, diversity, out_path, threads)
# }
