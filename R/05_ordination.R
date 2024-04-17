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
#' @param out_path
#' @param threads

#'------------------------------------------------------------------------------
#' Function
ordination <- function(beta, stand, type, out_path, threads = 1) {
  
  # Read distances
  distance <- readRDS(beta)
  
  # Select total distances
  distance <- distance$Btotal/10000
  gc()
  
  # Run NMDS
  NMDS <- metaMDS(comm = as.dist(distance),
                  k = 6,
                  try = 10, #15
                  trymax = 15, #30
                  maxit = 50, #50
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
                  maxdist = 1)
  
  # Ordination export
  saveRDS(NMDS, paste0(out_path, "/NMDS_", type, ".rds"))
  
}

#'------------------------------------------------------------------------------
#' @example 

#root_path <- "E:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "/panfs/jay/groups/17/cavender/guzman/DESIS_mapping"

# Taxonomic
beta <- paste0(root_path, "/03_diversity/beta_taxo_VOL.rds")
out_path <- paste0(root_path, "/05_ordination")
type <- "taxo"
threads <- 1
ordination(beta, stand, type, out_path, threads)

# Phylogenetic
beta <- paste0(root_path, "/03_diversity/beta_phylo_VOL.rds")
out_path <- paste0(root_path, "/05_ordination")
type <- "phylo"
threads <- 1
ordination(beta, stand, type, out_path, threads)

# Functional
beta <- paste0(root_path, "/03_diversity/beta_func_VOL.rds")
out_path <- paste0(root_path, "/05_ordination")
type <- "func"
threads <- 1
ordination(beta, stand, out_path, threads)
