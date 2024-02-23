#'------------------------------------------------------------------------------
#' @title Classical multidimensional scaling (PCoA) on beta diversity
#'------------------------------------------------------------------------------

#' @description It perform a MDS ordination from taxonomic, phylogenetic and 
#' functional beta diversity
#' @return .rds file with scores per communities and species

################################################################################
#' @warming If you have to many communities this need to be run using HPC
################################################################################

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(stats)
library(bit64)

#'------------------------------------------------------------------------------
#' Arguments
#' @param beta
#' @param stand
#' @param out_path
#' @param threads

#'------------------------------------------------------------------------------
#' Function
ordination <- function(beta, 
                       out_path, 
                       type) {

  # Read distances
  beta <- readRDS(beta)
  
  ### --------------------------------------------------------------------------
  # Btotal
  
  # Get beta 
  beta <- as.dist(beta$Btotal/10000)
  gc()
  
  # Run MDS
  MDS <- cmdscale(beta, 
                  k = 10, 
                  eig = TRUE, 
                  add = FALSE, 
                  x.ret = FALSE)
  
  # Ordination export
  saveRDS(MDS, paste0(paste0(out_path, "/PCoA_Btotal_", type,".rds")))
  
}

#'------------------------------------------------------------------------------
#' @example 

#root_path <- "F:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "/panfs/jay/groups/17/cavender/guzman/DESIS_mapping"

# Taxonomic
beta <- paste0(root_path, "/03_diversity/beta_taxo_VOL.rds")
out_path <- paste0(root_path, "/05_MDS")
type <- "taxo"
ordination(beta, out_path, type)

# Phylogenetic
beta <- paste0(root_path, "/03_diversity/beta_phylo_VOL.rds")
out_path <- paste0(root_path, "/05_MDS")
type <- "phylo"
ordination(beta, out_path, type)

# Functional
beta <- paste0(root_path, "/03_diversity/beta_func_VOL.rds")
out_path <- paste0(root_path, "/05_MDS")
type <- "func"
ordination(beta, out_path, type)
