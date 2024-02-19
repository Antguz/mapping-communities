#'------------------------------------------------------------------------------
#' @title Beta taxonomic diversity between communities 
#'------------------------------------------------------------------------------

#' @description It estimates beta taxonomic diversity between communities.
#' @return Return a distance matrices as .rds files.

################################################################################
#' @warming If you have to many communities this need to be run using HPC
################################################################################

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(bit64)

#'------------------------------------------------------------------------------
#' Source code
setwd("/panfs/jay/groups/17/cavender/guzman/DESIS_mapping/mapping-species")
source("R/beta_rcpp.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param out_path
#' @param threads
#' @param progress

#'------------------------------------------------------------------------------
#' Function
taxonomic_diversity <- function(fia, 
                                outpath, 
                                progress,
                                threads) {
  
  ### Get taxonomic beta diversity ---------------------------------------------
  
  # Compute taxonomic beta diversity
  beta <- beta_rcpp(comm = fia, 
                    tree = NULL, 
                    taxonomic = TRUE,
                    threads = threads,
                    progress = progress)
  
  saveRDS(beta, paste0(out_path, "/beta_taxo_VOL.rds"))
  
}

#'------------------------------------------------------------------------------
#' @example 

#root_path <- "E:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
root_path <- "/panfs/jay/groups/17/cavender/guzman/DESIS_mapping"
fia <- readRDS(paste0(root_path, "/03_diversity/master_matched_comm_phylo-VOL.rds"))
out_path <- paste0(root_path, "/03_diversity")
progress <- TRUE
threads <- 120

taxonomic_diversity(fia, 
                    outpath, 
                    progress,
                    threads)
