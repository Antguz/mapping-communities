#'------------------------------------------------------------------------------
#' @title Beta taxonomic diversity between communities 
#'------------------------------------------------------------------------------

#' @description It estimates beta taxonomic diversity between communities and its
#' partition between replacement and richness. Values are scaled to 10000 to reduce
#' the side of the files and keep significant digits.
#' 
#' @return Return a list of distance matrices (Btotal, Brepl, Brich) as .rds files

################################################################################
#' @warming If you have to many communities, you need to run it using using HPC
################################################################################

#'------------------------------------------------------------------------------
#' @Library

library(data.table)
library(bit64)

#'------------------------------------------------------------------------------
#' @Source-code

setwd("/panfs/jay/groups/17/cavender/guzman/MMDB")
source("R/beta_rcpp.R")

#'------------------------------------------------------------------------------
#' Arguments

#' @param fia A data.table with the forest communities (from 01_FIA-cleaning.R)
#' @param out_path A path to export estimations of beta diversity
#' @param progress Logic, provide progress bar
#' @param threads Number of threads for parallel processing

#'------------------------------------------------------------------------------
#' @Function

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

#root_path <- "/media/antonio/Work/MMDB"
root_path <- "/panfs/jay/groups/17/cavender/guzman/MMDB"

fia <- readRDS(paste0(root_path, "/Diversity/beta/master_matched_comm_phylo_VOL.rds"))
out_path <- paste0(root_path, "/Diversity/beta")
progress <- TRUE
threads <- 120

taxonomic_diversity(fia, 
                    outpath, 
                    progress,
                    threads)
