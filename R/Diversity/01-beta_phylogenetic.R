#'------------------------------------------------------------------------------
#' @title Beta phylogenetic diversity between communities 
#'------------------------------------------------------------------------------

#' @description It estimates beta phylogenetic diversity between communities and its
#' partition between replacement and richness. Values are scaled to 10000 to reduce
#' the size of the files and keep significant digits.
 
#' @return Return a list of distance matrices (Btotal, Brepl, Brich) as .rds files

################################################################################
#' @warning If you have too many communities, you need to run it using HPC
################################################################################

#'------------------------------------------------------------------------------
#' @Library

library(data.table)
library(ape)
library(picante)
library(V.PhyloMaker2)
library(bit64)

#'------------------------------------------------------------------------------
#' @Source-code

#setwd("/panfs/jay/groups/17/cavender/guzman/MMDB")
setwd("/n/netscratch/jcavender_lab/Everyone/guzman/MMDB")
source("R/auxiliary/beta_rcpp.R")
source("R/auxiliary/normalize_tree.R")

#'------------------------------------------------------------------------------
#' @Arguments

#' @param fia A data.table with the forest communities (from 01_FIA-cleaning.R)
#' @param SPCD A data.table of functional traits per species
#' @param out_path A path to export estimations of beta diversity
#' @param progress Logic, provide progress bar
#' @param threads Number of threads for parallel processing

#'------------------------------------------------------------------------------
#' @Function

phylogenetic_diversity <- function(fia, 
                                   SPCD, 
                                   outpath, 
                                   progress,
                                   threads) {
  
  ###---------------------------------------------------------------------------
  # Merge inventory with taxonomic information
  fia_summary <- merge(fia, SPCD[, c("SPCD", "Taxon2", "Genus", "Family")], 
                       by = "SPCD", all.x = TRUE, all.y = FALSE)
  
  ### Prepare phylo  -----------------------------------------------------------
  # Summary of species and families
  species <- unique(fia_summary[, c("Taxon2", "Genus", "Family")])
  colnames(species) <- c("species", "genus", "family")
  species <- na.exclude(species)
  hypotheses <- phylo.maker(species, scenarios = "S3")
  
  # Resolve Multichotomies
  phylo <- multi2di(hypotheses$scenario.3)
  is.binary.phylo(phylo) #Test for Binary Tree
  is.ultrametric(phylo) #Test if a Tree is Ultrametric
  
  # Normalize tree
  phylo <- normalize_tree(phylo)
  
  #plot(phylo, type = "fan")

  ### Community matrix ---------------------------------------------------------
  
  # dcast
  master_matrix <- sample2matrix(fia_summary[, c("PLT_CN", "VOL_TPA", "Taxon2")])
  
  ### Get phylogenetic beta diversity ------------------------------------------
  
  #Get match
  matched <- picante::match.phylo.comm(phy = phylo, comm = master_matrix)
  cat(paste0("Number of communities:", nrow(master_matrix)))
  
  # Export
  saveRDS(matched$comm, paste0(out_path, "/master_matched_comm_phylo_VOL.rds"))
  saveRDS(matched$phy, paste0(out_path, "/master_matched_phylo.rds"))
  
  # Compute functional beta diversity
  beta <- beta_rcpp(comm = as.matrix(matched$comm), 
                    tree = matched$phy, 
                    taxonomic = FALSE,
                    threads = threads,
                    progress = progress)
  
  saveRDS(beta, paste0(out_path, "/beta_phylo_VOL.rds"))
  
}

#'------------------------------------------------------------------------------
#' @example 

#root_path <- "/media/antonio/Work/MMDB"
# root_path <- "/panfs/jay/groups/17/cavender/guzman/MMDB"
root_path <- "/n/netscratch/jcavender_lab/Everyone/guzman/MMDB"

fia <- readRDS(paste0(root_path, "/FIA/01_FIA-reshape/fia_reshape.rds"))
SPCD <- fread(paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv"))
out_path <- paste0(root_path, "/Diversity/beta")
progress <- TRUE
threads <- 100

phylogenetic_diversity(fia, 
                       SPCD, 
                       outpath, 
                       progress,
                       threads)
