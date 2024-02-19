#'------------------------------------------------------------------------------
#' @title Beta phylogenetic diversity between communities 
#'------------------------------------------------------------------------------

#' @description It estimates beta phylogenetic diversity between communities.
#' @return Return a distance matrices as .rds files.

################################################################################
#' @warming If you have to many communities this need to be run using HPC
################################################################################

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(phytools)
library(V.PhyloMaker2)
library(picante)
library(bit64)

#'------------------------------------------------------------------------------
#' Source code
setwd("/panfs/jay/groups/17/cavender/guzman/DESIS_mapping/mapping-species")
source("R/beta_rcpp.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param SPCD
#' @param out_path
#' @param threads
#' @param progress

#'------------------------------------------------------------------------------
#' Function
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
  
  # Transform phylo
  tol = 1e-9
  phylo$edge.length[phylo$edge.length <= 0] <- tol
  is.ultrametric(phylo)
  phylo$edge.length <- (log10(phylo$edge.length)+(min(log10(phylo$edge.length))*-1))
  
  #plot(phylo, type = "fan")
  
  ### Community matrix ---------------------------------------------------------
  
  # dcast
  master_matrix <- sample2matrix(fia_summary[, c("PLT_CN", "VOL_TPA", "Taxon2")])
  
  # Standardization of communities
  master_matrix <- decostand(master_matrix, method = "total")
  master_matrix <- round(master_matrix, 6)
  
  ### Get phylogenetic beta diversity ------------------------------------------
  
  #Get match
  matched <- picante::match.phylo.comm(phy = phylo, comm = master_matrix)
  
  cat(paste0("Number of communities:", nrow(master_matrix)))
  
  # Export
  saveRDS(matched$comm, paste0(out_path, "/master_matched_comm_phylo-VOL.rds"))
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

#root_path <- "E:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
root_path <- "/panfs/jay/groups/17/cavender/guzman/DESIS_mapping"
fia <- readRDS(paste0(root_path, "/03_reshape/fia_reshape.rds"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD_meta.csv"), na.strings="NA")
out_path <- paste0(root_path, "/03_diversity")
progress <- TRUE
threads <- 120

phylogenetic_diversity(fia, 
                       SPCD, 
                       outpath, 
                       progress,
                       threads)
