#'------------------------------------------------------------------------------
#' @title Beta functional diversity between communities 
#'------------------------------------------------------------------------------

#' @description It estimates beta functional diversity between communities
#' @return Return a distance matrices as .rds files

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
source("R/MPDaw_rcpp.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param SPCD
#' @param out_path
#' @param threads
#' @param progress

#'------------------------------------------------------------------------------
#' Function
functional_diversity <- function(fia, 
                                 SPCD, 
                                 outpath, 
                                 threads) {
  
  ###---------------------------------------------------------------------------
  # Merge inventory with taxonomic information
  fia_summary <- merge(fia, SPCD[, c("SPCD", "Taxon2")], 
                       by = "SPCD", all.x = TRUE, all.y = FALSE)
  
  ### Prepare functional--------------------------------------------------------
  # Summary of species and families
  species <- fia[, .N, by = "SPCD"]
  species <- species[, 1]
  species <- merge(species, 
                   SPCD[, c("SPCD", "Taxon2", "lh", "HTq", "SLq", 
                            "CCLCD_mode", "WOOD_SPGR_GREENVOL_DRYWT", 
                            "MC_PCT_GREEN_WOOD", "CARBON_RATIO_LIVE", "shade", 
                            "drought", "waterlogging", "leaf_area", "LMA", "N")],
                   by = "SPCD",
                   all.x = TRUE,
                   all.y = FALSE)
  
  trait_distance <- dist(species[, c(3:15)], method="euclidean", 
                         diag = TRUE, upper = TRUE)
  
  # Rescale
  trait_distance <- (trait_distance-min(trait_distance))/
                    (max(trait_distance)-min(trait_distance))
  
  # Tree
  func <- hclust(trait_distance, "ave")
  func <- ape::as.phylo(func)
  func$tip.label <- species$Taxon2

  # Transform funct
  tol = 1e-9
  func$edge.length[func$edge.length <= 0] <- tol
  is.ultrametric(func)
  func$edge.length <- (log10(func$edge.length)+(min(log10(func$edge.length))*-1))
  
  #plot(func, type = "fan")
  
  ### Community matrix ---------------------------------------------------------
  
  # dcast
  master_matrix <- sample2matrix(fia_summary[, c("PLT_CN", "VOL_TPA", "Taxon2")])
  
  # Standardization of communities
  #master_matrix <- decostand(master_matrix, method = "total")
  #master_matrix <- round(master_matrix, 4)
  
  ### Get functional beta diversity ------------------------------------------
  
  #Get match
  matched <- picante::match.phylo.comm(phy = func, comm = master_matrix)
  
  cat(paste0("Number of communities:", nrow(master_matrix)))
  
  # Export
  saveRDS(matched$comm, paste0(out_path, "/master_matched_comm_func-VOL.rds"))
  saveRDS(matched$phy, paste0(out_path, "/master_matched_func-tree.rds"))
  
  # Compute functional beta diversity
  PDpw <- comdist(comm = matched$comm, 
                  tree = matched$phy, 
                  threads = threads)
  
  saveRDS(PDpw, paste0(out_path, "/beta_functional_VOL.rds"))
  
}

#'------------------------------------------------------------------------------
#' @example 

#root_path <- "E:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "/panfs/jay/groups/17/cavender/guzman/DESIS_mapping"
fia <- readRDS(paste0(root_path, "/03_reshape/fia_reshape.rds"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCP_meta.csv"), na.strings="NA")
out_path <- paste0(root_path, "/03_diversity")
threads <- 120

functional_diversity(fia, 
                     SPCD, 
                     outpath, 
                     threads)
