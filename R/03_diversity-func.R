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
library(ape)
library(picante)
library(scales)
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
#' @param progress
#' @param threads

#'------------------------------------------------------------------------------
#' Function
functional_diversity <- function(fia, 
                                 SPCD, 
                                 outpath, 
                                 progress,
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
                   SPCD[, c("SPCD", "Taxon2", "HTq", "SLq", 
                            "CCLCD_mode", "WOOD_SPGR_GREENVOL_DRYWT", 
                            "MC_PCT_GREEN_WOOD", "CARBON_RATIO_LIVE", "shade", 
                            "drought", "waterlogging", "leaf_area", "LMA", "N",
                            "lshape", "ltype", "lhabit")],
                   by = "SPCD",
                   all.x = TRUE,
                   all.y = FALSE)
  
  trait_distance <- dist(scale(species[, c(3:17)]), 
                         method="euclidean", 
                         diag = TRUE, 
                         upper = TRUE)
 
  # Tree
  func <- hclust(trait_distance)
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
  master_matrix <- decostand(master_matrix, method = "total")
  master_matrix <- round(master_matrix, 6)
  
  ### Get functional beta diversity ------------------------------------------
  
  #Get match
  matched <- picante::match.phylo.comm(phy = func, comm = master_matrix)
  
  cat(paste0("Number of communities:", nrow(master_matrix)))
  
  # Export
  saveRDS(matched$comm, paste0(out_path, "/master_matched_comm_func-VOL.rds"))
  saveRDS(matched$phy, paste0(out_path, "/master_matched_func.rds"))
  
  # Compute functional beta diversity
  beta <- beta_rcpp(comm = as.matrix(matched$comm), 
                    tree = matched$phy, 
                    taxonomic = FALSE,
                    threads = threads,
                    progress = progress)
  
  saveRDS(beta, paste0(out_path, "/beta_func_VOL.rds"))
  
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

functional_diversity(fia, 
                     SPCD, 
                     outpath, 
                     progress,
                     threads)

