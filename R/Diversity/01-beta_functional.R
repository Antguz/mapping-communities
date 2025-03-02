#'------------------------------------------------------------------------------
#' @title Beta functional diversity between communities 
#'------------------------------------------------------------------------------

#' @description It estimates beta functional diversity between communities and its
#' partition between replacement and richness. Values are scaled to 10000 to reduce
#' the side of the files and keep significant digits.

#' @return Return a list of distance matrices (Btotal, Brepl, Brich) as .rds files

################################################################################
#' @warming If you have to many communities, you need to run it using using HPC
################################################################################

#'------------------------------------------------------------------------------
#' @Library

library(data.table)
library(ape)
library(picante)
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
                   SPCD[, c("SPCD", "Taxon2", 
                            "HTq", "SLq", "WD", 
                            "shade", "drought", "waterlogging", 
                            "LMA", "N",
                            "lshape", "ltype", "lhabit")],
                   by = "SPCD",
                   all.x = TRUE,
                   all.y = FALSE)
  
  species <- as.data.frame(species)
  species_tips <- species$Taxon2
  
  # # Scale
  # species[, 3:10] <- scale(species[, 3:10])
  # 
  # # Get distance
  # trait_distance <- dist(species[, c(3:10)],
  #                        method= "euclidean",
  #                        diag = TRUE,
  #                        upper = TRUE)
  # 
  # # Tree
  # func <- hclust(trait_distance, "average")
  # func <- ape::as.phylo(func)
  # func$tip.label <- species_tips
  # 
  # # Normalize tree
  # func <- normalize_tree(func)
  
  # Run a PCA
  pca <- rda(species[, c(3:10)], scale = TRUE)
  val <- scores(pca, choices = c(1, 2, 3), display = "sites")
  
  # Get distance
  trait_distance <- dist(val, 
                         method= "euclidean", 
                         diag = TRUE, 
                         upper = TRUE)
  
  # Tree
  func <- hclust(trait_distance, "average")
  func <- ape::as.phylo(func)
  func$tip.label <- species_tips
  
  # Normalize tree
  func <- normalize_tree(func)
  
  ### Community matrix ---------------------------------------------------------
  
  # dcast
  master_matrix <- sample2matrix(fia_summary[, c("PLT_CN", "VOL_TPA", "Taxon2")])
  
  ### Get functional beta diversity ------------------------------------------
  
  #Get match
  matched <- picante::match.phylo.comm(phy = func, comm = master_matrix)
  cat(paste0("Number of communities:", nrow(master_matrix)))
  
  # Export
  saveRDS(matched$comm, paste0(out_path, "/master_matched_comm_func_VOL.rds"))
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

#root_path <- "/media/antonio/Work/MMDB"
#root_path <- "/panfs/jay/groups/17/cavender/guzman/MMDB" MSI
root_path <- "/n/netscratch/jcavender_lab/Everyone/guzman/MMDB"

fia <- readRDS(paste0(root_path, "/FIA/01_FIA-reshape/fia_reshape.rds"))
SPCD <- fread(paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv"))
out_path <- paste0(root_path, "/Diversity/beta")
progress <- TRUE
threads <- 100

functional_diversity(fia, 
                     SPCD, 
                     outpath, 
                     progress,
                     threads)

