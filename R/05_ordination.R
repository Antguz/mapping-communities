#'------------------------------------------------------------------------------
#' @title Non-metric Multidimensional Scaling
#'------------------------------------------------------------------------------

#' @description It creates a NMDS ordination of phylogenetic beta diversity
#' @return .rds file with scores per communities and species

################################################################################
#' @warming If you have to many communities this need to be run using HPC
################################################################################

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
library(bit64)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param beta
#' @param stand
#' @param out_path
#' @param threads

#'------------------------------------------------------------------------------
#' Function
ordination <- function(beta, stand, out_path, threads = 28) {
  
  # Read stands
  comm <- readRDS(stand)
  
  # Read distances
  distance <- readRDS(beta)
  
  # Select total distances
  distance <- distance$Btotal/10000
  
  # Run NMDS
  NMDS <- metaMDS(comm = distance,
                  k = 3,
                  try = 15,
                  trymax = 30,
                  maxit = 100,
                  engine = "monoMDS",
                  autotransform = FALSE,
                  plot = FALSE,
                  parallel = threads,
                  wascores = FALSE,
                  trace = 3,
                  model = "linear",
                  center = TRUE,
                  pc = FALSE,
                  weakties = TRUE,
                  stress = 1,
                  halfchange = TRUE,
                  threshold = 1.0)
  
  # Ordination export
  saveRDS(NMDS, paste0(out_path, "/NMDS_DRYBIO.rds"))
  #NMDS <- readRDS(paste0(out_path, "/NMDS_DRYBIO.rds"))
  
  # Get weighted averages scores of species
  species_scores <- wascores(NMDS$points, comm, expand = FALSE)
  species_names <- rownames(species_scores)
  species_scores <- cbind(Taxon2 = species_names, as.data.table(species_scores))
  saveRDS(species_scores, paste0(out_path, "/species_scores_DRYBIO.rds"))
  
  # Get dimensions
  PLT_CN <- as.integer64(rownames(NMDS$points))
  communities_scores <- as.data.table(NMDS$points)
  communities_scores$PLT_CN <- PLT_CN
  communities_scores$goodness <- goodness(NMDS)
  saveRDS(communities_scores, paste0(out_path, "/communities_scores_DRYBIO.rds"))
  
}

#'------------------------------------------------------------------------------
#' @example 

#root_path <- "F:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
root_path <- "/panfs/jay/groups/17/cavender/guzman/DESIS_mapping"
beta <- paste0(root_path, "/03_diversity/beta_DRYBIO.rds")
stand <- paste0(root_path, "/03_diversity/master_matched_comm-DRYBIO.rds")
out_path <- paste0(root_path, "/05_ordination")
threads <- 1
ordination(beta, stand, out_path, threads)
