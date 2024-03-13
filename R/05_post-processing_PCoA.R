#'------------------------------------------------------------------------------
#' @title A script for post processing the communities ordinations
#'------------------------------------------------------------------------------

#' @description It rotates ordinations according to environmental gradients to 
#' then extract information associated with the projection of communities and 
#' species

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
library(bit64)

#'------------------------------------------------------------------------------
#' Function

postprocessing <- function(ordination, species, type) {
  
  # Read ordination
  ord <- readRDS(ordination)
  comm <- as.integer64(rownames(ord$points))
  
  # Get community scores
  comm <- ord$points
  colnames(comm) <- paste0("PCoA", 1:ncol(comm))
  comm <- cbind(data.table(PLT_CN = as.integer64(rownames(comm))), comm)
  fwrite(comm, paste0(out_path, "/PCoA_scores_comm_", type, ".csv"))
  
  # Read species relative abundance
  spec <- readRDS(species)
  
  # Get species scores
  sp <- wascores(comm[, -1], spec, expand = FALSE)
  sp <- cbind(data.table(Taxon2 = rownames(sp)), sp)
  fwrite(sp, paste0(out_path, "/scores_sp_", type, ".csv"))
  
  sp <- wascores(comm[, -1], spec, expand = TRUE)
  sp <- cbind(data.table(Taxon2 = rownames(sp)), sp)
  fwrite(sp, paste0(out_path, "/scores-exp_sp_", type, ".csv"))
  
}

#'------------------------------------------------------------------------------
#' @example

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"

# Taxonomic
ordination <- paste0(root_path, "/05_MDS/PCoA_Btotal_taxo.rds")
species <- paste0(root_path, "/03_diversity/master_matched_comm_phylo-VOL.rds")
out_path <- paste0(root_path, "/05_MDS")
type <- "taxo"
postprocessing(ordination, species, type)

# Phylogenetic
ordination <- paste0(root_path, "/05_MDS/PCoA_Btotal_phylo-log.rds")
species <- paste0(root_path, "/03_diversity/master_matched_comm_phylo-log-VOL.rds")
out_path <- paste0(root_path, "/05_MDS")
type <- "phylo"
postprocessing(ordination, species, type)

# Phylogenetic
ordination <- paste0(root_path, "/05_MDS/PCoA_Btotal_func.rds")
species <- paste0(root_path, "/03_diversity/master_matched_comm_phylo-VOL.rds")
out_path <- paste0(root_path, "/05_MDS")
type <- "func"
postprocessing(ordination, species, type)
