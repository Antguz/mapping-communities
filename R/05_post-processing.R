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

postprocessing <- function(ordination, fia_environment, species, type) {
  
  # Read ordination
  ord <- readRDS(ordination)
  comm <- as.integer64(rownames(ord$points))
  
  # Read climate
  env <- fread(fia_environment)
  vector_match <- match(comm, env$PLT_CN)
  env <- env[vector_match,]
  
  # Rotate ordinations
  ord_new <- MDSrotate(ord, env[, c("Tavg", "Precipitation")], na.rm = TRUE)
  
  # Get community scores
  comm <- ord_new$points
  comm <- cbind(data.table(PLT_CN = as.integer64(rownames(comm))), comm)
  fwrite(comm, paste0(out_path, "/NDMS_scores_comm_", type, ".csv"))
  
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
ordination <- paste0(root_path, "/05_ordination/NMDS_taxo.rds")
fia_environment <- paste0(root_path, "/00_FIA/plots_climate.csv")
species <- paste0(root_path, "/03_diversity/master_matched_comm_phylo-VOL.rds")
out_path <- paste0(root_path, "/05_ordination")
type <- "taxo"
postprocessing(ordination, fia_environment, species, type)

# Phylogenetic
ordination <- paste0(root_path, "/05_ordination/NMDS_phylo.rds")
fia_environment <- paste0(root_path, "/00_FIA/plots_climate.csv")
species <- paste0(root_path, "/03_diversity/master_matched_comm_phylo-VOL.rds")
out_path <- paste0(root_path, "/05_ordination")
type <- "phylo"
postprocessing(ordination, fia_environment, species, type)
