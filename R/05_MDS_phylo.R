#'------------------------------------------------------------------------------
#' @title Classical multidimensional scaling (PCoA) on beta diversity
#'------------------------------------------------------------------------------

#' @description It perform a MDS ordination from taxonomic, phylogenetic and 
#' functional beta diversity
#' @return .rds file with scores per communities and species

################################################################################
#' @warming If you have to many communities this need to be run using HPC
################################################################################

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(hrms)
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
ordination <- function(beta, 
                       stand, 
                       out_path, 
                       type) {
  
  # Read stands
  comm <- readRDS(stand)
  
  # Read distances
  beta <- readRDS(beta)
  
  ### --------------------------------------------------------------------------
  # Btotal
  
  # Get beta 
  beta <- as.matrix(beta$Btotal/10000)
  gc()
  
  # Run MDS
  MDS <- hrms::pcoa(beta, k = 50)
  
  # Create name to export
  out_name_ordination <- paste0(paste0(out_path, "/PCoA_Btotal_", type,".rds"))
  out_name_species <- paste0(paste0(out_path, "/species_Btotal_", type,".rds"))
  out_name_species_exp <- paste0(paste0(out_path, "/species-exp_Btotal_", type,".rds"))
  
  # Ordination export
  saveRDS(MDS, out_name_ordination)
  
  # Get weighted averages scores of species
  species_scores <- vegan::wascores(x = MDS$pc,
                                    w = comm,
                                    expand = FALSE)
  saveRDS(species_scores, out_name_species)
  
  species_scores_expanded <- vegan::wascores(x = MDS$pc,
                                             w = comm,
                                             expand = TRUE)
  saveRDS(species_scores_expanded, out_name_species_exp)
  
  rm(list = c("beta", 
              "MDS", 
              "out_name_ordination", 
              "out_name_species",
              "out_name_species_exp", 
              "species_scores",
              "species_scores_expanded"))
  gc()
  
  ### --------------------------------------------------------------------------
  # Brepl
  
  # Get beta 
  #beta <- as.matrix(distance$Brepl/10000)
  #gc()
  
  # Run MDS
  #MDS <- cmdscale(beta, 
  #                k = 3, eig = TRUE, add = TRUE)
  
  # Create name to export
  #out_name_ordination <- paste0(paste0(out_path, "/MDS_Brepl_", type,".rds"))
  #out_name_species <- paste0(paste0(out_path, "/species_Brepl_", type,".rds"))
  
  # Ordination export
  #saveRDS(MDS, out_name_ordination)
  
  # Get weighted averages scores of species
  #species_scores <- wascores(MDS$points, comm, expand = FALSE)
  #species_names <- rownames(species_scores)
  #species_scores <- cbind(Taxon2 = species_names, as.data.table(species_scores))
  #colnames(species_scores) <- c("Taxon2", "pco1", "pco2", "pco3")
  #saveRDS(species_scores, out_name_species)
  
  #rm(list = c("beta", "MDS", "out_name_ordination", "out_name_species",
  #            "species_scores", "species_names"))
  #gc()
  
  ### --------------------------------------------------------------------------
  # Brich
  
  # Get beta 
  #beta <- as.matrix(distance$Brich/10000)
  #gc()
  
  # Run MDS
  #MDS <- cmdscale(beta, 
  #                k = 3, eig = TRUE, add = TRUE)
  
  # Create name to export
  #out_name_ordination <- paste0(paste0(out_path, "/MDS_Brich_", type,".rds"))
  #out_name_species <- paste0(paste0(out_path, "/species_Brich_", type,".rds"))
  
  # Ordination export
  #saveRDS(NMDS, out_name_ordination)
  
  # Get weighted averages scores of species
  #species_scores <- wascores(MDS$points, comm, expand = FALSE)
  #species_names <- rownames(species_scores)
  #species_scores <- cbind(Taxon2 = species_names, as.data.table(species_scores))
  #colnames(species_scores) <- c("Taxon2", "pco1", "pco2", "pco3")
  #saveRDS(species_scores, out_name_species)
  
  #rm(list = c("beta", "MDS", "out_name_ordination", "out_name_species",
  #            "species_scores", "species_names"))
  #gc()
}

#'------------------------------------------------------------------------------
#' @example 

#root_path <- "F:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "/panfs/jay/groups/17/cavender/guzman/DESIS_mapping"

# Phylogenetic
beta <- paste0(root_path, "/03_diversity/beta_phylo_VOL.rds")
stand <- paste0(root_path, "/03_diversity/master_matched_comm_phylo-VOL.rds")
out_path <- paste0(root_path, "/05_MDS")
type <- "phylo"
ordination(beta, stand, out_path, type)

