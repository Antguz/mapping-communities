#'------------------------------------------------------------------------------
#' @title Get correlations between MDS axes and attributes
#'------------------------------------------------------------------------------

#' @description It derives correlations between MDS axes and attributes associated
#' with climate, topography, plant community traits and abundance of plant
#' lineages.

#' @return .csv files with the correlations between MDS axes and species,
#' lineages, traits, and environment per diversity.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
library(bit64)
library(Rcpp)
library(RcppArmadillo)
library(dynutils)

#'------------------------------------------------------------------------------
#' Source code

source("R/auxiliary/MDS_rotate.R")
source("R/auxiliary/aggregate_plants.R")

#'------------------------------------------------------------------------------
#' Function

axes_correlations <- function(ordination,
                              fia_environment,
                              sp_meta,
                              species,
                              fia_traits,
                              out_path,
                              type) {
  
  # ----------------------------------------------------------------------------
  # Read ordination
  ord <- readRDS(ordination)
  comm_names <- as.integer64(rownames(ord))
  
  # ----------------------------------------------------------------------------
  # Read climate and rotate ordination
  
  # Climate
  env <- fread(fia_environment)
  vector_match <- match(comm_names, env$PLT_CN)
  env <- env[vector_match,]
  
  # Rotate
  ord_new <- MDS_rotate(ord, env[, c("MAT", "MAP")], na.rm = TRUE)
  rm(list = c("ord"))
  gc()
  
  #-----------------------------------------------------------------------------
  # Get community ordinations and goodness of fit
  
  # Get community ordinations
  comm <- ord_new
  colnames(comm) <- paste0("Axis", 1:ncol(comm))
  comm <- cbind(data.table(PLT_CN = as.integer64(rownames(comm))), comm)
  
  # Scale from -0.5 to 0.5
  comm[, 2:4] <- scale_uniform(comm[, 2:4], center = 0, max_range = 1)
  
  # ----------------------------------------------------------------------------
  # Get species scores
  
  # Read species
  sp <- readRDS(species)
  vector_match <- match(comm_names, as.integer64(rownames(sp)))
  sp <- sp[vector_match,]
  sp1 <- decostand(as.matrix(sp), method = "hellinger")
  
  # Initialize result matrices
  sp_cor_matrix <- matrix(NA, ncol = ncol(comm[,-1]), nrow = ncol(sp1))
  sp_pval_matrix <- matrix(NA, ncol = ncol(comm[,-1]), nrow = ncol(sp1))
  
  # Set row/col names for clarity
  rownames(sp_cor_matrix) <- colnames(sp1)
  colnames(sp_cor_matrix) <- colnames(comm[,-1])
  rownames(sp_pval_matrix) <- colnames(sp1)
  colnames(sp_pval_matrix) <- paste0("sig_", colnames(comm[,-1]))
  
  B <- as.matrix(comm[, -1])
  
  # Loop over columns of A and B
  for (i in 1:ncol(sp1)) {
    for (j in 1:ncol(B)) {
      test <- cor.test(sp1[, i], B[, j], method = "pearson")
      sp_cor_matrix[i, j] <- test$estimate
      sp_pval_matrix[i, j] <- test$p.value
    }
  }
  
  sp_cor <- cbind(data.table(plant = rownames(sp_pval_matrix)),
                  sp_cor_matrix,
                  sp_pval_matrix)
  
  fwrite(sp_cor, paste0(out_path, "/MDS_species-cor_", type, "_total.csv"))
  
  # ----------------------------------------------------------------------------
  # Get lineages scores
  
  # Aggregate species
  aggregation <- aggregate_batch(fread(sp_meta)[,c("Taxon2", "Family", "Genus", "growth_form", "mycorrhizal_association")], 
                                 decostand(sp, "total"))
  
  # Initialize result matrices
  linea_cor_matrix <- matrix(NA, ncol = ncol(comm[,-1]), nrow = ncol(aggregation[,-1]))
  linea_pval_matrix <- matrix(NA, ncol = ncol(comm[,-1]), nrow = ncol(aggregation[,-1]))
  
  # Set row/col names for clarity
  rownames(linea_cor_matrix) <- colnames(aggregation[, -1])
  colnames(linea_cor_matrix) <- colnames(comm[,-1])
  rownames(linea_pval_matrix) <- colnames(aggregation[, -1])
  colnames(linea_pval_matrix) <- paste0("sig_", colnames(comm[,-1]))
  
  B <- as.matrix(comm[, -1])
  A <- as.matrix(aggregation[, -1])
  
  # Loop over columns of A and B
  for (i in 1:ncol(A)) {
    for (j in 1:ncol(B)) {
      test <- cor.test(A[, i], B[, j], method = "pearson", na.action = na.exclude())
      linea_cor_matrix[i, j] <- test$estimate
      linea_pval_matrix[i, j] <- test$p.value
    }
  }
  
  linea_cor <- cbind(data.table(plant = rownames(linea_pval_matrix)),
                     linea_cor_matrix,
                     linea_pval_matrix)
  
  fwrite(linea_cor, paste0(out_path, "/MDS_lineages-cor_", type, "_total.csv"))
  
  # ----------------------------------------------------------------------------
  # Get traits scores
  
  # Read traits
  traits <- fread(fia_traits)
  vector_match <- match(colnames(sp), traits$Taxon2)
  traits <- traits[vector_match, ]
  meta <- traits[, c("Taxon2", 
                     "HTq", "SLq", "WD", 
                     "shade", "drought", "waterlogging", 
                     "LMA", "N")]
  communities <- decostand(sp, "total")
  meta <- aggregate_batch(meta, communities)
  
  # Initialize result matrices
  trait_cor_matrix <- matrix(NA, ncol = ncol(comm[,-1]), nrow = ncol(meta[,-1]))
  trait_pval_matrix <- matrix(NA, ncol = ncol(comm[,-1]), nrow = ncol(meta[,-1]))
  
  # Set row/col names for clarity
  rownames(trait_cor_matrix) <- colnames(meta[, -1])
  colnames(trait_cor_matrix) <- colnames(comm[,-1])
  rownames(trait_pval_matrix) <- colnames(meta[, -1])
  colnames(trait_pval_matrix) <- paste0("sig_", colnames(comm[,-1]))
  
  B <- as.matrix(comm[, -1])
  A <- as.matrix(meta[, -1])
  
  # Loop over columns of A and B
  for (i in 1:ncol(A)) {
    for (j in 1:ncol(B)) {
      test <- cor.test(A[, i], B[, j], method = "pearson", na.action = na.exclude())
      trait_cor_matrix[i, j] <- test$estimate
      trait_pval_matrix[i, j] <- test$p.value
    }
  }
  
  trait_cor <- cbind(data.table(plant = rownames(trait_pval_matrix)),
                     trait_cor_matrix,
                     trait_pval_matrix)
  
  fwrite(trait_cor, paste0(out_path, "/MDS_traits-cor_", type, "_total.csv"))
  
  # ----------------------------------------------------------------------------
  # Get environmental scores
  
  # Read traits
  vector_match <- match(env$PLT_CN, comm$PLT_CN)
  comm <- comm[vector_match,]
  
  # Initialize result matrices
  env_cor_matrix <- matrix(NA, ncol = ncol(comm[,-1]), nrow = ncol(env[,-1]))
  env_pval_matrix <- matrix(NA, ncol = ncol(comm[,-1]), nrow = ncol(env[,-1]))
  
  # Set row/col names for clarity
  rownames(env_cor_matrix) <- colnames(env[, -1])
  colnames(env_cor_matrix) <- colnames(comm[,-1])
  rownames(env_pval_matrix) <- colnames(env[, -1])
  colnames(env_pval_matrix) <- paste0("sig_", colnames(comm[,-1]))
  
  B <- as.matrix(comm[, -1])
  A <- as.matrix(env[, -1])
  
  # Loop over columns of A and B
  for (i in 1:ncol(A)) {
    for (j in 1:ncol(B)) {
      test <- cor.test(A[, i], B[, j], method = "pearson", na.action = na.exclude())
      env_cor_matrix[i, j] <- test$estimate
      env_pval_matrix[i, j] <- test$p.value
    }
  }
  
  env_cor <- cbind(data.table(plant = rownames(env_pval_matrix)),
                   env_cor_matrix,
                   env_pval_matrix)
  
  fwrite(env_cor, paste0(out_path, "/MDS_enviromental-cor_", type, "_total.csv"))

}

#'------------------------------------------------------------------------------
#' @example

root_path <- "/media/guzman/Work/MMDB"

# Taxonomic
# ordination <- paste0(root_path, "/Diversity/NMDS/NMDS_taxo_total.rds")
ordination <- paste0(root_path, "/Diversity/02_LMDS/lmds_taxo_VOL_hellinger.rds")
fia_environment <- paste0(root_path, "/FIA/03_FIA-env/fia-env.csv")
sp_meta <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
species <- paste0(root_path, "/Diversity/01_beta/master_matched_comm_phylo_VOL.rds")
fia_traits <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
out_path <- paste0(root_path, "/Diversity/04_MDS-corr")
type <- "taxo"
axes_correlations(ordination, fia_environment, sp_meta, species, fia_traits, out_path, type)

# Phylogenetic
ordination <- paste0(root_path, "/Diversity/02_LMDS/lmds_phylo_VOL_hellinger_norm.rds")
fia_environment <- paste0(root_path, "/FIA/03_FIA-env/fia-env.csv")
sp_meta <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
species <- paste0(root_path, "/Diversity/01_beta/master_matched_comm_phylo_VOL.rds")
fia_traits <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
out_path <- paste0(root_path, "/Diversity/04_MDS-corr")
type <- "phylo"
axes_correlations(ordination, fia_environment, sp_meta, species, fia_traits, out_path, type)

# Functional
ordination <- paste0(root_path, "/Diversity/02_LMDS/lmds_func_VOL_hellinger_norm.rds")
fia_environment <- paste0(root_path, "/FIA/03_FIA-env/fia-env.csv")
sp_meta <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
species <- paste0(root_path, "/Diversity/01_beta/master_matched_comm_func_VOL.rds")
fia_traits <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
out_path <- paste0(root_path, "/Diversity/04_MDS-corr")
type <- "func"
axes_correlations(ordination, fia_environment, sp_meta, species, fia_traits, out_path, type)

