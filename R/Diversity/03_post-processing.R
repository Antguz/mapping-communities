#'------------------------------------------------------------------------------
#' @title A script for post processing the NMDS ordinations
#'------------------------------------------------------------------------------

#' @description It rotates ordinations according to environmental gradients to 
#' then extract information associated with the projection of communities and 
#' species. It also extract the goodness of the communities in the ordinations 
#' and the data for the stressplot.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
library(bit64)
library(Rcpp)
library(RcppArmadillo)
library(rTLS)
library(dynutils)

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/goodness_rcpp.cpp")
source("R/auxiliary/aggregate_plants.R")

#'------------------------------------------------------------------------------
#' Function

postprocessing <- function(ordination, fia_environment, out_path, type, text) {
  
  # ----------------------------------------------------------------------------
  # Read ordination
  ord <- readRDS(ordination)
  comm_names <- as.integer64(rownames(ord))
  # comm_names <- as.integer64(rownames(ord$points))
  
  # ----------------------------------------------------------------------------
  # Read climate and rotate ordination
  
  # Climate
  env <- fread(fia_environment)
  vector_match <- match(comm_names, env$PLT_CN)
  env <- env[vector_match,]
  
  # Rotate
  # ord_new <- MDSrotate(ord, env[, c("MAT", "MAP")], na.rm = TRUE)
  # rm(list = c("ord"))
  # gc()
  
  #-----------------------------------------------------------------------------
  # Get community ordinations and goodness of fit
  
  # Get community ordinations
  # comm <- ord_new$points
  comm <- ord
  colnames(comm) <- paste0("Axis", 1:ncol(comm))
  comm <- cbind(data.table(PLT_CN = as.integer64(rownames(comm))), comm)
  
  # Scale from -0.5 to 0.5
  comm[, 2:4] <- scale_uniform(comm[, 2:4], center = 0, max_range = 1)
  # fwrite(comm, paste0(out_path, "/NMDS_scores_comm_", type, "_total.csv"))
  fwrite(comm, paste0(out_path, "/MDS_scores_comm_", type, "_total.csv"))
  
  # ----------------------------------------------------------------------------
  # Get species scores
  
  # Read species
  sp <- readRDS(species)
  vector_match <- match(comm_names, as.integer64(rownames(sp)))
  sp <- sp[vector_match,]
  
  # Get scores
  sp_scores <- wascores(x = comm[, 2:4], sp, expand = TRUE)
  sp_scores <- data.table(variable = rownames(sp_scores), sp_scores)
  # fwrite(sp_scores, paste0(out_path, "/NMDS_species-wascores_", type, "_total.csv"))
  fwrite(sp_scores, paste0(out_path, "/MDS_species-wascores_", type, "_total.csv"))
  rm(list = c("sp_scores"))
  
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
  traits_scores <- wascores(x = comm[, 2:4], meta[,-1], expand = TRUE)
  traits_scores <- data.table(variable = rownames(traits_scores), traits_scores)
  # fwrite(traits_scores, paste0(out_path, "/NMDS_traits-wascores_", type, "_total.csv"))
  fwrite(traits_scores, paste0(out_path, "/MDS_traits-wascores_", type, "_total.csv"))
  rm(list = c("traits_scores", "meta"))
  
  # ----------------------------------------------------------------------------
  # Get environmental scores
  
  # Read traits
  env <- na.exclude(env)
  vector_match <- match(env$PLT_CN, comm$PLT_CN)
  comm2 <- comm[vector_match,]
  env_scores <- wascores(x = comm2[, 2:4], env[,-1], expand = TRUE)
  env_scores <- data.table(variable = rownames(env_scores), env_scores)
  # fwrite(env_scores, paste0(out_path, "/NMDS_env-wascores_", type, "_total.csv"))
  fwrite(env_scores, paste0(out_path, "/MDS_env-wascores_", type, "_total.csv"))
  rm(list = c("env_scores", "comm2"))
  
  # ----------------------------------------------------------------------------
  # Get goodness
  
  gof <- goodness_rcpp(dist = ord_new$dist,
                      dhat = ord_new$dhat,
                      nobj = ord_new$nobj,
                      ndis = ord_new$ndis,
                      iidx = ord_new$iidx,
                      jidx = ord_new$jidx,
                      threads = 20)
  
  comm <- cbind(comm, gof = gof[,1])
  fwrite(comm, paste0(out_path, "/NMDS_scores_comm_", type, "_total.csv"))
  rm(list = c("comm", "gof", "env", "vector_match"))
  gc()

  #-----------------------------------------------------------------------------
  # Data for stressplot

  Inputs
  istart <- ord_new$istart
  ndis <- ord_new$ndis
  iregn <- ord_new$iregn
  ngrp <- ord_new$ngrp
  stress <- ord_new$stress
  dist <- ord_new$dist
  diss <- ord_new$diss
  dhat <- ord_new$dhat
  rm(ord_new)
  gc()
  
  # Collect values for 'linear fit'
  ralscal <- 0
  k <- NA
  
  # Monotone line
  if (iregn != 2) {
   ist <- c(istart, ndis + 1)
  
   for(j in 1:ngrp) {
     k <- seq(ist[j], ist[j+1]-1)
     ralscal <- ralscal + cor(dist[k], dhat[k])^2
   }
  }
  
  # Line
  b <- seq(from = 1, to = ndis, length.out = 10000)
  x_line <- diss[k][b]
  y_line <- dhat[k][b]
  
  rm(list = c("k", "b", "dhat"))
  gc()
  
  # Stress as R2
  rstress <- 1 - stress^2
  Rst <- format(rstress, digits = 2)
  Ral <- format(ralscal, digits = 2)
  
  # Vox
  frame <- data.table(x = diss, y = dist, z = 1)
  rm(list = c("diss", "dist"))
  gc()
  vox <- voxels(frame, edge_length = c(0.001, 0.001, 0.001), threads = 16, obj.voxels = FALSE)
  rm(list = c("frame"))
  gc()
  
  # List to export
  to_export <- list(max = max(vox$N),
                   vox = vox,
                   line = data.table(x_line = x_line,
                                     y_line = y_line),
                   rstress = rstress,
                   Rst = Rst,
                   Ral = Ral)
  
  saveRDS(to_export, paste0(out_path, "/data_stressplot_", type, ".rds"))

}

#'------------------------------------------------------------------------------
#' @example

root_path <- "/media/antonio/Work/MMDB"
#root_path <- "/panfs/jay/groups/17/cavender/guzman/DESIS_mapping"
# root_path <- "E:/backup/MMDB"

# Taxonomic
# ordination <- paste0(root_path, "/Diversity/NMDS/NMDS_taxo_total.rds")
ordination <- paste0(root_path, "/Diversity/02_LMDS/lmds_taxo_VOL.rds")
fia_environment <- paste0(root_path, "/FIA/03_FIA-env/fia-env.csv")
species <- paste0(root_path, "/Diversity/01_beta/master_matched_comm_phylo_VOL.rds")
fia_traits <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
out_path <- paste0(root_path, "/Diversity/03_MDS-post")
type <- "taxo"
text <- "Taxonomic"
postprocessing(ordination, fia_environment, out_path, type, text)

# Phylogenetic
ordination <- paste0(root_path, "/Diversity/02_NMDS/NMDS_phylo_total.rds")
fia_environment <- paste0(root_path, "/FIA/03_FIA-env/fia-env.csv")
species <- paste0(root_path, "/Diversity/beta/master_matched_comm_phylo_VOL.rds")
fia_traits <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
out_path <- paste0(root_path, "/Diversity/NMDS")
type <- "phylo"
text <- "Phylogenetic"
postprocessing(ordination, fia_environment, out_path, type, text)

# Functional
ordination <- paste0(root_path, "/Diversity/02_NMDS/NMDS_func_total.rds")
fia_environment <- paste0(root_path, "/FIA/03_FIA-env/fia-env.csv")
species <- paste0(root_path, "/Diversity/01_beta/master_matched_comm_func_VOL.rds")
fia_traits <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
out_path <- paste0(root_path, "/Diversity/02_NMDS")
type <- "func"
text <- "Functional"
postprocessing(ordination, fia_environment, out_path, type, text)

