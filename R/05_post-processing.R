#'------------------------------------------------------------------------------
#' @title A script for post processing the NMDS ordinations
#'------------------------------------------------------------------------------

#' @description It rotates ordinations according to environmental gradients to 
#' then extract information associated with the projection of communities and 
#' species. It also extract the goodness of the communities in the ordinations 
#' and a stressplot.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
library(bit64)
library(Rcpp)
library(RcppArmadillo)

#'------------------------------------------------------------------------------
#' Source code

#sourceCpp("src/goodness_rcpp.cpp")
root_path <- "/panfs/jay/groups/17/cavender/guzman/DESIS_mapping"
source(paste0(root_path, "/mapping-species/figures/figure_stressplot.R"))

#'------------------------------------------------------------------------------
#' Function

postprocessing <- function(ordination, fia_environment, out_path, type, text) {
  
  # Read ordination
  ord <- readRDS(ordination)
  comm_names <- as.integer64(rownames(ord$points))
  
  # Read climate
  env <- fread(fia_environment)
  vector_match <- match(comm_names, env$PLT_CN)
  env <- env[vector_match,]
  
  # Rotate ordinations
  ord_new <- MDSrotate(ord, env[, c("Tavg", "Precipitation")], na.rm = TRUE)
  
  # Get community scores
  comm <- ord_new$points
  colnames(comm) <- paste0("Axis", 1:ncol(comm))
  comm <- cbind(data.table(PLT_CN = as.integer64(rownames(comm))), comm)
  fwrite(comm, paste0(out_path, "/NMDS_scores_comm_", type, ".csv"))
  
  # Get goodness
  #gof <- goodness_rcpp(dist = ord_new$dist, 
  #                     dhat = ord_new$dhat, 
  #                     nobj = ord_new$nobj, 
  #                     ndis = ord_new$ndis, 
  #                     iidx = ord_new$iidx, 
  #                     jidx = ord_new$jidx,
  #                     threads = 128)
  #gof <- goodness(ord_new)
  #gof <- data.table(PLT_CN = comm_names, goodness = as.vector(gof))
  #fwrite(gof, paste0(out_path, "/NMDS_goodness_", type, ".csv"))

  # Get stressplot
  stress_figure(ord_new, text, out_path)
  
}

#'------------------------------------------------------------------------------
#' @example

#root_path <- "/media/antonio/Work/MMDB"
root_path <- "/panfs/jay/groups/17/cavender/guzman/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"

# Taxonomic
ordination <- paste0(root_path, "/04_NMDS/NMDS_taxo.rds")
fia_environment <- paste0(root_path, "/00_FIA/plots_climate.csv")
species <- paste0(root_path, "/03_diversity/master_matched_comm_phylo-VOL.rds")
out_path <- paste0(root_path, "/05_NMDS")
type <- "taxo"
text <- "Taxonomic"
postprocessing(ordination, fia_environment, out_path, type, text)

# Phylogenetic
ordination <- paste0(root_path, "/04_NMDS/NMDS_phylo-log.rds")
fia_environment <- paste0(root_path, "/00_FIA/plots_climate.csv")
species <- paste0(root_path, "/03_diversity/master_matched_comm_phylo-log-VOL.rds")
out_path <- paste0(root_path, "/05_NMDS")
type <- "phylo-log"
text <- "Phylogenetic"
postprocessing(ordination, fia_environment, out_path, type, text)

# Functional
ordination <- paste0(root_path, "/05_NMDS/NMDS_func.rds")
fia_environment <- paste0(root_path, "/00_FIA/plots_climate.csv")
species <- paste0(root_path, "/03_diversity/master_matched_comm_func-VOL.rds")
out_path <- paste0(root_path, "/05_NMDS")
type <- "func"
text <- "Functional"
postprocessing(ordination, fia_environment, out_path, type, text)

