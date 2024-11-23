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
options(scipen = 10000)

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/goodness_rcpp.cpp")

#'------------------------------------------------------------------------------
#' Function

postprocessing <- function(ordination, fia_environment, out_path, type, text) {
  
  # ----------------------------------------------------------------------------
  # Read ordination
  ord <- readRDS(ordination)
  comm_names <- as.integer64(rownames(ord$points))
  
  # ----------------------------------------------------------------------------
  # Read climate and rotate ordination
  
  # Climate
  env <- fread(fia_environment)
  vector_match <- match(comm_names, env$PLT_CN)
  env <- env[vector_match,]
  
  # Rotate
  ord_new <- MDSrotate(ord, env[, c("Tavg", "Precipitation")], na.rm = TRUE)
  rm(list = c("ord", "comm_names"))
  gc()
  
  #-----------------------------------------------------------------------------
  # Get community ordinations and goodness of fit
  
  # Get community ordinations
  comm <- ord_new$points
  colnames(comm) <- paste0("Axis", 1:ncol(comm))
  comm <- cbind(data.table(PLT_CN = as.integer64(rownames(comm))), comm)
  fwrite(comm, paste0(out_path, "/NMDS_scores_comm_", type, "_total.csv"))
  
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
  
  # Inputs
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
#root_path <- "E:/Projects/DESIS_mapping"

# Taxonomic
ordination <- paste0(root_path, "/04_NMDS/NMDS_taxo_total.rds")
fia_environment <- paste0(root_path, "/00_FIA/plots_climate.csv")
species <- paste0(root_path, "/03_diversity/master_matched_comm_phylo-VOL.rds")
out_path <- paste0(root_path, "/04_NMDS")
type <- "taxo"
text <- "Taxonomic"
postprocessing(ordination, fia_environment, out_path, type, text)

# Phylogenetic
ordination <- paste0(root_path, "/04_NMDS/NMDS_phylo_total.rds")
fia_environment <- paste0(root_path, "/00_FIA/plots_climate.csv")
species <- paste0(root_path, "/03_diversity/master_matched_comm_phylo-VOL.rds")
out_path <- paste0(root_path, "/04_NMDS")
type <- "phylo"
text <- "Phylogenetic"
postprocessing(ordination, fia_environment, out_path, type, text)

# Functional
ordination <- paste0(root_path, "/04_NMDS/NMDS_func_total.rds")
fia_environment <- paste0(root_path, "/00_FIA/plots_climate.csv")
species <- paste0(root_path, "/03_diversity/master_matched_comm_func-VOL.rds")
out_path <- paste0(root_path, "/04_NMDS")
type <- "func"
text <- "Functional"
postprocessing(ordination, fia_environment, out_path, type, text)

