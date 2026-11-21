#'------------------------------------------------------------------------------
#' @title A script for getting the goodness-of-fit MDS ordinations
#'------------------------------------------------------------------------------

#' @description The function helps to get the goodnness of fit from MDS ordinations
#' based on the correlation of distances.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
library(bit64)
library(Rcpp)
library(RcppArmadillo)
library(smatr)
library(dynutils)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' Source code

#setwd("/n/netscratch/jcavender_lab/Everyone/guzman/MMDB")
source("R/auxiliary/voxels.R")

#'------------------------------------------------------------------------------
#' Function

goodness_of_fit <- function(ordination_distances, 
                            beta, 
                            threads, 
                            out_path, 
                            type) {
  
  # Read beta
  beta_total <- sqrt(readRDS(beta)$Btotal/10000)
  gc()
  
  # Read distances
  lmds_distance <- readRDS(ordination_distances)
  gc()

  # Plot IDs (one per column of the matrices)
  plot_id <- rownames(lmds_distance)

  #-----------------------------------------------------------------------------
  # Slope estimation (per plot = per column)
  sma_application <- function(X, plot_id, beta_total, lmds_distance) {

    # Subset
    plot <- data.table(beta = beta_total[, X], mds = lmds_distance[, X])

    # Test
    test <- sma(mds ~ beta, data = plot)
    elevation <- test$coef[[1]][1,1]
    slope <- test$coef[[1]][2,1]

    # Results
    return(data.table(PLT_CN = plot_id[X],
                      elevation = elevation,
                      slope = slope,
                      rsquar = test$r2[[1]]))

  }

  sma <- mclapply(X = 1:length(plot_id),
                  FUN = sma_application,
                  plot_id = plot_id,
                  beta_total = beta_total,
                  lmds_distance = lmds_distance,
                  mc.cores = threads)

  sma <- do.call(rbind.data.frame, sma)

  # Make a frame (no PLT_CN column needed for the voxels)
  frame <- data.table(beta = c(beta_total),
                      mds = c(lmds_distance),
                      z = 1)
  rm(list = c("beta_total", "lmds_distance"))
  gc()

  # Correlation and variance explained
  r <- cor(frame$beta, frame$mds)
  rsquared <- r * r
  
  # Degrees of freedom
  freedom = nrow(frame) - 2
  Fvalue <- rsquared /((1 - rsquared)/freedom)
  Pvalue <- pf(Fvalue, 1, freedom, lower.tail=FALSE)
  
  # Return
  goodness_of_fit <- data.table(r = r,
                                rsquared = rsquared,
                                Fvalue = Fvalue,
                                Pvalue = Pvalue)
  
  fwrite(goodness_of_fit, paste0(out_path, "/goodness_of_fit_", type, "_total.csv"))
  gc()

  #-----------------------------------------------------------------------------
  # Data for stressplot

  vox <- voxels(frame, edge_length = c(0.001, 0.001, 0.001), threads = threads, obj.voxels = FALSE)
  rm(list = c("frame"))
  gc()
  
  # List to export
  to_export <- list(max = max(vox$N),
                    vox = vox,
                    linear = sma,
                    goodness_of_fit = goodness_of_fit)
  
  saveRDS(to_export, paste0(out_path, "/data_goodness_of_fit_", type, ".rds"))
  
}

#'------------------------------------------------------------------------------
#' @example

root_path <- "/media/guzman/Work/MMDB"
threads <- 16

# Taxonomic
ordination_distances <- paste0(root_path, "/Diversity/03_MDS-post/lmds_distance_taxo_total.rds")
beta <- paste0(root_path, "/Diversity/01_beta/beta_taxo_VOL_hellinger.rds")
out_path <- paste0(root_path, "/Diversity/03_MDS-post")
type <- "taxo"
goodness_of_fit(ordination_distances, beta, threads, out_path, type)
gc()

# Phylogenetic
ordination_distances <- paste0(root_path, "/Diversity/03_MDS-post/lmds_distance_phylo_total.rds")
beta <- paste0(root_path, "/Diversity/01_beta/beta_phylo_VOL_hellinger_norm.rds")
out_path <- paste0(root_path, "/Diversity/03_MDS-post")
type <- "phylo"
goodness_of_fit(ordination_distances, beta, threads, out_path, type)
gc()

# Functional
ordination_distances <- paste0(root_path, "/Diversity/03_MDS-post/lmds_distance_func_total.rds")
beta <- paste0(root_path, "/Diversity/01_beta/beta_func_VOL_hellinger_norm.rds")
out_path <- paste0(root_path, "/Diversity/03_MDS-post")
type <- "func"
goodness_of_fit(ordination_distances, beta, threads, out_path, type)
gc()
