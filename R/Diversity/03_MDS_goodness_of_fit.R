#'------------------------------------------------------------------------------
#' @title A script for getting the goodness of fit MDS ordinations
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

source("R/auxiliary/voxels.R")

#'------------------------------------------------------------------------------
#' Function

goodness_of_fit <- function(ordination_distances, 
                            beta, 
                            threads, 
                            out_path, 
                            type) {
  
  # Read beta
  beta_total <- readRDS(beta)$Btotal
  gc()
  
  # Read distances
  lmds_distance <- readRDS(ordination_distances)
  gc()
  
  # Make a frame
  frame <- data.table(PLT_CN = as.integer64(rownames(lmds_distance[1:100, 1:100])),
                      beta = c(beta_total[1:100, 1:100]), 
                      mds = c(lmds_distance[1:100, 1:100]),
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
  rm(list = c("goodness_of_fit"))
  gc()
  
  #-----------------------------------------------------------------------------
  # Slope estimation
  sma_application <- function(X, plot_id, frame) {
    
    # Subset
    plot <- frame[PLT_CN == plot_id[X], ]
    
    # Test
    test <- sma(mds ~ beta, data = plot)
    elevation <- test$coef[[1]][1,1]
    slope <- test$coef[[1]][2,1]
    
    # Results
    return(data.table(PLT_CN = X,
                      elevation = elevation,
                      slope = slope,
                      rsquar = test$r2[[1]]))
    
  }
  
  sma <- mclapply(X = 1:length(unique(frame$PLT_CN)), 
                  FUN = sma_application,
                  plot_id = unique(frame$PLT_CN),
                  frame = frame,
                  mc.cores = 5)
  
  sma <- do.call(rbind.data.frame, sma)
  
  #-----------------------------------------------------------------------------
  # Data for stressplot
  
  vox <- voxels(frame[,-1], edge_length = c(0.001, 0.001, 0.001), threads = threads, obj.voxels = FALSE)
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

root_path <- "/media/antonio/Work/MMDB"
#root_path <- "/n/netscratch/jcavender_lab/Everyone/guzman/MMDB"
threads <- 99

# Taxonomic
# ordination <- paste0(root_path, "/Diversity/NMDS/NMDS_taxo_total.rds")
ordination_distances <- paste0(root_path, "/Diversity/03_MDS-post/lmds_distance_taxo_total.rds")
beta <- paste0(root_path, "/Diversity/01_beta/beta_taxo_VOL_hellinger.rds")
out_path <- paste0(root_path, "/Diversity/03_MDS-post")
type <- "taxo"
goodness_of_fit(ordination_distances, beta, threads, out_path, type)

# Phylogenetic
ordination_distances <- paste0(root_path, "/Diversity/03_MDS-post/lmds_distance_phylo_total.rds")
beta <- paste0(root_path, "/Diversity/01_beta/beta_phylo_VOL_hellinger_norm.rds")
out_path <- paste0(root_path, "/Diversity/03_MDS-post")
type <- "phylo"
goodness_of_fit(ordination_distances, beta, threads, out_path, type)

# Functional
ordination_distances <- paste0(root_path, "/Diversity/03_MDS-post/lmds_distance_func_total.rds")
beta <- paste0(root_path, "/Diversity/01_beta/beta_func_VOL_hellinger_norm.rds")
out_path <- paste0(root_path, "/Diversity/03_MDS-post")
type <- "func"
goodness_of_fit(ordination_distances, beta, threads, out_path, type)

