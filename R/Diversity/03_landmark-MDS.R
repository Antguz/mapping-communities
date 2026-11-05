#'------------------------------------------------------------------------------
#' @title Perform Landmark Multi-Dimensional Scaling
#'------------------------------------------------------------------------------

#' Perform Landmark Multi-Dimensional Scaling on communities previously selected

#'------------------------------------------------------------------------------
#' @Library

library(data.table)
library(bit64)

#'------------------------------------------------------------------------------
#' @Source-code

#source("R/auxiliary/lmds.R")

#'------------------------------------------------------------------------------
#' @Function

comm_lMDS <- function(comm_landmark, beta, output, threads) {
  
  # Get PLT_CN
  PLT_CN <- as.integer64(rownames(beta$Btotal))
  
  # Match landmarks
  match_landmarks <- match(comm_landmark, PLT_CN)
  
  # Get landmarks
  beta_landmarks <- select_landmarks(x = sqrt(beta$Btotal/10000),
                                     landmarks_vector = 1:nrow(beta$Btotal),
                                     threads = threads)
  # Get MDS
  beta_lmds <- cmdscale_landmarks(beta_landmarks, ndim = 3, rescale = TRUE)
  
  # Export MDS
  saveRDS(beta_lmds, output)
  
}

#'------------------------------------------------------------------------------
#' @Application

# Define root path
#root_path <- "/media/antonio/Work/MMDB"
root_path <- "/n/netscratch/jcavender_lab/Everyone/guzman/MMDB"
source("/n/netscratch/jcavender_lab/Everyone/guzman/MMDB/R/auxiliary/lmds.R")
comm_landmark <- readRDS(paste0(root_path, "/Diversity/01_beta/comm_landmarks.rds"))
threads <- 100

# Taxo
beta <- readRDS(paste0(root_path, "/Diversity/01_beta/beta_taxo_VOL_hellinger.rds"))
output <- paste0(root_path, "/Diversity/02_LMDS/lmds_taxo_VOL_hellinger.rds")
comm_lMDS(comm_landmark, beta, output, threads)
rm(list = c('beta', 'output'))

# Phylo
beta <- readRDS(paste0(root_path, "/Diversity/01_beta/beta_phylo_VOL_hellinger_norm.rds"))
output <- paste0(root_path, "/Diversity/02_LMDS/lmds_phylo_VOL_hellinger_norm.rds")
comm_lMDS(comm_landmark, beta, output, threads)
rm(list = c('beta', 'output'))

# Func
beta <- readRDS(paste0(root_path, "/Diversity/01_beta/beta_func_VOL_hellinger_norm.rds"))
output <- paste0(root_path, "/Diversity/02_LMDS/lmds_func_VOL_hellinger_norm.rds")
comm_lMDS(comm_landmark, beta, output, threads)
rm(list = c('beta', 'output'))
