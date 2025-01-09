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

source("R/auxiliary/lmds.R")

#'------------------------------------------------------------------------------
#' @Function

comm_lMDS <- function(comm_landmark, beta, output, threads) {
  
  # Get PLT_CN
  PLT_CN <- as.integer64(rownames(beta$Btotal))
  #PLT_CN <- as.integer64(rownames(x))
  
  # Match landmarks
  match_landmarks <- match(comm_landmark, PLT_CN)
  
  # Get landmarks
  beta_landmarks <- select_landmarks(x = beta$Btotal/10000,
                                     landmarks_vector = match_landmarks,
                                     threads = 30)
  # Get MDS
  beta_lmds <- cmdscale_landmarks(beta_landmarks, ndim = 3, rescale = TRUE)
  
  # Export MDS
  saveRDS(beta_lmds, output)
  
}

#'------------------------------------------------------------------------------
#' @Application

# Define root path
root_path <- "/media/antonio/Work/MMDB"
comm_landmark <- readRDS(paste0(root_path, "/Diversity/01_beta/comm_landmarks.rds"))
threads <- 99

beta <- readRDS(paste0(root_path, "/Diversity/01_beta/beta_taxo_VOL.rds"))
output <- paste0(root_path, "/Diversity/02_LMDS/lmds_taxo_VOL.rds")
comm_lMDS(comm_landmark, beta, output, threads)