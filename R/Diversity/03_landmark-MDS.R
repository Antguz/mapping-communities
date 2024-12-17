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
#' @Perform_LMDS

# Define root path
root_path <- "/media/antonio/Work/MMDB"

# comm_landmark <- readRDS(paste0(root_path, "/Diversity/01_beta/comm_landmarks.rds"))
# # beta <- readRDS(paste0(root_path, "/Diversity/01_beta/beta_taxo_VOL.rds"))
# x <- readRDS(paste0(root_path, "/Diversity/01_beta/beta_taxo-total_VOL.rds"))

comm_landmark <- readRDS("data/comm_landmarks.rds")
x <- readRDS("data/beta_taxo-total_VOL.rds")

comm_lMDS <- function(comm_landmark, beta) {
  
  # Get PLT_CN
  PLT_CN <- as.integer64(rownames(beta$Btotal))
  PLT_CN <- as.integer64(rownames(x))
  
  # Match landmarks
  match_landmarks <- match(comm_landmark, PLT_CN)
  
  # Get MDS
  beta_landmarks <- select_landmarks(x = beta$Btotal/10000,
                                     landmarks_vector = match_landmarks,
                                     threads = 25)
  
  beta_lmds <- cmdscale_landmarks(beta_landmarks, ndim = 3, rescale = TRUE)
  
}