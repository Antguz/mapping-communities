#'------------------------------------------------------------------------------
#' @title Compute species approximation of the predicted values
#'------------------------------------------------------------------------------

#' @description It computes the approximation of the predicted values to species.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(terra)
library(Rcpp)
library(RcppArmadillo)

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/approximation_sp.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param communities_projections
#' @param species_projections

#'------------------------------------------------------------------------------
#' Function

sp_affinity <- function(communities_projections, 
                        species_projections) {
  
  #---------------------------------------------------------
  # Estimate distance
  distance <- approximation_sp(as.matrix(communities_projections[, 2:4]), 
                               as.matrix(species_projections[, 2:4]))
  
  #Get the inverse
  distance <- 1/distance
  total <- rowSums(distance, na.rm = FALSE)
  
  #Get fraction
  fraction <- distance/total
  #fraction <- 1 - (distance/total)
  fraction <- round(fraction*1000000, 0)
  fraction <- as.data.table(fraction)
  colnames(fraction) <- species_projections$Taxon2
  fraction <- cbind(communities_projections[, 1], fraction)
  
  return(fraction)
  
}

#'------------------------------------------------------------------------------
#' @example

# Run in a single scene

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"

sensor <- "DESIS"
type <- "taxonomic"
data <- "testing"

# Read and reshape community projections
communities_projections <- fread(paste0(root_path, "/06_model/SW_beta_PLS/", sensor, "/", 
                                        sensor, "_", type, "_predicted_", data, ".csv"))
communities_projections <- communities_projections[, c(1, 2, 4)]
communities_projections <- dcast(communities_projections, sample ~ Axes)

# Select just the dataset of interest
meta <- readRDS(paste0(root_path, "/06_model/SW_beta_PLS/", sensor, "/", sensor, "_", type, "_data.rds"))
split <- readRDS(paste0(root_path, "/06_model/SW_beta_PLS/", sensor, "/", sensor, "_split.rds"))
meta <- meta[!split, ]
communities_projections <- cbind(meta[,1], communities_projections[, 2:4])

# Read species projection
species_projections <- fread(paste0(root_path, "/05_MDS/scores_sp_", type, ".csv"))
species_projections <- species_projections[, 1:4]

# Estimate affinity
affinity <- sp_affinity(communities_projections, species_projections)
fwrite(affinity, paste0(root_path, "/07_affinity/", 
                        sensor, "_", type, "_", data, "_affinity.csv"))

