#'------------------------------------------------------------------------------
#' @title Model probabilities of plant groups based on NMDS and relative abundance
#'------------------------------------------------------------------------------

#' @description It models the probability of occurrence based on NMDS ordination 
#' and relative abundance.

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(bit64)
library(parallel)
library(caret)
library(Rcpp)
library(vegan)
library(CAST)
library(sf)
library(plotROC)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

sourceCpp("src/approximation_rcpp.cpp")
source("R/aggregate_plants.R")
source("R/data_split.R")
source("R/axes_segmentation.R") 
source("R/model_occurence.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"
type <- "func"
out_path <- "/media/antonio/Work/Projects/DESIS_mapping/07_model_occurrence"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Meta information
fia_meta <- fread(paste0(root_path, "/00_FIA/plots_eastern_20230427.csv"))

# NMDS ordination values
ordination <- fread(paste0(root_path, "/05_NMDS/NMDS_scores_comm_", type, ".csv"))

# FIA data for aggregation
inventory <- readRDS(paste0(root_path, "/03_diversity/master_matched_comm_", type, "-VOL.rds"))
sp_meta <- fread(paste0(root_path, "/00_FIA/SPCD_meta.csv"))
aggregation <- aggregate_batch(sp_meta[,c(9, 10)], inventory)
#aggregation <- aggregation[, 1:2]

# Test for match
all.equal(aggregation$PLT_CN, ordination$PLT_CN)
fia_meta_match <- match.integer64(ordination$PLT_CN, fia_meta$PLT_CN)
fia_meta <- fia_meta[fia_meta_match,]
all.equal(aggregation$PLT_CN, fia_meta$PLT_CN)

# Select the most abundant genus
aggregation <- aggregation[, names(aggregation)[-1] := lapply(.SD, function(x) as.integer(x!=0)), 
                                 .SDcols = 2:ncol(aggregation)]
presence <- colSums(aggregation[,-1])
rpresence <- round(presence/nrow(aggregation)*100, 4)
rpresence <- rpresence > 5.0
rpresence <- c(TRUE, rpresence)
aggregation <- aggregation[ , .SD, .SDcols = c(1:ncol(aggregation))[rpresence]]

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Applies just with one time using taxo data as reference
if(type == "taxo") {
  
  set.seed(1107)
  
  # Data-split based on State-County ID
  split <- data_split(fia_meta, p = .62)
  
  # Percentage of samples for training
  length(split)/nrow(fia_meta) # % For training
  length(split) # Training
  nrow(fia_meta) - length(split) # Testing
  
  # Export split and cross-validation
  saveRDS(split, paste0(out_path, "/data_split.rds"))
  
}

split <- readRDS(paste0(out_path, "/data_split.rds"))

#-------------------------------------------------------------------------------
#' @Cross-validation strategy
#-------------------------------------------------------------------------------

# Cross-validation strategy
if(type == "taxo") {
  
  # Create a vector file
  vec <- st_as_sf(x = fia_meta,                         
                  coords = c("LON", "LAT"),
                  crs = "EPSG:4269")
  
  # Project vector
  vec_projected <- st_transform(vec, 
                                crs = "+proj=lcc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")
  
  # Get vector for training and testing
  split_sf <- 1:nrow(fia_meta) %in% split
  train_vec <- vec_projected[split_sf,]
  test_vec <- vec_projected[!split_sf,]
  
  knndm_folds <- knndm(tpoints = vec_projected[split_sf,],
                       predpoints = vec_projected[!split_sf,],
                       space = "geographical",
                       k = 10,
                       maxp = 0.5,
                       clustering = "kmeans",
                       sampling = "regular")
  
  plot(knndm_folds)
  plot(knndm_folds, type = "simple") # For more accessible legend labels
  plot(knndm_folds, type = "simple", stat = "density") # To visualize densities rather than ECDFs
  
  saveRDS(knndm_folds, paste0(out_path, "/", "knndm_folds.rds"))
  
}

knndm_folds <- readRDS(paste0(out_path, "/", "knndm_folds.rds"))

#-------------------------------------------------------------------------------
#' @Model
#-------------------------------------------------------------------------------

# Create model
model <- model_occurence(fia_meta = fia_meta,
                         ordination = ordination,
                         aggregation = aggregation,
                         split = split,
                         knndm_folds = knndm_folds,
                         threads = 25)

# Export model
fwrite(model, paste0(out_path, "/", type, "_model-occurrence.csv"))

