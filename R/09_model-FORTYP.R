#'------------------------------------------------------------------------------
#' @title Model occurrence probabilities of forest types based on NMDS
#'------------------------------------------------------------------------------

#' @description It models the occurrence probability of forest types based on 
#' NMDS ordination.

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(bit64)
library(parallel)
library(pbmcapply)
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

#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "E:/Projects/DESIS_mapping"
type <- "phylo-log"
out_path <- "/media/antonio/Work/Projects/DESIS_mapping/07_model_occurrence"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Meta information
fia_meta <- fread(paste0(root_path, "/00_FIA/plots_eastern_20230427.csv"))

# NMDS ordination values
ordination <- fread(paste0(root_path, "/05_NMDS/NMDS_scores_comm_", type, ".csv"))

# Test for match
fia_meta_match <- match.integer64(ordination$PLT_CN, fia_meta$PLT_CN)
fia_meta <- fia_meta[fia_meta_match,]
all.equal(ordination$PLT_CN, fia_meta$PLT_CN)

# FIA forest type
forest_type <- fread(paste0(root_path, "/00_FIA/ENTIRE_COND.csv"))
forest_type <- forest_type[, c("PLT_CN", "FORTYPCD")]
FORTYPCD <- fread(paste0(root_path, "/00_FIA/ENTIRE_REF_FOREST_TYPE.csv"))
FORTYPCD <- FORTYPCD[, c("VALUE", "MEANING")]
forest_type <- merge(forest_type, FORTYPCD, 
                     by.x = "FORTYPCD", by.y = "VALUE",
                     all.x = TRUE, all.y = TRUE)
colnames(forest_type)[3] <- "FORTYP"

# Match plot order
FORTYP_match <- match.integer64(fia_meta$PLT_CN, forest_type$PLT_CN)
forest_type <- forest_type[FORTYP_match, ]
all.equal(fia_meta$PLT_CN, forest_type$PLT_CN)

# Clean forest type
to_remove <- c(181, 183, 223, 380, 382, 383, 803, 807, 925, 926, 
               932, 950, 951, 952, 953, 954, 955, 981, 999)
forest_type <- forest_type[!forest_type$FORTYPCD %in% to_remove, ]
forest_type <- forest_type[!is.na(FORTYP),]
forest_type <- forest_type[, c("PLT_CN", "FORTYP")]

# Match meta, NMDS, and forest type
fia_meta_match <- match.integer64(forest_type$PLT_CN, fia_meta$PLT_CN)
fia_meta <- fia_meta[fia_meta_match, ]
ordination_match <- match.integer64(forest_type$PLT_CN, ordination$PLT_CN)
ordination <- ordination[ordination_match, ]

all.equal(fia_meta$PLT_CN, forest_type$PLT_CN)
all.equal(fia_meta$PLT_CN, ordination$PLT_CN)
all.equal(ordination$PLT_CN, forest_type$PLT_CN)

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

