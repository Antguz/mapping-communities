#'------------------------------------------------------------------------------
#' @title Model occurrence probabilities of plant lineages based on NMDS
#'------------------------------------------------------------------------------

#' @description It models the probability of occurrence based on NMDS ordination 
#' and relative abundance.

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
library(MLmetrics)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

sourceCpp("src/approximation_rcpp.cpp")
source("R/aggregate_plants.R")
source("R/data_split_lineages.R")
source("R/axes_segmentation.R") 
source("R/model_lineages.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/DESIS_mapping"
type <- "phylo-log"
out_path <- "/media/antonio/Work/DESIS_mapping/07_model_occurrence"

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
aggregation <- aggregate_batch(sp_meta[,c(9, 3, 4, 10, 11, 12, 13, 14)], inventory)
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
rpresence <- rpresence > 10.0
rpresence <- c(TRUE, rpresence)
aggregation <- aggregation[ , .SD, .SDcols = c(1:ncol(aggregation))[rpresence]]

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Applies just with one time using taxo data as reference
if(type == "phylo-log") {
  
  set.seed(1107)
  
  # Data-split based on State-County ID
  data_split <- data_split_lineages(aggregation, 
                                    split = 0.6)
  
  # Export split and cross-validation
  saveRDS(data_split, paste0(out_path, "/split_lineages.rds"))
  
}

data_split <- readRDS(paste0(out_path, "/split_lineages.rds"))

#-------------------------------------------------------------------------------
#' @Model
#-------------------------------------------------------------------------------

# Create model
model <- model_lineages(fia_meta = fia_meta,
                        ordination = ordination,
                        aggregation = aggregation,
                        data_split = data_split,
                        iterations = 100,
                        threads = 32)

# Export model
fwrite(model, paste0(out_path, "/", type, "_model-occurrence.csv"))

