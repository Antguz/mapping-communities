#'------------------------------------------------------------------------------
#' @title Model the community weighted mean of plant traits
#'------------------------------------------------------------------------------

#' @description It models the CWM of plant traits based on NMDS ordination 
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

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("R/CWM.R")
source("R/auxiliary/CWM.R")
source("R/auxiliary/data_split_CWM.R")
source("R/auxiliary/model_CWM.R")
source("R/auxiliary/model_performance.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"
type <- "func"
out_path <- "/media/antonio/Work/MMDB/Models/03_model_CWM"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Meta information
fia_meta <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))

# NMDS ordination values
ordination <- fread(paste0(root_path, "/Diversity/02_NMDS/NMDS_scores_comm_", type, "_total.csv"))

# FIA data for aggregation
inventory <- readRDS(paste0(root_path, "/Diversity/01_beta/master_matched_comm_", type, "_VOL.rds"))
sp_meta <- fread(paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv"))
match <- match(colnames(inventory), sp_meta$Taxon2)
sp_meta <- sp_meta[match, ]
aggregation <- CWM_batch(sp_meta[,c("Taxon2", "HTq", "SLq", "WD", "shade", "drought", "waterlogging", "LMA", "N")], 
                         inventory)

# Test for match
all.equal(aggregation$PLT_CN, ordination$PLT_CN)
fia_meta_match <- match.integer64(ordination$PLT_CN, fia_meta$PLT_CN)
fia_meta <- fia_meta[fia_meta_match,]
all.equal(aggregation$PLT_CN, fia_meta$PLT_CN)

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Applies just with one time using taxo data as reference
if(type == "func") {
  
  set.seed(1107)
  
  # Data-split based on State-County ID
  data_split <- data_split_CWM(aggregation, 
                               split = .60)
  
  # Export split and cross-validation
  saveRDS(data_split, paste0(out_path, "/split_CWM.rds"))
  
}

data_split <- readRDS(paste0(out_path, "/split_CWM.rds"))

#-------------------------------------------------------------------------------
#' @Model
#-------------------------------------------------------------------------------

# Create model
model <- model_CWM(fia_meta = fia_meta,
                   ordination = ordination,
                   aggregation = aggregation,
                   data_split = data_split,
                   iterations = 100,
                   threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_model-CWM.csv"))

