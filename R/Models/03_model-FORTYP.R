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
library(ROSE)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("R/auxiliary/data_split_fortyp.R")
source("R/auxiliary/model_fortyp.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"
type <- "taxo"
out_path <- "/media/antonio/Work/MMDB/Models/03_models_fortyp"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Meta information
fia_meta <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))

# NMDS ordination values
ordination <- fread(paste0(root_path, "/Diversity/02_NMDS/NMDS_scores_comm_", type, "_total.csv"))

# Test for match
fia_meta_match <- match.integer64(ordination$PLT_CN, fia_meta$PLT_CN)
fia_meta <- fia_meta[fia_meta_match,]
all.equal(ordination$PLT_CN, fia_meta$PLT_CN)

# FIA forest type
forest_type <- fread(paste0(root_path, "/FIA/00_FIA/ENTIRE_COND.csv"))
forest_type <- forest_type[, c("PLT_CN", "FORTYPCD")]
FORTYPCD <- fread(paste0(root_path, "/FIA/00_FIA/ENTIRE_REF_FOREST_TYPE.csv"))
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
forest_type <- forest_type[forest_type$FORTYPCD %in% to_remove, FORTYP := NA]
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
  data_split <- data_split_fortyp(forest_type,
                                  split = 0.60)
  
  # Export split and cross-validation
  saveRDS(data_split, paste0(out_path, "/split_fortyp.rds"))
  
}

data_split <- readRDS(paste0(out_path, "/split_fortyp.rds"))

#-------------------------------------------------------------------------------
#' @Model
#-------------------------------------------------------------------------------

# Create model
model <- model_fortyp(forest_type = forest_type,
                      ordination = ordination[, 1:4],
                      data_split = data_split,
                      iterations = 100,
                      threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_model-FORTYP.csv"))

#-------------------------------------------------------------------------------
#' @END
#-------------------------------------------------------------------------------