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
# library(ROSE)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

sourceCpp("src/approximation_rcpp.cpp")
source("R/auxiliary/aggregate_plants.R")
source("R/auxiliary/data_split_lineages.R")
source("R/auxiliary/model_lineages.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"
type <- "phylo"
out_path <- "/media/antonio/Work/MMDB/Models/03_model_occurrence"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Meta information
fia_meta <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))

# NMDS ordination values
ordination <- fread(paste0(root_path, "/Diversity/02_LMDS/MDS_scores_comm_", type, "_total.csv"))

# Get enviromental data
fia_environment <- fread(paste0(root_path, "/FIA/03_FIA-env/fia-env.csv"))
ordination <- merge(ordination, fia_environment, by = "PLT_CN", all.x = FALSE, all.y = FALSE)
ordination <- na.exclude(ordination)

# FIA data for aggregation
inventory <- readRDS(paste0(root_path, "/Diversity/01_beta/master_matched_comm_", type, "_VOL.rds"))
sp_meta <- fread(paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv"))
aggregation <- aggregate_batch(sp_meta[,c("Taxon2", "Family", "Genus", "growth_form", 
                                          "mycorrhizal_association", "leaf_shape", 
                                          "leaf_type", "leaf_habit")], 
                               decostand(inventory, "total"))
match <- match.integer64(ordination$PLT_CN, aggregation$PLT_CN)
aggregation <- aggregation[match, ]
all.equal(aggregation$PLT_CN, ordination$PLT_CN)

# Test for match
fia_meta_match <- match.integer64(ordination$PLT_CN, fia_meta$PLT_CN)
fia_meta <- fia_meta[fia_meta_match,]
all.equal(aggregation$PLT_CN, fia_meta$PLT_CN)

# Select the most abundant groups
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
if(type == "phylo") {
  
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
model <- model_lineages(ordination = ordination,
                        aggregation = aggregation,
                        data_split = data_split,
                        iterations = 100,
                        threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_model-occurrence.csv"))

