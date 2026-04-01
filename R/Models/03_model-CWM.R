#'------------------------------------------------------------------------------
#' @title Model the community weighted mean of plant traits
#'------------------------------------------------------------------------------

#' @description It models the CWM of plant traits based on MDS ordination 
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
library(scales)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("R/auxiliary/CWM.R")
source("R/auxiliary/data_split_CWM.R")
source("R/auxiliary/model_CWM.R")
source("R/auxiliary/model_performance.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/guzman/Work/MMDB"
type <- "func"
out_path <- "/media/guzman/Work/MMDB/Models/03_model_CWM"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Meta information
fia_meta <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))

# MDS ordination values
ordination <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_scores_comm_", type, "_total.csv"))

# FIA data for aggregation
inventory <- readRDS(paste0(root_path, "/Diversity/01_beta/master_matched_comm_", type, "_VOL.rds"))
sp_meta <- fread(paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv"))
match <- match(colnames(inventory), sp_meta$Taxon2)
sp_meta <- sp_meta[match, ]
aggregation <- CWM_batch(sp_meta[,c("Taxon2", "HTq", "SLq", "WD", "shade", "drought", "waterlogging", "LMA", "N")], 
                         decostand(inventory, "total"))

# Get environmental data
fia_environment <- fread(paste0(root_path, "/FIA/03_FIA-env/fia-env.csv"))
ordination <- merge(ordination, fia_environment, by = "PLT_CN", all.x = FALSE, all.y = FALSE)
ordination <- na.exclude(ordination)

# Test for match
all.equal(aggregation$PLT_CN, ordination$PLT_CN)
match <- match.integer64(ordination$PLT_CN, aggregation$PLT_CN)
aggregation <- aggregation[match, ]
all.equal(aggregation$PLT_CN, ordination$PLT_CN)
fia_meta_match <- match.integer64(ordination$PLT_CN, fia_meta$PLT_CN)
fia_meta <- fia_meta[fia_meta_match,]
all.equal(aggregation$PLT_CN, fia_meta$PLT_CN)

# Data for ploting
saveRDS(aggregation, paste0(out_path, "/CWM_plant-traits.rds"))
saveRDS(ordination, paste0(out_path, "/CWM_variables.rds"))

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Applies just with one time using taxo data as reference
if(type == "func") {

  # Data-split based on State-County ID
  data_split <- data_split_CWM(aggregation, 
                               split = .60)
  
  # Export split and cross-validation
  saveRDS(data_split, paste0(out_path, "/split_CWM.rds"))
  
}

data_split <- readRDS(paste0(out_path, "/split_CWM.rds"))

#-------------------------------------------------------------------------------
#' @Scale-ordination
#-------------------------------------------------------------------------------

# Keep a copy of PLT_CN
PLT_CN <- ordination$PLT_CN

# Scale independent variables
ordination <- scale(ordination[,-1])

# Get scaling values (new_data - mean_values) / sd_values
mean_values <- attr(ordination, "scaled:center")  # Get the means
sd_values <- attr(ordination, "scaled:scale")  # Get the standard deviations

# Merge scaling values and export
scaling_values <- rbind(mean_values, sd_values)
fwrite(scaling_values, paste0(out_path, "/", type, "_scaling-CWM.csv"))

# Reconstruct independent variables
ordination <- cbind(data.table(PLT_CN = PLT_CN), 
                    as.data.table(ordination))

#-------------------------------------------------------------------------------
#' @Model
#-------------------------------------------------------------------------------

# ----------------------------------------------------------
# Model - just ordination axes as predictors

variables <- "axes"
model <- model_CWM(fia_meta = fia_meta,
                   ordination = ordination,
                   aggregation = aggregation,
                   variables = variables,
                   data_split = data_split,
                   iterations = 100,
                   threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_", variables,"_model-CWM.csv"))

# ----------------------------------------------------------
# Model - just climate and topography information as predictors

variables <- "climate-topography"
model <- model_CWM(fia_meta = fia_meta,
                   ordination = ordination,
                   aggregation = aggregation,
                   variables = variables,
                   data_split = data_split,
                   iterations = 100,
                   threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_", variables,"_model-CWM.csv"))

# ----------------------------------------------------------
# Model - axes, climate, and topography information as predictors

variables <- "altogether"
model <- model_CWM(fia_meta = fia_meta,
                   ordination = ordination,
                   aggregation = aggregation,
                   variables = variables,
                   data_split = data_split,
                   iterations = 100,
                   threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_", variables,"_model-CWM.csv"))

#-------------------------------------------------------------------------------
#' @END
#-------------------------------------------------------------------------------