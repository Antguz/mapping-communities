#'------------------------------------------------------------------------------
#' @title Model occurrence probabilities of forest types based on MDS - spatial CV
#'------------------------------------------------------------------------------

#' @description It models the occurrence probability of forest types based on 
#' MDS ordination.

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(bit64)
library(parallel)
library(pbmcapply)
library(caret)
library(ROSE)
library(sf)
library(blockCV)
library(MLmetrics)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("R/auxiliary/data_split_fortyp.R")
source("R/auxiliary/model_fortyp_spatial.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/guzman/Work/MMDB"
type <- "taxo"
out_path <- "/media/guzman/Work/MMDB/Models/03_models_fortyp"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# Meta information
fia_meta <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))

# NMDS ordination values
ordination <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_scores_comm_", type, "_total.csv"))

# Get enviromental data
fia_environment <- fread(paste0(root_path, "/FIA/03_FIA-env/fia-env.csv"))
ordination <- merge(ordination, fia_environment, by = "PLT_CN", all.x = FALSE, all.y = FALSE)
ordination <- na.exclude(ordination)

# Test for match
fia_meta_match <- match(ordination$PLT_CN, fia_meta$PLT_CN)
fia_meta <- fia_meta[fia_meta_match,]
all.equal(ordination$PLT_CN, fia_meta$PLT_CN)

# FIA forest type
COND <- fread(paste0(root_path, "/FIA/00_FIA/ENTIRE_COND.csv"))
COND <- COND[COND_STATUS_CD == 1]
setorder(COND, PLT_CN, -CONDPROP_UNADJ)
COND <- COND[, .SD[1], by = PLT_CN]
COND <- COND[, .(PLT_CN, FORTYPCD)]
COND <- na.exclude(COND)

# Read forest types
FORTYPCD <- fread(paste0(root_path, "/FIA/00_FIA/ENTIRE_REF_FOREST_TYPE.csv"))
FORTYPCD <- FORTYPCD[, c("VALUE", "MEANING")]

# Merge and clean
FORTYPCD <- merge(COND, FORTYPCD, 
                  by.x = "FORTYPCD", by.y = "VALUE",
                  all.x = FALSE, all.y = FALSE)
colnames(FORTYPCD)[3] <- "FORTYP"
to_remove <- c(181, 183, 223, 380, 382, 383, 803, 807, 925, 926, 
               932, 950, 951, 952, 953, 954, 955, 981, 999)
to_remove <- FORTYPCD$FORTYPCD %in% to_remove
FORTYPCD <- FORTYPCD[to_remove == FALSE,]
FORTYPCD <- FORTYPCD[, .(PLT_CN, FORTYP)]

# Match plot order
FORTYP_match <- match(fia_meta$PLT_CN, FORTYPCD$PLT_CN)
FORTYPCD <- FORTYPCD[FORTYP_match, ]
FORTYPCD <- na.exclude(FORTYPCD)

# Match meta, NMDS, and forest type
fia_meta_match <- match(FORTYPCD$PLT_CN, fia_meta$PLT_CN)
fia_meta <- fia_meta[fia_meta_match, ]
ordination_match <- match(FORTYPCD$PLT_CN, ordination$PLT_CN)
ordination <- ordination[ordination_match, ]

all.equal(fia_meta$PLT_CN, FORTYPCD$PLT_CN)
all.equal(fia_meta$PLT_CN, ordination$PLT_CN)
all.equal(ordination$PLT_CN, FORTYPCD$PLT_CN)

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# No data split
data_split <- NULL

#-------------------------------------------------------------------------------
#' @Spatial-blocks
#-------------------------------------------------------------------------------

# Residuals of the full models (axes, climate, and topography) per forest type.
fortyp_residuals <- residuals_fortyp(forest_type = FORTYPCD,
                                     ordination = ordination,
                                     variables = "altogether")
saveRDS(fortyp_residuals, paste0(out_path, "/residuals_fortyp.rds"))

# Plot coordinates projected to CONUS Albers (meters)
plots <- st_as_sf(cbind(fia_meta[, .(LON, LAT)], fortyp_residuals[, -1]),
                  coords = c("LON", "LAT"), 
                  crs = 4326)
plots <- st_transform(plots, 5070)

# Block size based on the spatial autocorrelation range of the residuals
autocorrelation <- autocor_fortyp(x = plots,
                                  columns = colnames(fortyp_residuals)[-1],
                                  num_sample = 5000,
                                  threads = 30)
autocorrelation$range_table
block_size <- autocorrelation$range
block_size <- 200000
saveRDS(autocorrelation, paste0(out_path, "/autocorrelation_fortyp.rds"))

# Square blocks of block_size
coordinates <- st_coordinates(plots)
blocks <- paste0(floor(coordinates[, 1]/block_size), "_", 
                 floor(coordinates[, 2]/block_size))
length(unique(blocks)) # Number of blocks
saveRDS(blocks, paste0(out_path, "/blocks_fortyp.rds"))

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
fwrite(scaling_values, paste0(out_path, "/", type, "_scaling-FORTYP.csv"))

# Reconstruct independent variables
ordination <- cbind(data.table(PLT_CN = PLT_CN), 
                    as.data.table(ordination))

#-------------------------------------------------------------------------------
#' @Model
#-------------------------------------------------------------------------------

# ----------------------------------------------------------
# Model - just ordination axes as predictors

variables <- "axes"
model <- model_fortyp_spatial(forest_type = FORTYPCD,
                              ordination = ordination,
                              data_split = data_split,
                              blocks = blocks,
                              variables = variables,
                              k = 5,
                              n_tries = 100,
                              iterations = 100,
                              threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_", variables, "_model-FORTYP.csv"))

# ----------------------------------------------------------
# Model - just climate and topography information as predictors

variables <- "climate-topography"
model <- model_fortyp_spatial(forest_type = FORTYPCD,
                              ordination = ordination,
                              data_split = data_split,
                              blocks = blocks,
                              variables = variables,
                              k = 5,
                              n_tries = 100,
                              iterations = 100,
                              threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_", variables, "_model-FORTYP.csv"))

# ----------------------------------------------------------
# Model - axes, climate, and topography information as predictors

variables <- "altogether"
model <- model_fortyp_spatial(forest_type = FORTYPCD,
                              ordination = ordination,
                              data_split = data_split,
                              blocks = blocks,
                              variables = variables,
                              k = 5,
                              n_tries = 100,
                              iterations = 100,
                              threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_", variables, "_model-FORTYP.csv"))

#-------------------------------------------------------------------------------
#' @END
#-------------------------------------------------------------------------------
