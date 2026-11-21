#'------------------------------------------------------------------------------
#' @title Model occurrence probabilities of plant lineages based on MDS - spatial CV
#'------------------------------------------------------------------------------

#' @description It models the occurrence probability of plant lineages (species,
#' genera, and families present in more than 10% of the plots) using binomial
#' GLMs based on three sets of predictors: MDS axes, climate and topography, and
#' altogether. Class imbalance is handled with ROSE within training folds, and
#' models are evaluated with repeated k-fold spatial block cross-validation
#' (200 km blocks).

#' @return .rds files with the presence/absence of lineages per plot, residuals,
#' spatial autocorrelation, and blocks, and .csv files with the scaling values
#' and the GLM coefficients and performance per set of predictors.

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(bit64)
library(parallel)
library(pbmcapply)
library(caret)
library(ROSE)
library(Rcpp)
library(vegan)
library(sf)
library(blockCV)
library(MLmetrics)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

sourceCpp("src/approximation_rcpp.cpp")
source("R/auxiliary/aggregate_plants.R")
source("R/auxiliary/model_lineages_spatial.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/guzman/Work/MMDB"
type <- "phylo"
out_path <- "/media/guzman/Work/MMDB/Models/03_model_occurrence"

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

# FIA data for aggregation
inventory <- readRDS(paste0(root_path, "/Diversity/01_beta/master_matched_comm_", type, "_VOL.rds"))
sp_meta <- fread(paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv"))
aggregation <- aggregate_batch(sp_meta[,c("Taxon2",
                                          "Family",
                                          "Genus")],
                               decostand(inventory, "total"))
aggregation_match <- match(ordination$PLT_CN, aggregation$PLT_CN)
aggregation <- aggregation[aggregation_match, ]
all.equal(aggregation$PLT_CN, ordination$PLT_CN)

# Test for match
fia_meta_match <- match(ordination$PLT_CN, fia_meta$PLT_CN)
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

# Data for testing occurrence
saveRDS(aggregation, paste0(out_path, "/LINEAGES_occurrance.rds"))
saveRDS(ordination, paste0(out_path, "/LINEAGES_variables.rds"))

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# No data split
data_split <- NULL

#-------------------------------------------------------------------------------
#' @Spatial-blocks
#-------------------------------------------------------------------------------

# Residuals of the full models (axes, climate, and topography) per plant lineage.
lineages_residuals <- residuals_lineages(aggregation = aggregation,
                                         ordination = ordination,
                                         variables = "altogether")
saveRDS(lineages_residuals, paste0(out_path, "/residuals_lineages.rds"))

# Plot coordinates projected to CONUS Albers (meters)
plots <- st_as_sf(cbind(fia_meta[, .(LON, LAT)], lineages_residuals[, -1]),
                  coords = c("LON", "LAT"),
                  crs = 4326)
plots <- st_transform(plots, 5070)

# Block size based on the spatial autocorrelation range of the residuals
autocorrelation <- autocor_lineages(x = plots,
                                    columns = colnames(lineages_residuals)[-1],
                                    num_sample = 5000,
                                    threads = 30)
autocorrelation$range_table
block_size <- autocorrelation$range
block_size <- 200000
saveRDS(autocorrelation, paste0(out_path, "/autocorrelation_lineages.rds"))

# Square blocks of block_size
coordinates <- st_coordinates(plots)
blocks <- paste0(floor(coordinates[, 1]/block_size), "_",
                 floor(coordinates[, 2]/block_size))
length(unique(blocks)) # Number of blocks
saveRDS(blocks, paste0(out_path, "/blocks_lineages.rds"))

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
fwrite(scaling_values, paste0(out_path, "/", type, "_scaling-occurrence.csv"))

# Reconstruct independent variables
ordination <- cbind(data.table(PLT_CN = PLT_CN),
                    as.data.table(ordination))

#-------------------------------------------------------------------------------
#' @Model
#-------------------------------------------------------------------------------

# ----------------------------------------------------------
# Model - just ordination axes as predictors

variables <- "axes"
model <- model_lineages_spatial(ordination = ordination,
                                aggregation = aggregation,
                                data_split = data_split,
                                blocks = blocks,
                                variables = variables,
                                k = 5,
                                n_tries = 100,
                                iterations = 100,
                                threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_", variables, "_model-occurrence.csv"))

# ----------------------------------------------------------
# Model - just climate and topography information as predictors

variables <- "climate-topography"
model <- model_lineages_spatial(ordination = ordination,
                                aggregation = aggregation,
                                data_split = data_split,
                                blocks = blocks,
                                variables = variables,
                                k = 5,
                                n_tries = 100,
                                iterations = 100,
                                threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_", variables, "_model-occurrence.csv"))

# ----------------------------------------------------------
# Model - axes, climate, and topography information as predictors

variables <- "altogether"
model <- model_lineages_spatial(ordination = ordination,
                                aggregation = aggregation,
                                data_split = data_split,
                                blocks = blocks,
                                variables = variables,
                                k = 5,
                                n_tries = 100,
                                iterations = 100,
                                threads = 30)

# Export model
fwrite(model, paste0(out_path, "/", type, "_", variables, "_model-occurrence.csv"))

#-------------------------------------------------------------------------------
#' @END
#-------------------------------------------------------------------------------
