#'------------------------------------------------------------------------------
#' @title Integration of spectra with ordination through machine learning models
#'------------------------------------------------------------------------------

#' @description It tune, model, and predict the ordinations based on spectra.

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(caret)
library(pls)
library(plsVarSel)
library(CAST)
library(sf)
library(parallel)
library(pbmcapply)
library(bit64)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("R/auxiliary/data_split.R")
source("R/auxiliary/axis_segmentation.R")
source("R/auxiliary/iscv.R") 
source("R/auxiliary/model_tune.R") 
source("figures/figure_optimal_inspection.R") 
source("R/auxiliary/model_build.R")
source("R/auxiliary/model_performance.R")
source("R/auxiliary/pls_vip.R")
source("R/auxiliary/pls_coefficients.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"
sensor <- "DESIS"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

type <- "taxo"
diversity <- "total"
meta <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))
ordination <- fread(paste0(root_path, "/Diversity/NMDS/NMDS_scores_comm_", type, "_", diversity,".csv"))
ncol_ordination <- ncol(ordination) - 1
match_meta <- fread(paste0(root_path, "/Pixel_extraction/03_match-spectra/", sensor, "/", sensor,"_match_usefull.csv"))
spectra <- fread(paste0(root_path, "/Pixel_extraction/03_match-spectra/", sensor, "/", sensor,"_spectra_SW.csv"))
iterations <- 100
out_path <- paste0(root_path, "/Models/01_model_beta/", sensor)

#-------------------------------------------------------------------------------
#' @Data_reshape  
#-------------------------------------------------------------------------------

# Get location values from meta files
match_meta <- merge(match_meta,
                    meta[, c("PLT_CN", "STATECD", "COUNTYCD", "ELEV")],
                    by = "PLT_CN", all.x = TRUE, all.y = FALSE)

# Years threshold
if(sensor == "DESIS") {
  match_meta <- match_meta[abs(dtime)/365.25 <= 6, ]
} else if(sensor == "EMIT") {
  match_meta <- match_meta[abs(dtime)/365.25 <= 6, ]
}

# Merge meta with ordination
ordination <- merge(match_meta, 
                    ordination,
                    by = "PLT_CN", all.x = TRUE, all.y = FALSE)

# Select bands to work with
if(sensor == "DESIS") {
  
  # Bands from 449.53 to 899.53.
  spectra <- spectra[, .SD, .SDcols = c(2, 50:200)]
  
} else if(sensor == "EMIT") {
  
  # Bands from 447.74 to 1265.17, from 1488.07 to 1710.97, from 2015.60 to 2424.25
  spectra <- spectra[, .SD, .SDcols = c(2, 
                                        41:139, # from 448.50 to 1266
                                        182:209, #from 1506 to 1708
                                        264:305)] #from 2015.60 to 2424.25
}

# Merge ordination with spectra 
data <- merge(ordination,
              spectra,
              by = "PLT_CN", all.x = TRUE, all.y = FALSE)

# New frames
data <- na.exclude(data)
data <- data[order(PLT_CN)]

# Save order for other simulations
if(type == "taxo") {
  
  PLT_CN <- data$PLT_CN
  saveRDS(PLT_CN, paste0(out_path, "/", sensor, "_PLT_CN-order.rds"))
  
} else {
  
  ref <- data$PLT_CN
  PLT_CN <- readRDS(paste0(out_path, "/", sensor, "_PLT_CN-order.rds"))
  all.equal(ref, PLT_CN)
  
}

# # Outliers
# outliers <- fread(paste0(out_path, "/", sensor, "_outliers.csv"))
# data <- data[outliers$remove != TRUE, ]

# Create frames
if(sensor == "DESIS") {
  
  saveRDS(data, paste0(out_path, "/DESIS_", type,"_data.rds"))
  meta <- data[, .SD, .SDcols = 1:28]
  ordination <- data[, .SD, .SDcols = 29:(28+ncol_ordination)]
  ordination <- ordination[, 1:3]
  spectra <- data[, .SD, .SDcols = (29+ncol_ordination):ncol(data)]
  weigths <- (data$gof - min(data$gof)) / (quantile(data$gof, 0.95) - min(data$gof))
  
} else if(sensor == "EMIT") {
  
  saveRDS(data, paste0(out_path, "/EMIT_data.rds"))
  meta <- data[, .SD, .SDcols = 1:27]
  ordination <- data[, .SD, .SDcols = 28:(27+ncol_ordination)]
  ordination <- ordination[, 1:3]
  spectra <- data[, .SD, .SDcols = (28+ncol_ordination):ncol(data)]
  weigths <- (data$gof - min(data$gof)) / (quantile(data$gof, 0.95) - min(data$gof))
  
}

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Applies just with one time using taxo data as reference
if(type == "taxo") {
  
  if(sensor == "DESIS") {
    p <- 0.64
    
  } else if(sensor == "EMIT") {
    p  <- 0.638
  }
  
  # Data-split based on State-County ID
  split <- data_split(meta, p = p)
  
  # Percentage of samples for training
  length(split)/nrow(meta) # % For training
  length(split) # Training
  nrow(meta) - length(split) # Testing
  
  # Export split and cross-validation
  saveRDS(split, paste0(out_path, "/", sensor, "_split.rds"))

}

split <- readRDS(paste0(out_path, "/", sensor, "_split.rds"))

#-------------------------------------------------------------------------------
#' @Segments and the iterative spatial cross-validation strategy
#-------------------------------------------------------------------------------

# Segments
segments <- axis_segmentation(ordination_split = ordination[split,], 
                              weigths_split = weigths[split],
                              fraction = 0.85,
                              threads = 20)

saveRDS(segments, paste0(out_path, "/", sensor, "_", type, "_segments.rds"))
segments <- readRDS(paste0(out_path, "/", sensor, "_", type, "_segments.rds"))

#-------------------------------------------------------------------------------
#' @Spatial cross-validation strategy
#-------------------------------------------------------------------------------

# Get strategy
iscv <- iscv_strategy(meta_split = meta[split, ], 
                      segments = segments,
                      threads = 30)

# Export strategy
saveRDS(iscv, paste0(out_path, "/", sensor, "_", type, "_iscv.rds"))
iscv <- readRDS(paste0(out_path, "/", sensor, "_", type, "_iscv.rds"))

#-------------------------------------------------------------------------------
#' @Model tune estimation
#-------------------------------------------------------------------------------

# PLSR
models <- model_tune(meta_training = meta[split,], 
                     ordination_training = ordination[split,], 
                     spectra_training = spectra[split,],
                     segments = segments,
                     iscv = iscv,
                     ncomp_max = 30,
                     ntimes = iterations,
                     threads = 16)

# Export for inspection
fwrite(models, paste0(out_path, "/", sensor, "_", type, "_optimal.csv"))
axis <- "Axis3"
optimal_figure(models, rAxes = axis, rMetric = "RMSEP", rEstimate = "train")
optimal_figure(models, rAxes = axis, rMetric = "RMSEP", rEstimate = "CV")
optimal_figure(models, rAxes = axis, rMetric = "R2", rEstimate = "train")
optimal_figure(models, rAxes = axis, rMetric = "R2", rEstimate = "CV")
optimal_figure(models, rAxes = axis, rMetric = "PRESS", rEstimate = "PRESS")

#' @note Do a manual checking of the optimal number of components

### DESIS 
# Taxonomic c(18, 12, 11)
# Phylogenetic c(XX, XX, XX)
# Functional c(XX, XX, XX)

### EMIT
# Taxonomic c(15, 11, 11)
# Phylogenetic c(XX, XX, XX)
# Functional c(XX, XX, XX)

# Manual checked
tune_optimal <- c(12, 10, 10)

#-------------------------------------------------------------------------------
#' @Final_model for prediction
#-------------------------------------------------------------------------------

final_models <- model_build(meta_training = meta[split,], 
                            ordination_training = ordination[split,], 
                            spectra_training = spectra[split,],
                            segments = segments,
                            iscv = iscv,
                            tune_optimal = tune_optimal,
                            ntimes = iterations,
                            threads = 5)

#saveRDS(final_models, paste0(out_path, "/", sensor, "_final_models.rds")) # Do it if you really need it

#-------------------------------------------------------------------------------
#' @Model_performance
#-------------------------------------------------------------------------------

training <- model_performance(meta_split = meta[split,],
                              ordination_split = ordination[split,],
                              spectra_split = spectra[split,],
                              models = final_models,
                              tune_optimal = tune_optimal,
                              threads = 25)

testing <- model_performance(meta_split = meta[!split,],
                             ordination_split = ordination[!split,],
                             spectra_split = spectra[!split,],
                             models = final_models,
                             tune_optimal = tune_optimal,
                             threads = 25)

# Export performance values
performance_training <- training$performance
performance_training$type <- "Training"
performance_testing <- testing$performance
performance_testing$type <- "Testing"
performance <- rbind(performance_training, performance_testing)
fwrite(performance, paste0(out_path, "/", sensor, "_", type, "_performance.csv"))

# Export predicted values
predicted_training <- training$predicted
predicted_training$type <- "Training"
fwrite(predicted_training, paste0(out_path, "/", sensor, "_", type, "_predicted_training.csv"))
predicted_testing <- testing$predicted
predicted_testing$type <- "Testing"
fwrite(predicted_testing, paste0(out_path, "/", sensor, "_", type, "_predicted_testing.csv"))

# Export three-dimensional residuals
residuals_3D_training <- training$residuals_3D
residuals_3D_training$type <- "Training"
fwrite(residuals_3D_training, paste0(out_path, "/", sensor, "_", type, "_residuals_3D_training.csv"))
residuals_3D_testing <- testing$residuals_3D
residuals_3D_testing$type <- "Testing"
fwrite(residuals_3D_testing, paste0(out_path, "/", sensor, "_", type, "_residuals_3D_testing.csv"))

#-------------------------------------------------------------------------------
#' @Model_coefficients
#-------------------------------------------------------------------------------

coefficients <- pls_coefficients(models = final_models,
                                 tune_optimal = tune_optimal)
fwrite(coefficients, paste0(out_path, "/", sensor, "_", type, "_coefficients.csv"))

#-------------------------------------------------------------------------------
#' @Model_VIP
#-------------------------------------------------------------------------------

vip <- pls_vip(models = final_models,
               tune_optimal = tune_optimal)
fwrite(vip, paste0(out_path, "/", sensor, "_", type, "_vip.csv"))

#-------------------------------------------------------------------------------
#' @END
#-------------------------------------------------------------------------------
