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
library(misc3d)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

source("R/data_split.R")
source("R/axes_segmentation.R")
source("R/iscv.R") 
source("R/model_tune.R") 
source("R/model_build.R")
source("R/model_performance.R")
source("R/pls_vip.R")
source("R/pls_coefficients.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
sensor <- "DESIS"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

type <- "func"
meta <- fread(paste0(root_path, "/00_FIA/plots_eastern_20230427.csv"))
ordination <- fread(paste0(root_path, "/05_NMDS/NMDS_scores_comm_", type, ".csv"))
ncol_ordination <- ncol(ordination) - 1
match_meta <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor, "/", sensor,"_match_usefull.csv"))
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor, "/", sensor,"_spectra_SW.csv"))
iterations <- 100
#outliers <-  fread(paste0(root_path, "/06_model/SW_beta_PLS/DESIS_outliers.csv"))
out_path <- paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor)

#-------------------------------------------------------------------------------
#' @Data_reshape  
#-------------------------------------------------------------------------------

# Get location values from meta files
match_meta <- merge(match_meta,
                    meta[, c("PLT_CN", "STATECD", "COUNTYCD", "ELEV")],
                    by = "PLT_CN", all.x = TRUE, all.y = FALSE)

# 400 days difference
if(sensor == "DESIS") {
  match_meta <- match_meta[abs(dtime) <= 455, ]
}

# Merge meta with ordination
ordination <- merge(match_meta, 
                    ordination,
                    by = "PLT_CN", all.x = TRUE, all.y = FALSE)

# Select bands to work with
if(sensor == "DESIS") {
  
  # Bands from 450 to 962.53
  spectra <- spectra[, .SD, .SDcols = c(2, 40:212)]
  
} else if(sensor == "EMIT") {
  
  # Bands from 447.74 to 1265.17, from 1488.07 to 1710.97, from 2015.60 to 2424.25
  spectra <- spectra[, .SD, .SDcols = c(2, 
                                        32:142, # from 447.74 to 1265.17
                                        172:202, #from 1488.07 to 1710.97
                                        243:298)] #from 2015.60 to 2424.25
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
  saveRDS(PLT_CN, paste0(out_path, "/", "PLT_CN-order.rds"))
  
} else {
  
  ref <- data$PLT_CN
  PLT_CN <- readRDS(paste0(out_path, "/", "PLT_CN-order.rds"))
  all.equal(ref, PLT_CN)
  
}

# Outliers
outliers <- fread(paste0(out_path, "/", sensor, "_outliers.csv"))
data <- data[outliers$remove != TRUE, ]

# Create frames
if(sensor == "DESIS") {
  
  saveRDS(data, paste0(out_path, "/DESIS_", type,"_data.rds"))
  meta <- data[, .SD, .SDcols = 1:18]
  ordination <- data[, .SD, .SDcols = 19:(18+ncol_ordination)]
  ordination <- ordination[, 1:3]
  spectra <- data[, .SD, .SDcols = (19+ncol_ordination):ncol(data)]
  #weigths <- data$goodness
  
} else if(sensor == "EMIT") {
  
  saveRDS(data, paste0(out_path, "/EMIT_data.rds"))
  meta <- data[, .SD, .SDcols = 1:17]
  ordination <- data[, .SD, .SDcols = 18:(17+ncol_ordination)]
  ordination <- ordination[, 1:3]
  spectra <- data[, .SD, .SDcols = (18+ncol_ordination):ncol(data)]
  #weigths <- data$goodness
  
}

#-------------------------------------------------------------------------------
#' @Data-split
#-------------------------------------------------------------------------------

# Applies just with one time using taxo data as reference
if(type == "taxo") {
  
  if(sensor == "DESIS") {
    p <- 0.665
    
  } else if(sensor == "EMIT") {
    p  <- 0.70
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
if(sensor == "DESIS") {
  
  segments <- pbmclapply(X = 1:iterations,
                         FUN = axes_segmentation,
                         ordination = ordination[split, ],
                         nsegments = 50,
                         lower.p = 0.01,
                         proportion = 0.8,
                         mc.set.seed = TRUE,
                         mc.cores = 20)
  
  saveRDS(segments, paste0(out_path, "/", sensor, "_", type, "_segments.rds"))
  segments <- readRDS(paste0(out_path, "/", sensor, "_", type, "_segments.rds"))
  
} else if(sensor == "EMIT") {
  
  segments <- pbmclapply(X = 1:iterations,
                         FUN = axes_segmentation,
                         ordination = ordination[split, ],
                         nsegments = 50,
                         lower.p = 0.01,
                         proportion = 0.8,
                         mc.set.seed = TRUE,
                         mc.cores = 20)
  
  saveRDS(segments, paste0(out_path, "/", sensor, "_", type, "_segments.rds"))
  segments <- readRDS(paste0(out_path, "/", sensor, "_", type, "_segments.rds"))

}

#-------------------------------------------------------------------------------
#' @Spatial cross-validation strategy
#-------------------------------------------------------------------------------

# Create a vector file
vec <- st_as_sf(x = meta,                         
                coords = c("LON", "LAT"),
                crs = "EPSG:4269")

# Project vector
vec_projected <- st_transform(vec, 
                              crs = "+proj=lcc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs")

# Get strategy
iscv <- iscv_strategy(meta = meta, 
                      vec_projected = vec_projected, 
                      split = split, 
                      segments = segments,
                      threads = 25)

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
                     threads = 20)

# Export for inspection
fwrite(models, paste0(out_path, "/", sensor, "_", type, "_optimal.csv"))
#plot(colMeans(models[type == "RMSEP" & estimate == "train" & Axes == "Axis2", 6:ncol(models)], na.rm = TRUE))
plot(colMeans(models[type == "RMSEP" & estimate == "CV" & Axes == "Axis1", 6:ncol(models)], na.rm = TRUE))
#plot(colMeans(models[type == "R2" & estimate == "train" & Axes == "Axis2", 6:ncol(models)], na.rm = TRUE))
#plot(colMeans(models[type == "PRESS" & estimate == "CV" & Axes == "Axis2", 6:ncol(models)], na.rm = TRUE))

#' @note Do a manual checking of the optimal number of components

### DESIS 
# Taxonomic c(23, 20, 17)
# Phylogenetic-log c(26, 12, 24)
# Functional c(25, 25, 25)

### EMIT
# Taxonomic c(10, 15, 12)
# Phylogenetic-log c(15, 16, 10)
# Functional c(12, 14, 9)

# Manual checked
tune_optimal <- c(25, 25, 25)

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
                            threads = 20)

#saveRDS(final_models, paste0(out_path, "/", sensor, "_final_models.rds")) # Do it if you really need it

#-------------------------------------------------------------------------------
#' @Model_performance
#-------------------------------------------------------------------------------

training <- model_performance(meta = meta[split,],
                              ordination = ordination[split,],
                              spectra = spectra[split,],
                              models = final_models,
                              tune_optimal = tune_optimal,
                              threads = 25)

testing <- model_performance(meta = meta[!split,],
                             ordination = ordination[!split,],
                             spectra = spectra[!split,],
                             models = final_models,
                             tune_optimal = tune_optimal,
                             threads = 25)

performance_training <- training$performance
performance_training$type <- "Training"
performance_testing <- testing$performance
performance_testing$type <- "Testing"
performance <- rbind(performance_training, performance_testing)
fwrite(performance, paste0(out_path, "/", sensor, "_", type, "_performance.csv"))

predicted_training <- training$predicted
predicted_training$type <- "Training"
fwrite(predicted_training, paste0(out_path, "/", sensor, "_", type, "_predicted_training.csv"))
predicted_testing <- testing$predicted
predicted_testing$type <- "Testing"
fwrite(predicted_testing, paste0(out_path, "/", sensor, "_", type, "_predicted_testing.csv"))

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
