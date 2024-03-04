#'------------------------------------------------------------------------------
#' @title Integration of spectra with ordination through machine learning models
#'------------------------------------------------------------------------------

#' @description It tune, model, and predict the ordinations based on spectra.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(caret)
library(pls)
library(plsVarSel)
library(CAST)
library(parallel)
library(misc3d)

#'------------------------------------------------------------------------------
#' Working directory

#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
sensor <- "EMIT"

#'------------------------------------------------------------------------------
#' Read information

type <- "phylo-log"
meta <- fread(paste0(root_path, "/00_FIA/plots_eastern_20230427.csv"))
ordination <- fread(paste0(root_path, "/05_MDS/PCoA_scores_comm_", type, ".csv"))
ncol_ordination <- ncol(ordination) - 1
match_meta <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor, "/", sensor,"_match_usefull.csv"))
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor, "/", sensor,"_spectra_SW.csv"))
iterations <- 100
#outliers <-  fread(paste0(root_path, "/06_model/SW_beta_PLS/DESIS_outliers.csv"))
out_path <- paste0(root_path, "/06_model/SW_beta_PLS/", sensor)

#'------------------------------------------------------------------------------
#' Load source code

source("R/data_split.R")
source("R/axes_segmentation.R") 
source("R/model_tune.R") 
source("R/model_build.R")
source("R/model_performance.R")
source("R/pls_vip.R")
source("R/pls_coefficients.R")

#' @Analysis  
#-------------------------------------------------------------------------------
# Data reshaping
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
#data <- data[outlier[remove == FALSE, sample],] #Outlier

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
#outliers <- fread(paste0(out_path, "/DESIS_outliers.csv"))
#data <- data[outliers$remove == FALSE, ]

# Create frames
if(sensor == "DESIS") {
  
  saveRDS(data, paste0(out_path, "/DESIS_", type,"_data.rds"))
  #data <- readRDS(paste0(out_path, "/DESIS_", type,"_data.rds"))
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
# Data-split and spatial cross-validation
#-------------------------------------------------------------------------------

# Applies just with one time using taxo data as reference
if(type == "taxo") {
  
  set.seed(1107)
  
  # Data-split based on State-County ID
  split <- data_split(meta, p = .665)
  
  # Percentage of samples for training
  length(split)/nrow(meta) # % For training
  length(split) # Training
  nrow(meta) - length(split) # Testing
  
  # Export split and cross-validation
  saveRDS(split, paste0(out_path, "/", sensor, "_split.rds"))

}

split <- readRDS(paste0(out_path, "/", sensor, "_split.rds"))

#-------------------------------------------------------------------------------
# Creates segments and the iterative spatial cross-validation strategy
#-------------------------------------------------------------------------------

# Segments
if(sensor == "DESIS") {
  
  segments <- mclapply(X = 1:iterations,
                       FUN = PCoA_segmentation,
                       ordination = ordination[split, ],
                       nsegments = 50,
                       lower.p = 0.01,
                       proportion = 0.8,
                       mc.set.seed = FALSE,
                       mc.cores = 28)
  
  saveRDS(segments, paste0(out_path, "/", sensor, "_", type, "_segments.rds"))
  segments <- readRDS(paste0(out_path, "/", sensor, "_", type, "_segments.rds"))
  
} else if(sensor == "EMIT") {
  
  segments <- list()
  
  for(i in 1:iterations) {
    segments[[i]] <- 1:nrow(ordination[split, ])
  }
  
  saveRDS(segments, paste0(out_path, "/", sensor, "_", type, "_segments.rds"))
  segments <- readRDS(paste0(out_path, "/", sensor, "_", type, "_segments.rds"))
  
}


# Cross-validation strategy
iscv <- list()

for(i in 1:iterations) {
  
  scv <- CreateSpacetimeFolds(x = meta[split, ][segments[[i]],],
                              spacevar = "STATECD",
                              k = 10,
                              #k = length(unique(meta$STATECD)),
                              class = "COUNTYCD")
  
  iscv[[i]] <- scv$index
  
}

saveRDS(iscv, paste0(out_path, "/", sensor, "_", type, "_iscv.rds"))
iscv <- readRDS(paste0(out_path, "/", sensor, "_", type, "_iscv.rds"))

#-------------------------------------------------------------------------------
# Model tune estimation
#-------------------------------------------------------------------------------

# PLSR
models <- model_tune(meta_training = meta[split,], 
                     ordination_training = ordination[split,], 
                     spectra_training = spectra[split,],
                     segments = segments,
                     iscv = iscv,
                     ncomp_max = 40,
                     ntimes = iterations,
                     threads = 24)

# Export for inspection
fwrite(models, paste0(out_path, "/", sensor, "_", type, "_optimal.csv"))
plot(colMeans(models[type == "RMSEP" & estimate == "train" & Axes == "PCoA3", 6:ncol(models)], na.rm = TRUE))
plot(colMeans(models[type == "RMSEP" & estimate == "CV" & Axes == "PCoA3", 6:ncol(models)], na.rm = TRUE))
plot(colMeans(models[type == "R2" & estimate == "train" & Axes == "PCoA3", 6:ncol(models)], na.rm = TRUE))

#' @note Do a manual checking of the optimal number of components

### DESIS 
# Taxonomic 23
# Phylogenetic 26
# Phylogenetic-log 26
# Functional 26

### EMIT
# Taxonomic 14
# Phylogenetic 20
# Phylogenetic-log 20
# Functional 18


# Manual checked
tune_optimal <- c(20)

#-------------------------------------------------------------------------------
# Final model for prediction
#-------------------------------------------------------------------------------

final_models <- model_build(meta_training = meta[split,], 
                            ordination_training = ordination[split,], 
                            spectra_training = spectra[split,],
                            segments = segments,
                            iscv = iscv,
                            tune_optimal = tune_optimal,
                            ntimes = iterations,
                            threads = 22)

#saveRDS(final_models, paste0(out_path, "/", sensor, "_final_models.rds")) # Do it if you really need it

#-------------------------------------------------------------------------------
# Model performance
#-------------------------------------------------------------------------------

training <- model_performance(meta = meta[split,],
                              ordination = ordination[split,],
                              spectra = spectra[split,],
                              models = final_models,
                              threads = 20)

testing <- model_performance(meta = meta[!split,],
                             ordination = ordination[!split,],
                             spectra = spectra[!split,],
                             models = final_models,
                             threads = 20)

performance_training <- training$performance
performance_training$type <- "Training"
performance_testing <- testing$performance
performance_testing$type <- "Testing"
performance <- rbind(performance_training, performance_testing)
fwrite(performance, paste0(out_path, "/", sensor, "_", type, "_performance.csv"))

predicted_training <- training$predicted
predicted_training$type <- "Training"
predicted_testing <- testing$predicted
predicted_testing$type <- "Testing"
fwrite(predicted_training, paste0(out_path, "/", sensor, "_", type, "_predicted_training.csv"))
fwrite(predicted_testing, paste0(out_path, "/", sensor, "_", type, "_predicted_testing.csv"))

#-------------------------------------------------------------------------------
# Model coefficients
#-------------------------------------------------------------------------------

coefficients <- pls_coefficients(models = final_models)
fwrite(coefficients, paste0(out_path, "/", sensor, "_", type, "_coefficients.csv"))

#-------------------------------------------------------------------------------
# Model VIP
#-------------------------------------------------------------------------------

vip <- pls_vip(models = final_models)
fwrite(vip, paste0(out_path, "/", sensor, "_", type, "_vip.csv"))

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
