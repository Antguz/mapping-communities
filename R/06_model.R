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
library(metrica)

#'------------------------------------------------------------------------------
#' Working directory

#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
sensor <- "DESIS"

#'------------------------------------------------------------------------------
#' Read information

meta <- fread(paste0(root_path, "/00_FIA/plots_eastern_20230427.csv"))
ordination <- fread(paste0(root_path, "/05_MDS/PCoA_scores_comm_func.csv"))
ncol_ordination <- ncol(ordination) - 1
match_meta <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor, "/", sensor,"_match_usefull.csv"))
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor, "/", sensor,"_spectra_SW.csv"))
type <- "func"
#outliers <-  fread(paste0(root_path, "/06_model/SW_beta_PLS/DESIS_outliers.csv"))
out_path <- paste0(root_path, "/06_model/SW_beta_PLS/", sensor)

#'------------------------------------------------------------------------------
#' Load source code

source("R/data_split.R")
source("R/model_tune.R") 
source("R/model_build.R")
source("R/model_predict.R")
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
                                        30:140, # from 447.74 to 1265.17
                                        170:200, #from 1488.07 to 1710.97
                                        241:296)] #from 2015.60 to 2424.25
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
if(type == "taxo ") {
  PLT_CN <- data$PLT_CN
  saveRDS(PLT_CN, paste0(out_path, "/", "PLT_CN-order.rds"))
} else {
  ref <- data$PLT_CN
  PLT_CN <- readRDS(paste0(out_path, "/", "PLT_CN-order.rds"))
  all.equal(ref, PLT_CN)
}

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
  meta <- data[, .SD, .SDcols = 1:15]
  ordination <- data[, .SD, .SDcols = 16:18]
  spectra <- data[, .SD, .SDcols = 20:ncol(data)]
  #weigths <- data$goodness
  
}

#-------------------------------------------------------------------------------
# Data-split and spatial cross-validation
#-------------------------------------------------------------------------------

set.seed(1107)

# Data-split based on State-County ID
split <- data_split(meta, p = .6)

# Percentage of samples for training
length(split)/nrow(meta) # % For training
length(split) # Training
nrow(meta) - length(split) # Testing

# Define iterative cross-validation strategy
iscv <- list()

for(i in 1:100) {
  
  scv <- CreateSpacetimeFolds(x = meta[split,],
                              spacevar = "STATECD",
                              timevar = NA,
                              k = length(unique(meta$STATECD)),
                              class = "COUNTYCD",
                              seed = sample(1:1000, 1))
  iteration <- scv$index
  
  iscv[[i]] <- iteration
  
}

# Export split and cross-validation
#saveRDS(split, paste0(out_path, "/", sensor, "_split.rds"))
#saveRDS(iscv, paste0(out_path, "/", sensor, "_iscv.rds"))

#------------------------------------------
# Model tune estimation

split <- readRDS(paste0(out_path, "/", sensor, "_split.rds"))
iscv <- readRDS(paste0(out_path, "/", sensor, "_iscv.rds"))

# PLSR
models <- model_tune(meta = meta,
                     ordination = ordination, 
                     spectra = spectra, 
                     split = split,
                     ncomp_max = 40,
                     ntimes = 100,
                     threads = 28)

# Export for inspection
fwrite(models, paste0(out_path, "/", sensor, "_", type, "_optimal.csv"))
plot(colMeans(models[type == "R2" & estimate == "train" & Axes == "PCoA2", 6:45]))

#' @note Do a manual checking of the optimal number of components

# DESIS 12, 14, 17
# EMIT 12, 14, 13

# Manual checked
ncomp_tuned <- c(20, 20, 30)

#------------------------------------------
# Final Model of prediction

final_models <- model_build(meta = meta,
                            ordination = ordination, 
                            spectra = spectra, 
                            split = split, 
                            iscv = iscv,
                            tune_optimal = ncomp_tuned,
                            ntimes = 10,
                            threads = 4)

#saveRDS(final_models, paste0(out_path, "/", sensor, "_final_models.rds")) # Do it if you really need it

#------------------------------------------
# Performance

training <- model_predict(meta = meta[split,],
                          ordination = ordination[split,],
                          spectra = spectra[split,],
                          models = final_models,
                          threads = 5)

testing <- model_predict(meta = meta[!split,],
                         ordination = ordination[!split,],
                         spectra = spectra[!split,],
                         models = final_models,
                         threads = 5)

performance_training <- training$performance
performance_training$type <- "Training"
performance_testing <- testing$performance
performance_testing$type <- "Testing"
performance <- rbind(performance_training, performance_testing)
fwrite(performance, paste0(out_path, "/", sensor, "_performance.csv"))

predicted_training <- training$predicted
predicted_training$type <- "Training"
predicted_testing <- testing$predicted
predicted_testing$type <- "Testing"
fwrite(predicted_training, paste0(out_path, "/", sensor, "_predicted_training.csv"))
fwrite(predicted_testing, paste0(out_path, "/", sensor, "_predicted_testing.csv"))

#------------------------------------------
# Get coefficients

coefficients <- pls_coefficients(models = final_models)
fwrite(coefficients, paste0(out_path, "/", sensor, "_coefficients.csv"))

#------------------------------------------
# Get VIP

vip <- pls_vip(models = final_models)
fwrite(vip, paste0(out_path, "/", sensor, "_vip.csv"))
