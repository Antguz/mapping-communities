#'------------------------------------------------------------------------------
#' @title Integration of spectra with ordination through machine learning models
#'------------------------------------------------------------------------------

#' @description It tune, model, and predict the ordination based on spectra.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(caret)
library(pls)
library(plsVarSel)
library(CAST)
library(sf)
library(doParallel)
library(metrica)
library(rTLS)
#library(trqwe)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
sensor <- "DESIS"

#'------------------------------------------------------------------------------
#' Read information

meta <- fread(paste0(root_path, "/00_FIA/plots_eastern_20230427.csv"))
ordination <- readRDS(paste0(root_path, "/05_ordination/communities_scores_VOL.rds"))
match_meta <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor,"_match_usefull.csv"))
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor,"_spectra_SW.csv"))
outliers <-  fread(paste0(root_path, "/06_model/SW_beta_PLS/DESIS_outliers.csv"))
out_path <- paste0(root_path, "/06_model/SW_beta_PLS/", sensor)

#'------------------------------------------------------------------------------
#' Load source code

source("R/data_split.R")
source("R/mds_segmentation.R") 
source("R/model_tune.R") 
source("R/model_build.R")
source("R/model_predict.R")
source("R/pls_vip.R")
source("R/pls_coefficients.R")

#'------------------------------------------------------------------------------
#' @Analysis  

# Outliers removal
match_meta <- merge(outliers, match_meta, 
                    by = c("PLT_CN"), 
                    all.x = TRUE, all.y = FALSE)

# Get location values from meta files
match_meta <- merge(match_meta,
                    meta[, c("PLT_CN", "STATECD", "COUNTYCD", "ELEV")],
                    by = "PLT_CN", all.x = TRUE, all.y = FALSE)

# Merge meta with ordination
ordination <- merge(match_meta, 
                    ordination,
                    by = "PLT_CN", all.x = TRUE, all.y = FALSE)

# Select bands to work with
if(sensor == "DESIS") {
  
  # Bands from 450 to 902
  spectra <- spectra[, .SD, .SDcols = c(2, 40:211)]
  
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
#data <- data[outlier[remove == FALSE, sample],] #Outlier

if(sensor == "DESIS") {
  
  saveRDS(data, paste0(out_path, "/DESIS_data.rds"))
  meta <- data[, .SD, .SDcols = 1:17]
  ordination <- data[, .SD, .SDcols = 18:20]
  spectra <- data[, .SD, .SDcols = 22:ncol(data)]
  weigths <- data$goodness
  
} else if(sensor == "EMIT") {
  
  saveRDS(data, paste0(out_path, "/EMIT_data.rds"))
  meta <- data[, .SD, .SDcols = 1:15]
  ordination <- data[, .SD, .SDcols = 16:18]
  spectra <- data[, .SD, .SDcols = 20:ncol(data)]
  weigths <- data$goodness
  
}

#------------------------------------------
# Data-split
set.seed(1107)

#DESIS 0.565 - EMIT 0.38
split <- data_split(ordination, ptraining = 0.565, steps = 15)

# Percentage of samples for training
length(split)/nrow(ordination)
nrow(ordination) - length(split)

# Export split
saveRDS(split, paste0(out_path, "/", sensor, "_split.rds"))
split <- readRDS(paste0(out_path, "/", sensor, "_split.rds"))

#------------------------------------------
# Model tune estimation

# PLSR
pls_tune <- expand.grid(ncomp = 1:30)
models <- model_tune(meta = meta,
                     ordination = ordination, 
                     spectra = spectra, 
                     split = split, 
                     weigths = weigths,
                     method = "pls",
                     tune_range = pls_tune,
                     ntimes = 10,
                     threads = 28)

fwrite(models, paste0(out_path, "/", sensor, "_optimal_models_test.csv"))

#' @note Do a manual checking of the optimal number of components

# DESIS 12, 14, 17
# EMIT 12, 14, 13

# Manual checked
ncomp_tuned <- list(expand.grid(ncomp = 24), #21
                    expand.grid(ncomp = 28), #20
                    expand.grid(ncomp = 21)) #13

#------------------------------------------
# Final Model of prediction

final_models <- model_build(meta = meta,
                            ordination = ordination, 
                            spectra = spectra, 
                            split = split, 
                            weigths = weigths,
                            method = "pls",
                            tune_optimal = ncomp_tuned,
                            ntimes = 10,
                            threads = 28)

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
