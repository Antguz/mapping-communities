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
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "E:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme Pro/Projects/DESIS_mapping"
meta <- fread(paste0(root_path, "/00_FIA/plots_WI_MN_20221102.csv"))
ordination <- readRDS(paste0(root_path, "/05_ordination/dimensions_plots_beta_DRYBIO.rds"))
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/spectra_SW.csv"))
PD <- fread(paste0(root_path, "/03_diversity/master_phylogenetic-diversity_DRYBIO.txt"))
out_path <- paste0(root_path, "/06_model/SW_beta_PLS")

#'------------------------------------------------------------------------------
#' Load source code

source("R/data_split.R")
source("R/model_tune.R") 
source("R/mds_segmentation.R") 
source("R/model_build.R")
source("R/performance.R")
source("R/pls_vip.R")
source("R/pls_coefficients.R")

#'------------------------------------------------------------------------------
#' @Analysis  

#Meta file
meta <- meta[, c("PLT_CN", "INVYR", "STATECD", "COUNTYCD", "LAT", "LON", "ELEV")]

#Merge phylogenetic diversity
PD <- PD[, c("PLT_CN", "ntaxa", "pd.obs", "PSR")]
PD <- merge(meta, 
            PD,
            by = "PLT_CN", all.x = FALSE, all.y = TRUE)

#Merge meta with ordination
ordination <- merge(PD, 
                    ordination,
                    by = "PLT_CN", all.x = FALSE, all.y = TRUE)

#Merge ordination with spectra 
data <- merge(ordination, 
              spectra[, .SD, .SDcols = c(1, 7:ncol(spectra))],
              by = "PLT_CN", all.x = FALSE, all.y = TRUE)

#New frames
data <- na.exclude(data)
meta <- data[, .SD, .SDcols = 1:10]
ordination <- data[, .SD, .SDcols = 11:13]
spectra <- data[, .SD, .SDcols = 15:ncol(data)]
weigths <- data$goodness

#------------------------------------------
# Data-split
set.seed(1107)
split <- data_split(ordination)
saveRDS(split, paste0(out_path, "/split.rds"))

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
                     ntimes = 2,
                     threads = 28)

optimal <- models[, .(Rsq = mean(Rsquared), Rsq_SD = sd(Rsquared), 
                      RMSE = mean(RMSE), RMSE_SD = sd(RMSE)), by = c("MDS", "ncomp")]
optimal_tune <- optimal[, which.min(RMSE), by = "MDS"]$V1
ncomp_tuned <- list(expand.grid(ncomp = optimal_tune[1]),
                    expand.grid(ncomp = optimal_tune[2]),
                    expand.grid(ncomp = optimal_tune[3]))

fwrite(models, paste0(out_path, "/optimal_models.rds"))

#------------------------------------------
# Final Model of prediction

final_models <- modul_build(meta = meta,
                            ordination = ordination, 
                            spectra = spectra, 
                            split = split, 
                            weigths = weigths,
                            method = "pls",
                            tune_optimal = ncomp_tuned,
                            ntimes = 100,
                            threads = 28)

saveRDS(final_models, paste0(out_path, "/final_models.rds"))
final_models <- readRDS(paste0(out_path, "/final_models.rds"))

#------------------------------------------
# Performance

training <- model_performance(meta = meta[split,],
                              ordination = ordination[split,], 
                              spectra = spectra[split,],
                              models = final_models)

testing <- model_performance(meta = meta[!split,],
                             ordination = ordination[!split,], 
                             spectra = spectra[!split,],
                             models = final_models)

performance_training <- training$performance
performance_training$type <- "Training"
performance_testing <- testing$performance
performance_testing$type <- "Testing"
performance <- rbind(performance_training, performance_testing)
fwrite(performance, paste0(out_path, "/performance.csv"))

predicted_training <- training$predicted
predicted_training$type <- "Training"
predicted_testing <- testing$predicted
predicted_testing$type <- "Testing"
fwrite(predicted_training, paste0(out_path, "/predicted_training.csv"))
fwrite(predicted_testing, paste0(out_path, "/predicted_testing.csv"))

#------------------------------------------
# Get coefficients

coefficients <- pls_coefficients(models = final_models)
fwrite(coefficients, paste0(out_path, "/coefficients.csv"))

#------------------------------------------
# Get VIP

vip <- pls_vip(models = final_models)
fwrite(vip, paste0(out_path, "/vip.csv"))

