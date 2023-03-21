#'------------------------------------------------------------------------------
#' @title Integration of PLSR and NMDS
#'------------------------------------------------------------------------------

#' @description It model the of species PLSR with NMDS.
#' @return .txt file with the stack of results

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
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param spectra
#' @param ordination
#' @param out_path

#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "/media/antonio/Extreme Pro/Projects/DESIS_mapping"
meta <- fread(paste0(root_path, "/00_FIA/plots_WI_MN_20221102.csv"))
ordination <- readRDS(paste0(root_path, "/05_ordination/dimensions_plots_taxonomic_DRYBIO.rds"))
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/spectra_VN.csv"))
PD <- fread(paste0(root_path, "/03_diversity/master_phylogenetic-diversity.txt"))
out_path <- paste0(root_path, "/06_model")

#'------------------------------------------------------------------------------
#' Load source code

source("R/FIA/optimal_sigma.R") 
source("R/FIA/gausspr_build.R")
source("R/FIA/performance.R")

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
meta <- data[, .SD, .SDcols = 1:10]
ordination <- data[, .SD, .SDcols = 11:13]
spectra <- data[, .SD, .SDcols = 14:ncol(data)]

#------------------------------------------
# Data-split
set.seed(1107)
split <- createDataPartition(meta$PLT_CN, 
                             p = 0.70,
                             list = FALSE, 
                             times = 1)

#------------------------------------------
# Sigma estimation model

model <- optimal_sigma(meta = meta,
                       ordination = ordination, 
                       spectra = spectra, 
                       split = split, 
                       threads = 28)

saveRDS(model, paste0(root_path, "/06_model/model_sigma.rds"))

#------------------------------------------
# Final Model of prediction

model <- gausspr_build(meta = meta,
                       ordination = ordination, 
                       spectra = spectra, 
                       split = split, 
                       sigma_components = c(0.01107595, 0.01359494, 0.01863291),
                       repeated = 100,
                       threads = 28)

saveRDS(model, paste0(root_path, "/06_model/model_final.rds"))

#------------------------------------------
# Performance

training <- model_performance(meta = meta[split[,1],],
                              ordination = ordination[split[,1],], 
                              spectra = spectra[split[,1],],
                              models = model)

testing <- model_performance(meta = meta[!split[,1],],
                             ordination = ordination[!split[,1],], 
                             spectra = spectra[!split[,1],],
                             models = model)

performance_training <- training$performance
performance_training$type <- "Training"
performance_testing <- testing$performance
performance_testing$type <- "Testing"
performance <- rbind(performance_training, performance_testing)
fwrite(performance, paste0(root_path, "/06_model/performance.csv"))

predicted_training <- training$predicted
predicted_training$type <- "Training"
predicted_testing <- testing$predicted
predicted_testing$type <- "Testing"
fwrite(predicted_training, paste0(root_path, "/06_model/predicted_training.csv"))
fwrite(predicted_testing, paste0(root_path, "/06_model/predicted_testing.csv"))
