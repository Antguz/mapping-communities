#'------------------------------------------------------------------------------
#' @title Performance of FORTYP's model probabilities based predicted axes
#'------------------------------------------------------------------------------

#' @description It evaluate the performance of FORTYP's models based on the probability 
#' predicted from NDMS axes using spaceborne spectroscopy.

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(bit64)
library(parallel)
library(pbmcapply)
library(Rcpp)
library(RcppArmadillo)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

sourceCpp("src/probability_rcpp.cpp")
source("R/auxiliary/model_performance.R")
source("R/auxiliary/performance_batch_FORTYR.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"
out_path <- "/media/antonio/Work/MMDB/Models/04_fortyp_performance"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

sensor <- "DESIS"

#### Taxonomic training #### ---------------------------------------------------
dataset <- "training"
type <- "taxo"

# Get glm coefficients
coefficients <- fread(paste0(root_path, "/Models/03_models_fortyp/", type, "_model-FORTYP.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 14:17)]

# Predict probabilities based on coefficients
predicted <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

# Predict and export
taxo_perf <- performance_batch(predicted, coefficients)
saveRDS(taxo_perf, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_probabilities.rds"))

#### Taxonomic testing ####
dataset <- "testing"
type <- "taxo"

# Get glm coefficients
coefficients <- fread(paste0(root_path, "/Models/03_models_fortyp/", type, "_model-FORTYP.csv"))
coefficients <- subset(coefficients, dataset == "testing")
coefficients <- coefficients[, .SD, .SDcols = c(1, 14:17)]

# Predict probabilities based on coefficients
predicted <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

taxo_perf <- performance_batch(predicted, coefficients)
saveRDS(taxo_perf, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_probabilities.rds"))

#-------------------------------------------------------------------------------
#' @Collect-performance
#-------------------------------------------------------------------------------

folder_path <- out_path

files <- list.files(path = folder_path,
                    pattern = ".rds",
                    recursive = TRUE,
                    include.dirs = TRUE)

data <- data.table()

for(i in 1:length(files)) {
  
  sensor <- strsplit(files[i], "_")[[1]][1]
  type <- strsplit(files[i], "_")[[1]][2]
  dataset <- strsplit(files[i], "_")[[1]][3]
  
  performance <- readRDS(paste0(folder_path, "/", files[i]))
  performance$sensor <- sensor
  performance$type <- type
  performance$dataset <- dataset
  
  data <- rbind(data, performance)
  
}

fwrite(data, paste0(out_path, "/performance_summary_fortyp.csv"))
