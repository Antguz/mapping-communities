#'------------------------------------------------------------------------------
#' @title Performance of CWM models based predicted axes
#'------------------------------------------------------------------------------

#' @description It evaluate the performance of models that predict CWM based on 
#' the predicted axes using spaceborne spectroscopy.

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(bit64)
library(parallel)
library(Rcpp)
library(RcppArmadillo)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

sourceCpp("src/cwm_rcpp.cpp")
source("R/auxiliary/model_performance.R")
source("R/auxiliary/performance_batch_CWM.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"
out_path <- "/media/antonio/Work/MMDB/Models/04_CWM_performance"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

sensor <- "DESIS"
type <- "func"

# Get glm coefficients
coefficients <- fread(paste0(root_path, "/Models/03_model_CWM/", 
                             type, "_model-CWM.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 11:22)]

# Get environmental variables
environment <- fread(paste0(root_path, "/FIA/03_FIA-env/fia-env.csv"))

# Training
dataset <- "training"

# Get predicted axes
predicted <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

# Predict CWM based on coefficients
func_perf <- performance_batch(predicted, coefficients, environment)
saveRDS(func_perf, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_CWM.rds"))

# Testing
dataset <- "testing"

# Get predicted axes
predicted <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

# Predict CWM based on coefficients
func_perf <- performance_batch(predicted, coefficients, environment)
saveRDS(func_perf, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_CWM.rds"))

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

fwrite(data, paste0(out_path, "/performance_summary_CWM.csv"))
