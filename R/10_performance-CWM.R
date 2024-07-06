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

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

sourceCpp("src/cwm_rcpp.cpp")
source("R/model_performance.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/DESIS_mapping"
out_path <- "/media/antonio/Work/DESIS_mapping/08_CWM_performance"

#-------------------------------------------------------------------------------
#' @Probabilities-function
#-------------------------------------------------------------------------------

#' @param predicted Axes predictions from spaceborne spectroscopy
#' @param coefficients Coefficients (and intercept) value from the occurrence models

get_CWM <- function(predicted, coefficients) {
  
  # Select columns
  obs <- predicted[, c("Axes", "sample", "observed")]
  pred <- predicted[, c("Axes", "sample", "predicted")]
  
  # Dcast columns
  obs <- dcast(obs, sample ~ Axes)
  pred <- dcast(pred, sample ~ Axes)
  
  # Get probabilities
  trait <- as.vector(as.matrix(coefficients[, 1]))
  obs_cwm <- multi_cwm(as.matrix(obs[,-1]), 
                       as.matrix(coefficients[, -1]))
  pred_cwm <- multi_cwm(as.matrix(pred[,-1]),
                        as.matrix(coefficients[, -1]))
  
  # Make it data.table
  obs_cwm <- as.data.table(obs_cwm)
  colnames(obs_cwm) <- trait
  pred_cwm <- as.data.table(pred_cwm)
  colnames(pred_cwm) <- trait
  
  # Return list
  list(obs = obs_cwm,
       pred = pred_cwm)
  
}

#-------------------------------------------------------------------------------
#' @Performance-function
#-------------------------------------------------------------------------------

#' @param all_probabilities A list with the observed and predicted probabilities

performance_batch <- function(all_cwm) {
  
  application_performance <- function(X,
                                      all_cwm) {
    
    # Predict new values
    per <- parameters(obs = as.matrix(all_cwm$obs)[, X], 
                      pred = as.matrix(all_cwm$pred)[,X])
    
    return(per)
  }
  
  # Performance in parallel
  performance <- pbmclapply(X = 1:ncol(all_cwm$obs),
                            FUN = application_performance,
                            all_cwm = all_cwm,
                            mc.cores = 28)
  
  performance <- do.call(rbind.data.frame, performance)
  performance <- cbind(plants = colnames(all_cwm$obs),
                       performance)
  
  return(performance)
  
}

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

sensor <- "EMIT"
dataset <- "training"

# Functional
type <- "func"
coefficients <- fread(paste0(root_path, "/07_model_CWM/", 
                             type, "_model-CWM.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 11:14)]

predicted <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

func_CWM <- get_CWM(predicted, coefficients)
func_perf <- performance_batch(func_CWM)
func_CWM$perf <- func_perf
saveRDS(func_CWM, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_CWM.rds"))


#-------------------------------------------------------------------------------
#' @Collect-performance
#-------------------------------------------------------------------------------

folder_path <- paste0(root_path, "/08_CWM_performance")

files <- list.files(path = folder_path, 
                    pattern = ".rds", 
                    recursive = TRUE,
                    include.dirs = TRUE)

data <- data.table()

for(i in 1:length(files)) {
  
  sensor <- strsplit(files[i], "_")[[1]][1]
  type <- strsplit(files[i], "_")[[1]][2]
  dataset <- strsplit(files[i], "_")[[1]][3]
  
  performance <- readRDS(paste0(folder_path, "/", files[i]))$perf
  performance$sensor <- sensor
  performance$type <- type
  performance$dataset <- dataset
  
  data <- rbind(data, performance)
  
}

fwrite(data, paste0(root_path, "/08_CWM_performance/performance_summary_CWM.csv"))
