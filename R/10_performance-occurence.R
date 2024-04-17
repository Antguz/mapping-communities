#'------------------------------------------------------------------------------
#' @title Performance of model probabilities based predicted axes
#'------------------------------------------------------------------------------

#' @description It evaluate the performance of models the probability based on the
#' predicted axes using spaceborne spectroscopy.

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

sourceCpp("src/probability_rcpp.cpp")
source("R/model_performance.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
out_path <- "/media/antonio/Work/Projects/DESIS_mapping/08_occurence_performance"

#-------------------------------------------------------------------------------
#' @Probabilities-function
#-------------------------------------------------------------------------------

#' @param predicted Axes predictions from spaceborne spectroscopy
#' @param coefficients Coefficients (and intercept) value from the occurrence models

get_probabilities <- function(predicted, coefficients) {
  
  # Select columns
  obs <- predicted[, c("Axes", "sample", "observed")]
  pred <- predicted[, c("Axes", "sample", "predicted")]
  
  # Dcast columns
  obs <- dcast(obs, sample ~ Axes)
  pred <- dcast(pred, sample ~ Axes)
  
  # Get probabilities
  plant <- as.vector(as.matrix(coefficients[, 1]))
  obs_prob <- multi_probability(as.matrix(obs[,-1]), 
                                as.matrix(coefficients[, -1]))
  pred_prob <- multi_probability(as.matrix(pred[,-1]),
                                 as.matrix(coefficients[, -1]))
  
  # Make it data.table
  obs_prob <- as.data.table(obs_prob)
  colnames(obs_prob) <- plant
  pred_prob <- as.data.table(pred_prob)
  colnames(pred_prob) <- plant
  
  # Return list
  list(obs = obs_prob/10000,
       pred = pred_prob/10000)
  
}

#-------------------------------------------------------------------------------
#' @Performance-function
#-------------------------------------------------------------------------------

#' @param all_probabilities A list with the observed and predicted probabilities

performance_batch <- function(all_probabilities) {
  
  application_performance <- function(X,
                                      obs,
                                      pred) {
    # Predict new values
    per <- parameters(obs = as.matrix(obs)[, X], pred = as.matrix(pred)[,X])
    return(per)
  }
  
  # Performance in parallel
  performance <- mclapply(X = 1:ncol(all_probabilities$obs),
                          FUN = application_performance,
                          obs = all_probabilities$obs,
                          pred = all_probabilities$pred,
                          mc.cores = 4)
  
  performance <- do.call(rbind.data.frame, performance)
  performance <- cbind(plants = colnames(all_probabilities$obs),
                       performance)
  
  return(performance)
  
}

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

sensor <- "EMIT"
dataset <- "testing"

# Taxonomic
type <- "taxo"
coefficients <- fread(paste0(root_path, "/07_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 14:17)]

predicted <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

taxo_probabilities <- get_probabilities(predicted, coefficients)
taxo_perf <- performance_batch(taxo_probabilities)
taxo_probabilities$perf <- taxo_perf
saveRDS(taxo_probabilities, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_probabilities.rds"))

# Phylogenetic
type <- "phylo-log"
coefficients <- fread(paste0(root_path, "/07_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 14:17)]

predicted <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

phylo_probabilities <- get_probabilities(predicted, coefficients)
phylo_perf <- performance_batch(phylo_probabilities)
phylo_probabilities$perf <- phylo_perf
saveRDS(phylo_probabilities, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_probabilities.rds"))


# Functional
type <- "func"
coefficients <- fread(paste0(root_path, "/07_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 14:17)]

predicted <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

func_probabilities <- get_probabilities(predicted, coefficients)
func_perf <- performance_batch(func_probabilities)
func_probabilities$perf <- func_perf
saveRDS(func_probabilities, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_probabilities.rds"))

