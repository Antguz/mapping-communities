#'------------------------------------------------------------------------------
#' @title Performance of model probabilities based predicted axes
#'------------------------------------------------------------------------------

#' @description It evaluate the performance of models based on the probability 
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
source("R/model_performance.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/DESIS_mapping"
out_path <- "/media/antonio/Work/DESIS_mapping/08_occurence_performance"

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

performance_batch <- function(all_probabilities, threads) {
  
  all_probabilities$obs[all_probabilities$obs == "NaN"] <- NA
  all_probabilities$pred[all_probabilities$pred == "NaN"] <- NA
  
  application_performance <- function(X,
                                      all_probabilities) {
    
    cat(X)
    
    name <- colnames(all_probabilities$obs)[X]
    
    if(max(as.matrix(all_probabilities$obs)[, X]) == 0 &
       max(as.matrix(all_probabilities$pred)[,X]) == 0) {
      
    } else {
      
      # Predict new values
      per <- try(parameters(obs = as.matrix(all_probabilities$obs)[, X], 
                            pred = as.matrix(all_probabilities$pred)[,X]),
                 silent = TRUE)
      
      if(class(per)[1] == "try-error") {
        per <- data.table(R2 = NA, BIAS = NA, RMSE = NA, 
                          perRMSE = NA, intercept = NA, slope = NA)
      }
    }
    
    per <- cbind(plant = name, per)
    return(per)
    
  }
  
  # Performance in parallel
  performance <- pbmclapply(X = 1:ncol(all_probabilities$obs),
                            FUN = application_performance,
                            all_probabilities = all_probabilities,
                            mc.cores = threads)
  
  performance <- do.call(rbind.data.frame, performance)
  
  return(performance)
  
}

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

sensor <- "DESIS"

#### Phylogenetic training ####
dataset <- "training"
type <- "phylo-log"
coefficients <- fread(paste0(root_path, "/07_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 15:18)]
coefficients <- coefficients[plant != "angiosperm" & plant != "broadleaf", ]

predicted <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

phylo_probabilities <- get_probabilities(predicted, coefficients)
phylo_perf <- performance_batch(all_probabilities = phylo_probabilities,
                                threads = 30)
phylo_probabilities$perf <- phylo_perf
saveRDS(phylo_probabilities, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_probabilities.rds"))

#### Phylogenetic testing ####
dataset <- "testing"
type <- "phylo-log"
coefficients <- fread(paste0(root_path, "/07_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 15:18)]
coefficients <- coefficients[plant != "angiosperm" & plant != "broadleaf", ]

predicted <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

phylo_probabilities <- get_probabilities(predicted, coefficients)
phylo_perf <- performance_batch(all_probabilities = phylo_probabilities,
                                threads = 30)
phylo_probabilities$perf <- phylo_perf
saveRDS(phylo_probabilities, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_probabilities.rds"))

#### Taxonomic training ####
dataset <- "training"
type <- "taxo"
coefficients <- fread(paste0(root_path, "/07_model_fortyp/", 
                             type, "_model-FORTYP.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 15:18)]

predicted <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

taxo_probabilities <- get_probabilities(predicted, coefficients)
taxo_perf <- performance_batch(all_probabilities = taxo_probabilities,
                               threads = 28)
taxo_probabilities$perf <- taxo_perf
saveRDS(taxo_probabilities, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_probabilities.rds"))

#### Taxonomic testing ####
dataset <- "testing"
type <- "taxo"
coefficients <- fread(paste0(root_path, "/07_model_fortyp/", 
                             type, "_model-FORTYP.csv"))
coefficients <- subset(coefficients, dataset == "testing")
coefficients <- coefficients[, .SD, .SDcols = c(1, 15:18)]

predicted <- fread(paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

taxo_probabilities <- get_probabilities(predicted, coefficients)
taxo_perf <- performance_batch(all_probabilities = taxo_probabilities,
                               threads = 30)
taxo_probabilities$perf <- taxo_perf
saveRDS(taxo_probabilities, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_probabilities.rds"))