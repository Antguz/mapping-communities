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
source("R/auxiliary/model_performance.R")

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"
out_path <- "/media/antonio/Work/MMDB/Models/04_occurence_performance"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

sensor <- "EMIT"

#### Phylogenetic training ####
dataset <- "training"
type <- "phylo"

# Get glm coefficients
coefficients <- fread(paste0(root_path, "/Models/03_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 14:17)]

# Get predicted axes
predicted <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

# Predict probabilities based on coefficients
phylo_perf <- performance_batch(predicted, coefficients)
saveRDS(phylo_perf, 
        paste0(out_path, "/", sensor, "_", type, "_", dataset, "_probabilities.rds"))

#### Phylogenetic testing ####
dataset <- "testing"
type <- "phylo"
coefficients <- fread(paste0(root_path, "/Models/03_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- subset(coefficients, dataset == "testing")
coefficients <- coefficients[, .SD, .SDcols = c(1, 14:17)]

# Predict probabilities based on coefficients
predicted <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                          sensor, "_", type, "_predicted_", dataset, ".csv"))

# Predict probabilities based on coefficients
phylo_perf <- performance_batch(predicted, coefficients)
saveRDS(phylo_perf, 
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

fwrite(data, paste0(out_path, "/performance_summary_lineages.csv"))

