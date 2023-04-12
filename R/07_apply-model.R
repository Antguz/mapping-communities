#'------------------------------------------------------------------------------
#' @title Apply PLS model to DESIS scenes
#'------------------------------------------------------------------------------

#' @description It applies the model to DESIS scenes.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(terra)
library(Rcpp)
library(RcppArmadillo)
library(parallel)

#'------------------------------------------------------------------------------
#' Source code

source("R/DESIS-transform.R")
sourceCpp("src/model-application.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param scene_main_path
#' @param coefficients
#' @param output_raster
#' @param threads

#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "E:/Projects/DESIS_mapping"
scene_main_path <- paste0(root_path, "/00_imagery/", "DESIS-HSI-L2A-DT0599346912 015-20210615T161539-V0215")
coefficients <- fread(paste0(root_path, "/06_model/SW_beta_PLS/", "coefficients.csv"))
output_raster <- paste0(scene_main_path, "_predicted.tif")

#'------------------------------------------------------------------------------
#' Function

apply_model <- function(scene_main_path, 
                        coefficients, 
                        output_raster,
                        threads = 28) {
  
  #Colnames
  bands <- colnames(coefficients)
  bands <- bands[4:length(bands)]
  
  #Transform scene
  sw <- desis_transform(scene_main_path, 
                        DESIS_bands,
                        threads = threads)
  
  sw <- sw[, .SD, .SDcols = 25:195]
  colnames(sw) <- bands
  
  sw <- fread("E:/Projects/DESIS_mapping/sw_test.csv")
  
  #--------------------------------------------------------
  # Makes coefficients matrices
  
  spectra <- as.matrix(sw)
  spectra <- as.matrix(spectra)
  MDS1 <- as.matrix(coefficients[MDS == "1", .SD, .SDcols = 3:ncol(coefficients)])
  MDS2 <- as.matrix(coefficients[MDS == "2", .SD, .SDcols = 3:ncol(coefficients)])
  MDS3 <- as.matrix(coefficients[MDS == "3", .SD, .SDcols = 3:ncol(coefficients)])
  
  #--------------------------------------------------------
  #Function for batch application
  
  predicted_ordination <- predict_pls(spectra,
                                      MDS1,
                                      MDS2,
                                      MDS3)
  
  predictions <- as.data.table(predictions)
  colnames(predictions) <- c("Healthy", "Symtomatic", "Dead", 
                             "amp-healthy", "amp-symtomatic", "amp-dead")
  
  
}