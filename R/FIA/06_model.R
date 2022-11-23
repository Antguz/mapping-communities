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
library(blockCV)
library(sf)
library(parallel)
library(doSNOW)

options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param spectra
#' @param NDMS
#' @param out_path

root_path <- "E:/DESIS"
#root_path <- "/media/antonio/Extreme SSD/DESIS"
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/spectra.csv"))
NDMS <- fread(paste0(root_path, "/05_NMDS/NMDS_sites.csv"))
out_path <- paste0(root_path, "/06_model")

#'------------------------------------------------------------------------------
#' Load source code

source("R/FIA/optimal_number.R") 
source("R/FIA/pls-ndms_model.R")
source("R/FIA/CV_spatialBlock.R")
source("R/FIA/performance.R")

#'------------------------------------------------------------------------------
#' @Analysis  

meta <- spectra[, .SD, .SDcols=c(1:5)]
spectra <- spectra[, .SD, .SDcols=c(15:ncol(spectra))]
NDMS <- NDMS[, .SD, .SDcols=c(6:ncol(NDMS))]

#------------------------------------------
# Data-split
set.seed(1107)
split <- createDataPartition(meta$PLT_CN, 
                             p = 0.75,
                             list = FALSE, 
                             times = 1)

#------------------------------------------
# Estimation of the optimal number of components

ncomponents <- optimal_number(meta = meta,
                              NDMS = NDMS,
                              spectra = spectra,
                              int_ncomp = 50, 
                              repeated = 10, 
                              k = 10, 
                              threads = 1)

#------------------------------------------
# Final Model of prediction

model <- pls_ndms_model(meta = meta,
                        NDMS = NDMS, 
                        spectra = spectra, 
                        split = split, 
                        folds = 10,
                        repeated = 4,
                        opt_components = c(1, 1, 1),
                        threads = 4)
