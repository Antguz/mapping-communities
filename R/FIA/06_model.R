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
#' @param ordination
#' @param out_path

#root_path <- "E:/DESIS"
root_path <- "/media/antonio/Extreme SSD/DESIS"
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/spectra.csv"))
ordination <- fread(paste0(root_path, "/05_NMDS/PCoA_plots.csv"))
out_path <- paste0(root_path, "/06_model")

#'------------------------------------------------------------------------------
#' Load source code

source("R/FIA/optimal_number.R") 
source("R/FIA/pls-ndms_model.R")
source("R/FIA/CV_spatialBlock.R")
source("R/FIA/performance.R")

#'------------------------------------------------------------------------------
#' @Analysis  

#Match files

spectra_anchor <- spectra[, 1]
spectra_anchor$spectra_row <- 1:nrow(spectra_anchor)
plot_anchor <- ordination[, 1]
plot_anchor$plot_row <- 1:nrow(plot_anchor)
keep <- merge(spectra_anchor, plot_anchor, by = "PLT_CN", all.x = TRUE, all.y = FALSE)

#observations

spectra <- spectra[keep$spectra_row,]
meta <- spectra[, c(1:5)]
spectra <- spectra[, .SD, .SDcols = c(13:ncol(spectra))]
ordination <- ordination[keep$plot_row,]
ordination <- ordination[, .SD, .SDcols = c(6:ncol(ordination))]

#Create spectra and meta



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
                              ordination = ordination,
                              spectra = spectra,
                              int_ncomp = 50, 
                              repeated = 10, 
                              k = 10, 
                              threads = 1)

#------------------------------------------
# Final Model of prediction

model <- pls_ordination_model(meta = meta,
                        ordination = ordination, 
                        spectra = spectra, 
                        split = split, 
                        folds = 10,
                        repeated = 4,
                        opt_components = c(1, 1, 1),
                        threads = 4)
