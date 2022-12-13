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
root_path <- "F:/DESIS"
root_path <- "/media/antonio/Extreme SSD/DESIS"
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/spectra.csv"))
ordination <- fread(paste0(root_path, "/05_NMDS/PCoA_plots_PD.csv"))
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
ordination <- ordination[, .SD, .SDcols = c(2:ncol(ordination))]

#Create spectra and meta
smooth <- as.data.frame(spectra)

for(i in 1:nrow(smooth)) {
  
  sample <- as.numeric(smooth[i, ])
  filter <- sgolayfilt(sample, p = 1, n = 5, m = 0, ts = 1)
  smooth[i, ] <- filter
  
}

smooth <- as.data.table(smooth)
spectra <- smooth[, .SD, .SDcols = c(3:(ncol(smooth)-2))]
smooth <- as.data.frame(spectra)
plot(as.numeric(smooth[1,]))

for(i in 1:nrow(smooth)) {
  
  sample <- as.numeric(smooth[i, ])
  norm <- sqrt(sum(sample^2))
  smooth[i, ] <- sample / norm
  
}

for(i in 1:nrow(smooth)) {
  
  sample <- as.numeric(smooth[i, ])
  sample_mean <- mean(sample)
  sample_sd <- sd(sample)
  smooth[i, ] <- (sample - sample_mean) / sample_sd
  
}

spectra <- as.data.table(smooth)
plot(as.numeric(spectra[1,]))

#------------------------------------------
# Data-split
set.seed(1107)
split <- createDataPartition(meta$PLT_CN, 
                             p = 0.6,
                             list = FALSE, 
                             times = 1)

#------------------------------------------
# Estimation of the optimal number of components

ncomponents <- optimal_number(meta = meta,
                              ordination = ordination,
                              spectra = spectra,
                              split = split,
                              int_ncomp = 20, 
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
