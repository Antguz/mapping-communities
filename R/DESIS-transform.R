#'------------------------------------------------------------------------------
#' @title Spectral transformation for DESIS
#'------------------------------------------------------------------------------

#' @description Transform DESIS imagery using summed wavelet spectra.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(CWT)

#'------------------------------------------------------------------------------
#' Arguments

#' @param scene_main_path A path to the scene folder
#' @param DESIS_bands A vector with the bands of DESIS
#' @param threads N threads to run in parallel

DESIS_bands <- readRDS("data/DESIS_bands.rds")

#'------------------------------------------------------------------------------
#' Function

DESIS_transform <- function(spectra, 
                            DESIS_bands,
                            cell,
                            new_bands = seq(min(DESIS_bands), max(DESIS_bands), by = 3.0),
                            threads = 28) {
  
  # Apply resampling
  resampled <- resampling_FWHM(spectra = spectra,
                               wavelengths = DESIS_bands,
                               new_wavelengths = new_bands,
                               FWHM = rep(3, length(new_bands)),
                               threads = threads)
  
  # Remove artifacts of resampling (no nice but works!)
  resampled$'401.53' <- resampled$'407.53'
  resampled$'404.53' <- resampled$'407.53'
  resampled$'998.53' <- resampled$'995.53'
  
  # Apply CWT transformation
  transformed <- cwt(t = resampled,
                     scales = c(3, 4, 5),
                     variance = 1,
                     summed_wavelet = TRUE,
                     threads = threads)
  
  transformed <- cbind(cell = cell, transformed) 
  transformed <- transformed[order(cell)]
  
  # Return transformation
  return(transformed)
  
}
