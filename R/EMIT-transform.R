#'------------------------------------------------------------------------------
#' @title Spectral transformation for EMIT
#'------------------------------------------------------------------------------

#' @description Transform EMIT imagery using summed wavelet.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(CWT)

#'------------------------------------------------------------------------------
#' Arguments

#' @param spectra A frame of spectral data
#' @param EMIT_bands A vector with the bands of EMIT
#' @param cell An ID for cell number
#' @param new_bands A new integer for the bands to compute
#' @param threads N threads to run in parallel

EMIT_bands <- readRDS("data/EMIT_bands.rds")

#'------------------------------------------------------------------------------
#' Function

EMIT_transform <- function(spectra, 
                           EMIT_bands,
                           cell,
                           new_bands = seq(min(EMIT_bands), max(EMIT_bands), by = 7.43),
                           threads = 28) {
  
  # Apply resampling
  resampled <- resampling_FWHM(spectra = spectra,
                               wavelengths = EMIT_bands,
                               new_wavelengths = new_bands,
                               FWHM = rep(7.43, length(new_bands)),
                               threads = threads)
  
  # Remove artifacts of resampling (no nice but works!)
  resampled$'381.005584716796' <- resampled$'395.865584716796'
  resampled$'388.435584716796' <- resampled$'395.865584716796'
  resampled$'2483.6955847168' <- resampled$'2476.2655847168'
  resampled$'2491.1255847168' <- resampled$'2476.2655847168'
  
  # Apply CWT transformation
  transformed <- cwt(t = resampled,
                     scales = c(2, 3),
                     variance = 1,
                     summed_wavelet = TRUE,
                     threads = threads)
  
  transformed <- cbind(cell = cell, transformed) 
  transformed <- transformed[order(cell)]
  
  # Return transformation
  return(transformed)

}
