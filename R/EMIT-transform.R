#'------------------------------------------------------------------------------
#' @title Spectral transformation for EMIT
#'------------------------------------------------------------------------------

#' @description Transform EMIT imagery using summed wavelet.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(wmtsa)
library(terra)
library(parallel)

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
  
  #--------------------------------------------------------
  #Function for batch transformation
  application <- function(i, 
                          spectra = spectra,
                          bands = EMIT_bands,
                          new_bands = new_bands,
                          cell = cell) {
    
    # Sample
    sample <- as.numeric(spectra[i, ]) 
    
    if(max(sample) != 0) {
      
      # Constant values
      scales <- c(2, 3)
      sampling <- 7.43
      
      # Cubic resample
      resample <- prospectr::resample2(X = sample,
                                       wav = bands,
                                       new.wav = new_bands,
                                       fwhm = sampling)
      
      # Create time series
      ts_resample <- ts(resample, 
                        start = min(new_bands), 
                        deltat = sampling)
      
      # Transformation
      cwt <- wavCWT(ts_resample, 
                    scale.range= c(sampling, sampling*max(scales)), 
                    variance= 1, 
                    wavelet= "gaussian2")
      
      summed_wavelet <- c(cell[i], rowSums(as.matrix(cwt)[, scales]))
      
    } else {
      
      summed_wavelet <- c(cell[i], rep(NA, 285))
      
    }
    
    return(summed_wavelet)
    
  }
  
  # Parallel
  transformed <- mclapply(1:nrow(spectra), 
                          FUN = application,
                          spectra = spectra,
                          bands = EMIT_bands,
                          cell = cell,
                          new_bands = new_bands,
                          mc.cores = threads,
                          mc.preschedule = TRUE,
                          mc.cleanup = TRUE)

  transformed <- data.table(do.call(rbind, transformed))  
  colnames(transformed) <- c("cell", as.character(new_bands)) 
  transformed <- transformed[order(cell)]
  
  return(transformed)
  
}
