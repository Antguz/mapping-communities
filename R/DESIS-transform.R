#'------------------------------------------------------------------------------
#' @title Spectral transformation for DESIS
#'------------------------------------------------------------------------------

#' @description Transform DESIS imagery using summed wavelet.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(wmtsa)
library(terra)
library(parallel)


#'------------------------------------------------------------------------------
#' Arguments

#' @param scene_main_path A path to the scene folder
#' @param DESIS_bands A vector with the bands of DESIS
#' @param threads N threads to run in parallel

#root_path <- "/media/antonio/Work/Projects/DESIS_mapping/00_imagery"
#scene_main_path <- paste0(root_path, "/", "DESIS-HSI-L2A-DT0599346912 015-20210615T161539-V0215")
#output_raster <- paste0(root_path, "/", scene, "_SW.tif")
DESIS_bands <- readRDS("data/DESIS_bands.rds")

#'------------------------------------------------------------------------------
#' Function

DESIS_transform <- function(spectra, 
                            DESIS_bands,
                            cell,
                            new_bands = seq(min(DESIS_bands), max(DESIS_bands), by = 3.0),
                            threads = 28) {
  
  #--------------------------------------------------------
  #Function for batch transformation
  application <- function(i, 
                          spectra = spectra,
                          bands = DESIS_bands,
                          new_bands = new_bands,
                          cell = cell) {
    
    # Sample
    sample <- as.numeric(spectra[i, ]) 
    
    if(max(sample) != 0) {
      
      # Constant values
      scales <- c(3, 4, 5)
      sampling <- 3.0
      
      # Cubic resample
      resample <- prospectr::resample2(X = sample,
                                       wav = bands,
                                       new.wav = new_bands,
                                       fwhm = sampling)
      
      resample <- resample[3:(length(resample)-1)]
      
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
      
      summed_wavelet <- c(cell[i], rep(NA, 197))
      
    }
    
    return(summed_wavelet)
    
  }
  
  # Parallel
  transformed <- mclapply(1:nrow(spectra), 
                          FUN = application,
                          spectra = spectra,
                          bands = DESIS_bands,
                          cell = cell,
                          new_bands = new_bands,
                          mc.cores = threads,
                          mc.preschedule = TRUE,
                          mc.cleanup = TRUE)
  
  transformed <- data.table(do.call(rbind, transformed))  
  colnames(transformed) <- c("cell", as.character(new_bands[3:(length(new_bands)-1)])) 
  transformed <- transformed[order(cell)]
  
  return(transformed)
  
}
