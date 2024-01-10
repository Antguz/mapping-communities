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

#' @param scene_main_path A path to the scene folder
#' @param EMIT_bands A vector with the bands of EMIT
#' @param threads N threads to run in parallel

#root_path <- "/media/antonio/Work/Projects/EMIT_mapping/00_imagery"
#scene_main_path <- paste0(root_path, "/", "EMIT-HSI-L2A-DT0599346912 015-20210615T161539-V0215")
#output_raster <- paste0(root_path, "/", scene, "_SW.tif")
EMIT_bands <- readRDS("data/EMIT_bands.rds")

#'------------------------------------------------------------------------------
#' Function

EMIT_transform <- function(spectra, 
                            EMIT_bands,
                            cell,
                            new_bands = seq(min(EMIT_bands), max(EMIT_bands), by = 3.0),
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
    
    # Constant values
    scales <- 5
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
                  scale.range= c(sampling, sampling*scales), 
                  variance= 1, 
                  wavelet= "gaussian2")
    
    summed_wavelet <- c(cell[i], rowSums(as.matrix(cwt)[,c(3, 4, 5)]))
    
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
  colnames(transformed) <- c("cell", as.character(new_bands[3:(length(new_bands)-1)])) 
  transformed <- transformed[order(cell)]
  
  return(transformed)
  
}
