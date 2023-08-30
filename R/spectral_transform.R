#'------------------------------------------------------------------------------
#' @tile Auxiliary functions for spectra transformation
#'------------------------------------------------------------------------------

#' @description Transform the spectra in different ways.
#' 
#' @return A frame with the transformed spectra.

#'------------------------------------------------------------------------------
#' Library
library(signal)
library(wmtsa)
library(prospectr)

#install.packages("wmtsa", repos="http://R-Forge.R-project.org")
#remotes::install_github("spkaluzny/ifultools")

#'------------------------------------------------------------------------------
#' Arguments
#' @param spectra A frame of spectra where rows are samples and columns bands

#'------------------------------------------------------------------------------
#' @RDS_load

DESIS_bands <- readRDS("data/DESIS_bands.rds")
DESIS_mask <- readRDS("data/DESIS_mask.rds")
EMIT_bands <- readRDS("data/EMIT_bands.rds")
EMIT_mask <- readRDS("data/EMIT_mask.rds")
EMIT_good_bands <- readRDS("data/EMIT_good_bands.rds")

#'------------------------------------------------------------------------------
#' Functions

# Vector normalization
vector_normalization <- function(spectra) {
  
  norm <- sqrt(apply(spectra^2, 1, sum, na.rm = TRUE))
  normalized <- spectra/norm
  return(normalized)

}

# Savitzky-Golay smoothing
smooth <- function(spectra, p = 1, n = 5) {
  
  smoothed <- data.table()
  
  for(s in 1:nrow(spectra)) {
    sample <- as.numeric(spectra[s, ])
    smoothed <- rbind(smoothed, matrix(sgolayfilt(sample, p = p, n = n), nrow = 1))
  }
  
  colnames(smoothed) <- colnames(spectra)
  return(smoothed)
}

# Continuous Wavelet Transform
CWT <- function(spectra, 
                sampling = 3.0, 
                scales = c(3, 4),
                sensor_bands = NULL) {
  
  #Copy
  transformed <- data.table()
  
  #Bands and new bands
  if(is.null(sensor_bands) == TRUE) {
    bands <- as.numeric(colnames(spectra))
    new_bands <- seq(min(bands), max(bands), by = sampling)
  } else {
    bands <- sensor_bands
    new_bands <- seq(min(bands), max(bands), by = sampling)
  }
  
  #Loop over samples
  for(t in 1:nrow(spectra)) {
    
    # Get the spectra
    sample <- as.numeric(spectra[t, ]) 
    
    # Resample
    resample <- prospectr::resample2(X = sample,
                                     wav = bands,
                                     new.wav = new_bands,
                                     fwhm = sampling)
    
    # Time series
    ts_resample <- ts(resample, 
                      start = min(new_bands), 
                      deltat = sampling)
    
    # Transformation
    cwt <- wavCWT(ts_resample, 
                  scale.range= c(sampling, sampling*max(scales)), 
                  variance= 1, 
                  wavelet= "gaussian2")
    
    if(length(scales) == 1) {
      #A single scales
      wavelet <- matrix(as.matrix(cwt)[,c(scales)], nrow = 1)
      transformed <- rbind(transformed, wavelet)
      
    } else if(length(scales) > 1) {
      #Summed wavelet spectra
      summed_wavelet <- matrix(rowSums(as.matrix(cwt)[,c(scales)]), nrow = 1)
      transformed <- rbind(transformed, summed_wavelet)
    }
  }
  
  colnames(transformed) <- as.character(new_bands)
  return(transformed)
  
}

# Mean extraction 
mean_extraction <- function(spectra) {
  
  # Get values
  mean_extracted <- apply(spectra, 1, mean, na.rm = TRUE)
  sd_extracted <- apply(spectra, 1, sd, na.rm = TRUE)
  
  # Copy
  mean_extracted <- (spectra-mean_extracted)/sd_extracted
  
  # Export 
  return(mean_extracted)
  
}

# Normalization
spectra_normalization <- function(spectra) {
  
  # Get values
  min_extracted <- apply(spectra, 1, min, na.rm = TRUE)
  max_extracted <- apply(spectra, 1, max, na.rm = TRUE)
  
  # Copy
  normalized <- (spectra - min_extracted)/(max_extracted - min_extracted)
  
  # Export 
  return(normalized)
  
}

# Linear interpolation
interpolation <- function(spectra)  {
  
  # Copy
  interpolated <- data.table()
  
  # Bands
  bands <- as.numeric(colnames(spectra))
  
  #Loop over samples
  for(t in 1:nrow(spectra)) {
    
    # Get the spectra
    int_sample <- matrix(approx(x = bands, 
                                y = as.numeric(spectra[t, ]), 
                                xout = bands, 
                                method = "linear", 
                                n = 1, rule = 2, f = 0)$y, nrow = 1)
    interpolated <- rbind(interpolated, int_sample)
    
  }
  
  colnames(interpolated) <- as.character(bands)
  return(interpolated)
  
}

