#'------------------------------------------------------------------------------
#' @tile Auxiliary functions for spectra transformation
#'------------------------------------------------------------------------------

#' @description Transform the spectra in different ways.
#' 
#' @return A frame with the transformed spectra.

#'------------------------------------------------------------------------------
#' Library
library(signal)
library(CWT)
#library(prospectr)

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
                scales = c(2, 3),
                sensor_bands = NULL,
                sensor = "EMIT") {
  
  if(sensor == "EMIT") {
    new_wavelengths <- seq(min(EMIT_bands), max(EMIT_bands), by = sampling)
  } else if(sensor == "DESIS") {
    
  }
  
  # Apply resampling
  resampled <- resampling_FWHM(spectra = spectra,
                               wavelengths = sensor_bands,
                               new_wavelengths = new_wavelengths,
                               FWHM = rep(sampling, length(new_wavelengths)),
                               threads = 1)
  
  if(sensor == "EMIT") {
    resampled$'381.005584716796' <- resampled$'396.005584716796'
    resampled$'388.505584716796' <- resampled$'396.005584716796'
    resampled$'2488.5055847168' <- resampled$'2481.0055847168'
  } else if(sensor == "DESIS") {
    
  }
  
  transformed <- cwt(t = resampled,
                     scales = scales,
                     variance = 1,
                     summed_wavelet = TRUE,
                     threads = 1)
  
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

