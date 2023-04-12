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

#'------------------------------------------------------------------------------
#' Arguments
#' @param spectra A frame of spectra where rows are samples and columns bands

#'------------------------------------------------------------------------------
#' @example 

#'------------------------------------------------------------------------------
#' Functions

# Vector normalization
vector_normalization <- function(spectra) {
  
  norm <- sqrt(apply(spectra^2, 1, sum, na.rm = TRUE))
  normalized <- spectra/norm
  return(normalized)

}

# Savitzky-Golay smoothing
smooth <- function(spectra) {
  
  smoothed <- spectra
  
  for(s in 1:nrow(spectra)) {
    sample <- as.numeric(spectra[s, ])
    smoothed[s, ] <- sgolayfilt(sample, p = 1, n = 5)
  }
  
  return(smoothed)
}

# Continous Wavelet Transform
CWT <- function(spectra, sampling = 2.55, scales = 5) {
  
  #Copy
  transformed <- data.table()
  
  #Bands and new bands
  bands <- as.numeric(colnames(spectra))
  new_bands <- seq(min(bands), max(bands), by = sampling)
  
  #Loop over samples
  for(t in 1:nrow(spectra)) {
    
    # Get the spectra
    sample <- as.numeric(spectra[t, ]) 
    
    # Cubic resample
    resample <- spline(x = bands, 
                       y = sample, 
                       xout = new_bands) 
    
    # Create time series
    ts_resample <- ts(resample$y, 
                      start = min(new_bands), 
                      deltat = sampling)
    
    # Transformation
    cwt <- wavCWT(ts_resample, 
                  scale.range= c(sampling, sampling*scales), 
                  variance= 1, 
                  wavelet= "gaussian2") 
    
    #Summed wavelet
    summed_wavelet <- matrix(rowSums(as.matrix(cwt)[,c(3, 4, 5)]), nrow = 1)
    transformed <- rbind(transformed, summed_wavelet)
    
  }
  
  colnames(transformed) <- as.character(new_bands)
  final_bands <- as.character(new_bands[13:(length(new_bands)-12)])
  transformed <- transformed[, ..final_bands]
  
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