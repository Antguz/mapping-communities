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

# Continuous Wavelet Transform
CWT <- function(spectra, sampling = 3.0, scales = 5) {
  
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
    
    #Summed wavelet
    summed_wavelet <- matrix(rowSums(as.matrix(cwt)[,c(3, 4, 5)]), nrow = 1)
    transformed <- rbind(transformed, summed_wavelet)
    
  }
  
  colnames(transformed) <- as.character(new_bands[3:(length(new_bands)-1)])
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

DESIS_bands <- c(401.53,
                 404.2,
                 406.82,
                 409.33,
                 411.81,
                 414.34,
                 416.91,
                 419.44,
                 422.01,
                 424.64,
                 427.29,
                 429.86,
                 432.45,
                 434.98,
                 437.49,
                 439.95,
                 442.49,
                 445.02,
                 447.76,
                 450.34,
                 452.89,
                 455.57,
                 458.18,
                 460.75,
                 463.28,
                 465.82,
                 468.3,
                 470.77,
                 473.36,
                 475.87,
                 478.57,
                 481.24,
                 483.8,
                 486.39,
                 488.99,
                 491.55,
                 493.97,
                 496.5,
                 499.24,
                 501.75,
                 504.32,
                 506.92,
                 509.54,
                 512.15,
                 514.69,
                 517.17,
                 519.66,
                 522.19,
                 524.69,
                 527.22,
                 529.76,
                 532.33,
                 534.93,
                 537.56,
                 540.15,
                 542.69,
                 545.21,
                 547.76,
                 550.38,
                 552.91,
                 555.56,
                 558.08,
                 560.62,
                 563.18,
                 565.79,
                 568.41,
                 571.0,
                 573.58,
                 576.1,
                 578.6,
                 581.13,
                 583.66,
                 586.23,
                 588.78,
                 591.34,
                 593.96,
                 596.58,
                 599.04,
                 601.61,
                 604.19,
                 606.75,
                 609.26,
                 611.8,
                 614.34,
                 616.9,
                 619.52,
                 622.11,
                 624.71,
                 627.18,
                 629.66,
                 632.22,
                 634.72,
                 637.24,
                 639.76,
                 642.32,
                 644.95,
                 647.57,
                 650.02,
                 652.59,
                 655.24,
                 657.72,
                 660.22,
                 662.76,
                 665.3,
                 667.85,
                 670.48,
                 673.16,
                 675.8,
                 678.36,
                 680.91,
                 683.47,
                 685.87,
                 688.42,
                 690.88,
                 693.47,
                 696.11,
                 698.82,
                 701.41,
                 703.92,
                 706.74,
                 709.37,
                 711.78,
                 714.02,
                 716.36,
                 718.89,
                 721.55,
                 724.21,
                 726.85,
                 729.45,
                 732.03,
                 734.37,
                 736.9,
                 739.5,
                 741.99,
                 744.52,
                 747.13,
                 749.74,
                 752.27,
                 755.17,
                 757.74,
                 760.33,
                 762.96,
                 764.91,
                 767.57,
                 770.33,
                 772.71,
                 775.34,
                 778.0,
                 780.59,
                 783.08,
                 785.58,
                 788.26,
                 790.57,
                 793.17,
                 795.9,
                 798.37,
                 801.16,
                 804.08,
                 806.67,
                 809.15,
                 811.7,
                 814.26,
                 816.88,
                 819.76,
                 822.81,
                 824.26,
                 827.17,
                 829.18,
                 832.2,
                 834.88,
                 836.69,
                 840.07,
                 841.97,
                 844.66,
                 847.73,
                 849.93,
                 852.46,
                 855.37,
                 857.88,
                 860.3,
                 862.87,
                 865.44,
                 867.97,
                 870.59,
                 873.18,
                 875.77,
                 878.76,
                 881.52,
                 883.08,
                 885.31,
                 888.08,
                 890.92,
                 894.05,
                 895.97,
                 898.36,
                 901.2,
                 903.76,
                 905.98,
                 908.63,
                 911.62,
                 914.77,
                 916.56,
                 918.38,
                 920.93,
                 923.84,
                 927.07,
                 929.66,
                 931.85,
                 934.49,
                 937.29,
                 939.4,
                 941.92,
                 944.74,
                 947.28,
                 949.58,
                 951.9,
                 954.22,
                 957.31,
                 959.57,
                 962.3,
                 965.42,
                 968.1,
                 970.43,
                 972.9,
                 976.03,
                 978.62,
                 980.0,
                 981.99,
                 984.84,
                 988.87,
                 991.66,
                 993.15,
                 995.64,
                 997.92,
                 1000.08)

fwhm  <- c(
  4.2,
  4.0,
  3.6,
  3.8,
  3.6,
  3.8,
  3.4,
  3.4,
  3.6,
  3.6,
  3.4,
  3.4,
  3.2,
  3.2,
  3.0,
  3.2,
  3.4,
  3.2,
  3.4,
  3.4,
  3.4,
  3.6,
  3.4,
  3.4,
  3.2,
  3.2,
  3.4,
  3.2,
  3.2,
  3.4,
  3.4,
  3.4,
  3.2,
  3.2,
  3.2,
  3.2,
  3.2,
  3.2,
  3.4,
  3.4,
  3.4,
  3.6,
  3.6,
  3.4,
  3.4,
  3.4,
  3.4,
  3.4,
  3.4,
  3.6,
  3.4,
  3.4,
  3.6,
  3.6,
  3.4,
  3.4,
  3.6,
  3.6,
  3.6,
  3.4,
  3.4,
  3.4,
  3.6,
  3.6,
  3.6,
  3.6,
  3.6,
  3.6,
  3.4,
  3.4,
  3.4,
  3.4,
  3.6,
  3.6,
  3.6,
  3.6,
  3.6,
  3.4,
  3.6,
  3.6,
  3.4,
  3.4,
  3.6,
  3.4,
  3.4,
  3.4,
  3.4,
  3.4,
  3.4,
  3.4,
  3.2,
  3.4,
  3.4,
  3.4,
  3.4,
  3.6,
  3.6,
  3.4,
  3.6,
  3.4,
  3.4,
  3.4,
  3.4,
  3.4,
  3.4,
  3.4,
  3.6,
  3.6,
  3.4,
  3.4,
  3.4,
  3.4,
  3.2,
  3.4,
  3.4,
  3.4,
  3.6,
  3.6,
  3.8,
  3.6,
  3.6,
  4.0,
  4.0,
  3.8,
  3.4,
  3.4,
  3.6,
  3.4,
  3.4,
  3.2,
  3.2,
  3.4,
  3.4,
  3.2,
  3.4,
  3.4,
  3.4,
  3.4,
  3.6,
  3.4,
  3.4,
  3.6,
  3.0,
  3.4,
  3.4,
  3.4,
  3.4,
  3.6,
  3.2,
  3.4,
  3.2,
  3.4,
  3.2,
  3.4,
  3.4,
  3.4,
  3.2,
  3.0,
  3.0,
  3.0,
  3.0,
  3.0,
  3.2,
  3.6,
  3.2,
  3.0,
  3.6,
  2.8,
  3.6,
  2.6,
  3.0,
  3.4,
  2.8,
  3.4,
  3.4,
  3.0,
  3.2,
  3.2,
  3.2,
  3.0,
  3.4,
  3.2,
  3.4,
  3.4,
  3.4,
  3.6,
  4.0,
  3.0,
  3.0,
  3.4,
  3.6,
  3.6,
  3.6,
  3.2,
  3.2,
  3.6,
  3.6,
  3.0,
  4.0,
  4.0,
  3.8,
  3.0,
  3.6,
  4.0,
  3.6,
  3.4,
  3.4,
  3.2,
  3.8,
  3.4,
  3.2,
  3.8,
  3.6,
  3.6,
  3.6,
  3.0,
  3.6,
  3.4,
  3.4,
  3.4,
  3.8,
  3.8,
  3.6,
  3.4,
  4.8,
  3.6,
  2.8,
  4.0,
  6.6,
  4.2,
  3.4,
  3.0,
  3.4,
  3.4,
  3.0)
