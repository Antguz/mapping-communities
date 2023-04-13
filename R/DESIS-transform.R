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

#'------------------------------------------------------------------------------
#' Function

desis_transform <- function(spectra, 
                            DESIS_bands,
                            cell,
                            threads = 28) {
  
  #--------------------------------------------------------
  #Function for batch transformation
  application <- function(i, 
                          spectra = spectra,
                          bands = DESIS_bands,
                          new_bands = seq(min(DESIS_bands), max(DESIS_bands), by = 2.55),
                          cell = cell) {
    
    # Sample
    sample <- as.numeric(spectra[i, ]) 
    
    # Constant values
    scales = 5
    sampling = 2.55
    
    # Cubic resample
    resample <- spline(x = DESIS_bands, 
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
    
    summed_wavelet <- c(cell[i], rowSums(as.matrix(cwt)[,c(3, 4, 5)]))
    
    return(summed_wavelet)
    
  }
  
  # Parallel
  transformed <- mclapply(1:nrow(spectra), 
                          FUN = application,
                          spectra = spectra,
                          bands = DESIS_bands,
                          cell = cell,
                          new_bands = seq(min(DESIS_bands), max(DESIS_bands), by = 2.55),
                          mc.cores = threads,
                          mc.preschedule = TRUE,
                          mc.cleanup = TRUE)
  
  transformed <- data.table(do.call(rbind, transformed))  
  colnames(transformed) <- c("cell", seq(min(DESIS_bands), 
                                        max(DESIS_bands), 
                                        by = 2.55)) 
  transformed <- transformed[order(cell)]
  
  return(transformed)
  
}

#Auxiliary bands
DESIS_bands <- c(401.9, 404.1, 406.7, 409.3, 411.7, 414.3, 416.8, 419.4, 
                 421.9, 424.5, 427.2, 429.8, 432.4, 434.9, 437.4, 439.9,
                 442.4, 444.9, 447.7, 450.3, 452.8, 455.5, 458,   460.7, 
                 463.1, 465.7, 468.2, 470.7, 473.3, 475.8, 478.5, 481.1,
                 483.7, 486.3, 488.9, 491.4, 493.8, 496.4, 499.2, 501.6, 
                 504.2, 506.8, 509.5, 512,   514.6, 517,   519.6, 522,  
                 524.6, 527.1, 529.7, 532.2, 534.8, 537.5, 540,   542.5, 
                 545.1, 547.7, 550.3, 552.9, 555.5, 558,   560.5, 563,  
                 565.7, 568.3, 570.9, 573.5, 576,   578.5, 581.1, 583.5, 
                 586.1, 588.7, 591.3, 593.9, 596.5, 598.9, 601.5, 604.1,
                 606.6, 609.2, 611.7, 614.3, 616.8, 619.4, 622,   624.6, 
                 627,   629.6, 632.1, 634.6, 637.2, 639.7, 642.2, 644.9,
                 647.5, 649.9, 652.5, 655.1, 657.6, 660.1, 662.6, 665.2, 
                 667.7, 670.3, 673.1, 675.7, 678.3, 680.8, 683.4, 685.8,
                 688.3, 690.8, 693.4, 696,   698.8, 701.3, 703.9, 706.7, 
                 709.3, 711.6, 713.9, 716.3, 718.8, 721.4, 724.1, 726.7,
                 729.4, 731.9, 734.2, 736.8, 739.4, 741.9, 744.4, 747,   
                 749.7, 752.2, 755,   757.7, 760.3, 762.9, 764.8, 767.5,
                 770.2, 772.6, 775.2, 777.9, 780.5, 783,   785.5, 788.2, 
                 790.4, 793.1, 795.8, 798.2, 801,   804,   806.5, 809.1,
                 811.6, 814.2, 816.8, 819.6, 822.7, 824.2, 827.1, 829.1, 
                 832.1, 834.8, 836.6, 839.9, 841.9, 844.6, 847.7, 849.8,
                 852.4, 855.2, 857.8, 860.2, 862.8, 865.4, 867.9, 870.5, 
                 873.1, 875.6, 878.6, 881.4, 883,   885.3, 888,   890.8,
                 893.9, 895.9, 898.3, 901.1, 903.7, 905.8, 908.5, 911.5, 
                 914.7, 916.4, 918.3, 920.8, 923.7, 927,   929.6, 931.7,
                 934.4, 937.2, 939.3, 941.8, 944.7, 947.2, 949.4, 951.8, 
                 954.2, 957.2, 959.5, 962.2, 965.3, 968.1, 970.3, 972.9,
                 975.9, 978.5, 979.9, 981.9, 984.8, 988.8, 991.6, 993.1, 
                 995.6, 997.8, 999.5)
