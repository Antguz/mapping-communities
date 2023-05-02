#'------------------------------------------------------------------------------
#' @tile Pixel Subtraction
#'------------------------------------------------------------------------------

#' @description Subtract pixels from Geotiffs of DESIS L2. During the extraction
#' use a file of coordinates and plot labels (PLT_CN, "LON", "LAT") and transform
#' the extracted bands using three methods.
#' 
#' @return several .txt files with names according the DESIS file of pixels 
#' transformed.

#'------------------------------------------------------------------------------
#' Library

library(terra)
library(data.table)
library(utils)
library(parallel)
options(scipen = 100000)

#'------------------------------------------------------------------------------
#' Load source code

source("R/spectral_transform.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param scene_path Path to the root folder where the scenes are located (compressed).
#' @param FIA_path Path where the FIA plot file with labels is located.
#' @param temp Path to a folder to storage temporal information during extraction.
#' @param out_path Path to export the extracted and transformed pixels.
#' @param threads An integer with the number of threads to use (no more than CPUs).

#'------------------------------------------------------------------------------
#' @example 

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
folder <- "Appalachians"
scene_path <- paste0("/media/antonio/Extreme SSD/Imagery", "/", folder)
temp <- paste0("/media/antonio/Extreme SSD/Imagery", "/temp")
FIA_path <- paste0(root_path, "/00_FIA/plots_eastern_20230427.csv")
out_path <- paste0("/media/antonio/Extreme SSD/Imagery", "/01_pixel-subtraction/", folder)
threads <- 5

rast_subst(scene_path, FIA_path, out_path, threads)

#'------------------------------------------------------------------------------
#' Function
rast_subst <- function(scene_path, FIA_path, out_path, threads) {
  
  ### Vector file processing 
  # Read the FI
  FIA <- fread(FIA_path)
  FIA <- FIA[, c("PLT_CN", "LAT", "LON")]
  
  # Create a virtual vector file
  FIA_vector <- vect(x = FIA, 
                     geom=c("LON", "LAT"), 
                     crs = "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs", 
                     keepgeom = FALSE)
  
  ### Raster files processing
  #List files
  files <- list.files(path = scene_path, 
                      pattern = ".zip", 
                      recursive = TRUE,
                      include.dirs = TRUE)
  
  # Arrange path in a frame
  frame <- data.table(matrix(unlist(strsplit(files, ".zip")), 
                             nrow= length(files), 
                             byrow=TRUE), 
                      stringsAsFactors=FALSE)
  colnames(frame) <- c("scene")
  
  # Unique scenes
  n_scenes <- nrow(frame)
  
  #Function to apply
  extraction <- function(X, 
                         frame = frame, 
                         scene_path = scene_path, 
                         out_path = out_path, 
                         temp = temp, 
                         FIA_vector = FIA_vector, 
                         DESIS_bands = DESIS_bands) {
    
    scene <- frame$scene[X]
    
    #Unzip files
    zip_scene <- paste0(scene_path, "/", scene, ".zip")
    out_zip <- paste0(temp, "/", scene)
    unzip(zip_scene, exdir = out_zip)

    # Select products
    spectra <- paste0(out_zip, "/", 
                      scene, "-", 
                      "SPECTRAL_IMAGE.tif")
    
    quality_cube <- paste0(out_zip, "/", 
                           scene, "-", 
                           "QL_QUALITY.tif")
    
    quality_flags <- paste0(out_zip, "/", 
                            scene, "-", 
                            "QL_QUALITY-2.tif")
    
    # Read products
    rast_spectra <- rast(spectra)
    rast_quality_cube <- rast(quality_cube)
    rast_quality_flags <- rast(quality_flags)
    
    # Test for projections
    if(crs(FIA_vector) != crs(rast_quality_flags)) {
      points <- project(FIA_vector, rast_quality_flags)
    } else {
      points <- FIA_vector
    }
    
    #Crop
    e <- ext(rast_quality_flags)
    points <- crop(points, e)
    
    # If there is overlap
    if(length(points) != 0) {
      
      # Pixel extraction
      point_spectra <- as.data.frame(extract(rast_spectra, points, cells=FALSE, xy=FALSE))
      point_quality_cube <- as.data.frame(extract(rast_quality_cube, points, cells=FALSE, xy=FALSE))
      point_quality_flags <-as.data.frame(extract(rast_quality_flags, points, cells=FALSE, xy=FALSE))
      
      # Colnames
      colnames(point_spectra) <- c("PLT_CN", DESIS_bands)
      colnames(point_quality_cube) <- c("PLT_CN", DESIS_bands)
      
      #Change ID
      point_spectra$PLT_CN <- points$PLT_CN
      point_quality_cube$PLT_CN <- points$PLT_CN
      colnames(point_quality_flags)[1] <- "PLT_CN"
      point_quality_flags$PLT_CN <- points$PLT_CN
      
      # To remove
      remove <- NA
      
      # Test values for bands
      for(j in 1:nrow(point_spectra)) {
        for(k in 2:ncol(point_spectra)) {
          if(is.na(point_quality_cube[j, k])) {
            remove[j] <- TRUE
          } else if(point_quality_cube[j, k] == 32) {
            point_spectra[j, k] <- NA
          }
        }
      }
      
      # NA observations
      remove[is.na(remove)] <- FALSE
      point_spectra <- point_spectra[remove != TRUE,]
      point_quality_flags <- point_quality_flags[remove != TRUE,]
      
      # If there are points at all
      if(nrow(point_spectra) > 0) {
        
        # Work with the spectra
        spectra <- point_spectra[, 13:207] #429 and 927.1 nm
        
        ### Perform transformations
        VN <- cbind(point_quality_flags, vector_normalization(spectra))
        SG <- cbind(point_quality_flags, smooth(spectra))
        SW <- cbind(point_quality_flags, CWT(spectra, sampling = 2.55, scales = 5))
        #ME <- cbind(point_quality_flags, mean_extraction(spectra))
        #SN <- cbind(point_quality_flags, spectra_normalization(spectra))
        
        ### Export 
        # Base file name
        export_name <- paste0(out_path, "/",
                              scene)
        
        # Export files
        fwrite(VN, paste0(export_name, "_VN.txt"), sep = "\t")
        fwrite(SG, paste0(export_name, "_SG.txt"), sep = "\t")
        fwrite(SW, paste0(export_name, "_SW.txt"), sep = "\t")
        #fwrite(ME, paste0(export_name, "_ME.txt"), sep = "\t")
        #fwrite(SN, paste0(export_name, "_SN.txt"), sep = "\t")
        
      }
    }
    
    #Remove temporal directory
    unlink(out_zip, recursive = TRUE)
    return(X)
    
  }
  
  #cl <- makeForkCluster(threads)
  #clusterExport(cl, varlist = c("rast", "crs", "crop", "ext"))
  #parLapply(cl, 
  #          X = 1:nrow(frame), 
  #          fun = extraction, 
  #          frame = frame, 
  #          scene_path = scene_path, 
  #          out_path = out_path, 
  #          temp = temp, 
  #          FIA_vector = FIA_vector, 
  #          DESIS_bands = DESIS_bands)
  #stopCluster(cl)
  
  #Add progress bar
  pb <- txtProgressBar(min = 1, max = n_scenes, style = 3)
  
  for(i in 1:n_scenes) {
    
    setTxtProgressBar(pb, i)
    
    extraction(X = i,
               frame = frame, 
               scene_path = scene_path, 
               out_path = out_path, 
               temp = temp, 
               FIA_vector = FIA_vector, 
               DESIS_bands = DESIS_bands)
  }
  
}


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

#Auxiliary bands
DESIS_bands <- c(401.53, 404.20, 406.82, 409.33, 411.81, 414.34, 416.91, 419.44, 
                 422.01, 424.64, 427.29, 429.86, 432.45, 434.98, 437.49, 439.95, 
                 442.49, 445.02, 447.76, 450.34, 452.89, 455.57, 458.18, 460.75, 
                 463.28, 465.82, 468.30, 475.87, 478.57, 481.24, 483.80, 486.39, 
                 488.99, 491.55, 493.97, 496.50, 499.24, 501.75, 504.32, 506.92, 
                 509.54, 512.15, 514.69, 517.17, 519.66, 522.19, 524.69, 527.22, 
                 529.76, 532.33, 534.93, 540.15, 542.69, 545.21, 547.76, 550.38, 
                 552.91, 555.56, 558.08, 560.62, 563.18, 565.79, 568.41, 571.00, 
                 573.58, 576.10, 578.60, 581.13, 583.66, 586.23, 588.78, 591.34, 
                 593.96, 596.58, 599.04, 601.61, 604.19, 606.75, 609.26, 611.80, 
                 614.34, 616.90, 619.52, 622.11, 624.71, 627.18, 629.66, 632.22, 
                 639.76, 642.32, 644.95, 647.57, 650.02, 652.59, 655.24, 657.72, 
                 660.22, 662.76, 665.30, 667.85, 670.48, 673.16, 675.80, 680.91, 
                 683.47, 685.87, 688.42, 690.88, 693.47, 696.11, 698.82, 701.41, 
                 703.92, 706.74, 709.37, 711.78, 714.02, 716.36, 718.89, 721.55, 
                 724.21, 726.85, 729.45, 732.03, 734.37, 736.90, 739.50, 741.99, 
                 744.52, 747.13, 749.74, 752.27, 755.17, 757.74, 760.33, 762.96, 
                 764.91, 767.57, 770.33, 772.71, 775.34, 778.00, 780.59, 783.08, 
                 785.58, 788.26, 790.57, 793.17, 795.90, 798.37, 801.16, 804.08, 
                 806.67, 809.15, 811.70, 814.26, 816.88, 819.76, 822.81, 824.26, 
                 827.17, 829.18, 832.20, 834.88, 836.69, 840.07, 841.97, 844.66, 
                 847.73, 849.93, 852.46, 855.37, 857.88, 860.30, 862.87, 865.44, 
                 867.97, 870.59, 873.18, 875.77, 878.76, 881.52, 883.08, 885.31, 
                 888.08, 890.92, 894.05, 895.97, 898.36, 901.20, 903.76, 905.98, 
                 908.63, 911.62, 914.77, 916.56, 918.38, 920.93, 923.84, 927.07, 
                 929.66, 931.85, 934.49, 937.29, 939.40, 941.92, 944.74, 947.28, 
                 949.58, 951.90, 954.22, 957.31, 959.57, 962.30, 965.42, 968.10, 
                 970.43, 972.90, 976.03, 978.62, 980.00, 981.99, 984.84, 988.87, 
                 991.66, 993.15, 995.64, 997.92, 1000.08)
