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
options(scipen = 100000)

#'------------------------------------------------------------------------------
#' Load source code

source("R/spectral_transform.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param scene_path Path to the root folder where the scenes are located
#' @param FIA_path Path where the FIA plot file with labels is located 
#' @param out_path Path to export the transformed pixels.

#'------------------------------------------------------------------------------
#' @example 

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
scene_path <- paste0(root_path, "/00_imagery")
FIA_path <- paste0(root_path, "/00_FIA/plots_WI_MN_20221102.csv")
out_path <- paste0(root_path, "/01_pixel-subtraction")

rast_subst(scene_path, FIA_path, out_path)

#'------------------------------------------------------------------------------
#' Function
rast_subst <- function(scene_path, FIA_path, out_path) {
  
  ### Vector file processing 
  # Read the FI
  FIA <- fread(FIA_path)
  FIA <- subset(FIA, STATECD == 27 | STATECD == 55) #Select specific states
  FIA <- FIA[, c("PLT_CN", "LAT", "LON")]
  
  # Create a virtual vector file
  FIA_vector <- vect(x = FIA, 
                     geom=c("LON", "LAT"), 
                     crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs", 
                     keepgeom = FALSE)
  
  ### Raster files processing
  #List files
  files <- list.files(path = scene_path, 
                      pattern = ".tif", 
                      recursive = TRUE,
                      include.dirs = TRUE)
  
  # Arrange path in a frame
  frame <- data.table(matrix(unlist(strsplit(files, "/")), 
                             nrow= length(files), 
                             byrow=TRUE), 
                      stringsAsFactors=FALSE)
  colnames(frame) <- c("folder", "scene")
  
  # Extract Product type
  frame[, type := unlist(strsplit(scene, "V0215-"))[2], 
        by = seq_along(1:nrow(frame))]
  frame[, type := unlist(strsplit(type, ".tif"))[1], 
        by = seq_along(1:nrow(frame))]
  
  # Unique scenes
  unique_scenes <- frame[, .N, by= c("folder")]
  n_scenes <- nrow(unique_scenes)
  unique_scenes$N <- 1:n_scenes
  
  # Add unique scenes to frame
  frame <- merge(unique_scenes, frame, by = c("folder"))
  
  #Add progress bar
  pb <- txtProgressBar(min = 1, max = n_scenes, style = 3)
  
  ### Loop for pixel extraction 
  for(i in 1:n_scenes) {
    
    # Track progress
    setTxtProgressBar(pb, i)
    
    # Select scene of interest
    sub_frame <- frame[N == i]
    
    # Select products
    spectra <- paste0(scene_path, "/", 
                      sub_frame$folder[1], "/", 
                      sub_frame[type == "SPECTRAL_IMAGE"]$scene[1])
    
    quality_cube <- paste0(scene_path, "/", 
                           sub_frame$folder[1], "/", 
                           sub_frame[type == "QL_QUALITY"]$scene[1])
    
    quality_flags <- paste0(scene_path, "/", 
                            sub_frame$folder[1], "/",
                            sub_frame[type == "QL_QUALITY-2"]$scene[1])
    
    # Read products
    rast_spectra <- rast(spectra)
    rast_quality_cube <- rast(quality_cube)
    rast_quality_flags <- rast(quality_flags)
    
    # Test for projections
    if(crs(FIA_vector) != crs(rast_quality_flags)) {
      FIA_vector <- project(FIA_vector, rast_quality_flags)
    }
    
    # Crop points by extend
    points <- FIA_vector 
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
        ME <- cbind(point_quality_flags, mean_extraction(spectra))
        SN <- cbind(point_quality_flags, spectra_normalization(spectra))
        
        ### Export 
        # Base file name
        export_name <- paste0(out_path, "/",
                              sub_frame$folder[1], "/",
                              unlist(strsplit(sub_frame[type == "QL_QUALITY-2"]$product[1], "-V0215-"))[1])
        
        # Export files
        fwrite(VN, paste0(export_name, "_VN.txt"), sep = "\t")
        fwrite(SG, paste0(export_name, "_SG.txt"), sep = "\t")
        fwrite(SW, paste0(export_name, "_SW.txt"), sep = "\t")
        fwrite(ME, paste0(export_name, "_ME.txt"), sep = "\t")
        fwrite(SN, paste0(export_name, "_SN.txt"), sep = "\t")
        
        
      }
    }
  }
  
  #Close progress
  close(pb)
  
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
