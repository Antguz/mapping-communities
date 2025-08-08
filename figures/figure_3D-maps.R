#'------------------------------------------------------------------------------
#' @title 3D topographic image for mapping biodiversity
#'------------------------------------------------------------------------------

#' @description 

#' @caution Run in Linux base machine due to parallel needs

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(terra)
library(rayshader)
library(rgl)
library(viridis)
library(png)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"
source("figures/RGB.R")

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# DTM plot ---------------------------------------------------------------------

# Area of interest
aoi <- vect(paste0(root_path, "/figures/Maps/limits.gpkg"))

# Read DTM data
dtm <- rast(paste0(root_path, "/figures/Maps/DTM/DEM.tif"))

# Resample DTM to match hyperspectral image
dtm_resampled <- resample(dtm, 
                          rast(paste0(root_path, "/figures/Maps/Rflec/reflectance.tif")), 
                          method = "bilinear")

# Crop 
dtm <- mask(crop(dtm_resampled, aoi), aoi)
dtm_3d <- raster_to_matrix(dtm)

# Initial ray-tracing
zscale_value <- 12
rayshade <- ray_shade(dtm_3d, 
                      zscale = zscale_value,
                      multicore = TRUE)

ambientshade <- ambient_shade(dtm_3d, 
                              zscale = zscale_value,
                              multicore = TRUE)

# Digital elevation model

dtm_3d %>%
  height_shade(texture = viridis(256, alpha = 1, begin = 0, end = 1, direction = 1, option = "turbo")) %>%
  add_shadow(rayshade, 0.5) %>%
  add_shadow(ambientshade, 0.8) %>%
  plot_3d(heightmap = dtm_3d,
          zscale = zscale_value,
          theta = 143.5, #323
          phi = 45,
          fov = 45,
          windowsize = c(640, 480), #c(6400, 3000), #c(1600, 750), #
          background = "transparent",
          solidcolor = "grey20",
          solidlinecolor = "grey30",
          zoom = 1, #
          solid = TRUE)
render_camera(zoom = 0.65)
render_snapshot(filename = paste0(outpath, "/DTM_3d.png"), 
                clear = FALSE, 
                width = 3480,#round(6400), 
                height = 2160,#round(3000),
                software_render = TRUE,
                parallel = TRUE)

# Stretch function -------------------------------------------------------------

stretch_percentile <- function(raster, percentile = c(0.01, 0.99), values = NULL, DTM = dtm) {
  
  if(is.null(values) == TRUE) {
    
    # Compute lower and upper percentiles
    p <- quantile(raster[], probs = c(percentile[1], percentile[2]), na.rm = TRUE)
    p_low <- p[1]
    p_high <- p[2]
    cat(c(p_low, p_high))
    
  } else {
    p_low <- values[1]
    p_high <- values[2]
  }
  
  # Clamp and apply stretching
  r_stretched <- clamp(raster, lower = p_low, upper = p_high)
  r_stretched <- ((r_stretched - p_low) / (p_high - p_low) * 254) + 1
  
  # Fill NA with 0
  dtm_mask <- !is.na(DTM)
  na_overlap <- is.na(r_stretched) & dtm_mask
  r_stretched[na_overlap] <- 0

  return(r_stretched)
  
}

# RGB plot ---------------------------------------------------------------------

refl <- rast(paste0(root_path, "/figures/Maps/Rflec/reflectance.tif"))
names(refl) <- readRDS("data/EMIT_bands.rds")
rgb <- c(refl$`648.730285644531`, refl$`551.86669921875`, refl$`447.742797851562`)
names(rgb) <- c("RED", "GREEN", "BLUE")
rgb <- mask(crop(rgb, aoi), aoi)

# Normalize per band using global max values
max_vals <- as.numeric(global(rgb, fun = max, na.rm = TRUE)$max)
min_vals <- as.numeric(global(rgb, fun = min, na.rm = TRUE)$min)

for (i in 1:nlyr(rgb)) {
  rgb[[i]] <- (rgb[[i]] - min_vals[i]) / (max_vals[i] - min_vals[i])
}

RGB_3d_plot <- function(rgb,
                        dtm_3d = get("dtm_3d", envir = globalenv()), 
                        rayshade = get("rayshade", envir = globalenv()), 
                        ambientshade = get("ambientshade", envir = globalenv()),
                        DTM = get("dtm", envir = globalenv()),
                        outpath = get("outpath", envir = globalenv())) {
  
  # Convert to grayscale using luminance weights
  gray <- rgb$RED + rgb$GREEN + rgb$BLUE
  
  # Scale to 0–255 digital number
  gray_dn <- round(gray * 768)
  
  # Clamp to make sure values are within 0–255
  gray_dn <- clamp(gray_dn, lower=0, upper=768)
  gray_dn_3d <- raster_to_matrix(gray_dn)
  
  # # Create RGB quantized color cube: 8x8x4 = 256 colors
  reds   <- seq(0, 255, length.out = 8)
  greens <- seq(0, 255, length.out = 8)
  blues  <- seq(0, 255, length.out = 8)
  palette_rgb <- expand.grid(R = reds, G = greens, B = blues)
   
  # # Convert to hex color strings
  rgb_palette <- rgb(palette_rgb$R, palette_rgb$G, palette_rgb$B, maxColorValue = 768)
    
  gray_dn_3d %>%
      height_shade(texture = rgb_palette) %>%
      add_shadow(rayshade, 0.5) %>%
      add_shadow(ambientshade, 0.5) %>%
      plot_map()
      
      plot_3d(heightmap = dtm_3d,
              zscale = 12,
              theta = 143.5,
              phi = 45,
              fov = 45,
              windowsize = c(640, 480), # c(640, 480)
              background = "transparent",
              solidcolor = "grey20",
              solidlinecolor = "grey30",
              zoom = 1, #
              solid = TRUE)
    render_camera(zoom = 0.65)
    out_name <- gsub(" / ", "_", layer_names[i])
    out_name <- gsub(" ", "-", out_name)
    render_snapshot(filename = paste0(outpath, "/", out_name, ".png"), 
                    clear = FALSE, 
                    width = 3480,
                    height = 2160,
                    software_render = TRUE,
                    parallel = TRUE)
    
    gc()

}

fortyp_3d_plot(fortyp)

# RGB plot ---------------------------------------------------------------------

MDS_3d_plot <- function(MDS,
                        dtm_3d = get("dtm_3d", envir = globalenv()), 
                        rayshade = get("rayshade", envir = globalenv()), 
                        ambientshade = get("ambientshade", envir = globalenv()),
                        DTM = get("dtm", envir = globalenv()),
                        type = get("type", envir = globalenv()),
                        outpath = get("outpath", envir = globalenv())) {
  
  # Clamp
  MDS_clamp <- (clamp(MDS, -0.5, 0.5) + 0.5) * 255
  
  # Fill NA
  dtm_mask <- !is.na(DTM)
  na_overlap <- is.na(MDS_clamp[[1]]) & dtm_mask
  MDS_clamp$Axis1[na_overlap] <- 51
  MDS_clamp$Axis2[na_overlap] <- 51
  MDS_clamp$Axis3[na_overlap] <- 51
  
  # Export
  MDS_clamp <- as.array(MDS_clamp) / 255
  writePNG(MDS_clamp, paste0(outpath, "/", type, "_file.png"))
  MDS_raster <- readPNG(paste0(outpath, "/", type, "_file.png"))
  
  # Create RGB
  
  dtm_3d %>%
    sphere_shade(texture ="desert") %>%
    add_overlay(MDS_raster, alphalayer = .99) %>%
    add_shadow(rayshade, 0.5) %>%
    add_shadow(ambientshade, 0.5) %>%
    plot_3d(dtm_3d,
            zscale = 12,
            theta = 143.5,
            phi = 45,
            fov = 45,
            windowsize = c(640, 480), # c(640, 480)
            background = "transparent",
            solidcolor = "grey20",
            solidlinecolor = "grey30",
            zoom = 1, #
            solid = TRUE)
  render_camera(zoom = 0.65)
  render_snapshot(filename = paste0(outpath, "/", type, "_3d.png"), 
                  clear = FALSE, 
                  width = 3480,
                  height = 2160,
                  software_render = TRUE,
                  parallel = TRUE)
  
  gc()
  
}


type <- "taxonomic"
level <- "L3A"
MDS1 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175417_2323212_006_MDS-",type,"-composition.tif"))
MDS2 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175429_2323212_007_MDS-",type,"-composition.tif"))
MDS <- merge(MDS1, MDS2)
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[1:3]]/10000
MDS_3d_plot(fortyp)

type <- "phylogenetic"
level <- "L3B"
MDS1 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175417_2323212_006_MDS-",type,"-composition.tif"))
MDS2 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175429_2323212_007_MDS-",type,"-composition.tif"))
MDS <- merge(MDS1, MDS2)
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[1:3]]/10000
MDS_3d_plot(MDS)

type <- "functional"
level <- "L3C"
MDS1 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175417_2323212_006_MDS-",type,"-composition.tif"))
MDS2 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175429_2323212_007_MDS-",type,"-composition.tif"))
MDS <- merge(MDS1, MDS2)
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[1:3]]/10000
MDS_3d_plot(MDS)

################################################################################
##### FORTYPE plots
################################################################################

fortyp1 <- rast(paste0(root_path, "/figures/Maps/FORTYP/EMIT-L4A_001_20230820T175417_2323212_006_FORTYP_mean.tif"))
fortyp2 <- rast(paste0(root_path, "/figures/Maps/FORTYP/EMIT-L4A_001_20230820T175429_2323212_007_FORTYP_mean.tif"))
fortyp <- merge(fortyp1, fortyp2)
fortyp <- mask(crop(fortyp, aoi), aoi)
outpath <- paste0(root_path, "/figures/Maps/3d_figures")

fortyp_3d_plot <- function(fortyp, 
                           dtm_3d = get("dtm_3d", envir = globalenv()), 
                           rayshade = get("rayshade", envir = globalenv()), 
                           ambientshade = get("ambientshade", envir = globalenv()),
                           DTM = get("dtm", envir = globalenv()),
                           outpath = get("outpath", envir = globalenv())) {
  
  layer_names <- names(fortyp)
  
  for(i in 1:length(layer_names)) {
    
    sub_fortyp <- fortyp[[i]]
    sub_fortyp <- stretch_percentile(sub_fortyp, values = c(0, 10000), DTM = dtm)
    sub_fortyp_3d <- raster_to_matrix(sub_fortyp)
    
    sub_fortyp_3d %>%
      height_shade(texture = c("grey20", viridis(255, alpha = 1, begin = 0, end = 1, direction = 1, option = "A"))) %>%
      add_shadow(rayshade, 0.5) %>%
      add_shadow(ambientshade, 0.5) %>%
      plot_3d(heightmap = dtm_3d,
              zscale = 12,
              theta = 143.5,
              phi = 45,
              fov = 45,
              windowsize = c(640, 480), # c(640, 480)
              background = "transparent",
              solidcolor = "grey20",
              solidlinecolor = "grey30",
              zoom = 1, #
              solid = TRUE)
    render_camera(zoom = 0.65)
    out_name <- gsub(" / ", "_", layer_names[i])
    out_name <- gsub(" ", "-", out_name)
    render_snapshot(filename = paste0(outpath, "/", out_name, ".png"), 
                    clear = FALSE, 
                    width = 3480,
                    height = 2160,
                    software_render = TRUE,
                    parallel = TRUE)
    
    gc()
    
  }
}

fortyp_3d_plot(fortyp)

################################################################################
##### Plant lineage plots
################################################################################

occurrence1 <- rast(paste0(root_path, "/figures/Maps/Occurrence/EMIT-L4B_001_20230820T175417_2323212_006_plant-occurrence_mean.tif"))
occurrence2 <- rast(paste0(root_path, "/figures/Maps/Occurrence/EMIT-L4B_001_20230820T175429_2323212_007_plant-occurrence_mean.tif"))
occurrence <- merge(occurrence1, occurrence2)
occurrence <- mask(crop(occurrence, aoi), aoi)
outpath <- paste0(root_path, "/figures/Maps/3d_figures")

occurrence_3d_plot <- function(occurrence, 
                               dtm_3d = get("dtm_3d", envir = globalenv()), 
                               rayshade = get("rayshade", envir = globalenv()), 
                               ambientshade = get("ambientshade", envir = globalenv()),
                               DTM = get("dtm", envir = globalenv()),
                               outpath = get("outpath", envir = globalenv())) {
  
  layer_names <- names(occurrence)
  
  for(i in 1:length(layer_names)) {
    
    sub_occurrence <- occurrence[[i]]
    sub_occurrence <- stretch_percentile(sub_occurrence, values = c(0, 10000), DTM = dtm)
    sub_occurrence_3d <- raster_to_matrix(sub_occurrence)
    
    sub_occurrence_3d %>%
      height_shade(texture = c("grey20", viridis(255, alpha = 1, begin = 0, end = 1, direction = 1, option = "G"))) %>%
      add_shadow(rayshade, 0.5) %>%
      add_shadow(ambientshade, 0.5) %>%
      plot_3d(heightmap = dtm_3d,
              zscale = 12,
              theta = 143.5,
              phi = 45,
              fov = 45,
              windowsize = c(640, 480), # c(640, 480)
              background = "transparent",
              solidcolor = "grey20",
              solidlinecolor = "grey30",
              zoom = 1, #
              solid = TRUE)
    render_camera(zoom = 0.65)
    out_name <- gsub(" / ", "_", layer_names[i])
    out_name <- gsub(" ", "-", out_name)
    render_snapshot(filename = paste0(outpath, "/", out_name, ".png"), 
                    clear = FALSE, 
                    width = 3480, # 3480
                    height = 2160, #2160
                    software_render = TRUE,
                    parallel = TRUE)
    
    gc()
    
  }
}

occurrence_3d_plot(occurrence)


# CWM plot ---------------------------------------------------------------------

func1 <- rast(paste0(root_path, "/figures/Maps/CWM/EMIT-L4C_001_20230820T175417_2323212_006_CWM_mean.tif"))
func2 <- rast(paste0(root_path, "/figures/Maps/CWM/EMIT-L4C_001_20230820T175429_2323212_007_CWM_mean.tif"))
func <- merge(func1, func2)
func <- mask(crop(func, aoi), aoi)

# Leaf mass per Area
LMA <- func$LMA
stretch_percentile(LMA) # 5707 8761
LMA <- stretch_percentile(LMA, values = c(5000, 9000))
LMA_3d <- raster_to_matrix(LMA)

LMA_3d %>%
  height_shade(texture = c("grey20", viridis(255, alpha = 1, begin = 0, end = 1, direction = 1, option = "D"))) %>%
  add_shadow(rayshade, 0.5) %>%
  add_shadow(ambientshade, 0.5) %>%
  plot_3d(heightmap = dtm_3d,
          zscale = 12,
          theta = 143.5,
          phi = 45,
          fov = 45,
          windowsize = c(640, 480),
          background = "transparent",
          solidcolor = "grey20",
          solidlinecolor = "grey30",
          zoom = 1, #
          solid = TRUE)
render_camera(zoom = 0.65)
render_snapshot(filename = "LMA_3d.png", 
                clear = FALSE, 
                width = 3480,
                height = 2160,
                software_render = TRUE,
                parallel = TRUE)

# Nitrogen
Nitrogen <- func$N
stretch_percentile(Nitrogen) # 5707 8761
Nitrogen_3d <- stretch_percentile(Nitrogen, values = c(900, 1900))
Nitrogen_3d <- raster_to_matrix(Nitrogen_3d)

Nitrogen_3d %>%
  height_shade(texture = c("grey20", viridis(255, alpha = 1, begin = 0, end = 1, direction = 1, option = "D"))) %>%
  add_shadow(rayshade, 0.5) %>%
  add_shadow(ambientshade, 0.5) %>%
  plot_3d(heightmap = dtm_3d,
          zscale = 12,
          theta = 143.5,
          phi = 45,
          fov = 45,
          windowsize = c(640, 480),
          background = "transparent",
          solidcolor = "grey20",
          solidlinecolor = "grey30",
          zoom = 1, #
          solid = TRUE)
render_camera(zoom = 0.65)
render_snapshot(filename = "Nitrogen_3d.png", 
                clear = FALSE, 
                width = 3480,
                height = 2160,
                software_render = TRUE,
                parallel = TRUE)

# Wood density
WD <- func$WD
stretch_percentile(WD) # 5707 8761
WD_3d <- stretch_percentile(WD, values = c(50, 60))
WD_3d <- raster_to_matrix(WD_3d)

WD_3d %>%
  height_shade(texture = c("grey20", viridis(255, alpha = 1, begin = 0, end = 1, direction = 1, option = "D"))) %>%
  add_shadow(rayshade, 0.5) %>%
  add_shadow(ambientshade, 0.5) %>%
  plot_3d(heightmap = dtm_3d,
          zscale = 12,
          theta = 143.5,
          phi = 45,
          fov = 45,
          windowsize = c(640, 480),
          background = "transparent",
          solidcolor = "grey20",
          solidlinecolor = "grey30",
          zoom = 1, #
          solid = TRUE)
render_camera(zoom = 0.65)
render_snapshot(filename = "WD_3d.png", 
                clear = FALSE, 
                width = 3480,
                height = 2160,
                software_render = TRUE,
                parallel = TRUE)
