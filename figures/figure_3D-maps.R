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
render_snapshot(filename = "DTM_3d.png", 
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

rgb <- c(refl$`648.730285644531`, refl$`551.86669921875`, refl$`447.742797851562`)
names(rgb) <- c("RED", "GREEN", "BLUE")

rgb <- mask(crop(rgb, aoi), aoi)

# Normalize per band using global max values
max_vals <- as.numeric(global(rgb, fun = max, na.rm = TRUE)$max)
min_vals <- as.numeric(global(rgb, fun = min, na.rm = TRUE)$min)

for (i in 1:nlyr(rgb)) {
  rgb[[i]] <- (rgb[[i]] - min_vals[i]) / (max_vals[i] - min_vals[i])
}

rgb_3d <- as.array(rgb)
rgb_3d <- aperm(rgb_3d, c(2, 1, 3))
rgb_3d <- rgb_3d / max(rgb_3d)
writePNG(rgb_3d, target = "reflectance.png")
rgb_3d <- png::readPNG("reflectance.png")

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
  height_shade(texture = c("black", viridis(255, alpha = 1, begin = 0, end = 1, direction = 1, option = "D"))) %>%
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
Nitrogen_3d <- stretch_percentile(Nitrogen, values = c(1250, 1900))
Nitrogen_3d <- raster_to_matrix(Nitrogen_3d)

Nitrogen_3d %>%
  height_shade(texture = c("black", viridis(255, alpha = 1, begin = 0, end = 1, direction = 1, option = "D"))) %>%
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
  height_shade(texture = c("black", viridis(255, alpha = 1, begin = 0, end = 1, direction = 1, option = "D"))) %>%
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
