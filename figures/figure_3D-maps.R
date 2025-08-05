#'------------------------------------------------------------------------------
#' @title 3D topographic image for mapping biodiversity
#'------------------------------------------------------------------------------

#' @description 

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(terra)
library(rayshader)
library(rgl)
library(viridis)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# 1. Get area of interest
#-------------------
aoi <- vect(paste0(root_path, "/figures/Maps/limits.gpkg"))

# 2. Get reflectance data
#----------------------
# Read
refl <- rast(paste0(root_path, "/figures/Maps/Rflec/reflectance.tif"))
emit_bands <- readRDS("data/EMIT_bands.rds")
names(refl) <- emit_bands

# Get RGB
rgb <- c(refl$`648.730285644531`, refl$`551.86669921875`, refl$`447.742797851562`)
names(rgb) <- c("RED", "GREEN", "BLUE")

# Crop 
rgb <- mask(crop(rgb, aoi), aoi)

# Normalize per band using global max values
max_vals <- as.numeric(global(rgb, fun = max, na.rm = TRUE)$max)
min_vals <- as.numeric(global(rgb, fun = min, na.rm = TRUE)$min)
for (i in 1:nlyr(rgb)) {
  rgb[[i]] <- (rgb[[i]] - min_vals[i]) / (max_vals[i] - min_vals[i])
}

# # Quantile application
# max_vals <- as.numeric(global(rgb, fun = quantile, probs = 0.975, na.rm = TRUE)$X97.5.)
# min_vals <- as.numeric(global(rgb, fun = quantile, probs = 0.025, na.rm = TRUE)$X2.5.)
# for (i in 1:nlyr(rgb)) {
#   rgb[[i]] <- (rgb[[i]] - min_vals[i]) / (max_vals[i] - min_vals[i])
# }

rgb_3d <- as.array(rgb)
rgb_3d <- aperm(rgb_3d, c(2, 1, 3))
rgb_3d <- rgb_3d / max(rgb_3d)
writePNG(rgb_3d, target = "reflectance.png")
rgb_3d <- png::readPNG("reflectance.png")

# 3. Get elevation data
#----------------------
# Read
dtm <- rast(paste0(root_path, "/figures/Maps/DTM/DEM.tif"))

# Resample DTM to match hyperspectral image
dtm_resampled <- resample(dtm, rgb, method = "bilinear")

# Crop 
dtm <- mask(crop(dtm_resampled, aoi), aoi)
dtm_3d <- raster_to_matrix(dtm)

#-------------------------------------------------------------------------------
# Plots
#-------------------------------------------------------------------------------

# Digital elevation model

dtm_3d %>%
  height_shade(texture = viridis(256, alpha = 1, begin = 0, end = 1, direction = 1, option = "turbo")) %>%
  add_shadow(ray_shade(dtm_3d, zscale = 10), 0.3) %>%
  plot_3d(heightmap = dtm_3d,
          zscale = 8,
          theta = 143.5, #323
          phi = 45,
          fov = 45,
          windowsize = c(1000, 500), #c(6400, 3000), #c(1600, 750), #
          background = "transparent",
          solidcolor = "black",
          solidlinecolor = "black",
          zoom = 1, #
          solid = TRUE)
render_camera(zoom = 0.65)
render_snapshot()
render_highquality(filename = "3d_dtm.png", 
                   clear = FALSE, 
                   samples = 128,
                   light = TRUE,
                   lightdirection = 45,
                   lightintensity = 1000,
                   #lightaltitude = 45,
                   width = 1000,#round(6400), 
                   height = 500,#round(3000),
                   parallel = TRUE)


# RGB

dtm_3d %>%
  add_overlay(rgb_3d) %>%
  plot_3d(heightmap = dtm_3d,
          zscale = 8,
          theta = 143.5, #323
          phi = 45,
          fov = 45,
          windowsize = c(1000, 500), #c(6400, 3000), #c(1600, 750), #
          background = "transparent",
          solidcolor = "black",
          solidlinecolor = "black",
          zoom = 1, #
          solid = FALSE)
render_camera(zoom = 0.65)
render_snapshot()
render_highquality(filename = "3d_dtm.png", 
                   clear = FALSE, 
                   samples = 128,
                   light = TRUE,
                   lightdirection = 45,
                   lightintensity = 1000,
                   #lightaltitude = 45,
                   width = 1000,#round(6400), 
                   height = 500,#round(3000),
                   parallel = TRUE)

dtm_3d %>%
  sphere_shade(texture = "desert") %>%
  render_image_overlay(rgb_3d)  
#add_overlay(rgb_3d, alphalayer = 0.8)
  #add_water(detect_water(elevation_matrix), color = "desert") %>%
  #add_shadow(ray_shade(dtm_3d), 0.1) %>%
  #add_shadow(ambient_shade(dtm_3d), 0.1) %>%
  

dtm_3d %>%
  sphere_shade(texture ="desert") %>%
  #add_water(detect_water(elevation_matrix), color = "desert") %>%
  add_shadow(ray_shade(dtm_3d), max_darken = 0.5) %>%
  add_shadow(ambient_shade(dtm_3d), max_darken = 0.5) %>%
  add_overlay(rgb_3d, alphalayer = 0.8) %>%
  plot_3d(heightmap = dtm_3d,
          zscale = 8,
          theta = 143.5, #323
          phi = 45,
          fov = 45,
          windowsize = c(1000, 500), #c(6400, 3000), #c(1600, 750), #
          background = "transparent",
          solidcolor = "black",
          solidlinecolor = "black",
          zoom = 1, #
          solid = FALSE)
Sys.sleep(0.2)
render_snapshot()


#render_snapshot()

