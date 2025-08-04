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
          zscale = 6,
          theta = 143.5, #323
          phi = 45,
          fov = 45,
          windowsize = c(4000, 2000), #c(6400, 3000), #c(1600, 750), #
          background = "transparent",
          solidcolor = "black",
          solidlinecolor = "black",
          zoom = 1, #
          solid = TRUE)
render_camera(zoom = 0.65)
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
  sphere_shade(texture="desert") %>%
  plot_3d(dtm_3d, 
          zscale = 35, 
          theta = 143.5, #323
          phi = 45, 
          fov = 45,
          windowsize = c(400, 200),#c(1600, 750), #6400 3000
          background = "transparent",
          solidcolor = "black",
          solidlinecolor = "black",
          #asp = 1/cospi(36.84534/180),
          zoom = 1, #
          solid = TRUE)
render_camera(zoom = 0.65)
render_snapshot(filename = "3d_dtm.png", 
                clear = FALSE, 
                background = "transparent",
                reset_scene_cache = TRUE,
                software_render = TRUE)
render_snapshot()

















aoi_sf <- st_read(paste0(root_path, "/figures/Maps/limits.gpkg"))
aoi_dtm <- get_elev_raster(
  locations = aoi_sf, 
  z = 7, 
  clip = "locations")
aoi_mat <- raster_to_matrix(aoi_dtm)

aoi_mat %>%
  height_shade(texture = viridis(256, alpha = 1, begin = 0, end = 1, direction = 1, option = "turbo")) %>%
  add_shadow(ray_shade(aoi_mat, zscale = 65), 0.3) %>%
  plot_3d(heightmap = aoi_mat,
          zscale = 65,
          theta = 143.5, #323
          phi = 45,
          fov = 45,
          windowsize = c(1600, 750), #c(6400, 3000), #c(1600, 750), #
          background = "transparent",
          solidcolor = "black",
          solidlinecolor = "black",
          zoom = 1, #
          solid = TRUE)
render_camera(zoom = 0.65)
render_snapshot(filename = "3d_dtm.png", 
                clear = FALSE, 
                background = "transparent",
                reset_scene_cache = TRUE,
                software_render = TRUE)
render_snapshot(parallel = TRUE)

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



render_highquality(
  filename = "3d_dtm.png", 
  clear = FALSE, 
  background = "transparent",
  parallel = TRUE, 
  samples = 300,
  light = FALSE, 
  interactive = FALSE,
  environment_light = "env/phalzer_forest_01_4k.hdr",
  intensity_env = 1.5,
  rotate_env = 180,
  width = round(6000), 
  height = round(6000)
)
