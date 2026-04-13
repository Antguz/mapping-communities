#'------------------------------------------------------------------------------
#' @title 3D topographic image for mapping biodiversity
#'------------------------------------------------------------------------------

#' @description Render 3d relief images for our biodiversity products

#' @caution Run in Linux base machine for better performance due to parallel needs

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(terra)
library(rayshader)
library(rgl)
library(viridis)
library(png)
library(ggplot2)
library(viridis)
library(ggpubr)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/guzman/Work/MMDB"
outpath <- paste0(root_path, "/figures/EMIT/For_images")

#'------------------------------------------------------------------------------
#' @Source-code
#-------------------------------------------------------------------------------

source("figures/RGB.R")

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

# DTM plot ---------------------------------------------------------------------

# Area of interest
aoi <- vect(paste0(root_path, "/figures/EMIT/For_images/limits.gpkg"))

# Read DTM data
dtm <- rast(paste0(root_path, "/figures/EMIT/For_images/DTM_new.tif"))

# Resample DTM to match hyperspectral image
dtm_resampled <- resample(dtm, 
                          rast(paste0(root_path, "/figures/EMIT/For_images/RGB_mosaic.tif")), 
                          method = "bilinear")

# Crop 
dtm <- mask(crop(dtm_resampled, aoi), aoi)
dtm_3d <- raster_to_matrix(dtm)

# Check for potential errors
anyNA(dtm_3d) # TRUE
any(is.nan(dtm_3d)) # FALSE
any(is.infinite(dtm_3d)) # FALSE

################################################################################
##### Arguments for preparation

zscale_value <- 12
zoom_value <- 0.27 #0.25
theta_value <- 326.5
phi_value <- 45
fov_value <- 45
windowsize_value <- c(640, 480) 
background_color <- "black" #"white"
solidcolor_color <- "grey40" #"black"
solidlinecolor_color <- "grey60" #"grey40"
shadowcolor_color <- "grey70" #"grey30"
title_color_color <- "grey99" #"grey99"
title_size_value <-  400
title_offset_value <- c(448, 88)
width_value <- 12000
height_value <- 4000

# Initial ray-tracing ----------------------------------------------------------

# rayshade <- ray_shade(dtm_3d, 
#                       zscale = zscale_value,
#                       multicore = TRUE)
# saveRDS(rayshade, file = paste0(root_path, "/figures/EMIT/For_images/rayshade.rds"))

rayshade <- readRDS(paste0(root_path, "/figures/EMIT/For_images/rayshade.rds"))

# ambientshade <- ambient_shade(dtm_3d, 
#                               zscale = zscale_value,
#                               multicore = TRUE)
# saveRDS(ambientshade, file = paste0(root_path, "/figures/EMIT/For_images/ambientshade.rds"))

ambientshade <- readRDS(paste0(root_path, "/figures/EMIT/For_images/ambientshade.rds"))

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

################################################################################
##### DTM plot

rgl::clear3d()
rgl::rgl.viewpoint(theta = 0, phi = 0, fov = 60)

dtm_3d %>%
  height_shade(texture = viridis(256, alpha = 1, begin = 0, end = 1, direction = 1, option = "turbo")) %>%
  add_shadow(rayshade, 0.3) %>%
  add_shadow(ambientshade, 0.9) %>%
  plot_3d(heightmap = dtm_3d,
          zscale = zscale_value,
          theta = theta_value,
          phi = phi_value,
          fov = fov_value,
          windowsize = windowsize_value,
          background = background_color,
          solidcolor = solidcolor_color,
          solidlinecolor = solidlinecolor_color,
          shadowcolor = shadowcolor_color,
          zoom = 1,
          solid = TRUE)
render_camera(zoom = zoom_value)
render_snapshot(filename = paste0(outpath, "/DTM_3d.png"), 
                clear = TRUE, 
                width = width_value, 
                height = height_value,
                software_render = TRUE,
                parallel = TRUE,
                title_text = "Digital Elevation Model",
                title_font = "Arial",
                title_color = title_color_color,
                title_size = title_size_value ,
                title_offset = title_offset_value)

################################################################################
##### RGB plot

refl <- rast(paste0(root_path, "/figures/EMIT/For_images/RGB_mosaic.tif"))
rgb <- refl
# names(refl) <- readRDS("data/EMIT_bands.rds")
# rgb <- c(refl$`648.730285644531`, refl$`551.86669921875`, refl$`447.742797851562`)
names(rgb) <- c("RED", "GREEN", "BLUE")
rgb <- mask(crop(rgb, aoi), aoi)
rgb <- stretch_percentile(rgb, values = c(0.0, 3000))
plotRGB(rgb, scale = 255)

# Export and read
rgb <- as.array(rgb) / 255
writePNG(rgb, paste0(outpath, "/RGB_file.png"))
rgb <- readPNG(paste0(outpath, "/RGB_file.png"))

dtm_3d %>%
  sphere_shade(texture ="desert") %>%
  add_overlay(rgb, alphalayer = .99) %>%
  add_shadow(rayshade, 0.3) %>%
  add_shadow(ambientshade, 0.9) %>%
  plot_3d(dtm_3d,
          zscale = zscale_value,
          theta = theta_value,
          phi = phi_value,
          fov = fov_value,
          windowsize = windowsize_value,
          background = background_color,
          solidcolor = solidcolor_color,
          solidlinecolor = solidlinecolor_color,
          shadowcolor = shadowcolor_color,
          zoom = 1,
          solid = TRUE)
render_camera(zoom = zoom_value)
render_snapshot(filename = paste0(outpath, "/RGB_3d.png"), 
                clear = FALSE, 
                width = width_value, 
                height = height_value,
                software_render = TRUE,
                parallel = TRUE,
                title_text = "EMIT (RGB)",
                title_font = "Arial",
                title_color = title_color_color,
                title_size = title_size_value,
                title_offset = title_offset_value)

################################################################################
##### MDS plots

# Ordinations
MDS_3d_plot <- function(MDS,
                        title_name,
                        dtm_3d = get("dtm_3d", envir = globalenv()), 
                        rayshade = get("rayshade", envir = globalenv()), 
                        ambientshade = get("ambientshade", envir = globalenv()),
                        DTM = get("dtm", envir = globalenv()),
                        type = get("type", envir = globalenv()),
                        outpath = get("outpath", envir = globalenv())) {
  
  # Clamp
  MDS_clamp <- clamp(MDS, lower = -0.5, upper =  0.5, values = TRUE)
  # MDS_clamp <- (clamp(MDS, -0.5, 0.5) + 0.5) * 255
  
  # Fill NA
  dtm_mask <- !is.na(DTM)
  na_overlap <- is.na(MDS_clamp[[1]]) & dtm_mask
  MDS_clamp$Axis1[na_overlap] <- 51
  MDS_clamp$Axis2[na_overlap] <- 51
  MDS_clamp$Axis3[na_overlap] <- 51

  # Export
  # MDS_clamp <- as.array(MDS_clamp) / 255
  arr <- as.array(MDS_clamp)
  arr <- (arr + 0.5) / 1.0
  arr[is.na(arr)] <- 0
  
  writePNG(arr, paste0(outpath, "/", type, "_file.png"))
  MDS_raster <- readPNG(paste0(outpath, "/", type, "_file.png"))
  
  # Create RGB
  
  dtm_3d %>%
    sphere_shade(texture ="desert") %>%
    add_overlay(MDS_raster, alphalayer = .99) %>%
    add_shadow(rayshade, 0.3) %>%
    add_shadow(ambientshade, 0.9) %>%
    plot_3d(dtm_3d,
            zscale = zscale_value,
            theta = theta_value,
            phi = phi_value,
            fov = fov_value,
            windowsize = windowsize_value,
            background = background_color,
            solidcolor = solidcolor_color,
            solidlinecolor = solidlinecolor_color,
            shadowcolor = shadowcolor_color,
            zoom = 1,
            solid = TRUE)
  render_camera(zoom = zoom_value)
  render_snapshot(filename = paste0(outpath, "/", type, "_3d.png"), 
                  clear = FALSE, 
                  width = width_value, 
                  height = height_value,
                  software_render = TRUE,
                  parallel = TRUE,
                  title_text = title_name,
                  title_font = "Arial",
                  title_color = title_color_color,
                  title_size = title_size_value,
                  title_offset = title_offset_value)
  
  gc()
  
}

type <- "taxonomic"
level <- "L3A"
title_name <- "Taxonomic composition"
MDS <- rast(paste0(root_path, "/figures/EMIT/For_images/L3A_mosaic.tif"))
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[1:3]]/10000
MDS_3d_plot(MDS, title_name)

type <- "phylogenetic"
level <- "L3B"
title_name <- "Phylogenetic composition"
MDS <- rast(paste0(root_path, "/figures/EMIT/For_images/L3B_mosaic.tif"))
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[1:3]]/10000
MDS_3d_plot(MDS, title_name)

type <- "functional"
level <- "L3C"
title_name <- "Functional composition"
MDS <- rast(paste0(root_path, "/figures/EMIT/For_images/L3C_mosaic.tif"))
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[1:3]]/10000
MDS_3d_plot(MDS, title_name)

################################################################################
# Ordination error

Error_3d_plot <- function(MDS,
                          title_name,
                          dtm_3d = get("dtm_3d", envir = globalenv()), 
                          rayshade = get("rayshade", envir = globalenv()), 
                          ambientshade = get("ambientshade", envir = globalenv()),
                          DTM = get("dtm", envir = globalenv()),
                          type = get("type", envir = globalenv()),
                          outpath = get("outpath", envir = globalenv())) {
  
  # Clamp
  MDS_clamp <- stretch_percentile(MDS$`amp-total`, values = c(0, 0.0015), DTM = dtm)
  MDS_clamp_3d <- raster_to_matrix(MDS_clamp)
  
  # Plot
  MDS_clamp_3d %>%
    height_shade(texture = c("grey20", viridis(255, 
                                               alpha = 1, 
                                               begin = 0, 
                                               end = 1, 
                                               direction = 1, 
                                               option = "B"))) %>%
    add_shadow(rayshade, 0.5) %>%
    add_shadow(ambientshade, 0.5) %>%
    plot_3d(heightmap = dtm_3d,
            zscale = zscale_value,
            theta = theta_value,
            phi = phi_value,
            fov = fov_value,
            windowsize = windowsize_value,
            background = background_color,
            solidcolor = solidcolor_color,
            solidlinecolor = solidlinecolor_color,
            shadowcolor = shadowcolor_color,
            zoom = 1,
            solid = TRUE)
  render_camera(zoom = zoom_value)
  render_snapshot(filename = paste0(outpath, "/Uncertainty_", type, "_3d.png"), 
                  clear = FALSE, 
                  width = width_value, 
                  height = height_value,
                  software_render = TRUE,
                  parallel = TRUE,
                  title_text = title_name,
                  title_font = "Arial",
                  title_color = title_color_color,
                  title_size = title_size_value ,
                  title_offset = title_offset_value)
  
  gc()
  
}


type <- "taxonomic"
level <- "L3A"
title_name <- "Uncertainty for taxonomic composition"
MDS <- rast(paste0(root_path, "/figures/EMIT/For_images/L3A_mosaic.tif"))
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[4:7]]/1000000
Error_3d_plot(MDS, title_name)

type <- "phylogenetic"
level <- "L3B"
title_name <- "Uncertainty for phylogenetic composition"
MDS <- rast(paste0(root_path, "/figures/EMIT/For_images/L3B_mosaic.tif"))
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[4:7]]/1000000
Error_3d_plot(MDS, title_name)

type <- "functional"
level <- "L3C"
title_name <- "Uncertainty for functional composition"
MDS <- rast(paste0(root_path, "/figures/EMIT/For_images/L3C_mosaic.tif"))
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[4:7]]/1000000
Error_3d_plot(MDS, title_name)

################################################################################
##### Points plots

# Ordinations
points_3d_plot <- function(MDS,
                           title_name,
                           dtm_3d = get("dtm_3d", envir = globalenv()), 
                           rayshade = get("rayshade", envir = globalenv()), 
                           ambientshade = get("ambientshade", envir = globalenv()),
                           DTM = get("dtm", envir = globalenv()),
                           type = get("type", envir = globalenv()),
                           outpath = get("outpath", envir = globalenv())) {
  
  # Clamp
  MDS$RGB <- RGB(MDS$Axis1, MDS$Axis2, MDS$Axis3)
  
  # 
  DTM <- DTM$file246665d20473b
  elevation <- extract(DTM, MDS) + 50
  MDS$elevation <- elevation[, 2]
  coords <- crds(MDS, df = TRUE)

  # Create RGB
  
  dtm_3d %>%
    sphere_shade(texture = create_texture("grey75", "grey80", "grey85", "grey90", "grey95")) %>%
    add_shadow(rayshade, 0.3) %>%
    add_shadow(ambientshade, 0.9) %>%
    plot_3d(heightmap = dtm_3d,
            zscale = zscale_value,
            theta = theta_value,
            phi = phi_value,
            fov = fov_value,
            windowsize = windowsize_value,
            background = background_color,
            solidcolor = solidcolor_color,
            solidlinecolor = solidlinecolor_color,
            shadowcolor = shadowcolor_color,
            zoom = 1,
            solid = TRUE)
  render_camera(zoom = zoom_value)
  render_points(extent = DTM,
                lat = coords$y,
                long = coords$x,
                altitude = MDS$elevation,
                zscale = zscale_value,
                size = 5,
                color = MDS$RGB,
                clear_previous = TRUE)
  render_snapshot(filename = paste0(outpath, "/", type, "_points_3d.png"), 
                  clear = FALSE, 
                  width = width_value, 
                  height = height_value,
                  software_render = TRUE,
                  parallel = TRUE,
                  title_text = title_name,
                  title_font = "Arial",
                  title_color = title_color_color,
                  title_size = title_size_value,
                  title_offset = title_offset_value,
                  point_radius = 35)
  
    gc()
  
}

type <- "taxonomic"
level <- "L3A"
title_name <- "Tβtotal MDS"
MDS <- vect(paste0(root_path, "/figures/EMIT/For_images/L3A_points.gpkg"))
MDS <- mask(crop(MDS, aoi), aoi)
points_3d_plot(MDS, title_name)

type <- "phylogenetic"
level <- "L3B"
title_name <- "Pβtotal MDS"
MDS <- vect(paste0(root_path, "/figures/EMIT/For_images/L3B_points.gpkg"))
MDS <- mask(crop(MDS, aoi), aoi)
points_3d_plot(MDS, title_name)

type <- "functional"
level <- "L3C"
title_name <- "Fβtotal MDS"
MDS <- vect(paste0(root_path, "/figures/EMIT/For_images/L3C_points.gpkg"))
MDS <- mask(crop(MDS, aoi), aoi)
points_3d_plot(MDS, title_name)