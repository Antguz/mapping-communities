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

root_path <- "/media/antonio/Work/MMDB"
outpath <- paste0(root_path, "/figures/Maps/3d_figures")
  
#'------------------------------------------------------------------------------
#' @Source-code
#-------------------------------------------------------------------------------

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

################################################################################
##### Arguments for preparation

zscale_value <- 12
zoom_value <- 0.25
theta_value <- 143.5
phi_value <- 45
fov_value <- 45
windowsize_value <- c(640, 480) 
background_color <- "white"
solidcolor_color <- "grey40"
solidlinecolor_color <- "grey40"
shadowcolor_color <- "grey30"
title_color_color <- "grey99"
title_size_value <-  100
title_offset_value <- c(448, 88)
width_value <- 3000
height_value <- 1000

# Initial ray-tracing ----------------------------------------------------------

rayshade <- ray_shade(dtm_3d, 
                      zscale = zscale_value,
                      multicore = TRUE)

ambientshade <- ambient_shade(dtm_3d, 
                              zscale = zscale_value,
                              multicore = TRUE)

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
                clear = FALSE, 
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

refl <- rast(paste0(root_path, "/figures/Maps/Rflec/reflectance.tif"))
names(refl) <- readRDS("data/EMIT_bands.rds")
rgb <- c(refl$`648.730285644531`, refl$`551.86669921875`, refl$`447.742797851562`)
names(rgb) <- c("RED", "GREEN", "BLUE")
rgb <- mask(crop(rgb, aoi), aoi)
rgb <- stretch_percentile(rgb, values = c(0.0, 0.3))
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
MDS1 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175417_2323212_006_MDS-",type,"-composition.tif"))
MDS2 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175429_2323212_007_MDS-",type,"-composition.tif"))
MDS <- merge(MDS1, MDS2)
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[1:3]]/10000
MDS_3d_plot(MDS, title_name)

type <- "phylogenetic"
level <- "L3B"
title_name <- "Phylogenetic composition"
MDS1 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175417_2323212_006_MDS-",type,"-composition.tif"))
MDS2 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175429_2323212_007_MDS-",type,"-composition.tif"))
MDS <- merge(MDS1, MDS2)
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[1:3]]/10000
MDS_3d_plot(MDS, title_name)

type <- "functional"
level <- "L3C"
title_name <- "Functional composition"
MDS1 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175417_2323212_006_MDS-",type,"-composition.tif"))
MDS2 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175429_2323212_007_MDS-",type,"-composition.tif"))
MDS <- merge(MDS1, MDS2)
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[1:3]]/10000
MDS_3d_plot(MDS, title_name)

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
MDS1 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175417_2323212_006_MDS-",type,"-composition.tif"))
MDS2 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175429_2323212_007_MDS-",type,"-composition.tif"))
MDS <- merge(MDS1, MDS2)
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[4:7]]/1000000
Error_3d_plot(MDS, title_name)

type <- "phylogenetic"
level <- "L3B"
title_name <- "Uncertainty for phylogenetic composition"
MDS1 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175417_2323212_006_MDS-",type,"-composition.tif"))
MDS2 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175429_2323212_007_MDS-",type,"-composition.tif"))
MDS <- merge(MDS1, MDS2)
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[4:7]]/1000000
Error_3d_plot(MDS, title_name)

type <- "functional"
level <- "L3C"
title_name <- "Uncertainty for functional composition"
MDS1 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175417_2323212_006_MDS-",type,"-composition.tif"))
MDS2 <- rast(paste0(root_path, "/figures/Maps/Diversity/EMIT-",level,"_001_20230820T175429_2323212_007_MDS-",type,"-composition.tif"))
MDS <- merge(MDS1, MDS2)
MDS <- mask(crop(MDS, aoi), aoi)
MDS <- MDS[[4:7]]/1000000
Error_3d_plot(MDS, title_name)

################################################################################
##### FORTYPE plots

fortyp1 <- rast(paste0(root_path, "/figures/Maps/FORTYP/EMIT-L4A_001_20230820T175417_2323212_006_FORTYP_mean.tif"))
fortyp2 <- rast(paste0(root_path, "/figures/Maps/FORTYP/EMIT-L4A_001_20230820T175429_2323212_007_FORTYP_mean.tif"))
fortyp <- merge(fortyp1, fortyp2)
fortyp <- mask(crop(fortyp, aoi), aoi)

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
      height_shade(texture = c("grey20", viridis(255, 
                                                 alpha = 1, 
                                                 begin = 0, 
                                                 end = 1, 
                                                 direction = 1, 
                                                 option = "A"))) %>%
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
    out_name <- gsub(" / ", "_", layer_names[i])
    out_name <- gsub(" ", "-", out_name)
    render_snapshot(filename = paste0(outpath, "/", out_name, ".png"), 
                    clear = FALSE, 
                    width = width_value, 
                    height = height_value,
                    software_render = TRUE,
                    parallel = TRUE,
                    title_text = layer_names[i],
                    title_font = "Arial",
                    title_color = title_color_color,
                    title_size = title_size_value ,
                    title_offset = title_offset_value)
    
    gc()
    
  }
}

fortyp_3d_plot(fortyp)

################################################################################
##### Plant lineage plots

occurrence1 <- rast(paste0(root_path, "/figures/Maps/Occurrence/EMIT-L4B_001_20230820T175417_2323212_006_plant-occurrence_mean.tif"))
occurrence2 <- rast(paste0(root_path, "/figures/Maps/Occurrence/EMIT-L4B_001_20230820T175429_2323212_007_plant-occurrence_mean.tif"))
occurrence <- merge(occurrence1, occurrence2)
occurrence <- mask(crop(occurrence, aoi), aoi)
occurrence <- c(occurrence$Quercus, occurrence$Acer, occurrence$Picea)

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
      height_shade(texture = c("grey20", viridis(255, 
                                                 alpha = 1, 
                                                 begin = 0, 
                                                 end = 1, 
                                                 direction = 1, 
                                                 option = "G"))) %>%
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
    out_name <- gsub(" / ", "_", layer_names[i])
    out_name <- gsub(" ", "-", out_name)
    render_snapshot(filename = paste0(outpath, "/", out_name, ".png"), 
                    clear = FALSE, 
                    width = width_value, 
                    height = height_value,
                    software_render = TRUE,
                    parallel = TRUE,
                    title_text = layer_names[i],
                    title_font = "Arial",
                    title_color = title_color_color,
                    title_size = title_size_value ,
                    title_offset = title_offset_value)
    
    gc()
    
  }
}

occurrence_3d_plot(occurrence)

################################################################################
##### CWM plot -----------------------------------------------------------------

func1 <- rast(paste0(root_path, "/figures/Maps/CWM/EMIT-L4C_001_20230820T175417_2323212_006_CWM_mean.tif"))
func2 <- rast(paste0(root_path, "/figures/Maps/CWM/EMIT-L4C_001_20230820T175429_2323212_007_CWM_mean.tif"))
func <- merge(func1, func2)
func <- mask(crop(func, aoi), aoi)

# Function
trait_3d_plot <- function(trait,
                          title_name, 
                          dtm_3d = get("dtm_3d", envir = globalenv()),
                          rayshade = get("rayshade", envir = globalenv()),
                          ambientshade = get("ambientshade", envir = globalenv()),
                          DTM = get("dtm", envir = globalenv()),
                          outpath = get("outpath", envir = globalenv())) {
  
  trait %>%
    height_shade(texture = c("grey20", viridis(255, 
                                               alpha = 1, 
                                               begin = 0, 
                                               end = 1, 
                                               direction = 1, 
                                               option = "D"))) %>%
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
  render_snapshot(filename = paste0(outpath, "/", title_name, ".png"), 
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

# Tree height
HTq <- func$HTq
stretch_percentile(HTq)
HTq <- stretch_percentile(HTq, values = c(2600, 2850))
HTq_3d <- raster_to_matrix(HTq)
title_name <- "Tree height"
trait_3d_plot(HTq_3d, title_name)

# Slenderness
SLq <- func$SLq
stretch_percentile(SLq)
SLq <- stretch_percentile(SLq, values = c(8000, 8700))
SLq_3d <- raster_to_matrix(SLq)
title_name <- "Slenderness"
trait_3d_plot(SLq_3d, title_name)

# Shade tolerance index 
shade <- func$shade
stretch_percentile(shade)
shade <- stretch_percentile(shade, values = c(100, 500))
shade_3d <- raster_to_matrix(shade)
title_name <- "Shade tolerance index"
trait_3d_plot(shade_3d, title_name)

# Drought tolerance index
drought <- func$drought
stretch_percentile(drought)
drought <- stretch_percentile(drought, values = c(100, 500))
drought_3d <- raster_to_matrix(drought)
title_name <- "Drought tolerance index"
trait_3d_plot(drought_3d, title_name)

# Waterlogging tolerance index"
waterlogging <- func$waterlogging
stretch_percentile(waterlogging)
waterlogging <- stretch_percentile(waterlogging, values = c(100, 500))
waterlogging_3d <- raster_to_matrix(waterlogging)
title_name <- "Waterlogging tolerance index"
trait_3d_plot(waterlogging_3d, title_name)

# Leaf mass per area
LMA <- func$LMA
stretch_percentile(LMA) # 5707 8761
LMA <- stretch_percentile(LMA, values = c(6000, 9500))
LMA_3d <- raster_to_matrix(LMA)
title_name <- "Leaf mass per area"
trait_3d_plot(LMA_3d, title_name)

# Nitrogen
N <- func$N
stretch_percentile(N)
N <- stretch_percentile(N, values = c(1500, 2000))
N_3d <- raster_to_matrix(N)
title_name <- "Nitrogen"
trait_3d_plot(N_3d, title_name)

# Wood density
WD <- func$WD
stretch_percentile(WD)
WD <- stretch_percentile(WD, values = c(45, 60))
WD_3d <- raster_to_matrix(WD)
title_name <- "Wood density"
trait_3d_plot(WD_3d, title_name)

################################################################################
##### Create legends   ---------------------------------------------------------

# DTM
dtm <- ggplot(diamonds, aes(x = clarity, y = depth)) + 
  geom_point(aes(fill = carat), shape = 21, color = "black", show.legend = TRUE) +  # Use fill shape
  scale_fill_viridis_c("Elevation (m a.s.l.)", 
                       option = "turbo",
                       limits = c(180, 2030),
                       n.breaks = 4,
                       guide = guide_colorbar(
                         title.position = "bottom",
                         title.hjust = 0.5, 
                         barwidth = 15,
                         barheight = 1)) +
  theme_minimal(18) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", 
                                     linewidth = 0.5, 
                                     linetype = "solid",
                                     colour = "transparent")
  )

# Uncertainty
error <- ggplot(diamonds, aes(x = clarity, y = depth)) + 
  geom_point(aes(fill = carat), shape = 21, color = "black", show.legend = TRUE) +  # Use fill shape
  scale_fill_viridis_c("Uncertainty", 
                       option = "B",
                       limits = c(0, 0.0015),
                       n.breaks = 2,
                       guide = guide_colorbar(
                         title.position = "bottom",
                         title.hjust = 0.5, 
                         barwidth = 15,
                         barheight = 1)) +
  theme_minimal(18) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", 
                                     linewidth = 0.5, 
                                     linetype = "solid",
                                     colour = "transparent")
  )

# Forest types
probability_forest <- ggplot(diamonds, aes(x = clarity, y = depth)) + 
  geom_point(aes(fill = carat), shape = 21, color = "black", show.legend = TRUE) +  # Use fill shape
  scale_fill_viridis_c("Probability", 
                       option = "A",
                       limits = c(0, 1),
                       n.breaks = 3,
                       guide = guide_colorbar(
                         title.position = "bottom",
                         title.hjust = 0.5, 
                         barwidth = 15,
                         barheight = 1)) +
  theme_minimal(18) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", 
                                     linewidth = 0.5, 
                                     linetype = "solid",
                                     colour = "transparent")
  )

# Plant occurance 
probability_lineage <- ggplot(diamonds, aes(x = clarity, y = depth)) + 
  geom_point(aes(fill = carat), shape = 21, color = "black", show.legend = TRUE) +  # Use fill shape
  scale_fill_viridis_c("Probability", 
                       option = "G",
                       limits = c(0, 1),
                       n.breaks = 3,
                       guide = guide_colorbar(
                         title.position = "bottom",
                         title.hjust = 0.5, 
                         barwidth = 15,
                         barheight = 1)) +
  theme_minimal(18) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", 
                                     linewidth = 0.5, 
                                     linetype = "solid",
                                     colour = "transparent")
  )

# LMA
lma <- ggplot(diamonds, aes(x = clarity, y = depth)) + 
  geom_point(aes(fill = carat), shape = 21, color = "black", show.legend = TRUE) +  # Use fill shape
  scale_fill_viridis_c("g/m2", 
                       option = "D",
                       limits = c(60, 95),
                       n.breaks = 4,
                       guide = guide_colorbar(
                         title.position = "bottom",
                         title.hjust = 0.5, 
                         barwidth = 15,
                         barheight = 1)) +
  theme_minimal(18) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", 
                                     linewidth = 0.5, 
                                     linetype = "solid",
                                     colour = "transparent")
  )

# N
N <- ggplot(diamonds, aes(x = clarity, y = depth)) + 
  geom_point(aes(fill = carat), shape = 21, color = "black", show.legend = TRUE) +  # Use fill shape
  scale_fill_viridis_c("mg/g", 
                       option = "D",
                       limits = c(15, 20),
                       n.breaks = 2,
                       guide = guide_colorbar(
                         title.position = "bottom",
                         title.hjust = 0.5, 
                         barwidth = 15,
                         barheight = 1)) +
  theme_minimal(18) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", 
                                     linewidth = 0.5, 
                                     linetype = "solid",
                                     colour = "transparent")
  )

# WD
WD <- ggplot(diamonds, aes(x = clarity, y = depth)) + 
  geom_point(aes(fill = carat), shape = 21, color = "black", show.legend = TRUE) +  # Use fill shape
  scale_fill_viridis_c("g/cm2", 
                       option = "D",
                       limits = c(0.45, 0.6),
                       breaks = c(0.45, 0.50, 0.55, 0.6),
                       guide = guide_colorbar(
                         title.position = "bottom",
                         title.hjust = 0.5, 
                         barwidth = 15,
                         barheight = 1)) +
  theme_minimal(18) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", 
                                     linewidth = 0.5, 
                                     linetype = "solid",
                                     colour = "transparent")
  )


# Extract legend
dtm_legend <- as_ggplot(get_legend(dtm))
error_legend <- as_ggplot(get_legend(error))
probability_forest_legend <- as_ggplot(get_legend(probability_forest))
probability_lineage_legend <- as_ggplot(get_legend(probability_lineage))
lma_legend <- as_ggplot(get_legend(lma))
N_legend <- as_ggplot(get_legend(N))
WD_legend <- as_ggplot(get_legend(WD))

# Plot legend
plot <- ggarrange(dtm_legend, error_legend, probability_lineage_legend,
                  lma_legend, N_legend, WD_legend, ncol = 1, nrow = 6, align = "v")

# Export for design
ggsave("legends_grid.png", plot, units = "mm", 
       width = 150, height = 200, dpi = 300, bg = "transparent")

