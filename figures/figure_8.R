#'------------------------------------------------------------------------------
#' @title Three-dimensional figure for MDS (Fig 8)
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(bit64)
library(plot3D)
library(data.table)
library(schemr)
library(stringr)
library(rcartocolor)
library(viridis)
library(colorspace)
library(Rcpp)
library(RcppArmadillo)
library(terra)
library(ggplot2)
library(ggdist)

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

root_path <- "/media/guzman/Work/MMDB/figures/EMIT/Land-change"

# Rasters
phylo2024 <- rast(paste0(root_path, "/2024_L3B.tif"))
phylo2025 <- rast(paste0(root_path, "/2025_L3B.tif"))
euclidean <- rast(paste0(root_path, "/euclidean_distance.tif"))

# Vectors
border <- vect(paste0(root_path, "/border.gpkg"))
aoi <- vect(paste0(root_path, "/aoi.gpkg"))       # Extent of the figure
covers <- vect(paste0(root_path, "/covers.gpkg")) # Field 'Cover': 'change' / 'stable'

# Mask to border and crop to the area of interest
phylo2024 <- crop(phylo2024, border, mask = TRUE)
phylo2025 <- crop(phylo2025, border, mask = TRUE)
phylo2024 <- crop(phylo2024, aoi, mask = TRUE)
phylo2025 <- crop(phylo2025, aoi, mask = TRUE)
euclidean <- crop(euclidean, border, mask = TRUE)
euclidean <- crop(euclidean, aoi, mask = TRUE)

phylodifference <- c(phylo2024, phylo2025)
phylodifference <- phylodifference[[c(1,2,3,8,9,10)]]
names(phylodifference) <- c("Axis1_2024", "Axis2_2024", "Axis3_2024", "Axis1_2025", "Axis2_2025", "Axis3_2025")

# Cover class per pixel (change / stable), same grid as the rasters
cover_class <- rasterize(covers, phylo2024, field = "Cover")
names(cover_class) <- "cover"

phylo2024 <- c(phylo2024, cover_class)
phylo2025 <- c(phylo2025, cover_class)
phylodifference <- c(phylodifference, cover_class)

# As frames
phylo2024 <- as.data.table(as.data.frame(phylo2024, cells = FALSE, na.rm = NA, xy = TRUE))
phylo2025 <- as.data.table(as.data.frame(phylo2025, cells = FALSE, na.rm = NA, xy = TRUE))
phylodifference <- as.data.table(as.data.frame(phylodifference, cells = FALSE, na.rm = NA, xy = TRUE))

# Keep only pixels with ordination values inside the change / stable covers
phylo2024 <- phylo2024[!is.na(Axis1) & !is.na(cover)]
phylo2025 <- phylo2025[!is.na(Axis1) & !is.na(cover)]
phylodifference <- phylodifference[!is.na(Axis1_2024) & !is.na(Axis1_2025) & !is.na(cover)]

phylo2024[, cover := as.character(cover)]
phylo2025[, cover := as.character(cover)]
phylodifference[, cover := as.character(cover)]

# Class colors
class_colors <- c(change = "#c51b7d", stable = "#4d9221")

#'------------------------------------------------------------------------------
#' @Plot-functions
#-------------------------------------------------------------------------------

# RGB color
RGB <- function(X, Y, Z, min_range = -0.5, max_range = 0.5) {
  
  if(is.null(min_range) == TRUE & is.null(max_range) == TRUE) {
    
    x_new <-  ((X - min(X)) / ((max(X)) - (min(X)))) * (255 - 0) + 0
    y_new <-  ((Y - min(Y)) / ((max(Y)) - (min(Y)))) * (255 - 0) + 0
    z_new <-  ((Z - min(Z)) / ((max(Z)) - (min(Z)))) * (255 - 0) + 0
    
  } else {
    
    x_new <-  ((X - min_range) / (max_range - min_range)) * (255 - 0) + 0
    y_new <-  ((Y - min_range) / (max_range - min_range)) * (255 - 0) + 0
    z_new <-  ((Z - min_range) / (max_range - min_range)) * (255 - 0) + 0
    
  }
  
  # Convert colour channels to hex
  hex_value <- paste0(
    "#",
    # Convert red
    str_pad(
      string = as.hexmode(round(x_new, 0)),
      width = 2, pad = "0", side = "left"
    ),
    # Convert green
    str_pad(
      string = as.hexmode(round(y_new, 0)),
      width = 2, pad = "0", side = "left"
    ),
    # Convert blue
    str_pad(
      string = as.hexmode(round(z_new, 0)),
      width = 2, pad = "0", side = "left"
    )
  )
  
  return(hex_value)
  
}

# Ordination figure
ordination_3d <- function(ordination, 
                          title = " ", 
                          phi = 20, 
                          theta = 45,
                          class_colors,
                          axis_range,
                          alpha = 1.0) {
  
  # Scale
  Axis1 <- clamp(ordination$Axis1 / 10000, -0.5, 0.5)
  Axis2 <- clamp(ordination$Axis2 / 10000, -0.5, 0.5)
  Axis3 <- clamp(ordination$Axis3 / 10000, -0.5, 0.5)
  
  # Plot features
  point_shape <- 20
  point_size_shade <- 0.1
  point_size <- 0.1
  xlab <- paste0("Axis 1")
  ylab <- paste0("Axis 2")
  zlab <- paste0("Axis 3")
  shade_color <- "grey75"
  alpha_color <- alpha
  side <- 1
  dist <- -0.05
  null_colkey <- list(plot = FALSE, 
                      side = side, 
                      side.clab = 1,
                      dist = dist)
  axis_leter <- 1
  axis_label <- 1.2
  clab <- NULL

  # Pixels per class
  stable <- ordination$cover == "stable"
  change <- ordination$cover == "change"

  scatter3D(x = rep(axis_range[1], length(Axis1)), y = Axis2, z = Axis3, 
            phi = phi,
            theta = theta,
            bty = "u", 
            pch = 16, 
            cex = point_size_shade, 
            ticktype = "detailed",
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            xlim = axis_range,
            ylim = axis_range,
            zlim = axis_range,
            nticks = 3,
            col = shade_color,
            colvar = NULL,
            colkey = null_colkey,
            clab = clab,
            alpha.col = 0.2,
            #col.panel = "black",
            col.grid = "grey75",
            col.axis = "white",
            main = title,
            cex.axis = axis_leter,
            cex.lab = axis_label,
            scale = FALSE,
            d = 1,
            zoom = 1)
  
  scatter3D(x = Axis1, y = rep(axis_range[2], length(Axis2)), z = Axis3, 
            phi = phi, 
            theta = theta,
            bty = "u",
            pch = 16, 
            cex = point_size_shade, 
            ticktype = "detailed",
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            xlim = axis_range,
            ylim = axis_range,
            zlim = axis_range,
            nticks = 3,
            col = shade_color,
            colvar = NULL,
            colkey = null_colkey,
            clab = clab,
            alpha.col = 0.2,
            add = TRUE)
  
  scatter3D(x = Axis1, y = Axis2, z = rep(axis_range[1], length(Axis3)), 
            phi = phi, 
            theta = theta,
            bty = "u", 
            pch = 16, 
            cex = point_size_shade, 
            ticktype = "detailed",
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            xlim = axis_range,
            ylim = axis_range,
            zlim = axis_range,
            nticks = 3,
            col = shade_color,
            colvar = NULL,
            colkey = null_colkey,
            clab = clab,
            alpha.col = 0.2,
            add = TRUE)
  

  # Stable pixels
  scatter3D(x = Axis1[stable], y = Axis2[stable], z = Axis3[stable],
            phi = phi,
            theta = theta,
            bty = "u",
            pch = point_shape,
            cex = 0.5,
            ticktype = "detailed",
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            xlim = axis_range,
            ylim = axis_range,
            zlim = axis_range,
            nticks = 3,
            col = class_colors["stable"],
            colvar = NULL,
            colkey = null_colkey,
            clab = clab,
            alpha.col = alpha_color,
            add = TRUE)

  # Change pixels
  scatter3D(x = Axis1[change], y = Axis2[change], z = Axis3[change],
            phi = phi,
            theta = theta,
            bty = "u",
            pch = point_shape,
            cex = 0.5,
            ticktype = "detailed",
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            xlim = axis_range,
            ylim = axis_range,
            zlim = axis_range,
            nticks = 3,
            col = class_colors["change"],
            colvar = NULL,
            colkey = null_colkey,
            clab = clab,
            alpha.col = alpha_color,
            add = TRUE)

}


#-------------------------------------------------------------------------------

ordination_arrow_3d <- function(ordination_difference,
                                title = "Change 2024–2025",
                                phi = 20,
                                theta = 45,
                                class_colors,
                                axis_range,
                                arrow_width = 0.3) {
  
  # Plot features
  point_size_shade <- 0.1
  xlab <- paste0("Axis 1")
  ylab <- paste0("Axis 2")
  zlab <- paste0("Axis 3")
  shade_color <- "grey75"
  side <- 1
  dist <- -0.05
  arrlength <- 0.10
  arrangle  <- 30
  null_colkey <- list(plot = FALSE, 
                      side = side, 
                      side.clab = 1,
                      dist = dist)
  axis_leter <- 1
  axis_label <- 1.2
  clab <- NULL
  
  # Pixels per class: stable first, change drawn on top
  dt <- ordination_difference[cover %in% c("change", "stable"),]
  dt <- dt[order(cover != "stable")]
  
  Axis1_2024 <- clamp(dt$Axis1_2024 / 10000, -0.5, 0.5)
  Axis2_2024 <- clamp(dt$Axis2_2024 / 10000, -0.5, 0.5)
  Axis3_2024 <- clamp(dt$Axis3_2024 / 10000, -0.5, 0.5)
  
  Axis1_2025 <- clamp(dt$Axis1_2025 / 10000, -0.5, 0.5)
  Axis2_2025 <- clamp(dt$Axis2_2025 / 10000, -0.5, 0.5)
  Axis3_2025 <- clamp(dt$Axis3_2025 / 10000, -0.5, 0.5)
  
  col <- class_colors[dt$cover]
  
  scatter3D(x = 0,
            y = 0,
            z = 0,
            phi = phi,
            theta = theta,
            bty = "u", 
            pch = 16, 
            cex = point_size_shade, 
            ticktype = "detailed",
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            xlim = axis_range,
            ylim = axis_range,
            zlim = axis_range,
            nticks = 3,
            col = "white",
            colvar = NULL,
            colkey = null_colkey,
            clab = clab,
            alpha.col = 0.2,
            col.grid = "grey75",
            col.axis = "white",
            main = title,
            cex.axis = axis_leter,
            cex.lab = axis_label,
            scale = FALSE,
            d = 1,
            zoom = 1)
  
  # Projections on the panels
  arrows3D(x0 = rep(axis_range[1], length(Axis1_2024)), y0 = Axis2_2024, z0 = Axis3_2024,
           x1 = rep(axis_range[1], length(Axis1_2025)), y1 = Axis2_2025, z1 = Axis3_2025,
           phi = phi,
           theta = theta,
           col = shade_color,
           add = TRUE,
           length = arrlength,
           angle  = arrangle,
           lwd = 0.1)
  
  arrows3D(x0 = Axis1_2024, y0 = rep(axis_range[2], length(Axis2_2024)), z0 = Axis3_2024,
           x1 = Axis1_2025, y1 = rep(axis_range[2], length(Axis2_2025)), z1 = Axis3_2025,
           phi = phi,
           theta = theta,
           col = shade_color,
           add = TRUE,
           length = arrlength,
           angle  = arrangle,
           lwd = 0.1)
  
  arrows3D(x0 = Axis1_2024, y0 = Axis2_2024, z0 = rep(axis_range[1], length(Axis3_2024)),
           x1 = Axis1_2025, y1 = Axis2_2025, z1 = rep(axis_range[1], length(Axis3_2025)),
           phi = phi,
           theta = theta,
           col = shade_color,
           add = TRUE,
           length = arrlength,
           angle  = arrangle,
           lwd = 0.1)
  
  # Arrows per pixel colored by class
  arrows3D(x0 = Axis1_2024, y0 = Axis2_2024, z0 = Axis3_2024,
           x1 = Axis1_2025, y1 = Axis2_2025, z1 = Axis3_2025,
           phi = phi,
           theta = theta,
           col = col,
           alpha = 0.6,
           add = TRUE,
           length = arrlength,
           angle  = arrangle,
           lwd = arrow_width)
  
}

#-------------------------------------------------------------------------------
# Figure

# Axis range
observed_range <- c(-0.6, 0.6)

png(paste0(root_path, "/", "Figure_8.png"),
    width = 90, 
    height = 210, 
    units = "mm", 
    pointsize = 12,
    bg = "white",  
    res = 600)

par(mfrow=c(3,1), mar = c(0,1.1,1,0), oma = c(0,0,0,0))
ordination_3d(phylo2024, title = "2024", class_colors = class_colors, axis_range = observed_range)
ordination_3d(phylo2025, title = "2025", class_colors = class_colors, axis_range = observed_range)
ordination_arrow_3d(phylodifference, title = "Euclidean distance between 2025 and 2024", class_colors = class_colors, axis_range = observed_range)

dev.off()

#-------------------------------------------------------------------------------
# Independent panels

# Panel 2024
png(paste0(root_path, "/", "Figure_8a_2024.png"),
    width = 90,
    height = 70,
    units = "mm",
    pointsize = 12,
    bg = "white",
    res = 600)

par(mar = c(0,1.1,1,0), oma = c(0,0,0,0), cex = 0.66) # cex as in mfrow = c(3,1)
ordination_3d(phylo2024, title = "2024", class_colors = class_colors, axis_range = observed_range)

dev.off()

# Panel 2025
png(paste0(root_path, "/", "Figure_8b_2025.png"),
    width = 90,
    height = 70,
    units = "mm",
    pointsize = 12,
    bg = "white",
    res = 600)

par(mar = c(0,1.1,1,0), oma = c(0,0,0,0), cex = 0.66) # cex as in mfrow = c(3,1)
ordination_3d(phylo2025, title = "2025", class_colors = class_colors, axis_range = observed_range)

dev.off()

# Panel change 2024-2025
png(paste0(root_path, "/", "Figure_8c_change.png"),
    width = 90,
    height = 70,
    units = "mm",
    pointsize = 12,
    bg = "white",
    res = 600)

par(mar = c(0,1.1,1,0), oma = c(0,0,0,0), cex = 0.66) # cex as in mfrow = c(3,1)
ordination_arrow_3d(phylodifference, title = "Euclidean distance between 2025 and 2024", class_colors = class_colors, axis_range = observed_range)

dev.off()

#-------------------------------------------------------------------------------
# Raincloud plot

# Euclidean distance per pixel from 'euclidean_distance.tif'
euclidean_class <- rasterize(covers, euclidean, field = "Cover")
names(euclidean_class) <- "cover"
names(euclidean) <- "distance"

distance <- as.data.table(as.data.frame(c(euclidean, euclidean_class), cells = FALSE, na.rm = NA, xy = FALSE))
distance <- distance[!is.na(distance) & !is.na(cover)]
distance[, distance := distance / 10000] # Same scale as the ordination panels
distance[, cover := factor(as.character(cover), levels = c("stable", "change"), labels = c("Stable", "Disturbed"))]
distance[, position := as.numeric(cover)]

# Summary per class
print(distance[, .(n = .N, 
                   median = median(distance), 
                   q25 = quantile(distance, 0.25), 
                   q75 = quantile(distance, 0.75)), 
               by = cover])

# Cloud: half density per class
cloud <- distance[, {
  d <- density(distance, from = 0)
  .(x = d$x, y = d$y / max(d$y) * 0.4)
}, by = .(cover, position)]



# Colors
raincloud_colors <- c(Stable = class_colors[["stable"]], Disturbed = class_colors[["change"]])

raincloud <- ggplot(data = distance, aes(x = cover, y = distance, fill = cover)) +
  
  # Points
  geom_jitter(#colour = "grey",
              aes(colour = cover),
              width = 0.03,
              height = 0,
              shape = 1,
              size = 0.3,
              alpha = 0.25) +
  # Box
  geom_boxplot(width = 0.06,
               outlier.shape = NA,
               alpha = 0.85,
               linewidth = 0.3,
               colour = "grey25") +
  # Rain
  #stat_dots(side = "left", binwidth = 0.0015, justification = 1.25, width = 0.4) +
  # Cloud
  stat_halfeye(adjust = .5, width = .4, justification = -.2, .width = 0, point_colour = NA, alpha = 0.5) +
  scale_fill_manual(name = "Stand trajectory", values = raincloud_colors) +
  scale_colour_manual(name = "Stand trajectory", values = raincloud_colors) +
  coord_cartesian(xlim = c(0.5, 2.5), ylim = c(0, 0.5), expand = FALSE) +
  scale_y_continuous(n.breaks = 3) +
  theme_bw(base_size = 10, base_family = "Arial") +
  theme(legend.position = "inside",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.background = element_rect(fill= "transparent"),
        legend.position.inside = c(0.90, 0.98),
        legend.justification.inside = c(1, 1),
        legend.background = element_rect(fill = "white", colour = NA)) +
  ylab("Euclidean distance (-)") +
  xlab("Stand trajectory")

ggsave(paste0(root_path, "/", "Figure_8d_rainclouda.png"),
       plot = raincloud,
       width = 60,
       height = 62,
       units = "mm",
       dpi = 600,
       bg = "white")
