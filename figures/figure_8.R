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

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

root_path <- "/media/guzman/Work/MMDB/figures/EMIT/Land-change"

# Rasters
phylo2024 <- rast(paste0(root_path, "/2024_L3B.tif"))
phylo2025 <- rast(paste0(root_path, "/2025_L3B.tif"))
phylodifference <- rast(paste0(root_path, "/difference_L3B.tif"))

# Vectors
border <- vect(paste0(root_path, "/border.gpkg"))
areas_change <- vect(paste0(root_path, "/areas_of_change.gpkg"))

# Mask
phylo2024 <- crop(phylo2024, border, mask = TRUE)
phylo2025 <- crop(phylo2025, border, mask = TRUE)
phylodifference <- crop(c(phylo2024, phylo2025), border, mask = TRUE)
phylodifference <- phylodifference[[c(1,2,3,8,9,10)]]
names(phylodifference) <- c("Axis1_2024", "Axis2_2024", "Axis3_2024", "Axis1_2025", "Axis2_2025", "Axis3_2025")

# Mask areas of interest
aoi2024 <- crop(phylo2024, areas_change, mask = TRUE)
aoi2025 <-  crop(phylo2025, areas_change, mask = TRUE)
aoidifference <- crop(phylodifference, areas_change, mask = TRUE)

# As frames
phylo2024 <- as.data.frame(phylo2024, cells = FALSE, na.rm = NA, xy = TRUE)
phylo2025 <- as.data.frame(phylo2025, cells = FALSE, na.rm = NA, xy = TRUE)
phylodifference <- as.data.frame(phylodifference, cells = FALSE, na.rm = NA, xy = TRUE)

aoi2024 <- as.data.frame(aoi2024, cells = FALSE, na.rm = NA, xy = TRUE)
aoi2025 <- as.data.frame(aoi2025, cells = FALSE, na.rm = NA, xy = TRUE)
aoidifference <- as.data.frame(aoidifference, cells = FALSE, na.rm = NA, xy = TRUE)

# Overlaps
aoi2024$stands <- "Yes"
aoi2025$stands <- "Yes"
aoidifference$stands <- "Yes"

phylo2024 <- merge(phylo2024, aoi2024, by = c("x", "y", "Axis1", "Axis2", "Axis3", "amp-Axis1", "amp-Axis2", "amp-Axis3", "amp-total"), all.x = TRUE, all.y = TRUE)
phylo2025 <- merge(phylo2025, aoi2025, by = c("x", "y", "Axis1", "Axis2", "Axis3", "amp-Axis1", "amp-Axis2", "amp-Axis3", "amp-total"), all.x = TRUE, all.y = TRUE)
phylodifference <- merge(phylodifference, aoidifference, by = c("x", "y", "Axis1_2024", "Axis2_2024", "Axis3_2024", "Axis1_2025", "Axis2_2025", "Axis3_2025"), all.x = TRUE, all.y = TRUE)

phylo2024 <- as.data.table(phylo2024)
phylo2025 <- as.data.table(phylo2025)
phylodifference <- as.data.table(phylodifference)

phylo2024[is.na(stands), stands := "No"]
phylo2025[is.na(stands), stands := "No"]
phylodifference[is.na(stands), stands := "No"]

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
                          cRGB = TRUE,
                          col = NULL,
                          colvar = NULL,
                          key = TRUE,
                          clab = " ",
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
  
  if(cRGB == TRUE) {
    col <- RGB(Axis1, 
               Axis2, 
               Axis3,
               min_range = axis_range[1], 
               max_range = axis_range[2])
    
    col[ordination$stands == "No"] <- "grey40"
    
    colvar <- NULL
    colkey <- list(plot = FALSE, 
                   side = side, 
                   side.clab = 1,
                   dist = dist)
    clab <- NULL
    
  } else if(cRGB == FALSE) {
    if(key == TRUE) {
      colkey <- list(plot = TRUE, 
                     side = side, 
                     side.clab = 1,
                     dist = dist)
    } else {
      colkey <- list(plot = FALSE, 
                     side = side, 
                     side.clab = 1,
                     dist = dist)
    }
  }

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
  
  scatter3D(x = Axis1[ordination$stands == "No"], y = Axis2[ordination$stands == "No"], z = Axis3[ordination$stands == "No"], 
            phi = phi, 
            theta = theta,
            bty = "u", 
            pch = 16, 
            cex = point_size, 
            ticktype = "detailed",
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            xlim = axis_range,
            ylim = axis_range,
            zlim = axis_range,
            nticks = 3,
            col = col[ordination$stands == "No"],
            colvar = colvar,
            colkey = colkey,
            clab = clab,
            alpha.col = 0.5,
            add = TRUE)
  
  scatter3D(x = Axis1[ordination$stands == "Yes"], y = Axis2[ordination$stands == "Yes"], z = Axis3[ordination$stands == "Yes"], 
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
            col = col[ordination$stands == "Yes"],
            colvar = colvar,
            colkey = colkey,
            clab = clab,
            alpha.col = alpha_color,
            add = TRUE)
  
}


#-------------------------------------------------------------------------------

ordination_arrow_3d <- function(ordination_difference,
                                title = "Change 2024–2025",
                                phi = 20,
                                theta = 45,
                                axis_range,
                                alpha = 0.8) {
  
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
  arrlength <- 0.05
  arrangle  <- 20
  null_colkey <- list(plot = FALSE, 
                      side = side, 
                      side.clab = 1,
                      dist = dist)
  axis_leter <- 1
  axis_label <- 1.2
  
  dt <- ordination_difference[stands == "Yes",]
  
  Axis1_2024 <- clamp(dt$Axis1_2024 / 10000, -0.5, 0.5)
  Axis2_2024 <- clamp(dt$Axis2_2024 / 10000, -0.5, 0.5)
  Axis3_2024 <- clamp(dt$Axis3_2024 / 10000, -0.5, 0.5)
  
  Axis1_2025 <- clamp(dt$Axis1_2025 / 10000, -0.5, 0.5)
  Axis2_2025 <- clamp(dt$Axis2_2025 / 10000, -0.5, 0.5)
  Axis3_2025 <- clamp(dt$Axis3_2025 / 10000, -0.5, 0.5)
  
  col <- RGB(Axis1_2025, 
             Axis2_2025, 
             Axis3_2025,
             min_range = axis_range[1], 
             max_range = axis_range[2])
  
  
  xlab <- "Axis 1"
  ylab <- "Axis 2"
  zlab <- "Axis 3"

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
            #col.panel = "black",
            col.grid = "grey75",
            col.axis = "white",
            main = title,
            cex.axis = axis_leter,
            cex.lab = axis_label,
            scale = FALSE,
            d = 1,
            zoom = 1)
  
  arrows3D(x0 = rep(axis_range[1], length(Axis1_2024)), y0 = Axis2_2024, z0 = Axis3_2024,
           x1 = rep(axis_range[1], length(Axis1_2025)), y1 = Axis2_2025, z1 = Axis3_2025,
           phi = phi,
           theta = theta,
           col = shade_color,
           alpha = 1,
           add = TRUE,
           length = arrlength,
           angle  = arrangle,
           lwd = 0.1)

  arrows3D(x0 = Axis1_2024, y0 = rep(axis_range[2], length(Axis2_2024)), z0 = Axis3_2024,
           x1 = Axis1_2025, y1 = rep(axis_range[2], length(Axis2_2025)), z1 = Axis3_2025,
           phi = phi,
           theta = theta,
           col = shade_color,
           alpha = 1,
           add = TRUE,
           length = arrlength,
           angle  = arrangle,
           lwd = 0.01)

  arrows3D(x0 = Axis1_2024, y0 = Axis2_2024, z0 = rep(axis_range[1], length(Axis3_2024)),
           x1 = Axis1_2025, y1 = Axis2_2025, z1 = rep(axis_range[1], length(Axis3_2025)),
           phi = phi,
           theta = theta,
           col = shade_color,
           alpha = 1,
           add = TRUE,
           length = arrlength,
           angle  = arrangle,
           lwd = 0.1)

  # Halo pass
  arrows3D(x0 = Axis1_2024, y0 = Axis2_2024, z0 = Axis3_2024,
           x1 = Axis1_2025, y1 = Axis2_2025, z1 = Axis3_2025,
           phi = phi,
           theta = theta,
           col = "black",
           alpha = 0.4,
           add = TRUE,
           length = arrlength * 1.2,
           angle  = arrangle,
           lwd = point_size + 0.4)

  # Colored pass on top
  arrows3D(x0 = Axis1_2024, y0 = Axis2_2024, z0 = Axis3_2024,
           x1 = Axis1_2025, y1 = Axis2_2025, z1 = Axis3_2025,
           phi = phi,
           theta = theta,
           col = col,
           alpha = 1,
           add = TRUE,
           length = arrlength,
           angle  = arrangle,
           lwd = point_size)
  
}

#-------------------------------------------------------------------------------
# RGB

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
ordination_3d(phylo2024, title = "2024", cRGB = TRUE, axis_range = observed_range)
ordination_3d(phylo2025, title = "2025", cRGB = TRUE, axis_range = observed_range)
ordination_arrow_3d(phylodifference, title = "Change 2024–2025", axis_range = observed_range)

dev.off()
