#'------------------------------------------------------------------------------
#' @title Plot occurrence of plant lineages based on MDS
#'------------------------------------------------------------------------------

#' @description It plot occurrence based on MDS ordinations.

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(bit64)
library(parallel)
library(pbmcapply)
library(scales)
library(stringr)
library(plot3D)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/guzman/Work/MMDB"
type <- "phylo"
out_path <- "/media/guzman/Work/MMDB/Models/03_model_occurrence"

#'------------------------------------------------------------------------------
#' @Read-Information-plant_lineages
#-------------------------------------------------------------------------------

aggregation <- readRDS(paste0(out_path, "/LINEAGES_occurrance.rds"))
ordination <- readRDS(paste0(out_path, "/LINEAGES_variables.rds"))

#'------------------------------------------------------------------------------
#' @Plot-function
#-------------------------------------------------------------------------------

# RGB
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
ordination_3d <- function(X = NULL,
                          ordination, 
                          aggregation,
                          phi = 20, 
                          theta = 45,
                          cRGB = TRUE,
                          col = NULL,
                          colvar = NULL,
                          key = TRUE,
                          clab = " ",
                          axis_range = c(-0.6, 0.6),
                          alpha = 1.0,
                          main_title = NULL) {
  
  # Scale
  Axis1 <- ordination$Axis1
  Axis2 <- ordination$Axis2
  Axis3 <- ordination$Axis3
  
  
  if(is.null(X)) {
    keep <- rep(TRUE, length(Axis1))
    name <- if (is.null(main_title)) "All points" else main_title
    value <- rep(1, length(Axis1))
  
  } else {
    value <- aggregation[, .SD, .SDcols = X][[1]]
    keep <- value != 0
    name <- if (is.null(main_title)) colnames(aggregation)[X] else main_title
  }
  
  # Plot features
  point_shape <- 20
  point_size_shade <- 0.1
  xlab <- "Axis 1"
  ylab <- "Axis 2"
  zlab <- "Axis 3"
  shade_color <- "grey85"
  alpha_color <- alpha
  side <- 1
  dist <- 0
  
  null_colkey <- list(plot = FALSE,
                      side = side, 
                      side.clab = 1,
                      dist = dist)
  
  axis_leter <- 1
  axis_label <- 1.2
  
  if (cRGB == TRUE) {
    col <- RGB(Axis1, Axis2, Axis3,
               min_range = axis_range[1], 
               max_range = axis_range[2]    )
    
    if (!is.null(X)) {
      col[!keep] <- "grey40"
    }
    
    colvar <- NULL
    colkey <- list(plot = FALSE, 
                   side = side, 
                   side.clab = 1,
                   dist = dist)
    clab <- NULL
    
  } else {
    if (key == TRUE) {
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
  
  par(mar = c(1, 1, 1, 0))
  
  scatter3D(x = rep(axis_range[1], sum(keep)),
            y = Axis2[keep],
            z = Axis3[keep], 
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
            col.grid = "grey75",
            col.axis = "white",
            main = name,
            cex.axis = axis_leter,
            cex.lab = axis_label,
            scale = FALSE,
            d = 1,
            zoom = 1)
  
  scatter3D(x = Axis1[keep],
            y = rep(axis_range[2], sum(keep)),
            z = Axis3[keep], 
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
  
  scatter3D(x = Axis1[keep],
            y = Axis2[keep],
            z = rep(axis_range[1], sum(keep)), 
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
  
  scatter3D(x = Axis1[keep],
            y = Axis2[keep],
            z = Axis3[keep], 
            phi = phi, 
            theta = theta,
            bty = "u", 
            pch = point_shape, 
            cex = 0.1, 
            ticktype = "detailed",
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            xlim = axis_range,
            ylim = axis_range,
            zlim = axis_range,
            nticks = 3,
            col = if (cRGB) col[keep] else col,
            colvar = colvar,
            colkey = colkey,
            clab = clab,
            alpha.col = alpha_color,
            add = TRUE)
  
}

# Create plots
export_ordination_grid <- function(out_dir,
                                   ord_data,
                                   aggregation,
                                   X = seq_len(ncol(aggregation)),
                                   phi = 20,
                                   theta = 45,
                                   cRGB = TRUE,
                                   key = FALSE,
                                   axis_range = c(-0.6, 0.6),
                                   alpha = 1,
                                   page_width_mm = 210,
                                   page_height_mm = 297,
                                   plot_size_mm = 70,
                                   res = 300) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  width_in  <- page_width_mm / 25.4
  height_in <- page_height_mm / 25.4
  
  ncol_page <- floor(page_width_mm / plot_size_mm)
  nrow_page <- floor(page_height_mm / plot_size_mm)
  plots_per_page <- ncol_page * nrow_page
  
  # Add the first overview plot
  X_all <- c(NA, X)
  
  pages <- split(X_all, ceiling(seq_along(X_all) / plots_per_page))
  
  for (p in seq_along(pages)) {
    
    file <- file.path(out_dir, sprintf("ordination_page_%02d.jpg", p))
    
    jpeg(
      filename = file,
      width = width_in,
      height = height_in,
      units = "in",
      res = res,
      quality = 100
    )
    
    par(
      mfrow = c(nrow_page, ncol_page),
      mar = c(0.2, 0.2, 0.8, 0.2),
      oma = c(0, 0, 0, 0),
      pty = "s"
    )
    
    idx <- pages[[p]]
    
    for (i in idx) {
      if (is.na(i)) {
        ordination_3d(X = NULL,
                      ordination = ord_data,
                      aggregation = aggregation,
                      phi = phi,
                      theta = theta,
                      cRGB = cRGB,
                      key = key,
                      axis_range = axis_range,
                      alpha = alpha,
                      main_title = "All observations")
      } else {
        ordination_3d(X = i,
                      ordination = ord_data,
                      aggregation = aggregation,
                      phi = phi,
                      theta = theta,
                      cRGB = cRGB,
                      key = key,
                      axis_range = axis_range,
                      alpha = alpha)
      }
    }
    
    n_empty <- plots_per_page - length(idx)
    if (n_empty > 0) {
      for (j in seq_len(n_empty)) plot.new()
    }
    
    dev.off()
  }
}

export_ordination_grid(
  out_dir = "plant_lineages",
  ord_data = ordination,
  aggregation = aggregation[, c(-1, -30, -31)],
  X = seq_len(ncol(aggregation[, c(-1, -30, -31)])),
  res = 300
)


#'------------------------------------------------------------------------------
#' @Read-Information-forest_type
#-------------------------------------------------------------------------------

type <- "taxo"
out_path <- "/media/guzman/Work/MMDB/Models/03_models_fortyp"

aggregation <- readRDS(paste0(out_path, "/FORTYPCD_occurrance.rds"))
ordination <- readRDS(paste0(out_path, "/FORTYPCD_variables.rds"))
models <- fread(paste0(out_path, "/taxo_altogether_model-FORTYP.csv"))
forestypes <- unique(models$forty)

export_ordination_grid(
  out_dir = "forest_types",
  ord_data = ordination,
  aggregation = aggregation[, ..forestypes],
  X = seq_len(ncol(aggregation[, ..forestypes])),
  res = 300
)

#'------------------------------------------------------------------------------
#' @Read-Information-forest_type
#-------------------------------------------------------------------------------

type <- "func"
out_path <- "/media/guzman/Work/MMDB/Models/03_model_CWM"

aggregation <- readRDS(paste0(out_path, "/CWM_plant-traits.rds"))
ordination <- readRDS(paste0(out_path, "/CWM_variables.rds"))

variable_names <- c("PLT_CN", "Tree height", "Slenderness", "Wood density", "Shade tolerance index", "Drought tolerance index", 
                    "Waterlogging tolerance index", "Leaf mass per area", "Nitrogen")
colnames(aggregation) <- variable_names
units <- c("m", "-", "g/cm3", "-", "-", "-", "g/m2", "mg/g")

#'------------------------------------------------------------------------------
#' @Plot-function
#-------------------------------------------------------------------------------

# Ordination figure
ordination_3dm <- function(X = NULL,
                           ordination,
                           aggregation,
                           phi = 20,
                           theta = 45,
                           cRGB = TRUE,
                           col = NULL,
                           colvar = NULL,
                           key = TRUE,
                           clab = " ",
                           axis_range = c(-0.6, 0.6),
                           alpha = 1.0,
                           main_title = NULL,
                           color_X = NULL,
                           color_units = NULL,
                           palette_fun = function(n) viridisLite::viridis(n),
                           zlim_col = NULL) {
  
  # Coordinates
  Axis1 <- ordination$Axis1
  Axis2 <- ordination$Axis2
  Axis3 <- ordination$Axis3
  
  # Which points to display
  if(is.null(X)) {
    keep <- rep(TRUE, length(Axis1))
    name <- if (is.null(main_title)) "All points" else main_title
    value <- rep(1, length(Axis1))
  } else {
    value <- aggregation[[X]]
    keep <- value != 0
    name <- if (is.null(main_title)) colnames(aggregation)[X] else main_title
  }
  
  # Plot features
  point_shape <- 20
  point_size_shade <- 0.1
  point_size <- 0.1
  xlab <- "Axis 1"
  ylab <- "Axis 2"
  zlab <- "Axis 3"
  shade_color <- "grey85"
  alpha_color <- alpha
  side <- 4
  dist <- 0
  
  null_colkey <- list(plot = FALSE,
                      side = side,
                      side.clab = 1,
                      dist = dist)
  
  axis_leter <- 1
  axis_label <- 1.2

  if (!is.null(color_X)) {
    color_values <- aggregation[[color_X]]
    
    if (is.null(zlim_col)) {
      zlim_col <- range(color_values[keep], na.rm = TRUE)
    }
    
    colvar <- color_values[keep]
    col <- NULL
    
    if (is.null(color_units)) {
      clab <- colnames(aggregation)[color_X]
    } else {
      clab <- paste0(colnames(aggregation)[color_X], " (", color_units, ")")
    }
    
    colkey <- list(plot = key,
                   side = side,
                   side.clab = 1,
                   dist = dist)
    
  } else if (cRGB == TRUE) {
    col <- RGB(Axis1,
               Axis2,
               Axis3,
               min_range = axis_range[1],
               max_range = axis_range[2])
    
    if (!is.null(X)) {
      col[!keep] <- "grey40"
    }
    
    colvar <- NULL
    colkey <- list(plot = FALSE,
                   side = side,
                   side.clab = 1,
                   dist = dist)
    clab <- NULL
    
  } else {
    if (key == TRUE) {
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
  
  par(mar = c(1, 1, 1, 3))
  
  # Background projections
  scatter3D(x = rep(axis_range[1], sum(keep)),
            y = Axis2[keep],
            z = Axis3[keep],
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
            col.grid = "grey75",
            col.axis = "white",
            main = name,
            cex.axis = axis_leter,
            cex.lab = axis_label,
            scale = FALSE,
            d = 1,
            zoom = 1)
  
  scatter3D(x = Axis1[keep],
            y = rep(axis_range[2], sum(keep)),
            z = Axis3[keep],
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
  
  scatter3D(x = Axis1[keep],
            y = Axis2[keep],
            z = rep(axis_range[1], sum(keep)),
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
  
  # Main points
  if (!is.null(color_X)) {
    scatter3D(x = Axis1[keep],
              y = Axis2[keep],
              z = Axis3[keep],
              phi = phi,
              theta = theta,
              bty = "u",
              pch = point_shape,
              cex = point_size,
              ticktype = "detailed",
              xlab = xlab,
              ylab = ylab,
              zlab = zlab,
              xlim = axis_range,
              ylim = axis_range,
              zlim = axis_range,
              nticks = 3,
              col = palette_fun(100),
              colvar = colvar,
              colkey = colkey,
              clab = clab,
              alpha.col = alpha_color,
              zlim = zlim_col,
              add = TRUE)
    
  } else {
    scatter3D(x = Axis1[keep],
              y = Axis2[keep],
              z = Axis3[keep],
              phi = phi,
              theta = theta,
              bty = "u",
              pch = point_shape,
              cex = point_size,
              ticktype = "detailed",
              xlab = xlab,
              ylab = ylab,
              zlab = zlab,
              xlim = axis_range,
              ylim = axis_range,
              zlim = axis_range,
              nticks = 3,
              col = if (cRGB) col[keep] else col,
              colvar = colvar,
              colkey = colkey,
              clab = clab,
              alpha.col = alpha_color,
              add = TRUE)
  }
}

export_ordination_grid <- function(out_dir,
                                   ord_data,
                                   aggregation,
                                   X = seq_len(ncol(aggregation)),
                                   units_vec = NULL,
                                   phi = 20,
                                   theta = 45,
                                   cRGB = TRUE,
                                   key = FALSE,
                                   axis_range = c(-0.6, 0.6),
                                   alpha = 1,
                                   page_width_mm = 210,
                                   page_height_mm = 297,
                                   plot_size_mm = 70,
                                   res = 300,
                                   add_overview = TRUE,
                                   add_continuous_overview = FALSE,
                                   overview_color_X = NULL,
                                   overview_color_units = NULL,
                                   zlim_col = NULL) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  if (!is.null(units_vec) && length(units_vec) != length(X)) {
    stop("length(units_vec) must match length(X)")
  }
  
  width_in  <- page_width_mm / 25.4
  height_in <- page_height_mm / 25.4
  
  ncol_page <- floor(page_width_mm / plot_size_mm)
  nrow_page <- floor(page_height_mm / plot_size_mm)
  plots_per_page <- ncol_page * nrow_page
  
  # Build plot plan
  plot_plan <- list()
  
  if (add_overview) {
    plot_plan[[length(plot_plan) + 1]] <- list(
      type = "all_points"
    )
  }
  
  if (add_continuous_overview) {
    plot_plan[[length(plot_plan) + 1]] <- list(
      type = "continuous_overview",
      color_X = overview_color_X,
      color_units = overview_color_units
    )
  }
  
  for (k in seq_along(X)) {
    plot_plan[[length(plot_plan) + 1]] <- list(
      type = "filtered",
      X = X[k],
      units = if (is.null(units_vec)) NULL else units_vec[k]
    )
  }
  
  pages <- split(plot_plan, ceiling(seq_along(plot_plan) / plots_per_page))
  
  for (p in seq_along(pages)) {
    file <- file.path(out_dir, sprintf("ordination_page_%02d.jpg", p))
    
    jpeg(
      filename = file,
      width = width_in,
      height = height_in,
      units = "in",
      res = res,
      quality = 100
    )
    
    par(
      mfrow = c(nrow_page, ncol_page),
      mar = c(0.2, 0.2, 0.8, 0.2),
      oma = c(0, 0, 0, 0),
      pty = "s"
    )
    
    page_plots <- pages[[p]]
    
    for (obj in page_plots) {
      if (obj$type == "all_points") {
        ordination_3dm(
          X = NULL,
          ordination = ord_data,
          aggregation = aggregation,
          phi = phi,
          theta = theta,
          cRGB = TRUE,
          key = FALSE,
          axis_range = axis_range,
          alpha = alpha,
          main_title = "All points"
        )
      }
      
      if (obj$type == "continuous_overview") {
        ordination_3dm(
          X = NULL,
          ordination = ord_data,
          aggregation = aggregation,
          phi = phi,
          theta = theta,
          cRGB = FALSE,
          key = TRUE,
          axis_range = axis_range,
          alpha = alpha,
          main_title = colnames(aggregation)[obj$color_X],
          color_X = obj$color_X,
          color_units = obj$color_units,
          zlim_col = zlim_col
        )
      }
      
      if (obj$type == "filtered") {
        ordination_3dm(
          X = obj$X,
          ordination = ord_data,
          aggregation = aggregation,
          phi = phi,
          theta = theta,
          cRGB = FALSE,
          key = TRUE,
          axis_range = axis_range,
          alpha = alpha,
          color_X = obj$X,
          color_units = obj$units,
          zlim_col = zlim_col
        )
      }
    }
    
    n_empty <- plots_per_page - length(page_plots)
    if (n_empty > 0) {
      for (j in seq_len(n_empty)) plot.new()
    }
    
    dev.off()
  }
}

export_ordination_grid(
  out_dir = "cwm_plant_traits",
  ord_data = ordination,
  aggregation = aggregation,
  X = 2:9,
  units_vec = units,
  res = 300
)
