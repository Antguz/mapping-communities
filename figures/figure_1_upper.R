#'------------------------------------------------------------------------------
#' @title Three-dimensional figure for NMDS
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

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"

# Communities
taxo_comm <- fread(paste0(root_path, "/Diversity/02_LMDS/MDS_scores_comm_taxo_total.csv"))
phylo_comm <- fread(paste0(root_path, "/Diversity/02_LMDS/MDS_scores_comm_phylo_total.csv"))
func_comm <- fread(paste0(root_path, "/Diversity/02_LMDS/MDS_scores_comm_func_total.csv"))

# Axis range
observed_range <- c(-0.6, 0.6)

# Get domains
domain <- fread(paste0(root_path, "/00_FIA/", "plot_domain.csv"))
mat <- match.integer64(taxo_comm$PLT_CN, domain$PLT_CN)
domain <- domain[mat, ]
all.equal(domain$PLT_CN, taxo_comm$PLT_CN)

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
    col <- RGB(ordination$Axis1, 
               ordination$Axis2, 
               ordination$Axis3,
               min_range = axis_range[1], 
               max_range = axis_range[2])
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
  
  
  scatter3D(x = rep(axis_range[1], length(ordination$Axis1)), y = ordination$Axis2, z = ordination$Axis3, 
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
            cex.lab = axis_label)
  
  scatter3D(x = ordination$Axis1, y = rep(axis_range[2], length(ordination$Axis2)), z = ordination$Axis3, 
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
  
  scatter3D(x = ordination$Axis1, y = ordination$Axis2, z = rep(axis_range[1], length(ordination$Axis3)), 
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
  
  scatter3D(x = ordination$Axis1, y = ordination$Axis2, z = ordination$Axis3, 
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
            col = col,
            colvar = colvar,
            colkey = colkey,
            clab = clab,
            alpha.col = alpha_color,
            add = TRUE)
  
}

#-------------------------------------------------------------------------------
# RGB

png(paste0(root_path, "/figures/Manuscript_white/", "MDS_RGB.png"),
    width = 210, 
    height = 90, 
    units = "mm", 
    pointsize = 12,
    bg = "white",  
    res = 600)

par(mfrow=c(1,3), mar = c(0,1.1,1,0), oma = c(0,0,0,0))
ordination_3d(taxo_comm, title = "Taxonomic", cRGB = TRUE, axis_range = observed_range)
ordination_3d(phylo_comm, title = "Phylogenetic", cRGB = TRUE, axis_range = observed_range)
ordination_3d(func_comm, title = "Functional", cRGB = TRUE, axis_range = observed_range)

dev.off()

#-------------------------------------------------------------------------------
# Proportion of angiosperms

par(mar = c(1,1,1,1))
par(oma = c(0,0,0,0))

png(paste0(root_path, "/figures/Manuscript_black/", "MDS_angiosperms.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "Fall"),
              colvar = proportions$proportion_angiosperm,
              key = FALSE,
              clab = "Proportion of angiosperms")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "Fall"),
              colvar = proportions$proportion_angiosperm,
              key = TRUE,
              clab = "Proportion of angiosperms")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "Fall"),
              colvar = proportions$proportion_angiosperm,
              key = FALSE,
              clab = "Proportion of angiosperms")

dev.off()
