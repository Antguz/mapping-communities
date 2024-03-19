#'------------------------------------------------------------------------------
#' @title Three-dimensional figure of Axis for RGB
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library
library(bit64)
library(plot3D)
library(data.table)
library(schemr)
library(stringr)
library(rcartocolor)
library(viridis)
library(RColorBrewer)

#'------------------------------------------------------------------------------
#' Read data

#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
root_path <- "E:/Projects/DESIS_mapping"

# Communities
taxo_comm <- fread(paste0(root_path, "/05_NMDS/NMDS_scores_comm_taxo.csv"))
#taxo_comm <- taxo_comm[order(Axis1, -Axis2, Axis3), ]
phylo_comm <- fread(paste0(root_path, "/05_NMDS/NMDS_scores_comm_phylo-log.csv"))
#phylo_comm <- phylo_comm[order(Axis1, -Axis2, Axis3), ]
func_comm <- fread(paste0(root_path, "/05_NMDS/NMDS_scores_comm_func.csv"))
#func_comm <- func_comm[order(Axis1, -Axis2, Axis3), ]

# Get proportions
proportions <- fread(paste0(root_path, "/00_FIA/proportions.csv"))
plot <- taxo_comm$PLT_CN
match_proportions <- match.integer64(plot, proportions$PLT_CN)
proportions <- proportions[match_proportions,]

# Species
#taxo_sp <- fread(paste0(root_path, "/05_MDS/scores_sp_taxo.csv"))
#taxo_sp <- taxo_sp[order(Axis1, -Axis2, Axis3), ]
#phylo_sp <- fread(paste0(root_path, "/05_MDS/scores_sp_phylo-log.csv"))
#phylo_sp <- phylo_sp[order(Axis1, -Axis2, Axis3), ]
#func_sp <- fread(paste0(root_path, "/05_MDS/scores_sp_func.csv"))
#func_sp <- func_sp[order(Axis1, -Axis2, Axis3), ]

#'------------------------------------------------------------------------------
#' Functions
RGB <- function(X, Y, Z, min_range = -0.66, max_range = 0.66) {
  
  x_new <-  ((X - min_range) / (max_range - min_range)) * (255 - 0) + 0
  y_new <-  ((Y - min_range) / (max_range - min_range)) * (255 - 0) + 0
  z_new <-  ((Z - min_range) / (max_range - min_range)) * (255 - 0) + 0
  
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

ordination_3d <- function(ordination, 
                          title = " ", 
                          phi = 20, 
                          theta = 45,
                          cRGB = TRUE,
                          col = NULL,
                          colvar = NULL,
                          key = TRUE,
                          clab = " ",
                          axis_range = c(-0.55, 0.55)) {
  
  # Plot features
  point_shape <- 20
  point_size_shade <- 0.1
  point_size <- 0.1
  xlab <- paste0("Axis 1")
  ylab <- paste0("Axis 2")
  zlab <- paste0("Axis 3")
  shade_color <- "grey75"
  alpha_color <- 0.5
  
  if(cRGB == TRUE) {
    col <- RGB(ordination$Axis1, 
               ordination$Axis2, 
               ordination$Axis3,
               min_range = axis_range[1], 
               max_range = axis_range[2])
    colvar <- NULL
    colkey <- list(plot = FALSE, side = 1)
    clab <- NULL
    
  } else if(cRGB == FALSE) {
    if(key == TRUE) {
      colkey <- list(plot = TRUE, side = 1)
    } else {
      colkey <- list(plot = FALSE, side = 1)
    }
  }
  
  
  scatter3D(x = rep(axis_range[1], length(ordination$Axis1)), y = ordination$Axis2, z = ordination$Axis3, 
            phi = phi,
            theta = theta,
            bty = "g", 
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
            colkey = list(plot = FALSE, side = 1),
            clab = clab,
            alpha.col = 0.2,
            #col.panel = "black",
            #col.grid = "white",
            #col.axis = "black",
            main = title,
            cex.axis = 0.65,
            cex.lab = 0.8)
  
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
            colkey = list(plot = FALSE, side = 1),
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
            colkey = list(plot = FALSE, side = 1),
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

# RGB
png(paste0(root_path, "/figures/", "ordination_RGB.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, title = "Taxonomic", cRGB = TRUE)
ordination_3d(phylo_comm, title = "Phylogenetic", cRGB = TRUE)
ordination_3d(func_comm, title = "Functional", cRGB = TRUE)

dev.off()

# Proportion of angiosperms
png(paste0(root_path, "/figures/", "proportion_angiosperms.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = viridis(100),
              colvar = proportions$proportion_angiosperm,
              key = FALSE,
              clab = "Proportion of angiosperms")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = viridis(100),
              colvar = proportions$proportion_angiosperm,
              key = TRUE,
              clab = "Proportion of angiosperms")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = viridis(100),
              colvar = proportions$proportion_angiosperm,
              key = FALSE,
              clab = "Proportion of angiosperms")

dev.off()








