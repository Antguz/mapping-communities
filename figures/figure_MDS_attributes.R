#'------------------------------------------------------------------------------
#' @title Three-dimensional figure for MDS post processing
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(bit64)
library(plot3D)
library(data.table)
library(stringr)

#'------------------------------------------------------------------------------
#' Source code

source("figures/RGB.R")

#'------------------------------------------------------------------------------
#' @Plot-functions
#-------------------------------------------------------------------------------


application <- function(ordination,
                        title = " ",
                        phi = 20, 
                        theta = 45,
                        cRGB = TRUE,
                        col = NULL,
                        colvar = NULL,
                        key = TRUE,
                        clab = " ",
                        #axis_range = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
                        axis_range = c(-1.06, 1.06),
                        alpha = 1.0) {
  
  # Plot features
  point_shape <- 20
  point_size_shade <- 0.8
  point_size <- 1.5
  xlab <- paste0("Axis 1")
  ylab <- paste0("Axis 2")
  zlab <- paste0("Axis 3")
  shade_color <- "grey75"
  alpha_color <- 0.75
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

  arrows3D(x0 = rep(axis_range[1], length(ordination$Axis1)), 
           y0 = rep(0, nrow(ordination)), 
           z0 = rep(0, nrow(ordination)), 
           x1 = rep(axis_range[1], length(ordination$Axis1)), 
           y1 = ordination$Axis2, 
           z1 = ordination$Axis3, 
           phi = phi,
           theta = theta,
           bty = "u", 
           pch = 16, 
           cex = point_size_shade, 
           #ticktype = "detailed",
           xlab = xlab,
           ylab = ylab,
           zlab = zlab,
           xlim = axis_range,
           ylim = axis_range,
           zlim = axis_range,
           nticks = 3,
           col = "grey75",
           colvar = NULL,
           colkey = null_colkey,
           clab = clab,
           alpha.col = 1.0,
           #col.panel = "black",
           col.grid = "white",
           col.axis = "black",
           main = title,
           cex.axis = axis_leter,
           cex.lab = axis_label)
  
  arrows3D(x0 = rep(0, nrow(ordination)), 
           y0 = rep(axis_range[2], length(ordination$Axis1)), 
           z0 = rep(0, nrow(ordination)), 
           x1 = ordination$Axis1, 
           y1 = rep(axis_range[2], length(ordination$Axis1)), 
           z1 = ordination$Axis3, 
           phi = phi,
           theta = theta,
           bty = "u", 
           pch = 16, 
           cex = point_size_shade, 
           #ticktype = "detailed",
           xlab = xlab,
           ylab = ylab,
           zlab = zlab,
           xlim = axis_range,
           ylim = axis_range,
           zlim = axis_range,
           nticks = 3,
           col = "grey75",
           colvar = NULL,
           colkey = null_colkey,
           clab = clab,
           alpha.col = 1.0,
           #col.panel = "black",
           col.grid = "white",
           col.axis = "black",
           main = title,
           cex.axis = axis_leter,
           cex.lab = axis_label,
           add = TRUE)
  
  arrows3D(x0 = rep(0, nrow(ordination)), 
           y0 = rep(0, nrow(ordination)), 
           z0 = rep(axis_range[1], length(ordination$Axis1)), 
           x1 = ordination$Axis1, 
           y1 = ordination$Axis2, 
           z1 = rep(axis_range[1], length(ordination$Axis1)), 
           phi = phi,
           theta = theta,
           bty = "u", 
           pch = 16, 
           cex = point_size_shade, 
           #ticktype = "detailed",
           xlab = xlab,
           ylab = ylab,
           zlab = zlab,
           xlim = axis_range,
           ylim = axis_range,
           zlim = axis_range,
           nticks = 3,
           col = "grey75",
           colvar = NULL,
           colkey = null_colkey,
           clab = clab,
           alpha.col = 1.0,
           #col.panel = "black",
           col.grid = "white",
           col.axis = "black",
           main = title,
           cex.axis = axis_leter,
           cex.lab = axis_label,
           add = TRUE)
  
  arrows3D(x0 = rep(0, nrow(ordination)), 
           y0 = rep(0, nrow(ordination)), 
           z0 = rep(0, nrow(ordination)), 
           x1 = ordination$Axis1, 
           y1 = ordination$Axis2, 
           z1 = ordination$Axis3, 
           phi = phi,
           theta = theta,
           bty = "u", 
           pch = 16, 
           cex = point_size_shade, 
           #ticktype = "detailed",
           xlab = xlab,
           ylab = ylab,
           zlab = zlab,
           xlim = axis_range,
           ylim = axis_range,
           zlim = axis_range,
           nticks = 3,
           col = col,
           colvar = NULL,
           colkey = null_colkey,
           clab = clab,
           alpha.col = 1.0,
           #col.panel = "black",
           col.grid = "white",
           col.axis = "black",
           main = title,
           cex.axis = axis_leter,
           cex.lab = axis_label,
           add = TRUE)
  
  text3D(x = ordination$Axis1, 
         y = ordination$Axis2, 
         z = ordination$Axis3+0.1,
         labels = ordination$variable,
         col = rep("red", nrow(ordination)),
         add=TRUE, 
         colkey = FALSE,
         cex = 0.8)
  
}
  

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"
root_path <- "E:/MMDB"

#'------------------------------------------------------------------------------
#' @Plots
#-------------------------------------------------------------------------------

observed_range <- c(-1.06, 1.06)

# Traits
taxo_traits <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_traits-wascores_", "taxo", "_total.csv"))
phylo_traits <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_traits-wascores_", "phylo", "_total.csv"))
func_traits <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_traits-wascores_", "func", "_total.csv"))
names <- c("H", "S", "WD", "STI", "DTI", "WTI", "LMA", "N")
taxo_traits$variable <- names
phylo_traits$variable <- names
func_traits$variable <- names

# Environment
taxo_env <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_env-wascores_", "taxo", "_total.csv"))
phylo_env <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_env-wascores_", "phylo", "_total.csv"))
func_env <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_env-wascores_", "func", "_total.csv"))

# Lineages
taxo_line <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_lineages-wascores_", "taxo", "_total.csv"))
phylo_line <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_lineages-wascores_", "phylo", "_total.csv"))
func_line <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_lineages-wascores_", "func", "_total.csv"))
taxo_line <- taxo_line[c(146, 145, 149, 147),]
phylo_line <- phylo_line[c(146, 145, 149, 147),]
func_line <- func_line[c(145, 146, 147, 148),]
names <- c("Gymno.", "Angio.", "ECM", "AM")
taxo_line$variable <- names
phylo_line$variable <- names
func_line$variable <- names

png(paste0(root_path, "/figures/Manuscript_black/", "MDS_environment.png"),
    width = 210, 
    height = 210, 
    units = "mm", 
    pointsize = 12,
    bg = "white",  
    res = 600)

par(mfrow=c(3,3), mar = c(0,1.1,1,0), oma = c(0,0,0,0))
application(taxo_env, title = "Taxonomic", cRGB = TRUE, axis_range = observed_range)
application(phylo_env, title = "Phylogenetic", cRGB = TRUE, axis_range = observed_range)
application(func_env, title = "Functional", cRGB = TRUE, axis_range = observed_range)

application(taxo_traits, title = " ", cRGB = TRUE, axis_range = observed_range)
application(phylo_traits, title = " ", cRGB = TRUE, axis_range = observed_range)
application(func_traits, title = " ", cRGB = TRUE, axis_range = observed_range)

application(taxo_line, title = " ", cRGB = TRUE, axis_range = c(-0.5, 0.5))
application(phylo_line, title = " ", cRGB = TRUE, axis_range = c(-0.5, 0.5))
application(func_line, title = " ", cRGB = TRUE, axis_range = c(-0.5, 0.5))

dev.off()
