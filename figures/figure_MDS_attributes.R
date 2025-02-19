#'------------------------------------------------------------------------------
#' @title Three-dimensional figure for MDS post processing
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(bit64)
library(plot3D)
library(data.table)

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"

type <- "taxo"

# Communities
env <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_env-wascores_", type, "_total.csv"))
phylo_comm <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_scores_comm_", type, "_total.csv"))
func_comm <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_scores_comm_", type, "_total.csv"))


application <- function(X,
                        ordination,
                        phi = 20, 
                        theta = 45,
                        axis_range = c(-0.6, 0.6)) {
  
  # Plot features
  point_shape <- 20
  point_size_shade <- 0.1
  point_size <- 0.1
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
  
  
  scatter3D(x = rep(-max(abs(range(ordination[,-1]))), length(ordination$Axis1)), 
            y = ordination$Axis2, 
            z = ordination$Axis3, 
            phi = phi,
            theta = theta,
            bty = "u",
            pch = 16,
            cex = 0.5, 
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            xlim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
            ylim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
            zlim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
            col = "grey75",
            colvar = NULL,
            colkey = null_colkey,
            alpha.col = 0.85)
  
  scatter3D(x = ordination$Axis1, 
            y = rep(max(abs(range(ordination[,-1]))), length(ordination$Axis2)), 
            z = ordination$Axis3, 
            phi = phi,
            theta = theta,
            bty = "u",
            pch = 16,
            cex = 0.5, 
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            xlim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
            ylim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
            zlim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
            col = "grey75",
            colvar = NULL,
            colkey = null_colkey,
            alpha.col = 0.85,
            add = TRUE)
  
  scatter3D(x = ordination$Axis1, 
            y = ordination$Axis2, 
            z = rep(-max(abs(range(ordination[,-1]))), length(ordination$Axis3)), 
            phi = phi,
            theta = theta,
            bty = "u",
            pch = 16,
            cex = 0.5, 
            xlab = xlab,
            ylab = ylab,
            zlab = zlab,
            xlim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
            ylim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
            zlim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
            col = "grey75",
            colvar = NULL,
            colkey = null_colkey,
            alpha.col = 0.85,
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
           cex = 1.5, 
           xlab = xlab,
           ylab = ylab,
           zlab = zlab,
           xlim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
           ylim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
           zlim = c(-max(abs(range(ordination[,-1]))), max(abs(range(ordination[,-1])))),
           col = "grey25",
           colvar = NULL,
           colkey = null_colkey,
           alpha.col = 0.85,
           add = TRUE)
  
  text3D(x = ordination$Axis1, 
         y = ordination$Axis2, 
         z = ordination$Axis3+0.1,
         labels = ordination$variable,
         col = rep("red", nrow(ordination)),
         add=TRUE, 
         colkey = FALSE)
  
  
  
  