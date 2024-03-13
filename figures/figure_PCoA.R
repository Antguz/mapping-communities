#'------------------------------------------------------------------------------
#' @title Three-dimensional figure of PCoA
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library
library(bit64)
library(plot3D)
library(data.table)
library(schemr)
library(stringr)

#'------------------------------------------------------------------------------
#' Read data

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"

taxo_comm <- fread(paste0(root_path, "/05_MDS/PCoA_scores_comm_taxo.csv"))
taxo_comm <- taxo_comm[order(PCoA1, -PCoA2, PCoA3), ]
phylo_comm <- fread(paste0(root_path, "/05_MDS/PCoA_scores_comm_phylo-log.csv"))
phylo_comm <- phylo_comm[order(PCoA1, -PCoA2, PCoA3), ]
func_comm <- fread(paste0(root_path, "/05_MDS/PCoA_scores_comm_func.csv"))
func_comm <- func_comm[order(PCoA1, -PCoA2, PCoA3), ]

PCoA_taxo <- readRDS(paste0(root_path, "/05_MDS/PCoA_Btotal_taxo.rds"))
eigen_taxo <- PCoA_taxo$eig
eigen_taxo <- round(eigen_taxo/sum(eigen_taxo), 2) 

PCoA_phylo <- readRDS(paste0(root_path, "/05_MDS/PCoA_Btotal_phylo-log.rds"))
eigen_phylo <- PCoA_phylo$eig
eigen_phylo <- round(eigen_phylo/sum(eigen_phylo) * 100, 2) 

PCoA_func <- readRDS(paste0(root_path, "/05_MDS/PCoA_Btotal_func.rds"))
eigen_func <- PCoA_func$eig
eigen_func <- round(eigen_func/sum(eigen_func) * 100, 2) 


taxo_sp <- fread(paste0(root_path, "/05_MDS/scores_sp_taxo.csv"))
taxo_sp <- taxo_sp[order(PCoA1, -PCoA2, PCoA3), ]
phylo_sp <- fread(paste0(root_path, "/05_MDS/scores_sp_phylo-log.csv"))
phylo_sp <- phylo_sp[order(PCoA1, -PCoA2, PCoA3), ]
func_sp <- fread(paste0(root_path, "/05_MDS/scores_sp_func.csv"))
func_sp <- func_sp[order(PCoA1, -PCoA2, PCoA3), ]

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
                          variance, 
                          title = "Taxonomic", 
                          phi = 20, 
                          theta = 45,
                          comm = TRUE,
                          axis_range = c(-0.66, 0.66)) {
  
  if(comm == TRUE) {
    
    point_shape <- 20
    point_size_shade <- 0.1
    point_size <- 0.1
    #xlab <- paste0("Axis 1 (", variance[1], "%)")
    #ylab <- paste0("Axis 2 (", variance[2], "%)")
    #zlab <- paste0("Axis 3 (", variance[3], "%)")
    xlab <- paste0("Axis 1")
    ylab <- paste0("Axis 2")
    zlab <- paste0("Axis 3")
    shade_color <- "grey75"
    col_RGB <- RGB(ordination$PCoA1, ordination$PCoA2, ordination$PCoA3)
    alpha_color <- 0.5
    
    
  } else if(comm == FALSE) {
    
    point_shape <- 18
    point_size_shade <- 0.5
    point_size <- 0.5
    xlab <- paste0("Axis 1")
    ylab <- paste0("Axis 2")
    zlab <- paste0("Axis 3")
    shade_color <- "grey75"
    col_RGB <- RGB(ordination$PCoA1, ordination$PCoA2, ordination$PCoA3)
    alpha_color <- 1.0
    
  }
  
  scatter3D(x = rep(-0.60, length(ordination$PCoA1)), y = ordination$PCoA2, z = ordination$PCoA3, 
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
            alpha.col = 0.2,
            col.panel = "black",
            col.grid = "white",
            col.axis = "black",
            main = title,
            cex.axis = 0.65,
            cex.lab = 0.8)
  
  scatter3D(x = ordination$PCoA1, y = rep(0.60, length(ordination$PCoA2)), z = ordination$PCoA3, 
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
            add = TRUE,
            alpha.col = 0.2,
            col.panel = "black",
            col.grid = "white",
            col.axis = "black")
  
  scatter3D(x = ordination$PCoA1, y = ordination$PCoA2, z = rep(-0.60, length(ordination$PCoA3)), 
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
            add = TRUE,
            alpha.col = 0.2,
            col.panel = "black",
            col.grid = "white",
            col.axis = "black")
  
  scatter3D(x = ordination$PCoA1, y = ordination$PCoA2, z = ordination$PCoA3, 
            colvar = NULL,
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
            col = col_RGB,
            add = TRUE,
            alpha = alpha_color,
            col.panel = "black",
            col.grid = "white",
            col.axis = "black")
  
}


jpeg(paste0(root_path, "/figures/", "3D-ordination.jpeg"),
    width = 210, 
    height = 160, 
    units = "mm", 
    pointsize = 12,
    bg = "white",  
    res = 600)

par(mfrow=c(2,3))
ordination_3d(taxo_comm, eigen_taxo, title = "Taxonomic")
ordination_3d(phylo_comm, eigen_phylo, title = "Phylogenetic")
ordination_3d(func_comm, eigen_func, title = "Functional")
ordination_3d(taxo_sp, eigen_func, title = "Communities", comm = FALSE)
ordination_3d(phylo_sp, eigen_func, title = "Species", comm = FALSE)
ordination_3d(func_sp, eigen_func, title = " ", comm = FALSE)
dev.off()



ordination_3d(func_sp, eigen_func, title = "Functional", comm = FALSE)



taxo_comm <- fread(paste0(root_path, "/05_MDS/PCoA_scores_comm_taxo.csv"))
phylo_comm <- fread(paste0(root_path, "/05_MDS/PCoA_scores_comm_phylo-log.csv"))
func_comm <- fread(paste0(root_path, "/05_MDS/PCoA_scores_comm_func.csv"))
plot <- taxo_comm$PLT_CN

match_phylo <- match.integer64(plot, phylo_comm$PLT_CN)
match_func <- match.integer64(plot, func_comm$PLT_CN)
phylo_comm <- phylo_comm[match_phylo, ]
func_comm <- func_comm[match_func, ]


scatter3D(x = taxo_comm$PCoA1, 
          y = phylo_comm$PCoA1, 
          z = func_comm$PCoA1, 
          colvar = NULL,
          phi = 20, 
          theta = 315,
          bty = "g", 
          pch = 20, 
          cex = 0.5, 
          ticktype = "detailed",
          xlab = "Taxonomimc",
          ylab = "Phylogenetic",
          zlab = "Functional",
          xlim = c(-0.60, 0.60),
          ylim = c(-0.60, 0.60),
          zlim = c(-0.60, 0.60),
          nticks = 3,
          col = RGB(taxo_comm$PCoA1, phylo_comm$PCoA1, func_comm$PCoA1),
          alpha = 1.0,
          col.panel = "black",
          col.grid = "white",
          col.axis = "black",
          shade = TRUE,
          cex.axis = 0.75,
          cex.lab = 1)


taxo_sp <- fread(paste0(root_path, "/05_MDS/scores_sp_taxo.csv"))
taxo_sp <- taxo_sp[order(Taxon2)]
phylo_sp <- fread(paste0(root_path, "/05_MDS/scores_sp_phylo-log.csv"))
phylo_sp <- phylo_sp[order(Taxon2)]
func_sp <- fread(paste0(root_path, "/05_MDS/scores_sp_func.csv"))
func_sp <- func_sp[order(Taxon2)]


scatter3D(x = taxo_sp$PCoA1, 
          y = phylo_sp$PCoA1, 
          z = func_sp$PCoA1, 
          colvar = NULL,
          phi = 20, 
          theta = 45,
          bty = "g", 
          pch = 20, 
          cex = 0.5, 
          ticktype = "detailed",
          xlab = "Taxonomimc",
          ylab = "Phylogenetic",
          zlab = "Functional",
          xlim = c(-0.60, 0.60),
          ylim = c(-0.60, 0.60),
          zlim = c(-0.60, 0.60),
          nticks = 3,
          col = RGB(taxo_sp$PCoA1, phylo_sp$PCoA1, func_sp$PCoA1),
          alpha = 1.0,
          col.panel = "black",
          col.grid = "white",
          col.axis = "black",
          shade = TRUE,
          cex.axis = 0.75,
          cex.lab = 1)
