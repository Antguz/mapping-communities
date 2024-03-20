#'------------------------------------------------------------------------------
#' @title Three-dimensional figure for NMDS
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
library(colorspace)

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
func_comm$Axis3 <- -func_comm$Axis3
#func_comm <- func_comm[order(Axis1, -Axis2, Axis3), ]
#fwrite(func_comm, paste0(root_path, "/05_NMDS/NMDS_scores_comm_func.csv"))


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
                          axis_range = c(-0.55, 0.55),
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
png(paste0(root_path, "/figures/", "NMDS_RGB.png"),
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

#-------------------------------------------------------------------------------
# Proportion of angiosperms
png(paste0(root_path, "/figures/", "NMDS_angiosperms.png"),
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


#-------------------------------------------------------------------------------
# Proportion of AM
png(paste0(root_path, "/figures/", "NMDS_AM.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "Tropic"),
              colvar = proportions$proportion_AM,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "Tropic"),
              colvar = proportions$proportion_AM,
              key = TRUE,
              clab = "Proportion of AM-associated trees")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "Tropic"),
              colvar = proportions$proportion_AM,
              key = FALSE,
              clab = " ")

dev.off()

#-------------------------------------------------------------------------------
# Proportion of broadleaf
png(paste0(root_path, "/figures/", "NMDS_broadleaf.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = rev(carto_pal(n = 100, "Earth")),
              colvar = proportions$proportion_broadleaf,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = rev(carto_pal(n = 100, "Earth")),
              colvar = proportions$proportion_broadleaf,
              key = TRUE,
              clab = "Proportion of broadleaf trees")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = rev(carto_pal(n = 100, "Earth")),
              colvar = proportions$proportion_broadleaf,
              key = FALSE,
              clab = " ")

dev.off()

#-------------------------------------------------------------------------------
# Proportion of compound leaf
png(paste0(root_path, "/figures/", "NMDS_compound.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "TealRose"),
              colvar = proportions$proportion_compound,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "TealRose"),
              colvar = proportions$proportion_compound,
              key = TRUE,
              clab = "Proportion of compound-leaf trees")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "TealRose"),
              colvar = proportions$proportion_compound,
              key = FALSE,
              clab = " ")

dev.off()


#-------------------------------------------------------------------------------
# Proportion of deciduous trees
png(paste0(root_path, "/figures/", "NMDS_deciduous.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "ArmyRose"),
              colvar = proportions$proportion_deciduous,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "ArmyRose"),
              colvar = proportions$proportion_deciduous,
              key = TRUE,
              clab = "Proportion of deciduous trees")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = carto_pal(n = 100, "ArmyRose"),
              colvar = proportions$proportion_deciduous,
              key = FALSE,
              clab = " ")

dev.off()

#-------------------------------------------------------------------------------
# CMW of tree height
png(paste0(root_path, "/figures/", "NMDS_HTq.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "RdYlGn"),
              colvar = proportions$CWM_HTq,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "RdYlGn"),
              colvar = proportions$CWM_HTq,
              key = TRUE,
              clab = "CWM of maximun tree height (m)")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "RdYlGn"),
              colvar = proportions$CWM_HTq,
              key = FALSE,
              clab = " ")

dev.off()

#-------------------------------------------------------------------------------
# CMW of tree slenderness
png(paste0(root_path, "/figures/", "NMDS_SL.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "RdYlBu"),
              colvar = proportions$CWM_SLq,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "RdYlBu"),
              colvar = proportions$CWM_SLq,
              key = TRUE,
              clab = "CWM tree slenderness")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "RdYlBu"),
              colvar = proportions$CWM_SLq,
              key = FALSE,
              clab = " ")

dev.off()

#-------------------------------------------------------------------------------
# CMW of CCLCD
png(paste0(root_path, "/figures/", "NMDS_CCLCD.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "RdGy"),
              colvar = proportions$CWM_CCLCD_mode,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "RdGy"),
              colvar = proportions$CWM_CCLCD_mode,
              key = TRUE,
              clab = "CWM of CCLCD")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "RdGy"),
              colvar = proportions$CWM_CCLCD_mode,
              key = FALSE,
              clab = " ")

dev.off()


#-------------------------------------------------------------------------------
# CMW of gravity
png(paste0(root_path, "/figures/", "NMDS_gravity.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "BrBG"),
              colvar = proportions$CWM_WOOD_SPGR_GREENVOL_DRYWT,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "BrBG"),
              colvar = proportions$CWM_WOOD_SPGR_GREENVOL_DRYWT,
              key = TRUE,
              clab = "CWM of specific gravity of wood")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "BrBG"),
              colvar = proportions$CWM_WOOD_SPGR_GREENVOL_DRYWT,
              key = FALSE,
              clab = " ")

dev.off()


#-------------------------------------------------------------------------------
# CMW of moisture content of green wood
png(paste0(root_path, "/figures/", "NMDS_moisture-green.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "PRGn"),
              colvar = proportions$CWM_MC_PCT_GREEN_WOOD,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "PRGn"),
              colvar = proportions$CWM_MC_PCT_GREEN_WOOD,
              key = TRUE,
              clab = "CWM of moisture content of green wood")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "PRGn"),
              colvar = proportions$CWM_MC_PCT_GREEN_WOOD,
              key = FALSE,
              clab = " ")

dev.off()


#-------------------------------------------------------------------------------
# CMW of wood carbon fraction
png(paste0(root_path, "/figures/", "NMDS_carbon-fraction.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "PuOr"),
              colvar = proportions$CWM_CARBON_RATIO_LIVE,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "PuOr"),
              colvar = proportions$CWM_CARBON_RATIO_LIVE,
              key = TRUE,
              clab = "CWM of wood carbon fraction")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = divergingx_hcl(100, palette = "PuOr"),
              colvar = proportions$CWM_CARBON_RATIO_LIVE,
              key = FALSE,
              clab = " ")

dev.off()

#-------------------------------------------------------------------------------
# CMW of shade tolerance
png(paste0(root_path, "/figures/", "NMDS_shade.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "magma", direction = 1),
              colvar = proportions$CWM_shade,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "magma", direction = 1),
              colvar = proportions$CWM_shade,
              key = TRUE,
              clab = "CWM of shade tolerance index")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "magma", direction = 1),
              colvar = proportions$CWM_shade,
              key = FALSE,
              clab = " ")

dev.off()


#-------------------------------------------------------------------------------
# CMW of drought tolerance
png(paste0(root_path, "/figures/", "NMDS_drought.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "inferno", direction = 1),
              colvar = proportions$CWM_drought,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "inferno", direction = 1),
              colvar = proportions$CWM_drought,
              key = TRUE,
              clab = "CWM of drought tolerance index")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "inferno", direction = 1),
              colvar = proportions$CWM_drought,
              key = FALSE,
              clab = " ")

dev.off()

#-------------------------------------------------------------------------------
# CMW of waterlogging tolerance
png(paste0(root_path, "/figures/", "NMDS_waterlogging.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "mako", direction = 1),
              colvar = proportions$CWM_waterlogging,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "mako", direction = 1),
              colvar = proportions$CWM_waterlogging,
              key = TRUE,
              clab = "CWM of waterlogging tolerance index")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "mako", direction = 1),
              colvar = proportions$CWM_waterlogging,
              key = FALSE,
              clab = " ")

dev.off()


#-------------------------------------------------------------------------------
# CMW of leaf area
png(paste0(root_path, "/figures/", "NMDS_leaf-area.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "cividis", direction = 1),
              colvar = proportions$CWM_leaf_area,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "cividis", direction = 1),
              colvar = proportions$CWM_leaf_area,
              key = TRUE,
              clab = "CWM of leaf area")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "cividis", direction = 1),
              colvar = proportions$CWM_leaf_area,
              key = FALSE,
              clab = " ")

dev.off()

#-------------------------------------------------------------------------------
# CMW of LMA
png(paste0(root_path, "/figures/", "NMDS_LMA.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "viridis", direction = 1),
              colvar = proportions$CWM_LMA,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "viridis", direction = 1),
              colvar = proportions$CWM_LMA,
              key = TRUE,
              clab = "CWM of LMA")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "viridis", direction = 1),
              colvar = proportions$CWM_LMA,
              key = FALSE,
              clab = " ")

dev.off()

#-------------------------------------------------------------------------------
# CMW of Genus
png(paste0(root_path, "/figures/", "NMDS_genus.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = viridis(length(unique(genus_max$Genus)), option = "turbo", direction = -1),
              colvar = as.integer(as.factor(genus_max$Genus)),
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = viridis(length(unique(genus_max$Genus)), option = "turbo", direction = -1),
              colvar = as.integer(as.factor(genus_max$Genus)),
              key = FALSE,
              clab = " ")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = viridis(length(unique(genus_max$Genus)), option = "turbo", direction = -1),
              colvar = as.integer(as.factor(genus_max$Genus)),
              key = FALSE,
              clab = " ")

dev.off()


#-------------------------------------------------------------------------------
# CMW of Families
png(paste0(root_path, "/figures/", "NMDS_families.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = viridis(length(unique(family_max$Family)), option = "turbo", direction = -1),
              colvar = as.integer(as.factor(family_max$Family)),
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = viridis(length(unique(family_max$Family)), option = "turbo", direction = -1),
              colvar = as.integer(as.factor(family_max$Family)),
              key = FALSE,
              clab = " ")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = viridis(length(unique(family_max$Family)), option = "turbo", direction = -1),
              colvar = as.integer(as.factor(family_max$Family)),
              key = FALSE,
              clab = " ")

dev.off()

#-------------------------------------------------------------------------------
# Domains

domain <- fread(paste0(root_path, "/00_FIA/", "plot_domain.csv"))
mat <- match.integer64(taxo_comm$PLT_CN, domain$PLT_CN)
domain <- domain[mat, ]

png(paste0(root_path, "/figures/", "NMDS_domains.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = qualitative_hcl(length(unique(domain$DomainName)), palette = "Dark3"),
              colvar = as.integer(as.factor(domain$DomainName)),
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = qualitative_hcl(length(unique(domain$DomainName)), palette = "Dark3"),
              colvar = as.integer(as.factor(domain$DomainName)),
              key = FALSE,
              clab = " ")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = qualitative_hcl(length(unique(domain$DomainName)), palette = "Dark3"),
              colvar = as.integer(as.factor(domain$DomainName)),
              key = FALSE,
              clab = " ")

dev.off()


#-------------------------------------------------------------------------------
# Elevation from sea level

domain <- fread(paste0(root_path, "/00_FIA/", "plot_domain.csv"))
mat <- match.integer64(taxo_comm$PLT_CN, domain$PLT_CN)
domain <- domain[mat, ]
all.equal(domain$PLT_CN, taxo_comm$PLT_CN)

png(paste0(root_path, "/figures/", "NMDS_elevation.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "turbo", direction = 1),
              colvar = domain$ELEV,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "turbo", direction = 1),
              colvar = domain$ELEV,
              key = TRUE,
              clab = "Elevation from sea level")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "turbo", direction = 1),
              colvar = domain$ELEV,
              key = FALSE,
              clab = " ")

dev.off()


#-------------------------------------------------------------------------------
# Latitude

domain <- fread(paste0(root_path, "/00_FIA/", "plot_domain.csv"))
mat <- match.integer64(taxo_comm$PLT_CN, domain$PLT_CN)
domain <- domain[mat, ]

png(paste0(root_path, "/figures/", "NMDS_latitude.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "turbo", direction = 1),
              colvar = domain$LAT,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "turbo", direction = 1),
              colvar = domain$LAT,
              key = TRUE,
              clab = "Latitude")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "turbo", direction = 1),
              colvar = domain$LAT,
              key = FALSE,
              clab = " ")

dev.off()


#-------------------------------------------------------------------------------
# Longitude

domain <- fread(paste0(root_path, "/00_FIA/", "plot_domain.csv"))
mat <- match.integer64(taxo_comm$PLT_CN, domain$PLT_CN)
domain <- domain[mat, ]

png(paste0(root_path, "/figures/", "NMDS_longitude.png"),
    width = 210, 
    height = 80, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

par(mfrow=c(1,3))
ordination_3d(taxo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "turbo", direction = 1),
              colvar = domain$LON,
              key = FALSE,
              clab = " ")
ordination_3d(phylo_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "turbo", direction = 1),
              colvar = domain$LON,
              key = TRUE,
              clab = "Longitude")
ordination_3d(func_comm, 
              cRGB = FALSE,
              col = viridis(100, option = "turbo", direction = 1),
              colvar = domain$LON,
              key = FALSE,
              clab = " ")

dev.off()
