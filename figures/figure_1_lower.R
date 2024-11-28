#'------------------------------------------------------------------------------
#' @title Three-dimensional figure for NMDS
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(bit64)
library(plot3D)
library(ggdark)
library(data.table)
library(schemr)
library(stringr)
library(rcartocolor)
library(viridis)
library(colorspace)
library(Rcpp)
library(RcppArmadillo)
library(sf)
library(ggpubr)
library(cowplot)

#'------------------------------------------------------------------------------
#' @Function
#-------------------------------------------------------------------------------

RGB <- function(X, Y, Z, min_range = -0.50, max_range = 0.50) {
  
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

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------


root_path <- "/media/antonio/Work/MMDB"

# Communities
taxo_comm <- fread(paste0(root_path, "/Diversity/02_NMDS/NMDS_scores_comm_taxo_total.csv"))
phylo_comm <- fread(paste0(root_path, "/Diversity/02_NMDS/NMDS_scores_comm_phylo_total.csv"))
func_comm <- fread(paste0(root_path, "/Diversity/02_NMDS/NMDS_scores_comm_func_total.csv"))

# Get domains
domain <- fread(paste0(root_path, "/FIA/00_FIA/", "plot_domain.csv"))
mat <- match.integer64(taxo_comm$PLT_CN, domain$PLT_CN)
domain <- domain[mat, ]
all.equal(domain$PLT_CN, taxo_comm$PLT_CN)

# Make a vector file
taxo_sf <- cbind(domain[, c("LON", "LAT")], taxo_comm)
phylo_sf <- cbind(domain[, c("LON", "LAT")], phylo_comm)
func_sf <- cbind(domain[, c("LON", "LAT")], func_comm)

taxo_sf <- st_as_sf(taxo_sf, coords = c("LON","LAT"), crs = 4269)
phylo_sf <- st_as_sf(phylo_sf, coords = c("LON","LAT"), crs = 4269)
func_sf <- st_as_sf(func_sf, coords = c("LON","LAT"), crs = 4269)

# Read polygon
states <- st_read(paste0(root_path, "/GIS/States.gpkg"))
neon <- st_read(paste0(root_path, "/GIS/NEON.gpkg"))
neon_aoi <- st_read(paste0(root_path, "/GIS/NEON_aoi.gpkg"))

# Make plot
taxo_col <- RGB(taxo_comm$Axis1, 
           taxo_comm$Axis2, 
           taxo_comm$Axis3,
           min_range = NULL, 
           max_range = NULL)

phylo_col <- RGB(phylo_comm$Axis1, 
                 phylo_comm$Axis2, 
                 phylo_comm$Axis3,
                 min_range = NULL, 
                 max_range = NULL)

func_col <- RGB(func_comm$Axis1, 
                func_comm$Axis2, 
                func_comm$Axis3,
                min_range = NULL, 
                max_range = NULL)

th <- theme(panel.grid.major = element_line(colour = gray(0.5), 
                                            linetype = "dotted", 
                                            linewidth = 0.1), 
            legend.background = element_rect(fill= "transparent"),
            panel.border = element_rect(fill = "transparent"),
            panel.background = element_rect(fill= "transparent"),
            rect = element_rect(fill = "transparent"),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"))

th_trans <- dark_theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        #legend.position.inside = c(0.90, 0.90), 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth=0),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"))

th_selected <- th
th_selected <- th_trans

taxo <- ggplot() +
  geom_sf(data = states, fill = "grey80", colour = "grey50", linewidth = 0.1) + 
  #geom_sf(data = neon, fill = "NA", colour = "grey25", linewidth = 0.1) +
  geom_sf(data = neon_aoi, fill = "black", colour = "white", linewidth = 0.1, alpha = 1.0) +
  geom_sf(data = taxo_sf, size = 0.1, color = taxo_col, shape = 16) +
  geom_sf(data = states, fill = "NA", colour = "grey50", linewidth = 0.1) +
  geom_sf(data = neon_aoi, fill = "NA", colour = "white", linewidth = 0.1) +  
  coord_sf(crs = '+proj=lcc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs',
           datum = sf::st_crs(4326),
           xlim = c(-290000, 2300000)) +
  scale_x_continuous(breaks = c(-100, -90, -80, -70)) +
  scale_y_continuous(breaks = c(50, 40, 30, 20)) +
  xlab(" ") + ylab("Latitude (°)") +
  theme_bw() + th_selected

phylo <- ggplot() +
  geom_sf(data = states, fill = "grey80", colour = "grey50", linewidth = 0.1) + 
  #geom_sf(data = neon, fill = "NA", colour = "grey25", linewidth = 0.1) +
  geom_sf(data = neon_aoi, fill = "black", colour = "white", linewidth = 0.1, alpha = 1.0) +
  geom_sf(data = phylo_sf, size = 0.1, color = phylo_col, shape = 16) +
  geom_sf(data = states, fill = "NA", colour = "grey50", linewidth = 0.1) +
  geom_sf(data = neon_aoi, fill = "NA", colour = "white", linewidth = 0.1) +  
  coord_sf(crs = '+proj=lcc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs',
           datum = sf::st_crs(4326),
           xlim = c(-290000, 2300000)) +
  scale_x_continuous(breaks = c(-100, -90, -80, -70)) +
  scale_y_continuous(breaks = c(50, 40, 30, 20)) +
  xlab("Longitude (°)") + ylab(" ") +
  theme_bw() + th_selected

func <- ggplot() +
  geom_sf(data = states, fill = "grey80", colour = "grey50", linewidth = 0.1) + 
  #geom_sf(data = neon, fill = "NA", colour = "grey25", linewidth = 0.1) +
  geom_sf(data = neon_aoi, fill = "black", colour = "white", linewidth = 0.1, alpha = 1.0) +
  geom_sf(data = func_sf, size = 0.1, color = func_col, shape = 16) +
  geom_sf(data = states, fill = "NA", colour = "grey50", linewidth = 0.1) +
  geom_sf(data = neon_aoi, fill = "NA", colour = "white", linewidth = 0.1) +  
  coord_sf(crs = '+proj=lcc +lat_0=39 +lon_0=-96 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs',
           datum = sf::st_crs(4326),
           xlim = c(-290000, 2300000)) +
  scale_x_continuous(breaks = c(-100, -90, -80, -70)) +
  scale_y_continuous(breaks = c(50, 40, 30, 20)) +
  xlab(" ") + ylab(" ") +
  theme_bw() + th_selected


plot <- ggarrange(taxo, phylo, func, ncol = 3)

# Black
jpeg(paste0(root_path, "/figures/Manuscript_black/", "NMDS_spatial.jpeg"),
    width = 210, 
    height = 90, 
    units = "mm", 
    pointsize = 12,
    bg = "white",  
    res = 600)

plot

dev.off()

# Transparent
png(paste0(root_path, "/figures/Manuscript_white/", "NMDS_spatial_trans.png"),
    width = 210, 
    height = 90, 
    units = "mm", 
    pointsize = 12,
    bg = "transparent",  
    res = 600)

plot

dev.off()
