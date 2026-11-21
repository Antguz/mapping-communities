#'------------------------------------------------------------------------------
#' @title Extract environmental data from WorldClim and HYDRO1k
#'------------------------------------------------------------------------------

#' @description It extracts climatic (WorldClim: MAT, TS, TAR, MAP, PS) and
#' topographic (HYDRO1k: DEM, SL, CTI) data based on plot locations.
#'
#' @return A .csv file (fia-env.csv) with the PLT_CN and columns of
#' environmental information.

#'------------------------------------------------------------------------------
#' @Library

library(terra)
library(data.table)
library(bit64)

#'------------------------------------------------------------------------------
#' @Paths-and-files

# Root path
root_path <- "/media/antonio/Work/MMDB"

# Read files
fia_plots <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))
fia_reshape <- readRDS(paste0(root_path, "/FIA/01_FIA-reshape/fia_reshape.rds"))

# Select clean plots with coordinates
PLT_CN <- unique(fia_reshape$PLT_CN)
sel <- match.integer64(PLT_CN, fia_plots$PLT_CN)
fia_plots <- fia_plots[sel,]
all.equal(PLT_CN, fia_plots$PLT_CN)

#'------------------------------------------------------------------------------
#' @Geo-information

# Create vector to extract pixels
FIA_vector <- vect(x = fia_plots[, c("PLT_CN", "LAT", "LON")], 
                   geom=c("LON", "LAT"), 
                   crs = "EPSG:4269", 
                   keepgeom = FALSE)

# Read WorlClim information
MAT <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_1.tif"))
TS <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_4.tif"))
TAR <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_7.tif"))
MAP <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_12.tif"))
PS <- rast(paste0(root_path, "/env/wc2.1_30s_bio/wc2.1_30s_bio_15.tif"))
WorlClim <- c(MAT, TS, TAR, MAP, PS)
names(WorlClim) <- c("MAT", "TS", "TAR", "MAP", "PS")

# Read DEM information
DEM <- rast(paste0(root_path, "/env/gt30h1kna/na_dem.bil"))
DEM[DEM == -9999] <- NA
names(DEM) <- "DEM"

# Read HYDRO1k information
slope <- rast(paste0(root_path, "/env/gt30h1kna/na_slope.bil"))
slope[slope == -9999] <- NA
CTI <- rast(paste0(root_path, "/env/gt30h1kna/na_cti.bil"))
CTI[CTI == -9999] <- NA
HYDRO1k <- c(slope, CTI)
names(HYDRO1k) <- c("SL", "CTI")
crs(HYDRO1k) <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"

#'------------------------------------------------------------------------------
#' @Pixel-extraction

# Get WorlCLim
if(crs(FIA_vector) != crs(WorlClim)) {
  points <- project(FIA_vector, WorlClim)
}

WorlClim_extract <- extract(WorlClim, points, cells=FALSE, xy=FALSE)

# Get DEM
if(crs(FIA_vector) != crs(DEM)) {
  points <- project(FIA_vector, DEM)
}

DEM_extract <- extract(DEM, points, cells=FALSE, xy=FALSE)

# Get HYDRO1k
if(crs(FIA_vector) != crs(HYDRO1k)) {
  points <- project(FIA_vector, HYDRO1k)
}

HYDRO1k_extract <- extract(HYDRO1k, points, cells=FALSE, xy=FALSE)

#'------------------------------------------------------------------------------
#' @Merge-information-and-export

# Merge
env <- cbind(fia_plots[, c("PLT_CN")], 
             WorlClim_extract[,-1], 
             DEM = DEM_extract[,-1], 
             HYDRO1k_extract[, -1])

# Clean
env[DEM == 55537, DEM := NA]
env[SL == 55537, SL := NA]
env[CTI == 55537, CTI := NA]

# Export
fwrite(env, paste0(root_path, "/FIA/03_FIA-env/fia-env.csv"))

# DO A MANUAL INSPECTION 


