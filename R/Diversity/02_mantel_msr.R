#'------------------------------------------------------------------------------
#' @title Spatially-constrained Mantel test on distance matrices (MSR)
#'------------------------------------------------------------------------------

#' @description Tests for dissimilarity between matrices of spectral
#' dissimilarity (Spectral Angle) and beta diversity accounting for spatial
#' autocorrelation with Moran Spectral Randomization (MSR; Crabot et al. 2019;
#' Wagner & Dray 2015). Replicates of the spectral dissimilarity preserve its
#' spatial structure on distance-based MEMs (dbMEM; Dray et al. 2006). The
#' randomized matrix must be Euclidean, so the test uses sqrt(SAM) against
#' sqrt(Btotal), sqrt(Brepl), and sqrt(Brich). The beta types of a diversity
#' are tested in one pass (same replicates). Results do not overwrite existing
#' files.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(bit64)
library(spdep)
library(adespatial)
library(sf)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/sam.cpp")
source("R/auxiliary/mantel_msr.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param sensor
#' @param spectra_disimilarity
#' @param location data.table with PLT_CN, LAT, LON
#' @param diversity
#' @param beta_path
#' @param beta_type One or more of "total", "replacement", "richness"
#' @param mem_file File to store (or read) the MEMs
#' @param output_path
#' @param nrepet
#' @param threads

#'------------------------------------------------------------------------------
#' Function to get the MEMs

#' Distance-based MEMs (dbMEM; Dray et al. 2006) with the default truncation
#' threshold (longest edge of the minimum spanning tree). Coordinates are
#' projected to USA Contiguous Albers Equal Area (EPSG:5070) in km.

get_mem <- function(location, PLOT_CN_spectra, mem_file) {

  if(file.exists(mem_file)) {
    mem <- readRDS(mem_file)
    stopifnot(identical(rownames(mem), as.character(PLOT_CN_spectra)))
    return(mem)
  }

  # Match locations and spectral dissimilarity
  match <- match(PLOT_CN_spectra, location$PLT_CN)
  stopifnot(!anyNA(match))
  location <- location[match, ]

  # Projected coordinates
  xy <- st_as_sf(location, coords = c("LON", "LAT"), crs = 4326)
  xy <- st_coordinates(st_transform(xy, 5070))/1000
  rownames(xy) <- as.character(location$PLT_CN)

  # dbMEM
  mem <- dbmem(xy, MEM.autocor = "all", silent = FALSE)

  saveRDS(mem, mem_file)
  return(mem)

}

#'------------------------------------------------------------------------------
#' Function

mantel_msr_test <- function(sensor,
                            spectra_disimilarity,
                            location,
                            diversity,
                            beta_path,
                            beta_type = c("total", "replacement", "richness"),
                            mem_file,
                            output_path,
                            nrepet = 999,
                            threads) {

  out_file <- paste0(output_path, "/", sensor, "_", diversity, "_mantel_msr.rds")

  if(file.exists(out_file)) {
    message("Exists, skipped: ", out_file)
    return(invisible(NULL))
  }

  #'----------------------------------------------------------------------------
  # Get info from spectral dissimilarity

  PLOT_CN_spectra <- colnames(spectra_disimilarity)
  PLOT_CN_spectra <- as.integer64(PLOT_CN_spectra)

  #'----------------------------------------------------------------------------
  # Get beta data

  beta <- readRDS(beta_path)
  PLOT_CN_beta <- colnames(beta$Btotal)
  PLOT_CN_beta <- as.integer64(PLOT_CN_beta)

  # Match spectral distance and beta diversity
  match <- match(PLOT_CN_spectra, PLOT_CN_beta)
  stopifnot(!anyNA(match))

  # sqrt of the beta types (Euclidean)
  component <- c(total = "Btotal", replacement = "Brepl", richness = "Brich")
  Beta <- lapply(beta_type, function(x) {
    sqrt(beta[[component[x]]][match, match]/10000)
  })
  names(Beta) <- beta_type

  rm("beta")
  gc()

  #'----------------------------------------------------------------------------
  # Get MEMs

  mem <- get_mem(location, PLOT_CN_spectra, mem_file)

  #'----------------------------------------------------------------------------
  # Mantel test with Moran Spectral Randomization

  set.seed(123)
  mantel_MSR <- mantel_msr(m1 = sqrt(spectra_disimilarity),
                           m2 = Beta,
                           mem = mem,
                           nrepet = nrepet,
                           threads = threads,
                           progress = TRUE)

  frame <- data.table(sensor = sensor,
                      diversity = diversity,
                      beta_type = beta_type,
                      n = nrow(spectra_disimilarity),
                      mantel = sapply(mantel_MSR, function(x) x$obs),
                      mantel_significance = sapply(mantel_MSR, function(x) x$pvalue),
                      mantel_msr = lapply(mantel_MSR, function(x) x$sim))

  #'----------------------------------------------------------------------------
  # Export

  saveRDS(frame, out_file)

  return(frame[, c("sensor", "diversity", "beta_type", "n", "mantel", "mantel_significance")])

}

#'------------------------------------------------------------------------------
#' Mantel test

#---------------------------------------------------------
#' @DESIS

#root_path <- "/media/guzman/Work/MMDB"
root_path <- "/n/netscratch/jcavender_lab/Everyone/guzman/MMDB"
sensor <- "DESIS"

# # Spectra from plots used
# data <- readRDS(paste0(root_path, "/Models/01_model_beta/DESIS/DESIS_taxo_data.rds"))
# data <- data[, c("scene", "PLT_CN")]
# 
# # Spectra avaiable
# spectra <- fread(paste0(root_path, "/Pixel_extraction/03_match-spectra/DESIS/DESIS_spectra_SG.csv"))
# spectra <- merge(data, spectra, by = c("scene", "PLT_CN"), all.x = TRUE, all.y = FALSE)
# 
# # Select bands for sam
# spectra <- spectra[, .SD, .SDcols = c(2, 53:228)] #"447.76" to "895.97"
# spectra_disimilarity <- sam(as.matrix(spectra[,-1]), threads = 26, progress = TRUE)
# colnames(spectra_disimilarity) <- as.character(spectra$PLT_CN)
# rownames(spectra_disimilarity) <- as.character(spectra$PLT_CN)
# saveRDS(spectra_disimilarity, paste0(root_path, "/Models/01_model_beta/DESIS/spectra_disimilarity.rds"))

# Read spectra disimilarity
spectra_disimilarity <- readRDS(paste0(root_path, "/Models/01_model_beta/DESIS/spectra_disimilarity.rds"))

# Plot locations
location <- fread(paste0(root_path, "/Models/01_model_beta/DESIS/DESIS_taxo_location.csv"))

# Select output
output_path <- paste0(root_path, "/Diversity/03_Mantel_MSR")
dir.create(output_path, showWarnings = FALSE)
mem_file <- paste0(output_path, "/DESIS_dbMEM.rds")
threads <- 16

#---------------------------------------------------------
# Taxonomic
diversity <- "taxo"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger.rds")

mantel_msr_test(sensor = "DESIS",
                spectra_disimilarity = spectra_disimilarity,
                location = location,
                diversity = "taxo",
                beta_path = beta_path,
                beta_type = c("total", "replacement", "richness"),
                mem_file = mem_file,
                output_path = output_path,
                threads = threads)
gc()

#---------------------------------------------------------
# Phylogenetic
diversity <- "phylo"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger_norm.rds")

mantel_msr_test(sensor = "DESIS",
                spectra_disimilarity = spectra_disimilarity,
                location = location,
                diversity = "phylo",
                beta_path = beta_path,
                beta_type = c("total", "replacement", "richness"),
                mem_file = mem_file,
                output_path = output_path,
                threads = threads)
gc()

#---------------------------------------------------------
# Functional
diversity <- "func"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger_norm.rds")

mantel_msr_test(sensor = "DESIS",
                spectra_disimilarity = spectra_disimilarity,
                location = location,
                diversity = "func",
                beta_path = beta_path,
                beta_type = c("total", "replacement", "richness"),
                mem_file = mem_file,
                output_path = output_path,
                threads = threads)
gc()

#---------------------------------------------------------
#' @EMIT

#root_path <- "/media/guzman/Work/MMDB"
root_path <- "/n/netscratch/jcavender_lab/Everyone/guzman/MMDB"
sensor <- "EMIT"

# # Spectra from plots used
# data <- readRDS(paste0(root_path, "/Models/01_model_beta/EMIT/EMIT_taxo_data.rds"))
# data <- data[, c("scene", "PLT_CN")]
# 
# # Spectra available
# spectra <- fread(paste0(root_path, "/Pixel_extraction/03_match-spectra/EMIT/EMIT_spectra_SG.csv"))
# spectra <- merge(data, spectra, by = c("scene", "PLT_CN"), all.x = TRUE, all.y = FALSE)
# 
# # Select bands for sam
# spectra <- spectra[, .SD, .SDcols = c(2,
#                                       42:152,
#                                       184:211,
#                                       252:308)] #447.74 - 2426.44
# spectra_disimilarity <- sam(as.matrix(spectra[, -1]), threads = 26, progress = TRUE)
# colnames(spectra_disimilarity) <- as.character(spectra$PLT_CN)
# rownames(spectra_disimilarity) <- as.character(spectra$PLT_CN)
# saveRDS(spectra_disimilarity, paste0(root_path, "/Models/01_model_beta/EMIT/spectra_disimilarity.rds"))

# Read the spectra dissimilarity
spectra_disimilarity <- readRDS(paste0(root_path, "/Models/01_model_beta/EMIT/spectra_disimilarity.rds"))

# Plot locations
location <- fread(paste0(root_path, "/Models/01_model_beta/EMIT/EMIT_taxo_location.csv"))

# Select output
output_path <- paste0(root_path, "/Diversity/03_Mantel_MSR")
dir.create(output_path, showWarnings = FALSE)
mem_file <- paste0(output_path, "/EMIT_dbMEM.rds")
threads <- 16

#---------------------------------------------------------
# Taxonomic
diversity <- "taxo"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger.rds")

mantel_msr_test(sensor = "EMIT",
                spectra_disimilarity = spectra_disimilarity,
                location = location,
                diversity = "taxo",
                beta_path = beta_path,
                beta_type = c("total", "replacement", "richness"),
                mem_file = mem_file,
                output_path = output_path,
                threads = threads)
gc()

#---------------------------------------------------------
# Phylogenetic
diversity <- "phylo"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger_norm.rds")

mantel_msr_test(sensor = "EMIT",
                spectra_disimilarity = spectra_disimilarity,
                location = location,
                diversity = "phylo",
                beta_path = beta_path,
                beta_type = c("total", "replacement", "richness"),
                mem_file = mem_file,
                output_path = output_path,
                threads = threads)
gc()

#---------------------------------------------------------
# Functional
diversity <- "func"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger_norm.rds")

mantel_msr_test(sensor = "EMIT",
                spectra_disimilarity = spectra_disimilarity,
                location = location,
                diversity = "func",
                beta_path = beta_path,
                beta_type = c("total", "replacement", "richness"),
                mem_file = mem_file,
                output_path = output_path,
                threads = threads)
gc()

#---------------------------------------------------------
#' @EMIT-short

#root_path <- "/media/guzman/Work/MMDB"
root_path <- "/n/netscratch/jcavender_lab/Everyone/guzman/MMDB"
sensor <- "EMITs"

# # Spectra from plots used
# data <- readRDS(paste0(root_path, "/Models/01_model_beta/EMIT/EMIT_taxo_data.rds"))
# data <- data[, c("scene", "PLT_CN")]
# 
# # Spectra available
# spectra <- fread(paste0(root_path, "/Pixel_extraction/03_match-spectra/EMIT/EMIT_spectra_SG.csv"))
# spectra <- merge(data, spectra, by = c("scene", "PLT_CN"), all.x = TRUE, all.y = FALSE)
# 
# # Select bands for sam
# spectra <- spectra[, .SD, .SDcols = c(2,
#                                      42:102)] # 447.74 - 894.9
# spectra_disimilarity <- sam(as.matrix(spectra[, -1]), threads = 26, progress = TRUE)
# colnames(spectra_disimilarity) <- as.character(spectra$PLT_CN)
# rownames(spectra_disimilarity) <- as.character(spectra$PLT_CN)
# saveRDS(spectra_disimilarity, paste0(root_path, "/Models/01_model_beta/EMIT/spectra_disimilarity_short.rds"))

# Read the spectra dissimilarity
spectra_disimilarity <- readRDS(paste0(root_path, "/Models/01_model_beta/EMIT/spectra_disimilarity_short.rds"))

# Plot locations
location <- fread(paste0(root_path, "/Models/01_model_beta/EMIT/EMIT_taxo_location.csv"))

# Select output
output_path <- paste0(root_path, "/Diversity/03_Mantel_MSR")
dir.create(output_path, showWarnings = FALSE)
# Same plots as EMIT: uses the same MEMs
mem_file <- paste0(output_path, "/EMIT_dbMEM.rds")
threads <- 16

#---------------------------------------------------------
# Taxonomic
diversity <- "taxo"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger.rds")

mantel_msr_test(sensor = "EMITs",
                spectra_disimilarity = spectra_disimilarity,
                location = location,
                diversity = "taxo",
                beta_path = beta_path,
                beta_type = c("total", "replacement", "richness"),
                mem_file = mem_file,
                output_path = output_path,
                threads = threads)
gc()

#---------------------------------------------------------
# Phylogenetic
diversity <- "phylo"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger_norm.rds")

mantel_msr_test(sensor = "EMITs",
                spectra_disimilarity = spectra_disimilarity,
                location = location,
                diversity = "phylo",
                beta_path = beta_path,
                beta_type = c("total", "replacement", "richness"),
                mem_file = mem_file,
                output_path = output_path,
                threads = threads)
gc()

#---------------------------------------------------------
# Functional
diversity <- "func"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger_norm.rds")

mantel_msr_test(sensor = "EMITs",
                spectra_disimilarity = spectra_disimilarity,
                location = location,
                diversity = "func",
                beta_path = beta_path,
                beta_type = c("total", "replacement", "richness"),
                mem_file = mem_file,
                output_path = output_path,
                threads = threads)
gc()
