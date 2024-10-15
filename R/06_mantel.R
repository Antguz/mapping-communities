#'------------------------------------------------------------------------------
#' @title Mantel test on distance matrices
#'------------------------------------------------------------------------------

#' @description Tests for dissimilarity between matrices of beta diversity and
#' spectral dissimilarity (Spectral Angle)

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
library(bit64)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
options(scipen = 100000)

#'------------------------------------------------------------------------------
#' Source code

setwd("/panfs/jay/groups/17/cavender/guzman/MMDB")
sourceCpp("src/sam.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' @param sensor
#' @param spectral_path
#' @param diversity
#' @param beta_path
#' @param beta_type
#' @param output_path
#' @param threads

#'------------------------------------------------------------------------------
#' Function
mantel_test <- function(sensor,
                        spectral_path,
                        diversity,
                        beta_path, 
                        beta_type,
                        output_path, 
                        threads = 6) {
  
  #'----------------------------------------------------------------------------
  # Get spectral data
  
  spectra <- fread(spectral_path)
  
  if(sensor == "DESIS") {
    
    spectra <- spectra[, .SD, .SDcols = c(2, 39:192)] # Bands from 450 to 901
    
  } else if(sensor == "EMIT") {
    
    spectra <- spectra[, .SD, .SDcols = c(2, 
                                          32:142, # from 447.74 to 1265.17
                                          174:204, #from 1488.07 to 1720
                                          243:298)] #from 2015.60 to 2424.25
    
  }
  
  # Clean spectra frame
  plot_cn <- spectra$PLT_CN
  spectra <- spectra[, -1]
  spectra <- as.data.frame(spectra)
  rownames(spectra) <- plot_cn
  
  # Get distance
  spectral_distance <-  sam(as.matrix(spectra), threads = 28) 
  colnames(spectral_distance) <- as.character(plot_cn)
  rownames(spectral_distance) <- as.character(plot_cn)
  
  # Remove and clean
  rm(list = c("spectra", "plot_cn")) 
  gc()
  
  #'----------------------------------------------------------------------------
  # Get beta data
  beta <- readRDS(beta_path)
  
  if(diversity == "total") {
    Beta <- beta$Btotal/10000
  } else if(diversity == "replacement") {
    Beta <- beta$Brepl/10000
  } else if(diversity == "richness") {
    Beta <- beta$Brich/10000
  }
  
  rm("beta")
  gc()
  
  # Get plot cn
  PLOT_CN_spectra <- colnames(spectral_distance)
  PLOT_CN_spectra <- as.integer64(PLOT_CN_spectra)
  
  PLOT_CN_beta <- colnames(Beta)
  PLOT_CN_beta <- as.integer64(PLOT_CN_beta)
  
  #'----------------------------------------------------------------------------
  #' Match spectral distance and beta diversity
  match <- match.integer64(PLOT_CN_spectra, PLOT_CN_beta)
  
  # Re-select beta plots
  Beta <- Beta[match, match]
  gc()
  
  #'----------------------------------------------------------------------------
  #' Compute mantel test
  mantel_Beta <- mantel(xdis = Beta, 
                        ydis = spectral_distance, 
                        method = "pearson", 
                        permutations = 999, 
                        strata = NULL,
                        na.rm = TRUE,
                        parallel = threads)
  
  #'----------------------------------------------------------------------------
  #' Melt frames
  
  # As data.table
  spectral_distance <- as.data.table(spectral_distance)
  spectral_distance$PLOT_CN <- colnames(spectral_distance)
  Beta <- as.data.table(Beta)
  Beta$PLOT_CN <- colnames(Beta)
  
  # Melt frames
  spectral_distance <- melt(spectral_distance, 
                            id.vars = "PLOT_CN",
                            variable.name = "PLOT_CN_col",
                            value.name = "spectral")
  spectral_distance$PLOT_CN <- as.integer64(spectral_distance$PLOT_CN)
  spectral_distance$PLOT_CN_col <- as.integer64(as.character(spectral_distance$PLOT_CN_col))
  
  Beta <- melt(Beta, 
               id.vars = "PLOT_CN",
               variable.name = "PLOT_CN_col",
               value.name = "Beta")
  Beta$PLOT_CN <- as.integer64(Beta$PLOT_CN)
  Beta$PLOT_CN_col <- as.integer64(as.character(Beta$PLOT_CN_col))
  
  frame <- merge(spectral_distance, 
                 Beta, 
                 by = c("PLOT_CN", "PLOT_CN_col"), 
                 all.x = TRUE,
                 all.y = TRUE)
  
  return_list <- list(mantel = mantel_Beta,
                      frame = frame)
  
  saveRDS(return, paste0(output_path, "/", 
                         sensor, "_", type, "_", diversity, "_mantel.rds"))
  
}

#'------------------------------------------------------------------------------
#' Mantel test

sensor,
                        spectral_path,
                        diversity,
                        beta_path, 
                        beta_type,
                        output_path, 
                        threads = 6

#' @DESIS

sensor <- "DESIS"
spectral_path <- paste0(root_path, 
                        "/04_match_spectra-stand/", 
                        sensor, "/", sensor, "_spectra_SW.csv")

#---------------------------------------------------------
# Taxonomic
type <- "taxo"
beta <- paste0(root_path, "/03_diversity/beta_", type, "_VOL.rds")

diversity <- "total"
test <- mantel_test(spectral_distance, beta, diversity)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
gc()

diversity <- "replacement"
test <- mantel_test(spectral_distance, beta, diversity)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
gc()

diversity <- "richness"
test <- mantel_test(spectral_distance, beta, diversity)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
gc()

# Phylogenetic
type <- "phylo"
beta <- readRDS(paste0(root_path, "/03_diversity/beta_", type, "_VOL.rds"))

diversity <- "total"
test <- mantel_test(spectral_distance, beta, diversity)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
gc()

diversity <- "replacement"
test <- mantel_test(spectral_distance, beta, diversity)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
gc()

diversity <- "richness"
test <- mantel_test(spectral_distance, beta, diversity)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
gc()

# Functional
type <- "func"
beta <- readRDS(paste0(root_path, "/03_diversity/beta_", type, "_VOL.rds"))

diversity <- "total"
test <- mantel_test(spectral_distance, beta, diversity)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
gc()

diversity <- "replacement"
test <- mantel_test(spectral_distance, beta, diversity)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
gc()

diversity <- "richness"
test <- mantel_test(spectral_distance, beta, diversity)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
gc()

