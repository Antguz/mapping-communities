#'------------------------------------------------------------------------------
#' @title Mantel test on distance matrices
#'------------------------------------------------------------------------------

#' @description It tune, model, and predict the ordinations based on spectra.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(parallelDist)
library(vegan)
library(bit64)
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)
options(scipen = 100000)

#'------------------------------------------------------------------------------
#' Source code
sourceCpp("src/sam.cpp")

#'------------------------------------------------------------------------------
#' Working directory

#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
sensor <- "DESIS"

#'------------------------------------------------------------------------------
#' Read information

spectra <- fread(paste0(root_path, "/04_match_spectra-stand/", sensor, "/", sensor,"_spectra_SW.csv"))

#'------------------------------------------------------------------------------
#' @Analysis  

# 400 days difference
if(sensor == "DESIS") {
  spectra <- spectra[abs(dtime) <= 455, ]
  
  # Bands from 450 to 962.53
  spectra <- spectra[, .SD, .SDcols = c(2, 39:211)]

} else if(sensor == "EMIT") {
  
  # Bands from 447.74 to 1265.17, from 1488.07 to 1710.97, from 2015.60 to 2424.25
  spectra <- spectra[, .SD, .SDcols = c(2, 
                                        32:142, # from 447.74 to 1265.17
                                        172:202, #from 1488.07 to 1710.97
                                        243:298)] #from 2015.60 to 2424.25
  
}

if(sensor == "DESIS") {
  outliers <- fread(paste0(root_path, "/06_model/SW_beta_PLS/DESIS_outliers.csv"))
  outliers <- outliers[outliers$remove != TRUE, ]
  outliers <- outliers[,1]
  spectra <- merge(outliers, spectra, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
}

plot_cn <- spectra$PLT_CN
spectra <- spectra[, -1]
spectra <- as.data.frame(spectra)
rownames(spectra) <- plot_cn

# Get distance
spectral_distance <-  sam(as.matrix(spectra), threads = 28) 
colnames(spectral_distance) <- as.character(plot_cn)
rownames(spectral_distance) <- as.character(plot_cn)

#'------------------------------------------------------------------------------
#' Function
mantel_test <- function(spectral_distance, beta) {

  #'----------------------------------------------------------------------------
  # Get beta data
  Btotal <- beta$Btotal/10000
  #Brepl <- beta$Brepl/10000
  #Brich <- beta$Brich/10000
  
  # Get plot cn
  PLOT_CN_spectra <- colnames(spectral_distance)
  PLOT_CN_spectra <- as.integer64(PLOT_CN_spectra)
  
  PLOT_CN_beta <- colnames(Btotal)
  PLOT_CN_beta <- as.integer64(PLOT_CN_beta)
  
  #'----------------------------------------------------------------------------
  #' Match spectral distance and beta diversity
  match <- match.integer64(PLOT_CN_spectra, PLOT_CN_beta)
  
  # Re-select beta plots
  Btotal <- Btotal[match, match]
  #Brepl <- Brepl[match, match]
  #Brich <- Brich[match, match]
  
  #'----------------------------------------------------------------------------
  #' Compute mantel test
  mantel_Btotal <- mantel(xdis = Btotal, 
                          ydis = spectral_distance, 
                          method = "pearson", 
                          permutations = 999, 
                          strata = NULL,
                          na.rm = TRUE)
  
  #mantel_Brepl <- mantel(xdis = Brepl[1:100, 1:100], 
  #                       ydis = spectral_distance[1:100, 1:100], 
  #                       method = "pearson", 
  #                       permutations = 999, 
  #                       strata = NULL,
  #                       na.rm = TRUE)
  
  #mantel_Brich <- mantel(xdis = Brich[1:100, 1:100], 
  #                       ydis = spectral_distance[1:100, 1:100], 
  #                       method = "pearson", 
  #                       permutations = 999, 
  #                       strata = NULL,
  #                       na.rm = TRUE)
  
  #'----------------------------------------------------------------------------
  #' Melt frames
  
  # As data.table
  spectral_distance <- as.data.table(spectral_distance)
  spectral_distance$PLOT_CN <- colnames(spectral_distance)
  Btotal <- as.data.table(Btotal)
  Btotal$PLOT_CN <- colnames(Btotal)
  #Brepl <- as.data.table(Brepl)
  #Brepl$PLOT_CN <- colnames(Brepl)
  #Brich <- as.data.table(Brich)
  #Brich$PLOT_CN <- colnames(Brich)
  
  # Melt frames
  spectral_distance <- melt(spectral_distance, 
                            id.vars = "PLOT_CN",
                            variable.name = "PLOT_CN_col",
                            value.name = "spectral")
  spectral_distance$PLOT_CN <- as.integer64(spectral_distance$PLOT_CN)
  spectral_distance$PLOT_CN_col <- as.integer64(as.character(spectral_distance$PLOT_CN_col))
  
  Btotal <- melt(Btotal, 
                 id.vars = "PLOT_CN",
                 variable.name = "PLOT_CN_col",
                 value.name = "Btotal")
  Btotal$PLOT_CN <- as.integer64(Btotal$PLOT_CN)
  Btotal$PLOT_CN_col <- as.integer64(as.character(Btotal$PLOT_CN_col))
  
  #Brepl <- melt(Brepl, 
  #              id.vars = "PLOT_CN",
  #              variable.name = "PLOT_CN_col",
  #              value.name = "Brepl")
  #Brepl$PLOT_CN <- as.integer64(Brepl$PLOT_CN)
  #Brepl$PLOT_CN_col <- as.integer64(as.character(Brepl$PLOT_CN_col))
  
  #Brich <- melt(Brich, 
  #              id.vars = "PLOT_CN",
  #              variable.name = "PLOT_CN_col",
  #              value.name = "Brich")
  #Brich$PLOT_CN <- as.integer64(Brich$PLOT_CN)
  #Brich$PLOT_CN_col <- as.integer64(as.character(Brich$PLOT_CN_col))
  
  frame <- merge(spectral_distance, 
                 Btotal, 
                 by = c("PLOT_CN", "PLOT_CN_col"), 
                 all.x = TRUE,
                 all.y = TRUE)
  
  #frame <- merge(frame, 
  #               Brepl, 
  #               by = c("PLOT_CN", "PLOT_CN_col"), 
  #               all.x = TRUE,
  #               all.y = TRUE)
  
  #frame <- merge(frame, 
  #               Brich, 
  #               by = c("PLOT_CN", "PLOT_CN_col"), 
  #               all.x = TRUE,
  #               all.y = TRUE)
  
  return(list(mantel = mantel_Btotal,
              frame = frame))
  
}

#'------------------------------------------------------------------------------
#' Mantel test

# Taxonomic
type <- "taxo"
beta <- readRDS(paste0(root_path, "/03_diversity/beta_", type, "_VOL.rds"))
test <- mantel_test(spectral_distance, beta)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_mantel.rds"))
rm(list = c("beta", "test"))
gc()

# Phylogenetic
type <- "phylo-log"
beta <- readRDS(paste0(root_path, "/03_diversity/beta_", type, "_VOL.rds"))
test <- mantel_test(spectral_distance, beta)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_mantel.rds"))
rm(list = c("beta", "test"))
gc()

# Functional
type <- "func"
beta <- readRDS(paste0(root_path, "/03_diversity/beta_", type, "_VOL.rds"))
test <- mantel_test(spectral_distance, beta)
saveRDS(test, paste0(root_path, "/06_mantel/", sensor, "_", type, "_mantel.rds"))
rm(list = c("beta", "test"))
gc()
