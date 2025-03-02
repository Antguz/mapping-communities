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

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/sam.cpp")
sourceCpp("src/mantel_rcpp.cpp")

#'------------------------------------------------------------------------------
#' Arguments
#' @param sensor
#' @param spectra_disimilarity
#' @param diversity
#' @param beta_path
#' @param beta_type
#' @param output_path
#' @param threads

#'------------------------------------------------------------------------------
#' Function
mantel_test <- function(sensor,
                        spectra_disimilarity,
                        beta_path, 
                        diversity,
                        beta_type,
                        output_path, 
                        threads) {
  
  #'----------------------------------------------------------------------------
  # Get info from spectral dissimilarity
  
  PLOT_CN_spectra <- colnames(spectra_disimilarity)
  PLOT_CN_spectra <- as.integer64(PLOT_CN_spectra)
  
  #'----------------------------------------------------------------------------
  # Get beta data
  beta <- readRDS(beta_path)
  PLOT_CN_beta <- colnames(beta$Btotal)
  PLOT_CN_beta <- as.integer64(PLOT_CN_beta)
  
  if(beta_type == "total") {
    Beta <- beta$Btotal/10000
  } else if(beta_type == "replacement") {
    Beta <- (beta$Brepl/10000) # / (beta$Btotal/10000)
    Beta[is.nan(Beta)] <- 0
  } else if(beta_type == "richness") {
    Beta <- (beta$Brich/10000) # / (beta$Btotal/10000)
    Beta[is.nan(Beta)] <- 0
  } else if(beta_type == "replacement_total") {
    Beta <- (beta$Brepl/beta$Btotal)
    Beta[is.nan(Beta)] <- 0
  } else if(beta_type == "richness_total") {
    Beta <- (beta$Brich/beta$Btotal)
    Beta[is.nan(Beta)] <- 0
  }
  
  rm("beta")
  gc()
  
  # Get plot cn
  colnames(Beta) <- as.character(PLOT_CN_beta)
  rownames(Beta) <- as.character(PLOT_CN_beta)
  
  #'----------------------------------------------------------------------------
  #' Match spectral distance and beta diversity
  match <- match.integer64(PLOT_CN_spectra, PLOT_CN_beta)
  
  # Re-select beta plots
  Beta <- Beta[match, match]
  gc()
  
  #'----------------------------------------------------------------------------
  #' Compute mantel test
  mantel_Beta <- mantel_rcpp(mat1 = Beta, 
                             mat2 = spectra_disimilarity,
                             nperm = 999,
                             threads = threads,
                             progress = TRUE)
  
  stadistics <- (sum(mantel_Beta$perm_stats >= mantel_Beta$observed_stat) + 1)/(999 + 1)
  
  #'----------------------------------------------------------------------------
  #' Melt frames
  
  # As data.table
  spectral_distance <- as.data.table(spectra_disimilarity)
  spectral_distance$PLOT_CN <- colnames(spectra_disimilarity)
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
  
  # Merge frames for results
  frame <- merge(spectral_distance, 
                 Beta, 
                 by = c("PLOT_CN", "PLOT_CN_col"), 
                 all.x = TRUE,
                 all.y = TRUE)
  
  return_list <- list(mantel = mantel_Beta$observed_stat,
                      mantel_permutation = mantel_Beta$perm_stats[,1],
                      mantel_significance = stadistics,
                      frame = frame)
  
  saveRDS(return_list, paste0(output_path, "/",
                              sensor, "_", diversity, 
                              "_", beta_type, "_mantel.rds"))
  
}

#'------------------------------------------------------------------------------
#' Mantel test

#---------------------------------------------------------
#' @DESIS

root_path <- "/media/antonio/Work/MMDB"
sensor <- "DESIS"

# Spectra from plots used
# data <- readRDS(paste0(root_path, "/Models/01_model_beta/DESIS/DESIS_taxo_data.rds"))
# data <- data[, c("scene", "PLT_CN")]

# # Spectra avaiable but contaminated
# spectra <- fread(paste0(root_path, "/Pixel_extraction/03_match-spectra/DESIS/DESIS_spectra_SG.csv"))
# spectra <- merge(data, spectra, by = c("scene", "PLT_CN"), all.x = TRUE, all.y = FALSE)

# Select bands for sam
# spectra <- spectra[, .SD, .SDcols = c(2, 53:228)] #"446.53" to "896.53"
# spectra_disimilarity <- sam(as.matrix(spectra[, -1]), threads = 26, progress = TRUE)
# PLT_CN <- spectra$PLT_CN

# colnames(spectra_disimilarity) <- as.character(PLT_CN)
# rownames(spectra_disimilarity) <- as.character(PLT_CN)
# spectra_disimilarity <- round(spectra_disimilarity, 4)
# 
# saveRDS(spectra_disimilarity, paste0(root_path, "/Models/01_model_beta/DESIS/spectra_disimilarity_SG.rds"))
spectra_disimilarity <- readRDS(paste0(root_path, "/Models/01_model_beta/DESIS/spectra_disimilarity_SG.rds"))

# Select output
output_path <- paste0(root_path, "/Diversity/03_Mantel")
threads <- 22

#---------------------------------------------------------
# Taxonomic
diversity <- "taxo"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger.rds")

mantel_test(sensor = "DESIS",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "taxo",
            beta_path = beta_path,
            beta_type = "total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "taxo",
            beta_path = beta_path,
            beta_type = "replacement",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "taxo",
            beta_path = beta_path,
            beta_type = "richness",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "taxo",
            beta_path = beta_path,  
            beta_type = "replacement_total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "taxo",
            beta_path = beta_path,  
            beta_type = "richness_total",
            output_path = output_path,
            threads = threads)
gc()

# Phylogenetic
diversity <- "phylo"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger_norm.rds")

mantel_test(sensor = "DESIS",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "phylo",
            beta_path = beta_path,
            beta_type = "total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "phylo",
            beta_path = beta_path,
            beta_type = "replacement",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "phylo",
            beta_path = beta_path,
            beta_type = "richness",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "phylo",
            beta_path = beta_path,  
            beta_type = "replacement_total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "phylo",
            beta_path = beta_path,  
            beta_type = "richness_total",
            output_path = output_path,
            threads = threads)
gc()


# Functional
diversity <- "func"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger_norm.rds")

mantel_test(sensor = "DESIS",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "func",
            beta_path = beta_path,
            beta_type = "total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "func",
            beta_path = beta_path,
            beta_type = "replacement",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "func",
            beta_path = beta_path,
            beta_type = "richness",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "func",
            beta_path = beta_path,  
            beta_type = "replacement_total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "DESIS", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "func",
            beta_path = beta_path,  
            beta_type = "richness_total",
            output_path = output_path,
            threads = threads)
gc()

#---------------------------------------------------------
#' @EMIT

root_path <- "/media/antonio/Work/MMDB"
sensor <- "EMIT"

# Spectra from plots used
data <- readRDS(paste0(root_path, "/Models/01_model_beta/EMIT/EMIT_taxo_data.rds"))
data <- data[, c("scene", "PLT_CN")]

# Spectra avaiable but contaminated
# spectra <- fread(paste0(root_path, "/Pixel_extraction/03_match-spectra/EMIT/EMIT_spectra_SG.csv"))
# spectra <- merge(data, spectra, by = c("scene", "PLT_CN"), all.x = TRUE, all.y = FALSE)
# 
# # Select bands for sam
# spectra <- spectra[, .SD, .SDcols = c(2,
#                                       42:139, #441 - 1176
#                                       184:210, #1498 - 1701
#                                       266:307)] #2113 - 2421
# 
# spectra_disimilarity <- sam(as.matrix(spectra[, -1]), threads = 26, progress = TRUE)
# 
# colnames(spectra_disimilarity) <- as.character(spectra$PLT_CN)
# rownames(spectra_disimilarity) <- as.character(spectra$PLT_CN)
# 
# spectra_disimilarity <- round(spectra_disimilarity, 4)
# 
# # 
# saveRDS(spectra_disimilarity, paste0(root_path, "/Models/01_model_beta/EMIT/spectra_disimilarity_SG.rds"))
spectra_disimilarity <- readRDS(paste0(root_path, "/Models/01_model_beta/EMIT/spectra_disimilarity_SG.rds"))

# Select output
output_path <- paste0(root_path, "/Diversity/03_Mantel")
threads <- 20

#---------------------------------------------------------
# Taxonomic
diversity <- "taxo"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger.rds")

mantel_test(sensor = "EMIT",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "taxo",
            beta_path = beta_path,
            beta_type = "total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "taxo",
            beta_path = beta_path,
            beta_type = "replacement",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "taxo",
            beta_path = beta_path,
            beta_type = "richness",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "taxo",
            beta_path = beta_path,  
            beta_type = "replacement_total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "taxo",
            beta_path = beta_path,  
            beta_type = "richness_total",
            output_path = output_path,
            threads = threads)
gc()

# Phylogenetic
diversity <- "phylo"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger_norm.rds")

mantel_test(sensor = "EMIT",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "phylo",
            beta_path = beta_path,
            beta_type = "total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "phylo",
            beta_path = beta_path,
            beta_type = "replacement",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "phylo",
            beta_path = beta_path,
            beta_type = "richness",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "phylo",
            beta_path = beta_path,  
            beta_type = "replacement_total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "phylo",
            beta_path = beta_path,  
            beta_type = "richness_total",
            output_path = output_path,
            threads = threads)
gc()

# Functional
diversity <- "func"
beta_path <- paste0(root_path, "/Diversity/01_beta/beta_", diversity, "_VOL_hellinger_norm.rds")

mantel_test(sensor = "EMIT",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "func",
            beta_path = beta_path,
            beta_type = "total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "func",
            beta_path = beta_path,
            beta_type = "replacement",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT",
            spectra_disimilarity = spectra_disimilarity,
            diversity = "func",
            beta_path = beta_path,
            beta_type = "richness",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "func",
            beta_path = beta_path,  
            beta_type = "replacement_total",
            output_path = output_path,
            threads = threads)
gc()

mantel_test(sensor = "EMIT", 
            spectra_disimilarity = spectra_disimilarity,
            diversity = "func",
            beta_path = beta_path,  
            beta_type = "richness_total",
            output_path = output_path,
            threads = threads)
gc()

