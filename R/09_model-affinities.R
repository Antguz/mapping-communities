#'------------------------------------------------------------------------------
#' @title Model probabilities of species affinities
#'------------------------------------------------------------------------------

#' @description It models the probability of occurence based on species affinities.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(bit64)

#'------------------------------------------------------------------------------
#' Arguments
#' 
#' @param fia
#' @param affinity

#'------------------------------------------------------------------------------
#' Function

affinity_prob <- function(fia,
                          affinity) {
  
  #---------------------------------------------------------
  # Match databases
  PLT_CN <- affinity[, "PLT_CN"]
  fia$PLT_CN <- as.integer64(rownames(fia))
  fia <- as.data.table(fia)
  fia <- merge(PLT_CN, fia, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
  
  # Get col names
  col_affinity <- colnames(affinity)[-1]
  col_fia <- colnames(fia)[-1]
  col_match <- match(col_fia, col_affinity)
  
  # Get binomial model of species occurrence
  binomial <- function(X, fia, affinity, col_match) {
    
    X <- 216
    i <- col_match[X] + 1
    
    # Occurrence
    species_occu <- fia[, .SD, .SDcols = i]
    sp_o <- colnames(species_occu)
    colnames(species_occu) <- "occu"
    
    nobs <- nrow(species_occu[occu > 0,])
    
    # Affinity
    species_affi <- affinity[, .SD, .SDcols = i]
    sp_a <- colnames(species_affi)
    colnames(species_affi) <- "affi"
    
    if(!all.equal(sp_o, sp_a)) {
      stop()
    }
    
    # Join
    frame <- cbind(species_occu, species_affinity)
    #plot(frame$affinity, frame$occurrence)
    #frame[occu > 0, occu := 1] # Presence/absence
    
    # Model
    model <- glm(occu ~ affi, data = frame, family = binomial(link='logit'))
    plot(frame$affi, frame$occu)
    points(frame$affi, predict(model, frame, type='response'), col = "red")
    
  }
  
}

#'------------------------------------------------------------------------------
#' @example

# Run in a single scene

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"

sensor <- "DESIS"
type <- "taxo"
data <- "training"

# Estimate affinity
affinity <- fread(paste0(root_path, "/07_affinity/",
                         sensor, "_", type, "_", data, "_affinity.csv"))

fia <- readRDS(paste0(root_path, "/03_diversity/", "master_matched_comm_", type, "-VOL.rds"))
fia[fia > 0] <- 1
colSums(fia)
