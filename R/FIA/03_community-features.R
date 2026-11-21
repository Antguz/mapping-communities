#'------------------------------------------------------------------------------
#' @title Get relative abundance or CWM of community features
#'------------------------------------------------------------------------------

#' @description It estimates the relative abundance (categorical features) or
#' the community-weighted mean (continuous features) of species features per
#' inventory plot, using wood volume (VOL_TPA) as weights.
#'
#' @return A .csv file (proportions.csv) with the relative abundances and CWMs
#' per plot (PLT_CN).

#'------------------------------------------------------------------------------
#' @Library

library(data.table)

#'------------------------------------------------------------------------------
#' @Paths-and-files

root_path <- "/media/antonio/Work/MMDB"
out_path <- paste0(root_path, "/FIA/00_FIA")
fia <- readRDS(paste0(root_path, "/FIA/01_FIA-reshape/fia_reshape.rds"))
SPCD <- fread(paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv"))
voi <- SPCD[, c("SPCD", "growth_form")]
features <- SPCD[, c(1, 10:26)]

#'------------------------------------------------------------------------------
#' Function
get_proportions <- function(fia, voi) {
  
  # Modify voi
  variable_name <- colnames(voi)[2]
  colnames(voi) <- c("SPCD", "variable")
  
  # Merge databases
  frame <- merge(fia, voi, by = "SPCD", all.x = TRUE, all.y = FALSE)
  
  if(class(frame$variable) == "character") {
    
    types <- unique(frame$variable)
    summary_frame <- frame[, sum(VOL_TPA), by = c("PLT_CN", "variable")]
    colnames(summary_frame)[3] <- "total"
    summary_frame <- dcast(summary_frame, PLT_CN ~ variable, value.var = "total")
    summary_frame$total <- apply(summary_frame[, .SD, .SDcols = 2:ncol(summary_frame)], 1, sum, na.rm = TRUE)
    summary_frame$proportion <- summary_frame[, .SD, .SDcols = 2]/summary_frame$total
    final <- summary_frame[, .SD, .SDcols = c(1, ncol(summary_frame))]
    final[is.na(proportion), proportion := 0]
    colnames(final)[2] <- paste0("rel-abund_", colnames(summary_frame)[2])
    
  } else {
    
    summary_frame <- frame[, total := sum(VOL_TPA), by = c("PLT_CN")]
    summary_frame <- summary_frame[order(PLT_CN)]
    summary_frame$relative <- summary_frame$VOL_TPA/summary_frame$total
    summary_frame$weighted <- summary_frame$relative*summary_frame$variable
    final <- summary_frame[, sum(weighted), by = c("PLT_CN")]
    colnames(final)[2] <- paste0("CWM_", variable_name)
    
  }
  
  return(final)
  
}

#'------------------------------------------------------------------------------
#' Function by trait
features_summary <- function(fia, features) {
  
  nfeatures <- ncol(features) - 1
  complete <- unique(fia[, "PLT_CN"])

  for(i in 2:nfeatures) {
    
    voi <- features[, .SD, .SDcols = c(1, i)]
    sub <- get_proportions(fia, voi)
    complete <- merge(complete, sub, by = "PLT_CN", all.x = TRUE, all.y = TRUE)
    
  }
  
  return(complete)
  
}

#'------------------------------------------------------------------------------
#' @example

results <- features_summary(fia, features)
fwrite(results, paste0(out_path, "/", "proportions.csv"))
