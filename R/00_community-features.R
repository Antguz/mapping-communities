#'------------------------------------------------------------------------------
#' @title Get proportions of angiosperms/gymnosperms and mycorrhizal
#'------------------------------------------------------------------------------

#' @description It gets the proportions per stand.
#' @return .txt file with the proportions

#'------------------------------------------------------------------------------
#' Library

library(data.table)

#root_path <- "E:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"

out_path <- paste0(root_path, "/00_FIA")
fia <- readRDS(paste0(root_path, "/03_reshape/fia_reshape.rds"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD_meta.csv"))
#voi <- SPCD[, c("SPCD", "growth_form")]
voi <- SPCD[, c("SPCD", "shade")]
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
    colnames(final)[2] <- paste0("proportion_", colnames(summary_frame)[2])
    
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
  complete <- data.table(fia[, "PLT_CN"])

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
