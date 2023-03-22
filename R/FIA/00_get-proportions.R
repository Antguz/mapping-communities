#'------------------------------------------------------------------------------
#' @title Get proportions of angiosperms/gymnosperms and mycorrhizal
#'------------------------------------------------------------------------------

#' @description It gets the proportions per stand.
#' @return .txt file with the proportions

#'------------------------------------------------------------------------------
#' Library

library(data.table)

root_path <- "E:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
fia <- fread(paste0(root_path, "/00_FIA/trees_WI_MN_20221102.csv"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"), na.strings= "NA")
out_path <- paste0(root_path, "/03_reshape")
SPCD <- na.exclude(SPCD)

# States of interest
fia <- subset(fia, STATECD == 27 | STATECD == 55)
fia <- fia[order(PLT_CN, SUBP)]

# Merge species code
fia <- merge(fia, SPCD[, c("SPCD", "PLANTS-Code", "A_G", "AM_EM")], by = "SPCD", all.x = TRUE, all.y = TRUE)

# Get Basal Area
fia$BA <- pi*(fia$DIA/2)^2

# Basal area multiplied by TPA
fia$BA_TPA <- fia$BA*fia$TPA_UNADJ
fia$DRYBIO_AG_TPA <- fia$DRYBIO_AG*fia$TPA_UNADJ

# Remove dead trees
fia <- fia[STATUSCD == 1]

# Get plots with complete information
complete_plots <- fia[, sum(is.na(DRYBIO_AG_TPA)), by = "PLT_CN"]
complete_plots <- complete_plots[V1 == 0]

# Select plots with complete information
fia <- merge(complete_plots[, 1], 
             fia, 
             by = "PLT_CN", 
             all.x = TRUE, 
             all.y = FALSE)

#Summarize
fia_summary <- fia[, .(DRYBIO_AG_TPA = sum(DRYBIO_AG_TPA, na.rm = TRUE)), 
                   by = c("PLOT", "PLT_CN", "A_G", "AM_EM")]
fia_summary <- fia_summary[order(PLT_CN)]

# Get unique
unique_PLT_CN <- unique(fia_summary$PLT_CN)

# Loop to get proportions
results <- data.table()

for(i in 1:length(unique_PLT_CN)) {
  
  plot <- subset(fia_summary, PLT_CN == unique_PLT_CN[i])
  
  A <- as.numeric(sum(plot[A_G == "A", "DRYBIO_AG_TPA"]))
  G <- as.numeric(sum(plot[A_G == "G", "DRYBIO_AG_TPA"]))
  
  if(is.na(A)) {
    proportion_A <- 0
  } else if(is.na(G)) {
    proportion_A <- 1
  } else {
    Total <- A + G
    proportion_A <- round(A/Total, 2)
  }
  
  
  AM <- as.numeric(sum(plot[AM_EM == "AM", "DRYBIO_AG_TPA"]))
  EM <- as.numeric(sum(plot[AM_EM == "EM", "DRYBIO_AG_TPA"]))
  
  if(is.na(AM)) {
    proportion_AM <- 0
  } else if(is.na(EM)) {
    proportion_AM <- 1
  } else {
    Total_M <- AM + EM
    proportion_AM <- round(AM/Total_M, 2)
  }
  
  
  to_results <- data.table(PLT_CN = unique_PLT_CN[i],
                           A_to_T = proportion_A,
                           AM_to_T = proportion_AM)
  
  results <- rbind(results, to_results)
  
}

fwrite(results, paste0(out_path, "/", "proportions_DRYBIO.csv"))
