#'------------------------------------------------------------------------------
#' @title Get proportions of angiosperms/gymnosperms and mycorrhizal
#'------------------------------------------------------------------------------

#' @description It gets the proportions per stand.
#' @return .txt file with the proportions

#'------------------------------------------------------------------------------
#' Library

library(data.table)

#root_path <- "E:/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
for_proportions <- readRDS(paste0(root_path, "/03_diversity/for_proportions_comm.rds"))
out_path <- paste0(root_path, "/03_diversity")

# Get unique
unique_PLT_CN <- unique(for_proportions$PLT_CN)

# Loop to get proportions
results <- data.table()

# Create progress bar
pb <- txtProgressBar(min = 0, max = length(unique_PLT_CN), style = 3)

for(i in 1:length(unique_PLT_CN)) {
  
  # Update progress bar
  setTxtProgressBar(pb, i)
  
  # Plot of interest
  plot <- subset(for_proportions, PLT_CN == unique_PLT_CN[i])
  
  A <- as.numeric(sum(plot[Gymnosperm  == "N", "BA_TPA"]))
  G <- as.numeric(sum(plot[Gymnosperm  == "Y", "BA_TPA"]))
  
  if(is.na(A)) {
    proportion_A <- 0
  } else if(is.na(G)) {
    proportion_A <- 1
  } else {
    Total <- A + G
    proportion_A <- round(A/Total, 2)
  }
  
  
  AM <- as.numeric(sum(plot[mycorrhizal_association == "AM", "BA_TPA"]))
  EM <- as.numeric(sum(plot[mycorrhizal_association == "ECM", "BA_TPA"]))
  
  if(is.na(AM)) {
    proportion_AM <- 0
  } else if(is.na(EM)) {
    proportion_AM <- 1
  } else {
    Total_M <- AM + EM
    proportion_AM <- round(AM/Total_M, 2)
  }
  
  nspecies <- nrow(plot)
  
  if(nspecies == 1) {
    species <- plot$Taxon2[1]
  } else {
    species <- NA
  }
  
  ntrees <- sum(plot$ntrees)
  BA_TPA <- sum(plot$BA_TPA)
  DRYBIO_AG_TPA <- sum(plot$DRYBIO_AG_TPA)
  VOL_TPA <- sum(plot$VOL_TPA)
  
  to_results <- data.table(PLT_CN = unique_PLT_CN[i],
                           SR = nspecies,
                           species = species,
                           ntrees = ntrees,
                           BA_TPA = BA_TPA,
                           DRYBIO_AG_TPA= DRYBIO_AG_TPA,
                           VOL_TPA= VOL_TPA,
                           A_to_T = proportion_A,
                           AM_to_T = proportion_AM)
  
  results <- rbind(results, to_results)
  
}

fwrite(results, paste0(out_path, "/", "proportions.csv"))
