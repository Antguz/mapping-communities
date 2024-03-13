#'------------------------------------------------------------------------------
#' @title Aggregate community values from species to genus
#'------------------------------------------------------------------------------

#' @description Aggregate community values from species to genus, family, or any
#' categorical variable

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(bit64)

#'------------------------------------------------------------------------------
#' Function

aggregate_species <- function(meta, community) {
  
  # Get PLT_CN
  PLT_CN <- as.integer64(rownames(community))
  community <- as.data.table(community)
  
  # Get taxom
  sp_meta <- meta$Taxon2
  sp_community <- colnames(community)
  
  # Match frames
  sp <- match(sp_community, sp_meta)
  
  # Create new meta
  meta_new <- meta[sp, ]
  meta_new$ID <- 1:nrow(meta_new)
  
  # Get species to aggregate
  colnames(meta_new)[2] <- "groups"
  group <- unique(meta_new$groups)
  
  # Collector
  complete <- data.table(PLT_CN = PLT_CN)
  
  # Apply loop to groups
  for(i in 1:length(group)) {
    
    # Subset groups to aggregate
    sub_group <- subset(meta_new, groups == group[i])
    sub_community <- community[, .SD, .SDcols = sub_group$ID]
    
    # Aggregate groups and create frame
    sub_community <- rowSums(sub_community)
    sub_community <- data.table(sub_community)
    colnames(sub_community) <- group[i]
    
    # Collect frame
    complete <- cbind(complete, sub_community)
    
  }
  
  # Provide 
  
  return(complete)
  
}

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"

meta <- fread(paste0(root_path, "/00_FIA/SPCD_meta.csv"))
meta <- meta[, c(9, 4)]
#meta <- meta[, c(9, 10)]
community <- readRDS(paste0(root_path, "/03_diversity/master_matched_comm_phylo-VOL-log.rds"))
