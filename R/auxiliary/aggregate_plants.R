#'------------------------------------------------------------------------------
#' @title Aggregate abundance based on plants groups
#'------------------------------------------------------------------------------

#' @description Aggregate (sum) community values from species to genus, family, 
#' or any categorical variable

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(bit64)

#'------------------------------------------------------------------------------
#' Arguments

#' @param meta A frame of two columns: Species names and plant group to aggregate
#' @param communities A frame of relative abundance of species

#'------------------------------------------------------------------------------
#' Function

aggregate_plants <- function(meta, communities) {
  
  # Get PLT_CN
  PLT_CN <- as.integer64(rownames(communities))
  communities <- as.data.table(communities)
  
  # Get taxom
  sp_meta <- meta$Taxon2
  sp_communities <- colnames(communities)
  
  # Match frames
  sp <- match(sp_communities, sp_meta)
  
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
    sub_community <- communities[, .SD, .SDcols = sub_group$ID]
    
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

# Apply in batch
aggregate_batch <- function(meta, communities) {
  
  n_cols <- ncol(meta)
  
  results <- data.table(PLT_CN = as.integer64(rownames(communities)))
  
  for(i in 2:n_cols) {
    
    sub_meta <- meta[, .SD, .SDcols = c(1, i)]
    sub_class <- aggregate_plants(sub_meta, communities)
    results <- cbind(results,  sub_class[, -1])
    
  }
  
  return(results)
  
}
