#'------------------------------------------------------------------------------
#' @title Compute the community-weighted mean
#'------------------------------------------------------------------------------

#' @description Estimate the community community-weighted mean

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

CWM <- function(meta, communities) {
  
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
  meta_trans <- as.data.table(t(meta_new[,2]))
  colnames(meta_trans) <- meta_new$Taxon2
  meta_trans <- meta_trans[rep(1, length(PLT_CN)),]
  
  # Get species to aggregate
  trait <- colnames(meta_new)[2]
  
  # Matrix multiplication
  cwm <- as.matrix(communities) * as.matrix(meta_trans)
  cwm <- rowSums(cwm)
  
  # Return
  complete <- data.table(PLT_CN = PLT_CN,
                         trait = cwm)
  colnames(complete)[2] <- trait
  
  return(complete)
  
}


# Apply in batch
CWM_batch <- function(meta, communities) {
  
  n_cols <- ncol(meta)
  
  results <- data.table(PLT_CN = as.integer64(rownames(communities)))
  
  for(i in 2:n_cols) {
    
    sub_meta <- meta[, .SD, .SDcols = c(1, i)]
    sub_class <- CWM(sub_meta, communities)
    results <- cbind(results,  sub_class[, -1])
    
  }
  
  return(results)
  
}
