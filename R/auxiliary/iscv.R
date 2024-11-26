# Cross-validation strategy
iscv_strategy <- function(meta_split, 
                          segments,
                          threads = 4) {
  
  meta_split$ID <- paste0(meta_split$STATECD, "_", meta_split$COUNTYCD)
  
  iscv_application <- function(X,
                               meta_split,
                               segments) {
    
    # Get vector for training and testing
    axis1_split_seg <- 1:nrow(meta_split) %in% segments[[X]]$axis1_seg
    axis2_split_seg <- 1:nrow(meta_split) %in% segments[[X]]$axis2_seg
    axis3_split_seg <- 1:nrow(meta_split) %in% segments[[X]]$axis3_seg
    
    # Split mete
    axis1_meta_seg <- meta_split[axis1_split_seg,]
    axis2_meta_seg <- meta_split[axis2_split_seg,]
    axis3_meta_seg <- meta_split[axis3_split_seg,]
    
    # Folds
    axis1folds <- CreateSpacetimeFolds(x = axis1_meta_seg,
                                       spacevar = "ID",
                                       k = 10)
    axis2folds <- CreateSpacetimeFolds(x = axis2_meta_seg,
                                       spacevar = "ID",
                                       k = 10)
    axis3folds <- CreateSpacetimeFolds(x = axis3_meta_seg,
                                       spacevar = "ID",
                                       k = 10)
    
    return(list(folds_axis1 = axis1folds$index,
                folds_axis2 = axis2folds$index,
                folds_axis3 = axis3folds$index))
    
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 1:length(segments), 
                         FUN = iscv_application, 
                         meta_split = meta_split, 
                         segments = segments,
                         mc.preschedule = TRUE, 
                         mc.set.seed = TRUE,
                         mc.silent = FALSE, 
                         mc.cores = threads,
                         mc.cleanup = TRUE, 
                         mc.allow.recursive = TRUE, 
                         affinity.list = NULL)
  
  return(complete)
  
}
