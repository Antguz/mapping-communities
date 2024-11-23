# Cross-validation strategy
iscv_strategy <- function(meta, 
                          vec_projected, 
                          split, 
                          segments,
                          threads = 4) {
  
  iscv_application <- function(X,
                               meta, 
                               vec_projected, 
                               split, 
                               segments) {
    
    # Get vector for training and testing
    split_split <- 1:nrow(meta) %in% split
    split_seg <- 1:nrow(meta[split,]) %in% segments[[X]]
    
    # Get training and testing vectors
    train_vec <- vec_projected[split_split,]
    train_vec <- train_vec[split_seg,]
    test_vec <- vec_projected[!split_split,] #[!split_split,]
    
    knndm_folds <- knndm(tpoints = train_vec,
                         predpoints = test_vec,
                         space = "geographical",
                         k = 10,
                         maxp = 0.5,
                         clustering = "kmeans",
                         sampling = "regular")
    
    return(knndm_folds$indx_train)
    
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 1:length(segments), 
                         FUN = iscv_application, 
                         meta = meta, 
                         vec_projected = vec_projected, 
                         split = split,
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
