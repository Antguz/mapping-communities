model_build <- function(meta_training, 
                        ordination_training, 
                        spectra_training,
                        segments,
                        iscv,
                        tune_optimal, 
                        ntimes,
                        threads = 26) {
  
  # Application loop
  application_build <- function(X, 
                                meta_training, 
                                ordination_training, 
                                spectra_training,
                                segments,
                                iscv,
                                tune_optimal) {
    
    # Get segments
    axis1_data <- cbind(ordination_training[segments[[X]]$axis1_seg, 1], 
                        spectra_training[segments[[X]]$axis1_seg,])
    
    axis2_data <- cbind(ordination_training[segments[[X]]$axis2_seg, 2], 
                        spectra_training[segments[[X]]$axis2_seg,])
    
    axis3_data <- cbind(ordination_training[segments[[X]]$axis3_seg, 3], 
                        spectra_training[segments[[X]]$axis3_seg,])
    
    # Get models
    axis1_model <- plsr(Axis1 ~ ., 
                        data = axis1_data,
                        validation = "CV",
                        segments = iscv[[X]]$folds_axis1,
                        ncomp = tune_optimal[1],
                        center = TRUE,
                        method = "oscorespls",
                        maxit = 10000)
    
    axis2_model <- plsr(Axis2 ~ ., 
                        data = axis2_data,
                        validation = "CV",
                        segments = iscv[[X]]$folds_axis2,
                        ncomp = tune_optimal[2],
                        center = TRUE,
                        method = "oscorespls",
                        maxit = 10000)
    
    axis3_model <- plsr(Axis3 ~ ., 
                        data = axis3_data,
                        validation = "CV",
                        segments = iscv[[X]]$folds_axis3,
                        ncomp = tune_optimal[3],
                        center = TRUE,
                        method = "oscorespls",
                        maxit = 10000)
    
    # Clean memory
    rm(list = c("axis1_data", "axis2_data", "axis3_data"))
    
    return(list(axis1_model = axis1_model,
                axis2_model = axis2_model,
                axis3_model = axis3_model))
    
  }
  
  # Parallel processing
  complete <- mclapply(X = 1:ntimes, 
                       FUN = application_build, 
                       meta_training = meta_training, 
                       ordination_training = ordination_training, 
                       spectra_training = spectra_training,
                       tune_optimal = tune_optimal,
                       segments = segments,
                       iscv = iscv,
                       mc.preschedule = TRUE, 
                       mc.set.seed = FALSE,
                       mc.silent = FALSE, 
                       mc.cores = threads,
                       mc.cleanup = TRUE, 
                       mc.allow.recursive = TRUE, 
                       affinity.list = NULL)
  
  return(complete)
  
}
