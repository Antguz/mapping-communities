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
    
    t <- X
    
    # Homogeneous segments
    sub_ordination <- ordination_training[segments[[t]],]
    ncol_ord <- ncol(sub_ordination)
    sub_spectra <- spectra_training[segments[[t]],]
    
    # Model
    model_frame <- data.frame(axes = I(as.matrix(sub_ordination)),
                              bands = I(as.matrix(sub_spectra)))
    ml_model <- plsr(axes ~ bands, 
                     data = model_frame,
                     validation = "CV",
                     segments = iscv[[t]],
                     ncomp = max(tune_optimal),
                     center = TRUE,
                     method = "oscorespls",
                     maxit = 10000)
    
    return(ml_model)
    
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
