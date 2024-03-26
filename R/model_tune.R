model_tune <- function(meta_training, 
                       ordination_training, 
                       spectra_training,
                       segments,
                       iscv,
                       ncomp_max,
                       ntimes,
                       threads) {
  
  # Application loop
  application_tune <- function(X, 
                               meta_training, 
                               ordination_training, 
                               spectra_training,
                               ncomp_max,
                               segments,
                               iscv) {
    
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
                     ncomp = ncomp_max,
                     center = TRUE,
                     method = "oscorespls",
                     maxit = 1000)

    # Performance 
    RMSEP_train <- pls::RMSEP(object = ml_model,
                              estimate = "train")
    
    RMSEP_CV <- pls::RMSEP(object = ml_model,
                           estimate = "CV")
    
    R2_train <- pls::R2(object = ml_model,
                        estimate = "train")
    
    R2_CV <- pls::R2(object = ml_model,
                     estimate = "CV")
    
    # Reshape
    RMSEP_train <- cbind(data.table(iteration = rep(t, ncol_ord),
                                    type = rep("RMSEP", ncol_ord),
                                    estimate = rep("train", ncol_ord),
                                    Axes = paste0("Axis", 1:ncol_ord)),
                         as.data.table(RMSEP_train$val[1,1:ncol_ord,]))
    
    RMSEP_CV <- cbind(data.table(iteration = rep(t, ncol_ord),
                                 type = rep("RMSEP", ncol_ord),
                                 estimate = rep("CV", ncol_ord),
                                 Axes = paste0("Axis", 1:ncol_ord)),
                      as.data.table(RMSEP_CV$val[1,1:ncol_ord,]))
    
    R2_train <- cbind(data.table(iteration = rep(t, ncol_ord),
                                 type = rep("R2", ncol_ord),
                                 estimate = rep("train", ncol_ord),
                                 Axes = paste0("Axis", 1:ncol_ord)),
                      as.data.table(R2_train$val[1,1:ncol_ord,]))
    
    R2_CV <- cbind(data.table(iteration = rep(t, ncol_ord),
                              type = rep("R2", ncol_ord),
                              estimate = rep("CV", ncol_ord),
                              Axes = paste0("Axis", 1:ncol_ord)),
                   as.data.table(R2_CV$val[1,1:ncol_ord,]))
    
    # Collect results
    results <- rbind(RMSEP_train, RMSEP_CV,
                     R2_train, R2_CV)
    
    # Clean memory
    rm(list = c("sub_ordination", "ncol_ord", "sub_spectra", 
                "ml_model", "RMSEP_train", "RMSEP_CV", "R2_train", "R2_CV"))
    
    return(results)
    
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 1:ntimes, 
                         FUN = application_tune, 
                         meta_training = meta_training, 
                         ordination_training = ordination_training, 
                         spectra_training = spectra_training,
                         ncomp_max = ncomp_max,
                         segments = segments,
                         iscv = iscv,
                         mc.preschedule = TRUE, 
                         mc.set.seed = FALSE,
                         mc.silent = FALSE, 
                         mc.cores = threads,
                         mc.cleanup = TRUE, 
                         mc.allow.recursive = TRUE, 
                         affinity.list = NULL)
  
  frame <- data.table()
  
  for(i in 1:length(complete)) {
    frame <- rbind(frame, complete[[i]], fill = TRUE)
  }
  
  return(frame)
  
}


