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
                        ncomp = ncomp_max,
                        center = TRUE,
                        method = "oscorespls",
                        maxit = 10000)
    
    axis2_model <- plsr(Axis2 ~ ., 
                        data = axis2_data,
                        validation = "CV",
                        segments = iscv[[X]]$folds_axis2,
                        ncomp = ncomp_max,
                        center = TRUE,
                        method = "oscorespls",
                        maxit = 10000)
    
    axis3_model <- plsr(Axis3 ~ ., 
                        data = axis3_data,
                        validation = "CV",
                        segments = iscv[[X]]$folds_axis3,
                        ncomp = ncomp_max,
                        center = TRUE,
                        method = "oscorespls",
                        maxit = 10000)
    
    # Clean memory
    rm(list = c("axis1_data", "axis2_data", "axis3_data"))

    # Performance 
    axis1_RMSEP_train <- pls::RMSEP(object = axis1_model, estimate = "train")$val[,1,]
    axis2_RMSEP_train <- pls::RMSEP(object = axis2_model, estimate = "train")$val[,1,]
    axis3_RMSEP_train <- pls::RMSEP(object = axis3_model, estimate = "train")$val[,1,]
    
    axis1_RMSEP_CV <- pls::RMSEP(object = axis1_model, estimate = "CV")$val[,1,]
    axis2_RMSEP_CV <- pls::RMSEP(object = axis2_model, estimate = "CV")$val[,1,]
    axis3_RMSEP_CV <- pls::RMSEP(object = axis3_model, estimate = "CV")$val[,1,]
    
    axis1_R2_train <- pls::R2(object = axis1_model, estimate = "train")$val[,1,]
    axis2_R2_train <- pls::R2(object = axis2_model, estimate = "train")$val[,1,]
    axis3_R2_train <- pls::R2(object = axis3_model, estimate = "train")$val[,1,]
    
    axis1_R2_CV <- pls::R2(object = axis1_model, estimate = "CV")$val[,1,]
    axis2_R2_CV <- pls::R2(object = axis2_model, estimate = "CV")$val[,1,]
    axis3_R2_CV <- pls::R2(object = axis3_model, estimate = "CV")$val[,1,]
    
    axis1_PRESS <- c(Intercept = 0, axis1_model$validation$PRESS)
    axis2_PRESS <- c(Intercept = 0, axis2_model$validation$PRESS)
    axis3_PRESS <- c(Intercept = 0, axis3_model$validation$PRESS)
    
    onasigma <- c(selectNcomp(axis1_model, "onesigma"),
                  selectNcomp(axis2_model, "onesigma"),
                  selectNcomp(axis3_model, "onesigma"))
    
    randomization <- c(selectNcomp(axis1_model, "randomization"),
                       selectNcomp(axis2_model, "randomization"),
                       selectNcomp(axis3_model, "randomization"))
    
    min_value <- c(which.min(axis1_RMSEP_CV[-1]),
                   which.min(axis2_RMSEP_CV[-1]),
                   which.min(axis3_RMSEP_CV[-1]))
    
    # Merge
    results <- rbind(axis1_RMSEP_train, axis2_RMSEP_train, axis3_RMSEP_train,
                     axis1_RMSEP_CV, axis2_RMSEP_CV, axis3_RMSEP_CV, 
                     axis1_R2_train, axis2_R2_train, axis3_R2_train,
                     axis1_R2_CV, axis2_R2_CV, axis3_R2_CV, 
                     axis1_PRESS, axis2_PRESS, axis3_PRESS)
    
    base <- data.table(Iteration = X,
                       Axes = rep(c("Axis1", "Axis2", "Axis3"), 5),
                       Metric = c(rep("RMSEP", 6), rep("R2", 6), rep("PRESS", 3)),
                       Estimate = c(rep("train", 3), rep("CV", 3), rep("train", 3), rep("CV", 3), rep("PRESS", 3)),
                       onasigma = rep(onasigma, 5),
                       randomization = rep(randomization, 5),
                       abs_min = rep(min_value, 5))
    
    results <- cbind(base, results)
    
    # Clean memory
    rm(list = c("axis1_model", "axis2_model", "axis3_model",
                "axis1_RMSEP_train", "axis2_RMSEP_train", "axis3_RMSEP_train",
                "axis1_RMSEP_CV", "axis2_RMSEP_CV", "axis3_RMSEP_CV",
                "axis1_R2_train", "axis2_R2_train", "axis3_R2_train",
                "axis1_R2_CV", "axis2_R2_CV", "axis3_R2_CV",
                "axis1_PRESS", "axis2_PRESS", "axis3_PRESS",
                "onasigma", "randomization", "min_value",
                "base"))
                
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
  
  complete <- do.call(rbind, complete)
  return(complete)
  
}

