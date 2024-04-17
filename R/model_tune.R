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
    
    PRESS <- ml_model$validation$PRESS
    
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
    
    PRESS <- cbind(data.table(iteration = rep(t, ncol_ord),
                              type = rep("PRESS", ncol_ord),
                              estimate = rep("CV", ncol_ord),
                              Axes = paste0("Axis", 1:ncol_ord),
                              '(Intercept)' = c(0, 0, 0)),
                   as.data.table(PRESS))
    
    # Collect results
    results <- rbind(RMSEP_train, RMSEP_CV,
                     R2_train, R2_CV, PRESS)
    
    # Clean memory
    rm(list = c("sub_ordination", "ncol_ord", "sub_spectra", 
                "ml_model", "RMSEP_train", "RMSEP_CV", "R2_train", "R2_CV", "PRESS"))
    
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

selectNcomp_axes <- function(object,
                             method = "onesigma",
                             plot = FALSE) {
  
  # Optimal
  optimal <- NA
  Rmseps <- RMSEP(object, "CV")$val
  
  for(i in 1:dim(Rmseps)[2]) {
    
    rmseps <- Rmseps[,i,]
    maxIdx <- object$ncomp + 1
    absBest <- which.min(rmseps[seq_len(maxIdx)])
    
    if (absBest > 0) {
      
      if (is.null(origResponse <- object$y)) {
        origResponse <- c(predict(object, ncomp = 1)[,i,] + residuals(object)[,i,1])
      }
      
      allresids <- cbind(origResponse - (sum(origResponse) - origResponse) / (length(origResponse) - 1),
                         object$validation$pred[,i,] - origResponse)
      
      residsds <- apply(allresids, 2, sd) / sqrt(nrow(allresids))
      uls <- rmseps - residsds
      selection <- min(which(uls < rmseps[absBest])) - 1
      
      
      if (isTRUE(plot)) {
        xvals <- seq_along(rmseps) - 1
        plot(xvals, rmseps, ylab = "RMSEP",
             xlab = "Number of components", type = "b")
        arrows(xvals, uls, xvals, rmseps + residsds,
               code = 3, col = "gray", angle = 90, length = .1)
        abline(h = rmseps[absBest], col = "gray", lty = 3)
        abline(v = absBest - 1, col = "gray", lty = 3)
        abline(v = selection, col = "blue", lty = 2)
        legend("topright", legend = c("Abs. minimum", "Selection"),
               lty = 3:2, col = c("gray", "blue"), bg = "white")
      }
    }
    
    optimal[i] <- selection
    
  }
  
  return(optimal)
  
}
  
