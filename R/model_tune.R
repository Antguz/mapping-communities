model_tune <- function(meta,
                       ordination, 
                       spectra, 
                       split,
                       ncomp_max,
                       ntimes,
                       threads) {
  
  # Creates partition for training and testing-------------
  meta_training <- meta[split,]
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  
  # Loop of repeated
  application <- function(X, 
                          meta_training, 
                          ordination_training, 
                          spectra_training,
                          ncomp_max) {
    
    t <- X
    
    # Segmentation
    segmentation <- data_split(meta = meta_training, ptraining = 0.9)
    
    # Homogeneous segmentation
    sub_meta <- meta_training[segmentation, ]
    sub_ordination <- ordination_training[segmentation,]
    ncol_ord <- ncol(sub_ordination)
    sub_spectra <- spectra_training[segmentation,]
    
    # Create spatial folds
    folds <- CreateSpacetimeFolds(x = sub_meta,
                                  spacevar = "STATECD", 
                                  class = "COUNTYCD",
                                  k = length(unique(sub_meta$STATECD)),
                                  seed = sample(1:10000, 1))
    folds <- folds$index
    
    # Model
    model_frame <- data.frame(axes = I(as.matrix(sub_ordination)),
                              bands = I(as.matrix(sub_spectra)))
    ml_model <- plsr(axes ~ bands, 
                     data = model_frame,
                     validation = "CV",
                     segments = folds,
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
                                    Axes = paste0("PCoA", 1:ncol_ord)),
                         as.data.table(RMSEP_train$val[1,1:ncol_ord,]))
    
    RMSEP_CV <- cbind(data.table(iteration = rep(t, ncol_ord),
                                 type = rep("RMSEP", ncol_ord),
                                 estimate = rep("CV", ncol_ord),
                                 Axes = paste0("PCoA", 1:ncol_ord)),
                      as.data.table(RMSEP_CV$val[1,1:ncol_ord,]))
    
    R2_train <- cbind(data.table(iteration = rep(t, ncol_ord),
                                 type = rep("R2", ncol_ord),
                                 estimate = rep("train", ncol_ord),
                                 Axes = paste0("PCoA", 1:ncol_ord)),
                      as.data.table(R2_train$val[1,1:ncol_ord,]))
    
    R2_CV <- cbind(data.table(iteration = rep(t, ncol_ord),
                              type = rep("R2", ncol_ord),
                              estimate = rep("CV", ncol_ord),
                              Axes = paste0("PCoA", 1:ncol_ord)),
                   as.data.table(R2_CV$val[1,1:ncol_ord,]))
    
    # Collect results
    results <- rbind(RMSEP_train, RMSEP_CV,
                     R2_train, R2_CV)
    
    # Clean memory
    #rm(list = c("segmentation", 
    #            "sub_meta", "sub_ordination", "ncol_ord", "sub_spectra", 
    #            "folds", "ml_model", 
    #            "RMSEP_train", "RMSEP_CV", "R2_train", "R2_CV"))
    
    return(results)
    
  }
  
  # Parallel processing
  complete <- mclapply(X = 1:ntimes, 
                       FUN = application, 
                       meta_training = meta_training, 
                       ordination_training = ordination_training, 
                       spectra_training = spectra_training,
                       ncomp_max = ncomp_max,
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


