pls_ndms_model <- function(meta, NDMS, spectra, split = split, folds = 10, repeated = 10, opt_components, threads = 1) {
  
  # Creates partition for training and testing
  NDMS_training <- NDMS[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  dimentions <- ncol(NDMS_training)
  
  #Create spatial points
  spatial <- sf::st_as_sf(meta_training, coords = c("xcoord", "ycoord"))
  st_crs(spatial) <- "ESRI:102004"
  
  # Progress bar
  pb <- txtProgressBar(min = 1, max = repeated, style = 3) 
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Parallel
  if(threads > 1) {
    cl <- makeCluster(threads)
    doSNOW::registerDoSNOW(cl)
  }
  
  #Loop over iterations
  final <- foreach(i = 1:repeated, .packages= c('base', 'data.table', 'foreach', 'caret', 'pls', 'plsVarSel', 'blockCV'), .export= c("performance", "CV_spatialBlock"), .options.snow = opts) %dopar% {
    
    #Spatial cross-validation 
    sub_split <- CV_spatialBlock(spatial, k = folds)
    
    #Collectors
    VIP_export <- data.table()
    coeff_export <- data.table()
    perf_export <- data.table()
    
    for(j in 1:folds) {
      
      #Train and test sets of folds
      trainSet <- unlist(sub_split$folds[[j]][1]) 
      testSet <- unlist(sub_split$folds[[j]][2])
      
      # Creates partition for training and testing
      NDMS_trainSet <- NDMS_training[trainSet,]
      NDMS_testSet <- NDMS_training[testSet,]
      spectra_trainSet <- spectra_training[trainSet,]
      spectra_testSet <- spectra_training[testSet, ]
      
      for(q in 1:dimentions) {
        
        #Basic information
        info <- data.table(NDMS = q,
                           repeated = i,
                           fold = j)
        
        #Train set
        train <- cbind(NDMS_trainSet[, .SD, .SDcols = c(q)],
                       spectra_trainSet)
        colnames(train)[1] <- "NDMS"
        
        #Testing set
        test <- cbind(NDMS_testSet[, .SD, .SDcols = c(q)],
                      spectra_testSet)
        colnames(test)[1] <- "NDMS"
        
        #Training
        pls_components <- plsr(NDMS ~ ., 
                               data= train, 
                               scale = FALSE,
                               center = TRUE,
                               ncomp= opt_components[q],
                               validation = "none", 
                               method = "oscorespls")
        
        #Predict
        predict_testing <- predict(pls_components, 
                                   test, 
                                   ncomp = opt_components[q], 
                                   type= "response")[,,1]
        perf <- performance(test$NDMS, predict_testing)
        perf <- as.data.table(matrix(perf, nrow = 1))
        perf <- cbind(info, perf)
        
        #Coefficients
        coefficients <- as.vector(coef(pls_components, 
                                       ncomp = opt_components[q], 
                                       intercept=TRUE))
        coefficients <- as.data.table(matrix(coefficients, nrow = 1))
        coefficients <- cbind(info, coefficients)
        
        #VIP
        vip <- as.vector(VIP(pls_components, opt.comp = opt_components[q]))
        vip <- as.data.table(matrix(vip, nrow = 1))
        vip <- cbind(info, vip)
        
        #Collect
        perf_export <- rbind(perf_export, perf)
        coeff_export <- rbind(coeff_export, coefficients)
        VIP_export <- rbind(VIP_export, vip)
        
        #Remove
        rm(list = c("info", "train", "test", "pls_components", 
                    "predict_testing", "perf", "coefficients", "vip"))
        
      }
      
      #Remove
      rm(list = c("trainSet", "testSet", "NDMS_trainSet", 
                  "NDMS_testSet", "spectra_trainSet", "spectra_testSet"))
      
    }
    
    #Remove
    rm(list = c("sub_split"))
    
    #Export
    return(list(Coefficients = coeff_export, 
                VIP = VIP_export, 
                Performance = perf_export))
    
  }
  
  if(threads > 1) {
    close(pb)
    stopCluster(cl)
  }
  
  new_final <- combine_custom(final)
  
  return(new_final)
  
}

#Combine-----------------------------------
combine_custom <- function(final) {
  n <- length(final)
  coefficients <- NULL
  vip <- NULL
  performance <- NULL
  
  for (i in seq(n)) {
    coefficients <- rbind(coefficients, final[[i]]$Coefficients)
    vip <- rbind(vip, final[[i]]$VIP)
    performance <- rbind(performance, final[[i]]$Performance)
    
  }
  
  return(list(Coefficients = coefficients, 
              VIP = vip, 
              Performance = performance))
  
}

