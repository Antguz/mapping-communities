pls_ordination_model <- function(meta, ordination, spectra, split, repeated = 100, opt_components, threads = 1) {
  
  # Creates partition for training and testing
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  dimentions <- ncol(ordination_training)
  
  #Names
  performance_names <- c("MDS", "repeated", "Model", "N", "Rsq", "Bias", "RMSE", "RMSE_P")
  coefficients_names <- c("MDS", "repeated", "intercept", colnames(spectra_training))
  VIP_names <- c("MDS", "repeated", colnames(spectra_training))
  
  # Progress bar
  pb <- txtProgressBar(min = 1, max = repeated*dimentions, style = 3) 
  
  #Collectors
  VIP_export <- data.table()
  coeff_export <- data.table()
  perf_export <- data.table()
  
  #Create partitions and folds----------------------------------
  #Model collectors
  data_split <- list()
  data_folds <- list()
  
  for(i in 1:repeated) {
    
    #Split
    split <- createDataPartition(meta_training$COUNTYCD, 
                                 p = 0.80,
                                 list = FALSE, 
                                 times = 1)
    
    #Folds
    folds <- CreateSpacetimeFolds(x = meta_training[split[,1],], 
                                  spacevar = "COUNTYCD",
                                  seed = sample(1:10000, 1))
    folds <- folds$index
    
    #Collect
    data_split[[i]] <- split[,1]
    data_folds[[i]] <- folds
    
    #Residuals
    rm(list = c("split", "folds"))
    
  }
  
  # Parallel
  if(threads > 1) {
    cl <- makeCluster(threads)
    doSNOW::registerDoSNOW(cl)
  }
  
  #Loop of repeated
  for(q in 1:dimentions) {
    
    # Train set
    train_master <- cbind(ordination_training[, .SD, .SDcols = c(q)],
                          spectra_training)
    colnames(train_master)[1] <- "MDS"
    
    # Progress bar
    pb <- txtProgressBar(min = 1, max = repeated, style = 3) 
    
    #Loop over ordination dimensions
    for(i in 1:repeated) {
      
      setTxtProgressBar(pb, i)
      
      #Basic information
      info <- data.table(MDS = q,
                         repeate = i)
      
      #Training
      pls_components <- plsr(MDS ~ ., 
                             data= train_master[data_split[[i]],], 
                             scale = FALSE,
                             center = TRUE,
                             ncomp= opt_components[q],
                             validation = "CV",
                             method = "oscorespls",
                             segments = data_folds[[i]],
                             na.action = na.omit)
      #Performance
      predict_training <- predict(pls_components, 
                                  train_master[, 2:ncol(train_master)], 
                                  ncomp = opt_components[q], 
                                  type= "response")[,,1]
      
      perf_training <- c("Training", performance(train_master$MDS, predict_training))
      perf_training <- as.data.table(matrix(perf_training, nrow = 1))
      perf_training <- cbind(info, perf_training)
      colnames(perf_training) <- performance_names
      
      #Coefficients
      coefficients <- as.vector(coef(pls_components, 
                                     ncomp = opt_components[q], 
                                     intercept=TRUE))
      coefficients <- as.data.table(matrix(coefficients, nrow = 1))
      coefficients <- cbind(info, coefficients)
      colnames(coefficients) <- coefficients_names
      
      #VIP
      vip <- as.vector(VIP(pls_components, opt.comp = opt_components[q]))
      vip <- as.data.table(matrix(vip, nrow = 1))
      vip <- cbind(info, vip)
      colnames(vip) <- VIP_names
      
      #Collect
      perf_export <- rbind(perf_export, perf_training)
      coeff_export <- rbind(coeff_export, coefficients)
      VIP_export <- rbind(VIP_export, vip)
      
    }
  }
  
  if(threads > 1) {
    close(pb)
    stopCluster(cl)
  }
  
  #Export
  return(list(Coefficients = coeff_export, 
              VIP = VIP_export, 
              Performance = perf_export))
    
}
