optimal_number <- function(meta, ordination, spectra, split, int_ncomp = 50, repeated = 10, k = 10, threads = 16) {
  
  # Creates partition for training and testing-------------
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  dimentions <- ncol(ordination_training)

  #Results collector
  results <- data.table()
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = repeated*dimentions, style = 3) 
  
  #Create partitions and folds--3---------------------------
  #Model collectors
  data_folds <- list()
  
  for(i in 1:repeated) {
    
    #Folds
    folds <- CreateSpacetimeFolds(x = meta_training, 
                                  spacevar = "COUNTYCD",
                                  seed = sample(1:10000, 1))
    folds <- folds$index
    
    data_folds[[i]] <- folds
    
    #Residuals
    rm(list = c("folds"))
    
  }
  
  #Parallel 
  if(threads > 1 ) {
    cl <- makeCluster(threads)
    pls.options(parallel = cl)
  }
  
  #Loop of repeated
  for(q in 1:dimentions) {
    
    #Train set
    train_master <- cbind(ordination_training[, .SD, .SDcols = c(q)],
                          spectra_training)
    colnames(train_master)[1] <- "MDS"
    
    #Loop over ordination dimensions
    for(i in 1:repeated) {
      
      #Training
      pls_components <- plsr(MDS ~ ., 
                             data= train_master[data_split[[i]],], 
                             scale = FALSE,
                             center = TRUE,
                             ncomp= int_ncomp,
                             validation = "CV",
                             method = "oscorespls",
                             segments = data_folds[[i]],
                             na.action = na.omit)
      
      #Get optimal component
      component <- selectNcomp(pls_components, "onesigma")
      PRESS_min <- which.min(as.vector(pls_components$validation$PRESS))
      CV <- which.min(as.vector(RMSEP(pls_components)$val[1,,]))
      plot(as.vector(pls_components$validation$PRESS))
      plot(pls_components, ncomp = CV)
      
      #Frame to collect
      prediction <- data.table(ordination = colnames(ordination_training)[q],
                               Repeated = i,
                               Onesigma = component,
                               PRESS = PRESS_min,
                               CV = CV,
                               samples = nrow(train))
      
      results <- rbind(results, prediction)
      
    }
  }
  
  if(threads > 1 ) {
    stopCluster(cl)
  }
  
  return(results)
  
}
