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
  
  #Loop over ordination dimensions
  for(i in 1:repeated) {
    
    train <- data.frame(MDS = I(as.data.frame(ordination_training[data_split,][[i]])),
                        spectra = I(as.data.frame(spectra_training[data_split,][[i]])))
    
    #Training
    pls_components <- plsr(as.matrix(ordination_training) ~ as.matrix(spectra_training), 
                           data= train,
                           scale = FALSE,
                           center = TRUE,
                           ncomp= int_ncomp,
                           validation = "CV",
                           method = "oscorespls",
                           #segments = data_folds[[i]],
                           na.action = na.omit)
    
    #Get optimal component
    press <- data.table(MDS = paste0("MDS", 1:ncol(ordination_training)),
                        iteration = i)
    press <- cbind(press, as.data.table(pls_components$validation$PRESS))
    
    #Export
    results <- rbind(results, prediction)
    
  }

  if(threads > 1 ) {
    stopCluster(cl)
  }
  
  return(results)
  
}
