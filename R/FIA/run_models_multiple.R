run_models_multiple <- function(meta,
                                ordination, 
                                spectra, 
                                split, 
                                int_ncomp = 19, 
                                repeated = 100,
                                threads = 26) {

  
  # Creates partition for training and testing-------------
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  weights_training <- weigths[split]
  
  #Results collector
  models <- list()
  
  #Create partitions and folds--3---------------------------
  #Model collectors
  data_folds <- list()
  
  for(i in 1:repeated) {
    
    #Folds
    folds <- grid_selection(ordination_training, steps = 30)
    
    data_folds[[i]] <- folds
    
    #Residuals
    rm(list = c("folds"))
    
  }
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = repeated, style = 3) 
  
  #Parallel 
  if(threads > 1 ) {
    cl <- makeCluster(threads)
    pls.options(parallel = cl)
  }
  
  #Loop over ordination dimensions
  for(i in 1:repeated) {
    
    #Progress
    setTxtProgressBar(pb, i)
    
    #Select
    MDS <- ordination_training[data_folds[[i]],]
    bands <- spectra_training[data_folds[[i]],]
    county <- meta_training[data_folds[[i]],]
    ncounties <- length(unique(county$COUNTYCD))
    
    #we <- weights_training[data_folds[[i]]]
    #Spatial folds
    folds <- CreateSpacetimeFolds(county,
                                  spacevar = "COUNTYCD",
                                  k = 10)
    
    #Training
    pls_components <- plsr(as.matrix(MDS) ~ as.matrix(bands), 
                           scale = FALSE,
                           center = TRUE,
                           ncomp= int_ncomp,
                           validation = "CV",
                           method = "oscorespls",
                           segments = folds$index,
                           na.action = na.omit,
                           maxit = 10000,
                           jackknife = TRUE)
    
    #Get optimal component
    models[[i]] <- pls_components
    
  }
  
  if(threads > 1 ) {
    stopCluster(cl)
  }
  
  return(models)
  
}
