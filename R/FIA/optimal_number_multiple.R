optimal_number <- function(meta, 
                           ordination, 
                           spectra, 
                           weigths,
                           split, 
                           int_ncomp = 50, 
                           repeated = 10,
                           threads = 4) {
  
  # Creates partition for training and testing-------------
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  weights_training <- weigths[split]
  
  #Results collector
  results <- data.table()
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = repeated, style = 3) 
  
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
  
  #Loop over ordination dimensions
  for(i in 1:repeated) {
    
    MDS <- ordination_training[data_folds[[i]],]
    bands <- spectra_training[data_folds[[i]],]
    we <- weights_training[data_folds[[i]]]
    
    #Training
    pls_components <- plsr(as.matrix(MDS) ~ as.matrix(bands), 
                           scale = FALSE,
                           center = TRUE,
                           ncomp= int_ncomp,
                           validation = "CV",
                           method = "oscorespls",
                           na.action = na.omit,
                           maxit = 1000)
    
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
