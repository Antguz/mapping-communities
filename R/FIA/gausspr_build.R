gausspr_build <- function(meta, ordination, spectra, split, sigma_components, repeated = 100, threads = 16) {
  
  # Creates partition for training and testing-------------
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  dimentions <- ncol(ordination_training)
  
  #Results collector
  results <- data.table()
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = repeated, style = 3) 
  
  #Create partitions and folds---------------------------------------
  #Model collectors
  data_split <- list()
  data_folds <- list()
  data_weigths <- list()
  
  #Create partitions
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
    
    #Weigths
    weigths <- 1/meta_training[split[,1],]$ntaxa
    weigths <- weigths/sum(weigths)
    
    #Collect
    data_split[[i]] <- split[,1]
    data_folds[[i]] <- folds
    data_weigths[[i]] <- weigths
    
    #Residuals
    rm(list = c("split", "folds", "weigths"))
    
  }
  
  #Final model collector
  final_collector <- list()
  
  #Loop of repeated
  for(q in 1:dimentions) {
    
    #Train set
    train_master <- cbind(ordination_training[, .SD, .SDcols = c(q)],
                          spectra_training)
    colnames(train_master)[1] <- "MDS"
    
    #Progress bar
    pb <- txtProgressBar(min = 1, max = repeated, style = 3) 
    cat(paste0("Model on dimension ", q))
    
    #Model collector
    model_collector <- list()
    
    if(threads > 1) {
      cl <- makePSOCKcluster(threads)
      registerDoParallel(cl)
    }
    
    #Loop over ordination dimensions
    for(i in 1:repeated) {
      
      #Bar
      setTxtProgressBar(pb, i)
      
      #Control
      cfitControl <- trainControl(method = "cv",
                                  index = data_folds[[i]],
                                  classProbs = FALSE,
                                  savePredictions = TRUE,
                                  allowParallel = TRUE)
      
      #Model
      ml_model <- train(MDS ~ ., 
                        data = train_master[data_split[[i]],], 
                        method = "gaussprRadial", 
                        trControl = cfitControl,
                        tune = expand.grid(sigma = sigma_components[q]),
                        weights = data_weigths[[i]])
      
      model_collector[[i]] <- ml_model
      
    }
    
    if(threads > 1) {
      stopCluster(cl)
      Sys.sleep(2)
      gc()
    }
    
    #Collect model
    final_collector[[q]] <- model_collector
    
  }
  
  return(final_collector)
  
}
