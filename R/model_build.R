modul_build <- function(meta, 
                        ordination, 
                        spectra, 
                        split, 
                        weigths, 
                        method = "gaussprRadial",
                        tune_optimal, 
                        ntimes = 100, 
                        threads = 26) {
  
  # Creates partition for training and testing-------------
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  dimentions <- ncol(ordination_training)
  data_weigths <- weigths[split]
  data_weigths <- data_weigths/sum(data_weigths)
  
  #Results collector
  final_collector <- list()
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = ntimes, style = 3) 
  
  #Loop of dimensions
  for(q in 1:dimentions) {
    
    #Train set
    train_master <- cbind(ordination_training[, .SD, .SDcols = c(q)],
                          spectra_training)
    colnames(train_master)[1] <- "MDS"
    
    #Progress bar
    print(paste0("MDS", q))
    pb <- txtProgressBar(min = 1, max = ntimes, style = 3) 
    
    #Model collector
    model_collector <- list()
    
    if(threads > 1) {
      cl <- makePSOCKcluster(threads)
      registerDoParallel(cl)
    }
    
    #Loop over repetitions
    for(t in 1:ntimes) {
      
      #Progress
      setTxtProgressBar(pb, t)
      
      #Homogeneous segmentation
      segmentation <- mds_segmentation(train_master$MDS)
      sub_master <- train_master[segmentation, ]
      sub_weigths <- data_weigths[segmentation]
      sub_meta <- meta_training[segmentation, ]
      
      #Create spatial folds
      folds <- CreateSpacetimeFolds(x = sub_meta,
                                    spacevar = "COUNTYCD", 
                                    seed = sample(1:10000, 1))
      folds <- folds$index
      
      #Control
      cfitControl <- trainControl(method = "cv",
                                  number = 10,
                                  index = folds,
                                  classProbs = FALSE,
                                  savePredictions = TRUE,
                                  allowParallel = TRUE)
      
      #Model
      ml_model <- train(MDS ~ ., 
                        data = sub_master, 
                        method = method, 
                        trControl = cfitControl,
                        weights = sub_weigths,
                        preProc = c("center"),
                        tuneGrid = tune_optimal[[q]])
      
      model_collector[[t]] <- ml_model
      
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
