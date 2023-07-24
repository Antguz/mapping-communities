model_tune <- function(meta, 
                       ordination, 
                       spectra, 
                       split, 
                       weigths,
                       method = "gaussprRadial",
                       tune_range,
                       ntimes = 10,
                       threads = 24) {
  
  # Creates partition for training and testing-------------
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  dimentions <- ncol(ordination_training)
  data_weigths <- weigths[split]
 
  #Results collector
  results <- data.table()
  
  #Loop of repeated
  for(q in 1:dimentions) {
    
    #Train set
    train_master <- cbind(ordination_training[, .SD, .SDcols = c(q)],
                          spectra_training)
    colnames(train_master)[1] <- "MDS"
    
    #Print progress
    print(paste0("MDS", q))
    pb <- txtProgressBar(min = 0, max = ntimes, style = 3)
    
    #Parallel implementation
    if(threads > 1) {
      cl <- makePSOCKcluster(threads)
      registerDoParallel(cl)
    }
    
    for(t in 1:ntimes) {
      
      #Progress
      setTxtProgressBar(pb, t)
      
      # Homogeneous segmentation
      segmentation <- mds_segmentation(train_master$MDS, proportions = 0.75)
      sub_master <- train_master[segmentation, ]
      sub_weigths <- data_weigths[segmentation]
      sub_meta <- meta_training[segmentation, ]
      
      # Create spatial folds
      folds <- CreateSpacetimeFolds(x = sub_meta,
                                    spacevar = "STATECD", 
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
                        metric = "RMSE",
                        tuneGrid = tune_range)
      
      ncomp <- ml_model$bestTune$ncomp
      ml_model <- as.data.table(ml_model$results)
      ml_model$MDS <- paste0("MDS", q)
      ml_model$iteration <- t
      ml_model$bestTune <- ncomp
      
      #Collect results
      results <- rbind(results, ml_model)
      
      rm(list = c("ml_model", "segmentation", "sub_master",
                  "sub_weigths", "sub_meta", "folds",
                  "cfitControl"))
      gc()
      
    }
    
    #Close parallel
    if(threads > 1) {
      stopCluster(cl)
      Sys.sleep(2)
      gc()
    }
    
  }
  
  return(results)
    
}
