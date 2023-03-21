optimal_sigma <- function(meta, ordination, spectra, split, threads = 24) {
  
  # Creates partition for training and testing-------------
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  dimentions <- ncol(ordination_training)
  data_weigths <- 1/meta_training$ntaxa
  data_weigths <- data_weigths/sum(data_weigths)
  
  #Results collector
  results <- list()
  
  #Loop of repeated
  for(q in 1:dimentions) {
    
    #Train set
    train_master <- cbind(ordination_training[, .SD, .SDcols = c(q)],
                          spectra_training)
    colnames(train_master)[1] <- "MDS"
    
    #Control
    tune <- expand.grid(sigma = seq(0.001, 0.2, length.out = 80))
    cfitControl <- trainControl(method = "cv",
                                classProbs = FALSE,
                                savePredictions = TRUE,
                                allowParallel = TRUE)
    
    if(threads > 1) {
      cl <- makePSOCKcluster(threads)
      registerDoParallel(cl)
    }
    
    #Model
    ml_model <- train(MDS ~ ., 
                      data = train_master, 
                      method = "gaussprRadial", 
                      trControl = cfitControl,
                      weights = data_weigths,
                      metric = "RMSE",
                      tuneGrid = tune)
    
    if(threads > 1) {
      stopCluster(cl)
      Sys.sleep(2)
      gc()
    }
    
    results[[q]] <- ml_model
    rm("ml_model")
    gc()
    
  }
    
  return(results)
    
}
