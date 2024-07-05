model_CWM <- function(fia_meta,
                      ordination,
                      aggregation,
                      data_split,
                      iterations = 100,
                      threads = 4) {
  
  # Get number of plant groups
  groups <- ncol(aggregation)
  
  # Get model application in batch
  apply_occurrence <- function(X,
                               fia_meta,
                               ordination,
                               aggregation,
                               data_split,
                               iterations) {
    
    #---------------------------------------------------------------------------
    # Match databases
    if(all.equal(aggregation$PLT_CN, ordination$PLT_CN) == FALSE) {
      stop("Files most match")
    }
    
    # Group name fir future reference
    group_name <- colnames(aggregation)[X]
    
    #---------------------------------------------------------------------------
    # Reshape data
    frame <- cbind(aggregation[, .SD, .SDcols = c(X)], ordination[, 2:4])
    colnames(frame) <- c("trait", "Axis1", "Axis2", "Axis3") 
    
    #---------------------------------------------------------------------------
    # Multiple iterations
    
    data_training <- frame[as.matrix(data_split)[, c(X-1)] == 1,]
    data_testing <- frame[as.matrix(data_split)[, c(X-1)] == 0,]
    cwm_collector <- data.table()
    
    for(ii in 1:iterations) {
      
      #---------------------------------------------------------------------------
      # Down sampling
      down_train <- data_training[sample(1:nrow(data_training), 
                                         ceiling(nrow(data_training)*0.8)),]
      
      #---------------------------------------------------------------------------
      # Perform model
      
      # Model features
      ctrl <- trainControl(method = "cv", 
                           number = 10,
                           savePredictions = TRUE)
      
      # Model
      plsFit <- train(trait ~ .,
                      data = down_train,
                      method = "glm",
                      trControl = ctrl)
      
      #---------------------------------------------------------------------------
      # Performance on training
      
      training <- data.frame(obs = data_training$trait,
                             pred = predict(plsFit, 
                                            newdata = data_training))
      
      training_performance <- parameters(training$obs, training$pred)
      
      #---------------------------------------------------------------------------
      # Performance on testing
      
      testing <- data.frame(obs = data_testing$trait,
                            pred = predict(plsFit, newdata = data_testing))
      
      testing_performance <- parameters(testing$obs, testing$pred)
      
      #---------------------------------------------------------------------------
      # Collect data
      collect <- rbind(training_performance, testing_performance)
      collect <- cbind(collect, matrix(plsFit$finalModel$coefficients, nrow = 1))
      colnames(collect)[7:10] <- names(plsFit$finalModel$coefficients)
      collect <- cbind(plant = c(group_name, group_name), 
                       iteration = c(ii, ii),
                       dataset = c("training", "testing"),
                       psamples = c(nrow(training)/(nrow(training)+nrow(testing))*100,
                                    nrow(testing)/(nrow(training)+nrow(testing))*100),
                       collect)
      
      cwm_collector <- rbind(cwm_collector, collect)
    
    }
    
    return(cwm_collector)
    
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 2:groups, 
                         FUN = apply_occurrence, 
                         fia_meta = fia_meta,
                         ordination = ordination,
                         aggregation = aggregation,
                         data_split = data_split,
                         iterations = iterations,
                         mc.preschedule = TRUE, 
                         mc.set.seed = FALSE,
                         mc.silent = FALSE, 
                         mc.cores = 25,
                         mc.cleanup = TRUE, 
                         mc.allow.recursive = TRUE, 
                         affinity.list = NULL)
  
  complete <- do.call(rbind.data.frame, complete)
  
  return(complete)
  
}
