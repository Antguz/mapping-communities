model_CWM <- function(fia_meta,
                      ordination,
                      aggregation,
                      data_split = data_split,
                      threads = 4) {
  
  # Get number of plant groups
  groups <- ncol(aggregation) - 1
  
  # Get model application in batch
  apply_occurrence <- function(X,
                               fia_meta,
                               ordination,
                               aggregation,
                               data_split) {
    
    #---------------------------------------------------------------------------
    # Match databases
    if(all.equal(aggregation$PLT_CN, ordination$PLT_CN) == FALSE) {
      stop("Files most match")
    }
    
    # Group name fir future reference
    group_name <- colnames(aggregation)[X]
    
    #---------------------------------------------------------------------------
    # Reshape data
    frame <- cbind(aggregation[, .SD, .SDcols = X], ordination[, 2:4])
    colnames(frame) <- c("trait", "Axis1", "Axis2", "Axis3") 
    
    #---------------------------------------------------------------------------
    # Perform model
    
    # Model features
    ctrl <- trainControl(method = "repeatedcv", 
                         number = 10,
                         repeats = 10,
                         savePredictions = TRUE)
    
    # Model
    plsFit <- train(trait ~ .,
                    data = frame[as.matrix(data_split)[,(X-1)] == 1,],
                    method = "glm",
                    trControl = ctrl)
    
    #---------------------------------------------------------------------------
    # Performance on training
    
    training <- data.frame(obs = frame$trait[as.matrix(data_split)[,(X-1)] == 1],
                           pred = predict(plsFit, 
                                          newdata = frame[as.matrix(data_split)[,(X-1)] == 1,]))
    
    training_performance <- parameters(training$obs, training$pred)
    
    #---------------------------------------------------------------------------
    # Performance on testing
    
    testing <- data.frame(obs = frame[as.matrix(data_split)[,(X-1)] == 0, ]$trait,
                          pred = predict(plsFit, newdata = frame[as.matrix(data_split)[,(X-1)] == 0,]))
    
    testing_performance <- parameters(testing$obs, testing$pred)
    
    #---------------------------------------------------------------------------
    # Collect data
    training <- as.data.table(training)
    testing <- as.data.table(testing)
    
    collect <- rbind(c(training_performance, plsFit$finalModel$coefficients),
                     c(testing_performance, plsFit$finalModel$coefficients))
    collect <- as.data.table(collect)
    collect <- cbind(plant = c(group_name, group_name), 
                     dataset = c("training", "testing"),
                     psamples = c(nrow(training)/(nrow(training)+nrow(testing))*100,
                                  nrow(testing)/(nrow(training)+nrow(testing))*100),
                     collect)
    
    return(collect)
    
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 2:groups, 
                         FUN = apply_occurrence, 
                         fia_meta = fia_meta,
                         ordination = ordination,
                         aggregation = aggregation,
                         data_split = data_split,
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
