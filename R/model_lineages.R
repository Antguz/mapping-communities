model_lineages <- function(fia_meta,
                           ordination,
                           aggregation,
                           data_split,
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
    
    # Group name for future reference
    group_name <- colnames(aggregation)[X]
    
    #---------------------------------------------------------------------------
    # Reshape data
    frame <- cbind(aggregation[, .SD, .SDcols = X], ordination[, 2:4])
    colnames(frame) <- c("occurrence", "Axis1", "Axis2", "Axis3") 
    frame[occurrence > 0, occurrence := 1]
    frame$occurrence <- as.character(frame$occurrence)
    frame[occurrence == "0", occurrence := "absence"]
    frame[occurrence == "1", occurrence := "presence"]
    frame$occurrence <- as.factor(frame$occurrence)
    
    #---------------------------------------------------------------------------
    # Perform model
    
    # Model features
    ctrl <- trainControl(method = "repeatedcv", 
                         number = 10,
                         repeats = 10,
                         classProbs = TRUE,
                         summaryFunction = twoClassSummary,
                         savePredictions = TRUE)
    
    # Model
    plsFit <- train(occurrence ~ .,
                    data = frame[as.matrix(data_split)[,(X-1)] == 1,],
                    method = "glm",
                    trControl = ctrl,
                    metric = "ROC")
    
    #---------------------------------------------------------------------------
    # Performance on training
    
    training <- data.frame(obs = frame$occurrence[as.matrix(data_split)[,(X-1)] == 1],
                           pred = predict(plsFit, 
                                          newdata = frame[as.matrix(data_split)[,(X-1)] == 1,]),
                           presence = predict(plsFit, 
                                              newdata = frame[as.matrix(data_split)[,(X-1)] == 1,], 
                                              type = "prob")[,2],
                           absence = predict(plsFit, 
                                             newdata = frame[as.matrix(data_split)[,(X-1)] == 1,], 
                                             type = "prob")[,1])
    
    training_performance <- c(defaultSummary(training, 
                                             lev = levels(frame$occurrence)),
                              twoClassSummary(training, 
                                              lev = levels(frame$occurrence)),
                              prSummary(training, 
                                        lev = levels(frame$occurrence)),
                              mnLogLoss(training, 
                                        lev = levels(frame$occurrence)))
    
    #---------------------------------------------------------------------------
    # Performance on testing
    
    testing <- data.frame(obs = frame[as.matrix(data_split)[,(X-1)] == 0, ]$occurrence,
                          pred = predict(plsFit, newdata = frame[as.matrix(data_split)[,(X-1)] == 0,]),
                          presence = predict(plsFit, 
                                             newdata = frame[as.matrix(data_split)[,(X-1)] == 0,], 
                                             type = "prob")[,2],
                          absence = predict(plsFit, 
                                            newdata = frame[as.matrix(data_split)[,(X-1)] == 0,], 
                                            type = "prob")[,1])
    
    testing_performance <- c(defaultSummary(testing, 
                                            lev = levels(frame$occurrence)),
                             twoClassSummary(testing, 
                                             lev = levels(frame$occurrence)),
                             prSummary(testing, 
                                       lev = levels(frame$occurrence)),
                             mnLogLoss(testing, 
                                       lev = levels(frame$occurrence)))
    
    #---------------------------------------------------------------------------
    # Collect data
    training <- as.data.table(training)
    testing <- as.data.table(testing)
    #training[obs == "presence",]
    #testing
    
    collect <- rbind(c(training_performance, plsFit$finalModel$coefficients),
                     c(testing_performance, plsFit$finalModel$coefficients))
    collect <- as.data.table(collect)
    collect <- cbind(plant = c(group_name, group_name), 
                     dataset = c("training", "testing"),
                     psamples = c(nrow(training[obs == "presence",])/nrow(training)*100,
                                  nrow(testing[obs == "presence",])/nrow(testing)*100),
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
