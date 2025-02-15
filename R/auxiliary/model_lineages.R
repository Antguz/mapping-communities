model_lineages <- function(ordination,
                           aggregation,
                           data_split,
                           iterations = 100,
                           threads = 4) {
  
  # Get number of plant groups
  groups <- ncol(aggregation) 
  
  # Get model application in batch
  apply_occurrence <- function(X,
                               ordination,
                               aggregation,
                               data_split,
                               iterations = 100) {
    
    #---------------------------------------------------------------------------
    # Match databases
    if(all.equal(aggregation$PLT_CN, ordination$PLT_CN) == FALSE) {
      stop("Files most match")
    }
    
    # Group name for future reference
    group_name <- colnames(aggregation)[X]
    
    #---------------------------------------------------------------------------
    # Reshape data
    frame <- cbind(aggregation[, .SD, .SDcols = c(X)], ordination[, .SD, .SDcols = c(2:ncol(ordination))])
    colnames(frame) <- c("occurrence", colnames(ordination[, .SD, .SDcols = c(2:ncol(ordination))])) 
    frame$occurrence <- as.character(frame$occurrence)
    frame[occurrence == "0", occurrence := "absence"]
    frame[occurrence == "1", occurrence := "presence"]
    frame$occurrence <- as.factor(frame$occurrence)
    frame$occurrence <- factor(frame$occurrence, levels = c("presence", "absence"))
    
    #---------------------------------------------------------------------------
    # Dealing with class imbalance
    
    data_training <- frame[as.matrix(data_split)[, c(X-1)] == 1,]
    data_testing <- frame[as.matrix(data_split)[, c(X-1)] == 0,]
    lineages_collector <- data.table()
    
    for(ii in 1:iterations) {
      
      #---------------------------------------------------------------------------
      # Down sampling
      split <- split_lineages(data_training$occurrence, split = 0.85, seed = ii)
      
      #---------------------------------------------------------------------------
      # Class balance
      down_train <- ROSE(occurrence ~ ., data = data_training[split,])$data
      down_train$occurrence <- factor(down_train$occurrence, levels = c("presence", "absence"))
      
      # down_train <- upSample(data_training[split, -1], data_training[split, ]$occurrence)
      # colnames(down_train)[ncol(down_train)] <- "occurrence"
      # down_train$occurrence <- factor(down_train$occurrence, levels = c("presence", "absence"))
      
      #---------------------------------------------------------------------------
      # Perform model
      
      # Model features
      ctrl <- trainControl(method = "cv", 
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                           savePredictions = TRUE,
                           seeds = 1:11)
      
      # Model
      glmFit <- train(occurrence ~ .,
                      data = down_train,
                      method = "glm",
                      trControl = ctrl,
                      metric = "ROC")
      
      #---------------------------------------------------------------------------
      # Performance on training
      
      training <- data.frame(obs = data_training$occurrence,
                             pred = predict(glmFit, 
                                            newdata = data_training),
                             presence = predict(glmFit, 
                                                newdata = data_training, 
                                                type = "prob")[,1],
                             absence = predict(glmFit, 
                                               newdata = data_training, 
                                               type = "prob")[,2])
      
      training_performance <- c(defaultSummary(training, 
                                               lev = levels(data_training$occurrence)),
                                twoClassSummary(training, 
                                                lev = levels(data_training$occurrence)),
                                prSummary(training, 
                                          lev = levels(data_training$occurrence)),
                                mnLogLoss(training, 
                                          lev = levels(data_testing$occurrence)))
      
      #---------------------------------------------------------------------------
      # Performance on testing
      
      testing <- data.frame(obs = data_testing$occurrence,
                            pred = predict(glmFit, newdata = data_testing),
                            presence = predict(glmFit, 
                                               newdata = data_testing, 
                                               type = "prob")[,1],
                            absence = predict(glmFit, 
                                              newdata = data_testing, 
                                              type = "prob")[,2])
      
      testing_performance <- c(defaultSummary(testing, 
                                              lev = levels(data_testing$occurrence)),
                               twoClassSummary(testing, 
                                               lev = levels(data_testing$occurrence)),
                               prSummary(testing, 
                                         lev = levels(data_testing$occurrence)),
                               mnLogLoss(testing, 
                                         lev = levels(data_testing$occurrence)))
      
      #---------------------------------------------------------------------------
      # Collect data
      training <- as.data.table(training)
      testing <- as.data.table(testing)
      
      # Coefficients
      collect <- rbind(c(training_performance, glmFit$finalModel$coefficients),
                       c(testing_performance, glmFit$finalModel$coefficients))
      collect <- as.data.table(collect)
      
      # varImp
      var <- rbind(as.matrix(t(varImp(glmFit$finalModel)), nrow = 1), as.matrix(t(varImp(glmFit$finalModel)), nrow = 1))
      colnames(var) <- paste0("imp_", colnames(var))
      collect <- cbind(collect, var)
      
      # Merge
      collect <- cbind(plant = c(group_name, group_name), 
                       iteration = c(ii, ii),
                       dataset = c("training", "testing"),
                       psamples = c(nrow(training[obs == "presence",])/nrow(training)*100,
                                    nrow(testing[obs == "presence",])/nrow(testing)*100),
                       collect)
      
      lineages_collector <- rbind(lineages_collector, collect)
      
    }
    
    return(lineages_collector)
    
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 2:groups, 
                         FUN = apply_occurrence, 
                         ordination = ordination,
                         aggregation = aggregation,
                         data_split = data_split,
                         iterations = iterations,
                         mc.preschedule = TRUE, 
                         mc.set.seed = FALSE,
                         mc.silent = FALSE, 
                         mc.cores = threads,
                         mc.cleanup = TRUE, 
                         mc.allow.recursive = TRUE, 
                         affinity.list = NULL)
  
  complete <- do.call(rbind, complete)
  
  return(complete)
  
}


split_lineages <- function(occurrence,
                         split = 0.85,
                         seed) {
  
  vec <- integer()
  
  split_frame <- data.table(sample = 1:length(occurrence), occurrence = occurrence)
  
  # Presence
  presence <- split_frame[occurrence == "presence",]
  set.seed(seed)
  train_presence <- sample(1:nrow(presence), floor(nrow(presence)*split))
  vec <- c(vec, presence$sample[train_presence])
  
  # Absence
  absence <- split_frame[occurrence == "absence",]
  set.seed(seed)
  train_absence <- sample(1:nrow(absence), floor(nrow(absence)*split))
  vec <- c(vec, absence$sample[train_absence])
  
  return(vec)
  
}
