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
    frame <- cbind(aggregation[, .SD, .SDcols = c(X)], ordination[, 2:4])
    colnames(frame) <- c("occurrence", "Axis1", "Axis2", "Axis3") 
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
      split <- split_fortyp(data_training$occurrence, split = 0.85)
      
      #---------------------------------------------------------------------------
      # Class balance
      down_train <- ROSE(occurrence ~ ., data = data_training[split,])$data
      down_train$occurrence <- factor(down_train$occurrence, levels = c("presence", "absence"))
      
      #---------------------------------------------------------------------------
      # Perform model
      
      # Model
      model <- glm(occurrence ~ Axis1 + Axis2 + Axis3, 
                   data = down_train, 
                   family = binomial(link = "logit"))
      
      #---------------------------------------------------------------------------
      # Performance on training
      
      train_pred <- ifelse(predict(model, data_training, type = "response") > 0.5, 
                           "presence", 
                           "absence")
      train_pred <- factor(train_pred, levels = c("presence", "absence"))
      
      training <- data.frame(obs = data_training$occurrence,
                             pred = train_pred,
                             presence = predict(model, data_training, type = "response"),
                             absence = 1 - predict(model, data_training, type = "response"))
      
      training_performance <- c(defaultSummary(training, 
                                               lev = levels(training$obs)),
                                twoClassSummary(training, 
                                                lev = levels(training$obs)),
                                prSummary(training, 
                                          lev = levels(training$obs)))
      
      #---------------------------------------------------------------------------
      # Performance on testing
      
      test_pred <- ifelse(predict(model, data_testing, type = "response") > 0.5, 
                          "presence", 
                          "absence")
      test_pred <- factor(test_pred, levels = c("presence", "absence"))
      
      testing <- data.frame(obs = data_testing$occurrence,
                            pred = test_pred,
                            presence = predict(model, data_testing, type = "response"),
                            absence = 1 - predict(model, data_testing, type = "response"))
      
      testing_performance <- c(defaultSummary(testing, 
                                              lev = levels(testing$obs)),
                               twoClassSummary(testing, 
                                               lev = levels(testing$obs)),
                               prSummary(testing, 
                                         lev = levels(testing$obs)))
      
      #---------------------------------------------------------------------------
      # Collect data
      training <- as.data.table(training)
      testing <- as.data.table(testing)
      
      collect <- rbind(c(training_performance, model$coefficients),
                       c(testing_performance, model$coefficients))
      collect <- as.data.table(collect)
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


split_fortyp <- function(occurrence,
                         split = 0.85) {
  
  vec <- integer()
  
  split_frame <- data.table(sample = 1:length(occurrence), occurrence = occurrence)
  
  # Presence
  presence <- split_frame[occurrence == "presence",]
  train_presence <- sample(1:nrow(presence), floor(nrow(presence)*split))
  vec <- c(vec, presence$sample[train_presence])
  
  # Absence
  absence <- split_frame[occurrence == "absence",]
  train_absence <- sample(1:nrow(absence), floor(nrow(absence)*split))
  vec <- c(vec, absence$sample[train_absence])
  
  return(vec)
  
}
