model_fortyp <- function(fia_meta,
                         forest_type,
                         ordination,
                         data_split,
                         iterations = 100,
                         threads = 25) {
  
  # Get number of plant groups
  groups <- ncol(data_split)
  
  # Get model application in batch
  apply_occurrence <- function(X,
                               fia_meta,
                               forest_type,
                               ordination,
                               iterations = 100,
                               data_split) {
    
    #---------------------------------------------------------------------------
    # Match databases
    if(all.equal(forest_type$PLT_CN, ordination$PLT_CN) == FALSE) {
      stop("Files most match")
    }
    
    # Group name for future reference
    group_name <- colnames(data_split)[X]
    ftype <- forest_type
    ftype$occurance <- "absence"
    ftype[FORTYP == group_name, occurance := "presence"]
    
    if(nrow(ftype[occurance == "presence",]) < 30) {
      
      return(data.table(plant = group_name,
                        iteration = NA,
                        dataset = NA,
                        psamples = NA,
                        Accuracy = NA,
                        Kappa = NA,
                        ROC = NA, 
                        Sens = NA, 
                        Spec = NA, 
                        AUC = NA, 
                        Precision = NA,
                        Recall = NA,
                        'F' = NA,
                        logLoss = NA,
                        '(Intercept)' = NA,
                        Axis1 = NA,
                        Axis2 = NA,
                        Axis3 = NA))
      
    } else {
      
      #---------------------------------------------------------------------------
      # Reshape data
      frame <- cbind(ftype[, .SD, .SDcols = 3], ordination[, 2:4])
      colnames(frame) <- c("occurrence", "Axis1", "Axis2", "Axis3") 
      frame$occurrence <- as.factor(frame$occurrence)
      frame$occurrence <- factor(frame$occurrence, levels = c("presence", "absence"))
      
      #---------------------------------------------------------------------------
      # Dealing with class imbalance
      
      data_training <- frame[as.matrix(data_split)[, X] == 1,]
      data_testing <- frame[as.matrix(data_split)[, X] == 0,]
      fortyp_collector <- data.table()
      
      for(ii in 1:iterations) {
        
        #---------------------------------------------------------------------------
        # Down sampling
        down_train <- ROSE(occurrence ~ ., data = data_training)$data
        down_train$occurrence <- factor(down_train$occurrence, levels = c("presence", "absence"))
        
        #---------------------------------------------------------------------------
        # Perform model
        
        # Model features
        ctrl <- trainControl(method = "cv", 
                             classProbs = TRUE,
                             summaryFunction = twoClassSummary,
                             savePredictions = TRUE)
        
        # Model
        plsFit <- train(occurrence ~ .,
                        data = down_train,
                        method = "glm",
                        trControl = ctrl,
                        metric = "ROC")
        
        #---------------------------------------------------------------------------
        # Performance on training
        
        training <- data.frame(obs = data_training$occurrence,
                               pred = predict(plsFit, 
                                              newdata = data_training),
                               presence = predict(plsFit, 
                                                  newdata = data_training, 
                                                  type = "prob")[,1],
                               absence = predict(plsFit, 
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
                              pred = predict(plsFit, newdata = data_testing),
                              presence = predict(plsFit, 
                                                 newdata = data_testing, 
                                                 type = "prob")[,1],
                              absence = predict(plsFit, 
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
        #training[obs == "presence",]
        #testing
        
        collect <- rbind(c(training_performance, plsFit$finalModel$coefficients),
                         c(testing_performance, plsFit$finalModel$coefficients))
        collect <- as.data.table(collect)
        collect <- cbind(plant = c(group_name, group_name), 
                         iteration = c(ii, ii),
                         dataset = c("training", "testing"),
                         psamples = c(nrow(training[obs == "presence",])/nrow(training)*100,
                                      nrow(testing[obs == "presence",])/nrow(testing)*100),
                         collect)
        
        fortyp_collector <- rbind(fortyp_collector, collect)
        
      }
      
      return(fortyp_collector)
      
    }
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 1:groups, 
                         FUN = apply_occurrence, 
                         fia_meta= fia_meta,
                         forest_type = forest_type,
                         ordination = ordination,
                         data_split = data_split,
                         iterations = iterations, 
                         mc.preschedule = TRUE, 
                         mc.set.seed = FALSE,
                         mc.silent = FALSE, 
                         mc.cores = threads,
                         mc.cleanup = TRUE, 
                         mc.allow.recursive = TRUE, 
                         affinity.list = NULL)
  
  complete <- do.call(rbind.data.frame, complete)
  complete <- complete[!is.na(dataset),]
  
  return(complete)
  
}