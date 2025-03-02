model_fortyp <- function(forest_type,
                         ordination,
                         data_split,
                         iterations = 100,
                         threads = 25) {
  
  # Get number of plant groups
  groups <- ncol(data_split)
  
  # Get model application in batch
  apply_occurrence <- function(X,
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
    
    if(group_name == "NA" | nrow(ftype[occurance == "presence",]) < 50) {
      
      return(data.table())
      
    } else {
      
      #---------------------------------------------------------------------------
      # Reshape data
      frame <- cbind(ftype[, .SD, .SDcols = 3], ordination[, .SD, .SDcols = c(2:ncol(ordination))])
      colnames(frame) <- c("occurrence", colnames(ordination[, .SD, .SDcols = c(2:ncol(ordination))])) 
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
        split <- split_fortyp(data_training$occurrence, split = 0.85, seed = ii)
        
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
        collect <- cbind(forty = c(group_name, group_name), 
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

  complete <- do.call(rbind, complete)
  complete <- complete[!is.na(dataset),]
  
  return(complete)
  
}


split_fortyp <- function(occurrence,
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
