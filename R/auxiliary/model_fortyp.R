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
      
      return(data.table(forty = group_name,
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
                        '(Intercept)' = NA,
                        Axis1 = NA,
                        Axis2 = NA,
                        Axis3 = NA))
      
    } else {
      
      #---------------------------------------------------------------------------
      # Reshape data
      frame <- cbind(ftype[, .SD, .SDcols = 3], ordination[, .SD, .SDcols = c(2:4)])
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
        split <- split_fortyp(data_training$occurrence, split = 0.85)
        
        #---------------------------------------------------------------------------
        # Class balance
        down_train <- ROSE(occurrence ~ ., data = data_training[split,])$data
        down_train$occurrence <- factor(down_train$occurrence, levels = c("presence", "absence"))
        
        #---------------------------------------------------------------------------
        # Perform model
        
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
