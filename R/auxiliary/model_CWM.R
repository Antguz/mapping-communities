model_CWM <- function(fia_meta,
                      ordination,
                      aggregation,
                      variables, 
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
                               variables, 
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
    frame <- cbind(aggregation[, .SD, .SDcols = c(X)], ordination[, .SD, .SDcols = c(2:ncol(ordination))])
    colnames(frame) <- c("trait", colnames(ordination[, .SD, .SDcols = c(2:ncol(ordination))]))
    
    #---------------------------------------------------------------------------
    # Multiple iterations
    
    data_training <- frame[as.matrix(data_split)[, c(X-1)] == 1,]
    data_testing <- frame[as.matrix(data_split)[, c(X-1)] == 0,]
    cwm_collector <- data.table()
    
    #---------------------------------------------------------------------------
    # Define variables
    
    if(variables == "axes") {
      coi <- c("trait", "Axis1", "Axis2", "Axis3")
      data_training <- data_training[, ..coi] 
      data_testing <- data_testing[, ..coi]
    
    } else if(variables == "climate-topography") {
      coi <- c("trait", "MAT", "TS", "TAR", "MAP", "PS", "DEM", "SL", "CTI")
      data_training <- data_training[, ..coi] 
      data_testing <- data_testing[, ..coi]
      
    } else if(variables == "altogether") {
      coi <- c("trait", "Axis1", "Axis2", "Axis3", "MAT", "TS", "TAR", "MAP", "PS", "DEM", "SL", "CTI")
      data_training <- data_training[, ..coi] 
      data_testing <- data_testing[, ..coi]
    }
    
    for(ii in 1:iterations) {
      
      #-------------------------------------------------------------------------
      # Down sampling
      sub_split <- cwm_subsplit(trait = data_training$trait, 
                                split = 0.85,
                                seed = ii)
      down_train <- data_training[sub_split,]
      
      #---------------------------------------------------------------------------
      # Perform model
      
      # Model features
      ctrl <- trainControl(method = "cv", 
                           number = 10,
                           savePredictions = TRUE,
                           seeds = 1:11)
      
      # Model
      set.seed(ii)
      glmFit <- train(trait ~ .,
                      data = down_train,
                      method = "glm",
                      trControl = ctrl)
      
      #---------------------------------------------------------------------------
      # Performance on training
      
      training <- data.frame(obs = data_training$trait,
                             pred = predict(glmFit, 
                                            newdata = data_training))
      
      training_performance <- parameters(training$obs, training$pred)
      
      #---------------------------------------------------------------------------
      # Performance on testing
      
      testing <- data.frame(obs = data_testing$trait,
                            pred = predict(glmFit, newdata = data_testing))
      
      testing_performance <- parameters(testing$obs, testing$pred)
      
      #---------------------------------------------------------------------------
      # Collect data
      # Coefficients
      collect <- rbind(training_performance, testing_performance)
      collect <- cbind(collect, matrix(glmFit$finalModel$coefficients, nrow = 1))
      colnames(collect)[7:ncol(collect)] <- names(glmFit$finalModel$coefficients)
      
      # varImp
      var <- rbind(as.matrix(t(varImp(glmFit$finalModel)), nrow = 1), as.matrix(t(varImp(glmFit$finalModel)), nrow = 1))
      colnames(var) <- paste0("imp_", colnames(var))
      collect <- cbind(collect, var)
      
      # Collect
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
                         variables = variables,
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
  
  return(complete)
  
}

cwm_subsplit <- function(trait, 
                         split = 0.85,
                         seed) {
  
  test <- data.table(sample = 1:length(trait),
                     trait = trait)
  test$cut <- cut(test$trait, breaks = 10)
  test$train <- 0
  cuts_unite <- unique(test$cut)
  all <- data.frame()
  
  for(k in 1:length(cuts_unite)) {
    
    sub <- subset(test, cut == cuts_unite[k])
    set.seed(seed)
    sample_test <- sample(1:nrow(sub), ceiling(nrow(sub)*split))
    sub[sample_test, train := 1]
    all <- rbind(all, sub)
  
  }
  
  all <- all[order(sample),]
  all <- all[train == 1,]
  
  return(all$sample)
  
}
