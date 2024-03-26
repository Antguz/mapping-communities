model_occurence <- function(fia_meta,
                            ordination,
                            aggregation,
                            split,
                            knndm_folds,
                            threads = 4) {
  
  # Get number of plant groups
  groups <- ncol(aggregation) - 1
  
  # Get model application in batch
  apply_occurrence <- function(X,
                               fia_meta,
                               ordination,
                               aggregation,
                               split,
                               knndm_folds) {
                               #out_path) {
    
    # Group name fir future reference
    group_name <- colnames(aggregation)[X]
    
    #---------------------------------------------------------------------------
    # Match databases
    if(all.equal(aggregation$PLT_CN, ordination$PLT_CN) == FALSE) {
      stop("Files most match")
    }
    
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
    ctrl <- trainControl(method = "cv", 
                         classProbs = TRUE,
                         summaryFunction = twoClassSummary,
                         savePredictions = TRUE,
                         index = knndm_folds$indx_train)
    
    # Model
    plsFit <- train(occurrence ~ .,
                    data = frame[split,],
                    method = "bayesglm",
                    trControl = ctrl,
                    metric = "ROC")
    
    #---------------------------------------------------------------------------
    # Performance on training
    
    training <- data.frame(obs = frame$occurrence[split],
                           pred = predict(plsFit, 
                                          newdata = frame[split,]),
                           presence = predict(plsFit, 
                                              newdata = frame[split,], 
                                              type = "prob")[,2],
                           absence = predict(plsFit, 
                                             newdata = frame[split,], 
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
    
    testing <- data.frame(obs = frame[!split, ]$occurrence,
                          pred = predict(plsFit, newdata = frame[!split,]),
                          presence = predict(plsFit, 
                                             newdata = frame[!split,], 
                                             type = "prob")[,2],
                          absence = predict(plsFit, 
                                            newdata = frame[!split,], 
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
    # Create plots
    #p <- ggplot(plsFit$pred, 
    #       aes(m = presence, d = factor(obs, levels = c("absence", "presence")))) + 
    #  geom_roc(hjust = -0.4, vjust = 1.5, n.cuts = 10) + coord_equal() + style_roc()
      
    #p <- p + annotate("text", x= 0.05, y = 1.0, label=paste("AUC =", round((calc_auc(p))$AUC, 4)))
    #p <- p + labs(subtitle = group_name)
    
    #---------------------------------------------------------------------------
    # Collect data
    training <- as.data.table(training)
    testing <- as.data.table(testing)
    training[obs == "presence",]
    testing
    
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
                         split = split,
                         knndm_folds = knndm_folds,
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
