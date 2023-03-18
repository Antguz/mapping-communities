model_performance <- function(meta, ordination, spectra, split, model) {
  
  # Creates partition for training and testing
  ordination_testing <- ordination[!split,]
  spectra_testing <- spectra[!split,]
  meta_testing <- meta[!split,]
  
  #Names
  predictions_names <- c("MDS", "observed", paste0("rep_", 1:(nrow(model$Performance)/3)))
  metrics <- c("R2","RSS", "RMSE", "RRMSE", "MAE", "RMAE", "RSE", "MBE", "B1", "B0")
  
  #Coefficients
  coefficients <- model$Coefficients
  dimentions <- unique(coefficients$MDS)
  
  #Collectors
  predictions <- data.table()
  perf <- data.table()
  names_metrics <- NA
  
  #Loop of repeated
  for(q in dimentions) {
    
    #Subset coefficients
    mds_coef <- subset(coefficients, MDS == dimentions[q])
    
    #Main coef
    main_coef <- as.matrix(mds_coef[, .SD, .SDcols = c(4:ncol(mds_coef))])
    
    #Perform prediction
    for(i in 1:nrow(spectra_testing)) {
      
      sample <- as.matrix(spectra_testing[i, ])
      
      #Matrix multiplication
      value <- matrix(rep(sample, nrow(main_coef)), nrow = nrow(main_coef), byrow = TRUE)
      value <- value * main_coef
      value <- rowSums(value)
      value <- value + mds_coef$intercept
      
      pred <- c(q, as.numeric(ordination_testing[i, .SD, .SDcols = c(q)]), value)
      predictions <- rbind(predictions, as.data.table(matrix(pred, nrow = 1)))
      
    }
  }
  
  colnames(predictions) <- predictions_names
    
  #Loop of repeated
  for(q in dimentions) {
    
    #Subset coefficients
    mds_predictions <- subset(predictions, MDS == dimentions[q])
    
    #True value
    observed <- mds_predictions$observed
    
    #Ncol
    rep <- ncol(predictions) - 2
    
    #Perform prediction
    for(i in 1:rep) {
      
      predicted <- as.vector(mds_predictions[, .SD, .SDcols = c(i+2)])[[1]]
      rep_perf <- metrics_summary(obs = observed, 
                                  pred = predicted,  
                                  type = "regression",
                                  metrics_list = metrics)
      score <- rep_perf$Score
      
      perf <- rbind(perf, as.data.table(matrix(c(q, i, score), nrow = 1)))
      
      if(i == 1) {
        names_metrics <- rep_perf$Metric
      }
      
    }
  }
  
  colnames(perf) <- c("MDS", "repeated", names_metrics)
  
  #Export
  return(list(Predicted = predictions, 
              Performance = perf))
  
}
