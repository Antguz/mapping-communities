#-------------------------------------------------------------------------------
# Model Predict Function
model_performance_multiple <- function(meta = meta,
                                       ordination = ordination, 
                                       spectra = spectra,
                                       models = models,
                                       ncomp = 19) {
  
  #Performance descriptors
  my.metrics <- c("R2", "MAE", "RMAE", "RMSE", "RRMSE")
  
  #Number of models
  nMDS <- length(models)
  
  #Collectors
  final_performance <- data.table()
  final_predicted <- data.table()
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = nMDS, style = 3) 
  
  for(i in 1:nMDS) {
    
    #Progress
    setTxtProgressBar(pb, i)
    
    mod <- models[[i]]
    predicted <- as.data.table(predict(mod,
                         newdata = as.matrix(spectra),
                         ncomp = ncomp,
                         type = "response")[,,1])
    
    
    for(q in 1:ncol(ordination)) {
      
      MDS <- paste0("MDS", q)
      
      metrics <- metrics_summary(obs = ordination[, ..MDS]$MDS, 
                                 pred = predicted[, ..MDS]$MDS, 
                                 type = "regression",
                                 metrics_list = my.metrics)
      
      names <- data.table(MDS = MDS, Iteration = i)
      metrics <- metrics$Score
      metrics <- as.data.table(matrix(metrics, nrow = 1, ncol = length(metrics)))
      metrics <- cbind(names, metrics)
      
      exp_predicted <- round(predicted[, ..MDS]$MDS, 4)
      exp_predicted <- as.data.table(matrix(exp_predicted, nrow = 1, ncol = length(exp_predicted)))
      exp_predicted <- cbind(names, exp_predicted)
      
      #Export
      final_performance <- rbind(final_performance, metrics)
      final_predicted <- rbind(final_predicted, exp_predicted)
      
    }
  }
  
  colnames(final_performance) <- c("MDS", "Iteration", "R2", "MAE", "RMAE", "RMSE", "RRMSE")
  colnames(final_predicted)[1:2] <- c("MDS", "Iteration")
  
  return(list(performance = final_performance,
              predicted = final_predicted))
  
}
