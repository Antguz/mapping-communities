#-------------------------------------------------------------------------------
# Model Predict Function
model_performance <- function(meta = meta,
                              ordination = ordination, 
                              spectra = spectra,
                              models = models) {
  
  #Performance descriptors
  my.metrics <- c("R2", "MAE", "RMAE", "RMSE", "RRMSE")
  
  #Number of models
  nMDS <- length(models)
  
  #Collectors
  final_performance <- data.table()
  final_predicted <- data.table()
  
  for(q in 1:nMDS) {
    
    #Repeated
    repeated <- length(models[[q]])
    
    #Observed
    observed <- ordination[, .SD, .SDcols = q]
    colnames(observed)[1] <- "MDS"
    observed <- observed$MDS
    
    
    #Collectors
    result_performance <- data.table()
    result_predicted <- data.table()
    
    for(i in 1:repeated) {
      
      grp <- models[[q]][[i]]
      predicted <- predict(grp, spectra)
      
      metrics <- metrics_summary(obs = observed, 
                                 pred = predicted, 
                                 type = "regression",
                                 metrics_list = my.metrics)
      metrics <- data.table(MDS = paste0("MDS", q),
                            repeated = i,
                            R2 = metrics[1,2],
                            MAE = metrics[2,2],
                            RMAE = metrics[3,2],
                            RMSE = metrics[4,2],
                            RRMSE = metrics[5,2])
      result_performance <- rbind(result_performance, metrics)
      
      predicted <- c(paste0("MDS", q), i, predicted)
      predicted <- as.data.table(matrix(predicted, nrow = 1))
      colnames(predicted)[1:2] <- c("MDS", "repeated")
      result_predicted <- rbind(result_predicted, predicted)
      
    }
    
    final_performance <- rbind(final_performance, result_performance)
    final_predicted <- rbind(final_predicted, result_predicted)
    
  }
  
  return(list(performance = final_performance,
              predicted = final_predicted))
  
}
