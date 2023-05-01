#-------------------------------------------------------------------------------
# Model Predict Function
model_performance <- function(meta = meta,
                              ordination = ordination, 
                              spectra = spectra,
                              models = models) {
  
  #Performance descriptors
  my.metrics <- c("R2", "MAE", "RMAE", "RMSE", "RRMSE", "B0", "B1", "MBE")
  
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
    
    #Print progress
    print(paste0("MDS", q))
    pb <- txtProgressBar(min = 0, max = repeated, style = 3)
    
    for(i in 1:repeated) {
      
      #Progress
      setTxtProgressBar(pb, i)
      
      #Select model
      grp <- models[[q]][[i]]
      
      #Predict model
      predicted <- predict(grp, newdata = spectra)
      
      #Get performance
      metrics <- metrics_summary(obs = observed, 
                                 pred = predicted, 
                                 type = "regression",
                                 metrics_list = my.metrics,
                                 orientation = "OP")
      metrics <- as.data.table(metrics)
      
      #Collect metrics
      metrics <- data.table(MDS = paste0("MDS", q),
                            repeated = i,
                            R2 = metrics[Metric == my.metrics[1], ]$Score,
                            MAE = metrics[Metric == my.metrics[2], ]$Score, 
                            RMAE = metrics[Metric == my.metrics[3], ]$Score,
                            RMSE = metrics[Metric == my.metrics[4], ]$Score,
                            RRMSE = metrics[Metric == my.metrics[5], ]$Score,
                            B0 = metrics[Metric == my.metrics[6], ]$Score,
                            B1 = metrics[Metric == my.metrics[7], ]$Score,
                            MBE = metrics[Metric == my.metrics[8], ]$Score)
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
