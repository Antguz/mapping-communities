#-------------------------------------------------------------------------------
# Model Predict Function
model_predict <- function(meta = meta,
                          ordination = ordination, 
                          spectra = spectra,
                          models = models,
                          threads = 28) {
  
  # Functions for loops --------------------------------------------------------
  # Apply prediction
  apply_predict <- function(i,
                            model_q = models[[q]],
                            spectra = spectra) {
    # Predict new values
    predicted <- predict(model_q[[i]], newdata = spectra)
    return(predicted)
  }
  
  # Apply metrics
  apply_performance <- function(i,
                                observed,
                                predicted,
                                metrics = metrics_feature) {
    
    # Get metrics
    metrics <- metrics_summary(obs = observed, 
                               pred = as.numeric(as.matrix(predicted[, .SD, .SDcols = i])), 
                               type = "regression",
                               metrics_list = metrics,
                               orientation = "OP")
    
    ite_performance <- metrics$Score
    names(ite_performance) <- metrics$Metric
    return(ite_performance)
  }
  
  # Parameters for predicting --------------------------------------------------
  # Number of models
  nMDS <- length(models)
  
  #Performance descriptors
  metrics_feature <- c("R2", "MAE", "RMAE", "RMSE", "RRMSE", "B0", "B1", "MBE")
  
  # Collectors
  final_predicted <- data.table()
  final_performance <- data.table()
  
  # Loop over axes
  for(q in 1:nMDS) {
    
    # Repeated
    repeated <- length(models[[q]])
    
    # Get predictions in parallel ----------------------------------------------
    predictions <- mclapply(X = 1:repeated, 
                            FUN = apply_predict,
                            model_q = models[[q]],
                            spectra = spectra,
                            mc.cores = threads,
                            mc.preschedule = TRUE,
                            mc.cleanup = TRUE)
    
    # Arrange the observed and predicted values
    predictions <- data.table(do.call(cbind, predictions))  
    colnames(predictions) <- paste0("iteration_", 1:ncol(predictions))
    frame <- ordination[, .SD, .SDcols = q]
    colnames(frame)[1] <- "observed"
    frame <- cbind(meta, 
                   data.table(MDS = paste0("MDS", q),
                              MDS_observed = frame$observed,
                              MDS_predicted_mean = rowMeans(predictions),
                              MDS_predicted_sd = apply(predictions, 1, sd)),
                   predictions)
    
    # Export from loop
    final_predicted <- rbind(final_predicted, frame)
    
    # Get performance in parallel ----------------------------------------------
    performance <- mclapply(X = 1:repeated, 
                            FUN = apply_performance,
                            observed = frame$MDS_observed,
                            predicted = frame[, .SD, .SDcols = c((ncol(frame) - repeated + 1):ncol(frame))],
                            metrics = metrics_feature,
                            mc.cores = threads,
                            mc.preschedule = TRUE,
                            mc.cleanup = TRUE)
    
    # Arrange the performance
    performance <- data.table(do.call(rbind, performance))
    performance <- cbind(data.table(MDS = paste0("MDS", q),
                                    iteration = 1:nrow(performance)),
                         performance)
    
    # Export from loop
    final_performance <- rbind(final_performance, performance)
    
    # Remove elements
    rm(list = c("repeated", "predictions", "frame"))
    
  }
  
  # Export from the function
  return(list(performance = final_performance,
              predicted = final_predicted))
  
}




