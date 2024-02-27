#-------------------------------------------------------------------------------
# Model Performance Function
model_performance <- function(meta,
                              ordination, 
                              spectra,
                              models,
                              tune_optimal,
                              threads = 28) {
  
  # Functions for loops --------------------------------------------------------
  # Apply prediction
  apply_predict <- function(i,
                            model_t = models[[i]],
                            spectra = spectra,
                            tune_optimal) {
    # Predict new values
    predicted <- predict(model_t, ncomp = tune_optimal, newdata = spectra)
    predicted <- as.data.table(cbind(MDS1 = predicted[,1,1],
                                     MDS2 = predicted[,2,2],
                                     MDS3 = predicted[,3,3]))
    return(predicted)
    
  }
  
  # Apply metrics
  apply_performance <- function(i,
                                obs,
                                pred) {
    
    # Get metrics
    MDS1 <- parameters(obs$MDS1, pred$MDS1)
    MDS2 <- parameters(obs$MDS2, pred$MDS2)
    MDS3 <- parameters(obs$MDS3, pred$MDS3)
    MDS <- rbind(MDS1, MDS2, MDS3)
    results <- cbind(data.table(iteration = c(i, i, i),
                                Axes = c("MDS1", "MDS2", "MDS3")),
                     MDS)
    
    return(results)
    
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

parameters <- function(obs, pred) {
  
  R2 <- 1 - (sum((obs-pred)^2)/sum((obs-mean(obs))^2))
  MAE <- sum(abs(obs-pred))/length(obs)
  RMSE <- sqrt(sum((obs-pred)^2)/length(obs))
  perRMSE <- sqrt(sum((obs-pred)^2)/length(obs))/(quantile(obs, 0.99) - quantile(obs, 0.01))
  names(perRMSE) <- NULL
  B1 <- sqrt(sum((obs - mean(obs))^2)/length(obs)) / sqrt(sum((pred - mean(pred))^2)/length(pred))
  B0 <- mean(obs) - (B1*mean(pred))
  
  return(round(c(R2 = R2,
           MAE = MAE,
           RMSE = RMSE,
           perRMSE = perRMSE,
           B1 = B1,
           B0 = B0), 4))
  
}
