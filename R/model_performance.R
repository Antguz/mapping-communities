#-------------------------------------------------------------------------------
# Model Performance Function
model_performance <- function(meta,
                              ordination, 
                              spectra,
                              models,
                              tune_optimal,
                              threads) {
  
  # Secure I()
  ord <- as.matrix(ordination)
  model_frame <- data.frame(bands = I(as.matrix(spectra)))
  
  #-----------------------------------------------------------------------------
  # Apply prediction function
  application_predict <- function(X,
                                  models,
                                  spectra,
                                  tune_optimal) {
    # Predict new values
    predicted <- predict(models[[X]], 
                         newdata = model_frame,
                         ncom = tune_optimal)
    predicted <- cbind(predicted[,1,1],
                       predicted[,2,2],
                       predicted[,3,3])
    predicted <- as.data.frame(predicted)
    colnames(predicted) <- c("Axis1", "Axis2", "Axis3")
    
    return(predicted)
    
  }
  
  # Predicted in parallel
  predicted <- pbmclapply(X = 1:length(models), 
                          FUN = application_predict, 
                          models = models,
                          spectra = spectra,
                          tune_optimal = tune_optimal,
                          mc.preschedule = TRUE, 
                          mc.set.seed = FALSE,
                          mc.silent = FALSE, 
                          mc.cores = threads,
                          mc.cleanup = TRUE, 
                          mc.allow.recursive = TRUE, 
                          affinity.list = NULL)

  #-----------------------------------------------------------------------------
  # Evaluate performance
  application_performance <- function(X,
                                      ordination,
                                      predicted) {
    # Predict new values
    iteration <- predicted[[X]]
    collector <- data.table()
    
    for(i in 1:ncol(ordination)) {
      per <- parameters(obs = ord[, i], pred = iteration[,i])
      per <- cbind(Axes = paste0("Axis", i),
                   iteration = X,
                   per)
      collector <- rbind(collector, per)
    }
    
    return(collector)
  }
  
  # Performance in parallel
  performance <- mclapply(X = 1:length(models),
                          FUN = application_performance,
                          ordination = ordination,
                          predicted = predicted,
                          mc.preschedule = TRUE, 
                          mc.set.seed = FALSE,
                          mc.silent = FALSE, 
                          mc.cores = threads,
                          mc.cleanup = TRUE, 
                          mc.allow.recursive = TRUE, 
                          affinity.list = NULL)
  
  frame <- data.table()
  
  for(i in 1:length(performance)) {
    frame <- rbind(frame, performance[[i]], fill = TRUE)
  }
  
  #-----------------------------------------------------------------------------
  # Get residuals
  # Apply prediction function
  predicted <- array(unlist(predicted), dim = c(dim(predicted[[1]])[1],
                                                dim(predicted[[1]])[2],
                                                length(predicted)))
  
  application_residuals <- function(ordination, predicted) {
    
    collector <- data.table()
    
    for(i in 1:ncol(ordination)) {
      
      axis <- as.matrix(ordination)[, i]
      predictions <- predicted[,i,]
      
      # Predictions
      predictions_means <- rowMeans(predictions)
      predictions_sd <- apply(predictions, 1, sd)
      
      # Residuals
      residuals <- predictions - matrix(axis, 
                                        ncol = dim(predictions)[2],
                                        nrow = dim(predictions)[1])
      residuals_means <- rowMeans(residuals)
      residuals_sd <- apply(residuals, 1, sd)
      
      results <- data.table(Axes = paste0("Axis", i),
                            sample = 1:length(predictions_means),
                            observed = axis,
                            predicted = predictions_means,
                            predicted_sd = predictions_sd,
                            residuals = residuals_means,
                            residuals_sd = residuals_sd)
      
      collector <- rbind(collector, results)
      
    }
    
    return(collector)
    
  }
    
  predicted_summary <- application_residuals(ordination, predicted)
  
  #-----------------------------------------------------------------------------
  # Return results
  return(list(performance = frame,
              predicted = predicted_summary))
  
}

#-------------------------------------------------------------------------------
# Parameters to estimate
parameters <- function(obs, pred) {
  
  linear <- lm(obs ~ pred, na.action = na.exclude)
  linear_summary <- summary(linear)
  
  R2 <- linear_summary$r.squared
  intercept <- linear_summary$coefficients[1,1]
  slope <- linear_summary$coefficients[2,1]
  BIAS <- mean(obs-pred, na.rm = TRUE)
  RMSE <- sqrt(sum(((obs-pred)^2))/length(obs))
  perRMSE <- sqrt(sum(((obs-pred)^2))/length(obs))/(quantile(obs, 0.99) - quantile(obs, 0.01))
  names(perRMSE) <- NULL
  
  
  return(round(data.table(R2 = R2,
                          BIAS = BIAS,
                          RMSE = RMSE,
                          perRMSE = perRMSE,
                          intercept = intercept,
                          slope = slope), 4))
  
}
