#-------------------------------------------------------------------------------
# Model Performance Function
model_performance <- function(meta_split,
                              ordination_split, 
                              spectra_split,
                              models,
                              tune_optimal,
                              threads) {
  
  #-----------------------------------------------------------------------------
  # Apply prediction function
  application_predict <- function(X,
                                  models,
                                  spectra_split,
                                  tune_optimal) {
    # Predict new values
    axis1_predicted <- predict(models[[X]]$axis1_model, 
                               newdata = spectra_split,
                               ncom = tune_optimal[1])[,1,]
    
    axis2_predicted <- predict(models[[X]]$axis2_model, 
                               newdata = spectra_split,
                               ncom = tune_optimal[2])[,1,]
    
    axis3_predicted <- predict(models[[X]]$axis3_model, 
                               newdata = spectra_split,
                               ncom = tune_optimal[3])[,1,]
    
    
    
    predicted <- cbind(axis1_predicted,
                       axis2_predicted,
                       axis3_predicted)
    predicted <- as.data.frame(predicted)
    colnames(predicted) <- c("Axis1", "Axis2", "Axis3")
    
    return(predicted)
    
  }
  
  # Predicted in parallel
  predicted <- pbmclapply(X = 1:length(models), 
                          FUN = application_predict, 
                          models = models,
                          spectra_split = spectra_split,
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
  ord <- as.data.frame(ordination_split)
  
  application_performance <- function(X,
                                      ord,
                                      predicted) {
    # Predict new values
    iteration <- predicted[[X]]
    collector <- data.table()
    
    # Performance at the Axis level
    for(i in 1:ncol(ordination_split)) {
      per <- parameters(obs = ord[, i], pred = iteration[,i])
      per <- cbind(Axes = paste0("Axis", i),
                   iteration = X,
                   per)
      collector <- rbind(collector, per)
    }
    
    # Performance at the Axes level
    per <- parameters_3d(observed = ordination_split, 
                         predicted_values = as.data.table(iteration))
    
    collector <- cbind(collector, per)
    
    return(collector)
    
  }
  
  # Performance in parallel
  performance <- pbmclapply(X = 1:length(models),
                            FUN = application_performance,
                            ord = ord,
                            predicted = predicted,
                            mc.preschedule = TRUE, 
                            mc.set.seed = FALSE,
                            mc.silent = FALSE, 
                            mc.cores = threads,
                            mc.cleanup = TRUE, 
                            mc.allow.recursive = TRUE, 
                            affinity.list = NULL)
  
  frame <- do.call(rbind, performance)
  
  #-----------------------------------------------------------------------------
  # Get residuals
  # Apply prediction function
  predicted <- array(unlist(predicted), dim = c(dim(predicted[[1]])[1],
                                                dim(predicted[[1]])[2],
                                                length(predicted)))
  
  application_residuals <- function(observed, predicted, PLT_CN = meta_split$PLT_CN) {
    
    collector <- data.table()
    
    for(i in 1:ncol(observed)) {
      
      axis <- as.matrix(observed)[, i]
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
      
      results <- data.table(PLT_CN = PLT_CN,
                            Axes = paste0("Axis", i),
                            observed = axis,
                            predicted = predictions_means,
                            predicted_sd = predictions_sd,
                            residuals = residuals_means,
                            residuals_sd = residuals_sd)
      
      collector <- rbind(collector, results)
      
    }
    
    return(collector)
    
  }
    
  predicted_summary <- application_residuals(ord, predicted)
  
  #-----------------------------------------------------------------------------
  # Get 3D residuals
  # Apply prediction function
  residuals_3D <- function(observed, predicted, PLT_CN = meta_split$PLT_CN) {
    
    collector <- data.table()
    
    for(i in 1:dim(predicted)[3]) {
      collector <- cbind(collector, rowSums((observed - predicted[,,i])^2))
    }
    
    residuals_means <- rowMeans(collector)
    residuals_sd <- apply(collector, 1, sd)
    
    results <- data.table(PLT_CN = PLT_CN,
                          residuals_3D = residuals_means,
                          residuals_3D_sd = residuals_sd)
    
  }
  
  residuals_3D_predicted <- residuals_3D(ord, predicted)
  
  #-----------------------------------------------------------------------------
  # Return results
  return(list(performance = frame,
              predicted = predicted_summary,
              residuals_3D = residuals_3D_predicted))
  
}

#-------------------------------------------------------------------------------
# Parameters of performance to estimate 1D
parameters <- function(obs, pred) {
  
  linear <- lm(obs ~ pred, na.action = na.exclude)
  linear_summary <- summary(linear)
  
  R2 <- 1 - sum((obs - pred)^2, na.rm = TRUE) / sum((obs - mean(obs, na.rm = TRUE))^2)
  intercept <- linear_summary$coefficients[1,1]
  slope <- linear_summary$coefficients[2,1]
  BIAS <- mean(obs-pred, na.rm = TRUE)
  RMSE <- sqrt(sum(((obs-pred)^2))/length(obs))
  perRMSE <- sqrt(sum(((obs-pred)^2))/length(obs))/(quantile(obs, 0.99) - quantile(obs, 0.01))
  names(perRMSE) <- NULL
  
  return(round(data.table(R2 = R2,
                          BIAS = BIAS,
                          RMSE = RMSE,
                          perRMSE = perRMSE*100,
                          intercept = intercept,
                          slope = slope), 4))
  
}

# Parameters of performance  estimate 3D
parameters_3d <- function(observed, predicted_values) {
  
  # Mean of observed vectors
  mean_observed <- colMeans(observed)
  
  # Residual Sum of Squares (RSS)
  rss <- sum(rowSums((observed - predicted_values)^2))
  
  # Total Sum of Squares (TSS)
  tss <- sum(rowSums((observed - mean_observed)^2))
  
  # R-squared
  r_squared <- 1 - (rss / tss)
  
  # RMSE (Root Mean Squared Error)
  n <- nrow(observed)
  rmse <- sqrt(sum(rowSums((observed - predicted_values)^2)) / n)
  
  # Overall bias magnitude (optional, as a scalar)
  bias <- colMeans(observed - predicted_values)
  bias_magnitude <- sqrt(sum(bias^2))
  
  return(data.table(R2_3D = r_squared,
                    RMSE_3D = rmse,
                    BIAS_3D = bias_magnitude))
  
}
