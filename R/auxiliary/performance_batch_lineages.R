#-------------------------------------------------------------------------------
#' @Performance-function
#-------------------------------------------------------------------------------

#' @param predicted Axes predictions from spaceborne spectroscopy
#' @param coefficients Coefficients (and intercept) value from the occurance models

performance_batch <- function(predicted, coefficients, environment, scaling) {
  
  # Reshape predicted
  # Select columns
  obs <- predicted[, c("Axes", "PLT_CN", "observed")]
  pred <- predicted[, c("Axes", "PLT_CN", "predicted")]
  
  # Dcast columns
  obs <- dcast(obs, PLT_CN ~ Axes)
  pred <- dcast(pred, PLT_CN ~ Axes)
  
  # Merge with environment 
  obs <- merge(obs, environment, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
  obs <- na.exclude(obs)
  pred <- merge(pred, environment, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
  pred <- na.exclude(pred)
  
  # Scale values
  obs[, 2:ncol(obs)] <- (obs[, 2:ncol(obs)] - as.list(as.numeric(scaling[1]))) / as.list(as.numeric(scaling[2]))
  pred[, 2:ncol(pred)] <- (pred[, 2:ncol(pred)] - as.list(as.numeric(scaling[1]))) / as.list(as.numeric(scaling[2]))
  
  # Reshape coefficients
  models <- unique(coefficients$plant)
  coef_array <- array(dim = c(100, ncol(coefficients)-2, length(models)))
  intercept_mat <- matrix(nrow = 100, ncol = length(models))
  
  for(i in 1:length(models)) {
    
    sub <- coefficients[plant == models[i]]
    coef_array[,,i] <- as.matrix(sub[, .SD, .SDcols = c(3:ncol(sub))])
    intercept_mat[,i] <- sub$`(Intercept)`
    
  }
  
  # Predict observed and predicted
  observed <- cube_prob_coeff_prediction(as.matrix(obs[,-1]),
                                         coef_array,
                                         intercept_mat)
  
  predicted <- cube_prob_coeff_prediction(as.matrix(pred[,-1]),
                                          coef_array,
                                          intercept_mat)
  
  # Get dimensions
  dimensions <- dim(observed)
  
  # Performance frame
  performance <- data.table()
  
  # Estimate performance
  for(i in 1:dimensions[3]) {
    for(ii in 1:dimensions[2]) {
      
      # Predict new values
      per <- parameters(obs = observed[, ii, i], 
                        pred = predicted[, ii, i])
      
      # Compile results
      res <- data.table(plant = models[i],
                        iteration = ii,
                        per)
      performance <- rbind(performance, res)
      
    }
  }
  
  return(performance)
  
}