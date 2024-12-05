#-------------------------------------------------------------------------------
#' @Performance-function
#-------------------------------------------------------------------------------

#' @param predicted Axes predictions from spaceborne spectroscopy
#' @param coefficients Coefficients (and intercept) value from the CWM models

performance_batch <- function(predicted, coefficients) {
  
  # Reshape predicted
  # Select columns
  obs <- predicted[, c("Axes", "sample", "observed")]
  pred <- predicted[, c("Axes", "sample", "predicted")]
  
  # Dcast columns
  obs <- dcast(obs, sample ~ Axes)
  pred <- dcast(pred, sample ~ Axes)
  
  # Reshape coefficients
  models <- unique(coefficients$plant)
  coef_array <- array(dim = c(100, 3, length(models)))
  intercept_mat <- matrix(nrow = 100, ncol = length(models))
  
  for(i in 1:length(models)) {
    
    sub <- coefficients[plant == models[i]]
    coef_array[,,i] <- as.matrix(sub[, 3:5])
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