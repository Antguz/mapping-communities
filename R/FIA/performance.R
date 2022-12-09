performance <- function(observed, predicted) {
  
  obser <- observed
  pred <- predicted
  
  rss <- sum((obser - pred) ^ 2)  ## residual sum of squares
  tss <- sum((obser - mean(obser))^2)  ## total sum of squares
  
  n <- length(predicted)
  Rsq <- 1 - (rss/tss)
  Bias <- mean(((obser - pred)/obser))
  RMSE <- sqrt(mean((obser - pred)^2))
  RMSE_P <- (RMSE/(max(obser)-min(obser)))*100
  
  return(round(c(N = n, Rsq = Rsq, Bias = Bias, RMSE = RMSE, RMSE_P = RMSE_P), 6))
  
}


