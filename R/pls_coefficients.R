pls_coefficients <- function(models = final_models,
                             tune_optimal) {
  
  # Collector
  coefficients <- data.table()
  
  # Number of models
  iterations <- length(models)
  
  # Coefficients
  for(i in 1:iterations) {
      
      #Select model
      object <- models[[i]]
      
      #B <- rowSums(coef(object, comps = ncomp), dims = 2)
      #B0 <- object$Ymeans - object$Xmeans %*% B
      
      coef <- coef(object, ncomp = tune_optimal, intercept=TRUE)
      coef <- rbind(coef[,1,1], coef[,2,2], coef[,3,3])
      coef <- cbind(data.table(iteration = i, 
                               Axes = c("Axis1", "Axis2", "Axis3")), coef)
      
      coefficients <- rbind(coefficients, 
                            coef)
      
    }
  
  return(coefficients)
  
}
