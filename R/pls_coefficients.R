pls_coefficients <- function(models = final_models) {
  
  # Collector
  coefficients <- data.table()
  
  # Number of models
  nMDS <- length(models)
  
  # Coefficients
  for(i in 1:nMDS) {
      
      #Select model
      object <- models[[i]]
      ncomp <- object$ncomp
      
      #B <- rowSums(coef(object, comps = ncomp), dims = 2)
      #B0 <- object$Ymeans - object$Xmeans %*% B
      coef <- coef(object, ncomp = ncomp, intercept=TRUE)[,,1]
      coef <- cbind(data.table(iteration = i, bands = rownames(coef)), coef)
      
      coefficients <- rbind(coefficients, 
                            coef)
      
    }
  
  return(coefficients)
  
}
