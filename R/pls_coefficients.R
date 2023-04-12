pls_coefficients <- function(models = final_models) {
  
  #Collector
  coefficients <- data.table()
  
  #Number of models
  nMDS <- length(models)
  
  #Name of bands
  bands_names <- numeric()
  
  for(q in 1:nMDS) {
    
    #Repeated
    repeated <- length(models[[q]])
    
    #Print progress
    print(paste0("MDS", q))
    pb <- txtProgressBar(min = 0, max = repeated, style = 3)
    
    for(i in 1:repeated) {
      
      #Progress
      setTxtProgressBar(pb, i)
      
      #Select model
      object <- models[[q]][[i]]
      object <- object$finalModel
      ncomp <- object$tuneValue$ncomp
      
      if(q == 1 & i == 1) {
        bands_names <- object$xNames
      }
      
      #B <- rowSums(coef(object, comps = ncomp), dims = 2)
      #B0 <- object$Ymeans - object$Xmeans %*% B
      coef <- c(q, i, as.numeric(coef(object, ncomp = ncomp, intercept=TRUE)[,,1]))
      
      coefficients <- rbind(coefficients, 
                            matrix(coef, nrow = 1, ncol = length(coef)))
      
    }
  }
  
  # Provide names
  colnames(coefficients)[1:3] <- c("MDS", "iteration", "Intercept")
  colnames(coefficients)[4:ncol(coefficients)] <- bands_names
  return(coefficients)
  
}
