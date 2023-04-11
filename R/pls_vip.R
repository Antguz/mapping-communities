pls_vip <- function(models = final_models) {
  
  #Collector
  vip <- data.table()
  
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
      
      #Get VIP
      vip_model <- VIP(object, opt.comp = ncomp)
      vip_model <- c(q, i, vip_model)
      vip <- rbind(vip, matrix(vip_model, nrow = 1, ncol = length(vip_model)))
      
    }
  }
  
  # Provide names
  colnames(vip)[1:2] <- c("MDS", "iteration")
  colnames(vip)[3:ncol(vip)] <- bands_names
  return(vip)
  
}