pls_vip <- function(models = final_models) {
  
  #Collector
  vip <- data.table()
  bands_names <- character()
  
  #Number of models
  nMDS <- length(models)
  
  # VIP
  for(i in 1:nMDS) {
    
    #Select model
    object <- models[[i]]
    ncomp <- object$ncomp
    
    #Get VIP
    vip_model <- plsVarSel::VIP(object, opt.comp = ncomp)
    vip_model <- c(iteration = i, vip_model)
    bands_names <- names(vip_model)
    vip <- rbind(vip, matrix(vip_model, nrow = 1, ncol = length(vip_model)))
    
  }
  
  # Provide names
  colnames(vip) <- bands_names
  return(vip)
  
}
