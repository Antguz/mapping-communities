pls_vip <- function(models,
                    tune_optimal,
                    threads = 5) {
  
  
  vip_application <- function(X,
                              models,
                              tune_optimal) {
    
    # Select model
    object <- models[[X]]
    
    # Get coef
    axis1_vip <- VIP(object$axis1_model, tune_optimal[1])
    axis2_vip <- VIP(object$axis2_model, tune_optimal[2])
    axis3_vip <- VIP(object$axis3_model, tune_optimal[3])
    
    # Merge
    vip_val <- rbind(axis1_vip, axis2_vip, axis3_vip)
    vip_val <- cbind(data.table(iteration = X, 
                             Axes = c("Axis1", "Axis2", "Axis3")), vip_val)
    colnames(vip_val) <- gsub("`", "", colnames(vip_val))
    
    return(vip_val)
    
  }
  
  # Coefficients
  vip <- pbmclapply(X = 1:length(models),
                    FUN = vip_application,
                    models = models,
                    tune_optimal = tune_optimal,
                    mc.set.seed = TRUE,
                    mc.cores = threads)
  
  vip <- do.call(rbind, vip)
  
  return(vip)
  
}
