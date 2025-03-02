pls_coefficients <- function(models,
                             tune_optimal,
                             threads = 10) {
  
  coeff_application <- function(X,
                                models,
                                tune_optimal) {
    
    # Select model
    object <- models[[X]]
    
    # Get coef
    axis1_coef <- coef(object$axis1_model, ncomp = tune_optimal[1], intercept=TRUE)[,1,]
    axis2_coef <- coef(object$axis2_model, ncomp = tune_optimal[2], intercept=TRUE)[,1,]
    axis3_coef <- coef(object$axis3_model, ncomp = tune_optimal[3], intercept=TRUE)[,1,]
    
    # Merge
    coef <- rbind(axis1_coef, axis2_coef, axis3_coef)
    coef <- cbind(data.table(iteration = X, 
                             Axes = c("Axis1", "Axis2", "Axis3")), coef)
    colnames(coef) <- gsub("`", "", colnames(coef))
    
    return(coef)
    
  }
  
  # Coefficients
  coefficients <- pbmclapply(X = 1:length(models),
                             FUN = coeff_application,
                             models = models,
                             tune_optimal = tune_optimal,
                             mc.set.seed = TRUE,
                             mc.cores = threads)
  
  coefficients <- do.call(rbind, coefficients)
  
  return(coefficients)
  
}
