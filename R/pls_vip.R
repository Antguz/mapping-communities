pls_vip <- function(models = final_models,
                    tune_optimal = tune_optimal) {
  
  #Collector
  vip <- data.table()
  
  #Number of models
  iterations <- length(models)
  
  # VIP
  for(i in 1:iterations) {
    
    #Select model
    object <- models[[i]]
    
    #Get VIP
    vip_model <- VIP(object, tune_optimal)
    vip_model <- cbind(data.table(iteration = i,
                                  Axes = c("Axis1", "Axis2", "Axis3")), vip_model)
    vip <- rbind(vip, vip_model)
    
  }
  
  return(vip)
  
}

VIP <- function(object, tune_optimal) {
  
  complete <- data.table()
  
  for(i in 1:length(tune_optimal)) {
    
    W <- object$loading.weights[, 1:tune_optimal[i]]
    WW <- W * W/apply(W, 2, function(x) sum(x * x))
    Q <- object$Yloadings[i, 1:tune_optimal[i]]
    Q2 <- as.numeric(Q) * as.numeric(Q)
    TT <- object$scores[, 1:tune_optimal[i]]
    Q2TT <- Q2[1:tune_optimal[i]] * diag(crossprod(TT))[1:tune_optimal[i]]
    vip_axis <- sqrt(dim(object$coef)[1] * apply(sweep(WW[, 1:tune_optimal[i], drop = FALSE], 
                                2, Q2TT, "*"), 1, sum)/sum(Q2TT))
    complete <- rbind(complete, matrix(vip_axis, nrow = 1))
    
  }
  
  colnames(complete) <- names(object$coef[,1,1])
  return(complete)
  
}
