model_build <- function(meta, 
                        ordination, 
                        spectra, 
                        iscv,
                        split, 
                        tune_optimal, 
                        ntimes = 100, 
                        threads = 26) {
  
  # Creates partition for training and testing-------------
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  
  #Results collector
  collector <- list()
  
  # Parallel implementation
  if(threads > 1) {
    if(.Platform$OS.type == "unix") {
      pls.options(parallel = threads)
    } else {
      pls.options(parallel = quote(makeCluster(threads, type = "PSOCK")))
    }
  }
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = ntimes, style = 3) 
  
  # Loop of repeated
  for(t in 1:ntimes) {
    
    # Progress
    setTxtProgressBar(pb, t)
    
    # Model
    ml_model <- plsr(cbind(MDS1 = ordination_training$MDS1, 
                           MDS2 = ordination_training$MDS2, 
                           MDS3 = ordination_training$MDS3) ~ ., 
                     data = spectra_training, 
                     validation = "CV",
                     segments = iscv[[t]],
                     ncomp = max(tune_optimal),
                     center = TRUE,
                     method = "oscorespls",
                     jackknife = TRUE,
                     maxit = 1000)
    
    collector[[t]] <- ml_model
    
  }
  
  return(collector)
  
}
