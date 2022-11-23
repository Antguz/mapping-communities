optimal_number <- function(meta, NDMS, spectra, split, int_ncomp = 50, repeated = 100, k = 10, threads = 6) {
  
  # Creates partition for training and testing
  NDMS_training <- NDMS[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  dimentions <- ncol(NDMS_training)
  
  #Create spatial points
  spatial <- sf::st_as_sf(meta_training, coords = c("xcoord", "ycoord"))
  st_crs(spatial) <- "ESRI:102004"
  
  #Results collector
  results <- data.table()
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = repeated, style = 3) 
  
  #Parallel 
  if(threads > 1 ) {
    cl <- makeCluster(threads, type = "PSOCK")
    pls.options(parallel = cl)
  }
  
  #Loop of repeated
  for(i in 1:repeated) {
    
    setTxtProgressBar(pb, i)
    
    #Loop over NDMS dimentions
    for(q in 1:dimentions) {
      
      #Create spatial bloks
      #cv <- CV_spatialBlock(spatial, k = 10)
      
      #Train set
      train <- cbind(NDMS_training[, .SD, .SDcols = c(q)],
                     spectra_training)
      colnames(train)[1] <- "NDMS"
      
      #Training
      pls_components <- plsr(NDMS ~ ., 
                             data= train, 
                             scale = FALSE,
                             center = TRUE,
                             ncomp= int_ncomp,
                             segments = 10,
                             validation = "CV", 
                             method = "oscorespls")
      
      #Get optimal component
      component <- selectNcomp(pls_components, "onesigma")
      PRESS_min <- which.min(as.vector(pls_components$validation$PRESS))
      CV <- which.min(as.vector(RMSEP(pls_components)$val[1,,]))
      
      #Frame to collect
      prediction <- data.table(NDMS = colnames(NDMS_training)[q],
                               Repeated = i,
                               Onesigma = component,
                               PRESS = PRESS_min,
                               CV = CV)
      
      results <- rbind(results, prediction)
      
    }
  }
  
  if(threads > 1 ) {
    stopCluster(cl)
  }
  
  return(results)
  
}
