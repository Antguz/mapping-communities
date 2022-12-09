optimal_number <- function(meta, ordination, spectra, split, int_ncomp = 50, repeated = 100, k = 10, threads = 16) {
  
  # Creates partition for training and testing
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  dimentions <- ncol(ordination_training)
  
  #Create spatial points
  #spatial <- sf::st_as_sf(meta_training, coords = c("xcoord", "ycoord"))
  #st_crs(spatial) <- "ESRI:102004"
  
  #Results collector
  results <- data.table()
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = repeated*10, style = 3) 
  
  #Parallel 
  if(threads > 1 ) {
    cl <- makeCluster(threads)
    pls.options(parallel = cl)
  }
  
  #Loop of repeated
  for(i in 1:repeated) {
    
    setTxtProgressBar(pb, i)
    
    #Loop over ordination dimentions
    for(q in 1:10) {
      
      setTxtProgressBar(pb, i)
      
      #Create spatial bloks
      #cv <- CV_spatialBlock(spatial, k = 10)
      
      #Train set
      train <- cbind(ordination_training[, .SD, .SDcols = c(q)],
                     spectra_training)
      colnames(train)[1] <- "ordination"
      
      #Training
      pls_components <- plsr(ordination ~ ., 
                             data= train, 
                             scale = FALSE,
                             center = TRUE,
                             ncomp= int_ncomp,
                             segments = 10,
                             validation = "CV",
                             method = "oscorespls",
                             na.action = na.omit)
      
      #Get optimal component
      component <- selectNcomp(pls_components, "onesigma")
      PRESS_min <- which.min(as.vector(pls_components$validation$PRESS))
      CV <- which.min(as.vector(RMSEP(pls_components)$val[1,,]))
      
      #Frame to collect
      prediction <- data.table(ordination = colnames(ordination_training)[q],
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
