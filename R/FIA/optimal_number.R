optimal_number <- function(meta, ordination, spectra, split, int_ncomp = 100, repeated = 10, k = 10, threads = 16) {
  
  # Creates partition for training and testing
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  dimentions <- ncol(ordination_training) - 1

  #Results collector
  results <- data.table()
  
  #Progress bar
  pb <- txtProgressBar(min = 1, max = repeated*dimentions, style = 3) 
  
  #Parallel 
  if(threads > 1 ) {
    cl <- makeCluster(threads)
    pls.options(parallel = cl)
  }
  
  #Loop of repeated
  for(q in 1:dimentions) {
    
    #Train set
    train_master <- cbind(ordination_training[, .SD, .SDcols = c(q+1)],
                   spectra_training)
    colnames(train_master)[1] <- "MDS"
    
    #Loop over ordination dimentions
    for(i in 1:repeated) {
      
      setTxtProgressBar(pb, i)
      
      #Freedman-Diaconis rule
      bw <- 2 * IQR(train_master$MDS) / length(train_master$MDS)^(1/3)
      bw <- floor(bw)
      
      #Select by gradient
      select <- break_continous(train_master$MDS, n = 5)
      train <- train_master[select,]
      
      print(paste0("Dimension:", q, "; repeated", i, "; samples", nrow(train)))
      
      #Training
      pls_components <- plsr(MDS ~ ., 
                             data= train, 
                             scale = FALSE,
                             center = TRUE,
                             ncomp= int_ncomp,
                             validation = "LOO",
                             method = "oscorespls",
                             na.action = na.omit)
      
      #Get optimal component
      component <- selectNcomp(pls_components, "onesigma")
      PRESS_min <- which.min(as.vector(pls_components$validation$PRESS))
      CV <- which.min(as.vector(RMSEP(pls_components)$val[1,,]))
      #plot(as.vector(pls_components$validation$PRESS))
      #plot(pls_components, ncomp = CV)
      
      #Frame to collect
      prediction <- data.table(ordination = colnames(ordination_training)[q+1],
                               Repeated = i,
                               Onesigma = component,
                               PRESS = PRESS_min,
                               CV = CV,
                               samples = nrow(train))
      
      results <- rbind(results, prediction)
      
    }
  }
  
  if(threads > 1 ) {
    stopCluster(cl)
  }
  
  return(results)
  
}

break_continous <- function(axis, n = 20) {
  
  range_axis <- range(axis)
  breaks <- seq((range_axis[1]-0.001), (range_axis[2]+0.001), length.out = n + 1)
  axis_cut <- cut(axis, breaks)
  frame_cut <- data.table(cut = axis_cut, row = 1:length(axis))
  
  cut_categories <- unique(axis_cut)
  nmin <- min(table(axis_cut))
  
  vector_rows <- NA
  
  for(i in 1:length(cut_categories)) {
    
    sub_frame <- subset(frame_cut, cut == cut_categories[i])
    sample_row <- sample(sub_frame$row, nmin)
    vector_rows <- c(vector_rows, sample_row)
  
  }
  
  vector_rows <- vector_rows[!is.na(vector_rows)]
  
  return(vector_rows)
  
}
