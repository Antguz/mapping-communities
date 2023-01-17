pls_ordination_model <- function(meta, ordination, spectra, split, repeated = 100, opt_components, threads = 1) {
  
  # Creates partition for training and testing
  ordination_training <- ordination[split,]
  spectra_training <- spectra[split,]
  meta_training <- meta[split,]
  dimentions <- ncol(ordination_training) - 1
  
  #Names
  performance_names <- c("MDS", "repeated", "Model", "N", "Rsq", "Bias", "RMSE", "RMSE_P")
  coefficients_names <- c("MDS", "repeated", "intercept", colnames(spectra_training))
  VIP_names <- c("MDS", "repeated", colnames(spectra_training))
  
  # Progress bar
  pb <- txtProgressBar(min = 1, max = repeated*dimentions, style = 3) 
  
  #Collectors
  VIP_export <- data.table()
  coeff_export <- data.table()
  perf_export <- data.table()
  
  # Parallel
  if(threads > 1) {
    cl <- makeCluster(threads)
    doSNOW::registerDoSNOW(cl)
  }
  
  #Loop of repeated
  for(q in 1:dimentions) {
    
    #Train set
    train_master <- cbind(ordination_training[, .SD, .SDcols = c(q+1)],
                          spectra_training)
    colnames(train_master)[1] <- "MDS"
    
    #Loop over ordination dimentions
    for(i in 1:repeated) {
      
      setTxtProgressBar(pb, i*q)
      
      #Basic information
      info <- data.table(MDS = q,
                         repeate = 1)
      
      #Freedman-Diaconis rule
      bw <- 2 * IQR(train_master$MDS) / length(train_master$MDS)^(1/3)
      bw <- floor(bw)
      
      #Select by gradient
      select <- break_continous(train_master$MDS, n = 10)
      train <- train_master[select,]
      
      #Training
      pls_components <- plsr(MDS ~ ., 
                             data= train, 
                             scale = FALSE,
                             center = TRUE,
                             ncomp= opt_components[q],
                             validation = "none", 
                             method = "oscorespls")
      #Performance
      predict_training <- predict(pls_components, 
                                  train_master[, 2:ncol(train_master)], 
                                  ncomp = opt_components[q], 
                                  type= "response")[,,1]
      
      perf_training <- c("Training", performance(train_master$MDS, predict_training))
      perf_training <- as.data.table(matrix(perf_training, nrow = 1))
      perf_training <- cbind(info, perf_training)
      colnames(perf_training) <- performance_names
      
      #Coefficients
      coefficients <- as.vector(coef(pls_components, 
                                     ncomp = opt_components[q], 
                                     intercept=TRUE))
      coefficients <- as.data.table(matrix(coefficients, nrow = 1))
      coefficients <- cbind(info, coefficients)
      colnames(coefficients) <- coefficients_names
      
      #VIP
      vip <- as.vector(VIP(pls_components, opt.comp = opt_components[q]))
      vip <- as.data.table(matrix(vip, nrow = 1))
      vip <- cbind(info, vip)
      colnames(vip) <- VIP_names
      
      #Collect
      perf_export <- rbind(perf_export, perf_training)
      coeff_export <- rbind(coeff_export, coefficients)
      VIP_export <- rbind(VIP_export, vip)
      
    }
  }
  
  if(threads > 1) {
    close(pb)
    stopCluster(cl)
  }
  
  #Export
  return(list(Coefficients = coeff_export, 
              VIP = VIP_export, 
              Performance = perf_export))
    
}

#Break contium-----------------------------------
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

