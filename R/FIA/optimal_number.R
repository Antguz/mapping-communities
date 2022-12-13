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
      
      #Select by gradient
      select <- break_continous(train$ordination, n = 20)
      train <- train[select, ]
      hist(train$ordination)
      dim(train)
      
      #Training
      pls_components <- plsr(ordination ~ ., 
                             data= train, 
                             scale = FALSE,
                             center = FALSE,
                             ncomp= 20,
                             validation = "LOO",
                             method = "oscorespls",
                             na.action = na.omit)
      
      saveRDS(pls_components, "pls3.rds")
      
      #Get optimal component
      component <- selectNcomp(pls_components, "onesigma")
      PRESS_min <- which.min(as.vector(pls_components$validation$PRESS))
      plot(as.vector(pls_components$validation$PRESS))
      plot(pls_components, ncomp = 8)
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

break_continous <- function(axis, n = 10) {
  
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

pls_components <- readRDS("pls2.rds")

ordination_testing <- ordination[!split,]
spectra_testing <- spectra[!split,]
meta_testing <- meta[!split,]
dimentions <- ncol(ordination_testing)
test <- cbind(ordination_testing[, .SD, .SDcols = c(q)],
              spectra_testing)
colnames(test)[1] <- "ordination"

plots <- fread("/media/antonio/Extreme SSD/DESIS/00_FIA/plots_WI_MN_20221102.csv")

predicted <- predict(pls_components, newdata = test, ncomp = 8)
observed <- test$ordination

frame <- data.table(predicted, observed, PLT_CN = meta_testing$PLT_CN)
frame <- merge(frame, plots, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
proportions <- fread("/media/antonio/antonio4/Github/mapping-species/proportions.csv")
frame <- merge(frame, proportions, by = "PLT_CN", all.x = TRUE, all.y = FALSE)

plot(predicted, observed)

jpeg("PLSR_PCoA1_M.jpg", width = 150, height = 150, units = "mm", res = 300)
ggplot(frame) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_point(aes(x = predicted, y = observed, fill = AM_to_T), shape = 21, colour = "grey", size = 2) +
  geom_smooth(aes(x = predicted, y = observed), method=lm, colour = "black", se = FALSE) +
  coord_fixed(xlim = c(-200, 400), ylim = c(-200, 400)) +
  theme_bw(base_size = 12) + ylab("Observed PCoA1") + xlab("Predicted PCoA1") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn")
dev.off()


jpeg("PLSR_PCoA2_M.jpg", width = 150, height = 150, units = "mm", res = 300)
ggplot(frame) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_point(aes(x = predicted, y = observed, fill = AM_to_T), shape = 21, colour = "grey", size = 2) +
  geom_smooth(aes(x = predicted, y = observed), method=lm, colour = "black", se = FALSE) +
  coord_fixed(xlim = c(-150, 150), ylim = c(-150, 150)) +
  theme_bw(base_size = 12) + ylab("Observed PCoA2") + xlab("Predicted PCoA2") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn")
dev.off()

jpeg("PLSR_PCoA3_M.jpg", width = 150, height = 150, units = "mm", res = 300)
ggplot(frame) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_point(aes(x = predicted, y = observed, fill = AM_to_T), shape = 21, colour = "grey", size = 2) +
  geom_smooth(aes(x = predicted, y = observed), method=lm, colour = "black", se = FALSE) +
  coord_fixed(xlim = c(-250, 150), ylim = c(-250, 150)) +
  theme_bw(base_size = 12) + ylab("Observed PCoA3") + xlab("Predicted PCoA3") +
  scale_fill_distiller("Proportion of trees with AM", palette = "PRGn")
dev.off()
