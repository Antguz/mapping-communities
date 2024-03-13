PCoA_segmentation <- function(ordination, 
                              nsegments = 50, 
                              lower.p = 0.01, 
                              proportion = 0.8,
                              ...) {
  
  densities <- kde3d(ordination$PCoA1, 
                     ordination$PCoA2,
                     ordination$PCoA3,
                     n = nsegments)

  # Divide data by quantiles
  x_bins <- cut(as.matrix(ordination)[,1], 
                breaks = densities$x, 
                include.lowest = TRUE, 
                labels = FALSE)
  
  y_bins <- cut(as.matrix(ordination)[,2], 
                breaks = densities$y, 
                include.lowest = TRUE, 
                labels = FALSE)
  
  z_bins <- cut(as.matrix(ordination)[,3], 
                breaks = densities$z, 
                include.lowest = TRUE, 
                labels = FALSE)
  
  # Togeter data and quantiles
  frame <- data.table(x = x_bins, 
                      y = y_bins, 
                      z = z_bins,
                      row = 1:nrow(ordination))
  
  frame$total <- 0
  
  for(i in 1:nrow(frame)) {
    frame$total[i] <- densities$d[frame$x[i],frame$y[i], frame$z[i]]
  }
  
  frame$total <- 1 - frame$total/max(frame$total)
  
  # Normalize
  frame$total <- (frame$total - min(frame$total)) / (max(frame$total) - min(frame$total))
  frame$total <- (1 - lower.p) * frame$total + lower.p
  
  # Get sample
  sub_sample <- sample(frame$row, 
                       ceiling(nrow(frame)*proportion), 
                       prob = frame$total)
  
  #plot(ordination[,c(2, 3)])
  #plot(ordination[sub_sample,c(2, 3)])
  #points(ordination[!sub_sample,c(2, 3)], col = "red")
  
  return(sub_sample)
  
}


