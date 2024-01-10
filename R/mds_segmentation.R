mds_segmentation <- function(MDS, proportions = 0.80, nsegments = 30) {
  
  # Get quantiles
  quantile_values <- quantile(MDS, 
                              probs = seq(0, 1, length.out = nsegments + 1), 
                              na.rm = TRUE)
  # Divide data by quantiles
  quantile_bins <- cut(MDS, 
                       breaks = quantile_values, 
                       include.lowest = TRUE, 
                       labels = FALSE)

  # Togeter data and quantiles
  frame <- data.table(MDS = MDS, 
                      segment = quantile_bins,
                      row = 1:length(MDS))
  
  # Count samples per segment
  n_frame <- frame[, .N, by = "segment"]
  min_sample <- min(n_frame$N)
  unique_segments <- unique(frame$segment)
  
  rows <- numeric()
  
  # Select rows
  for(i in 1:length(unique_segments)) {
    
    sub <- subset(frame, segment == unique_segments[i])
    sub_sample <- sample(sub$row, ceiling(proportions*nrow(sub)), replace = FALSE)
    rows <- c(rows, sub_sample)
    
  }
  
  return(rows)
  
}
