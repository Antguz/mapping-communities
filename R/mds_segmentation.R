mds_segmentation <- function(MDS, proportions = 0.80, nsegments = 30) {
  
  min_MDS <- min(MDS)
  segment_lengt <- (max(MDS) - min_MDS) / (nsegments - 1)
  segment <- NA
  
  for(i in 1:length(MDS)) {
    segment[i] <- floor(((MDS[i] - min_MDS)/segment_lengt));
  }

  frame <- data.table(MDS = MDS, 
                      segment = segment,
                      row = 1:length(MDS))  
  unique_segments <- unique(frame$segment)
  
  rows <- numeric()
  
  for(i in 1:length(unique_segments)) {
    
    sub <- subset(frame, segment == unique_segments[i])
    sub_sample <- sample(sub$row, ceiling(proportions*nrow(sub)), replace = FALSE)
    rows <- c(rows, sub_sample)
    
  }
  
  return(rows)
  
}
