mds_segmentation <- function(MDS, nsegments = 30) {
  
  min_MDS <- min(MDS)
  segment_lengt <- (max(MDS) - min_MDS) / (nsegments - 1)
  segment <- NA
  
  for(i in 1:length(MDS)) {
    segment[i] <- floor(((MDS[i] - min_MDS)/segment_lengt));
  }

  frame <- data.table(MDS = MDS, 
                      segment = segment,
                      row = 1:length(MDS))  
  quant <- quantile(frame[, .N, by = segment]$N, 0.25)
  unique_segments <- unique(frame$segment)
  
  rows <- 1
  rows <- rows[-1]
  
  for(i in 1:length(unique_segments)) {
    
    sub <- subset(frame, segment == unique_segments[i])
    
    if(quant >= nrow(sub)) {
      rows <- c(rows, sub$row)
      
    } else {
      sub_sample <- sample(sub$row, quant, replace = FALSE)
      rows <- c(rows, sub_sample)
    }
  }
  
  return(rows)
  
}
