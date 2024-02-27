data_split <- function(meta, ptraining = 0.6) {
 
  frame <- meta
  frame$site <- paste0(frame$STATECD, "_", frame$COUNTYCD)
  frame$ID <- 1:nrow(frame)
  
  sites <- unique(frame$site)
  training <- numeric()
  
  for(i in 1:length(sites)) {
    
    sub_site <- subset(frame, site == sites[i])
    selected <- sample(1:nrow(sub_site), floor(nrow(sub_site)*ptraining))
    training <- c(training, sub_site$ID[selected])
    
  }
  
  return(training)
  
}
