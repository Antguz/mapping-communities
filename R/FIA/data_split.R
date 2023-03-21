data_split <- function(frame, training = 0.6) {
  
  frame <- frame
  frame$sample <- 1:nrow(frame)
  county <- unique(frame$COUNTYCD)
  to_training <- as.numeric()
  
  for(i in 1:length(county)) {
    
    x_tile <- subset(frame, COUNTYCD == county[i])
    min_samples <- min(table(x_tile$COUNTYCD))
    min_samples <- floor(min_samples * training)
    
    healthy <- subset(x_tile, Condition == "Healthy")
    symptomatic <- subset(x_tile, Condition == "Symptomatic")
    dead <- subset(x_tile, Condition == "Dead")
    
    healthy <- healthy[sample(1:nrow(healthy), min_samples), sample]
    symptomatic <- symptomatic[sample(1:nrow(symptomatic), min_samples), sample]
    dead <- dead[sample(1:nrow(dead), min_samples), sample]
    
    to_training <- c(to_training, healthy, symptomatic, dead)
    
  }
  
  return(to_training)
  
}