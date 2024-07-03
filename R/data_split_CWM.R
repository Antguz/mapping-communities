data_split_CWM <- function(aggregation,
                           split = 0.6) {
  
  # Create data for training and testing
  PLT_CN <- aggregation$PLT_CN
  data_split <- matrix(NA, nrow = nrow(aggregation), ncol = (length(aggregation)-1))
  colnames(data_split) <- colnames(aggregation)[-1]
  
  for(i in 1:ncol(data_split)) {
    
    # Select trait
    trait <- aggregation[, .SD, .SDcols = c(1, (i+1))]
    colnames(trait)[2] <- "trait"
    
    # Select for training
    trait$cut <- cut(trait$trait, breaks = 20)
    trait$train <- 0
    cuts_unite <- unique(trait$cut)
    all <- data.frame()
    
    for(ii in 1:length(cuts_unite)) {
      
      sub <- subset(trait, cut == cuts_unite[ii])
      sample_trait <- sample(1:nrow(sub), ceiling(nrow(sub)*split))
      sub[sample_trait, train := 1]
      all <- rbind(all, sub)
      
    }
    
    all_match <- match.integer64(PLT_CN, all$PLT_CN)
    all <- all[all_match,]
    data_split[,i] <- all$train
    
  }
  
  return(data_split)
  
}
