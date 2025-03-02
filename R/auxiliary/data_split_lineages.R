data_split_lineages <- function(aggregation,
                                split = 0.6) {
  
  # Create data for training and testing
  PLT_CN <- aggregation$PLT_CN
  unique_lineages <- colnames(aggregation)[-1]
  data_split <- matrix(NA, nrow = nrow(aggregation), ncol = length(unique_lineages))
  colnames(data_split) <- unique_lineages
  
  for(i in 1:length(unique_lineages)) {
    
    spl <- aggregation[, .SD, .SDcols = c(1, (i+1))]
    colnames(spl)[2] <- "name"
    spl$split <- 0
    
    # Presence
    presence <- spl[name == 1,]
    set.seed(i)
    train_presence <- sample(1:nrow(presence), floor(nrow(presence)*split))
    presence[train_presence, split := 1] #Training
    
    # Absence
    absence <- spl[name == 0,]
    set.seed(i)
    train_absence <- sample(1:nrow(absence), floor(nrow(absence)*split))
    absence[train_absence, split := 1] #Training
    
    # Merge
    spl <- rbind(presence, absence)
    spl_match <- match.integer64(PLT_CN, spl$PLT_CN)
    spl <- spl[spl_match,]
    data_split[, i] <- spl$split
    
  }
  
  return(data_split)
  
}
