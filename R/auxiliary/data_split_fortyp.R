data_split_fortyp <- function(forest_type,
                              split = 0.6) {
  
  # Create data for training and testing
  PLT_CN <- forest_type$PLT_CN
  unique_fortyp <- unique(forest_type$FORTYP)
  data_split <- matrix(NA, nrow = nrow(forest_type), ncol = length(unique_fortyp))
  colnames(data_split) <- unique_fortyp
  
  for(i in 1:length(unique_fortyp)) {
    
    spl <- forest_type
    spl$name <- 0
    spl[FORTYP == unique_fortyp[i], name := 1]
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
