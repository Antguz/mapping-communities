grid_selection <- function(ordination, steps = 30) {
  
  rMDS1 <- round((max(ordination[,1]) - min(ordination[,1])) / steps, 3)
  rMDS2 <- round((max(ordination[,2]) - min(ordination[,2])) / steps, 3)
  rMDS3 <- round((max(ordination[,3]) - min(ordination[,3])) / steps, 3)
  
  grid <- voxels(ordination, 
                 edge_length = c(rMDS1, rMDS2, rMDS3), 
                 obj.voxels = FALSE)
  
  ordination$ID <- 1:nrow(ordination)
  
  samples <- NA
  
  for(i in 1:nrow(grid)) {
    
    sub_ordination <- ordination[MDS1 >= (grid$X[i] - rMDS1/2) & MDS1 < (grid$X[i] + rMDS1/2) &
                                 MDS2 >= (grid$Y[i] - rMDS2/2) & MDS2 < (grid$Y[i] + rMDS2/2) & 
                                 MDS3 >= (grid$Z[i] - rMDS3/2) & MDS3 < (grid$Z[i] + rMDS3/2), ]
    
    if(nrow(sub_ordination) == 1) {
      samples[i] <- sub_ordination$ID[1]
    } else {
      select <- sample(nrow(sub_ordination), 1)
      samples[i] <- sub_ordination$ID[select]
    }
  }
  
  return(samples)
  
}
