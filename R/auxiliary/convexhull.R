estimate_hull_volume <- function(group, scores) {
  
  # Clean group
  observed <- group
  group_name <- colnames(observed)[2]
  colnames(observed)[2] <- "occurrance"
  observed <- subset(observed, occurrance == 1)
  
  # Merge group with pca
  observed <- merge(observed, scores, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
  observed_scores <- observed[, 3:5]
  observed_scores <- as.matrix(na.exclude(observed_scores))
  
  # Compute convex hull
  hull <- tryCatch(
    convhulln(observed_scores, options = "FA"),
    error = function(e) NULL
  )
  
  return(data.table(group = group_name, hull_vol = hull$vol))

}