axis_segmentation <- function(ordination_split, 
                              weigths_split,
                              fraction = 0.85,
                              threads = 20) {
  # Frames
  axis1_frame <- data.table(axes = ordination_split$Axis1,
                            weigths = weigths_split,
                            groups = cut(ordination_split$Axis1,
                                         breaks = hist(ordination_split$Axis1, breaks = 10)$breaks,
                                         include.lowest = TRUE,
                                         labels = FALSE),
                            ID = 1:nrow(ordination_split))
  
  axis2_frame <- data.table(axes = ordination_split$Axis2,
                            weigths = weigths_split,
                            groups = cut(ordination_split$Axis2,
                                         breaks = hist(ordination_split$Axis2, breaks = 10)$breaks,
                                         include.lowest = TRUE,
                                         labels = FALSE),
                            ID = 1:nrow(ordination_split))
  
  axis3_frame <- data.table(axes = ordination_split$Axis3,
                            weigths = weigths_split,
                            groups = cut(ordination_split$Axis3,
                                         breaks = hist(ordination_split$Axis3, breaks = 10)$breaks,
                                         include.lowest = TRUE,
                                         labels = FALSE),
                            ID = 1:nrow(ordination_split))
  
  # Split the data by group
  axis1_grouped <- split(axis1_frame, axis1_frame$groups)
  axis2_grouped <- split(axis2_frame, axis2_frame$groups)
  axis3_grouped <- split(axis3_frame, axis3_frame$groups)
  
  # Application function
  segmentation_application <- function(X,
                                       axis1_grouped,
                                       axis2_grouped,
                                       axis3_grouped,
                                       fraction) {
    
    # Create segmentation
    proportional_axis1 <- lapply(axis1_grouped, function(axis, frac = fraction) {
      sample_size <- ceiling(frac * nrow(axis))  # 10% of group size
      if(nrow(axis) > sample_size) {
        axis$ID[sample(1:nrow(axis), sample_size, prob = axis$weigths)]
      } else {
        axis$ID
      }
    })
    proportional_axis1 <- do.call(c, proportional_axis1)
    
    proportional_axis2 <- lapply(axis2_grouped, function(axis, frac = fraction) {
      sample_size <- ceiling(frac * nrow(axis))  # 10% of group size
      if(nrow(axis) > sample_size) {
        axis$ID[sample(1:nrow(axis), sample_size, prob = axis$weigths)]
      } else {
        axis$ID
      }
    })
    proportional_axis2 <- do.call(c, proportional_axis2)
    
    proportional_axis3 <- lapply(axis3_grouped, function(axis, frac = fraction) {
      sample_size <- ceiling(frac * nrow(axis))  # 10% of group size
      if(nrow(axis) > sample_size) {
        axis$ID[sample(1:nrow(axis), sample_size, prob = axis$weigths)]
      } else {
        axis$ID
      }
    })
    proportional_axis3 <- do.call(c, proportional_axis3)
    
    return(list(axis1_seg = proportional_axis1,
                axis2_seg = proportional_axis2,
                axis3_seg = proportional_axis3))
    
  }
  
  seg <- pbmclapply(X = 1:iterations,
                    FUN = segmentation_application,
                    axis1_grouped = axis1_grouped,
                    axis2_grouped = axis2_grouped,
                    axis3_grouped = axis3_grouped,
                    fraction = 0.85,
                    mc.set.seed = TRUE,
                    mc.cores = threads)
  
  return(seg)
  
}
