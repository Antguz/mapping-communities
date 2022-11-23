CV_spatialBlock <- function(spatial, k = 10) {
  
  cv <- spatialBlock(speciesData = spatial,
                     theRange = 30000,
                     k = k,
                     selection = "random",
                     #iteration = 100,
                     numLimit = 10,
                     showBlocks = FALSE,
                     progress = FALSE)
  
  return(cv)
  
}
