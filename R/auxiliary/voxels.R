sourceCpp("src/voxelization_rcpp.cpp")
voxels <- function(cloud, edge_length, threads = 1L, obj.voxels = TRUE) {
  
  vox <- as.data.table(voxelization_rcpp(as.matrix(cloud), edge_length, threads))
  colnames(vox) <- c("X", "Y", "Z")
  
  vox <- vox[ , .N, by = .(X, Y, Z)] #Count the number of points per voxel
  
  if(obj.voxels == TRUE) {
    parameter <- edge_length
    names(parameter) <- c("X.size", "Y.size", "Z.size")
    
    final <- list(cloud = cloud, parameter = parameter, voxels = vox)
    class(final) <- "voxels"
    
  } else {
    final <- vox
  }
  return(final)
}