### This is a modified version of the MDSrotate function of vegan for ordinations ###

MDS_rotate <- function(ordination, vec, na.rm = TRUE) {
  
  x <- as.matrix(ordination)
  N <- NCOL(x)
  
  if (N < 2) {
    stop(gettextf("needs at least two dimensions"))
  }
  
  vec <- as.matrix(vec)
  NV <- NCOL(vec)
  
  if (NV >= N) {
    stop(gettextf("you can have max %d vectors, but you had %d", N-1, NV))
  }
  
  if (!is.numeric(vec)) {
    stop(gettextf("'vec' must be numeric"))
  }
  
  # vectorfit finds the direction cosine. We rotate first axis to
  # 'vec' which means that we make other axes orthogonal to 'vec'
  # one by one
  
  if (na.rm) {
    keep <- complete.cases(vec)
  } else {
    keep <- !logical(NROW(vec))
  }
  
  ## Rotation loop
  for(v in seq_len(NV)) {
    for (k in (v+1):N) {
      
      arrs <- vectorfit(x[keep,], vec[keep,v], permutations = 0)$arrows
      rot <- arrs[c(v,k)]/sqrt(sum(arrs[c(v,k)]^2))
      rot <- drop(rot)
      
      ## counterclockwise rotation matrix:
      ## [cos theta   -sin theta]
      ## [sin theta    cos theta]
      rot <- rbind(rot, rev(rot))
      rot[1,2] <- -rot[1,2]
      
      ## Rotation of points
      x[, c(v,k)] <- x[, c(v,k)] %*% rot
      
    }
  }
  
  x
  
}
