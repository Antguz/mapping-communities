################################################################################
### A modification from picante::comdist for speed using RcppArmadillo
################################################################################

library(Rcpp)
library(RcppArmadillo)
sourceCpp("src/distmat_rcpp.cpp") 

comdist_rcpp <- function(comm, dis, abundance.weighted= FALSE, threads = 1) {
  
  # Test match
  dat <- match.comm.dist(comm, dis)
  x <- dat$comm
  dis <- as.matrix(dat$dist)
  
  # Test weighted
  if (!abundance.weighted) {
    x <- decostand(x, method="pa")
  }
  
  x <- decostand(x, method="total", MARGIN=1)
  
  # Get distance matrix
  comdist <- distmat_rcpp(as.matrix(x), dis, threads = threads)
  
  # Provide com names
  row.names(comdist) <- row.names(x)
  colnames(comdist) <- row.names(x)
  
  return(as.dist(comdist))
  
}