#ifdef _OPENMP
#include <omp.h>
#endif
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace arma;

// [[Rcpp::export]]
arma::vec stresscomp(arma::vec y, arma::vec yf, int form) {
  arma::vec num = arma::pow(y - yf, 2);
  double den;
  
  if (form == 1) {
    den = arma::accu(arma::pow(y, 2));
  } else {
    den = arma::accu(arma::pow(y - mean(y), 2));
  }
  
  return num / den;
}

// [[Rcpp::export]]
arma::vec goodness_rcpp(arma::vec dist, 
                        arma::vec dhat, 
                        int nobj, 
                        int ndis, 
                        arma::vec iidx, 
                        arma::vec jidx,
                        int threads = 1,
                        bool progress = true) {
  
  arma::vec x = stresscomp(dist, dhat, 1);
  arma::mat mat = arma::zeros<arma::mat>(nobj, nobj);
  
  //create progress
  Progress p(ndis, progress);
  
  // Parallel for speed
#ifdef _OPENMP
  if (threads > 0) {
    omp_set_num_threads(threads);
  }
#endif
  
  //Loop for distance
#pragma omp parallel for
  for (int i = 0; i < (ndis - 1); ++i) {
    
    if (! Progress::check_abort() ) {
      p.increment(); // update progress
    }
    
    mat(iidx(i) - 1, jidx(i) - 1) = x(i);
    
  }
  
  arma::vec res = arma::sqrt(arma::sum(mat + mat.t(), 1) / 2);
  
  return res;
  
}