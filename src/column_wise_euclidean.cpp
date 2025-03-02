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
arma::mat column_wise_euclidean(const arma::mat X, 
                                const arma::mat Y,
                                int n_threads,
                                bool progress = true) {
  
  arma::uword n = X.n_rows;
  arma::uword ny = Y.n_rows;
  arma::mat dist_matrix(n, ny, arma::fill::zeros);
  
  //create progress
  Progress p(n*ny, progress);
  
  // Parallel for speed
#ifdef _OPENMP
  if ( n_threads > 0 ) {
    omp_set_num_threads( n_threads );
  } 
#endif
  
  //Loop for distance
#pragma omp parallel for collapse(2)
  for (arma::uword i = 0; i < n; i++) {
    for (arma::uword j = 0; j < ny; j++) {
      
      if (! Progress::check_abort() ) {
        p.increment(); // update progress
      }
      
      double dist = arma::norm(X.row(i) - Y.row(j), 2);
      dist_matrix(i, j) = dist;
      
    } 
  }
  
  return dist_matrix;
} 