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
double spectral_angle_mapper_cpp(const arma::vec& spectrum1, const arma::vec& spectrum2) {
  // Normalize spectral vectors
  arma::vec norm_spectrum1 = spectrum1 / arma::norm(spectrum1, 2);
  arma::vec norm_spectrum2 = spectrum2 / arma::norm(spectrum2, 2);
  
  // Compute dot product
  double dot_product = arma::dot(norm_spectrum1, norm_spectrum2);
  
  // Compute spectral angle
  double spectral_angle = (std::acos(dot_product) * 180/3.141593)/180;
  
  //return sam_value;
  return spectral_angle;
}

// [[Rcpp::export]]
arma::mat sam(arma::mat data, int threads = 1, bool progress = true) {
  
  //Matrix of data
  int n = data.n_rows;
  arma::mat distances(n, n);
  
  //create progress
  Progress p(n*n, progress);
  
  // Parallel for speed
#ifdef _OPENMP
  if ( threads > 0 ) {
    omp_set_num_threads( threads );
  }
#endif
  
  //Loop for distance
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      
      if (! Progress::check_abort() ) {
        p.increment(); // update progress
      } 
      
      if(i == j) {
        distances(i, j) = 0;
        
      } else {
        
        arma::vec a = arma::trans(data.row(i));
        arma::vec b = arma::trans(data.row(j));
        
        distances(i, j) = spectral_angle_mapper_cpp(a, b);
        
      }
    }
  }
  
  return distances;
  
}
