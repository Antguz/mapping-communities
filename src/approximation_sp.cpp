#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::mat approximation_sp(arma::mat predicted,
                           arma::mat wascores) {
  
  // Dimensions
  int npixels = predicted.n_rows; //n pixels
  int nspecies = wascores.n_rows; //n models
  
  //Mat to return
  arma::mat result(npixels, nspecies);
  
  for(int i = 0; i < npixels; i++) {
    for(int j = 0; j < nspecies; j++) {
    
    //Estimate distance
    result(i, j) = sqrt(pow((predicted(i, 0) - wascores(j, 0)), 2.0) + pow((predicted(i, 1) - wascores(j, 1)), 2.0) + pow((predicted(i, 2) - wascores(j, 2)), 2.0));
    
    }
  }
  
  return result;
  
}
