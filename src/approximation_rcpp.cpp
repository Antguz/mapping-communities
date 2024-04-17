#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::cube approximation_rcpp(arma::mat predicted,
                              arma::mat wascores) {
  
  // Dimensions
  int npixels = predicted.n_rows; //n pixels
  int nspecies = wascores.n_rows; //n models
  int naxes = wascores.n_cols; //
  
  //Mat to return
  arma::cube result(npixels, nspecies, naxes);
  
  for(int i = 0; i < npixels; i++) {
    for(int j = 0; j < nspecies; j++) {
      for(int k = 0; k < naxes; k++) {
        
        //Estimate distance
        result(i, j, k) = sqrt(pow((predicted(i, k) - wascores(j, k)), 2.0));
      
      }
    }
  }
  
  return result;
  
}
