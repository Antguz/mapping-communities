#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::vec cwm(arma::mat axes, 
              arma::rowvec coefficients) {
  
  // N observations and axes
  int nr = axes.n_rows;
  int nc = axes.n_cols;
  
  // Vector to fill
  arma:vec ans(nr);
  
  // Sum the intercept and get probability
  for (int r = 0; r < nr; r++) {
    
    double value = coefficients[0];
    
    for (int c = 0; c < nc; c++) {
      value += coefficients[(c+1)] * axes(r, c);
    }
    
    ans[r] = value;
    
  } 
  
  return ans;
  
} 

// [[Rcpp::export]]
arma::mat multi_cwm(arma::mat axes, 
                    arma::mat mcoefficients) {
  
  // N observations and axes
  int nr = axes.n_rows;
  int ncoef = mcoefficients.n_rows;
  
  // Vector to fill
  arma::mat prob(nr, ncoef);
  arma::mat pcol;
  
  // Sum the intercept and get probability
  for (int j = 0; j < ncoef; j++) {
    
    arma::rowvec coef = mcoefficients.row(j);
    arma::mat pcol = cwm(axes, coef);
    prob.col(j) = pcol.reshape(nr, 1);
    
  }  
  
  return prob;
  
}  