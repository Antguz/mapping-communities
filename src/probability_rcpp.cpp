#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec probability(arma::mat axes, 
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
    
    ans[r] = (exp(value)/(1+exp(value)));
    
  }
  
  return ans;
  
}

// [[Rcpp::export]]
List multi_probability(arma::mat axes, 
                       arma::mat mcoefficients,
                       arma::vec plants) {
  
  // N observations and axes
  int nr = axes.n_rows;
  int ncoef = mcoefficients.n_rows;
  
  // Vector to fill
  arma::mat prob(nr, ncoef);
  arma::mat pcol;
  
  // Sum the intercept and get probability
  for (int j = 0; j < ncoef; j++) {
    
    arma::rowvec coef = mcoefficients.row(j);
    arma::mat pcol = probability(axes, coef);
    prob.col(j) = pcol.reshape(nr, 1);
    
  }
  
  // Find the unique elements in the vector
  arma::vec unique_plants = arma::unique(plants);
  int nplants = unique_plants.n_elem;
  
  // Matrix to return
  arma::mat mean_prob(nr, nplants);
  arma::mat error_prob(nr, nplants);
  
  // Get mean and error
  for(int k = 0; k < nr; k++) {
    
    // Select pixel
    arma::rowvec pixel = prob.row(k);
    
    for(int i = 0; i < nplants; i++) {
      
      // Select plants within pixel for metrics
      arma::vec values = pixel(find(plants == unique_plants[i]));
      
      double n_root = sqrt(values.n_elem);
      double pixel_mean = mean(values);
      double pixel_sd = stddev(values);
      double pixel_z = 1.96*(pixel_sd/n_root);
      double pixel_error = (pixel_mean + pixel_z) - (pixel_mean - pixel_z);
      
      //Fill results
      mean_prob(k, i) = round(pixel_mean*10000);
      error_prob(k, i) = round(pixel_error*1000000);
    
    }
   }
  
  // Return a list containing two matrices
  return List::create(mean_prob, error_prob);
  
}