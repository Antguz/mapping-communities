#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// Prediction of probabilities using coefficients for a single model
// [[Rcpp::export]]
arma::mat prob_coeff_prediction(arma::mat axes,
                                arma::mat coefficients,
                                arma::vec intercept) {
  
  // axes: cols are axes, rows observations
  // coefficients: cols are coefficients, rows iterations
  // intercept: a vector of intercepts
  
  // N observations and axes
  int nsamples = axes.n_rows;
  int niterations = coefficients.n_rows;
  int naxes = axes.n_cols;
  
  // Vector to fill
  arma::mat results(nsamples, niterations);
  
  // Sum the intercept and get probability
  for (int r = 0; r < nsamples; r++) {
    
    //Pixel of interest
    arma::rowvec ipixel = axes.row(r);
    
    //N models multiplication (m)
    arma::mat coef_multi = coefficients.each_row() % ipixel;
    arma::vec coef_sum = arma::sum(coef_multi, 1) + intercept;
    arma::vec coef_prob = (1/(1+exp(coef_sum)));
    
    // Storage result
    results.row(r) = trans(coef_prob);
    
  }
  
  return results;
  
}  

//Prediction of probabilities using coefficients for several models
// [[Rcpp::export]]
arma::cube cube_prob_coeff_prediction(arma::mat axes,
                                      arma::cube cube_coefficients, 
                                      arma::mat mat_intercept) { 
  
  // axes: cols are axes, rows observations
  // cube_coefficients: cols are coefficients, rows iterations, slices are models
  // mat_intercept: cols are traits, rows iterations
   
  // Get n of traits coefficients to predict
  int samples = axes.n_rows;
  int nmodels = cube_coefficients.n_slices;
  int iterations = cube_coefficients.n_rows;
   
  // Vector to fill
  arma::cube results(samples, iterations, nmodels);
  
  // Sum the intercept and get probability
  for (int s = 0; s < nmodels; s++) {
    
    //Get coefficients and intercept of the trait of interest
    arma::mat coef = cube_coefficients.slice(s);
    arma::vec interp = mat_intercept.col(s);
    
    //Get prediction by coefficients
    arma::mat prob_predicted = prob_coeff_prediction(axes, coef, interp);
    
    //Storage
    results.slice(s) = prob_predicted;
  }
  
  return results;
   
}

// Summary of prediction of CWM using coefficients for several traits
// [[Rcpp::export]]
List summary_cwm(arma::mat axes,
                 arma::cube cube_coefficients, 
                 arma::mat mat_intercept) {
  
  // Get n of traits coefficients to predict
  int samples = axes.n_rows;
  int nmodels = cube_coefficients.n_slices;
  int iterations = cube_coefficients.n_rows;
  
  // Get results
  arma::cube results = cube_prob_coeff_prediction(axes,
                                                  cube_coefficients,
                                                  mat_intercept);
  
  // Mat to fill
  arma::mat summary_mean(samples, nmodels);
  arma::mat summary_error(samples, nmodels);
  
  // Sum the intercept and get probability
  for(int t = 0; t < nmodels; t++) {
    
    // Get a slice of results
    arma::mat sub = results.slice(t);
     
    // Get summary
    for(int p = 0; p < samples; ++p) {
      
      // Get metrics
      double n_root = sqrt(iterations);
      double pixel_mean = arma::mean(sub.row(p));
      double pixel_sd = arma::stddev(sub.row(p));
      double pixel_z = 1.96*(pixel_sd/n_root);
      double pixel_error = (pixel_mean + pixel_z) - (pixel_mean - pixel_z);
      
      // Store results
      summary_mean(p, t) = round(pixel_mean*10000);
      summary_error(p, t) = round(pixel_error*100000000);
      
    } 
  }
  
  // Return a list containing two matrices
  return List::create(summary_mean, summary_error);
  
} 