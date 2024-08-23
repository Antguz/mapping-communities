#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::mat summary_rcpp(arma::mat slide_results) {
  
  // Dimensions
  int npixels = slide_results.n_rows; //n pixels
  int niterations = slide_results.n_cols; //n iterations
  int nmodels = slide_results.n_elem; //n models
  
  //Root of nmodels for IC
  double n_root = sqrt(nmodels);
  
  //Mat to return
  arma::mat result(npixels, 2);
  
  for(int i = 0; i < npixels; i++) {
    
    //Pixel of interest
    arma::rowvec ipixel = slide_results.row(i);
    
    //Pixel of interest
    double mean_pixel = mean(ipixel);
    double sd_pixel = stddev(ipixel);
    double z_pixel = 1.96*(sd_pixel/n_root);
    double amplitude_pixel = (mean_pixel + z_pixel) - (mean_pixel - z_pixel);
    
    //Fill results
    result(i, 0) = mean_pixel;
    result(i, 1) = amplitude_pixel;
    
  }
  
  return result;
  
} 
