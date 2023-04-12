#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::mat product(arma::mat multiplication, arma::vec intercept) {
  
  int nr = multiplication.n_rows;
  int nc = multiplication.n_cols;
  
  arma:vec ans(nr);
  
  for (int r = 0; r < nr; r++) {
    
    double value = intercept[r];
    
    for (int c = 0; c < nc; c++) {
      value += multiplication(r, c);
    }
    
    ans[r] = exp(value);
    
  }
  
  return ans;
  
}

// [[Rcpp::export]]
arma::mat predicted_prob(arma::mat spectra, 
                         arma::mat MDS1, 
                         arma::mat MDS2,
                         arma::mat MDS3) {
  
  //Dim
  int an = spectra.n_rows; //n pixels
  int nmodels = MDS1.n_elem; //n models
  
  //Root of nmodels for IC
  double n_root = sqrt(nmodels);
  
  //Mat to return
  arma::mat result(an, 6);
  
  //Get coefficients
  arma::mat coef_MDS1 = MDS1.cols(1, 171);
  arma::mat coef_MDS2 = MDS2.cols(1, 171);
  arma::mat coef_MDS3 = MDS3.cols(1, 171);
  
  for(int i = 0; i < an; i++) {
    
    //Pixel of interest
    arma::rowvec ipixel = spectra.row(i);
    
    //N models multiplication (m)
    arma::mat mMDS1 = MDS1.each_row() % ipixel;
    arma::mat mMDS2 = MDS2.each_row() % ipixel;
    arma::mat mMDS3 = MDS3.each_row() % ipixel;
    
    //Get product (pro_)
    arma::vec pro_MDS1 = product(mMDS1, MDS1.col(0));
    arma::vec pro_MDS2 = product(mMDS1, MDS2.col(0));
    arma::vec pro_MDS3 = product(mMDS1, MDS3.col(0));
    
    //Get mean
    double mean_healthy = mean(prob.col(0));
    double mean_symptomatic = mean(prob.col(1));
    double mean_dead = mean(prob.col(2));
    
    //Get sd
    double sd_healthy = stddev(prob.col(0));
    double sd_symptomatic = stddev(prob.col(1));
    double sd_dead = stddev(prob.col(2));
    
    //Get limits
    double z_healthy = 1.96*(sd_healthy/n_root);
    double z_symptomatic = 1.96*(sd_symptomatic/n_root);
    double z_dead = 1.96*(sd_dead/n_root);
    
    //Difference
    double amplitude_healthy = (mean_healthy + z_healthy) - (mean_healthy - z_healthy);
    double amplitude_symptomatic = (mean_symptomatic + z_symptomatic) - (mean_symptomatic - z_symptomatic);
    double amplitude_dead = (mean_dead + z_dead) - (mean_dead - z_dead);
    
    //Fill results
    result(i, 0) = round(mean_healthy*10000);
    result(i, 1) = round(mean_symptomatic*10000);
    result(i, 2) = round(mean_dead*10000);
    result(i, 3) = round(amplitude_healthy*10000);
    result(i, 4) = round(amplitude_symptomatic*10000);
    result(i, 5) = round(amplitude_dead*10000);
    
  }
  
  return result;
  
}
