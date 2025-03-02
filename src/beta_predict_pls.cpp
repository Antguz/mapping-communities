#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::mat beta_predict_pls(arma::mat spectra, 
                           arma::mat Axis1, 
                           arma::mat Axis2,
                           arma::mat Axis3) {
  
  // Dimensions
  int npixels = spectra.n_rows; //n pixels
  int nbands = spectra.n_cols; //n bands
  int nmodels = Axis1.n_elem; //n models
  
  //Root of nmodels for IC
  double n_root = sqrt(nmodels);
  
  //Mat to return
  arma::mat result(npixels, 7);
  
  //Get coefficients
  arma::mat coef_Axis1 = Axis1.cols(1, nbands);
  arma::mat coef_Axis2 = Axis2.cols(1, nbands);
  arma::mat coef_Axis3 = Axis3.cols(1, nbands);
  
  for(int i = 0; i < npixels; i++) {
    
    //Pixel of interest
    arma::rowvec ipixel = spectra.row(i);
    
    //N models multiplication (m)
    arma::mat mAxis1 = coef_Axis1.each_row() % ipixel;
    arma::mat mAxis2 = coef_Axis2.each_row() % ipixel;
    arma::mat mAxis3 = coef_Axis3.each_row() % ipixel;
    
    //Dimensions
    int nr = mAxis1.n_rows;
    int nc = mAxis1.n_cols;
    
    arma::vec pAxis1(nr);
    arma::vec pAxis2(nr);
    arma::vec pAxis3(nr);
    
    for(int j = 0; j < nr; j++) {
      
      // Get intercepts
      double value_Axis1 = Axis1(j, 0);
      double value_Axis2 = Axis2(j, 0);
      double value_Axis3 = Axis3(j, 0);
      
      for(int k = 0; k < nc; k++) {
        
        // Sum intercepts and coefficients
        value_Axis1 = value_Axis1 + mAxis1(j, k);
        value_Axis2 = value_Axis2 + mAxis2(j, k);
        value_Axis3 = value_Axis3 + mAxis3(j, k);
        
      }
      
      // Predicted value per model
      pAxis1(j) = value_Axis1;
      pAxis2(j) = value_Axis2;
      pAxis3(j) = value_Axis3;
      
    }
    
    //Get mean
    double mean_Axis1 = mean(pAxis1);
    double mean_Axis2 = mean(pAxis2);
    double mean_Axis3 = mean(pAxis3);
    
    //Get sd
    double sd_Axis1 = stddev(pAxis1);
    double sd_Axis2 = stddev(pAxis2);
    double sd_Axis3 = stddev(pAxis3);
    
    //Get limits
    double z_Axis1 = 1.96*(sd_Axis1/n_root);
    double z_Axis2 = 1.96*(sd_Axis2/n_root);
    double z_Axis3 = 1.96*(sd_Axis3/n_root);
    
    //Difference
    double amplitude_Axis1 = (mean_Axis1 + z_Axis1) - (mean_Axis1 - z_Axis1);
    double amplitude_Axis2 = (mean_Axis2 + z_Axis2) - (mean_Axis2 - z_Axis2);
    double amplitude_Axis3 = (mean_Axis3 + z_Axis3) - (mean_Axis3 - z_Axis3);
    
    //Fill results
    result(i, 0) = round(mean_Axis1*10000);
    result(i, 1) = round(mean_Axis2*10000);
    result(i, 2) = round(mean_Axis3*10000);
    result(i, 3) = round(amplitude_Axis1*1000000);
    result(i, 4) = round(amplitude_Axis2*1000000);
    result(i, 5) = round(amplitude_Axis3*1000000);
    result(i, 6) = result(i, 3) + result(i, 4) + result(i, 5);

  }
  
  return result;
  
}
