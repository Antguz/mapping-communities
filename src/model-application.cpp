#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::mat predict_pls(arma::mat spectra, 
                      arma::mat MDS1, 
                      arma::mat MDS2,
                      arma::mat MDS3) {
  
  // Dimensions
  int npixels = spectra.n_rows; //n pixels
  int nbands = spectra.n_cols; //n bands
  int nmodels = MDS1.n_elem; //n models
  
  //Root of nmodels for IC
  double n_root = sqrt(nmodels);
  
  //Mat to return
  arma::mat result(npixels, 7);
  
  //Get coefficients
  arma::mat coef_MDS1 = MDS1.cols(1, nbands);
  arma::mat coef_MDS2 = MDS2.cols(1, nbands);
  arma::mat coef_MDS3 = MDS3.cols(1, nbands);
  
  for(int i = 0; i < npixels; i++) {
    
    //Pixel of interest
    arma::rowvec ipixel = spectra.row(i);
    
    //N models multiplication (m)
    arma::mat mMDS1 = coef_MDS1.each_row() % ipixel;
    arma::mat mMDS2 = coef_MDS2.each_row() % ipixel;
    arma::mat mMDS3 = coef_MDS3.each_row() % ipixel;
    
    //Dimensions
    int nr = mMDS1.n_rows;
    int nc = mMDS1.n_cols;
    
    arma::vec pMDS1(nr);
    arma::vec pMDS2(nr);
    arma::vec pMDS3(nr);
    
    for(int j = 0; j < nr; j++) {
      
      // Get intercepts
      double value_MDS1 = MDS1(j, 0);
      double value_MDS2 = MDS2(j, 0);
      double value_MDS3 = MDS3(j, 0);
      
      for(int k = 0; k < nc; k++) {
        
        // Sum intercepts and coefficients
        value_MDS1 = value_MDS1 + mMDS1(j, k);
        value_MDS2 = value_MDS2 + mMDS2(j, k);
        value_MDS3 = value_MDS3 + mMDS3(j, k);
        
      }
      
      // Predicted value per model
      pMDS1(j) = value_MDS1;
      pMDS2(j) = value_MDS2;
      pMDS3(j) = value_MDS3;
      
    }
    
    //Get mean
    double mean_MDS1 = mean(pMDS1);
    double mean_MDS2 = mean(pMDS2);
    double mean_MDS3 = mean(pMDS3);
    
    //Get sd
    double sd_MDS1 = stddev(pMDS1);
    double sd_MDS2 = stddev(pMDS2);
    double sd_MDS3 = stddev(pMDS3);
    
    //Get limits
    double z_MDS1 = 1.96*(sd_MDS1/n_root);
    double z_MDS2 = 1.96*(sd_MDS2/n_root);
    double z_MDS3 = 1.96*(sd_MDS3/n_root);
    
    //Difference
    double amplitude_MDS1 = (mean_MDS1 + z_MDS1) - (mean_MDS1 - z_MDS1);
    double amplitude_MDS2 = (mean_MDS2 + z_MDS2) - (mean_MDS2 - z_MDS2);
    double amplitude_MDS3 = (mean_MDS3 + z_MDS3) - (mean_MDS3 - z_MDS3);
    
    //Fill results
    result(i, 0) = round(mean_MDS1*1000);
    result(i, 1) = round(mean_MDS2*1000);
    result(i, 2) = round(mean_MDS3*1000);
    result(i, 3) = round(amplitude_MDS1*1000);
    result(i, 4) = round(amplitude_MDS2*1000);
    result(i, 5) = round(amplitude_MDS3*1000);
    result(i, 6) = result(i, 3) + result(i, 4) + result(i, 5);

  }
  
  return result;
  
}
