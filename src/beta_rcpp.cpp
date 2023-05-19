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
arma::vec betaObs_rcpp(arma::mat commBoth, arma::vec xtree_1, arma::mat xtree_2) {
  
  arma::mat BA = commBoth * xtree_2; // matrix samples X branches
  
  //Apply function
  int ncol = BA.n_cols;
  arma::vec minBA(ncol);
  
  for(int i = 0; i < ncol; i++) {
    minBA[i] = min(BA.col(i));
  }
  
  double a = sum(xtree_1 % minBA);
  
  arma::vec b_diff(ncol);
  arma::vec c_diff(ncol);
  
  for(int i = 0; i < ncol; i++) {
    
    double diff = (BA(0, i) - BA(1, i)) * xtree_1[i];
    
    if(diff < 0) {
      b_diff[i] = 0;
    } else {
      b_diff[i] = diff;
    }
    
    if(diff > 0) {
      c_diff[i] = 0;
    } else {
      c_diff[i] = diff * -1;
    }
  }
  
  double b = sum(b_diff);
  double c = sum(c_diff);
  
  double denominator = (a + b + c) + a;
  
  //Beta diversity
  arma::vec beta(3); 
  
  //Total
  beta(0) = (b+c)/denominator; 
  
  //Replacement
  arma::vec bc(2);
  bc(0) = b;
  bc(1) = c;
  beta(1) = 2*(min(bc)/denominator);
  
  //Richness
  double substraction = b - c;
  beta(2) = fabs(substraction)/denominator;
  
  return beta;
  
}

// [[Rcpp::export]]
arma::cube beta_arma(arma::mat comm, arma::vec xtree1, arma::mat xtree2, int threads = 1, bool progress = true) {
  
  //Matrix of beta
  int ncomm = comm.n_rows;
  arma:cube distance(ncomm, ncomm, 3);
  
  //create progress
  Progress p(ncomm*ncomm, progress);
  
  // Parallel for speed
#ifdef _OPENMP
  if ( threads > 0 ) {
    omp_set_num_threads( threads );
  }
#endif
  
  //Loop for distance
#pragma omp parallel for
  for(int k = 0; k < ncomm; k++) {
    for(int j = 0; j < ncomm; j++) {
      
      if (! Progress::check_abort() ) {
        p.increment(); // update progress
      }
      
      arma::uvec indices(2);
      indices[0] = k;
      indices[1] = j;
      
      arma::mat commBoth = comm.rows(indices);
      arma::vec beta_comm = betaObs_rcpp(commBoth, xtree1, xtree2);
      
      distance(k, j, 0) = beta_comm[0];
      distance(k, j, 1) = beta_comm[1];
      distance(k, j, 2) = beta_comm[2];
      
    }
  }
  
  return distance;
  
}