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
arma::vec beta_tree_rcpp(arma::mat commBoth, 
                         arma::vec xtree_1, 
                         arma::mat xtree_2) {
  
  // Get a
  arma::mat BA = commBoth * xtree_2; // matrix samples X branches
  
  int ncol = BA.n_cols;
  arma::vec minBA(ncol);
  
  for(int i = 0; i < ncol; ++i) {
    minBA[i] = min(BA.col(i));
  }
  
  double a = sum(xtree_1 % minBA);
  
  arma::vec comm_diff(ncol);
  
  for(int i = 0; i < ncol; ++i) {
    comm_diff[i] = (BA(0, i) - BA(1, i));
  }
  
  arma::vec diff = xtree_1 % comm_diff;
  arma::vec for_b = diff;
  arma::vec for_c = diff;
  
  for_b.elem( find(diff < 0) ).zeros();
  for_c.elem( find(diff > 0) ).zeros();
  
  double b = sum(for_b);
  double c = sum(for_c * -1);
  
  double denominator = (a + b + c);
  
  //Beta diversity
  arma::vec beta(3); 
  
  //Total
  beta[0] = (b+c)/denominator; 
  
  //Replacement
  arma::vec bc(2);
  bc[0] = b;
  bc[1] = c;
  beta[1] = (2*min(bc))/denominator;
  
  //Richness
  double substraction = b - c;
  beta[2] =  std::abs(substraction)/denominator;
  
  //return beta;
  return beta;
  
}

// [[Rcpp::export]]
arma::vec beta_abundance_rcpp(arma::mat commBoth) {
  
  double a = 0;
  double b = 0;
  double c = 0;
  
  int ncol = commBoth.n_cols;
  
  for(int i = 0; i < ncol; i++) {
    
    double minComm = min(commBoth.col(i));
    a = a + minComm;
    b = b + commBoth(0, i) - minComm;
    c = c + commBoth(1, i) - minComm;
    
  }
  
  double denominator = (a + b + c);
  
  //Beta diversity
  arma::vec beta(3); 
  
  //Total
  beta[0] = (b+c)/denominator;
  
  //Replacement
  arma::vec bc(2);
  bc[0] = b;
  bc[1] = c;
  beta[1] = (2*min(bc))/denominator;
  
  //Richness
  double substraction = b - c;
  beta[2] = std::abs(substraction)/denominator;
  
  //return beta;
  return beta;
  
}

// [[Rcpp::export]]
Rcpp::List beta_arma(arma::mat comm, 
                     arma::vec xtree1, 
                     arma::mat xtree2, 
                     bool taxonomic = true, 
                     int threads = 1, 
                     bool progress = true) {
  
  //Matrix of beta
  int ncomm = comm.n_rows;
  arma::mat Btotal(ncomm, ncomm);
  arma::mat Brepl(ncomm, ncomm);
  arma::mat Brich(ncomm, ncomm);
  
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
  for(int k = 0; k < ncomm; ++k) {
    for(int j = 0; j < ncomm; ++j) {
      
      if (! Progress::check_abort() ) {
        p.increment(); // update progress
      }
      
      arma::uvec indices(2);
      indices[0] = k;
      indices[1] = j;
      
      arma::mat commBoth = comm.rows(indices);
      
      if(!taxonomic) {
        arma::vec beta_comm = beta_tree_rcpp(commBoth = commBoth, 
                                             xtree1 = xtree1, 
                                             xtree2 = xtree2);
        //Btotal(k, j) = beta_comm[0];
        //Brepl(k, j) = beta_comm[1];
        //Brich(k, j) = beta_comm[2];
        Btotal(k, j) = round(beta_comm[0]*10000);
        Brepl(k, j) = round(beta_comm[1]*10000);
        Brich(k, j) = round(beta_comm[2]*10000);
        
      } else {
        arma::vec beta_comm = beta_abundance_rcpp(commBoth = commBoth);
        //Btotal(k, j) = beta_comm[0];
        //Brepl(k, j) = beta_comm[1];
        //Brich(k, j) = beta_comm[2];
        Btotal(k, j) = round(beta_comm[0]*10000);
        Brepl(k, j) = round(beta_comm[1]*10000);
        Brich(k, j) = round(beta_comm[2]*10000);
      }
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("Btotal") = Btotal,
                            Rcpp::Named("Brepl") = Brepl,
                            Rcpp::Named("Brich") = Brich);
  
}