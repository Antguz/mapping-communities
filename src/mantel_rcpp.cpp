#ifdef _OPENMP
#include <omp.h>
#endif
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <algorithm> // for std::random_shuffle
#include <numeric>   // for std::iota
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]

// [[Rcpp::depends(RcppArmadillo)]]

// Function to compute the Mantel statistic (Pearson correlation)
// between two distance matrices
double mantel_stat(const arma::mat& mat1, const arma::mat& mat2) {
  arma::uvec upper_tri_idx = arma::trimatu_ind(arma::size(mat1), 1);
  arma::vec vec1 = mat1(upper_tri_idx);
  arma::vec vec2 = mat2(upper_tri_idx);
  return arma::as_scalar(arma::cor(vec1, vec2));
}

// Function to perform the Mantel test with permutations
// [[Rcpp::export]]
Rcpp::List mantel_rcpp(const arma::mat& mat1, 
                       const arma::mat& mat2, 
                       int nperm = 999,
                       int threads = 1, 
                       bool progress = true) {
  
  // Compute the observed Mantel statistic
  double observed_stat = mantel_stat(mat1, mat2);
  
  // Create a vector to store permutation results
  arma::vec perm_stats(nperm);
  
  // Vector for indices to shuffle
  arma::uvec idx(mat1.n_rows);
  std::iota(idx.begin(), idx.end(), 0);  // Fill with 0, 1, 2, ..., n-1
  
  //create progress
  Progress p(nperm, progress);
  
  // Parallel for speed
#ifdef _OPENMP
  if ( threads > 1 ) {
    omp_set_num_threads( threads );
  }
#endif 

#pragma omp parallel for
  for (int i = 0; i < nperm; i++) {
    
    if (! Progress::check_abort() ) {
      p.increment(); // update progress
    } 
    
    std::random_shuffle(idx.begin(), idx.end());
    arma::mat permuted_mat2 = mat2(idx, idx);  // Permute the rows and columns
    perm_stats[i] = mantel_stat(mat1, permuted_mat2);
    
  }
  
  return Rcpp::List::create(
    Rcpp::Named("observed_stat") = observed_stat,
    Rcpp::Named("perm_stats") = perm_stats
  );
}