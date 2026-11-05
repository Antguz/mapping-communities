#ifdef _OPENMP
#include <omp.h>
#endif
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppProgress)]]

// Sum of mat1(i, j) * mat2(perm[i], perm[j]) over the upper triangle (i < j)
double cross_sum(const arma::mat& mat1,
                 const arma::mat& mat2,
                 const arma::uvec& perm) {

  int n = mat1.n_rows;
  double s = 0;

  for (int j = 1; j < n; ++j) {
    const double* b = mat2.colptr(perm[j]);
    const double* a = mat1.colptr(j);
    for (int i = 0; i < j; ++i) {
      s += a[i] * b[perm[i]];
    }
  }

  return s;
}

// Function to perform the Mantel test with permutations of mat2.
// The Mantel statistic is the Pearson correlation between the upper triangles.
// Permutations are drawn with the R random number generator (reproducible with
// set.seed), or given in perms (n x nperm, 1-based as in vegan).
// [[Rcpp::export]]
Rcpp::List mantel_rcpp(const arma::mat& mat1,
                       const arma::mat& mat2,
                       int nperm = 999,
                       int threads = 1,
                       bool progress = true,
                       Rcpp::Nullable<Rcpp::IntegerMatrix> perms = R_NilValue) {

  int n = mat1.n_rows;
  double N = 0.5 * (double)n * (double)(n - 1);

  // Permutations
  arma::umat P;
  if (perms.isNotNull()) {
    Rcpp::IntegerMatrix pm(perms);
    nperm = pm.ncol();
    P.set_size(n, nperm);
    for (int k = 0; k < nperm; ++k) {
      for (int i = 0; i < n; ++i) {
        P(i, k) = pm(i, k) - 1;
      }
    }
  } else {
    P.set_size(n, nperm);
    for (int k = 0; k < nperm; ++k) {
      arma::uvec idx = arma::regspace<arma::uvec>(0, n - 1);
      for (int i = n - 1; i > 0; --i) {
        int j = (int)(R::unif_rand() * (i + 1)); // Fisher-Yates
        std::swap(idx[i], idx[j]);
      }
      P.col(k) = idx;
    }
  }

  // Parallel for speed
#ifdef _OPENMP
  if ( threads > 1 ) {
    omp_set_num_threads( threads );
  }
#endif

  // Sums that do not change with the permutations
  double s1 = 0, s11 = 0, s2 = 0, s22 = 0;

#pragma omp parallel for reduction(+:s1,s11,s2,s22) schedule(dynamic, 64)
  for (int j = 1; j < n; ++j) {
    for (int i = 0; i < j; ++i) {
      s1 += mat1(i, j); s11 += mat1(i, j) * mat1(i, j);
      s2 += mat2(i, j); s22 += mat2(i, j) * mat2(i, j);
    }
  }

  double den = std::sqrt((s11 - s1 * s1 / N) * (s22 - s2 * s2 / N));

  // Observed Mantel statistic
  arma::uvec id = arma::regspace<arma::uvec>(0, n - 1);
  double observed_stat = (cross_sum(mat1, mat2, id) - s1 * s2 / N) / den;

  // Create a vector to store permutation results
  arma::vec perm_stats(nperm);

  //create progress
  Progress p(nperm, progress);

#pragma omp parallel for schedule(dynamic, 1)
  for (int k = 0; k < nperm; ++k) {

    if (! Progress::check_abort() ) {
      p.increment(); // update progress
    }

    arma::uvec perm = P.col(k);
    perm_stats[k] = (cross_sum(mat1, mat2, perm) - s1 * s2 / N) / den;

  }

  return Rcpp::List::create(
    Rcpp::Named("observed_stat") = observed_stat,
    Rcpp::Named("perm_stats") = perm_stats
  );
}
