#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>
#include <cmath>


// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Function to calculate Euclidean distance between two points
inline double euclidean_distance(const NumericVector& p1, const NumericVector& p2) {
  double sum = 0.0;
  int n = p1.size();
  for (int i = 0; i < n; ++i) {
    double diff = p1[i] - p2[i];
    sum += diff * diff;
  }
  return sqrt(sum);
}

// Function to calculate Euclidean distance matrix using OpenMP for parallelization
// [[Rcpp::export]]
NumericMatrix euclidean_distance_matrix(const NumericMatrix& mat, int threads = 1) {
  int n = mat.nrow();
  NumericMatrix dist_mat(n, n);
  
  // Parallel for speed
#ifdef _OPENMP
  if ( threads > 0 ) {
    omp_set_num_threads( threads );
  }
#endif
  
  // Parallelizing the computation of the distance matrix
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      double dist = euclidean_distance(mat.row(i), mat.row(j));
      dist_mat(i, j) = dist;
      dist_mat(j, i) = dist; // Symmetric matrix
    }
  }

  return dist_mat;
}
