#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppEigen.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <cmath>
#include <vector>
#include <algorithm>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppProgress)]]

// Mantel test with Moran Spectral Randomization (MSR, "pair" method).
//
// Numerically equivalent to adespatial::msr.mantelrtest (Crabot et al. 2019)
// using adespatial::msr.default (Wagner & Dray 2015), but avoids the PCoA and
// the n x n matrix of randomized distances:
//
//   adespatial: X = pcoscaled(d1); Xnew = E[, pairs] %*% a; sim = cor(dist(Xnew), d2)
//   here:       G = -0.5 J (d1^2) J; H = E' G E; Gnew = F H F'; d = sqrt(gii + gjj - 2gij)

typedef Eigen::Map<Eigen::MatrixXd> MapMat;

// Accumulate sums for the Pearson correlation between d (from the Gram matrix)
// and each d2 over the upper triangle of one block of rows (i) and columns (j)
template <typename Scalar>
void block_sums(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& Gb,
                const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& g,
                const std::vector<MapMat>& d2,
                int i0, int j0, double tol,
                double& sd, double& sdd, std::vector<double>& sdy, long& nneg, double& dmin) {

  int ni = Gb.rows();
  int nj = Gb.cols();
  int nk = d2.size();

  for (int jj = 0; jj < nj; ++jj) {
    int j = j0 + jj;
    int imax = std::min(ni, j - i0); // only i < j
    for (int ii = 0; ii < imax; ++ii) {
      int i = i0 + ii;
      double d = (double)g[i] + (double)g[j] - 2.0 * (double)Gb(ii, jj);
      if (d < -tol) {
        ++nneg;
      }
      if (d < dmin) {
        dmin = d;
      }
      d = d > 0 ? std::sqrt(d) : 0;
      sd += d;
      sdd += d * d;
      for (int k = 0; k < nk; ++k) {
        sdy[k] += d * d2[k](i, j);
      }
    }
  }
}

template <typename Scalar>
Rcpp::List mantel_msr_impl(const MapMat& d1,
                           const std::vector<MapMat>& d2,
                           const MapMat& mem,
                           const Rcpp::IntegerVector& idx0,
                           const MapMat& angles,
                           const Rcpp::NumericVector& signs,
                           int threads,
                           int block,
                           double tol_rel,
                           bool progress) {

  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Mat;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vec;

  int n = d1.rows();
  int m = mem.cols();
  int nk = d2.size();
  int nrepet = idx0.size();
  int npair = m / 2;

#ifdef _OPENMP
  if (threads > 0) {
    omp_set_num_threads(threads);
  }
#endif
  Eigen::setNbThreads(threads);

  //----------------------------------------------------------------------------
  // Observed statistics and fixed sums of each d2 (upper triangle)

  double N = 0.5 * (double)n * (double)(n - 1);
  double s1 = 0, s11 = 0;
  std::vector<double> s2(nk, 0), s22(nk, 0), s12(nk, 0);

#pragma omp parallel
{
  double t1 = 0, t11 = 0;
  std::vector<double> t2(nk, 0), t22(nk, 0), t12(nk, 0);

#pragma omp for schedule(dynamic, 64)
  for (int j = 1; j < n; ++j) {
    for (int i = 0; i < j; ++i) {
      double a = d1(i, j);
      t1 += a; t11 += a * a;
      for (int k = 0; k < nk; ++k) {
        double b = d2[k](i, j);
        t2[k] += b; t22[k] += b * b;
        t12[k] += a * b;
      }
    }
  }

#pragma omp critical
{
  s1 += t1; s11 += t11;
  for (int k = 0; k < nk; ++k) {
    s2[k] += t2[k]; s22[k] += t22[k]; s12[k] += t12[k];
  }
}
}

  Eigen::VectorXd var2(nk), obs(nk);
  for (int k = 0; k < nk; ++k) {
    var2[k] = s22[k] - s2[k] * s2[k] / N;
    obs[k] = (s12[k] - s1 * s2[k] / N) / std::sqrt((s11 - s1 * s1 / N) * var2[k]);
  }

  //----------------------------------------------------------------------------
  // Gower centred matrix of d1 and its projection on the MEMs (H = E'GE)

  Mat E = mem.cast<Scalar>();
  Mat H;
  double tol, gmean;
  {
    Eigen::MatrixXd A = -0.5 * d1.array().square().matrix();
    Eigen::VectorXd rm = A.rowwise().mean();
    double gm = rm.mean();
    Mat G(n, n);

#pragma omp parallel for
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        G(i, j) = (Scalar)(A(i, j) - rm[i] - rm[j] + gm);
      }
    }

    // Negative squared distances beyond this are not rounding errors
    gmean = G.diagonal().template cast<double>().mean();
    tol = tol_rel * gmean;

    A.resize(0, 0);
    Mat GE = G * E;
    G.resize(0, 0);
    H.noalias() = E.transpose() * GE;
  }

  //----------------------------------------------------------------------------
  // Blocks of the upper triangle

  std::vector<int> bi, bj;
  int nblock = (n + block - 1) / block;
  for (int b2 = 0; b2 < nblock; ++b2) {
    for (int b1 = 0; b1 <= b2; ++b1) {
      bi.push_back(b1);
      bj.push_back(b2);
    }
  }
  int npairs_block = bi.size();

  //----------------------------------------------------------------------------
  // Replicates

  Eigen::MatrixXd sim(nrepet, nk);
  Mat F(n, m);
  Mat T(n, m);
  Vec g(n);
  std::vector<int> order(m);
  double dmin_all = 0;

  Progress p(nrepet, progress);

  for (int r = 0; r < nrepet; ++r) {

    if (Progress::check_abort()) {
      Rcpp::stop("Interrupted");
    }

    // Pairs of MEMs as in adespatial:::cutbypair (1-based idx0 from R)
    int c0 = 0;
    for (int c = 0; c < m; ++c) {
      if ((m % 2 == 0) || (c != idx0[r] - 1)) {
        order[c0] = c;
        ++c0;
      }
    }

    // F: rotated pairs of MEMs (as in adespatial:::genbypair)
    for (int q = 0; q < npair; ++q) {
      int a = order[2 * q];
      int b = order[2 * q + 1];
      Scalar cs = (Scalar)std::cos(angles(q, r));
      Scalar sn = (Scalar)std::sin(angles(q, r));
      F.col(a) = cs * E.col(a) + sn * E.col(b);
      F.col(b) = -sn * E.col(a) + cs * E.col(b);
    }

    // Singleton MEM with random sign (as in adespatial:::genbysingleton)
    if (m % 2 == 1) {
      int a = idx0[r] - 1;
      F.col(a) = (Scalar)signs[r] * E.col(a);
    }

    // Gram matrix of the randomized coordinates: Gnew = F H F'
    T.noalias() = F * H;

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
      g[i] = T.row(i).dot(F.row(i));
    }

    double sd = 0, sdd = 0;
    std::vector<double> sdy(nk, 0);
    long nneg = 0;
    double dmin = 0;

#pragma omp parallel
{
    double td = 0, tdd = 0;
    std::vector<double> tdy(nk, 0);
    long tneg = 0;
    double tmin = 0;

#pragma omp for schedule(dynamic, 1)
    for (int b = 0; b < npairs_block; ++b) {
      int i0 = bi[b] * block;
      int j0 = bj[b] * block;
      int ni = std::min(block, n - i0);
      int nj = std::min(block, n - j0);
      Mat Gb = T.middleRows(i0, ni) * F.middleRows(j0, nj).transpose();
      block_sums<Scalar>(Gb, g, d2, i0, j0, tol, td, tdd, tdy, tneg, tmin);
    }

#pragma omp critical
{
    sd += td; sdd += tdd; nneg += tneg;
    dmin = std::min(dmin, tmin);
    for (int k = 0; k < nk; ++k) {
      sdy[k] += tdy[k];
    }
}
}

    dmin_all = std::min(dmin_all, dmin);

    if (nneg > 0) {
      Rcpp::stop("m1 is not Euclidean: replicate %d has %d negative squared distances", r + 1, (int)nneg);
    }

    for (int k = 0; k < nk; ++k) {
      sim(r, k) = (sdy[k] - sd * s2[k] / N) / std::sqrt((sdd - sd * sd / N) * var2[k]);
    }

    p.increment();
  }

  return Rcpp::List::create(
    Rcpp::Named("obs") = obs,
    Rcpp::Named("sim") = sim,
    Rcpp::Named("min_d2") = dmin_all / gmean
  );
}

// [[Rcpp::export]]
Rcpp::List mantel_msr_rcpp(const Eigen::Map<Eigen::MatrixXd> d1,
                           const Rcpp::List d2,
                           const Eigen::Map<Eigen::MatrixXd> mem,
                           const Rcpp::IntegerVector idx0,
                           const Eigen::Map<Eigen::MatrixXd> angles,
                           const Rcpp::NumericVector signs,
                           int threads = 1,
                           int block = 2048,
                           bool single = true,
                           double tol = 1e-3,
                           bool progress = true) {

  int n = d1.rows();
  std::vector<MapMat> d2maps;
  for (int k = 0; k < d2.size(); ++k) {
    Rcpp::NumericMatrix x = d2[k];
    if (x.nrow() != n || x.ncol() != n) {
      Rcpp::stop("All matrices must have the same dimensions");
    }
    d2maps.push_back(MapMat(x.begin(), n, n));
  }

  if (single) {
    return mantel_msr_impl<float>(d1, d2maps, mem, idx0, angles, signs, threads, block, tol, progress);
  } else {
    return mantel_msr_impl<double>(d1, d2maps, mem, idx0, angles, signs, threads, block, tol, progress);
  }
}
