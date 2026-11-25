
#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// Gaussian kernel
inline double gaussian_kernel(
    double x, double y, double bandwidth) {

  double z = (x - y) / bandwidth;

  return std::exp(-0.5 * z * z);
}

// [[Rcpp::export]]
double mmd_cpp(
    NumericVector x,
    NumericVector y,
    double bandwidth) {

  R_xlen_t n = x.size();
  R_xlen_t m = y.size();

  if (n < 2 || m < 2) {
    stop("Both samples must have at least two observations.");
  }

  if (!std::isfinite(bandwidth) || bandwidth <= 0) {
    stop("Bandwidth must be positive and finite.");
  }

  double sum_xx = 0.0;
  double sum_yy = 0.0;
  double sum_xy = 0.0;

  // Within-sample kernel sum for X
  for (R_xlen_t i = 0; i < n; ++i) {
    for (R_xlen_t j = i + 1; j < n; ++j) {
      sum_xx += gaussian_kernel(x[i], x[j], bandwidth);
    }
  }

  // Within-sample kernel sum for Y
  for (R_xlen_t i = 0; i < m; ++i) {
    for (R_xlen_t j = i + 1; j < m; ++j) {
      sum_yy += gaussian_kernel(y[i], y[j], bandwidth);
    }
  }

  // Between-sample kernel sum
  for (R_xlen_t i = 0; i < n; ++i) {
    for (R_xlen_t j = 0; j < m; ++j) {
      sum_xy += gaussian_kernel(x[i], y[j], bandwidth);
    }
  }

  double nd = static_cast<double>(n);
  double md = static_cast<double>(m);

  // Unbiased squared MMD
  return 2.0 * sum_xx / (nd * (nd - 1.0)) +
         2.0 * sum_yy / (md * (md - 1.0)) -
         2.0 * sum_xy / (nd * md);
}

