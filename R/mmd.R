
# Internal bandwidth selection function
rpmmd_bandwidth <- function(x, y) {

  distances <- stats::dist(c(x, y))
  positive_distances <- distances[distances > 0]

  if (length(positive_distances) == 0) {
    return(1)
  }

  stats::median(positive_distances)
}

# Internal R reference implementation
mmd_r <- function(x, y, bandwidth) {

  gaussian_kernel <- function(a, b) {
    exp(-0.5 * (outer(a, b, "-") / bandwidth)^2)
  }

  Kxx <- gaussian_kernel(x, x)
  Kyy <- gaussian_kernel(y, y)
  Kxy <- gaussian_kernel(x, y)

  n <- length(x)
  m <- length(y)

  (sum(Kxx) - sum(diag(Kxx))) / (n * (n - 1)) +
    (sum(Kyy) - sum(diag(Kyy))) / (m * (m - 1)) -
    2 * mean(Kxy)
}

#' Compute Gaussian-kernel MMD for two projected samples
#'
#' Computes the unbiased empirical squared maximum mean
#' discrepancy between two one-dimensional samples.
#'
#' @param x Numeric vector containing the first sample.
#' @param y Numeric vector containing the second sample.
#' @param bandwidth Positive finite Gaussian-kernel bandwidth.
#'   If NULL, the median of positive pooled pairwise distances is used.
#' @param engine Computational engine, either "cpp" or "r".
#'
#' @return A numeric value containing the unbiased squared
#'   MMD estimate, which may be negative.
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(30)
#' y <- stats::rnorm(30, mean = 0.5)
#' rpmmd_mmd(x, y)
#'
#' @export
rpmmd_mmd <- function(
    x,
    y,
    bandwidth = NULL,
    engine = c("cpp", "r")
) {

  engine <- match.arg(engine)

  if (!is.numeric(x) || !is.null(dim(x)) || length(x) < 2) {
    stop("`x` must be a numeric vector with at least two observations.")
  }

  if (!is.numeric(y) || !is.null(dim(y)) || length(y) < 2) {
    stop("`y` must be a numeric vector with at least two observations.")
  }

  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Both samples must contain only finite values.")
  }

  if (is.null(bandwidth)) {
    bandwidth <- rpmmd_bandwidth(x, y)
  }

  if (length(bandwidth) != 1 ||
      !is.numeric(bandwidth) ||
      !is.finite(bandwidth) ||
      bandwidth <= 0) {
    stop("`bandwidth` must be a positive finite number.")
  }

  if (engine == "cpp") {
    return(mmd_cpp(x, y, bandwidth))
  }

  mmd_r(x, y, bandwidth)
}

