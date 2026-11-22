#' Compute Gaussian-kernel MMD for two projected samples
#'
#' Computes the empirical squared maximum mean discrepancy (MMD) between
#' two one-dimensional samples using a Gaussian kernel.
#'
#' This function is intended to be used after the multivariate observations
#' have been projected onto a single random direction.
#'
#' @param x Numeric vector containing the first projected sample.
#' @param y Numeric vector containing the second projected sample.
#' @param bandwidth Positive numeric value giving the Gaussian-kernel
#'   bandwidth. If NULL, a median-distance rule is used.
#'
#' @return A non-negative numeric value containing the empirical
#'   squared MMD statistic.
#'
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(30)
#' y <- stats::rnorm(30, mean = 0.5)
#'
#' rpmmd_mmd(x, y)
#'
#' @export
rpmmd_mmd <- function(x, y, bandwidth = NULL) {
  
  # The projected samples should be numeric vectors with at least
  # two observations in each sample.
  if (!is.numeric(x) || length(x) < 2) {
    stop("`x` must be a numeric vector with at least two observations.")
  }
  
  if (!is.numeric(y) || length(y) < 2) {
    stop("`y` must be a numeric vector with at least two observations.")
  }
  
  # Missing or infinite values are not supported because they would make
  # the pairwise kernel calculations undefined.
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`x` must contain only finite values.")
  }
  
  if (anyNA(y) || any(!is.finite(y))) {
    stop("`y` must contain only finite values.")
  }
  
  # If no bandwidth is supplied, use the median of the positive pairwise
  # distances in the pooled projected sample.
  if (is.null(bandwidth)) {
    
    pooled <- c(x, y)
    distances <- stats::dist(pooled)
    positive_distances <- distances[distances > 0]
    
    # If every projected value is identical, all pairwise distances are zero.
    # In that special case, use bandwidth 1 to keep the kernel well-defined.
    if (length(positive_distances) == 0) {
      bandwidth <- 1
    } else {
      bandwidth <- stats::median(positive_distances)
    }
  }
  
  # The Gaussian bandwidth must be positive.
  if (length(bandwidth) != 1 ||
      !is.numeric(bandwidth) ||
      !is.finite(bandwidth) ||
      bandwidth <= 0) {
    stop("`bandwidth` must be a positive finite number.")
  }
  
  # Helper for computing a Gaussian kernel matrix between two vectors.
  gaussian_kernel <- function(a, b) {
    
    squared_distances <- outer(a, b, "-")^2
    
    exp(
      -squared_distances /
        (2 * bandwidth^2)
    )
  }
  
  # Kernel matrices within and between the two samples.
  Kxx <- gaussian_kernel(x, x)
  Kyy <- gaussian_kernel(y, y)
  Kxy <- gaussian_kernel(x, y)
  
  # Biased empirical squared MMD.
  mmd_squared <-
    mean(Kxx) +
    mean(Kyy) -
    2 * mean(Kxy)
  
  # Numerical roundoff can occasionally produce a tiny negative value,
  # although squared MMD is theoretically non-negative.
  max(mmd_squared, 0)
}