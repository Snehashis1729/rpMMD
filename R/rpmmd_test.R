
#' Random-projection MMD two-sample test
#'
#' Tests whether two multivariate samples come from the same distribution
#' using random projections and Gaussian-kernel MMD statistics.
#' The projection-specific statistics are aggregated using either a
#' Cramer-von Mises (CvM) or Kolmogorov-Smirnov (KS) type statistic.
#' A permutation test is used to calculate the p-value.
#'
#' @param X Numeric matrix containing the first multivariate sample.
#'   Observations are stored in rows and variables in columns.
#' @param Y Numeric matrix containing the second multivariate sample.
#'   Observations are stored in rows and variables in columns.
#' @param statistic Character string specifying the combined statistic.
#'   Either "cvm" or "ks".
#' @param nproj Positive integer giving the number of random projections.
#' @param bandwidth Optional positive bandwidth passed to \code{rpmmd_mmd()}.
#'   If NULL, the default bandwidth rule in \code{rpmmd_mmd()} is used.
#' @param nperm Positive integer giving the number of random permutations.
#'
#' @return A list containing the statistic type, number of projections,
#'   projection directions, projection-specific MMD values, combined
#'   statistic, permutation statistics, and permutation p-value.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(stats::rnorm(60), nrow = 20, ncol = 3)
#' Y <- matrix(stats::rnorm(60, mean = 0.3), nrow = 20, ncol = 3)
#' rpmmd_test(X, Y, nproj = 10, nperm = 19)
#'
#' @export
rpmmd_test <- function(
    X,
    Y,
    statistic = c("cvm", "ks"),
    nproj = 100,
    bandwidth = NULL,
    nperm = 999
) {
  
  # Validate the two multivariate samples.
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.")
  }
  
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("`Y` must be a numeric matrix.")
  }
  
  if (ncol(X) != ncol(Y)) {
    stop("`X` and `Y` must have the same number of columns.")
  }
  
  # The unbiased MMD estimator requires at least two observations
  # in each sample.
  if (nrow(X) < 2 || nrow(Y) < 2) {
    stop("Both samples must contain at least two observations.")
  }
  
  if (ncol(X) < 1) {
    stop("Both samples must contain at least one variable.")
  }
  
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("`X` and `Y` must contain only finite values.")
  }
  
  statistic <- match.arg(statistic)
  
  # Validate the number of random projections.
  if (length(nproj) != 1 ||
      !is.numeric(nproj) ||
      is.na(nproj) ||
      !is.finite(nproj) ||
      nproj < 1 ||
      nproj != floor(nproj)) {
    stop("`nproj` must be a positive integer.")
  }
  
  # Validate the number of permutations.
  if (length(nperm) != 1 ||
      !is.numeric(nperm) ||
      is.na(nperm) ||
      !is.finite(nperm) ||
      nperm < 1 ||
      nperm != floor(nperm)) {
    stop("`nperm` must be a positive integer.")
  }
  
  # Generate random unit directions.
  directions <- rpmmd_directions(
    p = ncol(X),
    nproj = nproj
  )
  
  # Project the original samples.
  projected <- rpmmd_project(
    X = X,
    Y = Y,
    directions = directions
  )
  
  # Compute MMD for each original projection.
  mmd_by_projection <- vapply(
    seq_len(nproj),
    function(j) {
      rpmmd_mmd(
        x = projected$X_proj[, j],
        y = projected$Y_proj[, j],
        bandwidth = bandwidth
      )
    },
    numeric(1)
  )
  
  # Combine projection-specific MMD values.
  if (statistic == "cvm") {
    combined_statistic <- mean(mmd_by_projection)
  } else {
    combined_statistic <- sqrt(max(0, max(mmd_by_projection)))
  }
  
  # Pool the observations for permutation testing.
  n <- nrow(X)
  m <- nrow(Y)
  Z <- rbind(X, Y)
  
  # Allocate storage for permutation statistics.
  perm_statistics <- numeric(nperm)
  
  # Perform permutation testing.
  for (b in seq_len(nperm)) {
    
    # Randomly permute the pooled observations.
    perm_indices <- sample.int(n + m)
    
    # Split the permuted observations into two groups.
    X_perm <- Z[perm_indices[seq_len(n)], , drop = FALSE]
    Y_perm <- Z[perm_indices[n + seq_len(m)], , drop = FALSE]
    
    # Use the same projection directions as the original test.
    projected_perm <- rpmmd_project(
      X = X_perm,
      Y = Y_perm,
      directions = directions
    )
    
    # Compute MMD for each permuted projection.
    mmd_perm <- vapply(
      seq_len(nproj),
      function(j) {
        rpmmd_mmd(
          x = projected_perm$X_proj[, j],
          y = projected_perm$Y_proj[, j],
          bandwidth = bandwidth
        )
      },
      numeric(1)
    )
    
    # Aggregate the permutation-specific MMD values.
    if (statistic == "cvm") {
      perm_statistics[b] <- mean(mmd_perm)
    } else {
      perm_statistics[b] <- sqrt(max(0, max(mmd_perm)))
    }
  }
  
  # Monte Carlo permutation p-value.
  p_value <- (
    1 + sum(perm_statistics >= combined_statistic)
  ) / (nperm + 1)
  
  # Return the test results.
  list(
    statistic = statistic,
    nproj = nproj,
    nperm = nperm,
    directions = directions,
    mmd_by_projection = mmd_by_projection,
    combined_statistic = combined_statistic,
    perm_statistics = perm_statistics,
    p_value = p_value
  )
}
