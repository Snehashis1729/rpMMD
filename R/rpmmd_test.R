
#' Random-projection MMD two-sample test
#'
#' Tests equality of two multivariate distributions using
#' random projections and Gaussian-kernel MMD statistics.
#'
#' @param X Numeric matrix containing the first sample.
#'   Observations are stored in rows.
#' @param Y Numeric matrix containing the second sample.
#'   Observations are stored in rows.
#' @param statistic Combined statistic: "cvm" or "ks".
#' @param nproj Positive integer giving the number of projections.
#' @param bandwidth Positive finite bandwidth or NULL.
#'   If NULL, a pooled median bandwidth is selected separately
#'   for each projection and reused across permutations.
#' @param nperm Positive integer giving the number of permutations.
#' @param engine Computational engine: "cpp" or "r".
#'
#' @details Projection directions and bandwidths remain fixed
#'   across permutations. A small floating-point tolerance is
#'   used when comparing permutation and observed statistics,
#'   on the squared scale for the KS statistic.
#'
#' @return An object of class rpmmd_result containing projection directions, bandwidths,
#'   projection-specific MMD values, combined statistic,
#'   permutation statistics and permutation p-value.
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
    nperm = 999,
    engine = c("cpp", "r")
) {

  statistic <- match.arg(statistic)
  engine <- match.arg(engine)

  # Validate the input matrices.
  if (!is.matrix(X) || !is.numeric(X) || is.complex(X)) {
    stop("`X` must be a real numeric matrix.")
  }

  if (!is.matrix(Y) || !is.numeric(Y) || is.complex(Y)) {
    stop("`Y` must be a real numeric matrix.")
  }

  if (ncol(X) != ncol(Y)) {
    stop("Both samples must have the same number of columns.")
  }

  if (nrow(X) < 2 || nrow(Y) < 2) {
    stop("Both samples must contain at least two observations.")
  }

  if (ncol(X) < 1) {
    stop("Both samples must contain at least one variable.")
  }

  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("Both samples must contain only finite values.")
  }

  # Validate the number of projections and permutations.
  positive_integer <- function(z) {
    length(z) == 1 &&
      is.numeric(z) &&
      !is.complex(z) &&
      is.finite(z) &&
      z >= 1 &&
      z <= .Machine$integer.max &&
      z == floor(z)
  }

  if (!positive_integer(nproj)) {
    stop("`nproj` must be a positive integer.")
  }

  if (!positive_integer(nperm)) {
    stop("`nperm` must be a positive integer.")
  }

  if (!is.null(bandwidth)) {

    if (length(bandwidth) != 1 ||
        !is.numeric(bandwidth) ||
        is.complex(bandwidth) ||
        !is.finite(bandwidth) ||
        bandwidth <= 0) {
      stop("`bandwidth` must be a positive finite number.")
    }
  }

  n <- nrow(X)
  m <- nrow(Y)

  # Generate random directions only once.
  directions <- rpmmd_directions(ncol(X), nproj)

  # Project all observations only once.
  projected <- rpmmd_project(X, Y, directions)

  Z_proj <- rbind(
    projected$X_proj,
    projected$Y_proj
  )

  if (any(!is.finite(Z_proj))) {
    stop("Projected observations must be finite.")
  }

  ix <- seq_len(n)
  iy <- n + seq_len(m)

  # Compute the bandwidth for each direction only once.
  bandwidths <- if (is.null(bandwidth)) {

    vapply(
      seq_len(nproj),
      function(j) {
        rpmmd_bandwidth(
          Z_proj[ix, j],
          Z_proj[iy, j]
        )
      },
      numeric(1)
    )

  } else {

    rep(bandwidth, nproj)
  }

  # Compute MMD for a given partition of the pooled sample.
  compute_mmd <- function(index_x, index_y) {

    vapply(
      seq_len(nproj),
      function(j) {
        rpmmd_mmd(
          Z_proj[index_x, j],
          Z_proj[index_y, j],
          bandwidth = bandwidths[j],
          engine = engine
        )
      },
      numeric(1)
    )
  }

  # Aggregate the projection-specific MMD estimates.
  aggregate_mmd <- function(values) {

    if (statistic == "cvm") {
      return(mean(values))
    }

    sqrt(max(0, max(values)))
  }

  # Observed statistic.
  mmd_by_projection <- compute_mmd(ix, iy)

  combined_statistic <- aggregate_mmd(
    mmd_by_projection
  )

  # Permutation calibration.
  perm_statistics <- numeric(nperm)

  for (b in seq_len(nperm)) {

    indices <- sample.int(n + m)

    index_x <- sort(indices[ix])
    index_y <- sort(indices[iy])

    perm_statistics[b] <- aggregate_mmd(
      compute_mmd(index_x, index_y)
    )
  }

  # Compare statistics with a small numerical tolerance.
  observed <- combined_statistic
  permuted <- perm_statistics

  if (statistic == "ks") {
    observed <- observed^2
    permuted <- permuted^2
  }

  tolerance <- 64 * .Machine$double.eps *
    max(1, abs(observed))

  p_value <- (
    1 + sum(permuted >= observed - tolerance)
  ) / (nperm + 1)

  # Return test results.
  new_rpmmd_result(list(
    statistic = statistic,
    nproj = nproj,
    nperm = nperm,
    directions = directions,
    bandwidths = bandwidths,
    mmd_by_projection = mmd_by_projection,
    combined_statistic = combined_statistic,
    perm_statistics = perm_statistics,
    p_value = p_value,
    engine = engine,
    n = n,
    m = m,
    dimension = ncol(X)
  ))
}

