#' Run the initial random-projection MMD workflow
#'
#' Applies the current rpMMD workflow to two multivariate samples.
#' Random projection directions are generated, both samples are projected
#' onto those directions, and an MMD value is computed for each projection.
#'
#' The final combined test statistic and permutation-based p-value will be
#' added in a later stage of the project.
#'
#' @param X Numeric matrix containing the first multivariate sample.
#'   Observations are stored in rows and variables in columns.
#' @param Y Numeric matrix containing the second multivariate sample.
#'   Observations are stored in rows and variables in columns.
#' @param statistic Character string specifying the planned combined
#'   statistic. Currently either "cvm" or "ks".
#' @param nproj Positive integer giving the number of random projections.
#' @param bandwidth Optional positive bandwidth passed to \code{rpmmd_mmd()}.
#'   If NULL, the default bandwidth rule in \code{rpmmd_mmd()} is used.
#'
#' @return A list containing the requested statistic type, random projection
#'   directions, projection-specific MMD values, and placeholders for the
#'   final combined statistic and p-value.
#'
#' @examples
#' set.seed(1)
#'
#' X <- matrix(stats::rnorm(60), nrow = 20, ncol = 3)
#' Y <- matrix(stats::rnorm(60, mean = 0.3), nrow = 20, ncol = 3)
#'
#' rpmmd_test(X, Y, nproj = 10)
#'
#' @export
rpmmd_test <- function(
    X,
    Y,
    statistic = c("cvm", "ks"),
    nproj = 100,
    bandwidth = NULL
) {
  
  # X and Y must be numeric matrices because the method works with
  # multivariate observations stored row-wise.
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.")
  }
  
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("`Y` must be a numeric matrix.")
  }
  
  # Both samples must have the same number of variables so that the same
  # random projection directions can be applied to each sample.
  if (ncol(X) != ncol(Y)) {
    stop("`X` and `Y` must have the same number of columns.")
  }
  
  # Match the requested statistic to one of the two planned choices.
  statistic <- match.arg(statistic)
  
  # Check the number of random projection directions.
  if (length(nproj) != 1 ||
      !is.numeric(nproj) ||
      nproj < 1 ||
      nproj != as.integer(nproj)) {
    stop("`nproj` must be a positive integer.")
  }
  
  # Generate random unit directions in the same dimension as the data.
  directions <- rpmmd_directions(
    p = ncol(X),
    nproj = nproj
  )
  
  # Project both multivariate samples onto all generated directions.
  projected <- rpmmd_project(
    X = X,
    Y = Y,
    directions = directions
  )
  
  # Compute one MMD value for each random projection.
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
  
  # At this stage, the projection-specific MMD values are available.
  # The exact CvM/KS aggregation and permutation p-value will be added later.
  list(
    statistic = statistic,
    nproj = nproj,
    directions = directions,
    mmd_by_projection = mmd_by_projection,
    combined_statistic = NA_real_,
    p_value = NA_real_
  )
}