#' Project two multivariate samples onto specified directions
#'
#' Projects the rows of two multivariate samples onto one or more supplied
#' projection directions. Each observation is stored in a row, and each
#' variable is stored in a column.
#'
#' @param X Numeric matrix for the first sample, with observations in rows.
#' @param Y Numeric matrix for the second sample, with observations in rows.
#' @param directions Numeric matrix whose rows contain the projection
#'   directions.
#'
#' @return A list with components \code{X_proj} and \code{Y_proj}.
#'   Each column corresponds to one projection direction.
#'
#' @examples
#' set.seed(1)
#'
#' X <- matrix(stats::rnorm(30), nrow = 10, ncol = 3)
#' Y <- matrix(stats::rnorm(24), nrow = 8, ncol = 3)
#'
#' U <- rpmmd_directions(p = 3, nproj = 5)
#' rpmmd_project(X, Y, U)
#'
#' @export
rpmmd_project <- function(X, Y, directions) {
  
  # X and Y should both be numeric matrices because the projection is based
  # on matrix multiplication.
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.")
  }
  
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("`Y` must be a numeric matrix.")
  }
  
  # The two samples must have the same number of variables so that they live
  # in the same p-dimensional space.
  if (ncol(X) != ncol(Y)) {
    stop("`X` and `Y` must have the same number of columns.")
  }
  
  # The projection directions are stored row-wise in a numeric matrix.
  if (!is.matrix(directions) || !is.numeric(directions)) {
    stop("`directions` must be a numeric matrix.")
  }
  
  # Each projection direction must have one coordinate for each variable in
  # X and Y.
  if (ncol(directions) != ncol(X)) {
    stop(
      "`directions` must have the same number of columns as `X` and `Y`."
    )
  }
  
  # If U contains the projection directions in its rows, then
  # X %*% t(U) gives all projections of X at once. Column j of X_proj
  # contains the projection of every observation onto direction j.
  X_proj <- X %*% t(directions)
  Y_proj <- Y %*% t(directions)
  
  # Return both projected samples so that later steps can compute
  # one-dimensional MMD statistics for each direction.
  list(
    X_proj = X_proj,
    Y_proj = Y_proj
  )
}