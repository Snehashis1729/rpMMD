#' Generate random projection directions
#'
#' Generates random unit-length directions in p dimensions.
#' Each row of the returned matrix represents one projection direction.
#'
#' @param p Positive integer giving the dimension of the data.
#' @param nproj Positive integer giving the number of random projections.
#'
#' @return A numeric matrix with nproj rows and p columns.
#'   Each row has Euclidean norm equal to one.
#'
#' @examples
#' set.seed(1)
#' rpmmd_directions(p = 3, nproj = 5)
#'
#' @export
rpmmd_directions <- function(p, nproj = 100) {
  
  # Check that p is a positive integer.
  if (length(p) != 1 || !is.numeric(p) || p < 1 || p != as.integer(p)) {
    stop("`p` must be a positive integer.")
  }
  
  # Check that nproj is a positive integer.
  if (length(nproj) != 1 ||
      !is.numeric(nproj) ||
      nproj < 1 ||
      nproj != as.integer(nproj)) {
    stop("`nproj` must be a positive integer.")
  }
  
  # Generate nproj random Gaussian vectors in p dimensions.
  directions <- matrix(
  stats::rnorm(nproj * p),
    nrow = nproj,
    ncol = p
  )
  
  # Normalize every row so that each direction has length one.
  norms <- sqrt(rowSums(directions^2))
  directions <- directions / norms
  
  directions
}