test_that("projected samples have the correct dimensions", {
  set.seed(1)
  
  X <- matrix(stats::rnorm(30), nrow = 10, ncol = 3)
  Y <- matrix(stats::rnorm(24), nrow = 8, ncol = 3)
  
  U <- rpmmd_directions(p = 3, nproj = 5)
  
  out <- rpmmd_project(X, Y, U)
  
  expect_equal(dim(out$X_proj), c(10, 5))
  expect_equal(dim(out$Y_proj), c(8, 5))
})

test_that("projection values agree with direct matrix multiplication", {
  set.seed(1)
  
  X <- matrix(stats::rnorm(12), nrow = 4, ncol = 3)
  Y <- matrix(stats::rnorm(9), nrow = 3, ncol = 3)
  
  U <- rpmmd_directions(p = 3, nproj = 2)
  
  out <- rpmmd_project(X, Y, U)
  
  expect_equal(out$X_proj, X %*% t(U))
  expect_equal(out$Y_proj, Y %*% t(U))
})