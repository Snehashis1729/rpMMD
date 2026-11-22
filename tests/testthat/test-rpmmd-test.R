test_that("rpmmd_test returns correct number of projections", {
  
  set.seed(123)
  
  X <- matrix(rnorm(60), nrow = 20, ncol = 3)
  Y <- matrix(rnorm(60), nrow = 20, ncol = 3)
  
  fit <- rpmmd_test(
    X,
    Y,
    nproj = 10
  )
  
  expect_equal(
    length(fit$mmd_by_projection),
    10
  )
  
})

test_that("projected MMD values are finite", {
  
  set.seed(123)
  
  X <- matrix(rnorm(60), nrow = 20, ncol = 3)
  Y <- matrix(rnorm(60), nrow = 20, ncol = 3)
  
  fit <- rpmmd_test(
    X,
    Y,
    nproj = 10
  )
  
  expect_true(
    all(is.finite(fit$mmd_by_projection))
  )
  
})
