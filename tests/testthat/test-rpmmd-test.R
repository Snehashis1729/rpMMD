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

test_that("permutation p-value is calculated correctly", {
  
  set.seed(123)
  
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  Y <- matrix(rnorm(30, mean = 0.5), nrow = 10, ncol = 3)
  
  fit <- rpmmd_test(
    X,
    Y,
    statistic = "cvm",
    nproj = 5,
    nperm = 19
  )
  
  expected_p <- (
    1 + sum(fit$perm_statistics >= fit$combined_statistic)
  ) / 20
  
  expect_equal(fit$p_value, expected_p)
  expect_length(fit$perm_statistics, 19)
  expect_gte(fit$p_value, 0)
  expect_lte(fit$p_value, 1)
  
})

test_that("PKS statistic is calculated correctly", {
  
  set.seed(123)
  
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  Y <- matrix(rnorm(30, mean = 0.5), nrow = 10, ncol = 3)
  
  fit <- rpmmd_test(
    X,
    Y,
    statistic = "ks",
    nproj = 5,
    nperm = 19
  )
  
  expected_statistic <- sqrt(
    max(0, max(fit$mmd_by_projection))
  )
  
  expect_equal(
    fit$combined_statistic,
    expected_statistic
  )
  
  expect_gte(fit$combined_statistic, 0)
  expect_length(fit$perm_statistics, 19)
  
})
