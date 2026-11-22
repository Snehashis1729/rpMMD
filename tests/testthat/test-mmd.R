test_that("MMD is zero for identical samples", {
  set.seed(1)
  
  x <- stats::rnorm(20)
  
  expect_equal(
    rpmmd_mmd(x, x),
    0,
    tolerance = 1e-12
  )
})

test_that("MMD is non-negative", {
  set.seed(1)
  
  x <- stats::rnorm(20)
  y <- stats::rnorm(20, mean = 0.5)
  
  expect_gte(
    rpmmd_mmd(x, y),
    0
  )
})
