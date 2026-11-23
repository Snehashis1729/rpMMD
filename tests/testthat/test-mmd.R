test_that("unbiased MMD matches a manually computed example", {
  x <- c(0, 1)
  y <- c(2, 3)
  bandwidth <- 1
  
  k <- function(a, b) {
    exp(-(a - b)^2 / (2 * bandwidth^2))
  }
  
  expected <-
    k(0, 1) +
    k(2, 3) -
    (k(0, 2) + k(0, 3) +
       k(1, 2) + k(1, 3)) / 2
  
  expect_equal(
    rpmmd_mmd(x, y, bandwidth = bandwidth),
    expected,
    tolerance = 1e-12
  )
})

test_that("unbiased squared MMD can be negative", {
  x <- c(0, 1, 2)
  
  result <- rpmmd_mmd(x, x, bandwidth = 1)
  
  expect_true(is.numeric(result))
  expect_length(result, 1)
  expect_lt(result, 0)
})

test_that("MMD is zero when all observations are identical", {
  x <- rep(1, 5)
  y <- rep(1, 5)
  
  expect_equal(
    rpmmd_mmd(x, y),
    0,
    tolerance = 1e-12
  )
})