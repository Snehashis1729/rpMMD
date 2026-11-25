
test_that("R and C++ MMD implementations agree", {

  set.seed(100)

  x <- rnorm(25)
  y <- rnorm(35, mean = 0.5)

  for (h in c(0.1, 0.5, 1, 2, 5)) {

    expect_equal(
      rpmmd_mmd(x, y, h, engine = "cpp"),
      rpmmd_mmd(x, y, h, engine = "r"),
      tolerance = 1e-12
    )
  }
})

test_that("Automatic bandwidth gives equivalent results", {

  set.seed(101)

  x <- rnorm(20)
  y <- rnorm(30)

  expect_equal(
    rpmmd_mmd(x, y, engine = "cpp"),
    rpmmd_mmd(x, y, engine = "r"),
    tolerance = 1e-12
  )
})

test_that("Bandwidth is invariant to group labels", {

  expect_equal(
    rpmmd_bandwidth(c(0, 1, 3), c(2, 4, 5)),
    rpmmd_bandwidth(c(0, 2, 4), c(1, 3, 5))
  )

  expect_equal(
    rpmmd_bandwidth(rep(1, 3), rep(1, 4)),
    1
  )
})

test_that("Complete R and C++ permutation tests agree", {

  set.seed(102)

  X <- matrix(rnorm(36), ncol = 3)
  Y <- matrix(rnorm(45, mean = 0.4), ncol = 3)

  for (method in c("cvm", "ks")) {

    set.seed(103)

    fit_r <- rpmmd_test(
      X, Y,
      statistic = method,
      nproj = 5,
      nperm = 19,
      engine = "r"
    )

    set.seed(103)

    fit_cpp <- rpmmd_test(
      X, Y,
      statistic = method,
      nproj = 5,
      nperm = 19,
      engine = "cpp"
    )

    expect_equal(
      fit_r$combined_statistic,
      fit_cpp$combined_statistic,
      tolerance = 1e-12
    )

    expect_equal(
      fit_r$perm_statistics,
      fit_cpp$perm_statistics,
      tolerance = 1e-12
    )

    expect_equal(
      fit_r$p_value,
      fit_cpp$p_value
    )

    expect_equal(
      fit_r$directions,
      fit_cpp$directions
    )
  }
})

test_that("Test results have valid S3 structure", {

  set.seed(104)

  X <- matrix(rnorm(30), ncol = 3)
  Y <- matrix(rnorm(30), ncol = 3)

  fit <- rpmmd_test(
    X, Y,
    nproj = 5,
    nperm = 19
  )

  expect_s3_class(fit, "rpmmd_result")

  expect_length(fit$bandwidths, 5)

  expect_true(all(fit$bandwidths > 0))

  expect_output(print(fit), "Permutation p-value")

  expect_type(summary(fit), "list")

  expect_equal(
    summary(fit)$p_value,
    fit$p_value
  )

  expect_equal(
    fit$combined_statistic,
    mean(fit$mmd_by_projection)
  )
})

test_that("Validator rejects invalid results", {

  X <- matrix(1:12, ncol = 3)
  Y <- matrix(2:13, ncol = 3)

  fit <- rpmmd_test(
    X, Y,
    nproj = 3,
    nperm = 9
  )

  fit$p_value <- -1

  expect_error(validate_rpmmd_result(fit))
})

test_that("Identical constant samples give p-value one", {

  X <- matrix(1, nrow = 5, ncol = 2)
  Y <- matrix(1, nrow = 6, ncol = 2)

  fit <- rpmmd_test(
    X, Y,
    nproj = 3,
    nperm = 19
  )

  expect_equal(fit$combined_statistic, 0)
  expect_equal(fit$p_value, 1)
})

test_that("MMD rejects invalid inputs", {

  expect_error(rpmmd_mmd(1, c(1, 2)))
  expect_error(rpmmd_mmd(c(1, NA), c(1, 2)))
  expect_error(rpmmd_mmd(c(1, 2), c(1, 2), -1))
  expect_error(rpmmd_mmd(c(1, 2), c(1, 2), Inf))
})

test_that("Main test rejects invalid arguments", {

  X <- matrix(1:12, ncol = 3)
  Y <- matrix(2:13, ncol = 3)

  expect_error(rpmmd_test(X, Y, nproj = 0))
  expect_error(rpmmd_test(X, Y, nperm = 0))
  expect_error(rpmmd_test(X, Y, nproj = NA))
  expect_error(rpmmd_test(X, Y, bandwidth = -1))
  expect_error(rpmmd_test(X, Y[, 1:2]))
  expect_error(rpmmd_test(X[1, , drop = FALSE], Y))
  expect_error(rpmmd_test(X, Y, statistic = "invalid"))
})

