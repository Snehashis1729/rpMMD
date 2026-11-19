test_that("random directions have the correct dimensions", {
  set.seed(1)
  
  U <- rpmmd_directions(p = 3, nproj = 5)
  
  expect_equal(dim(U), c(5, 3))
})

test_that("random directions have unit norm", {
  set.seed(1)
  
  U <- rpmmd_directions(p = 3, nproj = 5)
  
  expect_equal(
    rowSums(U^2),
    rep(1, 5),
    tolerance = 1e-12
  )
})