test_that("sumrast works", {
  # Set seed for reproducibility
  set.seed(42, "Mersenne-Twister", sample.kind="Rejection")
  
  x <- terra::rast(matrix(rbinom(100, 1, 0.2), nrow=10))
  y <- terra::rast(matrix(rbinom(100, 1, 0.8), nrow=10))
  # The above SpatRaster objects have values {0, 1, 2}, so max=2, min=0
  
  # Same name
  test_rast <- (x + y)/2
  names(test_rast) <- "sum"
  expect_true(terra::all.equal(prior3D::sumrast(list(x,y)), test_rast))
})