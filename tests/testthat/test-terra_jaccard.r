test_that("terra_jaccard works", {
  # Set seed for reproducibility
  set.seed(42, "Mersenne-Twister", sample.kind="Rejection")
  
  x <- terra::rast(matrix(rbinom(100, 1, 0.2), nrow=10))
  y <- terra::rast(matrix(rbinom(100, 1, 0.8), nrow=10))
  expect_equal(prior3D::terra_jaccard(x, y), 0.24705882352941177516)
})