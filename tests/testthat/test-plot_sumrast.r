test_that("plot_sumrast works", {
  # Set seed for reproducibility
  set.seed(42, "Mersenne-Twister", sample.kind="Rejection")

  x <- terra::rast(matrix(rbinom(100, 1, 0.2), nrow=10))
  y <- terra::rast(matrix(rbinom(100, 1, 0.8), nrow=10))
  # The above SpatRaster objects have values {0, 1, 2}, so max=2, min=0

  prior3D::plot_sumrast(list(x,y))
  prior3D_plot_sumrast <- recordPlot()
  test_rast <- (x + y)/2
  # Same name
  names(test_rast) <- "sum"
  terra::plot(test_rast)
  test_plot_sumrast <- recordPlot()
  dev.off()
  expect_true(
    all.equal(as.character(prior3D_plot_sumrast)[[3]],
              as.character(test_plot_sumrast)[[3]]
              )
  )
})