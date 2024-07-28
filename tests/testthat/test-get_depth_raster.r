test_that("get_depth_raster works", {
  expect_equal(names(prior3D::get_depth_raster()), "depth_raster")
})