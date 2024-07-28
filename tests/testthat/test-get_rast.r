test_that("get_rast works", {
  expect_equal(
    names(
      prior3D::get_rast(system.file("get_rast_example", package="prior3D"))
    ),
    c("aaptos_aaptos", "abietinaria_abietina")
  )
})