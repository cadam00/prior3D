test_that("get_biodiv_raster works", {
  expect_equal(
    head(names(prior3D::get_biodiv_raster())),
    c("aaptos_aaptos", "abietinaria_abietina", "abra_alba", "abralia_veranyi",
      "abraliopsis_morisii", "abraliopsis_pfefferi")
  )
})