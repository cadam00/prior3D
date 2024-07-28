test_that("biodiv_df works", {
  data(biodiv_df, package = "prior3D")
  expect_equal(
    head(biodiv_df),
    data.frame(
      species_name = c("acanthocybium_solandri", "acantholabrus_palloni",
                       "acanthomysis_longicornis", "abraliopsis_morisii",
                       "abralia_veranyi", "abraliopsis_pfefferi"),
      pelagic = c(1, 0, 0, 0, 0, 0),
      min_z = c(-20, -500, -100, -3660, -900, -750),
      max_z = c(0, -30, -2, 0, -1, -1)
    )
  )
})