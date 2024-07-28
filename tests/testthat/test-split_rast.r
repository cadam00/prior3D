test_that("split_rast works", {
  library(terra)
  
  biodiv_raster <- prior3D::get_biodiv_raster()
  depth_raster <- prior3D::get_depth_raster()
  data(biodiv_df)
  # You can split features' 2D distributions into
  # 3D ones and then run only 3D analysis
  split_features <- prior3D::split_rast(
                               biodiv_raster,
                               depth_raster,
                               breaks = c(0, -40, -200, -2000, -Inf),
                               biodiv_df
                    )
  expect_length(split_features,4)
  expect_named(
    split_features,
    c("(-40,0]", "(-200,-40]", "(-2e+03,-200]", "[-Inf,-2e+03]")
  )
  for (i in 1:length(split_features)){
    expect_equal(
      head(names(split_features[[i]])),
      c("aaptos_aaptos", "abietinaria_abietina", "abra_alba", "abralia_veranyi",
      "abraliopsis_morisii", "abraliopsis_pfefferi")
    )
  }
  
  expect_length(split_features,4)
  
  split_features2 <- prior3D::split_rast(
                             biodiv_raster,
                             depth_raster,
                             breaks = c(0, -40, -200, -2000, -Inf),
                             data.frame(
                               biodiv_df[,c("species_name", "pelagic")],
                               "min_z" = NA,
                               "max_z" = NA
                             ),
                             val_depth_range = FALSE
                     )
  expect_true(all.equal(split_features, split_features2))
  
  split_features2 <- prior3D::split_rast(
                             biodiv_raster,
                             depth_raster,
                             breaks = c(0, -40, -200, -2000, -Inf),
                             data.frame(
                               biodiv_df[,c("species_name", "pelagic")],
                               "min_z" = NA
                             ),
                             val_depth_range = FALSE
                     )
  expect_true(all.equal(split_features, split_features2))
  
  split_features2 <- prior3D::split_rast(
                             biodiv_raster,
                             depth_raster,
                             breaks = c(0, -40, -200, -2000, -Inf),
                             data.frame(
                               biodiv_df[,c("species_name", "pelagic")],
                               "max_z" = NA
                             ),
                             val_depth_range = FALSE
                     )
  expect_true(all.equal(split_features, split_features2))
  
  split_features2 <- prior3D::split_rast(
                                biodiv_raster,
                                depth_raster,
                             breaks = c(0, -40, -200, -2000, -Inf),
                             biodiv_df[,c("species_name", "pelagic")],
                             val_depth_range = FALSE
                     )
  expect_true(all.equal(split_features, split_features2))
  
  split_features2 <- prior3D::split_rast(
                       biodiv_raster = system.file("get_rast_example",
                                                          package="prior3D"),
                       depth_raster,
                       breaks = c(0, -40, -200, -2000, -Inf),
                       biodiv_df
                     )
  expect_equal(nlyr(split_features2[[1]]), 2)
  
  split_features2 <- prior3D::split_rast(
                              biodiv_raster = biodiv_raster,
                              depth_raster  =
                                system.file("external/depth_raster.tif",
                                       package="prior3D"),
                           breaks = c(0, -40, -200, -2000, -Inf),
                           biodiv_df =
                              system.file("test_datasets/biodiv_df.csv",
                                          package="prior3D")
                   )
  expect_true(all.equal(split_features, split_features2))
  
  split_features <- suppressWarnings(prior3D::split_rast(
    biodiv_raster,
    terra::project(depth_raster,"epsg:3857", method="average"),
    breaks = c(0, -40, -200, -2000, -Inf),
    biodiv_df)
  )
  
  expect_true(terra::crs(split_features[[1]]) == terra::crs(biodiv_raster))
  
  split_features <- suppressWarnings(prior3D::split_rast(
    terra::extend(biodiv_raster, c(1,1,0,0)),
    terra::extend(depth_raster, c(0,0,1,1)),
    breaks = c(0, -40, -200, -2000, -Inf),
    biodiv_df
  ))
  
  expect_true(terra::ext(split_features[[1]]) == terra::ext(biodiv_raster))
  
  split_features <- prior3D::split_rast(
                      biodiv_raster,
                      depth_raster,
                      breaks = c(0, -40, -200, -2000, -Inf),
                      biodiv_df,
                      val_depth_range = FALSE
                    )
  
  # aaptos_aaptos should have features only at [0,-40] level
  expect_false(all(is.na(terra::values(split_features[[1]]$aaptos_aaptos))))
  for (i in 2:4)
    expect_true(all(is.na(terra::values(split_features[[i]]$aaptos_aaptos))))
  
  # Test errors
  split_features <- try(prior3D::split_rast(
                          biodiv_raster = 1,
                          depth_raster,
                          breaks = c(0, -40, -200, -2000, -Inf),
                          biodiv_df
                    ),
                    silent = TRUE)
  
  expect_true(class(split_features) == "try-error")
  
  split_features <- try(prior3D::split_rast(
                          biodiv_raster,
                          depth_raster = 1,
                          breaks = c(0, -40, -200, -2000, -Inf),
                          biodiv_df
                  ),
                  silent = TRUE)
  
  # Suppose that we have only benthic with existence only at [0,-40] level
  split_features <- prior3D::split_rast(
                      biodiv_raster,
                      depth_raster,
                      breaks = c(0, -40, -200, -2000, -Inf),
                      data.frame(
                        species_name=biodiv_df$species_name,
                        pelagic = 0,
                        min_z = -30,
                        max_z = 0
                      ),
                      val_depth_range = FALSE
                    )
  
  expect_equal(suppressWarnings(sapply(1:4,
                 function(i, split_features)
                   all(is.na(values(split_features[[i]]))),
                 split_features=split_features)),
               c(FALSE, TRUE, TRUE, TRUE)
               )
  
  split_features <- try(prior3D::split_rast(
                          biodiv_raster,
                          depth_raster,
                          breaks = c(0, -40, -200, -2000, -Inf),
                          biodiv_df = 1
                  ),
                  silent = TRUE)
  
  expect_true(class(split_features) == "try-error")
  
  
  split_features <- try(prior3D::split_rast(
                          biodiv_raster,
                          depth_raster,
                          breaks = c(0, -40, -200, -2000, -Inf),
                          biodiv_df[,c("min_z", "max_z")]
                  ),
                  silent = TRUE)
  
  expect_true(class(split_features) == "try-error")
  
  split_features <- try(prior3D::split_rast(
                          biodiv_raster,
                          depth_raster,
                          breaks = c(0, -40, -200, -2000, -Inf),
                          biodiv_df[,c("pelagic", "min_z", "max_z")]
                  ),
                  silent = TRUE)
  
  expect_true(class(split_features) == "try-error")
  
  split_features <- try(prior3D::split_rast(
                          biodiv_raster,
                          depth_raster,
                          breaks = c(0, -40, -200, -2000, -Inf),
                          biodiv_df[,c("species_name", "min_z", "max_z")]
                  ),
                  silent = TRUE)
  
  expect_true(class(split_features) == "try-error")
  
  
    split_features <- try(prior3D::split_rast(
                                biodiv_raster,
                                depth_raster,
                             breaks = c(0, -40, -200, -2000, -Inf),
                             data.frame(
                               biodiv_df[,c("species_name", "pelagic")],
                               "min_z" = 100,
                               "max_z" = 0
                             ),
                             val_depth_range = FALSE
                     ),
                  silent = TRUE)
  
  expect_true(class(split_features) == "try-error")
  
  # No errors produced if only pelagic species
  split_features <- try(prior3D::split_rast(
                          biodiv_raster,
                          depth_raster,
                          breaks = c(0, -40, -200, -2000, -Inf),
                          data.frame(
                            biodiv_df[,-2],
                            pelagic = 1
                          ),
                          val_depth_range = TRUE
                        ))
  
  expect_false(class(split_features) == "try-error")
})