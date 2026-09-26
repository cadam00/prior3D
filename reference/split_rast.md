# Split 2D feature distributions into 3D ones

Split 2D feature distributions into 3D ones

## Usage

``` r
split_rast(biodiv_raster, depth_raster, breaks, biodiv_df, val_depth_range=TRUE,
sep_biodiv_df=",")
```

## Arguments

- biodiv_raster:

  SpatRaster object or folder path with 2D feature distributions as
  layers.

- depth_raster:

  SpatRaster object or file path with elevation/bathymetric map.

- breaks:

  Numeric vector defining the range of depth layers to use.

- biodiv_df:

  `data.frame` or a file path (CSV, TXT, XLS, or XLSX) containing
  additional information about biodiversity features.

- val_depth_range:

  No correction of the splitted 3D distributions based on depth range of
  the biodiversity features (`"min_z"` and `"max_z"` from `biodiv_df`)
  is needed.

- sep_biodiv_df:

  The separator used in biodiv_df file, if biodiv_df is in path format.

## Details

This function is used to convert 2D distributions of biodiversity
features (rasters) into a 3D format. Here the `biodiv_df` can have the
following column names (independently of their order and any other names
are ignored):

- `"species_name"`: **Mandatory** column with the feature names, which
  must be the same with biodiv_raster.

- `"pelagic"`: **Mandatory** column about the features' behaviour.
  `TRUE` means that this feature is pelagic and `FALSE` means that this
  feature is benthic.

- `"min_z"`: **Optional** column about the minimum vertical range of
  features. `NA` values are translated as unlimited upward feature
  movement.

- `"max_z"`: **Optional** column about the maximum vertical range of
  features. `NA` values are translated as unlimited downward feature
  movement.

`breaks` must be in correspondence to depth_raster file. For example, if
depth_raster has range \\\[10, -3000\]\\, then a breaks vector of
`c(0,-40,-200,-2000,-Inf)` will create depth levels \\\[0,-40\],\\
(-40,200\], (-200, -2000\], (-2000, -\infty)\\ and set to NA cells with
values greater than \\0\\.

If `val_depth_range = TRUE` (default), then no correction is done and
the depth range of the biodiversity features is derived from the
corresponding feature distribution raster and so `"min_z"` and `"max_z"`
are ignored. If `val_depth_range = FALSE`, then the function uses the
minimum and maximum depth information provided in the biodiv_df, so as
to remove feature occurrences outside their expected range.

## Value

A list containing species distributions for each bathymetric layer, that
are necessary for further 3D analysis. List names are indicating the
depth levels.

## Examples

``` r
biodiv_raster <- get_biodiv_raster()
depth_raster <- get_depth_raster()
data(biodiv_df)

# You can split features' 2D distributions into 3D ones and then run only 3D analysis
split_features <- split_rast(biodiv_raster,
                             depth_raster,
                             breaks = c(0, -40, -200, -2000, -Inf),
                             biodiv_df)
```
