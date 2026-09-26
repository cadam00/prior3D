# Example depth raster

Example depth raster

## Usage

``` r
get_depth_raster()
```

## Details

Example of input `depth_df` object for functions.

## Value

SpatRaster object with depth levels for Mediterranean.

## References

GEBCO Compilation Group. (2021). GEBCO 2021 Grid.
[doi:10.5285/c6612cbe-50b3-0cff-e053-6c86abc09f8f](https://doi.org/10.5285/c6612cbe-50b3-0cff-e053-6c86abc09f8f)
.

## Examples

``` r
depth_raster <- get_depth_raster()
terra::plot(depth_raster)
```
