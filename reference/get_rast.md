# Read multiple rast files

Read multiple rast files contained in a folder path. Raster files must
have either .asc or .tif extension.

## Usage

``` r
get_rast(path)
```

## Arguments

- path:

  Path string of folder containing rast files.

## Value

A SpatRaster object.

## Examples

``` r
feature_folder <- system.file("get_rast_example", package="prior3D")
get_rast(feature_folder)
#> class       : SpatRaster
#> size        : 31, 83, 2  (nrow, ncol, nlyr)
#> resolution  : 0.5, 0.5  (x, y)
#> extent      : -5.5, 36, 30.5, 46  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> sources     : aaptos_aaptos.tif
#>               abietinaria_abietina.tif
#> names       : aaptos_aaptos, abietinaria_abietina
#> min values  :          0.01,                 0.01
#> max values  :             1,                 0.63
```
