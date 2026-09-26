# Sum list of SpatRaster objects.

Sum list of SpatRaster objects.

## Usage

``` r
sumrast(x, normalize = TRUE)
```

## Arguments

- x:

  List of SpatRaster objects.

- normalize:

  If `TRUE`, then sum of solutions is normalized at a \\\[0,1\]\\ scale.

## Value

A SpatRaster object.

## See also

[`plot_sumrast`](https://cadam00.github.io/prior3D/reference/plot_sumrast.md)

## Examples

``` r
set.seed(42)
x <- terra::rast(matrix(rbinom(100, 1, 0.2), nrow=10))
y <- terra::rast(matrix(rbinom(100, 1, 0.8), nrow=10))
sumrast(list(x, y))
#> class       : SpatRaster
#> size        : 10, 10, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 10, 0, 10  (xmin, xmax, ymin, ymax)
#> coord. ref. : 
#> source(s)   : memory
#> name        : sum
#> min value   :   0
#> max value   :   1
```
