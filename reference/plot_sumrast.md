# Plot sum list of SpatRaster objects.

Plot sum list of SpatRaster objects.

## Usage

``` r
plot_sumrast(x, normalize = TRUE, add_lines = TRUE, ...)
```

## Arguments

- x:

  List of SpatRaster objects.

- normalize:

  If `TRUE`, then sum of solutions is normalized at a \\\[0,1\]\\ scale.

- add_lines:

  If `TRUE`, then border lines from **maps::map** are ploted as well.

- ...:

  Further arguments passed to **terra::plot**

## Value

A plot.

## See also

[`sumrast`](https://cadam00.github.io/prior3D/reference/sumrast.md)

## Examples

``` r
set.seed(42)
x <- terra::rast(matrix(rbinom(100, 1, 0.2), nrow=10))
y <- terra::rast(matrix(rbinom(100, 1, 0.8), nrow=10))
plot_sumrast(list(x, y), add_lines = FALSE)
```
