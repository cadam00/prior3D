# Jaccard similarity coefficient among two SpatRaster objects

Jaccard coefficient among two SpatRaster objects

## Usage

``` r
terra_jaccard(x, y)
```

## Arguments

- x:

  SpatRaster object with binary values.

- y:

  SpatRaster object with binary values.

## Details

Jaccard similarity coefficient evaluates the percentage number equal to
the intersection between two sets, divided by the size of the union of
these sets.

## Value

A numeric value \\\[0,1\]\\.

## Examples

``` r
set.seed(42)
x <- terra::rast(matrix(rbinom(100, 1, 0.2), nrow=10))
y <- terra::rast(matrix(rbinom(100, 1, 0.8), nrow=10))
terra_jaccard(x, y)
#> [1] 0.2470588
```
