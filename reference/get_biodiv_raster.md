# Example biodiversity raster

Example biodiversity raster

## Usage

``` r
get_biodiv_raster()
```

## Details

Example of input `biodiv_raster` used for functions.

## Value

SpatRaster object with distribution of features.

## References

Kaschner, K., Kesner-Reyes, K., Garilao, C., Segschneider, J.,
Rius-Barile, J., Rees, T., & Froese, R. (2019). AquaMaps: Predicted
range maps for aquatic species. <https://www.aquamaps.org>

## Examples

``` r
biodiv_raster <- get_biodiv_raster()
terra::plot(biodiv_raster[[1:4]])
```
