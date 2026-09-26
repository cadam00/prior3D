# Plot output of [Compare_2D_3D](https://cadam00.github.io/prior3D/reference/Compare_2D_3D.md)

Plot summarized output of
[Compare_2D_3D](https://cadam00.github.io/prior3D/reference/Compare_2D_3D.md)

## Usage

``` r
plot_Compare_2D_3D(x, to_plot = "all", add_lines = TRUE)
```

## Arguments

- x:

  Output of
  [Compare_2D_3D](https://cadam00.github.io/prior3D/reference/Compare_2D_3D.md).

- to_plot:

  Any of `"maps"`, `"relative_held"` or `"all"`. The default is `"all"`.
  See more in Details.

- add_lines:

  If `TRUE`, then border lines from **maps::map** are ploted as well.

## Details

This function plots the summarized output of
[Compare_2D_3D](https://cadam00.github.io/prior3D/reference/Compare_2D_3D.md)
for all selected budgets. The produced plot can contain information
about:

- `"maps"`: produced maps normalized at a \\\[0,1\]\\ scale.

- `"relative_held"`: percentage of protection for all features per depth
  level.

- `"all"`: both `"maps"` and `"relative_held"`.

## Value

A plot.

## References

Becker, R. A., Wilks, A. R., Brownrigg, R., & Minka, T. P. (2023). maps:
Draw Geographical Maps. R package version 3.4.2,
<https://CRAN.R-project.org/package=maps>

## See also

[`Compare_2D_3D`](https://cadam00.github.io/prior3D/reference/Compare_2D_3D.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## This example requires commercial solver from 'gurobi' package for
## portfolio = "gap". Else replace it with e.g. portfolio = "shuffle" for using
## a free solver like the one from 'highs' package.

biodiv_raster <- get_biodiv_raster()
depth_raster <- get_depth_raster()
data(biodiv_df)

out_2D_3D <- Compare_2D_3D(biodiv_raster = biodiv_raster,
                           depth_raster = depth_raster,
                           breaks = c(0, -40, -200, -2000, -Inf),
                           biodiv_df = biodiv_df,
                           budget_percents = seq(0, 1, 0.1),
                           budget_weights = "richness",
                           threads = parallel::detectCores(),
                           portfolio = "gap",
                           portfolio_opts = list(number_solutions = 10))

plot_Compare_2D_3D(out_2D_3D, to_plot="all", add_lines=FALSE)
plot_Compare_2D_3D(out_2D_3D, to_plot="all", add_lines=TRUE)
plot_Compare_2D_3D(out_2D_3D, to_plot="maps", add_lines=TRUE)
plot_Compare_2D_3D(out_2D_3D, to_plot="relative_held", add_lines=TRUE)
} # }
```
