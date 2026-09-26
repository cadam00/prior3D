# Plot output of [prioritize_3D](https://cadam00.github.io/prior3D/reference/prioritize_3D.md)

Plot summarized output of
[prioritize_3D](https://cadam00.github.io/prior3D/reference/prioritize_3D.md)

## Usage

``` r
plot_3D(x, to_plot = "all", add_lines = TRUE)
```

## Arguments

- x:

  Output of
  [prioritize_3D](https://cadam00.github.io/prior3D/reference/prioritize_3D.md).

- to_plot:

  Any of `"maps"`, `"relative_held"` or `"all"`. The default is `"all"`.
  See more in Details.

- add_lines:

  If `TRUE`, then border lines from **maps::map** are ploted as well.

## Details

This function plots the summarized output of
[prioritize_3D](https://cadam00.github.io/prior3D/reference/prioritize_3D.md)
for all the selected budgets. The produced plot can contain information
about:

- `"maps"`: produced maps normalized at a \\\[0,1\]\\ scale.

- `"relative_held"`: percentage of protection for all features per depth
  level.

- `"all"`: both `"maps"` and `"relative_held"`.

## Value

A plot.

## References

Becker, R.A., Wilks, A.R., Brownrigg, R., & Minka, T.P. (2023). maps:
Draw Geographical Maps. R package version 3.4.2,
<https://CRAN.R-project.org/package=maps>

## See also

[`prioritize_3D`](https://cadam00.github.io/prior3D/reference/prioritize_3D.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## This example requires commercial solver from 'gurobi' package for
## portfolio = "gap". Else replace it with e.g. portfolio = "shuffle" for using
## a free solver like the one from 'highs' package.

biodiv_raster <- get_biodiv_raster()
depth_raster <- get_depth_raster()
data(biodiv_df)

# You can split features' 2D distributions into 3D ones and then run only 3D analysis
split_features <- split_rast(biodiv_raster,
                             depth_raster,
                             breaks = c(0, -40, -200, -2000, -Inf),
                             biodiv_df,
                             val_depth_range=TRUE)

out_3D <- prioritize_3D(split_features = split_features,
                        depth_raster = depth_raster,
                        breaks = c(0, -40, -200, -2000, -Inf),
                        biodiv_df = biodiv_df,
                        priority_weights = NULL,#priority_weights,
                        budget_percents = seq(0, 1, 0.1),
                        budget_weights = "equal",
                        penalty = 0,
                        edge_factor = 0.5,
                        gap = 0.1,
                        threads = parallel::detectCores(),
                        sep_priority_weights = ",",
                        portfolio = "gap",
                        portfolio_opts = list(number_solutions = 10),
                        sep_biodiv_df = ",",
                        locked_in_raster = NULL,
                        locked_out_raster = NULL)

plot_3D(out_3D, to_plot="all", add_lines=FALSE)
plot_3D(out_3D, to_plot="all", add_lines=TRUE)
plot_3D(out_3D, to_plot="maps", add_lines=TRUE)
plot_3D(out_3D, to_plot="relative_held", add_lines=TRUE)
} # }
```
