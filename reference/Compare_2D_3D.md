# Compare 2D vs 3D prioritization algorithms

Compare 2D vs 3D prioritization algorithms

## Usage

``` r
Compare_2D_3D(biodiv_raster, depth_raster, breaks, biodiv_df, val_depth_range = TRUE,
priority_weights = NULL, budget_percents = seq(0,1,0.1), budget_weights = "equal",
penalty = 0, edge_factor = 0.5, gap = 0.1, threads = 1L, sep_priority_weights = ",",
portfolio = "gap", portfolio_opts = list(number_solutions = 10, pool_gap = 0.1),
sep_biodiv_df = ",", locked_in_raster = NULL, locked_out_raster = NULL, verbose = FALSE)
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

- priority_weights:

  `data.frame` object or file path (CSV, TXT, XLS, or XLSX) containing
  group names of biodiversity features in the first column and
  corresponding group weights in the second column. This data.frame
  attributes distinct prioritization weights to different biodiversity
  features or groups of features.

- budget_percents:

  Numeric value \\\[0,1\]\\ or vector containing budget percentages to
  use. The default is `seq(0,1,0.1)`.

- budget_weights:

  Numeric weight vector for budget_percents allocation among depth
  levels. Otherwise it can be a string with one of the choices
  `"equal"`, `"area"` or `"richness"`. Alternatively, it can be a
  numerical vector with custom weights corresponding to each depth
  layer, where the first value corresponds to the surface and last one
  corresponds to the bottom of the sea. The weights are normalized if
  their sum exceeds 1. If not specified, an equal distribution of budget
  among depth levels is used, as the default.

- penalty:

  Numeric penalty applied to each depth zone, as defined in the  
  [`prioritizr::add_boundary_penalties`](https://prioritizr.net/reference/add_boundary_penalties.html).

- edge_factor:

  Numeric edge factor applied to each depth zone, as defined in the  
  [`prioritizr::add_boundary_penalties`](https://prioritizr.net/reference/add_boundary_penalties.html).

- gap:

  The optimality gap for the solver, as defined in the prioritizr
  package. The default gap is 0.1.

- threads:

  The number of solver threads to be used. The default is 1.

- sep_priority_weights:

  Separator used in priority_weights file, if priority_weights is in
  path format.

- portfolio:

  The portfolio to be used, choosing between `"extra"`, `"gap"`,
  `"cuts"` and `"shuffle"` portfolios. The default is `"gap"`.
  `portfolio=""` indicates that no portfolio is used. For more about
  `portfolios` see prioritizr.

- portfolio_opts:

  The prioritizr portfolio options to be used.

- sep_biodiv_df:

  Separator  
  used in biodiv_df file, if biodiv_df is in path format.

- locked_in_raster:

  An optional `locked_in_raster` SpatRaster to be used. Note that these
  areas are considered as zero-cost.

- locked_out_raster:

  An optional `locked_out_raster` SpatRaster to be used. Note that these
  areas are excluded from the solution.

- verbose:

  If `verbose = TRUE`, then solver messages are printed as well. The
  default is `FALSE`.

## Details

To facilitate comparisons between 3D and 2D approaches, the
`compare_2D_3D()` function is provided in the package. This function
enables users to conduct all steps of the analysis (data generation,
setting and solving the optimization problem and producing outputs), by
executing both 2D and 3D approaches, with similar settings, that
facilitate comparisons. The function generates corresponding maps and
graphs for both approaches.

The
[`split_rast`](https://cadam00.github.io/prior3D/reference/split_rast.md)
function is used to convert 2D distributions of biodiversity features
(rasters) into a 3D format.

Here the `biodiv_df` can have the following column names (independently
of their order and any other names are ignored):

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

- `"group"`: **Optional** column with the group weights names.

Except from `biodiv_df`, an additional `data.frame` object can also be
used for defining group weights, named `priority_weights`. If used, this
`data.frame` object must have two columns:

- `"group"`: **Mandatory** column with the group weights names.

- `"weight"`: **Mandatory** column with the group weights.

In case that no feature weights are desired, then `priority_weights` can
be kept to `NULL`.

`breaks` must be in correspondence to depth_raster file. For example, if
depth_raster has range \\\[10, -3000\]\\, then a breaks vector of
`c(0,-40,-200,-2000,-Inf)` will create depth levels
\\\[0,-40\],\\(-40,-200\], (-200, -2000\], (-2000, -\infty)\\ and set to
NA cells with values greater than \\0\\.

If `val_depth_range = TRUE` (default), then no correction is done and
the depth range of the biodiversity features is derived from the
corresponding feature distribution raster and so `"min_z"` and `"max_z"`
are ignored. If `val_depth_range = FALSE`, then the function uses the
minimum and maximum depth information provided in the biodiv_df, so as
to remove feature occurrences outside their expected range.

`budget_percents`: Budget reflects the desired level of protection to be
modeled. It ranges from 0 to 1, with 0 indicating no resources available
for protection, while 1 signifies resources sufficient to protect the
entire study area. Typically, setting a budget of 0.3 corresponds to the
30% conservation target (i.e. 30% of the total area set aside for
conservation). Users also have the flexibility to define multiple budget
levels using a vector, allowing for the exploration of various
protection scenarios. For instance, a vector like `c(0.1, 0.3, 0.5)`
represents three scenarios where 10%, 30%, and 50% of the study area are
designated for protection.

`budget_weights`: The **Compare_2D_3D** function allows users to specify
how the budget is distributed among depth levels. Three allocation
methods are available:

1.  Equal Distribution: Allocates an equal share of the budget to each
    depth level  
    (`budget_weights = "equal"`).

2.  Proportional to Area: Allocates budget based on the spatial extent
    of each depth level  
    (`budget_weights = "area"`).

3.  Proportional to Species Richness: Prioritizes budget allocation to
    depth levels with higher species diversity (number of species).
    (`budget_weights = "richness"`)

Otherwise, it can be a numeric vector with length equal to the number of
depth levels, where each number indicates the budget share per depth
level.

The solver used for solving the prioritization problems is the best
available on the computer, following the solver hierarchy of prioritizr.

## Value

A list containing the following objects (non-referenced are identical to
the input ones):

- split_features: output of
  [`split_rast`](https://cadam00.github.io/prior3D/reference/split_rast.md)

- solution3D: list with 3D solution per budget percentage

- absolute_held3D: `absolute_held` for 3D solutions (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- overall_available3D: `overall_available` for 3D solutions (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- overall_held3D: `overall_held` for 3D solutions (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- relative_helds3D: `relative_held` for 3D solutions (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- mean_overall_helds3D: [`base::mean`](https://rdrr.io/r/base/mean.html)
  of `overall_held` for 3D solution (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))
  per budget

- sd_overall_helds3D: [`stats::sd`](https://rdrr.io/r/stats/sd.html) of
  `overall_held` for 3D solution (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))
  per budget

- depth_overall_available3D: `depth_overall_available` for 3D solutions
  (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- solution2D: list with 2D solution per budget percentage

- absolute_held2D: `absolute_held` for 2D solutions (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- overall_available2D: `overall_available` for 2D solutions (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- overall_held2D: `overall_held` for 2D solutions (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- relative_helds2D: `relative_held` for 2D solutions (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- mean_overall_helds2D: [`base::mean`](https://rdrr.io/r/base/mean.html)
  of `overall_held` for 2D solution (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))
  per budget

- sd_overall_helds2D: [`stats::sd`](https://rdrr.io/r/stats/sd.html) of
  `overall_held` for 2D solution (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))
  per budget

- depth_overall_available2D: `depth_overall_available` for 2D solutions
  (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- names_features: names of features used

- total_amount: `total_amount` of features used (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- overall_total_amount: `overal_total_amount` of names of features used
  (see
  [`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md))

- jaccard_coef:
  [`terra_jaccard`](https://cadam00.github.io/prior3D/reference/terra_jaccard.md)
  per pair of 2D and 3D solutions, given each budget

- depth_levels_names: Depth levels names

- biodiv_raster: `biodiv_raster` used, after cleaning

- biodiv_df: `biodiv_df` used after cleaning

## References

Hanson, Jeffrey O, Richard Schuster, Nina Morrell, Matthew
Strimas-Mackey, Brandon P M Edwards, Matthew E Watts, Peter Arcese,
Joseph Bennett, and Hugh P Possingham. 2024. prioritizr: Systematic
Conservation Prioritization in R. <https://prioritizr.net>.

Lehtomäki, Joona (2016). Comparing prioritization methods, 21 June.  
Available at: <https://rpubs.com/jlehtoma/priocomp> (Accessed 1 June
2024).

## See also

` `[`evaluate_3D`](https://cadam00.github.io/prior3D/reference/evaluate_3D.md)`, `[`terra_jaccard`](https://cadam00.github.io/prior3D/reference/terra_jaccard.md)`, `[`split_rast`](https://cadam00.github.io/prior3D/reference/split_rast.md)`, `[`plot_Compare_2D_3D`](https://cadam00.github.io/prior3D/reference/plot_Compare_2D_3D.md)` `

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

plot_Compare_2D_3D(out_2D_3D, to_plot = "all", add_lines=TRUE)

# Arbitrary random weights
priority_weights <- data.frame(c("A", "B", "C"), c(0.001, 1000, 1))
names(priority_weights) <- c("group", "weight")
biodiv_df$group <- rep(c("A", "B", "C"), length.out=20)
out_2D_3D <- Compare_2D_3D(biodiv_raster = biodiv_raster,
                           depth_raster = depth_raster,
                           breaks = c(0, -40, -200, -2000, -Inf),
                           biodiv_df = biodiv_df,
                           priority_weights = priority_weights,
                           budget_percents = seq(0, 1, 0.1),
                           budget_weights = "richness",
                           threads = parallel::detectCores(),
                           portfolio = "gap",
                           portfolio_opts = list(number_solutions = 10))

plot_Compare_2D_3D(out_2D_3D, to_plot = "all", add_lines=TRUE)
} # }
```
