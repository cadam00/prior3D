# Coherence metrics

Coherence metrics

## Usage

``` r
coherence(x, w, metric = "sa", normalize = TRUE, plot = TRUE, addlines = TRUE, ...)
```

## Arguments

- x:

  Output of **Compare_2D_3D** function.

- w:

  Integer indicating square window dimensions.

- metric:

  Cohension metric to use. It must be one of `"sa"`, `"sku"` and `"rao"`

- normalize:

  If `TRUE`, then sum of solutions is normalized at a \\\[0,1\]\\ scale.

- plot:

  If `TRUE`, then coherence maps are ploted.

- addlines:

  If `TRUE`, then border lines from
  [`maps::map`](https://rdrr.io/pkg/maps/man/map.html) are ploted.

- ...:

  Further arguments passed in function, based on metric choice. See more
  in Details.

## Details

`"sa"` and `"sku"` are derived from
[`geodiv::focal_metrics`](https://rdrr.io/pkg/geodiv/man/focal_metrics.html)
and in ellipsis (`...`) further arguments of
[`geodiv::focal_metrics`](https://rdrr.io/pkg/geodiv/man/focal_metrics.html)
are passed.

`metric = "rao"` is derived from
[`rasterdiv::paRao`](https://mattmar.github.io/rasterdiv/reference/paRao.html)
and in ellipsis (`...`) further arguments of
[`rasterdiv::paRao`](https://mattmar.github.io/rasterdiv/reference/paRao.html)
are passed.

## Value

`numeric` vector with 2D and 3D metrics. The result maps are optionally
ploted.

## References

Rocchini, Duccio, Matteo Marcantonio, Daniele Da Re, Giovanni Bacaro,
Enrico Feoli, Giles Foody, Reinhard Furrer, et al. 2021. "From zero to
infinity: Minimum to maximum diversity of the planet by
spatio-parametric Rao’s quadratic entropy." Global Ecology and
Biogeography 30 (5): 2315.
[doi:10.1111/geb.13270](https://doi.org/10.1111/geb.13270) .

Rocchini, Duccio, Elisa Thouverai, Matteo Marcantonio, Martina
Iannacito, Daniele Da Re, Michele Torresani, Giovanni Bacaro, et al.
2021. "rasterdiv - An Information Theory tailored R package for
measuring ecosystem heterogeneity from space: To the origin and back."
Methods in Ecology and Evolution 12 (6): 2195.
[doi:10.1111/2041-210X.13583](https://doi.org/10.1111/2041-210X.13583) .

Smith, Annie C., Phoebe Zarnetske, Kyla Dahlin, Adam Wilson, and Andrew
Latimer. 2023. Geodiv: Methods for Calculating Gradient Surface Metrics.
<https://CRAN.R-project.org/package=geodiv>.

Becker OScbRA, Minka ARWRvbRBEbTP, Deckmyn. A (2023). maps: Draw
Geographical Maps. R package version 3.4.2,
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

coherence(out_2D_3D, w = 3, metric = "sa")
coherence(out_2D_3D, w = 3, metric = "sku")
coherence(out_2D_3D, w = 3, metric = "rao")
} # }
```
