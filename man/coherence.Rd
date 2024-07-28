\name{coherence}
\alias{coherence}
\title{Coherence metrics}
\description{
Coherence metrics
}
\usage{
coherence(x, w, metric = "sa", normalize = TRUE, plot = TRUE, addlines = TRUE, ...)
}
\arguments{
\item{x}{
Output of \bold{Compare_2D_3D} function.
}
\item{w}{
Integer indicating square window dimensions.
}
\item{metric}{
Cohension metric to use. It must be one of \code{"sa"}, \code{"sku"} and \code{"rao"}
}
\item{normalize}{
If \code{TRUE}, then sum of solutions is normalized at a \eqn{[0,1]} scale.
}
\item{plot}{
If \code{TRUE}, then coherence maps are ploted.
}
\item{addlines}{
If \code{TRUE}, then border lines from \code{maps::map} are ploted.
}
\item{...}{
Further arguments passed in function, based on metric choice. See more in Details.
}
}
\details{
\code{"sa"} and \code{"sku"} are derived from \code{geodiv::focal_metrics} and
in ellipsis (\code{...}) further arguments of
\code{geodiv::focal_metrics}  are passed.

\code{metric = "rao"} is derived from \code{rasterdiv::paRao} and in ellipsis
(\code{...}) further arguments of \code{rasterdiv::paRao} are passed.
}

\references{
Rocchini, Duccio, Matteo Marcantonio, Daniele Da Re, Giovanni Bacaro, Enrico Feoli, Giles
Foody, Reinhard Furrer, et al. 2021. "From zero to infinity: Minimum to maximum diversity
of the planet by spatio-parametric Rao’s quadratic entropy." Global Ecology and
Biogeography 30 (5): 2315.
\doi{10.1111/geb.13270}.
 
Rocchini, Duccio, Elisa Thouverai, Matteo Marcantonio, Martina Iannacito, Daniele Da Re,
Michele Torresani, Giovanni Bacaro, et al. 2021. "rasterdiv - An Information Theory
tailored R package for measuring ecosystem heterogeneity from space: To the origin and
back." Methods in Ecology and Evolution 12 (6): 2195.
\doi{10.1111/2041-210X.13583}.

Smith, Annie C., Phoebe Zarnetske, Kyla Dahlin, Adam Wilson, and Andrew Latimer. 2023.
Geodiv: Methods for Calculating Gradient Surface Metrics.
\href{https://CRAN.R-project.org/package=geodiv}{
https://CRAN.R-project.org/package=geodiv}.

Becker OScbRA, Minka ARWRvbRBEbTP, Deckmyn. A (2023). maps: Draw
Geographical Maps. R package version 3.4.2,
\href{https://CRAN.R-project.org/package=maps}{
https://CRAN.R-project.org/package=maps}
}

\seealso{
\code{\link[prior3D]{Compare_2D_3D}}
}

\examples{
\dontrun{
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
}
}

