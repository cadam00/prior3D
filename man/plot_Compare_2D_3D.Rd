\name{plot_Compare_2D_3D}

\alias{plot_Compare_2D_3D}

\title{Plot output of \link[prior3D]{Compare_2D_3D}}

\description{
Plot summarized output of \link[prior3D]{Compare_2D_3D}
}

\usage{
plot_Compare_2D_3D(x, to_plot = "all", add_lines = TRUE)
}

\arguments{
\item{x}{
Output of \link[prior3D]{Compare_2D_3D}.
}
\item{to_plot}{
Any of \code{"maps"}, \code{"relative_held"} or \code{"all"}. The default is \code{"all"}.
See more in Details.
}
\item{add_lines}{
If \code{TRUE}, then border lines from \strong{maps::map} are ploted as well.
}
}

\details{
This function plots the summarized output of \link[prior3D]{Compare_2D_3D}
for all selected budgets.
The produced plot can contain information about:
\itemize{
\item{\code{"maps"}}: produced maps normalized at a \eqn{[0,1]} scale.
\item{\code{"relative_held"}}: percentage of protection for all features per depth level.
\item{\code{"all"}}: both \code{"maps"} and \code{"relative_held"}.
}
}

\value{
A plot.
}

\references{
Becker, R. A., Wilks, A. R., Brownrigg, R., & Minka, T. P. (2023). maps: Draw
Geographical Maps. R package version 3.4.2,
\href{https://CRAN.R-project.org/package=maps}{
https://CRAN.R-project.org/package=maps}
}

\seealso{
\code{\link[prior3D]{Compare_2D_3D}}
}

\examples{
\donttest{
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
}
}
