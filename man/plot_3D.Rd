\name{plot_3D}

\alias{plot_3D}

\title{Plot output of \link[prior3D]{prioritize_3D}}

\description{
Plot summarized output of \link[prior3D]{prioritize_3D}
}

\usage{
plot_3D(x, to_plot = "all", add_lines = TRUE)
}

\arguments{
\item{x}{
Output of \link[prior3D]{prioritize_3D}.
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
This function plots the summarized output of \link[prior3D]{prioritize_3D}
for all the selected budgets.
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
Becker, R.A., Wilks, A.R., Brownrigg, R., & Minka, T.P. (2023). maps: Draw
Geographical Maps. R package version 3.4.2,
\href{https://CRAN.R-project.org/package=maps}{
https://CRAN.R-project.org/package=maps}
}

\seealso{
\code{\link[prior3D]{prioritize_3D}}
}

\examples{
\dontrun{
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
}
}