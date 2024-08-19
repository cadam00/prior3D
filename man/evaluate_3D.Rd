\name{evaluate_3D}

\alias{evaluate_3D}

\title{Evaluate prioritization solution over 3D feature distributions}

\description{
Evaluate prioritization solution over 3D feature distributions.
}

\usage{
evaluate_3D(solution, split_features)
}

\arguments{
\item{solution}{
prioritization solution SpatRaster object.
}
\item{split_features}{
A list of SpatRaster objects representing depth zones, where each element
corresponds to a different depth level, ranging from surface to the bottom of the sea.
The function uses the output of the split_rast function, but other multilevel (3D)
distribution data that conform to this structure is also acceptable.
First list element corresponds to the shallowest distribution and the last list element to
the deepest one.
}
}


\details{
This function evaluates a prioritization solution over 3D feature distributions, by
estimating the relative and overall percentages of features under protection, as
designated by the prioritization solution.
}

\value{
A list containing the following evaluation metrics:
\itemize{
\item{relative_held_raw}: relative percentage under protection per feature and
per depth level
\item{relative_held}: mean percentage under protection of all features per depth level
\item{overall_held}: overall percentage under protection per feature
\item{overall_available}: relative percentage under protection per feature and per
depth level, over total amount of each feature
\item{depth_overall_available}: averages of overall_available per depth
\item{absolute_held}: absolute value per feature and per depth under protection
\item{total_amount}: absolute value per feature and depth level
}
}

\examples{
\dontrun{
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
                             biodiv_df)

out_3D <- prioritize_3D(split_features = split_features,
                        depth_raster = depth_raster,
                        breaks = c(0, -40, -200, -2000, -Inf),
                        budget_percents = 0.3,
                        budget_weights = "richness",
                        threads = parallel::detectCores(),
                        portfolio = "gap",
                        portfolio_opts = list(number_solutions = 10))

evaluate_3D(out_3D$solution3D$budget0.3, split_features)
}
}
