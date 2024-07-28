\name{get_depth_raster}

\alias{get_depth_raster}

\title{Example depth raster}

\description{
Example depth raster
}

\usage{
get_depth_raster()
}

\details{
Example of input \code{depth_df} object for functions.
}

\value{
SpatRaster object with depth levels for Mediterranean.
}

\references{
GEBCO Compilation Group. (2021). GEBCO 2021 Grid.
\doi{10.5285/c6612cbe-50b3-0cff-e053-6c86abc09f8f}.
}

\examples{
\dontrun{
depth_raster <- get_depth_raster()
terra::plot(depth_raster)
}
}