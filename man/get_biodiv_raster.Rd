\name{get_biodiv_raster}

\alias{get_biodiv_raster}

\title{Example biodiversity raster}

\description{
Example biodiversity raster
}

\usage{
get_biodiv_raster()
}

\details{
Example of input \code{biodiv_raster} used for functions.
}

\value{
SpatRaster object with distribution of features.
}

\references{
Kaschner, K., Kesner-Reyes, K., Garilao, C., Segschneider, J., Rius-Barile,
J., Rees, T., & Froese, R. (2019). AquaMaps: Predicted range maps
for aquatic species. \href{https://www.aquamaps.org}{https://www.aquamaps.org}
}

\examples{
biodiv_raster <- get_biodiv_raster()
terra::plot(biodiv_raster[[1:4]])
}