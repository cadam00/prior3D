\name{terra_jaccard}

\alias{terra_jaccard}

\title{Jaccard similarity coefficient among two SpatRaster objects}

\description{
Jaccard coefficient among two SpatRaster objects
}

\usage{
terra_jaccard(x, y)
}

\arguments{
\item{x}{
SpatRast object with binary values.
}
\item{y}{
SpatRast object with binary values.
}
}

\details{
Jaccard similarity coefficient evaluates the percentage number equal to the intersection
between two sets, divided by the size of the union of these sets.
}

\value{
A numeric value \eqn{[0,1]}.
}

\examples{
\dontrun{
set.seed(42)
x <- terra::rast(matrix(rbinom(100, 1, 0.2), nrow=10))
y <- terra::rast(matrix(rbinom(100, 1, 0.8), nrow=10))
terra_jaccard(x, y)
}
}