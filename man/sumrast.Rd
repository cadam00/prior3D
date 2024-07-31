\name{sumrast}

\alias{sumrast}

\title{Sum list of SpatRaster objects.}

\description{
Sum list of SpatRaster objects.
}

\usage{
sumrast(x, normalize = TRUE)
}

\arguments{
\item{x}{
List of SpatRaster objects.
}
\item{normalize}{
If \code{TRUE}, then sum of solutions is normalized at a \eqn{[0,1]} scale.
}
}

\value{
A SpatRaster object.
}

\seealso{
\code{\link[prior3D]{plot_sumrast}}
}

\examples{
set.seed(42)
x <- terra::rast(matrix(rbinom(100, 1, 0.2), nrow=10))
y <- terra::rast(matrix(rbinom(100, 1, 0.8), nrow=10))
sumrast(list(x, y))
}