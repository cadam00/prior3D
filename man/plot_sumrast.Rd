\name{plot_sumrast}

\alias{plot_sumrast}

\title{Plot sum list of SpatRaster objects.}

\description{
Plot sum list of SpatRaster objects.
}

\usage{
plot_sumrast(x, normalize = TRUE, add_lines = TRUE, ...)
}

\arguments{
\item{x}{
List of SpatRaster objects.
}
\item{normalize}{
If \code{TRUE}, then sum of solutions is normalized at a \eqn{[0,1]} scale.
}
\item{add_lines}{
If \code{TRUE}, then border lines from \strong{maps::map} are ploted as well.
}
\item{...}{
Further arguments passed to \strong{terra::plot}
}
}

\value{
A plot.
}

\seealso{
\code{\link[prior3D]{sumrast}}
}


\examples{
set.seed(42)
x <- terra::rast(matrix(rbinom(100, 1, 0.2), nrow=10))
y <- terra::rast(matrix(rbinom(100, 1, 0.8), nrow=10))
plot_sumrast(list(x, y), add_lines = FALSE)
}
