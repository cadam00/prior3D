\name{get_rast}

\alias{get_rast}

\title{Read multiple rast files}

\description{
Read multiple rast files contained in a folder path. Raster files must have either
.asc or .tif extension.
}

\usage{
get_rast(path)
}

\arguments{
\item{path}{
Path string of folder containing rast files.
}
}

\value{
A SpatRaster object.
}

\examples{
\dontrun{
feature_folder <- system.file("get_rast_example", package="prior3D")
get_rast(feature_folder)
}
}