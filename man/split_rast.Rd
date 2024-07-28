\name{split_rast}

\alias{split_rast}

\title{Split 2D feature distributions into 3D ones}

\description{
Split 2D feature distributions into 3D ones
}

\usage{
split_rast(biodiv_raster, depth_raster, breaks, biodiv_df, val_depth_range=TRUE,
sep_biodiv_df=",")
}

\arguments{
\item{biodiv_raster}{
SpatRaster object or folder path with 2D feature distributions as layers.
}
\item{depth_raster}{
SpatRaster object or file path with elevation/bathymetric map.
}
\item{breaks}{
Numeric vector defining the range of depth layers to use.
}
\item{biodiv_df}{
\code{data.frame} or a file path (CSV, TXT, XLS, or XLSX) containing additional
information about biodiversity features.
}
\item{val_depth_range}{
No correction of the splitted 3D distributions based on depth range of the biodiversity
features (\code{"min_z"} and \code{"max_z"} from \code{biodiv_df}) is needed.
}
\item{sep_biodiv_df}{
The separator used in biodiv_df file, if biodiv_df is in path format.
}
}



\details{
This function is used to convert 2D distributions of biodiversity
features (rasters) into a 3D format. Here the \code{biodiv_df} can have the
following column names (independently of their order and any other names are ignored):
\itemize{
\item{\code{"species_name"}}: \bold{Mandatory} column with the feature names, which must
be the same with biodiv_raster.
\item{\code{"pelagic"}}: \bold{Mandatory} column about the features' behaviour.
\code{TRUE} means that this feature is pelagic and \code{FALSE} means that this feature is
benthic.
\item{\code{"min_z"}}: \bold{Optional} column about the minimum vertical range of
features. \code{NA} values are translated as unlimited upward feature movement.
\item{\code{"max_z"}}: \bold{Optional} column about the maximum vertical range of
features. \code{NA} values are translated as unlimited downward feature movement.
}

\code{breaks} must be in correspondance to depth_raster file.
For example, if depth_raster has range \eqn{[10, -3000]}, then a breaks vector of
\code{c(0,-40,-200,-2000,-Inf)} will create depth levels \eqn{[0,-40],\\
(-40,200], (-200, -2000], (-2000, -\infty)}
and set to NA cells with values greater than \eqn{0}.

If \code{val_depth_range = TRUE} (default), then no correction is done and the depth range
of the biodiversity features is derived from the corresponding feature distribution
raster and so \code{"min_z"} and \code{"max_z"} are ignored.
If \code{val_depth_range = FALSE}, then the function uses the minimum and maximum depth
information provided in the biodiv_df, so as to remove feature occurances outside their
expected range.

}

\value{
A list containing species distributions for each bathymetric layer, that are necessary for
further 3D analysis. List names are indicating the depth levels.
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
                             biodiv_df)
}
}