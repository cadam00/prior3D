\name{Compare_2D_3D}
\alias{Compare_2D_3D}
\title{Compare 2D vs 3D prioritization algorithms}
\usage{
Compare_2D_3D(biodiv_raster, depth_raster, breaks, biodiv_df, val_depth_range = TRUE,
priority_weights = NULL, budget_percents = seq(0,1,0.1), budget_weights = "equal",
penalty = 0, edge_factor = 0.5, gap = 0.1, threads = 1L, sep_priority_weights = ",",
portfolio = "gap", portfolio_opts = list(number_solutions = 10, pool_gap = 0.1),
sep_biodiv_df = ",", locked_in_raster = NULL, locked_out_raster = NULL, verbose = FALSE)
}

\description{
Compare 2D vs 3D prioritization algorithms
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
\item{priority_weights}{
\code{data.frame} object or file path (CSV, TXT, XLS, or XLSX) containing group names of
biodiversity features in the first column and corresponding group weights in the second
column. This data.frame attributes distinct prioritization weights to different
biodiversity features or groups of features.
}
\item{budget_percents}{
Numeric value \eqn{[0,1]} or vector containing budget percentages to use. The default is
\code{seq(0,1,0.1)}.
}
\item{budget_weights}{
Numeric weight vector for budget_percents allocation among depth levels.
Otherwise it can be a string with one of the choices \code{"equal"}, \code{"area"} or
\code{"richness"}.
Alternatively, it can be a numerical vector with custom weights corresponding to each
depth layer, where the first value corresponds to the surface and last one corresponds to
the bottom of the sea. The weights are normalized if their sum exceeds 1. If not
specified, an equal distribution of budget among depth levels is used, as the default.
}
\item{penalty}{
Numeric penalty applied to each depth zone, as defined in the\cr
\code{prioritizr::add_boundary_penalties}.
}
\item{edge_factor}{
Numeric edge factor applied to each depth zone, as defined in the\cr
\code{prioritizr::add_boundary_penalties}.
}
\item{gap}{
The optimality gap for the solver, as defined in the \pkg{prioritizr} package.
The default gap is 0.1.
}
\item{threads}{
The number of solver threads to be used. The default is 1.
}
\item{sep_priority_weights}{
Separator used in priority_weights file, if priority_weights is in path format.
}
\item{portfolio}{
The portfolio to be used, choosing between \code{"extra"},
\code{"gap"}, \code{"cuts"} and \code{"shuffle"} portfolios.
The default is \code{"gap"}. \code{portfolio=""} indicates that no portfolio is used.
For more about \code{portfolios} see \pkg{prioritizr}.
}
\item{portfolio_opts}{
The \pkg{prioritizr} portfolio options to be used.
}
\item{sep_biodiv_df}{
Separator\cr
used in biodiv_df file, if biodiv_df is in path format.
}
\item{locked_in_raster}{
An optional \code{locked_in_raster} SpatRaster to be used.
Note that these areas are considered as zero-cost.
}
\item{locked_out_raster}{
An optional \code{locked_out_raster} SpatRaster to be used.
Note that these areas are excluded from the solution.
}
\item{verbose}{
If \code{verbose = TRUE}, then solver messages are printed as well. The default is
\code{FALSE}.
}
}

\details{
To facilitate comparisons between 3D and 2D approaches, the \code{compare_2D_3D()}
function is provided in the package. This function enables users to conduct all steps of
the analysis (data generation, setting and solving the optimization problem
and producing outputs), by executing both 2D and 3D approaches, with similar settings,
that facilitate comparisons.  The function generates corresponding maps and graphs for
both approaches.

The \code{\link[prior3D]{split_rast}} function is used to convert 2D distributions of
biodiversity features (rasters) into a 3D format.

Here the \code{biodiv_df} can have the following
column names (independently of their order and any other names are ignored):
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
\item{\code{"group"}}: \bold{Optional} column with the group weights names.
}

Except from \code{biodiv_df}, an additional \code{data.frame} object can also be used for
defining group weights, named \code{priority_weights}. If used, this \code{data.frame}
object must have two columns:
\itemize{
\item{\code{"group"}}: \bold{Mandatory} column with the group weights names.
\item{\code{"weight"}}: \bold{Mandatory} column with the group weights.
}

In case that no feature weights are desired, then \code{priority_weights} can be kept
to \code{NULL}.

\code{breaks} must be in correspondence to depth_raster file.
For example, if depth_raster has range \eqn{[10, -3000]}, then a breaks vector of
\code{c(0,-40,-200,-2000,-Inf)} will create depth levels
\eqn{[0,-40],\\(-40,-200], (-200, -2000], (-2000, -\infty)}
and set to NA cells with values greater than \eqn{0}.

If \code{val_depth_range = TRUE} (default), then no correction is done and the depth range
of the biodiversity features is derived from the corresponding feature distribution
raster and so \code{"min_z"} and \code{"max_z"} are ignored.
If \code{val_depth_range = FALSE}, then the function uses the minimum and maximum depth
information provided in the biodiv_df, so as to remove feature occurrences outside their
expected range.

\code{budget_percents}: Budget reflects the desired level of protection to be modeled.
It ranges from 0 to 1, with 0 indicating no resources available for protection,
while 1 signifies resources sufficient to protect the entire study area. Typically,
setting a budget of 0.3 corresponds to the 30\% conservation target (i.e. 30\% of the
total area set aside for conservation).
Users also have the flexibility to define multiple budget levels using a vector,
allowing for the exploration of various protection scenarios. For instance, a vector like
\code{c(0.1, 0.3, 0.5)} represents three scenarios where 10\%, 30\%, and 50\% of the
study area are designated for protection.

\code{budget_weights}: The \bold{Compare_2D_3D} function allows users to specify how the
budget is distributed among depth levels. Three allocation methods are available:
\enumerate{
\item{Equal Distribution: Allocates an equal share of the budget to each depth level\cr
(\code{budget_weights = "equal"}).}
\item{Proportional to Area: Allocates budget based on the spatial extent of each depth
level\cr
(\code{budget_weights = "area"}).}
\item{Proportional to Species Richness: Prioritizes budget allocation to depth levels with
higher species diversity (number of species). (\code{budget_weights = "richness"})}
}
Otherwise, it can be a numeric vector with length equal to the number of depth levels,
where each number indicates the budget share per depth level.

The solver used for solving the prioritization problems is the best available on the
computer, following the solver hierarchy of \pkg{prioritizr}.
}

\value{
A list containing the following objects (non-referenced are identical to the input ones):
\itemize{
\item{split_features}: output of \code{\link[prior3D]{split_rast}}
\item{solution3D}: list with 3D solution per budget percentage
\item{absolute_held3D}: \code{absolute_held} for 3D solutions (see
\code{\link[prior3D]{evaluate_3D}})
\item{overall_held3D}: \code{overall_held} for 3D solutions (see
\code{\link[prior3D]{evaluate_3D}})
\item{relative_helds3D}: \code{relative_held} for 3D solutions (see
\code{\link[prior3D]{evaluate_3D}})
\item{mean_overall_helds3D}: \code{base::mean} of \code{overall_held}
for 3D solution (see \code{\link[prior3D]{evaluate_3D}}) per budget
\item{sd_overall_helds3D}: \code{stats::sd} of \code{overall_held}
for 3D solution (see \code{\link[prior3D]{evaluate_3D}}) per budget
\item{depth_overall_available3D}: \code{depth_overall_available} for 3D solutions (see
\code{\link[prior3D]{evaluate_3D}})
\item{solution2D}: list with 2D solution per budget percentage
\item{absolute_held2D}: \code{absolute_held} for 3D solutions (see
\code{\link[prior3D]{evaluate_3D}})
\item{overall_held2D}: \code{overall_held} for 2D solutions (see
\code{\link[prior3D]{evaluate_3D}})
\item{relative_helds2D}: \code{relative_held} for 2D solutions (see
\code{\link[prior3D]{evaluate_3D}})
\item{mean_overall_helds2D}: \code{base::mean} of \code{overall_held}
for 2D solution (see \code{\link[prior3D]{evaluate_3D}}) per budget
\item{sd_overall_helds2D}: \code{stats::sd} of \code{overall_held}
for 2D solution (see \code{\link[prior3D]{evaluate_3D}}) per budget
\item{depth_overall_available2D}: \code{depth_overall_available} for 2D solutions (see
\code{\link[prior3D]{evaluate_3D}})
\item{names_features}: names of features used
\item{total_amount}: \code{total_amount} of features used
(see \code{\link[prior3D]{evaluate_3D}})
\item{overall_total_amount}: \code{overal_total_amount} of names of features used
(see \code{\link[prior3D]{evaluate_3D}})
\item{jaccard_coef}: \code{\link[prior3D]{terra_jaccard}} per pair of 2D and 3D solutions,
given each budget
\item{depth_levels_names}: Depth levels names
\item{biodiv_raster}: \code{biodiv_raster} used, after cleaning
\item{biodiv_df}: \code{biodiv_df} used after cleaning
}
}

\references{
Hanson, Jeffrey O, Richard Schuster, Nina Morrell, Matthew Strimas-Mackey, Brandon
P M Edwards, Matthew E Watts, Peter Arcese, Joseph Bennett, and Hugh P Possingham. 2024.
prioritizr: Systematic Conservation Prioritization in R. \href{https://prioritizr.net}{
https://prioritizr.net}.

Lehtomäki, Joona (2016). Comparing prioritization methods, 21 June.\cr
Available at: \href{https://rpubs.com/jlehtoma/priocomp}{
https://rpubs.com/jlehtoma/priocomp}
(Accessed 1 June 2024).
}

\seealso{
\code{
\link[prior3D]{evaluate_3D},
\link[prior3D]{terra_jaccard},
\link[prior3D]{split_rast},
\link[prior3D]{plot_Compare_2D_3D}
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

out_2D_3D <- Compare_2D_3D(biodiv_raster = biodiv_raster,
                           depth_raster = depth_raster,
                           breaks = c(0, -40, -200, -2000, -Inf),
                           biodiv_df = biodiv_df,
                           budget_percents = seq(0, 1, 0.1),
                           budget_weights = "richness",
                           threads = parallel::detectCores(),
                           portfolio = "gap",
                           portfolio_opts = list(number_solutions = 10))

plot_Compare_2D_3D(out_2D_3D, to_plot = "all", add_lines=TRUE)

# Arbitrary random weights
priority_weights <- data.frame(c("A", "B", "C"), c(0.001, 1000, 1))
names(priority_weights) <- c("group", "weight")
biodiv_df$group <- rep(c("A", "B", "C"), length.out=20)
out_2D_3D <- Compare_2D_3D(biodiv_raster = biodiv_raster,
                           depth_raster = depth_raster,
                           breaks = c(0, -40, -200, -2000, -Inf),
                           biodiv_df = biodiv_df,
                           priority_weights = priority_weights,
                           budget_percents = seq(0, 1, 0.1),
                           budget_weights = "richness",
                           threads = parallel::detectCores(),
                           portfolio = "gap",
                           portfolio_opts = list(number_solutions = 10))

plot_Compare_2D_3D(out_2D_3D, to_plot = "all", add_lines=TRUE)
}
}
