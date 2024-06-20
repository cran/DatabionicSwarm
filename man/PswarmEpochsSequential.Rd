\name{PswarmEpochsSequential}
\alias{PswarmEpochsSequential}
\title{Intern function, do not use yourself}
\description{
Finds the weak Nash equilibirium of the data bots for one epoch depending on a
radius, which requires the setting of constants, grid, and so on in, see
\code{\link{Pswarm}}.
}
\usage{
PswarmEpochsSequential(AllDataBotsPos, MyDistanceMatrix, IndPossibleDBPosR,
AllFreePosR0, NumAllDB, Lines, Columns, Origin, Happiness, GridRadii, GridAngle,
QuadOrHexa, RadiusVector, Rmin, Rmax, Cls, Debug, pp, PlotIt = FALSE,
Verbose = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{AllDataBotsPos}{Complex vector [1:n] of the current positions for the databots on a 2d and real plane in complex numbers.}
\item{MyDistanceMatrix}{Numeric vector with vectorized distance matrix of the datapoints in the original (high-dimensional) data space}
\item{IndPossibleDBPosR}{Numeric vector containing the possible positions around a databot dependent on the radius.}
\item{AllFreePosR0}{NumericMatrix, see \code{AllallowedDBPosR0} in \code{\link{setPolarGrid}}.}
\item{NumAllDB}{Integer total number of databots}
\item{Lines}{Integer stating the number of Lines the polar grid consists of.}
\item{Columns}{Integer stating the number of columns the polar grid consists of.}
\item{Origin}{Numeric origin of the positions of grid in two dimensions}
\item{Happiness}{Numeric value indicating the global happiness over all databots}
\item{GridRadii}{Numeric matrix with radius information of polar transformation for each grid position}
\item{GridAngle}{Numeric matrix with angle information of polar transformation for each grid position}
\item{QuadOrHexa}{optional, bool: If TRUE prints status every 100 iterations}
\item{RadiusVector}{Numeric vector stating all moving radius in a descending order (cooling down scheme).}
\item{Rmin}{Integer stating minimum radius.}
\item{Rmax}{Integer stating maximum radius.}
\item{Cls}{Integer vector stating the classification vector for each datapoints/databots.}
\item{Debug}{optional, bool: If TRUE prints information for debugging.}
\item{pp}{Numeric vector stating ratio of number of jumping simultaneously DataBots of one eppoch (per nash-equilibirum), this vector is linearly monotonically decreasing.}
\item{PlotIt}{optional, bool: If TRUE creates plot of projection after each epoch.}
\item{Verbose}{optional, integer stating degree of textual feedback. 0 = no output, 1 = basic notifications, 2 = progress bar, 3 = details.}
}

\value{
list of
\item{AllDataBotsPosRe}{Numeric vector [1:n] of the current positions for the databots on first of two dimensions.}
\item{AllDataBotsPosIm}{Numeric vector [1:n] of the current positions for the databots on second of two dimensions.}
\item{CourseOfHappiness}{NumericVector, states the global happiness value per epoch.}
\item{RadiusPerEpoch}{NumericVector, stating the radius used per epoch in order of computation.}
}
\details{
Algorithm is described in [Thrun, 2018, p. 95, Listing 8.1].
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \doi{10.1007/978-3-658-20540-9}, 2018. 
}
\author{
Quirin Stier
}
\keyword{swarm}
\keyword{nash}
\keyword{polar}
