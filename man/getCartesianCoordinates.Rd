\name{getCartesianCoordinates}
\alias{getCartesianCoordinates}
\title{Intern function: Transformation of Databot indizes to coordinates}
\description{
Transforms Databot indizes to exact cartesian coordinates on an toroid two dimensional grid.
}
\usage{
getCartesianCoordinates(DataBotsPosRe, DataBotsPosIm, GridRadius, GridAngle,
QuadOrHexa = TRUE)
}
\arguments{
\item{DataBotsPosRe}{[1:N] real part of complex vector Two Indizes per Databot
describing its positions in an two dimensional grid}
\item{DataBotsPosIm}{[1:N] imaginary part of complex vector Two Indizes per
Databot describing its positions in an two dimensional grid}
\item{GridRadius}{[Columns, Lines] Radii Matrix of all possible Positions of
DataBots in Grid, see also documentation of \code{\link{setPolarGrid}}}
\item{GridAngle}{[Columns, Lines] Angle Matrix of all possible Positions of
DataBots in Grid, see also documentation of \code{\link{setPolarGrid}}}
\item{QuadOrHexa}{Optional, FALSE=If DataPos on hexadiagonal grid, round to 2
decimals after value, Default=TRUE}
}
\value{
\item{BestMatchingUnits}{[1:N,2] coordinates on an two dimensional grid for
each databot excluding unique key, such that by using
\code{\link{GeneratePswarmVisualization}} a visualization of the Pswarm
projection is possible}
}
\author{Michael Thrun}
\details{Transformation is described in [Thrun, 2018, p. 93].}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \doi{10.1007/978-3-658-20540-9}, 2018. 
}
