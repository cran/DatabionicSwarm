\name{PswarmRadiusParallel}
\alias{PswarmRadiusParallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Intern function, do not use yourself}
\description{
Finds the weak Nash equilibirium of the data bots for one epoch depending on a
radius, which requires the setting of constants, grid, and so on in, see
\code{\link{Pswarm}}.
}
\usage{
PswarmRadiusParallel(DataBotsPos, DataDists, AllallowedDBPosR0, IndPossibleDBPosRe,
IndPossibleDBPosIm, Lines, Columns, Radius, NumAllDB, NumChoDB, NumFreeShape1,
NumJumps, Origin1, Origin2, Happiness, MinIterations, HappinessInclination, Eps, debug)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{DataBotsPos}{Numeric vector [1:NumJumps*n*2] containing the current
positions and all positions for considered/possible jumps which can be
computed (depending on number of jumps parameter NumJumps) for the databots
on two dimensions.}
\item{DataDists}{Numeric vector with vectorized distance matrix of the datapoints in the original (high-dimensional) data space}
\item{AllallowedDBPosR0}{NumericMatrix, see \code{AllallowedDBPosR0} in \code{\link{setPolarGrid}}}
\item{IndPossibleDBPosRe}{Numeric Vector of possible positions of the 1st coordinate.}
\item{IndPossibleDBPosIm}{Numeric Vector of possible positions of the 2nd coordinate.}
\item{Lines}{Integer stating the number of Lines the polar grid consists of.}
\item{Columns}{Integer stating the number of columns the polar grid consists of.}
\item{Radius}{Numeric (Integer) stating the moving radius of the databots}
\item{NumAllDB}{Integer total number of databots}
\item{NumChoDB}{Integer number of databots chosen for moving/jumps.}
\item{NumFreeShape1}{Integer stating the first dimension of the numeric matrix book keeping the possible position grid}
\item{NumJumps}{Integer number of jumps}
\item{Origin1}{Numeric origin coordinate 1}
\item{Origin2}{Numeric origin coordinate 2}
\item{Happiness}{Numeric value indicating the global happiness over all databots}
\item{MinIterations}{asdf}
\item{HappinessInclination}{asdf}
\item{Eps}{optional, double: Stop criterion for convergence of each epoche.}
\item{debug}{optional, bool: If TRUE prints status every 100 iterations}
}

\value{
list of
\item{AllDataBotsPos}{ComplexVector, indizes of DataBot Positions after a weak Nash equlibrium is found}
\item{stressverlauf}{NumericVector, intern result, for debugging only}
\item{fokussiertlaufind}{NumericVector, intern result, for debugging only}
}
\details{
Algorithm is described in [Thrun, 2018, p. 95, Listing 8.1].
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through
Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer,
Heidelberg, ISBN: 978-3-658-20539-3, \doi{10.1007/978-3-658-20540-9}, 2018.

[Thrun/Ultsch, 2021]  Thrun, M. C., and Ultsch, A.: Swarm Intelligence for
Self-Organized Clustering, Artificial Intelligence, Vol. 290, pp. 103237,
\doi{10.1016/j.artint.2020.103237}, 2021.

[Stier/Thrun, 2024] Stier, Q. and Thrun, M. C.: An efficient multicore CPU
implementation of the DatabionicSwarm, 18th conference of the International
Federation of Classification Societies (IFCS), San José, Costa Rica, July 14-19,
2024.
}
\author{
Quirin Stier
}
\keyword{swarm}
\keyword{nash}
\keyword{polar}
