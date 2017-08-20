\name{ProjectedPoints2Grid}
\alias{ProjectedPoints2Grid}

\title{
Transforms ProjectedPoints to a grid
}
\description{
quantized xy cartesianncoordinates of ProjectedPoints
}
\usage{
ProjectedPoints2Grid(ProjectedPoints, Lines, Columns,PlotIt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ProjectedPoints}{[1:n,1:2] matrix of cartesian xy coordinates}
  \item{Lines}{double, length of small side of the rectangular grid}
  \item{Columns}{double, length of big side of the rectangular grid}
  \item{PlotIt}{optional, bool, show the result?}
}
\details{
intern function
}
\value{
BestMatches[1:n,1:3] columns in order: Key,Lines,Columns
}

\author{
Michael Thrunb

}
\seealso{\code{\link{GeneratePswarmVisualization}}
}

\keyword{GeneratePswarmVisualization}% use one of  RShowDoc("KEYWORDS")
\keyword{generlized Umatrix}% __ONLY ONE__ keyword per line
