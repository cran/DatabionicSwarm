\name{Delta3DWeightsC}
\alias{Delta3DWeightsC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
intern function, do not use yourself
}
\description{
Delta3DWeightsC
}
\usage{
Delta3DWeightsC(vx, Datasample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{vx}{Array [1:n,1:m,1:d] of neuron weights on a nxm grid with d dimensional
weights.}
\item{Datasample}{One observation of a d-dimensional datapoint.}
}
\value{
\item{vx}{Array [1:n,1:m,1:l]}
}
\details{
Algorithm is described in [Thrun, 2018, p. 95, Listing 8.1].
}
\references{
[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \doi{10.1007/978-3-658-20540-9}, 2018. 
}
\author{
Michael Thrun
}
\keyword{swarm}
\keyword{nash}
\keyword{polar}
