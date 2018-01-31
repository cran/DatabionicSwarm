\name{RobustNormalization}
\alias{RobustNormalization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RobustNormalization
}
\description{
RobustNormalization as described in [Milligan/Cooper, 1988].
}
\usage{
RobustNormalization(Data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{[1:n,1:d]

}
}
\details{
Range of most of values for each feature is not between zero and one. Fpr a more precise description please read [Thrun, 2018, p.17].
}
\value{
TransformedData[1:n,1:d]
}
\references{
[Milligan/Cooper, 1988]  Milligan, G. W., & Cooper, M. C.: A study of standardization of variables in cluster analysis, Journal of Classification, Vol. 5(2), pp. 181-204. 1988.

[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 

}
\author{
Michael Thrun
}


\keyword{RobustNormalization}% use one of  RShowDoc("KEYWORDS")
