\name{rDistanceToroidCsingle}
\alias{rDistanceToroidCsingle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
intern function for \code{\link{pswarmCpp}}
}
\description{
toroid distance calculation
}
\usage{
rDistanceToroidCsingle( AllDataBotsPosX, 
 AllDataBotsPosY, AllallowedDBPosR0,
 Lines, Columns,  Nullpunkt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AllDataBotsPosX}{NumericVector [1:n,1], positions of on grid}
  \item{AllDataBotsPosY}{NumericVector [1:n,1], positions of on grid}
  \item{AllallowedDBPosR0}{NumericMatrix}
  \item{Lines}{double}
  \item{Columns}{double}
  \item{Nullpunkt}{NumericVector}
}
\details{do not use yourself}
\value{
numeric matrix of toroid Distances[1:n,1:n]
}

\author{
Michael Thrun
}

\seealso{
\code{\link{pswarmCpp}}
}

\keyword{toroid}% use one of  RShowDoc("KEYWORDS")
\keyword{distances}% __ONLY ONE__ keyword per line
