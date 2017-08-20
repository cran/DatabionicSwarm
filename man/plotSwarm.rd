\name{plotSwarm}
\alias{plotSwarm}

\title{
intern function for plotting during pswarm process
}
\description{
intern function, plots the progess of the Pswarm algorithm after every nash equlibirum
}
\usage{
plotSwarm(Points,Cls,xlab,ylab,main)
}

\arguments{
  \item{Points}{ProjectedPoints or DataBot positions in cartesian coordinates}
   \item{Cls}{optional, Classification as a numeric vector, if given}
    \item{xlab}{='X', optional, string}
     \item{ylab}{='Y', optional, string}
     \item{main}{="DataBots", optional, string}
}

\author{
Michael Thrun
}

\seealso{
\code{\link{pswarmCpp}} with \code{PlotIt}=TRUE
}

\keyword{pswarmCpp}
\keyword{swarm}
