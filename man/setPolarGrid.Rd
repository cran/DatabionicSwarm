\name{setPolarGrid}
\alias{setPolarGrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sets the polar grid
}
\description{
sets a polar grid for a swarm in an rectangular shape
}
\usage{
setPolarGrid(Lines,Columns,QuadOrHexa,PlotIt,global)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lines}{Integer, hast to be able to be divided by 2}
  \item{Columns}{Integer, with Columns>=Lines}
  \item{QuadOrHexa}{bool, default(TRUE) If False Hexagonal grid, default quad grid}
  \item{PlotIt}{bool, default(FALSE)}
  \item{global}{bool,  default(TRUE), intern parameter, how shall the radii be calculated?}
}

\value{
list of 
\item{GridRadii}{matrix [1:Lines,1:Columns], Radii Matrix of all possible Positions of DataBots in Grid}
\item{GridAngle}{matrix [1:Lines,1:Columns], Angle Matrix of all possible Positions of DataBots in Grid}
\item{AllallowedDBPosR0}{matrix [1:Lines+1,1:Columns+1],  Matrix of radii in polar coordinates respecting origin (0,0) of all allowed DataBots Positions in one jump}
\item{AllallowedDBPosPhi0}{matrix [1:Lines+1,1:Columns+1], 
  # V$AllallowedDBPosPhi0[Lines+1,Lines+1]        Matrix of angle in polar coordinates respecting origin (0,0) of all allowed DataBots Positions in one jump}
}
\author{
Michael Thrun
}
\seealso{
\code{\link{pswarmCpp}}
}

\keyword{hexagonal}% use one of  RShowDoc("KEYWORDS")
\keyword{grid}
\keyword{swarm}