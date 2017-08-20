\name{getCartesianCoordinates}
\alias{getCartesianCoordinates}
\title{ Transform DataBot Indizes two exac cartesian coordinates on an toroid two dimensional grid. }
\arguments{
\item{DataBotsPos[1:AnzData]}{complex vector Two Indizes per Databot describing its positions in an two dimensional grid}
\item{GridRadius[Lines,Columns]}{Radii Matrix of all possible Positions of DataBots in Grid, see als Doku of setPolarGrid()}
\item{GridAngle[Lines,Columns]}{Angle Matrix of all possible Positions of DataBots in Grid, see als Doku of setPolarGrid()}
\item{Lines,Columns}{Size of planar toroid two dimensional grid}
\item{QuadOrHexa}{FALSE=If DataPos on hexadiagonal grid, runde auf 2 nachkommestellen, Default=TRUE}
}
\value{BestMatchingUnits[1:AnzData,2]         coordinates on an two dimensional grid for each databot excluding unique key, sucht that}
\description{
 Transform DataBot Indizes two exac cartesian coordinates on an toroid two dimensional grid. }
\author{
MT 01/2015}
