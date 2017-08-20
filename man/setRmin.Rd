\name{setRmin}
\alias{setRmin}
\title{ estimates the minimal radius}
\arguments{
\item{Lines}{x-value determining the size of the map, i.e. how many open places for DataBots will be available  on the 2-dimensional grid BEWARE: has to be able to be divided by 2}
\item{Columns}{y-value  determining the size of the map, i.e. how many open places for DataBots will be available on the 2-dimensional grid Columns>Lines}
\item{AllallowedDBPosR0}{[1:Lines+1,1:Lines+1]Matrix of radii in polar coordinates respecting origin (0,0) of all allowed DataBots Positions in one jumpp                                   Prozentangabe als Zahl kleiner 1}
\item{p}{percent of gitterpositions, which should be considered}
}
\value{Rmin                                     Minimum Radius}
\description{
estimates the minimal radius in the annealing process of Pswarm}
\author{Michael Thrun
}
