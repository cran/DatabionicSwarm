\name{DatabionicSwarm-package}
\alias{DatabionicSwarm-package}
\alias{Databonic swarm}
\alias{DBS}
\docType{package}
\title{
Databionic Swarm (DBS)
}
\description{
Algorithms implementing populations of agents which interact with one 
another and sense their environment may exhibit emergent behavior such 
as self-organization and swarm intelligence. Here a swarm system, called
 databionic swarm (DBS), is introduced which is able to adapt itself to 
 structures of high-dimensional data such as natural clusters characterized
 by distance and/or density based structures in the data space [Thrun, 2017].
 The first module is the parameter-free projection method Pswarm,
 which exploits the concepts of self-organization and emergence,
 game theory, swarm intelligence and symmetry considerations. 
 The second module is a parameter-free high-dimensional data 
 visualization technique, which generates projected points on 
 a topographic map with hypsometric colors [Thrun et al., 2016]
 based on the generalized U-matrix [Ultsch/Thrun, 2017]. 
 The third module is the clustering method itself with non-critical parameters.
 The clustering can be verified by the visualization and vice versa. 
 The term DBS refers to the method as a whole. 
 DBS enables even a non-professional in the field of data mining 
 to apply its algorithms for visualization and/or clustering to 
 data sets with completely different structures drawn from diverse research fields.
}
\details{
\tabular{ll}{
Package: \tab Databonic swarm\cr
Type: \tab Package\cr
Version: \tab 0.9.6\cr
Date: \tab 2017-05-18\cr
License: \tab CC BY-NC-SA 4.0\cr
}
DBS is a flexible and robust clustering framework that consists
 of three independent modules. The first module is the parameter-free
 projection method Pswarm, which exploits the concepts of self-organization
 and emergence, game theory, swarm intelligence and symmetry considerations. 
 The second module is a parameter-free high-dimensional data visualization technique, 
 which generates projected points on a topographic map with hypsometric colors, 
 called the generalized U-matrix. The third module is a clustering method with no 
 sensitive parameters. The clustering can be verified by the visualization and vice versa. 
 The term DBS refers to the method as a whole. 
 
For further details, see Databionic swarm in the Phd thesis, chapter 8.
}
\author{
Michal Thrun

Maintainer: Michael Thrun <m.thrun@gmx.net>
}
\references{
[Thrun, 2017]  Thrun, M. C.:A System for Projection Based Clustering through Self-Organization and Swarm Intelligence, (Doctoral dissertation), Philipps-Universität Marburg, Marburg, 2017.

[Ultsch/Thrun, 2017]  Ultsch, A., & Thrun, M. C.: Credible Visualizations for Planar Projections, in Cottrell, M. (Ed.), 12th International Workshop on Self-Organizing Maps and Learning Vector Quantization, Clustering and Data Visualization (WSOM), IEEE Xplore, France, 2017.

[Thrun et al., 2016]  Thrun, M. C., Lerch, F., Loetsch, J., & Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision (WSCG), Vol. 24, Plzen, http://wscg.zcu.cz/wscg2016/short/A43-full.pdf, 2016.

}
\examples{
data('Lsun3D')
##2d projection, without instant visualization of steps
#DistanceMatrix hast to be defined by the user.
InputDistances=as.matrix(dist(Lsun3D$Data))
\donttest{
projection=Pswarm(InputDistances)
#2d projection, with instant visualization 
#of steps and DataMatrix (Distance is Euclidean per default)

projection=Pswarm(Lsun3D$Data,Cls=Lsun3D$Cls,PlotIt=T)
#
##Computation of Generalized Umatrix
# If Non Euclidean Distances are used, Please Use \code{SammonsMapping}
# from the ProjectionBasedClustering package with the correct OutputDimension
# to generate a new DataMatrix from the distances (see SheppardDiagram
# or KruskalStress)
visualization=GeneratePswarmVisualization(Data = Lsun3D$Data,

projection$ProjectedPoints,projection$LC)
## Visualizuation of GenerelizedUmatrix, 
# Estimation of the Number of Clusters=Number of valleys
library(GeneralizedUmatrix)
plotTopographicMap(visualization$Umatrix,visualization$Bestmatches)
## Automatic Clustering
# number of Cluster from dendrogram (PlotIt=TRUE) or visualization 
Cls=DBSclustering(k=3, Lsun3D$Data, 

visualization$Bestmatches, visualization$LC,PlotIt=FALSE)
# Verification, often its better to mark Outliers manually
library(GeneralizedUmatrix)
plotTopographicMap(visualization$Umatrix,visualization$Bestmatches,Cls)
}
\dontrun{
library(ProjectionBasedClustering)
Cls2=interactiveClustering(visualization$Umatrix, 
visualization$Bestmatches, Cls)
}
\dontshow{
data2=matrix(runif(n = 100),10,10)
distance=as.matrix(dist(data2))
res=Pswarm(distance,LC = c(10,12))
Cls=DBSclustering(k=2, data2, 

res$ProjectedPoints, res$LC,PlotIt=FALSE)

resUmatrix=GeneratePswarmVisualization(data2,res$ProjectedPoints,res$LC)
}

}
\keyword{DBS}
\keyword{DataBionic}
\keyword{swarm}
\keyword{swarm intelligence}
\keyword{swarm intelligence}
\keyword{Databionic}
\keyword{cluster analysis}
\keyword{clustering}
\keyword{visualization}
\keyword{DR}
\keyword{Dimensionality Reduction}
\keyword{projection}
\keyword{Pswarm}
\keyword{self-organization}
\keyword{emergence}
\keyword{nash}
\keyword{equilibrium}
\keyword{game theory}
\keyword{projection}
\keyword{projection method}