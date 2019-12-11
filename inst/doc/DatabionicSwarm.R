## ----setup, include=FALSE-----------------------------------------------------
library(rgl)
#library(rglwidget) #does not show any pictures in Rmarkdown if rglwidget()  is called after rgl
setupKnitr()
knitr::opts_chunk$set(echo = TRUE,
                      fig.align = "center",
                      message=FALSE,
                      warning = FALSE,
                      webgl = TRUE,
                      dpi=50,
                      fig.width = 7, 
                      fig.height = 7,
                      fig.keep = "all"
                      )

## ----results = "hide"---------------------------------------------------------
library(DatabionicSwarm)
data('Hepta')
InputDistances = as.matrix(dist(Hepta$Data))
projection = Pswarm(InputDistances)

## ----results = "hide",webGL = TRUE,fig.keep="none"----------------------------
library(DatabionicSwarm)
library(GeneralizedUmatrix)
visualization = GeneratePswarmVisualization(
  Data = Hepta$Data, 
  projection$ProjectedPoints, 
  projection$LC)

GeneralizedUmatrix::plotTopographicMap(
  visualization$Umatrix, 
  visualization$Bestmatches, 
  NoLevels = 10)

## ----webGL = TRUE,fig.keep="none"---------------------------------------------
library(DatabionicSwarm)
library(GeneralizedUmatrix)
Cls = DBSclustering(
  k = 7,
  Hepta$Data,
  visualization$Bestmatches,
  visualization$LC,
  PlotIt = FALSE
)
GeneralizedUmatrix::plotTopographicMap(
  visualization$Umatrix,
  visualization$Bestmatches,
  Cls,
  NoLevels = 10)
