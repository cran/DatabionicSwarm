## ----setup, include=FALSE-----------------------------------------------------
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

## ----echo = FALSE-------------------------------------------------------------
if (!requireNamespace("rmarkdown") || !rmarkdown::pandoc_available("1.12.3")) {
  warning("This vignette requires pandoc version 1.12.3; code will not run in older versions.")
  knitr::opts_chunk$set(eval = FALSE)
}

## ----echo = FALSE-------------------------------------------------------------
if (!requireNamespace("rgl")) {
  warning("This vignette requires the package rgl; code will not run without this package.")
  knitr::opts_chunk$set(eval = FALSE)
}else{
  library(rgl)
  setupKnitr()
}

## ----results = "hide"---------------------------------------------------------
library(DatabionicSwarm)
data('Hepta')
InputDistances = as.matrix(dist(Hepta$Data))
projection = Pswarm(InputDistances)

## ----results = "hide",webGL = TRUE,fig.keep="none"----------------------------
library(DatabionicSwarm)
library(GeneralizedUmatrix)
genUmatrixList = GeneratePswarmVisualization(
  Data = Hepta$Data, 
  projection$ProjectedPoints, 
  projection$LC)

GeneralizedUmatrix::plotTopographicMap(
  genUmatrixList$Umatrix, 
  genUmatrixList$Bestmatches, 
  NoLevels = 10)

## ----webGL = TRUE,fig.keep="none"---------------------------------------------
library(DatabionicSwarm)
library(GeneralizedUmatrix)
Cls = DBSclustering(
  k = 7,
  Hepta$Data,
  genUmatrixList$Bestmatches,
  genUmatrixList$LC,
  PlotIt = FALSE
)
GeneralizedUmatrix::plotTopographicMap(
  genUmatrixList$Umatrix,
  genUmatrixList$Bestmatches,
  Cls,
  NoLevels = 10)

