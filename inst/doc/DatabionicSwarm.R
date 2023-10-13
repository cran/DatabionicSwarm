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
library(GeneralizedUmatrix)
data('Hepta')
InputDistances = as.matrix(dist(Hepta$Data))
projection = Pswarm(InputDistances)

## ----results = "hide",webGL = TRUE,fig.keep="none"----------------------------
genUmatrixList = GeneratePswarmVisualization(
  Data = Hepta$Data, 
  projection$ProjectedPoints, 
  projection$LC,
  Parallel=FALSE)#CRAN guidelines do not allow =TRUE in vignette

GeneralizedUmatrix::plotTopographicMap(
  genUmatrixList$Umatrix, 
  genUmatrixList$Bestmatches, 
  NoLevels = 10)

## ----webGL = TRUE,fig.keep="none"---------------------------------------------
Cls = DBSclustering(k = 7,
                    DataOrDistance = Hepta$Data,
                    BestMatches = genUmatrixList$Bestmatches,
                    LC = genUmatrixList$LC,
                    PlotIt = FALSE)
GeneralizedUmatrix::plotTopographicMap(genUmatrixList$Umatrix,
                                       genUmatrixList$Bestmatches,
                                       Cls,
                                       NoLevels = 10)

