RobustNormalization = function (Data) 
{
  if (is.vector(Data)) {
    quants = quantile(Data, c(0.01, 0.5, 0.99))
    minX = quants[1]
    maxX = quants[3]
    return((Data - minX)/(maxX - minX))
  }
  else if (is.matrix(Data)) {
    cols = ncol(Data)
    xtrans = Data
    for (i in 1:cols) {
      xtrans[, i] = RobustNormalization(as.vector(Data[, i]))
    }
    return(xtrans)
  }
  else {
    tryCatch({
      warning("Data is not a vector or a matrix. Trying as.matrix")
      return(RobustNormalization(as.matrix(Data)))
    }, error = function(e) {
      stop("It did not work")
    })
  }
}