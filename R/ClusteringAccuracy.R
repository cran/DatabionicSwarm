ClusteringAccuracy=function(PriorCls,CurrentCls,K=9){
  if(length(unique(PriorCls))>9){
    warning('Too many clusters in PriorCls for RAM of single PC. Please use Cloud Computing, e.g. SparkR')
  }
  #author: 04/2018 MT
  #Note: symmetric ClsToTrueCls() which always works
  
  NormalizeCls <- function(Cls) {

    uniqueClasses <- sort(na.last = T, unique(Cls))
    numberOfClasses <- length(uniqueClasses)
    unique2Cls <- NULL #  initializing the vector
    
    for (i in 1:length(Cls)) {
      # calculating the indexes of elements of Cls in uniqueClasses
      unique2Cls <- c(unique2Cls, which(uniqueClasses == Cls[i]))
    }
    
    if (numberOfClasses > 0) {
      normalizedClasses <- c(1:numberOfClasses)
      normalizedCls <- normalizedClasses[unique2Cls]
    }
    else {
      normalizedClasses <- Cls
    }
    
    return(normalizedCls)
  }
  ####################################################
  ClassCount <- function(Cls) {
    # calulates statistics for   points in each group of the data
    # C <-ClassCount(Cls);
    # UniqueClasses <-C$uniqueClasses
    # CountPerClass <-C$countPerClass
    # NrOfClasses   <-C$numberOfClasses
    # ClassPercentages <-C$classPercentages
    #
    # INPUT
    # Cls(d)                          Cls(i) == ClusterNumber of data(i,:)
    #
    # OUTPUT list with:
    # UniqueClasses[1:NrOfClasses]      the  NrOfClasses unique classes in Cls
    # CountPerClass(NrOfClasses,n)    CountPerClass(i) is the Count of the data points in UniqueClasses(i)
    # NumberOfClasses                     the number of classes
    # ClassPercentages                the percentages of the classes
    # All2UniqInd                     UniqueClasses == Cls[All2UniqInd]
    # Uniq2AllInd                     Cls == UniqueClasses[Uniq2AllInd]
    
    uniqueClasses <- sort(na.last = T, unique(Cls))
    numberOfClasses <- length(uniqueClasses)
    countPerClass <-
      rep(0, numberOfClasses) # just initializing the vector, all values are replaced later.
    
    for (i in 1:numberOfClasses) {
      inClassI <-
        sum(Cls == uniqueClasses[i]) # counts all occurances of uniqueClass[i] in Cls
      countPerClass[i] = inClassI # updates countPerClass[i] to the number of occurances of uniqueClasses[i] in Cls.
    }
    
    classPercentages <- rep(0, numberOfClasses)
    
    for (i in 1:numberOfClasses) {
      classPercentages[i] <- (countPerClass[i] / sum(countPerClass)) * 100
    }
    
    return(
      list(
        UniqueClasses = uniqueClasses,
        CountPerClass = countPerClass,
        NumberOfClasses = numberOfClasses,
        ClassPercentages = classPercentages,
        All2UniqInd = vapply(uniqueClasses,function(i){max(which(Cls==i))},0),
        Uniq2AllInd = vapply(uniqueClasses, function(i){which(uniqueClasses == i)},0)
      )
    )
  }
  
  ReduceClassesToK=function(Cls,K=9){
    m=length(unique(Cls))
    ClsRes1=Cls
    if(m>K){
      p=length(seq(from=m,to=K,by=-1))-1
      for(i in 1:p){
        ClsTmp=ClsRes1
        V=ClassCount(ClsTmp)
        ind=order(V$CountPerClass)
        u=head(V$UniqueClasses[ind],2)
        ClsTmp[c(which(ClsTmp==u[1]),which(ClsTmp==u[2]))]=u[1]
        ClsRes1=ClsTmp
      }
    }
    return(ClsRes1)
  }
  ######################################################################
  
  
  standardCls <- NormalizeCls(PriorCls)
  givenCls <- NormalizeCls(CurrentCls)
  if(length(unique(givenCls))>K){
    warning('Too many clusters in CurrentCls for RAM of single PC. Combining clusters of smalest size.
            Alternativly, please use Cloud Computing, e.g. SparkR')
    givenCls=ReduceClassesToK(givenCls,K=K)
  }
  givenCls<- NormalizeCls(givenCls)
  
  uniqueClasses <- sort(na.last = T, unique(c(standardCls,givenCls)))
  requireNamespace('pracma')
  allPossiblePermutations <- pracma::perms(uniqueClasses)
  nrOfPermutations <- nrow(allPossiblePermutations)
  nrOfStdClasses <- ncol(allPossiblePermutations)
  givenClasses <- sort(na.last = T, unique(givenCls))
  nrOfGivenClasses <- length(givenClasses)
  renamedCls <- givenCls
  bestAccuracy <- 0
  #For every permutation
  for (i in 1:nrOfPermutations) {
    #set ground truth
    tryRenameCls <- givenCls
    
    #set a permutation of cls to be inspected
    newClassNames <- c(1:nrOfGivenClasses)
    newClassNames[1:nrOfStdClasses] <- allPossiblePermutations[i,]
    
    for (j in 1:nrOfGivenClasses) { #for every point
      #search
      tryRenameCls[which(givenCls == givenClasses[j])] <- newClassNames[j]
    }
    #true positives
    accuracy <- sum(tryRenameCls == standardCls)
    
    if (accuracy > bestAccuracy) {
      renamedCls <- tryRenameCls
      bestAccuracy <- accuracy
    }
  }
  
  bestAccuracy <- bestAccuracy / length(standardCls)
  
  return(Accuracy = bestAccuracy)
}