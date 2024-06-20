Pswarm = function(DataOrDistance, Cls = NULL, QuadOrHexa = "Hexa", NumJumps = 4,
                  LC = NULL, Parallel = FALSE, NCores = "max",
                  Verbose = 1, PlotIt = FALSE, Debug = FALSE,
                  DistanceMeasure = "euclidean", Eps = 0.001){
  # V = PSwarm(DataOrDistance)
  # 
  # DESCRIPTION
  # ???
  # 
  # INPUT
  # DataOrDistance[1:n, 1:d]    Numeric matrix nxd. Two cases here:
  #                                 d=n  => assuming distance matrix
  #                                 d!=n => assuming data matrix with n cases
  #                                         and d features implying the need
  #                                         to compute the distance matrix
  #                                         internally
  # 
  # OPTIONAL
  # Cls[1:n]                    Numeric vector with class labels for each
  #                             observation in DataOrDistance.
  # QuadOrHexa                  Character indicating the geometry of tiles the 2D
  #                             projection plane is built with.
  #                             "Hexa" for hexagonal grid
  #                             "Quad" for quadratix grid
  # LC                          Numeric vector with two entries: Number of lines
  #                             and columns of the two-dimensional projection
  #                             plane
  # Parallel                    Boolean: TRUE  = parallel execution,
  #                                      FALSE = single thread execution.
  # NCores                      Character or integer: choice of number of cores
  #                             of CPU (in case). Can be 'max' or a number. The 
  #                             max will always be 'all available cores - 1', to
  #                             avoid core overload.
  # Verbose                     Boolean TRUE: no messages. FALSE: messages.
  #                             (Default: Verbose=TRUE)
  # PlotIt                      Boolean TRUE: show plot. FALSE: do not show plot.
  #                             (Default: PlotIt=FALSE)
  # Debug                       Boolean TRUE: give debug messages. FALSE: do not
  #                             show debug messages. (Default: Debug=FALSE)
  # DistanceMeasure             Character giving the name of distance measure.
  #                             (Default: DistanceMeasure=Euclidean)
  # 
  # OUTPUT
  # ProjectedPoints[1:n, 1:2]    Numeric matrix with the n projected 
  #                              observations from the data on the 2D plane.
  # LC[1:2]                      Numeric vector with number of Lines and Columns
  #                              of the 2D projection plane (Lines, Columns)
  # Control                      List of three items:
  #                                  CourseOfHappiness
  #                                  RadiusPerEpoch
  #                                  HappinessLastEpoch
  # 
  # 1st Author: Michael C. Thrun
  # 2nd Author: Quirin Stier
  
  #----------------------------------------------------------------------------#
  # Check if input is correct, otherwise return
  #----------------------------------------------------------------------------#
  T_init_1 = Sys.time()
  # DataOrDistance is numerical matrix
  if(!is.matrix(DataOrDistance) | !is.numeric(DataOrDistance)){
    stop("DataOrDistance is not a numeric matrix.")
  }
  
  # QuadOrHexa is one of "Option1" or "Option2"
  #if(!is.logical(QuadOrHexa)){
  #  stop("QuadOrHexa is not a logical (TRUE/FALSE.")
  #}
  QuadOrHexa = "Hexa"
  QuadOrHexa = FALSE # == "Hexa"
  
  # LC is NULL or numerical vector of size 2 (if not, set LC = c(NULL,NULL) # Global variable: Lines and Columns of grid)
  if(!is.null(LC)){
    if(!is.vector(LC) | (length(LC) != 2)){
      stop("LC is not a numeric vector of length 2.")
    }
  }
  
  if(isFALSE(Parallel)){
    if(NumJumps != 4){
      print("Parameters NumJumps is fixed to 4 in the sequential version!")
    }
    NumJumps = 4
  }else{
    if(is.character(NCores)){
      NCores = defaultNumThreads() - 1
    }else if(is.integer(as.integer(NCores))){
      if((NCores < 1) & (NCores >= defaultNumThreads())){
        NCores = defaultNumThreads() - 1
      }
    }else{
      NCores = defaultNumThreads() - 1
    }
    RcppParallel::setThreadOptions(numThreads = NCores)
    if(Verbose == 1){
      print(paste0("Parallel version is executed on ", NCores, " threads."))
    }
  }
  
  if(!is.numeric(NumJumps)){
    NumJumps = 4
  }
  
  if((NumJumps < 1)){
    NumJumps = 4
  }
  
  #----------------------------------------------------------------------------#
  # Initialization phase
  #----------------------------------------------------------------------------#
  # Distance matrix of input/data space
  if(dim(DataOrDistance)[1] != dim(DataOrDistance)[2]){
    tmpRes         = parallelDist::parallelDist(DataOrDistance,                 # Compute distances in parallel
                                                method = DistanceMeasure)       # Choice of distance measure
    tmpDM          = as.matrix(tmpRes,                                          # The simple matrix call is able to transform it correctly nowadays!
                               nrow = dim(DataOrDistance)[1])
    MyDistanceMatrix = tmpDM
  }else{
    warning("Data matrix is square. Assuming a distance matrix!")
    MyDistanceMatrix = DataOrDistance
  }
  
  NumAllDB = NumCases = dim(MyDistanceMatrix)[1]      # Number of cases/observations/...
  
  #----------------------------------------------------------------------------#
  # Grid size: Columns (x-axis) x Lines (y-axis)
  if(is.null(LC)){                # Global variable: Lines and Columns of grid; LC=[X,Y]
    LC = setGridSize(InputDistances = MyDistanceMatrix, minp = 0.01, maxp = 0.99, Verbose = Verbose)
  }
  
  if(is.null(Cls)){
    Cls = rep(1, NumCases)                  # Classes for each case
  }
  
  Lines   = LC[1]             # Lines of grid
  Columns = LC[2]             # Columns of grid
  Rmax    = Lines/2           # Max radius for databots to move

  assertCorrectDataDrivenApproach(Lines, Columns, NumAllDB)
  
  # Y Construct Databionicswarm setting for given data
  # Binary matrix of all free DataBotsPositions(=0), Positions of DataBots will be defined by 1
  V1 = setPolarGrid(Lines = Lines, Columns = Columns, QuadOrHexa = QuadOrHexa, PlotIt = FALSE)

  # Polar grid
  GridRadii      = V1$GridRadii               # Radii Matrix of all possible Positions of DataBots in Grid
  GridAngle      = V1$GridAngle               # Angle Matrix of all possible Positions of DataBots in Grid
  AllFreePosR0   = V1$AllallowedDBPosR0       # Radius-Matrix in polar coordinates respecting origin (0,0) of all allowed DataBots Positions in one jump
  AllFreePosPhi0 = V1$AllallowedDBPosPhi0     # Angle-Matrix in polar coordinates respecting origin (0,0) of all allowed DataBots Positions in one jump
  ## Geht schief wenn LC zu klein
  # Linear array index of 2 dimensional field/grid/matrix of possible positions of the databot for initialization
  init = sample(x = Lines * Columns, replace = F, size = NumAllDB) # Get NumCases many random indices
  #InitialLinearIndex
  
  binary      = matrix(0, nrow = Lines, ncol = Columns)    # Get zero filled grid field
  dballind    = which(binary == 0, arr.ind = T)            # Get index of empty fields
  DataBotsPos = dballind[init, ]                           # 
  
  # NOTA BENE
  # Es gibt eine Schablone um die Polarkoordinaten relativ fuer databot zu bestimmen
  # Es gibt eine tatsaechliche Position auf Matrix (Lines x Columns)
  # Dynamische Berechnung von IST Zustand der Databots aus aus beiden Datenstrukturen!
  
  AllDataBotsPosRe = DataBotsPos[, 1]
  AllDataBotsPosIm = DataBotsPos[, 2]
  
  Origin = which(AllFreePosR0 == 0, arr.ind = T) # PointZero within Polar grid
  
  # Rmin: minimum radius depending on grid size (LC) and p=0.05 (assuming uniform distribution of databots)
  if(QuadOrHexa){ #Bei quadratischen Gitter ist minimaler Radius groesser 3 noetig
    Rmin = setRmin(AllFreePosR0, Lines, Columns, NumAllDB, p = 0.05)
  } else{ # !QuadOrHexa #Rmin=1
    Rmin = setRmin(AllFreePosR0, Lines, Columns, NumAllDB, p = 0.05)
  }
  m                 = (0.5 - 0.05) / (Rmax - Rmin)
  b                 = 0.5 - m * Rmax
  pp                = m * 1:Rmax + b                                     # Number of bots which are allowed to move simultaneously?
  jumpthreshold     = 0
  fokussiertlaufind = 1
  RadiusVector      = seq(from = Rmax, by = -1, to = Rmin)
  Stress            = Inf
  alpha             = Rmax * 0.01
  
  #----------------------------------------------------------------------------#
  # End initialization phase
  #----------------------------------------------------------------------------#
  # Possible positions around a databot depending on radius range
  if(isTRUE(Parallel)){
    JumpsPerRadius = lapply(RadiusVector, function(Radius, AllPos, alpha, Lines){
      findPossiblePositionsCsingle(AllPos, Radius, alpha, Lines)
    }, AllFreePosR0, alpha, Lines)
  }else{
    IndPossibleDBPosR = findPossiblePositionsCsingle(AllFreePosR0, Rmax, alpha, Lines)
  }
  
  # Distance Matrix of output/projection
  RadiusPerEpoch = c()
  if(isTRUE(Parallel)){
    OutputDistance = rcppPar_DistanceToroid(AllDataBotsPosRe, AllDataBotsPosIm,
                                            AllFreePosR0, Lines, Columns, as.vector(Origin))
  }else{
    OutputDistance = rDistanceToroidCsingle(AllDataBotsPosRe, AllDataBotsPosIm,
                                            AllFreePosR0, Lines, Columns, as.vector(Origin))
  }
  
  # Happiness measure
  # NbghFUN = 1 - OutputDistance ^ 2 / (pi * Rmax ^ 2)
  # all(1 - OutputDistance ^ 2 / (pi * Rmax ^ 2) == 1 - ((OutputDistance ^ 2) / (pi * Rmax ^ 2)))
  NbghFUN   = 1 - ((OutputDistance^2)/(pi * Rmax ^ 2))
  NbghFUN[NbghFUN < 0] = 0
  N         = sum(NbghFUN)
  Happiness = sum(NbghFUN * MyDistanceMatrix) / N
  
  T_init_2 = Sys.time()
  
  if(Verbose >= 1){
    t_diff_1 = T_init_2 - T_init_1
    print(paste0("Initialization phase successful (Time: ", round(t_diff_1, 2), " ", units(t_diff_1), ")."))
    print("Starting swarm game ...")
  }
  
  T_DBS_1 = Sys.time()
  if(isTRUE(Parallel)){
    tmpResPar = PswarmEpochsParallel(AllDataBotsPosRe, AllDataBotsPosIm, MyDistanceMatrix, AllFreePosR0,
                                     GridRadii, GridAngle, JumpsPerRadius,
                                     NumJumps, NumAllDB, Lines, Columns, Origin, Happiness,
                                     QuadOrHexa, RadiusVector, Rmin, Rmax,
                                     Cls, Debug, pp, PlotIt = PlotIt, Verbose = Verbose, Eps = Eps)
    AllDataBotsPosRe  = tmpResPar$AllDataBotsPosRe
    AllDataBotsPosIm  = tmpResPar$AllDataBotsPosIm
    CourseOfHappiness = tmpResPar$CourseOfHappiness
    RadiusPerEpoch    = tmpResPar$RadiusPerEpoch
  }else{
    AllDataBotsPos = AllDataBotsPosRe + 1i * AllDataBotsPosIm
    tmpResSeq = PswarmEpochsSequential(AllDataBotsPos, MyDistanceMatrix, IndPossibleDBPosR, AllFreePosR0,
                                       NumAllDB, Lines, Columns, Origin, Happiness,
                                       GridRadii, GridAngle, QuadOrHexa, RadiusVector, Rmin, Rmax,
                                       Cls, Debug, pp, PlotIt = PlotIt, Verbose = Verbose)
    AllDataBotsPos    = tmpResSeq$AllDataBotsPos
    AllDataBotsPosRe  = Re(AllDataBotsPos)
    AllDataBotsPosIm  = Im(AllDataBotsPos)
    CourseOfHappiness = tmpResSeq$CourseOfHappiness
    RadiusPerEpoch    = tmpResSeq$RadiusPerEpoch
  }
  
  # Annomalies: zero indices are not allowed!
  # Take care of zeros
  ShiftReIdx = which(AllDataBotsPosRe == 0)
  ShiftImIdx = which(AllDataBotsPosIm == 0)
  if(length(ShiftReIdx) > 0){
    AllDataBotsPosRe[ShiftReIdx] = Lines
  }
  if(length(ShiftImIdx) > 0){
    AllDataBotsPosIm[ShiftImIdx] = Columns
  }
  
  T_DBS_2 = Sys.time()
  
  if(Verbose >= 1){
    t_diff_2 = T_DBS_2 - T_DBS_1
    print(paste0("Swarm computation ended (Time: ", round(t_diff_2, 2), " ", units(t_diff_2), ")."))
    print("Adjusting ...")
  }
  
  #----------------------------------------------------------------------------#
  # End of simulated annealing scheme
  #----------------------------------------------------------------------------#
  bmu = getCartesianCoordinates(DataBotsPosRe = AllDataBotsPosRe,
                                DataBotsPosIm = AllDataBotsPosIm,
                                GridRadius    = GridRadii,
                                GridAngle     = GridAngle,
                                QuadOrHexa    = QuadOrHexa)
  ShiftReIdx = which(bmu[,1] == 0)
  ShiftImIdx = which(bmu[,2] == 0)
  
  if(length(ShiftReIdx) > 0){
    bmu[ShiftReIdx, 1] = Lines # Columns
  }
  if(length(ShiftImIdx) > 0){
    bmu[ShiftImIdx, 2] = Columns # Lines
  }
  
  #----------------------------------------------------------------------------#
  # Possible Minor Rounding error
  if(max(bmu[,1])>Lines){
    Lines=Lines+1
  }
  if(max(bmu[,2])>Columns){
    Columns=Columns+1
  }
  #----------------------------------------------------------------------------#
  return(list("ProjectedPoints" = round(bmu, 0),
              "LC"              = c(Lines, Columns),
              "Control"         = list("CourseOfHappiness"  = CourseOfHappiness,
                                       "RadiusPerEpoch"     = RadiusPerEpoch,
                                       "HappinessLastEpoch" = Happiness)))
}

