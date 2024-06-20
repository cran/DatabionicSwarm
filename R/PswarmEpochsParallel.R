PswarmEpochsParallel = function(AllDataBotsPosRe, AllDataBotsPosIm, MyDistanceMatrix, AllFreePosR0,
                                GridRadii, GridAngle, JumpsPerRadius,
                                NumJumps, NumAllDB, Lines, Columns, Origin, Happiness,
                                QuadOrHexa, RadiusVector, Rmin, Rmax,
                                Cls, Debug, pp, PlotIt = FALSE, Verbose = 1, Eps = 0.0001){
  
  NumFreeShape1      = dim(AllFreePosR0)[1]
  DataDists          = as.vector(MyDistanceMatrix)
  DataBotsPos        = c(rep(AllDataBotsPosRe, NumJumps + 1), rep(AllDataBotsPosIm, NumJumps + 1))
  DataBotsPos        = c(DataBotsPos, rep(0, NumAllDB))
  AllallowedDBPosR0  = as.vector(AllFreePosR0)
  #IndPossibleDBPosRe = Re(IndPossibleDBPosR)
  #IndPossibleDBPosIm = Im(IndPossibleDBPosR)
  
  IdxRe1 = (NumJumps * NumAllDB) + 1
  IdxRe2 = ((NumJumps + 1) * NumAllDB)
  IdxIm1 = (NumJumps * NumAllDB) + 1 + (NumJumps + 1) * NumAllDB
  IdxIm2 = ((NumJumps + 1) * NumAllDB) + (NumJumps + 1) * NumAllDB
  
  if(PlotIt){
    AllDataBotsPosRe = DataBotsPos[IdxRe1:IdxRe2]
    AllDataBotsPosIm = DataBotsPos[IdxIm1:IdxIm2]
    Title = "Random initialization"
    my_plot_function(AllDataBotsPosRe, AllDataBotsPosIm, GridRadii, GridAngle, QuadOrHexa, Title, Cls)
  }
  
  CourseOfHappiness = c()
  RadiusPerEpoch    = c()
  CurrEpoch         = 0
  NOfSteps          = length(RadiusVector)
  dummy             = 0
  
  # Z Start algorithm
  if(Verbose >= 2){
    progress = txtProgressBar(min = dummy, max = NOfSteps+1, style = 3)
    #ProzentualeZeitfolge=round(sort((RadiusVector-Rmin)/(Rmax-Rmin),decreasing=F)*100,0)
    #ProzentualeZeitfolge[NOfSteps]=99
  }
  #----------------------------------------------------------------------------#
  # Simulated Annealing Scheme
  #----------------------------------------------------------------------------#
  for(Radius in RadiusVector){
    CurrEpoch = CurrEpoch + 1
    dummy = dummy+1
    if(Verbose >= 2){
      progressm = setTxtProgressBar(progress, dummy)
      #print(paste0('Operator: ', ProzentualeZeitfolge[dummy],'% calculated.'))
      #print(paste0('Operator: Current focus: ', Radius))
    }
    Jumping = TRUE
    MinIterations = ceiling(1 / pp[Radius])    # Ab welcher Eppoche wird Abbruchbedingung geprueft
    #MinIterations = 20
    steigungsverlaufind = 20                   # limit #wieviele zurueckliegende eppochen werden maximal geprueft
    
    #Zeitfaktor: Je naeher DataBots springen, desto schneller riechen sie erneut entspricht, weniger DBs springen pro Eppoche
    NumChoDB = round(pp[Radius] * NumAllDB)
    # nBots = round(0.05*NumDataBots)    # s. Fast and reliable ESOM learning
    
    IndPossibleDBPosRe = Re(JumpsPerRadius[[CurrEpoch]])
    IndPossibleDBPosIm = Im(JumpsPerRadius[[CurrEpoch]])
    
    List = PswarmRadiusParallel(DataBotsPos, DataDists, AllallowedDBPosR0,
                                IndPossibleDBPosRe, IndPossibleDBPosIm,
                                Lines, Columns, Radius, NumAllDB, NumChoDB,
                                NumFreeShape1, NumJumps, Origin[1], Origin[2],
                                Happiness, MinIterations, steigungsverlaufind, Eps, Debug)
    DataBotsPos       = List$AllDataBotsPos
    CoursePerEpoch    = List$HappinessCourse
    Improvement       = signif(RelativeDifference(CoursePerEpoch[1], tail(CoursePerEpoch, 1)), 2)
    CourseOfHappiness = c(CourseOfHappiness, CoursePerEpoch)
    RadiusPerEpoch    = c(RadiusPerEpoch, List$fokussiertlaufind)
    
    if(PlotIt){
      AllDataBotsPosRe = DataBotsPos[IdxRe1:IdxRe2]
      AllDataBotsPosIm = DataBotsPos[IdxIm1:IdxIm2]
      Title = paste0("Epoch ", CurrEpoch, " - Radius ", Radius, " - Improvement: ", Improvement, "% - NumIter: ", RadiusPerEpoch[CurrEpoch])
      my_plot_function(AllDataBotsPosRe, AllDataBotsPosIm, GridRadii, GridAngle,
                       QuadOrHexa, Title, Cls)
    }
    
    if(Verbose >= 3){
      print(paste0('Operator: ', tail(RadiusPerEpoch, 1), '.iteration'))
      print(paste0('Operator: weak Nash equilibrium found. Payoff maximized with ', Improvement,' %'))
    }
  }
  #----------------------------------------------------------------------------#
  if(Verbose >= 2){
    progressm = setTxtProgressBar(progress, dummy+1)
    close(progress)
    #print(paste0('Operator: 100 % calculated.'))
  }
  
  AllDataBotsPosRe = DataBotsPos[IdxRe1:IdxRe2]
  AllDataBotsPosIm = DataBotsPos[IdxIm1:IdxIm2]
  
  return(list("AllDataBotsPosRe"  = AllDataBotsPosRe,
              "AllDataBotsPosIm"  = AllDataBotsPosIm,
              "CourseOfHappiness" = CourseOfHappiness,
              "RadiusPerEpoch"    = RadiusPerEpoch))
}

my_plot_function = function(AllDataBotsPosRe, AllDataBotsPosIm, GridRadii,
                            GridAngle, QuadOrHexa, Title, Cls){
  bmu = getCartesianCoordinates(DataBotsPosRe = AllDataBotsPosRe,
                                DataBotsPosIm = AllDataBotsPosIm,
                                GridRadius       = GridRadii,
                                GridAngle        = GridAngle,
                                QuadOrHexa       = QuadOrHexa)
  plotSwarm(bmu, Cls, main = Title)
}
