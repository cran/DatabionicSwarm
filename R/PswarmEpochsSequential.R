PswarmEpochsSequential = function(AllDataBotsPos, MyDistanceMatrix, IndPossibleDBPosR, AllFreePosR0,
                                  NumAllDB, Lines, Columns, Origin, Happiness,
                                  GridRadii, GridAngle, QuadOrHexa, RadiusVector, Rmin, Rmax,
                                  Cls, Debug, pp, PlotIt = FALSE, Verbose = 1){
  
  NOfSteps  = length(RadiusVector)
  dummy = 0
  # Z Start algorithm
  if(Verbose >= 2){
    progress = txtProgressBar(min = dummy, max = NOfSteps+1, style = 3)
    #ProzentualeZeitfolge=round(sort((RadiusVector-Rmin)/(Rmax-Rmin),decreasing=F)*100,0)
    #ProzentualeZeitfolge[NOfSteps]=99
  }
  NumIter = 0
  CourseOfStress = c()
  RadiusPerEpoch = c()
  #----------------------------------------------------------------------------#
  # Simulated Annealing Scheme
  #----------------------------------------------------------------------------#
  for(Radius in RadiusVector){
    dummy = dummy+1
    if(Verbose >= 2){
      progressm = setTxtProgressBar(progress, dummy)
      #print(paste0('Operator: ', ProzentualeZeitfolge[dummy],'% calculated.'))
      #print(paste0('Operator: Current focus: ', Radius))
    }
    Jumping = TRUE
    limit = ceiling(1 / pp[Radius])    # Ab welcher Eppoche wird Abbruchbedingung geprueft
    steigungsverlaufind = 20           # limit #wieviele zurueckliegende eppochen werden maximal geprueft
    if(PlotIt){
      AllDataBotsPosRe  = Re(AllDataBotsPos)
      AllDataBotsPosIm  = Im(AllDataBotsPos)
      if(dummy == 1){
        Title = "Random initialization"
      }else{
        Title = paste0("Epoch ", dummy-1, " - Radius ", Radius, " - ", "% - NumIter: ", NumIter)
        #Title = paste0("Epoch ", dummy-1, " - Radius ", Radius, " - Improvement: ", Improvement, "% - NumIter: ", NumIter)
      }
      my_plot_function(AllDataBotsPosRe, AllDataBotsPosIm, GridRadii, GridAngle,
                       QuadOrHexa, Title, Cls)
    }
    #Zeitfaktor: Je naeher DataBots springen, desto schneller riechen sie erneut entspricht, weniger DBs springen pro Eppoche
    nBots = round(pp[Radius] * NumAllDB)
    # nBots = round(0.05*NumDataBots)    # s. Fast and reliable ESOM learning

    List = PswarmRadiusSequential(AllDataBotsPos, Radius, MyDistanceMatrix,
                                  IndPossibleDBPosR, AllFreePosR0,
                                  pp, Origin, Lines, Columns, nBots, limit,
                                  steigungsverlaufind, Happiness, Debug)
    
    AllDataBotsPos = List$AllDataBotsPos
    CourseOfStress = c(CourseOfStress, List$stressverlauf)
    RadiusPerEpoch = c(RadiusPerEpoch, List$fokussiertlaufind)
    NumIter        = length(List$fokussiertlaufind)
    
    if(Verbose >= 3){
      print(paste0('Operator: ', tail(RadiusPerEpoch, 1), '.iteration'))
      print(paste0('Operator: weak Nash equilibrium found. Payoff maximized with ',
                   signif(RelativeDifference(List$stressverlauf[1],tail(List$stressverlauf, 1)),2),' %'))
    }
  }
  
  if(Verbose >= 2){
    progressm = setTxtProgressBar(progress, dummy+1)
    close(progress)
    #print(paste0('Operator: 100 % calculated.'))
  }
  
  return(list("AllDataBotsPos"    = AllDataBotsPos,
              "CourseOfHappiness" = CourseOfStress,
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
