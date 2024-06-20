getCartesianCoordinates = function(DataBotsPosRe, DataBotsPosIm,
                                   GridRadius, GridAngle, QuadOrHexa = TRUE){
  # BMUs = getCartesianCoordinates(DataBotsPos,GridRadius,GridAngle)
  # Transform DataBot Indizes two exac cartesian coordinates on an toroid two dimensional grid.
  #
  # INPUT
  # DataBotsPos[1:AnzData]               complex vector Two Indizes per Databot describing its positions in an two dimensional grid
  # GridRadius[Lines,Columns]              Radii Matrix of all possible Positions of DataBots in Grid, see als Doku of setPolarGrid()
  # GridAngle[Lines,Columns]               Angle Matrix of all possible Positions of DataBots in Grid, see als Doku of setPolarGrid()
  # Lines,Columns                          Size of planar toroid two dimensional grid
  # QuadOrHexa                           FALSE=If DataPos on hexadiagonal grid, runde auf 2 nachkommestellen, Default=TRUE
  #
  # OUTPUT
  # BestMatchingUnits[1:AnzData,2]         coordinates on an two dimensional grid for each databot excluding unique key, sucht that
  #                                        with gUmatrix(inputs=Data, ,projectionPoints=BestMatchingUnits,cls=cls,toroid=TRUE), see Doku there
  #                                        a visualization of the pswarm projection is possible
  # author: MT 01/2015  
  # Example:
  # If Classification cls of Data is availible
  #
  # bmu=getCartesianCoordinates(DataBotsPos,GridRadius,GridAngle)
  # ClassPlot(bmu[1,],bmu[2,],cls)
  DataBotsPosInd = cbind(DataBotsPosRe, DataBotsPosIm)
  bmR            = GridRadius[DataBotsPosInd]
  bmPhi          = GridAngle[DataBotsPosInd]*pi/180    
  if(!QuadOrHexa){
    bmY = round(bmR*cos(bmPhi),2)
    bmX = round(bmR*sin(bmPhi),2)
  }else{
    bmY = bmR*cos(bmPhi)
    bmX = bmR*sin(bmPhi)
  }
  return(BestMatchingUnits=cbind(bmX, bmY))
}
