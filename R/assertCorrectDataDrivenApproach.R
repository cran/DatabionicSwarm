assertCorrectDataDrivenApproach = function(Lines, Columns, NumDataBots){
  # Assertions: Helper function
  # Assertion2 = assertCorrectDataDrivenApproach(Lines, Columns, LinesColumns)
  if(Lines %% 2 != 0){
    warning('Lines has to be even')
  }
  if(Columns %% 4 != 0){
    stop('Columns has to be dividable by four')
  }
  #achtung unter umstaenden darf auch Lines==Columns nicht gelten, muss ich noch pruefen
  # in setPolarGrid
  if(Lines > Columns){
    stop('Number of Columns has to be higher than number of Lines')
  }
  if(Lines * Columns < NumDataBots){
    stop('Map is too small for the Number of DataBots, please chose higher Lines and Columns')
  }
  return()
}
