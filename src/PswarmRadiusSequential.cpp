// [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadilloExtensions/sample.h>
#include <tuple>
#include "sampleC.h"

using namespace Rcpp;
using namespace std;
using namespace sugar;


// // [[Rcpp::depends(RcppArmadillo)]]
//NumericVector rcppPar_sampleC(NumericVector x, double len) {
//  bool replace=0;
//  return(RcppArmadillo::sample(x,len,replace));
//}

// [[Rcpp::depends(RcppArmadillo)]]
NumericVector lmC(NumericVector x,NumericVector yr) {
  // doku see lm() in R 
  NumericMatrix Xr(yr.length(),2);
  for(int i=0;i<yr.length();i++){
    Xr(i,0)=1;
    Xr(i,1)=x(i);
  }
  int n = Xr.nrow(), k = Xr.ncol();
  
  arma::mat X(Xr.begin(), n, k, false);       // reuses memory and avoids extra copy
  arma::colvec y(yr.begin(), yr.size(), false);
  
  arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
  
  NumericVector coefs(coef.begin(),coef.end());
  
  return(coefs);
}


// [[Rcpp::depends(RcppArmadillo)]]
Rcomplex modComplexC(Rcomplex x,Rcomplex y) {
  // Berechnet den Modulo zweier komplexer Zahlen
  Rcomplex res;
  //res.i=x.i-((int)x.i/y.i)*y.i;
  //res.r=x.r-((int)x.r/y.r)*y.i;
  res.i=(int)abs(x.i)%(int)y.i;
  res.r=(int)abs(x.r)%(int)y.r;
  return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
Rcomplex makePolarPositionToroidC(Rcomplex DBpositionInd,Rcomplex IndPossibleDBPosR,double Lines, double Columns) {
  //INPUT
  // DBpositionInd                 complex number, Current row and column Indize of the DataBot  
  // IndPossibleDBPosR[        complex number of Indice of radius of AllallowedDBPosR0 of one possible DataBots Position in one radial jump  
  // Lines         Integer, hast to be able to be divided by 2
  // Columns       Integer, with Columns>Lines  
  // Optional
  // ToroidPositionInd       complex number of Indice of radius of all one DataBot Position in one radial jump in Packman-Universe
  // Output:
  
  // author: MT 02/2016  
  //Verschiebe alle moeglichen Positionen B bis zum DB Ort A);
  //ComplexVector ToroidPositionInds(n);
  Rcomplex ToroidPositionInd=IndPossibleDBPosR+DBpositionInd;
  
  // IndPossibleDBPosR[,1]=IndPossibleDBPosR[,1]+DBpositionInd[1]//-(Lines/2+1)# And den Ursprung des anderen Bezug-Systemes (BS) anpassen
  // IndPossibleDBPosR[,2]=IndPossibleDBPosR[,2]+DBpositionInd[2]//-(Lines/2+1)# And den Ursprung des anderen Bezug-Systemes anpassen
  //anderes BS, s.  ind=which(AllallowedDBPosR0<=Radius & AllallowedDBPosR0!=0,arr.ind=T) in getOriginPositionsByRadius()
  //  Toroides Feld anpassen
  //Trick: u.U. sind Positionen negativ und werden es mit Modulo positiv angepasst
  //cout <<db<<endl;
  //cout<<DBpositionInd(db)<<endl;
  //cout<<ToroidPositionInds(0)<<endl;
  Rcomplex eins;
  eins.i=1;
  eins.r=1;
  Rcomplex LC;
  LC.i=Columns;
  LC.r=Lines;
  
  //Der Eins Trick hat etwas mit dem Modulo bei negativen Zahlen zu tun, allerdings weis ich nichtmehr wieso
  ToroidPositionInd=modComplexC(ToroidPositionInd-eins,LC)+eins;
  
  //cout<<ToroidPositionInds<<endl;
  //IndPossibleDBPosR[,1]=(IndPossibleDBPosR[,1]-1)%Lines+1
  //IndPossibleDBPosR[,2]=(IndPossibleDBPosR[,2]-1)%Columns+1
  return ToroidPositionInd;
}


// [[Rcpp::depends(RcppArmadillo)]]
double Happiness4BotC(NumericVector DistInput,NumericVector DistOutput,double Radius,NumericVector Nachbahrschaftsfunktion, double DBAnzahl, double Happiness){
  // HappinessVector=stresskriterium(DistsInput,DistsOutput,Radius)
  // Stresskriterium fuer den sprung eines Databots berechnet fuer alle Databots seperat
  //
  // INPUT
  // DistsInput       vector: Pairwise distance between pairs of objects
  // DistsOutput      vector of polar projected points: Pairwise distance between pairs of objects
  // Radius          Radius der Umgebung fuer Nachbahrschaftsfunktion 
  //
  // OUTPUT
  // HappinessVector         One value
  //  
  // Autor: MT 01/2015
  // 1. Editor: CL 01/15
  // 2. Editor: MT 02/15
  
  
  //Nachbahrschaftsfunktion=1-DistOutput*DistOutput/(3.14159265*Radius*Radius);
  // double x;
  // double N=0;
  // for(int i=0;i<DBAnzahl;i++){
  //   x=1-DistOutput[i]*DistOutput[i]/(3.14159265*Radius*Radius);
  //   if(x<0){
  //     Nachbahrschaftsfunktion[i]=0;
  //   }else{
  //     Nachbahrschaftsfunktion[i]=x;
  //     N=N+x;
  //   }
  // }
  Nachbahrschaftsfunktion=1-(DistOutput*DistOutput)/(3.14159265*Radius*Radius);
  //Nachbahrschaftsfunktion[Nachbahrschaftsfunktion<0]=0;
  for(int i=0;i<DBAnzahl;i++){
    if(Nachbahrschaftsfunktion[i]<0)
      Nachbahrschaftsfunktion[i]=0;
  }
  //double N=std::accumulate(Nachbahrschaftsfunktion.begin(),Nachbahrschaftsfunktion.end());
  
  double N=sum(Nachbahrschaftsfunktion);
  
  if(N<=0.0000001){return(Happiness);}
  
  
  return(Happiness-sum(Nachbahrschaftsfunktion*DistInput)/N);
  
}

// [[Rcpp::depends(RcppArmadillo)]]
double vecminInd(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference 
  return it - x.begin();
}

// [[Rcpp::depends(RcppArmadillo)]]
double vecmaxInd(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::max_element(x.begin(), x.end());
  // we want the value so dereference 
  return it - x.begin();
}

// [[Rcpp::depends(RcppArmadillo)]]
NumericMatrix calcHappinessC(NumericMatrix DataDists,NumericMatrix OutputDistance,NumericMatrix OutputDistanceNeu,NumericMatrix OutputDistanceNeu2,NumericMatrix OutputDistanceNeu3,NumericMatrix OutputDistanceNeu4,double Radius,double Happiness, double DBAnzahl,NumericVector Nachbahrschaftsfunktion,NumericMatrix xx){
  
  //double DBAnzahl=DataDists.nrow();
  //NumericVector Nachbahrschaftsfunktion(DBAnzahl);
  //NumericMatrix xx(DBAnzahl, 3);
  double phi;
  //double phiNeu;
  //double phiNeu2;
  double phiTest;
  //double pos;
  for(int db =0;db<DBAnzahl;db++){  //espringen nur 15% der DataBots durch BotsJumping
    // aber es wird trotzdem fuer alle der HappinessVector neu berechnet//
    //Achtung Distanz zu sich selber auslassen!
    // phi=stress4Bot(DataDists[db,-db],OutputDistance[db,-db],Radius)  
    // phiNeu=stress4Bot(DataDists[db,-db],OutputDistanceNeu[db,-db],Radius)
    
    
    phi=Happiness4BotC(DataDists.row(db),OutputDistance.row(db),Radius,Nachbahrschaftsfunktion,DBAnzahl,Happiness) ; 
    NumericVector PhiNeu(4);
    PhiNeu[0]=Happiness4BotC(DataDists.row(db),OutputDistanceNeu.row(db),Radius,Nachbahrschaftsfunktion,DBAnzahl,Happiness);
    PhiNeu[1]=Happiness4BotC(DataDists.row(db),OutputDistanceNeu2.row(db),Radius,Nachbahrschaftsfunktion,DBAnzahl,Happiness);
    PhiNeu[2]=Happiness4BotC(DataDists.row(db),OutputDistanceNeu3.row(db),Radius,Nachbahrschaftsfunktion,DBAnzahl,Happiness);
    PhiNeu[3]=Happiness4BotC(DataDists.row(db),OutputDistanceNeu4.row(db),Radius,Nachbahrschaftsfunktion,DBAnzahl,Happiness);
    
    double ind=vecmaxInd(PhiNeu);
    // if(phiNeu<phiNeu2){
    //   phiTest=phiNeu2;
    //   pos=2;
    // }else{
    //     phiTest=phiNeu;
    // pos=1;
    // }  
    phiTest=PhiNeu[ind];
    if(phiTest>phi)
    {
      xx(db,2)=ind;
      xx(db,1)=db;
      xx(db,0)=phiTest;
      
    } //end if phiNeu<phi
    else
    {
      xx(db,2)=0;
      xx(db,1)=-1;
      xx(db,0)=phi;
      
    }
    //xx(db,2)=phi;
    //xx(db,3)=phiNeu;
  } // end for 1:DBAnzahl
  return(xx);
}


// [[Rcpp::depends(RcppArmadillo)]]
ComplexVector calcPolarPositionsC(ComplexVector DataBotsPos,NumericVector ChosenForJump,ComplexVector PossiblePositions,double Radius,double Lines,double Columns,ComplexVector ToroidPosition, double db, int n,ComplexVector DataBotsPosNeu){
  // calcPolarPositionsGauss(DataBotsPos,RadiusPositionsschablone,Radius,Lines,Columns)
  // Position random generation from the normal distribution in an two dimensional toroid grid defined by Lines and Columns
  //INPUT
  // DataBotsPos[1:AnzData,2]        complex vector of Two Indizes per Databot such that getCartesianCoordinates(DataBotPos,GridRadius,GridAngle) gets the cartesian positions on the grid
  // ChosenForJump                              numeric vector of DataBots, which where chosen to jump, s. Fast and reliable ESOM learning
  // PossiblePositions[m,2]                    complex vector of two indizes of possible jum positions
  // Radius                                        Jump radius of DataBot, around this DataBot, for each DataBots equal
  // Lines                           Integer, hast to be able to be divided by 2, sets Size of planar grid
  // Columns                         Integer, with Columns>Lines  sets Size of planar grid
  // Output:
  // DataBotsPos[1:AnzData]          random new Indizes respective to polar dataPoints positions on free grid places within a radius
  // author: MT 02/2016
  
  //ComplexVector ClosedPositions=clone(DataBotsPos) ;//Hier kommen alle Positionen rein, welche belegt sind, d.h. dorthin darf nicht gesprungen werdem
  //DataBotsPosNeu=clone(DataBotsPos);
  //ist anscheinend schneller als clone
  std::copy( DataBotsPos.begin(), DataBotsPos.end(), DataBotsPosNeu.begin() ) ;
  //IntegerVector DataBotsPosRealInt=as<IntegerVector>(DataBotsPosReal);
  //IntegerVector DataBotsPosImagInt=as<IntegerVector>(DataBotsPosImag);
  
  // deep copy, sonst diverse fehler
  // ComplexVector DataBotsPosClone = clone(DataBotsPos);
  // Muss abhaengig vom gesamtradius sein => bei kleinen Radien stehen viel weniger meogliche positionen zur verfuegung
  //int n=ChosenForJump.length();
  ///// wird nun einmalig vorher initialisiert
  //ComplexVector ToroidPosition(n);
  //double db;
  ////
  for(int i=0;i<n;i++){
    db=ChosenForJump[i];//R zaehlt um 1 anders wie C++
    Rcomplex PosOld=DataBotsPos[db];
    Rcomplex PosNew=PossiblePositions[i];
    //Nullpunkt Verschiebung der erlaubten Indizes bezueglich eines polaren gitters von Mitte zum linken unteren Ecke
    // die indizes dieser erlaubten positionen muessen nun an ein toroides gitter angepasst werden
    //ToroidPositionIndsCurrent=makePolarPositionsToroid(cbind(Re(DataBotsPos[db]),Im(DataBotsPos[db])),cbind(Re(IndPossibleDBPosR),Im(IndPossibleDBPosR)),Lines,Columns)
    ////ToroidPositionIndsCurrent=makePolarPositionsToroidC(DataBotsPos,db-1,PossiblePositions[i],Lines,Columns);
    Rcomplex ToroidPositionIndsCurrent= makePolarPositionToroidC(PosOld,PosNew,Lines,Columns);
    // besetzte Positionen muessen entfernt werden
    //OpenPositions=setdiffComplexVectors(ToroidPositionIndsCurrent,ClosedPositions)
    ;
    ToroidPosition[i]=ToroidPositionIndsCurrent;
  }
  //aehnlich zu R setdiff(a,b)
  int i;int j;
  IntegerVector Res(ToroidPosition.length());
  for(i=0; i < ToroidPosition.size();i++){
    for(j=0; j < DataBotsPos.size();j++){
      if((ToroidPosition[i].i == DataBotsPos[j].i) &&
         (ToroidPosition[i].r == DataBotsPos[j].r))
        Res[i] = 1;
    }
  }
  //Achtung es koennten theoretisch 2 DataBots auf die gleiche neue Position springen, wenn
  // a) die Gleiche Position in 2 verscheidenen samplen gewuerfelt worden ist
  // b) diese Position fuer beide DataBots einen besseren HappinessVector hat als jeweils 3 andere Positionen
  for(int i=0;i<n;i++){
    db=ChosenForJump[i];//R zaehlt um 1 anders wie C++
    // ComplexVector tzz(1);
    // tzz=ToroidPosition[i];
    //ComplexVector OpenPositions=setdiffComplexVectorC(tzz,DataBotsPosRealInt,DataBotsPosImagInt);
    //aus den offenen positionen wird eine position gezogen und die outputvarible dort veraendert
    //int lengthDbsindtest=OpenPositions.length();
    if(Res[i]==0){
      //DataBotsPosClone[db]=OpenPositions[0];
      
      DataBotsPosNeu[db]=ToroidPosition[i];//DataBot darf weder auf alte Positions eines anderen DBs springen
      //ClosedPositions.push_back(OpenPositions[0]);
    } // end iflengthDbsindtest>0
  } // end for db in vec
  return(DataBotsPosNeu);
}


// [[Rcpp::depends(RcppArmadillo)]]
NumericMatrix rDistanceToroidC(NumericVector AllDataBotsPosX, NumericVector AllDataBotsPosY,NumericMatrix AllallowedDBPosR0,double Lines,double Columns, NumericVector Nullpunkt,double DBanzahl,NumericMatrix Distances,NumericVector Dx,NumericVector Dy,NumericVector D1,NumericVector D2){
  
  //double DBanzahl=AllDataBotsPosX.length();
  //NumericMatrix Distances(DBanzahl,DBanzahl); //
  //NumericVector Dx(DBanzahl);
  //NumericVector Dy(DBanzahl);
  //NumericVector D1(DBanzahl);
  //NumericVector D2(DBanzahl);
  for(int i=0;i<DBanzahl;i++){
    Dx=abs(AllDataBotsPosX-AllDataBotsPosX[i]);
    Dy=abs(AllDataBotsPosY-AllDataBotsPosY[i]);  
    
    D1=Lines-Dx+1;
    D2=Columns-Dy+1;
    // -1 da im Gegensatz zu R die MAtrixenindizes bei 0 beginnen!
    Dx = pmin(Dx,D1)+Nullpunkt[0]-1; //Toroid machen und an Nullpunkt der Schablone anpassen
    Dy = pmin(Dy,D2)+Nullpunkt[1]-1; //Toroid machen und an Nullpunkt der Schablone anpassen
    
    for(int j=0;j<DBanzahl;j++){
      Distances(i,j)=AllallowedDBPosR0(Dx[j],Dy[j]);
    }
  }
  return Distances;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List PswarmRadiusSequential(ComplexVector AllDataBotsPosOld,
                                       double Radius,
                                       NumericMatrix DataDists,
                                       ComplexVector IndPossibleDBPosR,
                                       NumericMatrix RadiusPositionsschablone,
                                       NumericVector pp,
                                       NumericVector Nullpunkt,
                                       double Lines,
                                       double Columns,
                                       double nBots,
                                       int limit,
                                       int steigungsverlaufind, 
                                       double Happiness, 
                                       bool debug){
  // PswarmRadiusSequential( AllDataBotsPosOld,
  //                                      Radius, DataDists,
  //                                      IndPossibleDBPosR,
  //                                      RadiusPositionsschablone,  pp,
  //                                      Nullpunkt, Lines,  Columns,
  //                                      nBots,  limit, steigungsverlaufind,  Happiness)
  //   intern function, do not use yourself
  //   Finds the weak Nash equilibirium for DataBots in one epoch(Radius), requires the setting of constants, grid, and so on in \code{\link{pswarmCpp}}
  // INPUT
  //   AllDataBotsPosOld              ComplexVector [1:n,1], DataBots position in the last Nash-Equlibriuum}
  //   Radius                         double, Radius of payoff function, neighborhood, where other DatsBots can be smelled}
  //   DataDists                      NumericMatrix, Inputdistances[1:n,1:n]}
  //   IndPossibleDBPosR              ComplexVector, see output of \code{\link{findPossiblePositionsCsingle}}}
  //   RadiusPositionsschablone       NumericMatrix, see \code{AllallowedDBPosR0} in \code{\link{setPolarGrid}}}
  //   pp                             NumericVector, number of jumping simultaneously DataBots of one eppoch (per nash-equilibirum), this vector is linearly monotonically decreasing}
  //   Nullpunkt                      NumericVector, equals \code{which(AllallowedDBPosR0==0,arr.ind=T)}, see see \code{AllallowedDBPosR0} in \code{\link{setPolarGrid}}}
  //   Lines                          double, small edge length of rectangulare grid}
  //   Columns                        double, big edge length of rectangulare grid}
  //   nBots                          double, intern constant, equals \code{round(pp[Radius]*DBAnzahl)}}
  //   limit                          int, intern constant, equals \code{ceiling(1/pp[Radius])}}
  //   steigungsverlaufind            int, intern constant}
  //   Happiness              double, intern constant, sum of payoff of all databots in random condition before the algorithm starts}
  // OUTPUT:
  //   list V of
  //   V$AllDataBotsPos           ComplexVector, indizes of DataBot Positions after a weak Nash equlibrium is found}
  //   V$HappinessCourse            NumericVector, intern result, for debugging only}
  //   V$fokussiertlaufind        NumericVector, intern result, for debugging only}
  // 
  // Author: Michael Thrun, 04/16
  // Edit: Quirin Stier, 03/24
  
  ComplexVector  AllDataBotsPos=clone(AllDataBotsPosOld);
  bool Jumping          = 1;
  int fokussiertlaufind = 0;
  int leng              = AllDataBotsPos.length();
  int Iteration         = 0;
  //double DBanzahl=leng;
  NumericVector HappinessVector(leng);         // Numeric vector bookkeeping the current best local happiness for each databot
  NumericVector slopeVec(2);                   // Bookkeeping variable for inclination (regression) of global happiness course
  NumericVector KeyBot(leng);                  //dummy
  NumericMatrix PosAndHappy(leng,2);           // Bookkeeping matrix containing index of best choice (index: current state + 4 possible jumps = 5 choices) and its positioning
  ComplexVector DataBotsPosNeu(leng);          // New positions in complex vector
  ComplexVector DataBotsPosNeu2(leng);         // New possible positions for each
  ComplexVector DataBotsPosNeu3(leng);         // possible jump 
  ComplexVector DataBotsPosNeu4(leng);         // (total 4 jumps fix)
  NumericMatrix OutputDistanceNeu(leng,leng);  // New distance matrix of databots
  NumericMatrix OutputDistanceNeu2(leng,leng); // on polar projected grid
  NumericMatrix OutputDistanceNeu3(leng,leng); // One for each possible jump => 4 jumps
  NumericMatrix OutputDistanceNeu4(leng,leng); // => 4 new distance matrices
  NumericMatrix OutputDistance(leng,leng);     // Distance matrix of current state - distance matrix of projection
  
  //Initialisierung calcPolarPositions und calcSressC, calcdistancetoroid
  ComplexVector ToroidPosition(nBots);
  double db=-1;
  int nBotsalsInt=nBots;
  
  double DBAnzahl=DataDists.nrow();
  NumericVector Nachbahrschaftsfunktion(DBAnzahl);
  NumericMatrix xxVergleich(DBAnzahl, 3);
  
  NumericMatrix Distances(DBAnzahl,DBAnzahl); //Alternativ muesste in toroiddistances eine deep copy machen!
  NumericMatrix Distances1(DBAnzahl,DBAnzahl); // These objects are used to initialize an object for Rcpp computation methods
  NumericMatrix Distances2(DBAnzahl,DBAnzahl); // Here, its avoided to initialize objects within computational methods ...
  NumericMatrix Distances3(DBAnzahl,DBAnzahl); // This could be changed
  NumericMatrix Distances4(DBAnzahl,DBAnzahl);
  
  NumericVector Dx(DBAnzahl);
  NumericVector Dy(DBAnzahl);
  NumericVector D1(DBAnzahl);
  NumericVector D2(DBAnzahl);

  //NumericVector CurrentKeyBot(leng);
  NumericVector ChosenForJump(nBots);
  NumericVector AllDataBotsPosReal(DBAnzahl);
  NumericVector AllDataBotsPosImag(DBAnzahl);
  
  NumericVector PossiblePositionsIndi(nBots);
  NumericVector PossiblePositionsIndi2(nBots);
  NumericVector PossiblePositionsIndi3(nBots);
  NumericVector PossiblePositionsIndi4(nBots);
  for(int i=0;i<leng;i++){
    KeyBot(i)=i;
    //CurrentKeyBot(i)=i;
  }
  
  int KeyPossiblePositionLen=IndPossibleDBPosR.length();
  NumericVector KeyPossiblePosition(KeyPossiblePositionLen);
  ComplexVector PossiblePositions(nBots);
  ComplexVector PossiblePositions2(nBots);
  ComplexVector PossiblePositions3(nBots);
  ComplexVector PossiblePositions4(nBots);
  for(int i=0;i<KeyPossiblePositionLen;i++){
    KeyPossiblePosition(i)=i;
  }
  
  NumericVector HappinessCourse;
  NumericVector KeySteigung(steigungsverlaufind);
  NumericVector HappinessTail(steigungsverlaufind);
  double epsilon=0.0001;//Steigungsgenauigkeit, steigung nimmt kaum mehr ab
  for(int i=0;i<steigungsverlaufind;i++)
    KeySteigung(i)=i;
  
  while(Jumping){
    fokussiertlaufind=fokussiertlaufind+1;
    //NumericVector AllDataBotsPosReal(leng);
    //NumericVector AllDataBotsPosImag(leng);
    // for(int i=0;i<leng;i++){
    //    AllDataBotsPosReal(i)=AllDataBotsPos(i).r;
    //    AllDataBotsPosImag(i)=AllDataBotsPos(i).i;
    // }
    AllDataBotsPosReal=Re(AllDataBotsPos);
    AllDataBotsPosImag=Im(AllDataBotsPos);
    
    // Computation of the current distance between databots on 2d projected grid
    OutputDistance=rDistanceToroidC(AllDataBotsPosReal,AllDataBotsPosImag,RadiusPositionsschablone,Lines,Columns,Nullpunkt,DBAnzahl,Distances,Dx,Dy,D1,D2);
    ////// Bestimme Anzahl Position und Distanz Springender Databots   // - //
    //Zeitfaktor: Je naeher DataBots springen, desto schneller riechen sie erneut entspricht, weniger DBs springen pro Iteration
    
    //nBots=round(0.05*DBAnzahl)//s. Fast and reliable ESOM learning#
    // Sample als ziehen ohne zuruecklegen!!
    //if(CurrentKeyBot.length()<nBots){
    //NumericVector CurrentKeyBot(leng);
    //CurrentKeyBot=clone(KeyBot);
    // }
    
    // Sample a few databots to allow these to look for a new position which
    // might be better
    ChosenForJump=sampleC(KeyBot,nBots); //15% Wahrscheinlichkeit das der DataBot ueberhaupt versucht zu springen
    
    //NumericVector CurrentKeyBot2=DeleteAll(CurrentKeyBot,ChosenForJump);
    //NumericVector CurrentKeyBot(CurrentKeyBot2.length());
    // CurrentKeyBot=clone(CurrentKeyBot2);
    
    //Normalerweise koennte man direkt aus de Vektor ein sample ziehen, allerdings geht die
    //sample funktion nur fuer numerischeVektoren, hier ist aber ein ComplexerVektor vorhanden
    
    // For each chosen databot which is allowed to jump, sample random movements
    PossiblePositionsIndi=sampleC(KeyPossiblePosition,nBots);
    PossiblePositionsIndi2=sampleC(KeyPossiblePosition,nBots);
    PossiblePositionsIndi3=sampleC(KeyPossiblePosition,nBots);
    PossiblePositionsIndi4=sampleC(KeyPossiblePosition,nBots);
    for(int i=0;i<nBots;i++){
      PossiblePositions(i)=IndPossibleDBPosR[PossiblePositionsIndi(i)];
      PossiblePositions2(i)=IndPossibleDBPosR[PossiblePositionsIndi2(i)];
      PossiblePositions3(i)=IndPossibleDBPosR[PossiblePositionsIndi3(i)];
      PossiblePositions4(i)=IndPossibleDBPosR[PossiblePositionsIndi4(i)];
    }
    //DataBotsPosNeu=calcPolarPositionsV3(AllDataBotsPos,ChosenForJump,PossiblePositions,Radius,Lines,Columns)
    // Indize abzug von 1 in ChosenForJump, da R von 1 und C++ von 0 zaehlt
    
    // Computation of new random positions for a subset of the total number of
    // databots. For each moving databots, 4 new positions are randomly drawn
    DataBotsPosNeu=calcPolarPositionsC(AllDataBotsPos,ChosenForJump,PossiblePositions,Radius,Lines,Columns, ToroidPosition, db,nBotsalsInt,DataBotsPosNeu); //DataBotsPosNeu als pointen hinten uebergen
    DataBotsPosNeu2=calcPolarPositionsC(AllDataBotsPos,ChosenForJump,PossiblePositions2,Radius,Lines,Columns, ToroidPosition, db,nBotsalsInt,DataBotsPosNeu2);
    DataBotsPosNeu3=calcPolarPositionsC(AllDataBotsPos,ChosenForJump,PossiblePositions3,Radius,Lines,Columns, ToroidPosition, db,nBotsalsInt,DataBotsPosNeu3);
    DataBotsPosNeu4=calcPolarPositionsC(AllDataBotsPos,ChosenForJump,PossiblePositions4,Radius,Lines,Columns, ToroidPosition, db,nBotsalsInt,DataBotsPosNeu4);
    
    // Computation of the distances of the new possible positions on the 2d projection grid
    OutputDistanceNeu=rDistanceToroidC(Re(DataBotsPosNeu),Im(DataBotsPosNeu),RadiusPositionsschablone,Lines,Columns,Nullpunkt,DBAnzahl,Distances1,Dx,Dy,D1,D2);
    OutputDistanceNeu2=rDistanceToroidC(Re(DataBotsPosNeu2),Im(DataBotsPosNeu2),RadiusPositionsschablone,Lines,Columns,Nullpunkt,DBAnzahl,Distances2,Dx,Dy,D1,D2);
    OutputDistanceNeu3=rDistanceToroidC(Re(DataBotsPosNeu3),Im(DataBotsPosNeu3),RadiusPositionsschablone,Lines,Columns,Nullpunkt,DBAnzahl,Distances3,Dx,Dy,D1,D2);
    OutputDistanceNeu4=rDistanceToroidC(Re(DataBotsPosNeu4),Im(DataBotsPosNeu4),RadiusPositionsschablone,Lines,Columns,Nullpunkt,DBAnzahl,Distances4,Dx,Dy,D1,D2);
    
    //////                                                             // - //
    ////// Vegleiche Phi                                               // - //
    // Achtung: Es werden aber als "Jumps", Alle DBs gezaehlt, welche sich in einer
    // neuen Position wohl fuehlen wuerden, selbst wenn sie keine moeglichkeit haben
    // auf diese neue Position zu springen
    
    // Compute the happiness for each position of current state and the four new
    // proposed positions where a few bots (everywhere the same) were allowed to change
    PosAndHappy=calcHappinessC(DataDists,OutputDistance,OutputDistanceNeu,OutputDistanceNeu2,OutputDistanceNeu3,OutputDistanceNeu4,Radius,Happiness,DBAnzahl,Nachbahrschaftsfunktion, xxVergleich);
    
    // Choose best happiness value and set the databot positioning
    int j=0;
    double SumHappiness=0;
    for(int i=0;i<leng ;i++){
      HappinessVector(i)=PosAndHappy(i,0);
      if(PosAndHappy(i,1)!=-1){
        if(PosAndHappy(i,2)==0)
          AllDataBotsPos(PosAndHappy(i,1))=DataBotsPosNeu(PosAndHappy(i,1));
        if(PosAndHappy(i,2)==1)
          AllDataBotsPos(PosAndHappy(i,1))=DataBotsPosNeu2(PosAndHappy(i,1));
        if(PosAndHappy(i,2)==2)
          AllDataBotsPos(PosAndHappy(i,1))=DataBotsPosNeu3(PosAndHappy(i,1));
        if(PosAndHappy(i,2)==3)
          AllDataBotsPos(PosAndHappy(i,1))=DataBotsPosNeu4(PosAndHappy(i,1));
        j++;
      }
      if(HappinessVector(i)!=R_PosInf){
        SumHappiness=SumHappiness+HappinessVector(i);
      }
    }
    
    HappinessCourse.push_back(SumHappiness);///sqrt(sqrt(Nomierung)))///(pi*Radius^2)*DBAnzahl)
    Iteration=Iteration+1;
    //////                                                               // - //
    ////// Abbruchbedingung                                              // - //
    if(j==0){
      //Jumping=0; //macht den Algorithmus nur 4mal schneller sonst nichts
    }// end if 
    
    if(Iteration>limit){ //Pruefe ab 20. Iteration in einem Radius
      HappinessTail=tail(HappinessCourse,min(steigungsverlaufind,Iteration)); //Innerhalb einer Radius oder nur letzten Iterationn?
      //slope=lm(formula=HappinessTail~x)$coefficients[2]
      slopeVec=lmC(KeySteigung,HappinessTail);
      
      if(slopeVec[1]<epsilon){ //globale HappinessVector nicht mehr besonders zu
        Jumping=0;
        //vielleicht stattdessen random walk mit 1% der dbs?
      } // end if slope>=0
    }// end of Iteration>limit
    //////                                                               // - //
    if(debug){
      if(Iteration%100==(limit+1)){
        //CRAN limits cout,therefore this is depricated
        //std::cout<<"Steigung:"<<slopeVec[1]<<"; Iteration:"<<Iteration<<"; HappinessVector:"<<SumHappiness<<";"<<std::endl;
      }
    }
  }// end While Jumping
  return Rcpp::List::create(
    Rcpp::Named("AllDataBotsPos")    = AllDataBotsPos,
    Rcpp::Named("HappinessCourse")     = HappinessCourse,
    Rcpp::Named("fokussiertlaufind") = fokussiertlaufind
  ) ;
}
