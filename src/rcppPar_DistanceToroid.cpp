// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>

using namespace RcppParallel;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppParallel)]]
// Worker for parallelization
struct PolarDistance : public Worker {
  // inputs to read from
  const RVector<double> AllDataBotsPosX;
  const RVector<double> AllDataBotsPosY;
  const RMatrix<double> AllallowedDBPosR0;
  const double Lines;
  const double Columns;
  const RVector<double> Nullpunkt;
  const int DBanzahl;
  
  // output to write to
  RMatrix<double> Distances;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to form the Rcpp matrix type)
  PolarDistance(const NumericVector AllDataBotsPosX,
                const NumericVector AllDataBotsPosY,
                const NumericMatrix AllallowedDBPosR0,
                const double Lines,
                const double Columns,
                const NumericVector Nullpunkt,
                const int DBanzahl,
                NumericMatrix Distances):
    AllDataBotsPosX(AllDataBotsPosX),
    AllDataBotsPosY(AllDataBotsPosY),
    AllallowedDBPosR0(AllallowedDBPosR0),
    Lines(Lines),
    Columns(Columns),
    Nullpunkt(Nullpunkt),
    DBanzahl(DBanzahl),
    Distances(Distances) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      for(std::size_t j = 0; j < i; j++){
        double Dx = abs(AllDataBotsPosX[j]-AllDataBotsPosX[i]);
        double Dy = abs(AllDataBotsPosY[j]-AllDataBotsPosY[i]);
        double D1 = Lines-Dx+1;
        double D2 = Columns-Dy+1;
        //cout << D2;
        Dx = min(Dx, D1) + Nullpunkt[0]-1;
        Dy = min(Dy, D2) + Nullpunkt[1]-1;
        Distances(i,j) = AllallowedDBPosR0(Dx,Dy);
        Distances(j,i) = AllallowedDBPosR0(Dx,Dy);
        //Distances(i,j) = AllallowedDBPosR0(Dx[j],Dy[j]);
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix rcppPar_DistanceToroid(NumericVector AllDataBotsPosX,
                                     NumericVector AllDataBotsPosY,
                                     NumericMatrix AllallowedDBPosR0,
                                     double Lines,
                                     double Columns,
                                     NumericVector Nullpunkt){
  // Fuer aufrufe in PswarmCPP (einmalig) und pswarm() damit es sich nicht mit
  // der internen CPP fkt von PswarmCPP ueberschneidet
  int DBanzahl = AllDataBotsPosX.length();
  NumericMatrix Distances(DBanzahl,DBanzahl);
  //std::size_t end = 10;//AllDataBotsPosX.length();
  // create the worker
  PolarDistance polarDistance(AllDataBotsPosX,
                              AllDataBotsPosY,
                              AllallowedDBPosR0,
                              Lines,
                              Columns,
                              Nullpunkt,
                              DBanzahl,
                              Distances);
  // call it with parallelFor
  parallelFor(0, DBanzahl, polarDistance);
  return Distances;
}

