// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
#include <RcppParallel.h>
#include <tuple>
#include "sampleC.h"

using namespace RcppParallel;
using namespace Rcpp;
using namespace sugar;
using namespace std;

// // [[Rcpp::depends(RcppArmadillo)]]
//NumericVector rcppPar_sampleC(NumericVector x, double len) {
//  bool replace=0;
//  return(RcppArmadillo::sample(x,len,replace));
//}

// // [[Rcpp::depends(RcppArmadillo)]]
//arma::vec sampleC2(NumericVector x, double len) {
//  bool replace=0;
//  return(RcppArmadillo::sample(x,len,replace));
//} 


// [[Rcpp::depends(RcppArmadillo)]]
Rcpp::NumericVector rcppPar_lmC2(Rcpp::NumericVector x, Rcpp::NumericVector yr){  // doku see lm() in R 
  Rcpp::NumericMatrix Xr(yr.length(),2);
  for(int i=0;i<yr.length();i++){
    Xr(i,0)=1;
    Xr(i,1)=x(i);
    }
  int n = Xr.nrow(), k = Xr.ncol();
  arma::mat X(Xr.begin(), n, k, false);       // reuses memory and avoids extra copy
  arma::colvec y(yr.begin(), yr.size(), false);
  arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
  Rcpp::NumericVector coefs(coef.begin(),coef.end());
  return(coefs);
}

// [[Rcpp::depends(RcppParallel)]]
int vecmaxInd2(std::vector<double> MyRVector) {
  int End    = MyRVector.end() - MyRVector.begin();
  int IdxMax = 0;
  for(int i = 1; i < End; i++){
    if(MyRVector[IdxMax] <= MyRVector[i]){
      IdxMax = i;
    }
  }
  return IdxMax;
}

// [[Rcpp::depends(RcppParallel)]]
struct rcppPar_DataBotsPosNeu : public Worker{                            // Worker for parallelization
  const RVector<double> CoordsMoveRe, CoordsMoveIm, ChosenForJump;
  int Lines, Columns, NumBK1, NumBK2, NumAllDB, NumChoDB;
  RVector<double> DataBotsPos;
  
  rcppPar_DataBotsPosNeu(const NumericVector CoordsMoveRe, const NumericVector CoordsMoveIm,
                         const NumericVector ChosenForJump,
                         int Lines, int Columns, int NumBK1, int NumBK2, int NumAllDB, int NumChoDB,
                         const NumericVector DataBotsPos):
    CoordsMoveRe(CoordsMoveRe),
    CoordsMoveIm(CoordsMoveIm),
    ChosenForJump(ChosenForJump),
    Lines(Lines),
    Columns(Columns),
    NumBK1(NumBK1),
    NumBK2(NumBK2),
    NumAllDB(NumAllDB),
    NumChoDB(NumChoDB),
    DataBotsPos(DataBotsPos) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      
      int Counter = 0;
      int Check   = i;
      
      while(Check > NumChoDB){
        Check   = Check - NumChoDB;
        Counter = Counter + 1;
      }
      
      int db       = ChosenForJump[i];
      
      int PosOldRe = DataBotsPos[db + NumBK2];
      int PosOldIm = DataBotsPos[db + NumBK2 + NumBK1];
      
      int PosNewRe = PosOldRe + CoordsMoveRe[i];
      int PosNewIm = PosOldIm + CoordsMoveIm[i];
      
      if(PosNewRe > Lines){
        PosNewRe = PosNewRe - Lines;
      }
      if(PosNewRe < 0){
        PosNewRe = PosNewRe + Lines;
      }
      
      if(PosNewIm > Columns){
        PosNewIm = PosNewIm - Columns;
      }
      if(PosNewIm < 0){
        PosNewIm = PosNewIm + Columns;
      }
      
      int Res = 0;
      
      for(int j = 0; j < NumAllDB; j++){
        if((PosNewRe == DataBotsPos[j + NumBK2]) && (PosNewIm == DataBotsPos[j + NumBK2 + NumBK1])){
          Res = 1;
        }
      }
      
      if(Res == 0){
        DataBotsPos[db + Counter * NumAllDB]             = PosNewRe;
        DataBotsPos[db + Counter * NumAllDB + NumBK1]    = PosNewIm;
      }else{
        DataBotsPos[db + Counter * NumAllDB]             = PosOldRe;
        DataBotsPos[db + Counter * NumAllDB + NumBK1]    = PosOldIm;
      }
    }
  }
};

// [[Rcpp::depends(RcppParallel)]]
NumericVector NewPositions(NumericVector DataBotsPos, NumericVector CoordsMoveRe, NumericVector CoordsMoveIm, NumericVector ChosenForJump,
                           int Lines, int Columns, int NumBK1, int NumBK2, int NumAllDB, int NumChoDB, int NumBK4){
  rcppPar_DataBotsPosNeu dbsPosNew(CoordsMoveRe, CoordsMoveIm, ChosenForJump,
                                   Lines, Columns, NumBK1, NumBK2, NumAllDB, NumChoDB,
                                   DataBotsPos);
  parallelFor(0, NumBK4, dbsPosNew);
  return DataBotsPos;
}

// [[Rcpp::depends(RcppParallel)]]
struct GetHappiness : public Worker{
  const RVector<double> DataBotsPos, DataDists, AllallowedDBPosR0;
  int Lines, Columns, Origin1, Origin2, Radius, NumAllDB, NumFreeShape1, NumBK1;
  double Happiness, Eps;
  RVector<double> LocalHappiness;
  
  GetHappiness(const NumericVector DataBotsPos, const NumericVector DataDists, const NumericVector AllallowedDBPosR0,
             int Lines, int Columns, int Origin1, int Origin2, int Radius, int NumAllDB, int NumFreeShape1, int NumBK1,
             double Happiness, double Eps,
             const NumericVector LocalHappiness):
    DataBotsPos(DataBotsPos), DataDists(DataDists), AllallowedDBPosR0(AllallowedDBPosR0),
    Lines(Lines), Columns(Columns), Origin1(Origin1), Origin2(Origin2), Radius(Radius), NumAllDB(NumAllDB),
    NumFreeShape1(NumFreeShape1), NumBK1(NumBK1), Happiness(Happiness), Eps(Eps), LocalHappiness(LocalHappiness) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      
      int CurPosRe = DataBotsPos[i];
      int CurPosIm = DataBotsPos[i + NumBK1];
      int Counter = 0;
      int Check   = i;
      
      while(Check >= NumAllDB){
        Check   = Check - NumAllDB;
        Counter = Counter + 1;
      }
      
      double Sum1     = 0;
      double Sum2     = 0;
      double Happiness_i = 0;
      double Inv      = 1/(3.14159265*Radius*Radius);
      
      for(int j = 0; j < NumAllDB; j++){
        int Dx  = abs(CurPosRe - DataBotsPos[j + Counter * NumAllDB]);
        int Dy  = abs(CurPosIm - DataBotsPos[j + Counter * NumAllDB + NumBK1]);
        int D1  = Lines-Dx+1;
        int D2  = Columns-Dy+1;
        Dx      = min(Dx, D1) + Origin1 - 1;
        Dy      = min(Dy, D2) + Origin2 - 1;
        double tmpResD = AllallowedDBPosR0[Dx + Dy * NumFreeShape1];
        double tmpRes1 = max(0.0, (1 - ((tmpResD*tmpResD)*Inv)));
        Sum1    = Sum1 + tmpRes1;
        Sum2    = Sum2 + tmpRes1 * DataDists[Check + j * NumAllDB];
      }
      
      if(Sum1 <= Eps){
        Happiness_i = Happiness;
      }else{
        Happiness_i = Happiness - (Sum2/Sum1);
      }
      
      LocalHappiness[i] = Happiness_i;
    }
  }
};

// [[Rcpp::depends(RcppParallel)]]
NumericVector CalcHappiness(NumericVector DataBotsPos, NumericVector DataDists, NumericVector AllallowedDBPosR0, NumericVector LocalHappiness,
                         int Lines, int Columns, int Origin1, int Origin2, int Radius, int NumAllDB, int NumFreeShape1, int NumBK1, int NumBK2,
                         double Happiness, double Eps){
  GetHappiness dbsPosNew(DataBotsPos, DataDists, AllallowedDBPosR0,
                       Lines, Columns, Origin1, Origin2, Radius, NumAllDB, NumFreeShape1, NumBK1, Happiness, Eps, LocalHappiness);
  parallelFor(0, NumBK1, dbsPosNew);
  return LocalHappiness;
}

// [[Rcpp::depends(RcppParallel)]]
struct SelectBest : public Worker{
  const RVector<double> LocalHappiness;
  int NumAllDB, NumBK1, NumJumps;
  RVector<double> DataBotsPos;
  
  SelectBest(NumericVector LocalHappiness,
             int NumAllDB, int NumBK1, int NumJumps,
             NumericVector DataBotsPos):
    LocalHappiness(LocalHappiness), NumAllDB(NumAllDB), NumBK1(NumBK1), NumJumps(NumJumps), DataBotsPos(DataBotsPos) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      
      std::vector<double> ChooseHappiness(NumJumps, 0);
      
      for(int j = 0; j < NumJumps; j++){
        ChooseHappiness[j] = LocalHappiness[i + j * NumAllDB];
      }
      
      double HappinessOld = LocalHappiness[i + NumJumps * NumAllDB];
      
      int ind = vecmaxInd2(ChooseHappiness);
      double HappinessNew = ChooseHappiness[ind];
      
      if(HappinessNew > HappinessOld){
        for(int j = 0; j < (NumJumps+1); j++){
          if(j != ind){
            DataBotsPos[i + j * NumAllDB]          = DataBotsPos[i + ind * NumAllDB];
            DataBotsPos[i + j * NumAllDB + NumBK1] = DataBotsPos[i + ind * NumAllDB + NumBK1];
          }
        }
        DataBotsPos[2*NumBK1 + i] = HappinessNew;
      }else{
        for(int j = 0; j < NumJumps; j++){
          DataBotsPos[i + j * NumAllDB]          = DataBotsPos[i + NumJumps * NumAllDB];
          DataBotsPos[i + j * NumAllDB + NumBK1] = DataBotsPos[i + NumJumps * NumAllDB + NumBK1];
        }
        DataBotsPos[2*NumBK1 + i] = HappinessOld;
      }
    }
  }
};

// [[Rcpp::depends(RcppParallel)]]
NumericVector SelectBestHappiness(NumericVector DataBotsPos, NumericVector LocalHappiness, int NumAllDB, int NumBK1, int NumJumps){
  SelectBest dbsPosNew(LocalHappiness, NumAllDB, NumBK1, NumJumps, DataBotsPos);
  parallelFor(0, NumAllDB, dbsPosNew);
  return DataBotsPos;
}



// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
List PswarmRadiusParallel(NumericVector DataBotsPos, NumericVector DataDists, NumericVector AllallowedDBPosR0,
                          NumericVector IndPossibleDBPosRe, NumericVector IndPossibleDBPosIm,
                          int Lines,  int Columns, int Radius, int NumAllDB, int NumChoDB, int NumFreeShape1,
                          int NumJumps, int Origin1, int Origin2,
                          double Happiness, int MinIterations, int HappinessInclination, double Eps,
                          bool debug){
  // 
  // Author: Quirin Stier 2022
  
  bool Jumping = 1;
  int Iteration = 0;
  //double Eps = 0.0001;
  //double Eps = 0.000001;
  
  int NumBK1 = (NumJumps + 1) * NumAllDB;
  int NumBK2 = NumJumps * NumAllDB;
  //int NumBK3 = (NumJumps + 1) * NumChoDB;
  int NumBK4 = NumJumps * NumChoDB;
  
  NumericVector ChosenForJump(NumAllDB);
  NumericVector ChosenForJumpEnlarged(NumBK4);
  NumericVector CoordsMoveRe(NumBK4);
  NumericVector CoordsMoveIm(NumBK4);
  
  NumericVector tmpVec(NumChoDB);
  
  NumericVector LocalHappiness(NumBK1);
  
  NumericVector slopeVec(2);
  
  NumericVector KeyBot(NumAllDB); //dummy
  for(int i = 0; i < NumAllDB; i++){
    KeyBot(i) = i;
  }
  
  int KeyPossiblePositionLen = IndPossibleDBPosRe.length();
  NumericVector KeyPossiblePosition(KeyPossiblePositionLen);
  for(int i = 0; i < KeyPossiblePositionLen; i++){
    KeyPossiblePosition(i) = i;
  }
  
  NumericVector HappinessCourse;
  NumericVector KeySteigung(HappinessInclination);
  NumericVector HappinessTail(HappinessInclination);

  for(int i=0;i<HappinessInclination;i++){
    KeySteigung(i)=i;
  }
  
  while(Jumping){

    ChosenForJump = sampleC(KeyBot, NumChoDB);
    for(int i = 0; i < NumJumps; i++){
      for(int j = 0; j < NumChoDB; j++){
        ChosenForJumpEnlarged[j + i * NumChoDB] = ChosenForJump[j];
      }
    }
    
    for(int i = 0; i < NumJumps; i++){
      tmpVec = sampleC(KeyPossiblePosition, NumChoDB);
      for(int j = 0; j < NumChoDB; j++){
        CoordsMoveRe[j + i * NumChoDB] = IndPossibleDBPosRe[tmpVec[j]];
        CoordsMoveIm[j + i * NumChoDB] = IndPossibleDBPosIm[tmpVec[j]];
      }
    }
    
    DataBotsPos = NewPositions(DataBotsPos, CoordsMoveRe, CoordsMoveIm, ChosenForJumpEnlarged,
                               Lines, Columns, NumBK1, NumBK2, NumAllDB, NumChoDB, NumBK4);
    LocalHappiness = CalcHappiness(DataBotsPos, DataDists, AllallowedDBPosR0, LocalHappiness,
                             Lines, Columns, Origin1, Origin2, Radius, NumAllDB, NumFreeShape1, NumBK1, NumBK2, Happiness, Eps);
    DataBotsPos = SelectBestHappiness(DataBotsPos, LocalHappiness, NumAllDB, NumBK1, NumJumps);
    
    double HappinessSum = 0;
    for(int j = 0; j < NumAllDB; j++){
      HappinessSum = HappinessSum + DataBotsPos[2*NumBK1 + j];
    }
    
    HappinessCourse.push_back(HappinessSum);
    Iteration = Iteration+1;
    
    if(Iteration>MinIterations){
      int MinTail = min(HappinessInclination, Iteration);
      HappinessTail  = tail(HappinessCourse, MinTail);
      slopeVec=rcppPar_lmC2(KeySteigung,HappinessTail);
      //double tmpSum = 0;
      //for(int j = 1; j < MinTail; j++){
      //  tmpSum = tmpSum + abs(HappinessTail[j] - HappinessTail[j-1]);
      //}
      //slopeVec[1] = tmpSum;
      
      if(slopeVec[1] < Eps){
        Jumping=0;
      }
    }
    
    if(debug){
      if(Iteration%100==(MinIterations+1)){
        //CRAN limits cout,therefore this is deprecated
        //std::cout<<"Steigung:"<<slopeVec[1]<<"; Iteration:"<<Iteration<<"; Happiness:"<<SumHappiness<<";"<<std::endl;
      }
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("AllDataBotsPos")    = DataBotsPos,
                            Rcpp::Named("HappinessCourse")   = HappinessCourse,
                            Rcpp::Named("fokussiertlaufind") = Iteration
  );
}
