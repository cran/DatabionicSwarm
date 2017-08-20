// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Delta3DWeightsC
arma::cube Delta3DWeightsC(Rcpp::NumericVector vx, Rcpp::NumericVector Datasample);
RcppExport SEXP _DatabionicSwarm_Delta3DWeightsC(SEXP vxSEXP, SEXP DatasampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vx(vxSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Datasample(DatasampleSEXP);
    rcpp_result_gen = Rcpp::wrap(Delta3DWeightsC(vx, Datasample));
    return rcpp_result_gen;
END_RCPP
}
// DijkstraSSSP
NumericVector DijkstraSSSP(NumericMatrix Adj, NumericMatrix Costs, int source);
RcppExport SEXP _DatabionicSwarm_DijkstraSSSP(SEXP AdjSEXP, SEXP CostsSEXP, SEXP sourceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Adj(AdjSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Costs(CostsSEXP);
    Rcpp::traits::input_parameter< int >::type source(sourceSEXP);
    rcpp_result_gen = Rcpp::wrap(DijkstraSSSP(Adj, Costs, source));
    return rcpp_result_gen;
END_RCPP
}
// findPossiblePositionsCsingle
ComplexVector findPossiblePositionsCsingle(NumericMatrix RadiusPositionsschablone, double jumplength, double alpha, double Lines);
RcppExport SEXP _DatabionicSwarm_findPossiblePositionsCsingle(SEXP RadiusPositionsschabloneSEXP, SEXP jumplengthSEXP, SEXP alphaSEXP, SEXP LinesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type RadiusPositionsschablone(RadiusPositionsschabloneSEXP);
    Rcpp::traits::input_parameter< double >::type jumplength(jumplengthSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type Lines(LinesSEXP);
    rcpp_result_gen = Rcpp::wrap(findPossiblePositionsCsingle(RadiusPositionsschablone, jumplength, alpha, Lines));
    return rcpp_result_gen;
END_RCPP
}
// PswarmCurrentRadiusC2botsPositive
List PswarmCurrentRadiusC2botsPositive(ComplexVector AllDataBotsPosOld, double Radius, NumericMatrix DataDists, ComplexVector IndPossibleDBPosR, NumericMatrix RadiusPositionsschablone, NumericVector pp, NumericVector Nullpunkt, double Lines, double Columns, double nBots, int limit, int steigungsverlaufind, double StressConstAditiv, bool debug);
RcppExport SEXP _DatabionicSwarm_PswarmCurrentRadiusC2botsPositive(SEXP AllDataBotsPosOldSEXP, SEXP RadiusSEXP, SEXP DataDistsSEXP, SEXP IndPossibleDBPosRSEXP, SEXP RadiusPositionsschabloneSEXP, SEXP ppSEXP, SEXP NullpunktSEXP, SEXP LinesSEXP, SEXP ColumnsSEXP, SEXP nBotsSEXP, SEXP limitSEXP, SEXP steigungsverlaufindSEXP, SEXP StressConstAditivSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ComplexVector >::type AllDataBotsPosOld(AllDataBotsPosOldSEXP);
    Rcpp::traits::input_parameter< double >::type Radius(RadiusSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DataDists(DataDistsSEXP);
    Rcpp::traits::input_parameter< ComplexVector >::type IndPossibleDBPosR(IndPossibleDBPosRSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type RadiusPositionsschablone(RadiusPositionsschabloneSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nullpunkt(NullpunktSEXP);
    Rcpp::traits::input_parameter< double >::type Lines(LinesSEXP);
    Rcpp::traits::input_parameter< double >::type Columns(ColumnsSEXP);
    Rcpp::traits::input_parameter< double >::type nBots(nBotsSEXP);
    Rcpp::traits::input_parameter< int >::type limit(limitSEXP);
    Rcpp::traits::input_parameter< int >::type steigungsverlaufind(steigungsverlaufindSEXP);
    Rcpp::traits::input_parameter< double >::type StressConstAditiv(StressConstAditivSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(PswarmCurrentRadiusC2botsPositive(AllDataBotsPosOld, Radius, DataDists, IndPossibleDBPosR, RadiusPositionsschablone, pp, Nullpunkt, Lines, Columns, nBots, limit, steigungsverlaufind, StressConstAditiv, debug));
    return rcpp_result_gen;
END_RCPP
}
// rDistanceToroidCsingle
NumericMatrix rDistanceToroidCsingle(NumericVector AllDataBotsPosX, NumericVector AllDataBotsPosY, NumericMatrix AllallowedDBPosR0, double Lines, double Columns, NumericVector Nullpunkt);
RcppExport SEXP _DatabionicSwarm_rDistanceToroidCsingle(SEXP AllDataBotsPosXSEXP, SEXP AllDataBotsPosYSEXP, SEXP AllallowedDBPosR0SEXP, SEXP LinesSEXP, SEXP ColumnsSEXP, SEXP NullpunktSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type AllDataBotsPosX(AllDataBotsPosXSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type AllDataBotsPosY(AllDataBotsPosYSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type AllallowedDBPosR0(AllallowedDBPosR0SEXP);
    Rcpp::traits::input_parameter< double >::type Lines(LinesSEXP);
    Rcpp::traits::input_parameter< double >::type Columns(ColumnsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Nullpunkt(NullpunktSEXP);
    rcpp_result_gen = Rcpp::wrap(rDistanceToroidCsingle(AllDataBotsPosX, AllDataBotsPosY, AllallowedDBPosR0, Lines, Columns, Nullpunkt));
    return rcpp_result_gen;
END_RCPP
}
// trainstepC
arma::cube trainstepC(Rcpp::NumericVector vx, Rcpp::NumericVector vy, Rcpp::NumericMatrix DataSampled, Rcpp::NumericMatrix BMUsampled, double Lines, double Columns, double Radius, bool toroid);
RcppExport SEXP _DatabionicSwarm_trainstepC(SEXP vxSEXP, SEXP vySEXP, SEXP DataSampledSEXP, SEXP BMUsampledSEXP, SEXP LinesSEXP, SEXP ColumnsSEXP, SEXP RadiusSEXP, SEXP toroidSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vx(vxSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vy(vySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type DataSampled(DataSampledSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type BMUsampled(BMUsampledSEXP);
    Rcpp::traits::input_parameter< double >::type Lines(LinesSEXP);
    Rcpp::traits::input_parameter< double >::type Columns(ColumnsSEXP);
    Rcpp::traits::input_parameter< double >::type Radius(RadiusSEXP);
    Rcpp::traits::input_parameter< bool >::type toroid(toroidSEXP);
    rcpp_result_gen = Rcpp::wrap(trainstepC(vx, vy, DataSampled, BMUsampled, Lines, Columns, Radius, toroid));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DatabionicSwarm_Delta3DWeightsC", (DL_FUNC) &_DatabionicSwarm_Delta3DWeightsC, 2},
    {"_DatabionicSwarm_DijkstraSSSP", (DL_FUNC) &_DatabionicSwarm_DijkstraSSSP, 3},
    {"_DatabionicSwarm_findPossiblePositionsCsingle", (DL_FUNC) &_DatabionicSwarm_findPossiblePositionsCsingle, 4},
    {"_DatabionicSwarm_PswarmCurrentRadiusC2botsPositive", (DL_FUNC) &_DatabionicSwarm_PswarmCurrentRadiusC2botsPositive, 14},
    {"_DatabionicSwarm_rDistanceToroidCsingle", (DL_FUNC) &_DatabionicSwarm_rDistanceToroidCsingle, 6},
    {"_DatabionicSwarm_trainstepC", (DL_FUNC) &_DatabionicSwarm_trainstepC, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_DatabionicSwarm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
