// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;
using namespace sugar;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec sampleC(NumericVector x, double len) {
  bool replace=0;
  return(RcppArmadillo::sample(x,len,replace));
} 
