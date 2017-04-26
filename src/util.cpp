// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include "RxODE_types.h"

using namespace Rcpp;
using namespace R;
using namespace arma;

//' Echo cout to console for a number
//'
//' @param number number to output
//'
//'  @return nothing.
//'
//' @export
// [[Rcpp::export]]
void rxCoutEcho(NumericVector number){
  Rcout << number[0];
}
