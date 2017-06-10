// [[Rcpp::depends(lamW)]]

// Taken from
// https://github.com/cran/LambertW/blob/master/src/W_Cpp.cpp and
// modified to allow accessing in C, so that lambertW can be accessed in RxODE as a function
//
#include <Rcpp.h>
#include <lamW.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector W_Cpp(const NumericVector& z, int branch) {
  NumericVector Wz(z.length());
  if (branch == 0) {
    Wz = lamW::lambertW0_C(z);
  } else if (branch == -1) {
    Wz = lamW::lambertWm1_C(z);
  } else {
    stop("Only principal (0) and non-principal branch (-1) are implemented.");
  }
  return Wz;
}

extern "C" double W(double z, int branch);
double W(double z, int branch){
  NumericVector Wz;
  NumericVector zz(1);
  Wz = W_Cpp(zz, branch);
  double ret = Wz[0];
  return ret;
}
