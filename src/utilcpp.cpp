// [[Rcpp::interfaces(r, cpp)]]
//#undef NDEBUG
#include <RcppArmadillo.h>
#include "../inst/include/RxODE.h"
#include <R.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif
using namespace Rcpp;
using namespace arma;

extern "C" SEXP _RxODE_isNullZero(SEXP); 

//[[Rcpp::export]]
LogicalVector isNullZero(RObject obj) {
  if (Rf_isNull(obj)) {
    return true;
  }
  int t = TYPEOF(obj);
  if (t == INTSXP || t == REALSXP) {
    // This tests for thetaMat
    if (obj.hasAttribute("dim")) {
      mat cur = as<arma::mat>(obj);
      if (cur.is_zero()) {
	return true;
      }
    }
  }
  if (t == VECSXP) {
    List cur = as<List>(obj);
    bool allZero  = true;
    for (int i = cur.size(); i--;) {
      if (!_RxODE_isNullZero(cur[i])) {
	return false;
      }
    }
    return true;
  }
  return false;
}
