#include <RcppArmadillo.h>
#include <algorithm>
#include "../inst/include/RxODE.h"
using namespace Rcpp;

// expm::expm method="PadeRBS" is faster than rexpokit::expm BUT
// expm::expm method="Higham08" is "better" and has a C interface.  It
// is slower, though.

//' Enhance the Matrix for expm
//' This addst the columns indicating infusions if needed.
void enhanceMatrix(const arma::mat& m0, const arma::vec& InfusionRate,
		   const arma::mat& yp,
		   arma::mat &mout, arma::mat &ypout){
  unsigned int nrow = m0.n_rows;
  if (nrow != m0.n_cols)
    stop("m0 needs to be a square matrix");
  if (yp.n_elem != nrow)
    stop("yp needs to be the same dimension as m0.");
  if (InfusionRate.n_elem != nrow)
    stop("InfusionRate needs to be the same dimension as m0.");
  unsigned int i, nInf=0;
  arma::vec ypExtra(nrow);
  arma::mat m0extra(nrow, nrow);
  for (i = 0; i < nrow; i++){
    if (InfusionRate[i] != 0.0){
      nInf++;
      m0extra.resize(nrow*nInf,0);
      m0extra[nrow*(nInf-1)+i]=1;
      ypExtra[i] = InfusionRate[i];
    }
  }
  if (nInf == 0){
    mout = m0;
    ypout=yp;
    return;
  }
  mout = join_cols(join_rows(m0, m0extra.cols(0,nInf)),
		   arma::mat(nInf, nInf+nrow, arma::fill::zeros));
  ypout= join_cols(yp, ypExtra.head(nInf));
}

//' Enhance the Matrix for expm
//' This addst the columns indicating infusions if needed.
//' @param m0 Initial matrix
//' @param InfusionRate is a vector of infusion rates
//' @param yp is the last known state concentrations
//'
//' This is mostly for testing
//' @noRd
//[[Rcpp::export]]
List rxExpmMat(const arma::mat& m0, const arma::vec& InfusionRate,
	       const arma::mat& yp){
  arma::mat mout;
  arma::vec ypout;
  enhanceMatrix(m0, InfusionRate, yp, mout, ypout);
  List ret(2);
  ret[0] = wrap(mout);
  ret[1] = wrap(ypout);
  return ret;
}

//' Armadillo interface to R package expm
//'
//' @param inMat is the in matrix
//' @param t is the time for the calculation
//'
//' @inheritParams expm::expm
//'
//' @return expm(t*X)
//' @noRd
//[[Rcpp::export]]
arma::mat rxExpm(const arma::mat& inMat, double t = 1,
		 std::string method="Higham08.b"){
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment expmNS = loadNamespace("expm");
  Function expm = expmNS["expm"];
  arma::mat out0 =t*inMat;
  return as<arma::mat>(expm(out0,_["method"]=method));
}
