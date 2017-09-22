// [[Rcpp::depends(RcppArmadillo)]]
#include <stdarg.h>
#include <RcppArmadillo.h>
#include <R.h>
using namespace Rcpp;
using namespace R;
using namespace arma;

//' Invert matrix using Rcpp Armadilo.  
//'
//' @param matrix matrix to be inverted.
//' @return inverse or pseudo inverse of matrix.
//' @export
// [[Rcpp::export]]
NumericVector rxInv(SEXP matrix){
  mat smatrix= as<mat>(matrix);
  mat imat;
  try{
    imat = inv(smatrix);
  } catch(...){
    Rprintf("Warning: matrix seems singular; Using pseudo-inverse\n");
    imat = pinv(smatrix);
  }
  NumericVector ret;
  ret = wrap(imat);
  return(ret);
}

// [[Rcpp::export]]
void RxODE_finalize_focei_omega(SEXP rho){
  Environment e = as<Environment>(rho);
  List dOmega = as<List>(e["dOmega"]);
  mat omegaInv = as<mat>(e["omegaInv"]);
  mat c;
  vec diag;
  int ntheta = dOmega.length();
  int neta = omegaInv.n_rows;
  mat cEta = zeros(neta,1);
  NumericVector trInv(ntheta);
  List prod1(ntheta);
  List prod2(ntheta);
  int i,j;
  for (i = 0; i < ntheta; i++){
    c = omegaInv * as<mat>(dOmega[i]);
    diag = c.diag();
    trInv[i] = 0.5*sum(diag);
    c = c * omegaInv;
    prod1[i] = c;
    List prodI(neta);
    for (j = 0; j < neta; j++){
      cEta(j,0) = 1;
      prodI[j] = c*cEta;
      cEta(j,0) = 0;
    }
    prod2[i] = prodI;
  }
  e["tr.omegaInv.dOmega.0.5"] = trInv;
  e["omegaInv.dOmega.omegaInv"] = prod1;
  e["omegaInv.dOmega.omegaInv.dEta"] = prod2;
}
// [[Rcpp::export]]
NumericVector RxODE_finalize_log_det_OMGAinv_5(SEXP rho){
  // log(det(omegaInv^1/2)) = 1/2*log(det(omegaInv))
  Environment e = as<Environment>(rho);
  NumericVector reset;
  mat c;
  try{
    c = chol(as<mat>(e["omegaInv"]));
  } catch(...){
    e["reset"] = 0;
    c = as<mat>(e["omegaInv"]);
    Function nearPD = as<Function>(e["nearPD"]);
    c = as<mat>(nearPD(c, rho));
    reset = as<NumericVector>(e["reset"]);
    if (reset[0] != 1){
      Rprintf("Warning: The Omega^-1 is non-positive definite, correcting with nearPD\n");
    }
    c = chol(-as<mat>(e["H"]));
  }
  vec diag = c.diag();
  vec ldiag = log(diag);
  NumericVector ret = as<NumericVector>(wrap(sum(ldiag)));
  e["log.det.OMGAinv.5"] = ret;
  return ret;
}

