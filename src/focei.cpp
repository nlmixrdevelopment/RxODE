// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include "RxODE_types.h"

using namespace Rcpp;
using namespace R;
using namespace arma;

extern "C" SEXP RxODE_ode_solver_focei_eta (SEXP sexp_eta, SEXP sexp_rho);
extern "C" SEXP RxODE_ode_solver_focei_hessian(SEXP sexp_rho);

// [[Rcpp::export]]
NumericVector RxODE_focei_eta_lik(SEXP sexp_eta, SEXP sexp_rho){
  SEXP solve_env = RxODE_ode_solver_focei_eta(sexp_eta, sexp_rho);
  Environment e = as<Environment>(solve_env);
  mat omegaInv = as<mat>(e["omegaInv"]);
  vec eta = as<vec>(e["eta"]);
  vec aret = -(as<vec>(e["llik"])-0.5*(eta.t() * omegaInv * eta));
  e.assign("llik2",aret);
  NumericVector ret = as<NumericVector>(wrap(aret));
  return ret;
}

// [[Rcpp::export]]
NumericVector RxODE_focei_eta_lp(SEXP sexp_eta, SEXP sexp_rho){
  SEXP solve_env = RxODE_ode_solver_focei_eta(sexp_eta, sexp_rho);
  Environment e = as<Environment>(solve_env);
  mat omegaInv = as<mat>(e["omegaInv"]);
  vec eta = as<vec>(e["eta"]);
  vec aret = -(as<vec>(e["lp"])- omegaInv * eta);
  e.assign("ep2",aret);
  NumericVector ret = as<NumericVector>(wrap(aret));
  return ret;
}

// [[Rcpp::export]]
XPtr<rxFn2> RxODE_focei_eta(std::string fstr){
  if (fstr == "lik")
    return(XPtr<rxFn2>(new rxFn2(&RxODE_focei_eta_lik)));
  else if (fstr == "lp")
    return(XPtr<rxFn2>(new rxFn2(&RxODE_focei_eta_lp)));
  else 
    return XPtr<rxFn2>(R_NilValue); // runtime error as NULL no XPtr
}

// [[Rcpp::export]]
NumericVector RxODE_focei_finalize_llik(SEXP rho){
  RxODE_ode_solver_focei_hessian(rho);
  Environment e = as<Environment>(rho);
  // Calculate -1/2 log(det(-H)) by chol.
  mat c = chol(-as<mat>(e["H"]));
  vec diag = c.diag();
  vec ldiag = log(diag);
  NumericVector ret = -as<NumericVector>(e["llik2"]);
  // log(det(omegaInv^1/2)) = 1/2*log(det(omegaInv))
  ret += as<NumericVector>(e["log.det.OMGAinv.5"]);
  ret += -as<NumericVector>(wrap(sum(ldiag)));
  ret.attr("fitted") = as<NumericVector>(e["f"]);
  ret.attr("posthoc") = as<NumericVector>(e["eta"]);
  return ret;
}

// [[Rcpp::export]]
NumericVector RxODE_finalize_log_det_OMGAinv_5(SEXP rho){
  // log(det(omegaInv^1/2)) = 1/2*log(det(omegaInv))
  Environment e = as<Environment>(rho);
  mat c = chol(as<mat>(e["omegaInv"]));
  vec diag = c.diag();
  vec ldiag = log(diag);
  NumericVector ret = as<NumericVector>(wrap(sum(ldiag)));
  e["log.det.OMGAinv.5"] = ret;
  return ret;
}


// [[Rcpp::export]]
void RxODE_finalize_focei_omega(SEXP rho){
  Environment e = as<Environment>(rho);
  List dOmega = as<List>(e["dOmega"]);
  mat omegaInv = as<mat>(e["omegaInv"]);
  mat c;
  vec diag;
  int ntheta = dOmega.length();
  NumericVector trInv(ntheta);
  List prod1(ntheta);
  int i;
  for (i = 0; i < ntheta; i++){
    c = omegaInv * as<mat>(dOmega[i]);
    diag = c.diag();
    trInv[i] = 0.5*sum(diag);
    c = c * omegaInv;
    prod1[i] = c;
  }
  e["tr.omegaInv.dOmega.0.5"] = trInv;
  e["omegaInv.dOmega.omegaInv"] = prod1;
}

//' Calculate d(eta)/d(omega)
//'
//' @param eta the eta to caluclate the differential for.
//'
//' @param rho environment where omegaInv.dOmega.omegaInv and tr.omegaInv.dOmega.0.5
//' are calculated.  This is done with the rxSymInv function.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericVector rxDetaDomega(SEXP rho, SEXP eta_sexp){
  Environment e = as<Environment>(rho);
  List dOmega = as<List>(e["omegaInv.dOmega.omegaInv"]);
  NumericVector omegaInv = as<NumericVector>(e["tr.omegaInv.dOmega.0.5"]);
  mat eta = as<mat>(eta_sexp);
  mat c;
  vec ret;
  int ntheta = dOmega.length();
  NumericVector dEta(ntheta);
  int i;
  for (i = 0; i < ntheta; i++){
    c = 0.5*(eta.t() * as<mat>(dOmega[i]) * eta);
    ret = c.diag()-omegaInv[i];
    dEta[i] = sum(ret);
  }
  return(dEta);
}
