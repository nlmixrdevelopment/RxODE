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
  // I'm not sure why/if this is more efficient, it is equavialent to -1/2 log(det(-H)).
  mat c = chol(-as<mat>(e["H"]));
  vec diag = c.diag();
  vec ldiag = log(diag);
  // I don't understand the log.det.OMGAinv.5 piece, but we can add it in later?
  NumericVector ret = -as<NumericVector>(e["llik2"]);
  ret += as<NumericVector>(e["log.det.OMGAinv.5"]);
  ret += -as<NumericVector>(wrap(sum(ldiag)));
  ret.attr("fitted") = as<NumericVector>(e["f"]);
  ret.attr("posthoc") = as<NumericVector>(e["eta"]);
  return ret;
}
