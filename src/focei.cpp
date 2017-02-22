// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include "RxODE_types.hpp"

using namespace Rcpp;
using namespace R;
using namespace arma;

extern "C" SEXP RxODE_ode_solver_focei_eta (SEXP sexp_eta, SEXP sexp_rho);

// [[Rcpp::export]]
NumericVector RxODE_focei_eta_lik(SEXP sexp_eta, SEXP sexp_rho){
  SEXP solve_env = RxODE_ode_solver_focei_eta(sexp_eta, sexp_rho);
  Environment e = as<Environment>(solve_env);
  mat omegaInv = as<mat>(e["omegaInv"]);
  vec eta = as<vec>(e["eta"]);
  vec aret = -(as<vec>(e["llik"])-0.5*(eta.t() * omegaInv * eta));
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
  NumericVector ret = as<NumericVector>(wrap(aret));
  return ret;
}

// [[Rcpp::export]]
XPtr<rxFn2> RxODE_focei_eta(std::string fstr) {
  if (fstr == "lik")
    return(XPtr<rxFn2>(new rxFn2(&RxODE_focei_eta_lik)));
  else if (fstr == "lp")
    return(XPtr<rxFn2>(new rxFn2(&RxODE_focei_eta_lp)));
  else 
    return XPtr<rxFn2>(R_NilValue); // runtime error as NULL no XPtr
}
