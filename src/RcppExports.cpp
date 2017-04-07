// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "RxODE_types.h"
#include "../inst/include/RxODE_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rxInner
void rxInner(SEXP etanews, SEXP rho);
RcppExport SEXP RxODE_rxInner(SEXP etanewsSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type etanews(etanewsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    rxInner(etanews, rho);
    return R_NilValue;
END_RCPP
}
// rxHessian
void rxHessian(SEXP rho);
RcppExport SEXP RxODE_rxHessian(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    rxHessian(rho);
    return R_NilValue;
END_RCPP
}
// RxODE_focei_eta_lik
NumericVector RxODE_focei_eta_lik(SEXP sexp_eta, SEXP sexp_rho);
RcppExport SEXP RxODE_RxODE_focei_eta_lik(SEXP sexp_etaSEXP, SEXP sexp_rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sexp_eta(sexp_etaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexp_rho(sexp_rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(RxODE_focei_eta_lik(sexp_eta, sexp_rho));
    return rcpp_result_gen;
END_RCPP
}
// RxODE_focei_eta_lp
NumericVector RxODE_focei_eta_lp(SEXP sexp_eta, SEXP sexp_rho);
RcppExport SEXP RxODE_RxODE_focei_eta_lp(SEXP sexp_etaSEXP, SEXP sexp_rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sexp_eta(sexp_etaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sexp_rho(sexp_rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(RxODE_focei_eta_lp(sexp_eta, sexp_rho));
    return rcpp_result_gen;
END_RCPP
}
// RxODE_focei_eta
XPtr<rxFn2> RxODE_focei_eta(std::string fstr);
RcppExport SEXP RxODE_RxODE_focei_eta(SEXP fstrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fstr(fstrSEXP);
    rcpp_result_gen = Rcpp::wrap(RxODE_focei_eta(fstr));
    return rcpp_result_gen;
END_RCPP
}
// RxODE_focei_finalize_llik
NumericVector RxODE_focei_finalize_llik(SEXP rho);
RcppExport SEXP RxODE_RxODE_focei_finalize_llik(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(RxODE_focei_finalize_llik(rho));
    return rcpp_result_gen;
END_RCPP
}
// RxODE_finalize_log_det_OMGAinv_5
NumericVector RxODE_finalize_log_det_OMGAinv_5(SEXP rho);
RcppExport SEXP RxODE_RxODE_finalize_log_det_OMGAinv_5(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(RxODE_finalize_log_det_OMGAinv_5(rho));
    return rcpp_result_gen;
END_RCPP
}
// RxODE_finalize_focei_omega
void RxODE_finalize_focei_omega(SEXP rho);
RcppExport SEXP RxODE_RxODE_finalize_focei_omega(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    RxODE_finalize_focei_omega(rho);
    return R_NilValue;
END_RCPP
}
// rxDetaDomega
void rxDetaDomega(SEXP rho);
RcppExport SEXP RxODE_rxDetaDomega(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    rxDetaDomega(rho);
    return R_NilValue;
END_RCPP
}
// rxOuter_
void rxOuter_(SEXP rho);
RcppExport SEXP RxODE_rxOuter_(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    rxOuter_(rho);
    return R_NilValue;
END_RCPP
}
// rxDetaDtheta
void rxDetaDtheta(SEXP rho);
RcppExport SEXP RxODE_rxDetaDtheta(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    rxDetaDtheta(rho);
    return R_NilValue;
END_RCPP
}
// rxOuter
NumericVector rxOuter(SEXP rho);
RcppExport SEXP RxODE_rxOuter(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(rxOuter(rho));
    return rcpp_result_gen;
END_RCPP
}
// rxUpdateEtas
NumericVector rxUpdateEtas(SEXP DnDhS, SEXP DhS, SEXP initS, SEXP acceptNS);
RcppExport SEXP RxODE_rxUpdateEtas(SEXP DnDhSSEXP, SEXP DhSSEXP, SEXP initSSEXP, SEXP acceptNSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type DnDhS(DnDhSSEXP);
    Rcpp::traits::input_parameter< SEXP >::type DhS(DhSSEXP);
    Rcpp::traits::input_parameter< SEXP >::type initS(initSSEXP);
    Rcpp::traits::input_parameter< SEXP >::type acceptNS(acceptNSSEXP);
    rcpp_result_gen = Rcpp::wrap(rxUpdateEtas(DnDhS, DhS, initS, acceptNS));
    return rcpp_result_gen;
END_RCPP
}
// getMacroConstants
void getMacroConstants(SEXP rho);
RcppExport SEXP RxODE_getMacroConstants(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    getMacroConstants(rho);
    return R_NilValue;
END_RCPP
}
// getLinDerivs
void getLinDerivs(SEXP rho);
RcppExport SEXP RxODE_getLinDerivs(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    getLinDerivs(rho);
    return R_NilValue;
END_RCPP
}
