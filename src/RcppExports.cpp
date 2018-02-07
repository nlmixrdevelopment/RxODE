// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "RxODE_types.h"
#include "../inst/include/RxODE_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rxIs
bool rxIs(const RObject& obj, std::string cls);
RcppExport SEXP _RxODE_rxIs(SEXP objSEXP, SEXP clsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type obj(objSEXP);
    Rcpp::traits::input_parameter< std::string >::type cls(clsSEXP);
    rcpp_result_gen = Rcpp::wrap(rxIs(obj, cls));
    return rcpp_result_gen;
END_RCPP
}
// rxDataSetup
List rxDataSetup(const RObject& ro, const RObject& covNames, const RObject& sigma, const RObject& df, const int& ncoresRV, const bool& isChol, const int& nDisplayProgress, const StringVector& amountUnits, const StringVector& timeUnits);
RcppExport SEXP _RxODE_rxDataSetup(SEXP roSEXP, SEXP covNamesSEXP, SEXP sigmaSEXP, SEXP dfSEXP, SEXP ncoresRVSEXP, SEXP isCholSEXP, SEXP nDisplayProgressSEXP, SEXP amountUnitsSEXP, SEXP timeUnitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type ro(roSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type covNames(covNamesSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const int& >::type ncoresRV(ncoresRVSEXP);
    Rcpp::traits::input_parameter< const bool& >::type isChol(isCholSEXP);
    Rcpp::traits::input_parameter< const int& >::type nDisplayProgress(nDisplayProgressSEXP);
    Rcpp::traits::input_parameter< const StringVector& >::type amountUnits(amountUnitsSEXP);
    Rcpp::traits::input_parameter< const StringVector& >::type timeUnits(timeUnitsSEXP);
    rcpp_result_gen = Rcpp::wrap(rxDataSetup(ro, covNames, sigma, df, ncoresRV, isChol, nDisplayProgress, amountUnits, timeUnits));
    return rcpp_result_gen;
END_RCPP
}
// rxModelVars_
List rxModelVars_(const RObject& obj);
RcppExport SEXP _RxODE_rxModelVars_(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxModelVars_(obj));
    return rcpp_result_gen;
END_RCPP
}
// rxState
RObject rxState(const RObject& obj, RObject state);
RcppExport SEXP _RxODE_rxState(SEXP objSEXP, SEXP stateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type obj(objSEXP);
    Rcpp::traits::input_parameter< RObject >::type state(stateSEXP);
    rcpp_result_gen = Rcpp::wrap(rxState(obj, state));
    return rcpp_result_gen;
END_RCPP
}
// rxParams
CharacterVector rxParams(const RObject& obj);
RcppExport SEXP _RxODE_rxParams(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxParams(obj));
    return rcpp_result_gen;
END_RCPP
}
// rxDfdy
CharacterVector rxDfdy(const RObject& obj);
RcppExport SEXP _RxODE_rxDfdy(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxDfdy(obj));
    return rcpp_result_gen;
END_RCPP
}
// rxLhs
CharacterVector rxLhs(const RObject& obj);
RcppExport SEXP _RxODE_rxLhs(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxLhs(obj));
    return rcpp_result_gen;
END_RCPP
}
// rxInits
NumericVector rxInits(const RObject& obj, RObject vec, Nullable<CharacterVector> req, double defaultValue, bool noerror, bool noini);
RcppExport SEXP _RxODE_rxInits(SEXP objSEXP, SEXP vecSEXP, SEXP reqSEXP, SEXP defaultValueSEXP, SEXP noerrorSEXP, SEXP noiniSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type obj(objSEXP);
    Rcpp::traits::input_parameter< RObject >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< Nullable<CharacterVector> >::type req(reqSEXP);
    Rcpp::traits::input_parameter< double >::type defaultValue(defaultValueSEXP);
    Rcpp::traits::input_parameter< bool >::type noerror(noerrorSEXP);
    Rcpp::traits::input_parameter< bool >::type noini(noiniSEXP);
    rcpp_result_gen = Rcpp::wrap(rxInits(obj, vec, req, defaultValue, noerror, noini));
    return rcpp_result_gen;
END_RCPP
}
// rxSetupIni
NumericVector rxSetupIni(const RObject& obj, RObject inits);
RcppExport SEXP _RxODE_rxSetupIni(SEXP objSEXP, SEXP initsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type obj(objSEXP);
    Rcpp::traits::input_parameter< RObject >::type inits(initsSEXP);
    rcpp_result_gen = Rcpp::wrap(rxSetupIni(obj, inits));
    return rcpp_result_gen;
END_RCPP
}
// rxSetupScale
NumericVector rxSetupScale(const RObject& obj, RObject scale, Nullable<List> extraArgs);
RcppExport SEXP _RxODE_rxSetupScale(SEXP objSEXP, SEXP scaleSEXP, SEXP extraArgsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type obj(objSEXP);
    Rcpp::traits::input_parameter< RObject >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type extraArgs(extraArgsSEXP);
    rcpp_result_gen = Rcpp::wrap(rxSetupScale(obj, scale, extraArgs));
    return rcpp_result_gen;
END_RCPP
}
// rxDataParSetup
List rxDataParSetup(const RObject& object, const RObject& params, const RObject& events, const RObject& inits, const RObject& covs, const RObject& sigma, const RObject& sigmaDf, const int& nCoresRV, const bool& sigmaIsChol, const int& nDisplayProgress, const StringVector& amountUnits, const StringVector& timeUnits, const RObject& theta, const RObject& eta, const RObject& scale, const Nullable<List>& extraArgs);
RcppExport SEXP _RxODE_rxDataParSetup(SEXP objectSEXP, SEXP paramsSEXP, SEXP eventsSEXP, SEXP initsSEXP, SEXP covsSEXP, SEXP sigmaSEXP, SEXP sigmaDfSEXP, SEXP nCoresRVSEXP, SEXP sigmaIsCholSEXP, SEXP nDisplayProgressSEXP, SEXP amountUnitsSEXP, SEXP timeUnitsSEXP, SEXP thetaSEXP, SEXP etaSEXP, SEXP scaleSEXP, SEXP extraArgsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type object(objectSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type inits(initsSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type covs(covsSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type sigmaDf(sigmaDfSEXP);
    Rcpp::traits::input_parameter< const int& >::type nCoresRV(nCoresRVSEXP);
    Rcpp::traits::input_parameter< const bool& >::type sigmaIsChol(sigmaIsCholSEXP);
    Rcpp::traits::input_parameter< const int& >::type nDisplayProgress(nDisplayProgressSEXP);
    Rcpp::traits::input_parameter< const StringVector& >::type amountUnits(amountUnitsSEXP);
    Rcpp::traits::input_parameter< const StringVector& >::type timeUnits(timeUnitsSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const Nullable<List>& >::type extraArgs(extraArgsSEXP);
    rcpp_result_gen = Rcpp::wrap(rxDataParSetup(object, params, events, inits, covs, sigma, sigmaDf, nCoresRV, sigmaIsChol, nDisplayProgress, amountUnits, timeUnits, theta, eta, scale, extraArgs));
    return rcpp_result_gen;
END_RCPP
}
// rxSolveCsmall
SEXP rxSolveCsmall(const RObject& object, const Nullable<CharacterVector>& specParams, const Nullable<List>& extraArgs, const RObject& params, const RObject& events, const RObject& inits, const RObject& scale, const RObject& covs, const Nullable<List>& optsL);
RcppExport SEXP _RxODE_rxSolveCsmall(SEXP objectSEXP, SEXP specParamsSEXP, SEXP extraArgsSEXP, SEXP paramsSEXP, SEXP eventsSEXP, SEXP initsSEXP, SEXP scaleSEXP, SEXP covsSEXP, SEXP optsLSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const RObject& >::type object(objectSEXP);
    Rcpp::traits::input_parameter< const Nullable<CharacterVector>& >::type specParams(specParamsSEXP);
    Rcpp::traits::input_parameter< const Nullable<List>& >::type extraArgs(extraArgsSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type inits(initsSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const RObject& >::type covs(covsSEXP);
    Rcpp::traits::input_parameter< const Nullable<List>& >::type optsL(optsLSEXP);
    rcpp_result_gen = Rcpp::wrap(rxSolveCsmall(object, specParams, extraArgs, params, events, inits, scale, covs, optsL));
    return rcpp_result_gen;
END_RCPP
}
// rxSolveGet
RObject rxSolveGet(RObject obj, RObject arg, LogicalVector exact);
RcppExport SEXP _RxODE_rxSolveGet(SEXP objSEXP, SEXP argSEXP, SEXP exactSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type obj(objSEXP);
    Rcpp::traits::input_parameter< RObject >::type arg(argSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type exact(exactSEXP);
    rcpp_result_gen = Rcpp::wrap(rxSolveGet(obj, arg, exact));
    return rcpp_result_gen;
END_RCPP
}
// rxSolveUpdate
RObject rxSolveUpdate(RObject obj, RObject arg, RObject value);
RcppExport SEXP _RxODE_rxSolveUpdate(SEXP objSEXP, SEXP argSEXP, SEXP valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type obj(objSEXP);
    Rcpp::traits::input_parameter< RObject >::type arg(argSEXP);
    Rcpp::traits::input_parameter< RObject >::type value(valueSEXP);
    rcpp_result_gen = Rcpp::wrap(rxSolveUpdate(obj, arg, value));
    return rcpp_result_gen;
END_RCPP
}
// rxRmModelLib_
void rxRmModelLib_(std::string str);
RcppExport SEXP _RxODE_rxRmModelLib_(SEXP strSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type str(strSEXP);
    rxRmModelLib_(str);
    return R_NilValue;
END_RCPP
}
// rxGetRxODE
RObject rxGetRxODE(RObject obj);
RcppExport SEXP _RxODE_rxGetRxODE(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxGetRxODE(obj));
    return rcpp_result_gen;
END_RCPP
}
// rxIsCurrent
bool rxIsCurrent(RObject obj);
RcppExport SEXP _RxODE_rxIsCurrent(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxIsCurrent(obj));
    return rcpp_result_gen;
END_RCPP
}
// rxAssignPtr
void rxAssignPtr(SEXP object);
RcppExport SEXP _RxODE_rxAssignPtr(SEXP objectSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type object(objectSEXP);
    rxAssignPtr(object);
    return R_NilValue;
END_RCPP
}
// rxCores
IntegerVector rxCores();
RcppExport SEXP _RxODE_rxCores() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rxCores());
    return rcpp_result_gen;
END_RCPP
}
// rxDll
std::string rxDll(RObject obj);
RcppExport SEXP _RxODE_rxDll(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxDll(obj));
    return rcpp_result_gen;
END_RCPP
}
// rxC
std::string rxC(RObject obj);
RcppExport SEXP _RxODE_rxC(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxC(obj));
    return rcpp_result_gen;
END_RCPP
}
// rxIsLoaded
bool rxIsLoaded(RObject obj);
RcppExport SEXP _RxODE_rxIsLoaded(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxIsLoaded(obj));
    return rcpp_result_gen;
END_RCPP
}
// rxDynLoad
bool rxDynLoad(RObject obj);
RcppExport SEXP _RxODE_rxDynLoad(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxDynLoad(obj));
    return rcpp_result_gen;
END_RCPP
}
// rxDynUnload
bool rxDynUnload(RObject obj);
RcppExport SEXP _RxODE_rxDynUnload(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxDynUnload(obj));
    return rcpp_result_gen;
END_RCPP
}
// rxDelete
bool rxDelete(RObject obj);
RcppExport SEXP _RxODE_rxDelete(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(rxDelete(obj));
    return rcpp_result_gen;
END_RCPP
}
// cvPost
NumericMatrix cvPost(const double& nu, const NumericMatrix& Omega, const bool& omegaIsChol);
RcppExport SEXP _RxODE_cvPost(SEXP nuSEXP, SEXP OmegaSEXP, SEXP omegaIsCholSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type omegaIsChol(omegaIsCholSEXP);
    rcpp_result_gen = Rcpp::wrap(cvPost(nu, Omega, omegaIsChol));
    return rcpp_result_gen;
END_RCPP
}
// rwish
NumericMatrix rwish(double nu, NumericMatrix Omega);
RcppExport SEXP _RxODE_rwish(SEXP nuSEXP, SEXP OmegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Omega(OmegaSEXP);
    rcpp_result_gen = Rcpp::wrap(rwish(nu, Omega));
    return rcpp_result_gen;
END_RCPP
}
// riwish
NumericMatrix riwish(double nu, NumericMatrix Omega);
RcppExport SEXP _RxODE_riwish(SEXP nuSEXP, SEXP OmegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Omega(OmegaSEXP);
    rcpp_result_gen = Rcpp::wrap(riwish(nu, Omega));
    return rcpp_result_gen;
END_RCPP
}
// riwishDf
NumericMatrix riwishDf(double nu, NumericMatrix Omega);
RcppExport SEXP _RxODE_riwishDf(SEXP nuSEXP, SEXP OmegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Omega(OmegaSEXP);
    rcpp_result_gen = Rcpp::wrap(riwishDf(nu, Omega));
    return rcpp_result_gen;
END_RCPP
}
// rsiwish
NumericMatrix rsiwish(double nu, NumericMatrix Omega, NumericVector mu, NumericVector delta);
RcppExport SEXP _RxODE_rsiwish(SEXP nuSEXP, SEXP OmegaSEXP, SEXP muSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(rsiwish(nu, Omega, mu, delta));
    return rcpp_result_gen;
END_RCPP
}
// rxSimThetaOmega
List rxSimThetaOmega(const Nullable<NumericVector>& params, const Nullable<NumericMatrix>& omega, const Nullable<NumericMatrix>& omegaDf, const bool& omegaIsChol, int nSub, const Nullable<NumericMatrix>& thetaMat, const Nullable<NumericMatrix>& thetaDf, const bool& thetaIsChol, int nStud, const Nullable<NumericVector> sigma, int nCoresRV, bool simVariability, int nObs);
RcppExport SEXP _RxODE_rxSimThetaOmega(SEXP paramsSEXP, SEXP omegaSEXP, SEXP omegaDfSEXP, SEXP omegaIsCholSEXP, SEXP nSubSEXP, SEXP thetaMatSEXP, SEXP thetaDfSEXP, SEXP thetaIsCholSEXP, SEXP nStudSEXP, SEXP sigmaSEXP, SEXP nCoresRVSEXP, SEXP simVariabilitySEXP, SEXP nObsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Nullable<NumericVector>& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericMatrix>& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericMatrix>& >::type omegaDf(omegaDfSEXP);
    Rcpp::traits::input_parameter< const bool& >::type omegaIsChol(omegaIsCholSEXP);
    Rcpp::traits::input_parameter< int >::type nSub(nSubSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericMatrix>& >::type thetaMat(thetaMatSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericMatrix>& >::type thetaDf(thetaDfSEXP);
    Rcpp::traits::input_parameter< const bool& >::type thetaIsChol(thetaIsCholSEXP);
    Rcpp::traits::input_parameter< int >::type nStud(nStudSEXP);
    Rcpp::traits::input_parameter< const Nullable<NumericVector> >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nCoresRV(nCoresRVSEXP);
    Rcpp::traits::input_parameter< bool >::type simVariability(simVariabilitySEXP);
    Rcpp::traits::input_parameter< int >::type nObs(nObsSEXP);
    rcpp_result_gen = Rcpp::wrap(rxSimThetaOmega(params, omega, omegaDf, omegaIsChol, nSub, thetaMat, thetaDf, thetaIsChol, nStud, sigma, nCoresRV, simVariability, nObs));
    return rcpp_result_gen;
END_RCPP
}
// rxInv
NumericVector rxInv(SEXP matrix);
RcppExport SEXP _RxODE_rxInv(SEXP matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type matrix(matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(rxInv(matrix));
    return rcpp_result_gen;
END_RCPP
}
// rxToOmega
arma::mat rxToOmega(arma::mat cholMat);
RcppExport SEXP _RxODE_rxToOmega(SEXP cholMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type cholMat(cholMatSEXP);
    rcpp_result_gen = Rcpp::wrap(rxToOmega(cholMat));
    return rcpp_result_gen;
END_RCPP
}
// rxSymInvChol
RObject rxSymInvChol(RObject invObjOrMatrix, Nullable<NumericVector> theta, std::string type, int thetaNumber);
RcppExport SEXP _RxODE_rxSymInvChol(SEXP invObjOrMatrixSEXP, SEXP thetaSEXP, SEXP typeSEXP, SEXP thetaNumberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type invObjOrMatrix(invObjOrMatrixSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type thetaNumber(thetaNumberSEXP);
    rcpp_result_gen = Rcpp::wrap(rxSymInvChol(invObjOrMatrix, theta, type, thetaNumber));
    return rcpp_result_gen;
END_RCPP
}
// rxSymInvCholEnvCalculate
RObject rxSymInvCholEnvCalculate(List obj, std::string what, Nullable<NumericVector> theta);
RcppExport SEXP _RxODE_rxSymInvCholEnvCalculate(SEXP objSEXP, SEXP whatSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type obj(objSEXP);
    Rcpp::traits::input_parameter< std::string >::type what(whatSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(rxSymInvCholEnvCalculate(obj, what, theta));
    return rcpp_result_gen;
END_RCPP
}
// rxInvWishartVar
arma::mat rxInvWishartVar(arma::mat Omega, double nu);
RcppExport SEXP _RxODE_rxInvWishartVar(SEXP OmegaSEXP, SEXP nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    rcpp_result_gen = Rcpp::wrap(rxInvWishartVar(Omega, nu));
    return rcpp_result_gen;
END_RCPP
}
// RxODE_finalize_focei_omega
void RxODE_finalize_focei_omega(RObject rho);
RcppExport SEXP _RxODE_RxODE_finalize_focei_omega(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type rho(rhoSEXP);
    RxODE_finalize_focei_omega(rho);
    return R_NilValue;
END_RCPP
}
// RxODE_finalize_log_det_OMGAinv_5
NumericVector RxODE_finalize_log_det_OMGAinv_5(SEXP rho);
RcppExport SEXP _RxODE_RxODE_finalize_log_det_OMGAinv_5(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(RxODE_finalize_log_det_OMGAinv_5(rho));
    return rcpp_result_gen;
END_RCPP
}
// rxCoutEcho
void rxCoutEcho(NumericVector number);
RcppExport SEXP _RxODE_rxCoutEcho(SEXP numberSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type number(numberSEXP);
    rxCoutEcho(number);
    return R_NilValue;
END_RCPP
}
// removableDrive
bool removableDrive(std::string driveRoot);
RcppExport SEXP _RxODE_removableDrive(SEXP driveRootSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type driveRoot(driveRootSEXP);
    rcpp_result_gen = Rcpp::wrap(removableDrive(driveRoot));
    return rcpp_result_gen;
END_RCPP
}
