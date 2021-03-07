#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "../inst/include/RxODE.h"

#include "tran.h"
#include "threefry.h"

SEXP _rxHasOpenMp();

SEXP _vecDF(SEXP cv, SEXP n_);
SEXP _RxODE_dropUnitsRxSolve(SEXP);
SEXP _RxODE_atolRtolFactor_(SEXP);
SEXP _RxODE_etRep_(SEXP, SEXP, SEXP, SEXP, SEXP,
		   SEXP, SEXP);
SEXP _RxODE_etSeq_(SEXP, SEXP, SEXP, SEXP, SEXP,
		   SEXP, SEXP, SEXP, SEXP, SEXP,
		   SEXP);
SEXP _RxODE_rxSolveSEXP(SEXP, SEXP, SEXP, SEXP, SEXP,
			SEXP, SEXP, SEXP);
SEXP _RxODE_rxStack(SEXP, SEXP);
SEXP _RxODE_etUpdate(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_et_(SEXP, SEXP);
SEXP _RxODE_etTrans(SEXP, SEXP, SEXP, SEXP, SEXP,
		    SEXP, SEXP, SEXP);
SEXP _RxODE_rxUpdateTrans_(SEXP, SEXP, SEXP);
double powerDi(double x, double lambda, int yj);
double powerD(double x, double lambda, int yj);
double powerDD(double x, double lambda, int yj);
double powerDDD(double x, double lambda, int yj);
double powerL(double x, double lambda, int yj);
double powerDL(double x, double lambda, int yj);

SEXP _rxProgress(SEXP num, SEXP core);
SEXP _rxTick();
SEXP _rxProgressStop(SEXP);
SEXP _rxProgressAbort(SEXP);
SEXP _RxODE_codeLoaded();

SEXP _RxODE_codegen(SEXP c_file, SEXP prefix, SEXP libname, SEXP pMd5, SEXP timeId, SEXP lastMv);
SEXP _RxODE_parseModel(SEXP type);
SEXP _RxODE_isLinCmt();
SEXP _RxODE_RcppExport_registerCCallable();
SEXP _RxODE_setRstudio(SEXP);
SEXP _RxODE_rxSolveFree();
SEXP _RxODE_linCmtEnv(SEXP rho);
SEXP _RxODE_rxInv(SEXP matrix);
SEXP _RxODE_removableDrive(SEXP letter);
SEXP _rxCholInv(SEXP dms, SEXP theta, SEXP tn);
SEXP _RxODE_rxSymInvCholEnvCalculate(SEXP, SEXP, SEXP);
SEXP _RxODE_rxSymInvChol(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxIs(SEXP,SEXP);
SEXP _RxODE_rxModelVars_(SEXP);
SEXP _RxODE_rxState(SEXP, SEXP);
SEXP _RxODE_rxParams_(SEXP);
SEXP _RxODE_rxDfdy(SEXP);
SEXP _RxODE_rxLhs(SEXP);
SEXP _RxODE_rxInits(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxSetupIni(SEXP, SEXP);
SEXP _RxODE_rxSetupScale(SEXP,SEXP,SEXP);
SEXP _RxODE_rxSolveGet(SEXP, SEXP, SEXP);
SEXP _RxODE_rxSolveUpdate(SEXP, SEXP, SEXP);
SEXP _RxODE_rxAssignPtr(SEXP);
SEXP _RxODE_rxAssignPtr(SEXP objectSEXP);
SEXP _RxODE_dynLoad(SEXP dllSEXP);
SEXP _RxODE_rxOptRep_(SEXP);
SEXP _RxODE_rxIndLin_(SEXP);
SEXP _rxParProgress(SEXP);
SEXP _RxODE_forderForceBase(SEXP);
SEXP _RxODE_rxRmvn_(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxCholperm(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxGradpsi(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxNleq(SEXP, SEXP, SEXP);
SEXP _RxODE_rxMvnrnd(SEXP, SEXP, SEXP, SEXP, SEXP,
		     SEXP, SEXP);
SEXP _RxODE_rxMvrandn_(SEXP, SEXP, SEXP, SEXP, SEXP,
		       SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _RxODE_rxSolveDollarNames(SEXP);
SEXP _RxODE_etDollarNames(SEXP);
SEXP _RxODE_rxExpandNesting(SEXP, SEXP, SEXP);
SEXP _RxODE_expandTheta_(SEXP, SEXP, SEXP, SEXP, SEXP,
			 SEXP);

SEXP _RxODE_rxRmvnSEXP(SEXP, SEXP, SEXP, SEXP, SEXP,
		       SEXP, SEXP, SEXP, SEXP, SEXP,
		       SEXP, SEXP);

SEXP _cbindOme(SEXP et_, SEXP mat_, SEXP n_);

SEXP _RxODE_nestingInfo_(SEXP omega, SEXP data);

SEXP _RxODE_isNullZero(SEXP in);

SEXP RxODE_get_mv();

SEXP _gammap(SEXP, SEXP);
SEXP _gammaq(SEXP, SEXP);
SEXP _lowergamma(SEXP, SEXP);
SEXP _uppergamma(SEXP, SEXP);
SEXP _gammapDer(SEXP, SEXP);

SEXP _gammapInv(SEXP, SEXP);
SEXP _gammaqInv(SEXP, SEXP);

SEXP _gammapInva(SEXP, SEXP);
SEXP _gammaqInva(SEXP, SEXP);

double expit(double, double, double);
double logit(double, double, double);
SEXP _expit(SEXP, SEXP, SEXP);
SEXP _logit(SEXP, SEXP, SEXP);
SEXP _linCmtParse(SEXP vars, SEXP inStr, SEXP verbose);
SEXP _RxODE_linCmtGen(SEXP linCmt, SEXP vars, SEXP linCmtSens, SEXP verbose);

static R_NativePrimitiveArgType RxODE_Sum_t[] = {
  REALSXP, INTSXP
};

extern int RxODE_current_fn_pointer_id();
extern double RxODE_sum(double *input, int len);
extern double RxODE_prod(double *input, int len);

extern void RxODE_assign_fn_pointers(SEXP mv);


// Need to change to remove global variables
extern void RxODE_ode_free();

// Changed for Parallel
extern void RxODE_ode_freeP(rx_solve *rx, unsigned int id);

extern void rxRmModelLib(const char* s);
extern SEXP rxGetModelLib(const char *s);

extern SEXP _RxODE_rxRmModelLib_(SEXP);
extern SEXP _RxODE_rxDll(SEXP);
extern SEXP _RxODE_rxIsLoaded(SEXP);
extern SEXP _RxODE_rxDynUnload(SEXP);
extern SEXP _RxODE_rxDynLoad(SEXP);
extern SEXP _RxODE_rxDelete(SEXP);
extern SEXP _RxODE_rxGetRxODE(SEXP);
extern SEXP _RxODE_rxC(SEXP);

extern SEXP _RxODE_rxIsCurrent(SEXP);

SEXP _RxODE_cvPost_(SEXP, SEXP, SEXP, SEXP, SEXP,
		    SEXP, SEXP);
SEXP _RxODE_expandPars_(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rinvchisq(SEXP, SEXP, SEXP);

SEXP _RxODE_getRxFn(SEXP);
SEXP _RxODE_setProgSupported(SEXP);
SEXP _RxODE_getProgSupported();
SEXP _RxODE_rxSetIni0(SEXP);
SEXP _RxODE_rxSetSilentErr(SEXP silentSEXP);
SEXP _RxODE_rxUnloadAll_();
SEXP _RxODE_rxLock(SEXP);
SEXP _RxODE_rxUnlock(SEXP);
SEXP _RxODE_rxAllowUnload(SEXP);

SEXP _RxODE_rxExpandGrid_(SEXP, SEXP, SEXP);
SEXP _RxODE_rxExpandSens_(SEXP, SEXP);
SEXP _RxODE_rxExpandSens2_(SEXP, SEXP, SEXP);
SEXP _RxODE_rxExpandFEta_(SEXP, SEXP, SEXP);
SEXP _RxODE_rxRepR0_(SEXP);
SEXP _RxODE_rLKJ1(SEXP, SEXP, SEXP);
SEXP _RxODE_rLKJcv1(SEXP, SEXP);
SEXP _RxODE_rLKJcvLsd1(SEXP, SEXP, SEXP);
SEXP _RxODE_rcvC1(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _RxODE_rxSeedEng(SEXP);
SEXP _RxODE_rxnorm_(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxnormV_(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxpois_(SEXP, SEXP, SEXP);
SEXP _RxODE_rxt__(SEXP, SEXP, SEXP);
SEXP _RxODE_rxunif_(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxweibull_(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxgeom_(SEXP, SEXP, SEXP);
SEXP _RxODE_rxbeta_(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxgamma_(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxf_(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxexp_(SEXP, SEXP, SEXP);
SEXP _RxODE_rxchisq_(SEXP, SEXP, SEXP);
SEXP _RxODE_rxcauchy_(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxbinom_(SEXP, SEXP, SEXP, SEXP);

SEXP _RxODE_rxRmvn0(SEXP, SEXP, SEXP, SEXP, SEXP,
		    SEXP, SEXP, SEXP, SEXP, SEXP,
		    SEXP);

SEXP _hasStan();

SEXP _RxODE_invWR1d(SEXP dSEXP, SEXP nuSEXP, SEXP omegaIsCholSEXP);

SEXP _RxODE_rxSimThetaOmega(SEXP paramsSEXP, SEXP omegaSEXP, SEXP omegaDfSEXP, SEXP omegaLowerSEXP, SEXP omegaUpperSEXP,
			    SEXP omegaIsCholSEXP, SEXP omegaSeparationSEXP, SEXP omegaXformSEXP, SEXP nSubSEXP, SEXP thetaMatSEXP,
			    SEXP thetaLowerSEXP, SEXP thetaUpperSEXP, SEXP thetaDfSEXP, SEXP thetaIsCholSEXP, SEXP nStudSEXP,
			    SEXP sigmaSEXP, SEXP sigmaLowerSEXP, SEXP sigmaUpperSEXP, SEXP sigmaDfSEXP, SEXP sigmaIsCholSEXP,
			    SEXP sigmaSeparationSEXP, SEXP sigmaXformSEXP, SEXP nCoresRVSEXP, SEXP nObsSEXP, SEXP dfSubSEXP,
			    SEXP dfObsSEXP, SEXP simSubjectsSEXP);

SEXP _RxODE_convertId_(SEXP);

SEXP _RxODE_rpp_(SEXP nS, SEXP lambdaS, SEXP gammaS, SEXP probS, SEXP t0S,
		 SEXP tmaxS, SEXP randomOrderS);

SEXP _RxODE_rxQs(SEXP);

SEXP _RxODE_rxQr(SEXP);

extern int rxIsCurrentC(SEXP obj);

rx_solve *getRxSolve_();

void avoid_openmp_hang_within_fork();
void initRxThreads();

// Remove these functions later...

void rxOptionsIni();
/* void rxOptionsIniFocei(); */

double linCmtA(rx_solve *rx, unsigned int id, double t, int linCmt,
	       int ncmt, int trans, double d_ka,
	       double p1, double v1,
	       double p2, double p3,
	       double p4, double p5,
	       double d_tlag, double d_tlag2, double d_F, double d_F2,
	       double d_rate, double d_dur, double d_rate2, double d_dur2);

double linCmtC(rx_solve *rx, unsigned int id, double t, int linCmt,
	       int ncmt, int trans, double d_ka,
	       double p1, double v1,
	       double p2, double p3,
	       double p4, double p5,
	       double d_tlag, double d_tlag2, double d_F, double d_F2,
	       double d_rate, double d_dur, double d_rate2, double d_dur2);

double linCmtB(rx_solve *rx, unsigned int id, double t, int linCmt,
	       int i_cmt, int trans, int val,
	       double dd_p1, double dd_v1,
	       double dd_p2, double dd_p3,
	       double dd_p4, double dd_p5,
	       double dd_ka,
	       double dd_tlag, double dd_tlag2,
	       double dd_F, double dd_F2,
	       double dd_rate, double dd_dur,
	       double dd_rate2, double dd_dur2);

void _update_par_ptr(double t, unsigned int id, rx_solve *rx, int idx);
double _getParCov(unsigned int id, rx_solve *rx, int parNo, int idx);

int par_progress(int c, int n, int d, int cores, clock_t t0, int stop);
void ind_solve(rx_solve *rx, unsigned int cid, t_dydt_liblsoda dydt_lls, 
	       t_dydt_lsoda_dum dydt_lsoda, t_jdum_lsoda jdum,
	       t_dydt c_dydt, t_update_inis u_inis, int jt);
void par_solve(rx_solve *rx);
int isRstudio();

const char *rxGetId(int id);


double phi(double q);
SEXP _phi(SEXP q);
SEXP _calcDerived(SEXP transSXP, SEXP ncmt, SEXP inp, SEXP dig);

double gamma_p(double, double z);
double gamma_q(double, double z);
double tgamma_lower(double a, double z);
double tgamma_upper(double a, double z);
double gamma_p_derivative(double a, double x);
double gamma_q_inv(double a, double q);
double gamma_p_inv(double a, double p);
double gamma_q_inva(double a, double q);
double gamma_p_inva(double a, double p);

int compareFactorVal(int val, const char *factor, const char *value);
SEXP _RxODE_rxSolve_(SEXP, SEXP, SEXP, SEXP, SEXP,
		     SEXP, SEXP, SEXP);

SEXP getRxThreads_R(SEXP verbose);
SEXP setRxthreads(SEXP threads, SEXP percent, SEXP throttle);

int getSilentErr();

int iniSubjectE(int solveid, int inLhs, rx_solving_options_ind *ind, rx_solving_options *op, rx_solve *rx,
		t_update_inis u_inis);

t_calc_lhs getRxLhs();
t_update_inis getUpdateInis();

void sortIds(rx_solve* rx, int ini);

void handleTlast(double *time, rx_solving_options_ind *ind);

SEXP _probit(SEXP xS, SEXP lowS, SEXP highS);
SEXP _probitInv(SEXP xS, SEXP lowS, SEXP highS);
SEXP _RxODE_rxrandnV(SEXP, SEXP);

void simeps(int id);
void simeta(int id);

double getTime(int idx, rx_solving_options_ind *ind);
void transIniNull();
void nullGlobals();
void R_init_RxODE(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"_rxProgress", (DL_FUNC) &_rxProgress, 2},
    {"_rxTick", (DL_FUNC) &_rxTick, 0},
    {"_rxProgressStop", (DL_FUNC) &_rxProgressStop, 1},
    {"_rxProgressAbort", (DL_FUNC) &_rxProgressAbort, 1},
    {"_RxODE_trans", (DL_FUNC) &_RxODE_trans, 7},
    {"_RxODE_codegen", (DL_FUNC) &_RxODE_codegen, 6},
    {"_RxODE_codeLoaded", (DL_FUNC) &_RxODE_codeLoaded, 0},
    {"_RxODE_parseModel", (DL_FUNC) &_RxODE_parseModel, 1},
    {"_RxODE_isLinCmt", (DL_FUNC) &_RxODE_isLinCmt, 0},
    {"RxODE_get_mv", (DL_FUNC) &RxODE_get_mv, 0},
    {"_RxODE_rxInv", (DL_FUNC) &_RxODE_rxInv, 1},
    {"_RxODE_removableDrive", (DL_FUNC) &_RxODE_removableDrive, 1},
    {"_rxCholInv", (DL_FUNC) &_rxCholInv, 3},
    {"_RxODE_rxSymInvCholEnvCalculate", (DL_FUNC) &_RxODE_rxSymInvCholEnvCalculate, 3},
    {"_RxODE_rxSymInvChol", (DL_FUNC) &_RxODE_rxSymInvChol, 4},
    {"_RxODE_rxIs", (DL_FUNC) &_RxODE_rxIs, 2},
    {"_RxODE_rxModelVars_", (DL_FUNC) &_RxODE_rxModelVars_, 1},
    {"_RxODE_rxState", (DL_FUNC) &_RxODE_rxState, 2},
    {"_RxODE_rxParams_", (DL_FUNC) &_RxODE_rxParams_, 1},
    {"_RxODE_rxDfdy", (DL_FUNC) &_RxODE_rxDfdy, 1},
    {"_RxODE_rxLhs", (DL_FUNC) &_RxODE_rxLhs, 1},
    {"_RxODE_rxInits", (DL_FUNC) &_RxODE_rxInits, 7},
    {"_RxODE_rxSetupIni", (DL_FUNC) &_RxODE_rxSetupIni, 2},
    {"_RxODE_rxSetupScale", (DL_FUNC) &_RxODE_rxSetupScale, 3},
    {"_RxODE_rxSolveGet", (DL_FUNC) &_RxODE_rxSolveGet, 3},
    {"_RxODE_rxSolveUpdate", (DL_FUNC) &_RxODE_rxSolveUpdate, 3},
    {"_RxODE_rxAssignPtr", (DL_FUNC) &_RxODE_rxAssignPtr, 1},
    {"_RxODE_rxRmModelLib_",(DL_FUNC) &_RxODE_rxRmModelLib_, 1},
    {"_RxODE_rxDll",(DL_FUNC) &_RxODE_rxDll, 1},
    {"_RxODE_rxC",(DL_FUNC) &_RxODE_rxC, 1},
    {"_RxODE_rxIsLoaded", (DL_FUNC) &_RxODE_rxIsLoaded, 1},
    {"_RxODE_rxDynUnload", (DL_FUNC) &_RxODE_rxDynUnload, 1},
    {"_RxODE_rxDynLoad", (DL_FUNC) &_RxODE_rxDynLoad, 1},
    {"_RxODE_rxDelete", (DL_FUNC) &_RxODE_rxDelete, 1},
    {"_RxODE_rxGetRxODE", (DL_FUNC) &_RxODE_rxGetRxODE, 1},
    {"_RxODE_rxIsCurrent", (DL_FUNC) &_RxODE_rxIsCurrent, 1},
    {"_RxODE_cvPost_", (DL_FUNC) &_RxODE_cvPost_, 7},
    {"_RxODE_rinvchisq", (DL_FUNC) &_RxODE_rinvchisq, 3},
    {"_RxODE_dynLoad", (DL_FUNC) &_RxODE_dynLoad, 1},
    {"_RxODE_rxSolveFree", (DL_FUNC) &_RxODE_rxSolveFree, 0},
    {"_RxODE_setRstudio", (DL_FUNC) &_RxODE_setRstudio, 1},
    {"_RxODE_RcppExport_registerCCallable", (DL_FUNC) &_RxODE_RcppExport_registerCCallable, 0},
    {"_RxODE_getRxFn", (DL_FUNC) &_RxODE_getRxFn, 1},
    {"_RxODE_setProgSupported", (DL_FUNC) &_RxODE_setProgSupported, 1},
    {"_RxODE_getProgSupported", (DL_FUNC) &_RxODE_getProgSupported, 0},
    {"_RxODE_rxUpdateTrans_", (DL_FUNC) &_RxODE_rxUpdateTrans_, 3},
    {"_RxODE_etTrans", (DL_FUNC) &_RxODE_etTrans, 8},
    {"_RxODE_et_", (DL_FUNC) &_RxODE_et_, 2},
    {"_RxODE_etUpdate", (DL_FUNC) &_RxODE_etUpdate, 4},
    {"_RxODE_rxStack", (DL_FUNC) &_RxODE_rxStack, 2},
    {"_RxODE_etSeq_", (DL_FUNC) &_RxODE_etSeq_, 11},
    {"_RxODE_etRep_", (DL_FUNC) &_RxODE_etRep_, 7},
    {"_RxODE_rxSolveSEXP", (DL_FUNC) &_RxODE_rxSolveSEXP, 8},
    {"_RxODE_dropUnitsRxSolve", (DL_FUNC) &_RxODE_dropUnitsRxSolve, 1},
    {"_RxODE_atolRtolFactor_", (DL_FUNC) &_RxODE_atolRtolFactor_, 1},
    {"_RxODE_rxExpandGrid_", (DL_FUNC) &_RxODE_rxExpandGrid_, 3},
    {"_RxODE_rxExpandSens_", (DL_FUNC) &_RxODE_rxExpandSens_, 2},
    {"_RxODE_rxExpandSens2_",(DL_FUNC) &_RxODE_rxExpandSens2_, 3},
    {"_RxODE_rxExpandFEta_", (DL_FUNC) &_RxODE_rxExpandFEta_, 3},
    {"_RxODE_rxRepR0_", (DL_FUNC) &_RxODE_rxRepR0_, 1},
    {"_RxODE_rxOptRep_", (DL_FUNC) &_RxODE_rxOptRep_, 1},
    {"_RxODE_rxSetIni0", (DL_FUNC) &_RxODE_rxSetIni0, 1},
    {"_RxODE_rxSetSilentErr", (DL_FUNC) &_RxODE_rxSetSilentErr, 1},
    {"_RxODE_rxIndLin_",(DL_FUNC) &_RxODE_rxIndLin_, 1},
    {"_RxODE_rxUnloadAll_", (DL_FUNC) &_RxODE_rxUnloadAll_, 0},
    {"_RxODE_rxLock", (DL_FUNC) &_RxODE_rxLock, 1},
    {"_RxODE_rxUnlock", (DL_FUNC) &_RxODE_rxUnlock, 1},
    {"_RxODE_rxAllowUnload", (DL_FUNC) &_RxODE_rxAllowUnload, 1},
    {"_rxParProgress", (DL_FUNC) &_rxParProgress, 1},
    {"_RxODE_forderForceBase", (DL_FUNC) &_RxODE_forderForceBase, 1},
    {"_RxODE_rLKJ1", (DL_FUNC) &_RxODE_rLKJ1, 3},
    {"_RxODE_rLKJcv1", (DL_FUNC) &_RxODE_rLKJcv1, 2},
    {"_RxODE_rLKJcvLsd1", (DL_FUNC) &_RxODE_rLKJcvLsd1, 3},
    {"_RxODE_rcvC1", (DL_FUNC) &_RxODE_rcvC1, 5},
    {"_RxODE_rxRmvn_", (DL_FUNC) &_RxODE_rxRmvn_, 5},
    {"_RxODE_rxCholperm", (DL_FUNC) &_RxODE_rxCholperm, 4},
    {"_RxODE_rxGradpsi", (DL_FUNC) &_RxODE_rxGradpsi, 4},
    {"_RxODE_rxNleq", (DL_FUNC) &_RxODE_rxNleq, 3},
    {"_RxODE_rxMvnrnd", (DL_FUNC) &_RxODE_rxMvnrnd, 7},
    {"_RxODE_rxMvrandn_", (DL_FUNC) &_RxODE_rxMvrandn_, 10},
    {"_RxODE_rxSeedEng", (DL_FUNC) &_RxODE_rxSeedEng, 1},
    {"_RxODE_rxnorm_", (DL_FUNC) &_RxODE_rxnorm_, 4},
    {"_RxODE_rxnormV_", (DL_FUNC) &_RxODE_rxnormV_, 4},
    {"_RxODE_rxpois_", (DL_FUNC) &_RxODE_rxpois_, 3},
    {"_RxODE_rxt__", (DL_FUNC) &_RxODE_rxt__, 3},
    {"_RxODE_rxunif_", (DL_FUNC) &_RxODE_rxunif_, 4},
    {"_RxODE_rxweibull_", (DL_FUNC) &_RxODE_rxweibull_, 4},
    {"_RxODE_rxgeom_", (DL_FUNC) &_RxODE_rxgeom_, 3},
    {"_RxODE_rxbeta_", (DL_FUNC) &_RxODE_rxbeta_, 4},
    {"_RxODE_rxgamma_", (DL_FUNC) &_RxODE_rxgamma_, 4},
    {"_RxODE_rxf_", (DL_FUNC) &_RxODE_rxf_, 4},
    {"_RxODE_rxexp_", (DL_FUNC) &_RxODE_rxexp_, 3},
    {"_RxODE_rxchisq_", (DL_FUNC) &_RxODE_rxchisq_, 3},
    {"_RxODE_rxcauchy_", (DL_FUNC) &_RxODE_rxcauchy_, 4},
    {"_RxODE_rxbinom_", (DL_FUNC) &_RxODE_rxbinom_, 4},
    {"_RxODE_rxSolveDollarNames", (DL_FUNC) _RxODE_rxSolveDollarNames, 1},
    {"_RxODE_etDollarNames", (DL_FUNC) _RxODE_etDollarNames, 1},
    {"_RxODE_rxExpandNesting", (DL_FUNC) _RxODE_rxExpandNesting, 3},
    {"_RxODE_rxRmvn0", (DL_FUNC) _RxODE_rxRmvn0, 11},
    {"_vecDF", (DL_FUNC) _vecDF, 2},
    {"_cbindOme", (DL_FUNC) _cbindOme, 3},
    {"_RxODE_rxRmvnSEXP", (DL_FUNC) _RxODE_rxRmvnSEXP, 12},
    {"_RxODE_expandPars_", (DL_FUNC) _RxODE_expandPars_, 4},
    {"_RxODE_convertId_", (DL_FUNC) _RxODE_convertId_, 1},
    {"_RxODE_nestingInfo_", (DL_FUNC) _RxODE_nestingInfo_, 2},
    {"_phi", (DL_FUNC) _phi, 1},
    {"_RxODE_expandTheta_", (DL_FUNC) _RxODE_expandTheta_, 6},
    {"_gammap", (DL_FUNC) _gammap, 2},
    {"_gammaq", (DL_FUNC) _gammaq, 2},
    {"_uppergamma", (DL_FUNC) _uppergamma, 2},
    {"_lowergamma", (DL_FUNC) _lowergamma, 2},
    {"_gammapDer", (DL_FUNC) _gammapDer, 2},
    {"_gammapInv", (DL_FUNC) _gammapInv, 2},
    {"_gammaqInv", (DL_FUNC) _gammaqInv, 2},
    {"_gammapInva", (DL_FUNC) _gammapInv, 2},
    {"_gammaqInva", (DL_FUNC) _gammaqInv, 2},
    {"_expit", (DL_FUNC) _expit, 3},
    {"_logit", (DL_FUNC) _logit, 3},
    {"_calcDerived", (DL_FUNC) _calcDerived, 4},
    {"_linCmtParse", (DL_FUNC) _linCmtParse, 3},
    {"_RxODE_linCmtGen", (DL_FUNC) _RxODE_linCmtGen, 4},
    {"_RxODE_rpp_", (DL_FUNC) _RxODE_rpp_, 7},
    {"_RxODE_rxSolve_", (DL_FUNC) _RxODE_rxSolve_, 8},
    {"getRxThreads_R", (DL_FUNC) getRxThreads_R, 1},
    {"setRxthreads", (DL_FUNC) setRxthreads, 3},
    {"_rxHasOpenMp", (DL_FUNC) _rxHasOpenMp, 0},
    {"_probit", (DL_FUNC) _probit, 3},
    {"_probitInv", (DL_FUNC) _probitInv, 3},
    {"_RxODE_rxrandnV", (DL_FUNC) _RxODE_rxrandnV, 2},
    {"_RxODE_rxQs", (DL_FUNC) _RxODE_rxQs, 1},
    {"_RxODE_rxQr", (DL_FUNC) _RxODE_rxQr, 1},
    {"_RxODE_isNullZero", (DL_FUNC) _RxODE_isNullZero, 1},
    {"_RxODE_invWR1d", (DL_FUNC) _RxODE_invWR1d, 3},
    {"_RxODE_rxSimThetaOmega", (DL_FUNC) _RxODE_rxSimThetaOmega, 27},
    {"_hasStan", (DL_FUNC) _hasStan, 0},
    {NULL, NULL, 0}
  };
  // C callable to assign environments.
  R_RegisterCCallable("RxODE","_RxODE_rxQr", (DL_FUNC) &_RxODE_rxQr);
  R_RegisterCCallable("RxODE", "simeps", (DL_FUNC) &simeps);
  R_RegisterCCallable("RxODE", "simeta", (DL_FUNC) &simeta);
  R_RegisterCCallable("RxODE", "getSilentErr", (DL_FUNC) &getSilentErr);
  R_RegisterCCallable("RxODE", "logit", (DL_FUNC) &logit);
  R_RegisterCCallable("RxODE", "expit", (DL_FUNC) &expit);
  R_RegisterCCallable("RxODE", "phi", (DL_FUNC) &phi);
  
  R_RegisterCCallable("RxODE", "rxnormV", (DL_FUNC) &rxnormV);
  R_RegisterCCallable("RxODE", "rinormV", (DL_FUNC) &rinormV);
  
  R_RegisterCCallable("RxODE", "rxgamma", (DL_FUNC) &rxgamma);
  R_RegisterCCallable("RxODE", "rxbeta", (DL_FUNC) &rxbeta);
  R_RegisterCCallable("RxODE", "rxbinom", (DL_FUNC) &rxbinom);
  R_RegisterCCallable("RxODE", "rxcauchy", (DL_FUNC) &rxcauchy);
  R_RegisterCCallable("RxODE", "rxchisq", (DL_FUNC) &rxchisq);
  R_RegisterCCallable("RxODE", "rxexp", (DL_FUNC) &rxexp);
  R_RegisterCCallable("RxODE", "rxf", (DL_FUNC) &rxf);
  R_RegisterCCallable("RxODE", "rxgeom", (DL_FUNC) &rxgeom);
  R_RegisterCCallable("RxODE", "rxnorm", (DL_FUNC) &rxnorm);
  R_RegisterCCallable("RxODE", "rxpois", (DL_FUNC) &rxpois);
  R_RegisterCCallable("RxODE", "rxt_", (DL_FUNC) &rxt_);
  R_RegisterCCallable("RxODE", "rxunif", (DL_FUNC) &rxunif);
  R_RegisterCCallable("RxODE", "rxweibull", (DL_FUNC) &rxweibull);

  R_RegisterCCallable("RxODE", "rigamma", (DL_FUNC) &rigamma);
  R_RegisterCCallable("RxODE", "ribeta", (DL_FUNC) &ribeta);
  R_RegisterCCallable("RxODE", "ribinom", (DL_FUNC) &ribinom);
  R_RegisterCCallable("RxODE", "ricauchy", (DL_FUNC) &ricauchy);
  R_RegisterCCallable("RxODE", "richisq", (DL_FUNC) &richisq);
  R_RegisterCCallable("RxODE", "riexp", (DL_FUNC) &riexp);
  R_RegisterCCallable("RxODE", "rif", (DL_FUNC) &rif);
  R_RegisterCCallable("RxODE", "rigeom", (DL_FUNC) &rigeom);
  R_RegisterCCallable("RxODE", "rinorm", (DL_FUNC) &rinorm);
  R_RegisterCCallable("RxODE", "ripois", (DL_FUNC) &ripois);
  R_RegisterCCallable("RxODE", "rit_", (DL_FUNC) &rit_);
  R_RegisterCCallable("RxODE", "riunif", (DL_FUNC) &riunif);
  R_RegisterCCallable("RxODE", "riweibull", (DL_FUNC) &riweibull);
  
  R_RegisterCCallable("RxODE", "powerDi", (DL_FUNC) &powerDi);
  R_RegisterCCallable("RxODE", "powerD", (DL_FUNC) &powerD);
  R_RegisterCCallable("RxODE", "powerDD", (DL_FUNC) &powerDD);
  R_RegisterCCallable("RxODE", "powerDDD", (DL_FUNC) &powerDDD);
  R_RegisterCCallable("RxODE", "powerL", (DL_FUNC) &powerL);
  R_RegisterCCallable("RxODE", "powerDL", (DL_FUNC) &powerDL);
  R_RegisterCCallable("RxODE", "par_progress", (DL_FUNC) &par_progress);
  R_RegisterCCallable("RxODE", "isRstudio", (DL_FUNC) &isRstudio);
  R_RegisterCCallable("RxODE", "ind_solve", (DL_FUNC) &ind_solve);
  R_RegisterCCallable("RxODE", "par_solve", (DL_FUNC) &par_solve);
  R_RegisterCCallable("RxODE", "linCmtA", (DL_FUNC) &linCmtA);
  R_RegisterCCallable("RxODE", "linCmtC", (DL_FUNC) &linCmtC);
  R_RegisterCCallable("RxODE", "linCmtB", (DL_FUNC) &linCmtB);
  R_RegisterCCallable("RxODE", "_update_par_ptr", (DL_FUNC) &_update_par_ptr);
  R_RegisterCCallable("RxODE", "_getParCov", (DL_FUNC) &_getParCov);
  R_RegisterCCallable("RxODE","rxRmModelLib", (DL_FUNC) &rxRmModelLib);
  R_RegisterCCallable("RxODE","rxGetModelLib", (DL_FUNC) &rxGetModelLib);
  
  R_RegisterCCallable("RxODE","RxODE_ode_free",           (DL_FUNC) &RxODE_ode_free);
  
  //Functions
  
  R_RegisterCCallable("RxODE","RxODE_sum",                (DL_FUNC) &RxODE_sum);
  R_RegisterCCallable("RxODE","RxODE_prod",               (DL_FUNC) &RxODE_prod);

  R_RegisterCCallable("RxODE","RxODE_assign_fn_pointers", (DL_FUNC) &RxODE_assign_fn_pointers);

  R_RegisterCCallable("RxODE","_RxODE_rxAssignPtr",       (DL_FUNC) _RxODE_rxAssignPtr);
  R_RegisterCCallable("RxODE", "rxIsCurrentC", (DL_FUNC) rxIsCurrentC);
  R_RegisterCCallable("RxODE","RxODE_current_fn_pointer_id", (DL_FUNC) &RxODE_current_fn_pointer_id);
  R_RegisterCCallable("RxODE","getRxSolve_", (DL_FUNC) &getRxSolve_);

  R_RegisterCCallable("RxODE", "gammap", (DL_FUNC) &gamma_p);
  R_RegisterCCallable("RxODE", "gammaq", (DL_FUNC) &gamma_q);
  R_RegisterCCallable("RxODE", "lowergamma", (DL_FUNC) &tgamma_lower);
  R_RegisterCCallable("RxODE", "uppergamma", (DL_FUNC) &tgamma_upper);
  R_RegisterCCallable("RxODE", "gammapDer", (DL_FUNC) &gamma_p_derivative);
  R_RegisterCCallable("RxODE", "gammapInv", (DL_FUNC) &gamma_p_inv);
  R_RegisterCCallable("RxODE", "gammapInva", (DL_FUNC) &gamma_p_inva);
  R_RegisterCCallable("RxODE", "gammaqInv", (DL_FUNC) &gamma_q_inv);
  R_RegisterCCallable("RxODE", "gammaqInva", (DL_FUNC) &gamma_q_inva);
  R_RegisterCCallable("RxODE", "compareFactorVal", (DL_FUNC) &compareFactorVal);
  R_RegisterCCallable("RxODE", "iniSubjectE", (DL_FUNC) &iniSubjectE);
  R_RegisterCCallable("RxODE", "getRxLhs", (DL_FUNC) &getRxLhs);
  R_RegisterCCallable("RxODE", "getUpdateInis", (DL_FUNC) &getUpdateInis);
  R_RegisterCCallable("RxODE", "sortIds", (DL_FUNC) &sortIds);
  R_RegisterCCallable("RxODE", "handleTlast", (DL_FUNC) &handleTlast);
  R_RegisterCCallable("RxODE", "rxGetId", (DL_FUNC) &rxGetId);
  R_RegisterCCallable("RxODE", "getTime", (DL_FUNC) &getTime);
  
  static const R_CMethodDef cMethods[] = {
    {"RxODE_sum",               (DL_FUNC) &RxODE_sum, 2, RxODE_Sum_t},
    {"RxODE_prod",              (DL_FUNC) &RxODE_prod, 2, RxODE_Sum_t},
    {NULL, NULL, 0, NULL}
  };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  rxOptionsIni();
  initRxThreads();
  avoid_openmp_hang_within_fork();
  transIniNull();
  nullGlobals();
  /* rxOptionsIniFocei(); */
}

void parseFree(int last);
void rxOptionsFree();
void gFree();
void rxFreeLast();
void R_unload_RxODE(DllInfo *info){
  gFree();
  rxOptionsFree();
  rxOptionsIni();
  parseFree(1);
  rxFreeLast();
}
