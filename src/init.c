#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "solve.h"

static R_NativePrimitiveArgType RxODE_sign_exp_t[] = {
  REALSXP, REALSXP
};

static R_NativePrimitiveArgType RxODE_one_int_t[] = {
  INTSXP
};

static R_NativePrimitiveArgType RxODE_one_dbl_t[] = {
  REALSXP
};

SEXP trans(SEXP orig_file, SEXP parse_file, SEXP c_file, SEXP extra_c, SEXP prefix, SEXP model_md5, SEXP parse_model,SEXP parse_model3);
SEXP _RxODE_linCmtEnv(SEXP rho);
SEXP _RxODE_rxInv(SEXP matrix);
SEXP _RxODE_removableDrive(SEXP letter);
SEXP _RxODE_rxCoutEcho(SEXP number);
SEXP _RxODE_RxODE_finalize_focei_omega(SEXP);
SEXP _RxODE_RxODE_finalize_log_det_OMGAinv_5(SEXP);
SEXP _rxCholInv(SEXP dms, SEXP theta, SEXP tn);
SEXP _RxODE_rxSymInvCholEnvCalculate(SEXP, SEXP, SEXP);
SEXP _RxODE_rxInvWishartVar(SEXP, SEXP);
SEXP _RxODE_rxSymInvChol(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxDataSetup(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxIs(SEXP,SEXP);
SEXP _RxODE_rxModelVars_(SEXP);
SEXP _RxODE_rxState(SEXP, SEXP);
SEXP _RxODE_rxParams(SEXP);
SEXP _RxODE_rxDfdy(SEXP);
SEXP _RxODE_rxLhs(SEXP);
SEXP _RxODE_rxInits(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxSetupIni(SEXP, SEXP);
SEXP _RxODE_rxSetupScale(SEXP,SEXP,SEXP);
SEXP _RxODE_rxDataParSetup(SEXP, SEXP, SEXP, SEXP, SEXP,
                           SEXP, SEXP, SEXP, SEXP, SEXP,
                           SEXP, SEXP, SEXP, SEXP, SEXP,
			   SEXP);
/* SEXP _RxODE_rxSolvingData(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, */
/* 			  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); */
/* SEXP _RxODE_rxData(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP, */
/* 		   SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP, */
/* 		   SEXP,SEXP,SEXP,SEXP,SEXP); */
/* SEXP _RxODE_rxSolveC(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP, */
/* 		     SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP, */
/* 		     SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP, */
/* 		     SEXP); */

SEXP _RxODE_rxSolveCsmall(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP _RxODE_rxSolveGet(SEXP, SEXP, SEXP);
SEXP _RxODE_rxSolveUpdate(SEXP, SEXP, SEXP);
SEXP _RxODE_rxAssignPtr(SEXP);
SEXP _RxODE_rxCores();
SEXP _RxODE_rxAssignPtr(SEXP objectSEXP);
SEXP RxODE_get_mv();

double RxODE_solveLinB(rx_solve *rx, unsigned int id,double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);

SEXP _RxODE_rxToOmega(SEXP cholInv);

static R_NativePrimitiveArgType RxODE_Sum_t[] = {
  REALSXP, INTSXP
};

extern int RxODE_current_fn_pointer_id();
extern double RxODE_as_zero(double x);
extern double RxODE_safe_log(double x);
extern double RxODE_safe_zero(double x);
extern double RxODE_pow(double x, double y);
extern double RxODE_pow_di(double x, int i);
extern double RxODE_sign_exp(double sgn, double x);
extern double RxODE_abs_log(double x);
extern double RxODE_abs_log1p(double x);
extern double RxODE_factorial(double x);
extern double RxODE_sum(double *input, int len);
extern double RxODE_prod(double *input, int len);
extern void RxODE_ode_solve_env(SEXP sexp_rho);
extern int nEq ();
extern unsigned int nObs();
extern unsigned int nLhs ();
extern double RxODE_as_zero(double x);
extern double rxLhs(int i);
extern void rxCalcLhs(int i);
extern unsigned int nAllTimes ();
extern int rxEvid(int i);

// Changed for Parallel
extern int nEqP (rx_solve *rx, unsigned int id);
extern unsigned int nObsP(rx_solve *rx, unsigned int id);
extern unsigned int nLhsP (rx_solve *rx, unsigned int id);
extern double rxLhsP(int i, rx_solve *rx, unsigned int id);
extern void rxCalcLhsP(int i, rx_solve *rx, unsigned int id);
extern unsigned int nAllTimesP (rx_solve *rx, unsigned int id);
extern int rxEvidP(int i, rx_solve *rx, unsigned int id);

extern void RxODE_assign_fn_pointers(SEXP mv);

extern void rxSolveOldC(SEXP object, 
			int *neqa,
			double *theta,  //order:
			double *timep,
			int *evidp,
			int *ntime,
			double *initsp,
			double *dosep,
			double *retp,
			double *atol,
			double *rtol,
			int *stiffa,
			int *transit_abs,
			int *nlhsa,
			double *lhsp,
			int *rc);

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

extern SEXP _RxODE_rxSimThetaOmega(SEXP, SEXP, SEXP, SEXP, SEXP,
                                   SEXP, SEXP, SEXP, SEXP, SEXP,
                                   SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _RxODE_cvPost(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _RxODE_rinvchisq(SEXP, SEXP, SEXP);

extern int rxIsCurrentC(SEXP obj);


// Remove these functions later...

void rxOptionsIni();
void R_init_RxODE(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"trans", (DL_FUNC) &trans, 8},
    {"RxODE_get_mv", (DL_FUNC) &RxODE_get_mv, 0},
    {"_RxODE_rxInv", (DL_FUNC) &_RxODE_rxInv, 1},
    {"_RxODE_RxODE_finalize_focei_omega",(DL_FUNC) &_RxODE_RxODE_finalize_focei_omega, 1},
    {"_RxODE_RxODE_finalize_log_det_OMGAinv_5",(DL_FUNC) &_RxODE_RxODE_finalize_log_det_OMGAinv_5, 1},
    {"_RxODE_rxCoutEcho", (DL_FUNC) &_RxODE_rxCoutEcho, 1},
    {"_RxODE_removableDrive", (DL_FUNC) &_RxODE_removableDrive, 1},
    {"_rxCholInv", (DL_FUNC) &_rxCholInv, 3},
    {"_RxODE_rxToOmega", (DL_FUNC) &_RxODE_rxToOmega, 1},
    {"_RxODE_rxSymInvCholEnvCalculate", (DL_FUNC) &_RxODE_rxSymInvCholEnvCalculate, 3},
    {"_RxODE_rxInvWishartVar", (DL_FUNC) &_RxODE_rxInvWishartVar, 2},
    {"_RxODE_rxSymInvChol", (DL_FUNC) &_RxODE_rxSymInvChol, 4},
    {"_RxODE_rxDataSetup", (DL_FUNC) &_RxODE_rxDataSetup, 9},
    {"_RxODE_rxIs", (DL_FUNC) &_RxODE_rxIs, 2},
    {"_RxODE_rxModelVars_", (DL_FUNC) &_RxODE_rxModelVars_, 1},
    {"_RxODE_rxState", (DL_FUNC) &_RxODE_rxState, 2},
    {"_RxODE_rxParams", (DL_FUNC) &_RxODE_rxParams, 1},
    {"_RxODE_rxDfdy", (DL_FUNC) &_RxODE_rxDfdy, 1},
    {"_RxODE_rxLhs", (DL_FUNC) &_RxODE_rxLhs, 1},
    {"_RxODE_rxInits", (DL_FUNC) &_RxODE_rxInits, 6},
    {"_RxODE_rxSetupIni", (DL_FUNC) &_RxODE_rxSetupIni, 2},
    {"_RxODE_rxSetupScale", (DL_FUNC) &_RxODE_rxSetupScale, 3},
    {"_RxODE_rxDataParSetup", (DL_FUNC) &_RxODE_rxDataParSetup, 16},
    // Solaris needs 23 args; fix me...
    /* {"_RxODE_rxSolveC", (DL_FUNC) &_RxODE_rxSolveC, 31}, */
    {"_RxODE_rxSolveCsmall", (DL_FUNC) &_RxODE_rxSolveCsmall, 9},
    {"_RxODE_rxSolveGet", (DL_FUNC) &_RxODE_rxSolveGet, 3},
    {"_RxODE_rxSolveUpdate", (DL_FUNC) &_RxODE_rxSolveUpdate, 3},
    {"_RxODE_rxCores",(DL_FUNC) &_RxODE_rxCores, 0},
    {"_RxODE_rxAssignPtr", (DL_FUNC) &_RxODE_rxAssignPtr, 1},
    {"_RxODE_rxRmModelLib_",(DL_FUNC) &_RxODE_rxRmModelLib_, 1},
    {"_RxODE_rxDll",(DL_FUNC) &_RxODE_rxDll, 1},
    {"_RxODE_rxC",(DL_FUNC) &_RxODE_rxC, 1},
    {"_RxODE_rxIsLoaded", (DL_FUNC) &_RxODE_rxIsLoaded, 1},
    {"_RxODE_rxDynUnload", (DL_FUNC) &_RxODE_rxDynUnload, 1},
    {"_RxODE_rxDynLoad", (DL_FUNC) &_RxODE_rxDynLoad, 1},
    {"_RxODE_rxDelete", (DL_FUNC) &_RxODE_rxDelete, 1},
    {"_RxODE_rxGetRxODE", (DL_FUNC) &_RxODE_rxGetRxODE, 1},
    {"_RxODE_rxSimThetaOmega", (DL_FUNC) &_RxODE_rxSimThetaOmega, 15},
    {"_RxODE_rxIsCurrent", (DL_FUNC) &_RxODE_rxIsCurrent, 1},
    {"_RxODE_cvPost", (DL_FUNC) &_RxODE_cvPost, 5},
    {"_RxODE_rinvchisq", (DL_FUNC) &_RxODE_rinvchisq, 3},
    {NULL, NULL, 0}
  };
  // C callable to assign environments.
  R_RegisterCCallable("RxODE","rxRmModelLib", (DL_FUNC) rxRmModelLib);
  R_RegisterCCallable("RxODE","rxGetModelLib", (DL_FUNC) rxGetModelLib);
  // C callables needed in FOCEi
  R_RegisterCCallable("RxODE","nEq",                      (DL_FUNC) nEq);
  R_RegisterCCallable("RxODE","nLhs",                     (DL_FUNC) nLhs);
  R_RegisterCCallable("RxODE","rxLhs",                    (DL_FUNC) rxLhs);
  R_RegisterCCallable("RxODE","nAllTimes",                (DL_FUNC) nAllTimes);
  R_RegisterCCallable("RxODE","rxEvid",                   (DL_FUNC) rxEvid);
  R_RegisterCCallable("RxODE","rxCalcLhs",                (DL_FUNC) rxCalcLhs);
  R_RegisterCCallable("RxODE","nObs",                     (DL_FUNC) nObs);

  R_RegisterCCallable("RxODE","nEqP",                      (DL_FUNC) nEqP);
  R_RegisterCCallable("RxODE","nLhsP",                     (DL_FUNC) nLhsP);
  R_RegisterCCallable("RxODE","rxLhsP",                    (DL_FUNC) rxLhsP);
  R_RegisterCCallable("RxODE","nAllTimesP",                (DL_FUNC) nAllTimesP);
  R_RegisterCCallable("RxODE","rxEvidP",                   (DL_FUNC) rxEvidP);
  R_RegisterCCallable("RxODE","rxCalcLhsP",                (DL_FUNC) rxCalcLhsP);
  R_RegisterCCallable("RxODE","nObsP",                     (DL_FUNC) nObsP);

  R_RegisterCCallable("RxODE","RxODE_ode_solve_env",      (DL_FUNC) RxODE_ode_solve_env);
  R_RegisterCCallable("RxODE","RxODE_ode_free",           (DL_FUNC) RxODE_ode_free);
  R_RegisterCCallable("RxODE","RxODE_safe_zero",          (DL_FUNC) RxODE_safe_zero);
  R_RegisterCCallable("RxODE","RxODE_safe_log",           (DL_FUNC) RxODE_safe_log);
  R_RegisterCCallable("RxODE","RxODE_sign_exp",           (DL_FUNC) RxODE_sign_exp);
  R_RegisterCCallable("RxODE","RxODE_abs_log",            (DL_FUNC) RxODE_abs_log);
  
  //Functions
  R_RegisterCCallable("RxODE","rxSolveOldC",              (DL_FUNC) rxSolveOldC);
  
  // tranit compartment models
  R_RegisterCCallable("RxODE","RxODE_factorial",          (DL_FUNC) RxODE_factorial);
  R_RegisterCCallable("RxODE","RxODE_safe_log",           (DL_FUNC) RxODE_safe_log);
  R_RegisterCCallable("RxODE","RxODE_safe_zero",          (DL_FUNC) RxODE_safe_zero);
  R_RegisterCCallable("RxODE","RxODE_as_zero",            (DL_FUNC) RxODE_as_zero);
  R_RegisterCCallable("RxODE","RxODE_sign_exp",           (DL_FUNC) RxODE_sign_exp);
  R_RegisterCCallable("RxODE","RxODE_abs_log",            (DL_FUNC) RxODE_abs_log);
  R_RegisterCCallable("RxODE","RxODE_abs_log1p",          (DL_FUNC) RxODE_abs_log1p);
  R_RegisterCCallable("RxODE","RxODE_solveLinB",          (DL_FUNC) RxODE_solveLinB);

  R_RegisterCCallable("RxODE","RxODE_sum",                (DL_FUNC) RxODE_sum);
  R_RegisterCCallable("RxODE","RxODE_prod",               (DL_FUNC) RxODE_prod);

  R_RegisterCCallable("RxODE","RxODE_pow",                (DL_FUNC) RxODE_pow);
  R_RegisterCCallable("RxODE","RxODE_pow_di",             (DL_FUNC) RxODE_pow_di);
  R_RegisterCCallable("RxODE","RxODE_assign_fn_pointers", (DL_FUNC) &RxODE_assign_fn_pointers);

  R_RegisterCCallable("RxODE","_RxODE_rxAssignPtr",       (DL_FUNC) _RxODE_rxAssignPtr);
  R_RegisterCCallable("RxODE", "rxIsCurrentC", (DL_FUNC) rxIsCurrentC);
  R_RegisterCCallable("RxODE","RxODE_current_fn_pointer_id", (DL_FUNC) &RxODE_current_fn_pointer_id);

  
  static const R_CMethodDef cMethods[] = {
    {"RxODE_factorial",         (DL_FUNC) &RxODE_factorial, 1, RxODE_one_dbl_t},
    {"RxODE_safe_log",          (DL_FUNC) &RxODE_safe_log, 1, RxODE_one_dbl_t},
    {"RxODE_safe_zero",         (DL_FUNC) &RxODE_safe_zero, 1, RxODE_one_dbl_t},
    {"RxODE_as_zero",           (DL_FUNC) &RxODE_as_zero, 1, RxODE_one_dbl_t},
    {"RxODE_sign_exp",          (DL_FUNC) &RxODE_sign_exp, 2, RxODE_sign_exp_t},
    {"RxODE_abs_log",           (DL_FUNC) &RxODE_abs_log, 1, RxODE_one_dbl_t},
    {"RxODE_abs_log1p",         (DL_FUNC) &RxODE_abs_log1p, 1, RxODE_one_dbl_t},
    {"RxODE_sum",               (DL_FUNC) &RxODE_sum, 2, RxODE_Sum_t},
    {"RxODE_prod",              (DL_FUNC) &RxODE_prod, 2, RxODE_Sum_t},
    {NULL, NULL, 0, NULL}
  };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  rxOptionsIni();
}

void rxOptionsFree();
void R_unload_RxODE(DllInfo *info){
  rxOptionsFree();
}
