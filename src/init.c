#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "solve.h"

SEXP trans(SEXP orig_file, SEXP parse_file, SEXP c_file, SEXP extra_c, SEXP prefix, SEXP model_md5, SEXP parse_model,SEXP parse_model3);
SEXP _RxODE_linCmtEnv(SEXP rho);
SEXP _RxODE_rxInv(SEXP matrix);
SEXP _RxODE_removableDrive(SEXP letter);
SEXP _RxODE_RxODE_finalize_focei_omega(SEXP);
SEXP _RxODE_RxODE_finalize_log_det_OMGAinv_5(SEXP);
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
SEXP _RxODE_rxSolveCsmall(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP _RxODE_rxSolveGet(SEXP, SEXP, SEXP);
SEXP _RxODE_rxSolveUpdate(SEXP, SEXP, SEXP);
SEXP _RxODE_rxAssignPtr(SEXP);
SEXP _RxODE_rxCores();
SEXP _RxODE_rxAssignPtr(SEXP objectSEXP);
SEXP _RxODE_dynLoad(SEXP dllSEXP);
SEXP RxODE_get_mv();

SEXP _RxODE_rxToOmega(SEXP cholInv);

static R_NativePrimitiveArgType RxODE_Sum_t[] = {
  REALSXP, INTSXP
};

extern int RxODE_current_fn_pointer_id();
extern double RxODE_sum(double *input, int len);
extern double RxODE_prod(double *input, int len);
extern void RxODE_ode_solve_env(SEXP sexp_rho);
extern int nEq ();
extern unsigned int nObs();
extern unsigned int nLhs ();
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
                                   SEXP, SEXP, SEXP, SEXP, SEXP,
				   SEXP, SEXP);

SEXP _RxODE_cvPost(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _RxODE_rinvchisq(SEXP, SEXP, SEXP);
SEXP _RxODE_add_dosing_(SEXP, SEXP, SEXP, SEXP, SEXP,
                        SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_add_sampling_(SEXP, SEXP, SEXP);

extern int rxIsCurrentC(SEXP obj);


// Remove these functions later...

void rxOptionsIni();
void rxOptionsIniData();
void R_init_RxODE(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"trans", (DL_FUNC) &trans, 8},
    {"RxODE_get_mv", (DL_FUNC) &RxODE_get_mv, 0},
    {"_RxODE_rxInv", (DL_FUNC) &_RxODE_rxInv, 1},
    {"_RxODE_RxODE_finalize_focei_omega",(DL_FUNC) &_RxODE_RxODE_finalize_focei_omega, 1},
    {"_RxODE_RxODE_finalize_log_det_OMGAinv_5",(DL_FUNC) &_RxODE_RxODE_finalize_log_det_OMGAinv_5, 1},
    {"_RxODE_removableDrive", (DL_FUNC) &_RxODE_removableDrive, 1},
    {"_rxCholInv", (DL_FUNC) &_rxCholInv, 3},
    {"_RxODE_rxToOmega", (DL_FUNC) &_RxODE_rxToOmega, 1},
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
    {"_RxODE_rxSimThetaOmega", (DL_FUNC) &_RxODE_rxSimThetaOmega, 17},
    {"_RxODE_rxIsCurrent", (DL_FUNC) &_RxODE_rxIsCurrent, 1},
    {"_RxODE_cvPost", (DL_FUNC) &_RxODE_cvPost, 5},
    {"_RxODE_rinvchisq", (DL_FUNC) &_RxODE_rinvchisq, 3},
    {"_RxODE_add_dosing_", (DL_FUNC) &_RxODE_add_dosing_,10},
    {"_RxODE_add_sampling_", (DL_FUNC) &_RxODE_add_sampling_, 3},
    {"_RxODE_dynLoad", (DL_FUNC) &_RxODE_dynLoad, 1},
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
  
  //Functions
  R_RegisterCCallable("RxODE","rxSolveOldC",              (DL_FUNC) rxSolveOldC);
  
  R_RegisterCCallable("RxODE","RxODE_sum",                (DL_FUNC) RxODE_sum);
  R_RegisterCCallable("RxODE","RxODE_prod",               (DL_FUNC) RxODE_prod);

  R_RegisterCCallable("RxODE","RxODE_assign_fn_pointers", (DL_FUNC) &RxODE_assign_fn_pointers);

  R_RegisterCCallable("RxODE","_RxODE_rxAssignPtr",       (DL_FUNC) _RxODE_rxAssignPtr);
  R_RegisterCCallable("RxODE", "rxIsCurrentC", (DL_FUNC) rxIsCurrentC);
  R_RegisterCCallable("RxODE","RxODE_current_fn_pointer_id", (DL_FUNC) &RxODE_current_fn_pointer_id);

  
  static const R_CMethodDef cMethods[] = {
    {"RxODE_sum",               (DL_FUNC) &RxODE_sum, 2, RxODE_Sum_t},
    {"RxODE_prod",              (DL_FUNC) &RxODE_prod, 2, RxODE_Sum_t},
    {NULL, NULL, 0, NULL}
  };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  rxOptionsIni();
  rxOptionsIniData();
}

void rxOptionsFree();
void gFree();
void R_unload_RxODE(DllInfo *info){
  rxOptionsFree();
  gFree();
}
