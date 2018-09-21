#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "RxODE.h"

SEXP _rxProgress(SEXP num, SEXP core);
SEXP _rxTick();
SEXP _rxProgressStop(SEXP);
SEXP _rxProgressAbort();

SEXP trans(SEXP orig_file, SEXP parse_file, SEXP c_file, SEXP extra_c, SEXP prefix, SEXP model_md5, SEXP parse_model,SEXP parse_model3);
SEXP _RxODE_sqrtm(SEXP);
SEXP _RxODE_foceiOuterG(SEXP);
SEXP _RxODE_foceiOuterF(SEXP);
SEXP _RxODE_cholSE_(SEXP, SEXP);
SEXP _RxODE_cholSE0(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_setRstudio(SEXP);
SEXP _RxODE_coxBox_(SEXP, SEXP, SEXP);
SEXP _RxODE_foceiFitCpp_(SEXP);
SEXP _RxODE_foceiCalcCov(SEXP);
SEXP _RxODE_foceiOuter(SEXP);
SEXP _RxODE_foceiEtas();
SEXP _RxODE_foceiNumericGrad(SEXP);
SEXP _RxODE_foceiLik(SEXP);
SEXP _RxODE_foceiOfv(SEXP);
SEXP _RxODE_likInner(SEXP, SEXP);
SEXP _RxODE_foceiInnerLp(SEXP, SEXP);
SEXP _RxODE_foceiSetup_(SEXP,SEXP, SEXP, SEXP, SEXP,
                        SEXP,SEXP, SEXP, SEXP, SEXP);
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

extern double powerD(double x, double lambda, int yj);
extern double powerDD(double x, double lambda, int yj);
extern double powerDDD(double x, double lambda, int yj);
extern double powerDi(double x, double lambda, int yj);

SEXP _RxODE_cvPost(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _RxODE_rinvchisq(SEXP, SEXP, SEXP);
SEXP _RxODE_add_dosing_(SEXP, SEXP, SEXP, SEXP, SEXP,
                        SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_add_sampling_(SEXP, SEXP, SEXP);

SEXP _RxODE_rxSolveFree();

extern int rxIsCurrentC(SEXP obj);


// Remove these functions later...

void rxOptionsIni();
void rxOptionsIniData();
void rxOptionsIniFocei();

double solveLinB(rx_solve *rx, unsigned int id, double t, int linCmt, int diff1, int diff2, double d_A, double d_alpha, double d_B, double d_beta, double d_C, double d_gamma, double d_ka, double d_tlag);

void _update_par_ptr(double t, unsigned int id, rx_solve *rx, int idx);

void R_init_RxODE(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"_rxProgress", (DL_FUNC) &_rxProgress, 2},
    {"_rxTick", (DL_FUNC) &_rxTick, 0},
    {"_rxProgressStop", (DL_FUNC) &_rxProgressStop, 1},
    {"_rxProgressAbort", (DL_FUNC) &_rxProgressAbort, 0},
    {"trans", (DL_FUNC) &trans, 8},
    {"RxODE_get_mv", (DL_FUNC) &RxODE_get_mv, 0},
    {"_RxODE_rxInv", (DL_FUNC) &_RxODE_rxInv, 1},
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
    {"_RxODE_rxSolveFree", (DL_FUNC) &_RxODE_rxSolveFree, 0},
    {"_RxODE_foceiSetup_", (DL_FUNC) &_RxODE_foceiSetup_, 10},
    {"_RxODE_foceiLik", (DL_FUNC) &_RxODE_foceiLik, 1},
    {"_RxODE_foceiOfv", (DL_FUNC) &_RxODE_foceiOfv, 1},
    {"_RxODE_likInner", (DL_FUNC) &_RxODE_likInner, 2},
    {"_RxODE_foceiInnerLp", (DL_FUNC) &_RxODE_foceiInnerLp, 2},
    {"_RxODE_foceiNumericGrad", (DL_FUNC) &_RxODE_foceiNumericGrad, 1},
    {"_RxODE_foceiEtas", (DL_FUNC) &_RxODE_foceiEtas, 0},
    {"_RxODE_foceiOuter", (DL_FUNC) &_RxODE_foceiOuter, 1},
    {"_RxODE_foceiCalcCov", (DL_FUNC) &_RxODE_foceiCalcCov, 1},
    {"_RxODE_foceiFitCpp_", (DL_FUNC) &_RxODE_foceiFitCpp_, 1},
    {"_RxODE_coxBox_", (DL_FUNC) &_RxODE_coxBox_, 3},
    {"_RxODE_setRstudio", (DL_FUNC) &_RxODE_setRstudio, 1},
    {"_RxODE_cholSE_", (DL_FUNC) &_RxODE_cholSE_, 2},
    {"_RxODE_cholSE0", (DL_FUNC) &_RxODE_cholSE0, 4},
    {"_RxODE_foceiOuterG", (DL_FUNC) &_RxODE_foceiOuterG, 1},
    {"_RxODE_foceiOuterF", (DL_FUNC) &_RxODE_foceiOuterF, 1},
    {"_RxODE_sqrtm", (DL_FUNC) &_RxODE_sqrtm, 1},
    {NULL, NULL, 0}
  };
  // C callable to assign environments.
  R_RegisterCCallable("RxODE", "solveLinB", (DL_FUNC) solveLinB);
  R_RegisterCCallable("RxODE", "_update_par_ptr", (DL_FUNC) _update_par_ptr);
  R_RegisterCCallable("RxODE","rxRmModelLib", (DL_FUNC) rxRmModelLib);
  R_RegisterCCallable("RxODE","rxGetModelLib", (DL_FUNC) rxGetModelLib);
  
  R_RegisterCCallable("RxODE","RxODE_ode_free",           (DL_FUNC) RxODE_ode_free);
  
  //Functions
  R_RegisterCCallable("RxODE","rxSolveOldC",              (DL_FUNC) rxSolveOldC);
  
  R_RegisterCCallable("RxODE","RxODE_sum",                (DL_FUNC) RxODE_sum);
  R_RegisterCCallable("RxODE","RxODE_prod",               (DL_FUNC) RxODE_prod);

  R_RegisterCCallable("RxODE","RxODE_assign_fn_pointers", (DL_FUNC) &RxODE_assign_fn_pointers);

  R_RegisterCCallable("RxODE","_RxODE_rxAssignPtr",       (DL_FUNC) _RxODE_rxAssignPtr);
  R_RegisterCCallable("RxODE", "rxIsCurrentC", (DL_FUNC) rxIsCurrentC);
  R_RegisterCCallable("RxODE","RxODE_current_fn_pointer_id", (DL_FUNC) &RxODE_current_fn_pointer_id);
  R_RegisterCCallable("RxODE", "powerD", (DL_FUNC) &powerD);
  R_RegisterCCallable("RxODE", "powerDD", (DL_FUNC) &powerDD);
  R_RegisterCCallable("RxODE", "powerDDD", (DL_FUNC) &powerDDD);
  R_RegisterCCallable("RxODE", "powerDi", (DL_FUNC) &powerDi);

  
  static const R_CMethodDef cMethods[] = {
    {"RxODE_sum",               (DL_FUNC) &RxODE_sum, 2, RxODE_Sum_t},
    {"RxODE_prod",              (DL_FUNC) &RxODE_prod, 2, RxODE_Sum_t},
    {NULL, NULL, 0, NULL}
  };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  rxOptionsIni();
  rxOptionsIniData();
  rxOptionsIniFocei();
}


void rxOptionsFree();
void gFree();
void rxOptionsFreeFocei();
void R_unload_RxODE(DllInfo *info){
  rxOptionsFree();
  gFree();
  rxOptionsFreeFocei();
}
