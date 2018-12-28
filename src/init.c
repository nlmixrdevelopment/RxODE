#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "../inst/include/RxODE.h"

double powerDi(double x, double lambda, int yj);
double powerD(double x, double lambda, int yj);
double powerDD(double x, double lambda, int yj);
double powerDDD(double x, double lambda, int yj);
double powerL(double x, double lambda, int yj);
double powerDL(double x, double lambda, int yj);

SEXP _rxProgress(SEXP num, SEXP core);
SEXP _rxTick();
SEXP _rxProgressStop(SEXP);
SEXP _rxProgressAbort();
SEXP _RxODE_codeLoaded();

SEXP _RxODE_trans(SEXP parse_file, SEXP extra_c, SEXP prefix, SEXP model_md5, SEXP isStr);
SEXP _RxODE_codegen(SEXP c_file, SEXP prefix, SEXP libname, SEXP pMd5, SEXP timeId,
		    SEXP fixInis);
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
SEXP _RxODE_rxSolveCsmall(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP _RxODE_rxSolveGet(SEXP, SEXP, SEXP);
SEXP _RxODE_rxSolveUpdate(SEXP, SEXP, SEXP);
SEXP _RxODE_rxAssignPtr(SEXP);
SEXP _RxODE_rxCores();
SEXP _RxODE_rxAssignPtr(SEXP objectSEXP);
SEXP _RxODE_dynLoad(SEXP dllSEXP);
SEXP _RxODE_rxSolveC(SEXP, SEXP, SEXP, SEXP, SEXP, //5
		      SEXP, SEXP, SEXP, SEXP, SEXP, //10
		      SEXP, SEXP, SEXP, SEXP, SEXP, //15
		      SEXP, SEXP, SEXP, SEXP, SEXP, //20
		      SEXP, SEXP, SEXP, SEXP, SEXP, //25
		      SEXP, SEXP, SEXP, SEXP, SEXP, //30
		      SEXP, SEXP, SEXP, SEXP, SEXP, //35
		      SEXP, SEXP, SEXP, SEXP, SEXP, //40
		      SEXP, SEXP, SEXP, SEXP, SEXP, //45
		      SEXP);
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

SEXP _RxODE_cvPost(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _RxODE_rinvchisq(SEXP, SEXP, SEXP);
SEXP _RxODE_add_dosing_(SEXP, SEXP, SEXP, SEXP, SEXP,
                        SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_add_sampling_(SEXP, SEXP, SEXP);

SEXP _RxODE_getRxFn(SEXP);
SEXP _RxODE_setProgSupported(SEXP);
SEXP _RxODE_getProgSupported();

extern int rxIsCurrentC(SEXP obj);

rx_solve *getRxSolve_();

// Remove these functions later...

void rxOptionsIni();
void rxOptionsIniData();
/* void rxOptionsIniFocei(); */

double solveLinB(rx_solve *rx, unsigned int id, double t, int linCmt, int diff1, int diff2, double d_A, double d_alpha, double d_B, double d_beta, double d_C, double d_gamma, double d_ka, double d_tlag);

void _update_par_ptr(double t, unsigned int id, rx_solve *rx, int idx);

int par_progress(int c, int n, int d, int cores, clock_t t0, int stop);
void ind_solve(rx_solve *rx, unsigned int cid, t_dydt_liblsoda dydt_lls, 
	       t_dydt_lsoda_dum dydt_lsoda, t_jdum_lsoda jdum,
	       t_dydt c_dydt, t_update_inis u_inis, int jt);
int isRstudio();

void R_init_RxODE(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"_rxProgress", (DL_FUNC) &_rxProgress, 2},
    {"_rxTick", (DL_FUNC) &_rxTick, 0},
    {"_rxProgressStop", (DL_FUNC) &_rxProgressStop, 1},
    {"_rxProgressAbort", (DL_FUNC) &_rxProgressAbort, 0},
    {"_RxODE_trans", (DL_FUNC) &_RxODE_trans, 5},
    {"_RxODE_codegen", (DL_FUNC) &_RxODE_codegen, 6},
    {"_RxODE_codeLoaded", (DL_FUNC) &_RxODE_codeLoaded, 0},
    {"_RxODE_parseModel", (DL_FUNC) &_RxODE_parseModel, 1},
    {"_RxODE_isLinCmt", (DL_FUNC) &_RxODE_isLinCmt, 0},
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
    {"_RxODE_setRstudio", (DL_FUNC) &_RxODE_setRstudio, 1},
    {"_RxODE_RcppExport_registerCCallable", (DL_FUNC) &_RxODE_RcppExport_registerCCallable, 0},
    // Solaris needs 23 args; fix me...
    {"_RxODE_rxSolveC", (DL_FUNC) &_RxODE_rxSolveC, 46},
    {"_RxODE_getRxFn", (DL_FUNC) &_RxODE_getRxFn, 1},
    {"_RxODE_setProgSupported", (DL_FUNC) &_RxODE_setProgSupported, 1},
    {"_RxODE_getProgSupported", (DL_FUNC) &_RxODE_getProgSupported, 0},
    {NULL, NULL, 0}
  };
  // C callable to assign environments.
  R_RegisterCCallable("RxODE", "powerDi", (DL_FUNC) powerDi);
  R_RegisterCCallable("RxODE", "powerD", (DL_FUNC) powerD);
  R_RegisterCCallable("RxODE", "powerDD", (DL_FUNC) powerDD);
  R_RegisterCCallable("RxODE", "powerDDD", (DL_FUNC) powerDDD);
  R_RegisterCCallable("RxODE", "powerL", (DL_FUNC) powerL);
  R_RegisterCCallable("RxODE", "powerDL", (DL_FUNC) powerDL);
  R_RegisterCCallable("RxODE", "par_progress", (DL_FUNC) par_progress);
  R_RegisterCCallable("RxODE", "isRstudio", (DL_FUNC) isRstudio);
  R_RegisterCCallable("RxODE", "ind_solve", (DL_FUNC) ind_solve);
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
  R_RegisterCCallable("RxODE","getRxSolve_", (DL_FUNC) &getRxSolve_);
  
  static const R_CMethodDef cMethods[] = {
    {"RxODE_sum",               (DL_FUNC) &RxODE_sum, 2, RxODE_Sum_t},
    {"RxODE_prod",              (DL_FUNC) &RxODE_prod, 2, RxODE_Sum_t},
    {NULL, NULL, 0, NULL}
  };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  rxOptionsIni();
  rxOptionsIniData();
  /* rxOptionsIniFocei(); */
}


void rxOptionsFree();
void gFree();
/* void rxOptionsFreeFocei(); */
void R_unload_RxODE(DllInfo *info){
  rxOptionsFree();
  gFree();
  /* rxOptionsFreeFocei(); */
}
