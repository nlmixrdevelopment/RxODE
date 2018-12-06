#include <RxODE_model.h>

extern void __ODE_SOLVER_SOLVEDATA__ (rx_solve *solve){
  _solveData = solve;
}

extern rx_solve *__ODE_SOLVER_GET_SOLVEDATA__(){
  return _solveData;
}

SEXP __MODEL_VARS__();
extern void __ODE_SOLVER__(int *neq,
			   double *theta,      //order:
			   double *time,
			   int *evid,
			   int *ntime,
			   double *inits,
			   double *dose,
			   double *ret,
			   double *atol,
			   double *rtol,
			   int *stiff,
			   int *transit_abs,
			   int *nlhs,
			   double *lhs,
			   int *rc){
  // Backward compatible ode solver for 0.5* C interface
  //if (_ptrid() != __TIMEID__ ){ _assign_ptr(__MODEL_VARS__());}
  __FIX_INIS__
  _old_c(neq, _theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc);
}

static R_NativePrimitiveArgType __ODE_SOLVER__rx_t[] = {
  //*neq, *theta, *time,  *evid, *ntime, *inits,   *dose,   *ret,     *atol,  *rtol,   *stiff, *transit_abs, *nlhs, *lhs, *rc
  INTSXP,REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP
};

// CODE HERE

extern void __DYDT_LSODA__(int *neq, double *t, double *A, double *DADT)
{
  __DYDT__(neq, *t, A, DADT);
}

extern int __DYDT_LIBLSODA__(double t, double *y, double *ydot, void *data)
{
  int *neq = (int*)(data);
  __DYDT__(neq, t, y, ydot);
  return(0);
}

extern void __CALC_JAC_LSODA__(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){
  // Update all covariate parameters
  __CALC_JAC__(neq, *t, A, JAC, *nrowpd);
}

//Initilize the dll to match RxODE's calls
void __R_INIT__ (DllInfo *info){
  // Get C callables on load; Otherwise it isn't thread safe
  _assign_ptr=(RxODE_assign_ptr) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
  _rxRmModelLib=(_rxRmModelLibType) R_GetCCallable("RxODE","rxRmModelLib");
  _rxGetModelLib=(_rxGetModelLibType) R_GetCCallable("RxODE","rxGetModelLib");
  _old_c = (RxODE_ode_solver_old_c) R_GetCCallable("RxODE","rxSolveOldC");
  _RxODE_rxAssignPtr=(_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxAssignPtr");
  _rxIsCurrentC = (_rxIsCurrentC_type)R_GetCCallable("RxODE","rxIsCurrentC");
  _sumPS  = (_rxSumType) R_GetCCallable("PreciseSums","PreciseSums_sum_r");
  _prodPS = (_rxProdType) R_GetCCallable("PreciseSums","PreciseSums_prod_r");
  _prodType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_prod_get");
  _sumType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_sum_get");
  _ptrid=(RxODE_fn0i)R_GetCCallable("RxODE", "RxODE_current_fn_pointer_id");
  _powerD=(RxODE_fn3i)R_GetCCallable("RxODE", "powerD");
  _powerDi=(RxODE_fn3i)R_GetCCallable("RxODE", "powerDi");
  _powerDD=(RxODE_fn3i)R_GetCCallable("RxODE", "powerDD");
  _powerDDD=(RxODE_fn3i)R_GetCCallable("RxODE", "powerDDD");
  solveLinB=(solveLinB_p)R_GetCCallable("RxODE", "solveLinB");
  _update_par_ptr=(_update_par_ptr_p) R_GetCCallable("RxODE","_update_par_ptr");
  // Register the outside functions
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_STR__,       (DL_FUNC) __ODE_SOLVER__);
  R_RegisterCCallable(__LIB_STR__,__INIS_STR__, (DL_FUNC) __INIS__);
  R_RegisterCCallable(__LIB_STR__,__INIS_STR__, (DL_FUNC) __INIS__);
  R_RegisterCCallable(__LIB_STR__,__DYDT_STR__, (DL_FUNC) __DYDT__);
  R_RegisterCCallable(__LIB_STR__,__CALC_LHS_STR__, (DL_FUNC) __CALC_LHS__);
  R_RegisterCCallable(__LIB_STR__,__CALC_JAC_STR__, (DL_FUNC) __CALC_JAC__);
  R_RegisterCCallable(__LIB_STR__,__DYDT_LSODA_STR__, (DL_FUNC) __DYDT_LSODA__);
  R_RegisterCCallable(__LIB_STR__,__CALC_JAC_LSODA_STR__, (DL_FUNC) __CALC_JAC_LSODA__);
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_SOLVEDATA_STR__, (DL_FUNC) __ODE_SOLVER_SOLVEDATA__);
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_GET_SOLVEDATA_STR__, (DL_FUNC) __ODE_SOLVER_GET_SOLVEDATA__);
  R_RegisterCCallable(__LIB_STR__,__DYDT_LIBLSODA_STR__, (DL_FUNC) __DYDT_LIBLSODA__);
  
  static const R_CMethodDef cMethods[] = {
    {__ODE_SOLVER_STR__, (DL_FUNC) &__ODE_SOLVER__, 15, __ODE_SOLVER__rx_t},
    {NULL, NULL, 0, NULL}
  };
  
  R_CallMethodDef callMethods[]  = {
    {__MODEL_VARS_STR__, (DL_FUNC) &__MODEL_VARS__, 0},
    {NULL, NULL, 0}
  };
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
}

void __R_UNLOAD__ (DllInfo *info){
  // Free resources required for single subject solve.
  SEXP _mv = PROTECT(_rxGetModelLib(__MODEL_VARS_STR__));
  if (!isNull(_mv)){
    _rxRmModelLib(__MODEL_VARS_STR__);
  }
  UNPROTECT(1);
}
