#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#define JAC_Rprintf Rprintf
#define JAC0_Rprintf if (_jac_counter_val() == 0) Rprintf
#define ODE_Rprintf Rprintf
#define ODE0_Rprintf if (_dadt_counter_val() == 0) Rprintf
#define LHS_Rprintf Rprintf
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#define R_pow Rx_pow
#define R_pow_di Rx_pow_di

// Types for par pointers.r
typedef void (*RxODE_update_par_ptr)(double t, rx_solve *rx, unsigned int id);
typedef double (*RxODE_transit3P)(double t, rx_solve *rx, unsigned int id, double n, double mtt);
typedef double (*RxODE_fn) (double x);
typedef double (*RxODE_fn2) (double x, double y);
typedef double (*RxODE_fn2i) (double x, int i);
typedef double (*RxODE_transit4P)(double t, rx_solve *rx, unsigned int id, double n, double mtt, double bio);
typedef double (*RxODE_vec) (int val, rx_solve *rx, unsigned int id);
typedef long (*RxODE_cnt) (rx_solve *rx, unsigned int id);
typedef void (*RxODE_inc) (rx_solve *rx, unsigned int id);
typedef double (*RxODE_val) (rx_solve *rx, unsigned int id);
typedef SEXP (*RxODE_ode_solver) (SEXP sexp_theta, SEXP sexp_inits, SEXP sexp_lhs, SEXP sexp_time, SEXP sexp_evid,SEXP sexp_dose, SEXP sexp_pcov, SEXP sexp_cov, SEXP sexp_locf, SEXP sexp_atol, SEXP sexp_rtol, SEXP sexp_hmin, SEXP sexp_hmax, SEXP sexp_h0, SEXP sexp_mxordn, SEXP sexp_mxords, SEXP sexp_mx,SEXP sexp_stiff, SEXP sexp_transit_abs, SEXP sexp_object, SEXP sexp_extra_args, SEXP sexp_matrix, SEXP sexp_add_cov);
typedef void (*RxODE_assign_ptr)(SEXP);
typedef void (*RxODE_ode_solver_old_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc);
typedef double (*RxODE_solveLinB)(rx_solve *rx, unsigned int id, double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);
typedef double (*RxODE_sum_prod)(double *input, int n);

double solveLinB(rx_solve *rx, unsigned int id, double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag){
  static RxODE_solveLinB fun = NULL;
  if (fun == NULL) fun = (RxODE_solveLinB) R_GetCCallable("RxODE","RxODE_solveLinB");
  return fun(rx, id, t, linCmt, diff1, diff2, A, alpha, B, beta, C, gamma, ka, tlag);
}

void _assign_ptr(SEXP x){
  static RxODE_assign_ptr fun = NULL;
  if (fun == NULL) fun = (RxODE_assign_ptr) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
  fun(x);
} 

double _sum1(double *input, int n){
  static RxODE_sum_prod fun = NULL;
  if (fun == NULL) fun = (RxODE_sum_prod) R_GetCCallable("RxODE","RxODE_sum");
  return fun(input, n);
}

double _prod1(double *input, int n){
  static RxODE_sum_prod fun = NULL;
  if (fun == NULL) fun = (RxODE_sum_prod) R_GetCCallable("RxODE","RxODE_prod");
  return fun(input, n);
}

typedef void (*_rxRmModelLibType)(const char *inp);
void _rxRmModelLib(const char *inp){
  static _rxRmModelLibType fun = NULL;
  if (fun == NULL) fun = (_rxRmModelLibType) R_GetCCallable("RxODE","rxRmModelLib");
  fun(inp);
}

typedef SEXP (*_rxGetModelLibType)(const char *s);
SEXP _rxGetModelLib(const char *inp){
  static _rxGetModelLibType fun = NULL;
  if (fun == NULL) fun = (_rxGetModelLibType) R_GetCCallable("RxODE","rxGetModelLib");
  return fun(inp);
}

void _old_c(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc){
  static RxODE_ode_solver_old_c fun = NULL;
  if (fun == NULL) fun = (RxODE_ode_solver_old_c) R_GetCCallable("RxODE","rxSolveOldC");
  return fun(neq, theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, 
	     transit_abs, nlhs, lhs, rc);
}

double Rx_pow_di(double x, double y){
  static RxODE_fn2i fun = NULL;
  if (fun == NULL) fun = (RxODE_fn2i) R_GetCCallable("RxODE","RxODE_pow_di");
  return fun(x, y);
}

double Rx_pow(double x, double y){
  static RxODE_fn2 fun = NULL;
  if (fun == NULL) fun = (RxODE_fn2) R_GetCCallable("RxODE","RxODE_pow");
  return fun(x, y);
}

double sign_exp(double x, double y){
  static RxODE_fn2 fun = NULL;
  if (fun == NULL) fun = (RxODE_fn2) R_GetCCallable("RxODE","RxODE_sign_exp");
  return fun(x, y);
}

double _as_zero(double x){
  static RxODE_fn fun = NULL;
  if (fun == NULL) fun = (RxODE_fn) R_GetCCallable("RxODE","RxODE_as_zero");
  return fun(x);
}

double _safe_log(double x){
  static RxODE_fn fun = NULL;
  if (fun == NULL) fun = (RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_log");
  return fun(x);
}

double safe_zero(double x){
  static RxODE_fn fun = NULL;
  if (fun == NULL) fun = (RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_zero");
  return fun(x);
}

double abs_log(double x){
  static RxODE_fn fun = NULL;
  if (fun == NULL) fun = (RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log");
  return fun(x);
}

double abs_log1p(double x){
  static RxODE_fn fun = NULL;
  if (fun == NULL) fun = (RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log1p");
  return fun(x);
}

double factorial(double x){
  static RxODE_fn fun = NULL;
  if (fun == NULL) fun = (RxODE_fn) R_GetCCallable("RxODE","RxODE_factorial");
  return fun(x);
}

double _transit4P(double t, rx_solve *rx, unsigned int id, double n, double mtt, double bio){
  static RxODE_transit4P fun = NULL;
  if (fun == NULL) fun = (RxODE_transit4P) R_GetCCallable("RxODE","RxODE_transit4P");
  return fun(t, rx, id, n, mtt, bio);
}

double _transit3P(double t, rx_solve *rx, unsigned int id, double n, double mtt){
  static RxODE_transit3P fun = NULL;
  if (fun == NULL) fun = (RxODE_transit3P) R_GetCCallable("RxODE","RxODE_transit3P");
  return fun(t, rx, id, n, mtt);
}

double podo(rx_solve *rx, unsigned int id){
  static RxODE_val fun = NULL;
  if (fun == NULL) fun = (RxODE_val) R_GetCCallable("RxODE","RxODE_podoP");
  return fun(rx, id);
}

double tlast(rx_solve *rx, unsigned int id){
  static RxODE_val fun = NULL;
  if (fun == NULL) fun = (RxODE_val) R_GetCCallable("RxODE","RxODE_tlastP");
  return fun(rx, id);
}

void _dadt_counter_inc(rx_solve *rx, unsigned int id){
  static RxODE_inc fun = NULL;
  if (fun == NULL) fun = (RxODE_inc) R_GetCCallable("RxODE","RxODE_dadt_counter_incP");
  fun(rx, id);
}

void _jac_counter_inc(rx_solve *rx, unsigned int id){
  static RxODE_inc fun = NULL;
  if (fun == NULL) fun = (RxODE_inc) R_GetCCallable("RxODE","RxODE_jac_counter_incP");
  fun(rx, id);
}

long _dadt_counter_val(rx_solve *rx, unsigned int id){
  static RxODE_cnt fun = NULL;
  if (fun == NULL) fun = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_dadt_counter_valP");
  return fun(rx, id);
}

long _jac_counter_val(rx_solve *rx, unsigned int id){
  static RxODE_cnt fun = NULL;
  if (fun == NULL) fun = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_jac_counter_valP");
  return fun(rx, id);
}

void _update_par_ptr(double t, rx_solve *rx, unsigned int id){
  static RxODE_update_par_ptr fun = NULL;
  if (fun == NULL) fun = (RxODE_update_par_ptr) R_GetCCallable("RxODE","RxODE_update_par_ptrP");
  return fun(t, rx, id);
}

double _par_ptr(int val, rx_solve *rx, unsigned int id){
  static RxODE_vec fun = NULL ;
  if (fun == NULL) fun = (RxODE_vec) R_GetCCallable("RxODE","RxODE_par_ptrP");
  return fun(val, rx, id);
}

double _InfusionRate(int val, rx_solve *rx, unsigned int id){
  static RxODE_vec fun = NULL ;
  if (fun == NULL) fun = (RxODE_vec) R_GetCCallable("RxODE","RxODE_InfusionRateP");
  return fun(val, rx, id);
}

typedef  SEXP (*_rx_asgn) (SEXP objectSEXP);
SEXP _RxODE_rxAssignPtr(SEXP objectSEXP){
  static _rx_asgn fun = NULL;
  if (fun==NULL) fun = (_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxAssignPtr");
  return fun(objectSEXP);
}

double _sum(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = _sum1(p, n);
  Free(p);
  return s;
}

double _prod(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = _prod1(p, n);
  Free(p);
  return s;
}

extern double _sign(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double s = 1;
  for (unsigned int i = 0; i < n; i++){
    s = sign(va_arg(valist, double))*s;
    if (s == 0){
      break;
    }
  }
  va_end(valist);
  return s;
}

SEXP _rxModels;
SEXP __MODEL_VARS__0();
extern SEXP __MODEL_VARS__(){
  SEXP _mv = _rxGetModelLib(__ODE_SOLVER_PTR_STR__);
  if (isNull(_mv)){
    _mv = __MODEL_VARS__0();
    _assign_ptr(_mv);
    return _mv;
  } else {
    return _mv;
  }
}

rx_solve *_solveData = NULL;
extern SEXP __ODE_SOLVER_XPTR__ ();

extern void __ODE_SOLVER_SOLVEDATA__ (rx_solve *solve){
  _solveData = solve;
}

extern rx_solve *__ODE_SOLVER_GET_SOLVEDATA__(){
  return _solveData;
}

extern void __ODE_SOLVER_PTR__();
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
  _assign_ptr(__MODEL_VARS__());
  _old_c(neq, theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc);
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
  // Register the outside functions
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_STR__,       (DL_FUNC) __ODE_SOLVER__);
  R_RegisterCCallable(__LIB_STR__,"__INIS__", (DL_FUNC) __INIS__);
  R_RegisterCCallable(__LIB_STR__,"__INIS__", (DL_FUNC) __INIS__);
  R_RegisterCCallable(__LIB_STR__,"__DYDT__", (DL_FUNC) __DYDT__);
  R_RegisterCCallable(__LIB_STR__,"__CALC_LHS__", (DL_FUNC) __CALC_LHS__);
  R_RegisterCCallable(__LIB_STR__,"__CALC_JAC__", (DL_FUNC) __CALC_JAC__);
  R_RegisterCCallable(__LIB_STR__,"__DYDT_LSODA__", (DL_FUNC) __DYDT_LSODA__);
  R_RegisterCCallable(__LIB_STR__,"__CALC_JAC_LSODA__", (DL_FUNC) __CALC_JAC_LSODA__);
  R_RegisterCCallable(__LIB_STR__,"__ODE_SOLVER_SOLVEDATA__", (DL_FUNC) __ODE_SOLVER_SOLVEDATA__);
  R_RegisterCCallable(__LIB_STR__,"__ODE_SOLVER_GET_SOLVEDATA__", (DL_FUNC) __ODE_SOLVER_GET_SOLVEDATA__);
  R_RegisterCCallable(__LIB_STR__,"__DYDT_LIBLSODA__", (DL_FUNC) __DYDT_LIBLSODA__);
  
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
  Rprintf("Unloading __ODE_SOLVER__ ....");
  SEXP _mv = _rxGetModelLib(__ODE_SOLVER_PTR_STR__);
  if (!isNull(_mv)){
    _rxRmModelLib(__ODE_SOLVER_PTR_STR__);
  }
  Rprintf("done\\n");
}
