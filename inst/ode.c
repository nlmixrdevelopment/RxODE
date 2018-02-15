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
#define R_pow Rx_pow
#define R_pow_di Rx_pow_di

// Types for par pointers.r
typedef void (*RxODE_update_par_ptr)(double t, rx_solve *rx, unsigned int id);
typedef double (*RxODE_transit3P)(double t, rx_solve *rx, unsigned int id, double n, double mtt);
typedef double (*RxODE_fn) (double x);
typedef double (*RxODE_fn2) (double x, double y);
typedef double (*RxODE_fn2i) (double x, int i);
typedef int (*RxODE_fn0i) ();
typedef double (*RxODE_transit4P)(double t, rx_solve *rx, unsigned int id, double n, double mtt, double bio);
typedef double (*RxODE_vec) (int val, rx_solve *rx, unsigned int id);
typedef long (*RxODE_cnt) (rx_solve *rx, unsigned int id);
typedef void (*RxODE_inc) (rx_solve *rx, unsigned int id);
typedef double (*RxODE_val) (rx_solve *rx, unsigned int id);
typedef void (*RxODE_assign_ptr)(SEXP);
typedef void (*RxODE_ode_solver_old_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc);
typedef double (*RxODE_solveLinB)(rx_solve *rx, unsigned int id, double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);

RxODE_solveLinB solveLinB    = NULL;
RxODE_assign_ptr _assign_ptr = NULL;

typedef void (*_rxRmModelLibType)(const char *inp);
_rxRmModelLibType _rxRmModelLib = NULL;

typedef SEXP (*_rxGetModelLibType)(const char *s);
_rxGetModelLibType _rxGetModelLib = NULL;

RxODE_ode_solver_old_c _old_c = NULL;

RxODE_fn2i Rx_pow_di = NULL;
RxODE_fn2 Rx_pow = NULL;
RxODE_fn2 sign_exp = NULL;

RxODE_fn _as_zero = NULL;
RxODE_fn _safe_log = NULL;
RxODE_fn safe_zero = NULL;
RxODE_fn abs_log = NULL;
RxODE_fn abs_log1p = NULL;
RxODE_fn factorial = NULL;
RxODE_transit4P _transit4P = NULL;
RxODE_transit3P _transit3P =NULL;
RxODE_val podo = NULL;
RxODE_val tlast = NULL;
RxODE_inc _dadt_counter_inc=NULL;
RxODE_inc _jac_counter_inc =NULL;
RxODE_cnt _dadt_counter_val = NULL;
RxODE_cnt _jac_counter_val = NULL;
RxODE_vec _par_ptr = NULL;
RxODE_update_par_ptr _update_par_ptr=NULL;
RxODE_vec _InfusionRate= NULL;

typedef  SEXP (*_rx_asgn) (SEXP objectSEXP);
_rx_asgn _RxODE_rxAssignPtr =NULL;

typedef int(*_rxIsCurrentC_type)(SEXP);
_rxIsCurrentC_type _rxIsCurrentC=NULL;

typedef double(*_rxSumType)(double *, int, double *, int, int);
_rxSumType _sumPS = NULL;

double _sum(double *input, double *pld, int m, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  return _sumPS(input, n, pld, m, type);
}

typedef double(*_rxProdType)(double*, double*, int, int);
_rxProdType _prodPS = NULL;

extern double _prod(double *input, double *p, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  return _prodPS(input, p, n, type);
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

extern double _max(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double mx = NA_REAL;
  double tmp = 0;
  if (n >= 1){
    mx = va_arg(valist, double);
    for (unsigned int i = 1; i < n; i++){
      tmp = va_arg(valist, double);
      if (tmp>mx) mx=tmp;
    }
    va_end(valist);
  }
  return mx;
}

extern double _min(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double mn = NA_REAL;
  double tmp = 0;
  if (n >= 1){
    mn = va_arg(valist, double);
    for (unsigned int i = 1; i < n; i++){
      tmp = va_arg(valist, double);
      if (tmp<mn) mn=tmp;
    }
    va_end(valist);
  }
  return mn;
}

rx_solve *_solveData = NULL;

RxODE_fn0i _prodType = NULL;
RxODE_fn0i _sumType = NULL;

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
  // Get C callables on load; Otherwise it isn't thread safe
  solveLinB = (RxODE_solveLinB) R_GetCCallable("RxODE","RxODE_solveLinB");
  _assign_ptr=(RxODE_assign_ptr) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
  _rxRmModelLib=(_rxRmModelLibType) R_GetCCallable("RxODE","rxRmModelLib");
  _rxGetModelLib=(_rxGetModelLibType) R_GetCCallable("RxODE","rxGetModelLib");
  _old_c = (RxODE_ode_solver_old_c) R_GetCCallable("RxODE","rxSolveOldC");
  Rx_pow_di=(RxODE_fn2i) R_GetCCallable("RxODE","RxODE_pow_di");
  Rx_pow = (RxODE_fn2) R_GetCCallable("RxODE","RxODE_pow");
  sign_exp = (RxODE_fn2) R_GetCCallable("RxODE","RxODE_sign_exp");
  _as_zero = (RxODE_fn) R_GetCCallable("RxODE","RxODE_as_zero");
  _safe_log=(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_log");
  safe_zero=(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_zero");
  abs_log = (RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log");
  abs_log1p=(RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log1p");
  factorial = (RxODE_fn) R_GetCCallable("RxODE","RxODE_factorial");
  _transit4P= (RxODE_transit4P) R_GetCCallable("RxODE","RxODE_transit4P");
  _transit3P=(RxODE_transit3P) R_GetCCallable("RxODE","RxODE_transit3P");
  podo = (RxODE_val) R_GetCCallable("RxODE","RxODE_podoP");
  tlast = (RxODE_val) R_GetCCallable("RxODE","RxODE_tlastP");
  _dadt_counter_inc=(RxODE_inc) R_GetCCallable("RxODE","RxODE_dadt_counter_incP");
  _jac_counter_inc=(RxODE_inc) R_GetCCallable("RxODE","RxODE_jac_counter_incP");
  _dadt_counter_val= (RxODE_cnt) R_GetCCallable("RxODE","RxODE_dadt_counter_valP");
  _jac_counter_val=(RxODE_cnt) R_GetCCallable("RxODE","RxODE_jac_counter_valP");
  _update_par_ptr=(RxODE_update_par_ptr) R_GetCCallable("RxODE","RxODE_update_par_ptrP");
  _par_ptr=(RxODE_vec) R_GetCCallable("RxODE","RxODE_par_ptrP");
  _InfusionRate=(RxODE_vec) R_GetCCallable("RxODE","RxODE_InfusionRateP");
  _RxODE_rxAssignPtr=(_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxAssignPtr");
  _rxIsCurrentC = (_rxIsCurrentC_type)R_GetCCallable("RxODE","rxIsCurrentC");
  _sumPS  = (_rxSumType) R_GetCCallable("PreciseSums","PreciseSums_sum_r");
  _prodPS = (_rxProdType) R_GetCCallable("PreciseSums","PreciseSums_prod_r");
  _prodType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_prod_get");
  _sumType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_sum_get");
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
  SEXP _mv = PROTECT(_rxGetModelLib(__MODEL_VARS_STR__));
  if (!isNull(_mv)){
    _rxRmModelLib(__MODEL_VARS_STR__);
  }
  UNPROTECT(1);
}
