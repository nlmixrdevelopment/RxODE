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
typedef double (*RxODE_transit4P)(double t, rx_solve *rx, unsigned int id, double n, double mtt, double bio);
typedef double (*RxODE_vec) (int val, rx_solve *rx, unsigned int id);
typedef long (*RxODE_cnt) (rx_solve *rx, unsigned int id);
typedef void (*RxODE_inc) (rx_solve *rx, unsigned int id);
typedef double (*RxODE_val) (rx_solve *rx, unsigned int id);
typedef void (*RxODE_assign_ptr)(SEXP);
typedef void (*RxODE_ode_solver_old_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc);
typedef double (*RxODE_solveLinB)(rx_solve *rx, unsigned int id, double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);


RxODE_solveLinB _RxODE_solveLinB = NULL;
double solveLinB(rx_solve *rx, unsigned int id, double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag){
  return _RxODE_solveLinB(rx, id, t, linCmt, diff1, diff2, A, alpha, B, beta, C, gamma, ka, tlag);
}

RxODE_assign_ptr _RxODE_assign_ptr = NULL;
void _assign_ptr(SEXP x){
  _RxODE_assign_ptr(x);
}

typedef void (*_rxRmModelLibType)(const char *inp);
_rxRmModelLibType _rxRmModelLibType_ = NULL;
void _rxRmModelLib(const char *inp){
  _rxRmModelLibType_(inp);
}

typedef SEXP (*_rxGetModelLibType)(const char *s);
_rxGetModelLibType _rxGetModelLibType_ = NULL;
SEXP _rxGetModelLib(const char *inp){
  return _rxGetModelLibType_(inp);
}

RxODE_ode_solver_old_c _RxODE_ode_solver_old_c = NULL;
void _old_c(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc){
  return _RxODE_ode_solver_old_c(neq, theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, 
	     transit_abs, nlhs, lhs, rc);
}

RxODE_fn2i _Rx_pow_di = NULL;
double Rx_pow_di(double x, double y){
  return _Rx_pow_di(x, y);
}

RxODE_fn2 _Rx_pow = NULL;
double Rx_pow(double x, double y){
  return _Rx_pow(x, y);
}

RxODE_fn2 _sign_exp = NULL;
double sign_exp(double x, double y){
  return _sign_exp(x, y);
}

RxODE_fn _as_zero_ = NULL;
double _as_zero(double x){
  return _as_zero_(x);
}

RxODE_fn _safe_log_ = NULL;
double _safe_log(double x){
  return _safe_log_(x);
}

RxODE_fn _safe_zero_=NULL;
double safe_zero(double x){
  return _safe_zero_(x);
}

RxODE_fn _abs_log_=NULL;
double abs_log(double x){
  return _abs_log_(x);
}

RxODE_fn _abs_log1p_=NULL;
double abs_log1p(double x){
  return _abs_log1p_(x);
}

RxODE_fn _factorial_=NULL;
double factorial(double x){
  return _factorial_(x);
}

RxODE_transit4P _transit4P_=NULL;
double _transit4P(double t, rx_solve *rx, unsigned int id, double n, double mtt, double bio){
  return _transit4P_(t, rx, id, n, mtt, bio);
}

RxODE_transit3P _transit3P_=NULL;
double _transit3P(double t, rx_solve *rx, unsigned int id, double n, double mtt){
  return _transit3P_(t, rx, id, n, mtt);
}

RxODE_val _podo_=NULL;
double podo(rx_solve *rx, unsigned int id){
  return _podo_(rx, id);
}

RxODE_val _tlast_=NULL;
double tlast(rx_solve *rx, unsigned int id){
  return _tlast_(rx, id);
}

RxODE_inc _dadt_counter_inc_=NULL;
void _dadt_counter_inc(rx_solve *rx, unsigned int id){
  _dadt_counter_inc_(rx, id);
}

RxODE_inc _jac_counter_inc_ = NULL;
void _jac_counter_inc(rx_solve *rx, unsigned int id){
  _jac_counter_inc_(rx, id);
}

RxODE_cnt _dadt_counter_val_=NULL;
long _dadt_counter_val(rx_solve *rx, unsigned int id){
  return _dadt_counter_val_(rx, id);
}

RxODE_cnt _jac_counter_val_=NULL;
long _jac_counter_val(rx_solve *rx, unsigned int id){
  return _jac_counter_val_(rx, id);
}

RxODE_update_par_ptr _update_par_ptr_=NULL;
void _update_par_ptr(double t, rx_solve *rx, unsigned int id){
  return _update_par_ptr_(t, rx, id);
}

RxODE_vec _par_ptr_=NULL;
double _par_ptr(int val, rx_solve *rx, unsigned int id){
  return _par_ptr_(val, rx, id);
}

RxODE_vec _InfusionRate_=NULL;
double _InfusionRate(int val, rx_solve *rx, unsigned int id){
  _InfusionRate_(val, rx, id);
}

typedef  SEXP (*_rx_asgn) (SEXP objectSEXP);
_rx_asgn _rx_asgn_ = NULL;
SEXP _RxODE_rxAssignPtr(SEXP objectSEXP){
  return _rx_asgn_(objectSEXP);
}

int(*_rxIsCurrentC_)(SEXP)= NULL;
int _rxIsCurrentC(SEXP obj){
  return _rxIsCurrentC_(obj);
}

double (*_sum_)(double *, int, long double *, int, int)=NULL;

double _sum(double *p, long double *pld, int m, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  return _sum_(p, n, pld, m, type);
}

static double (*_prod_)(double*, double*, int, int)=NULL;
extern double _prod(double *input, double *p, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  return _prod_(input, p, n, type);
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

int _prodType(){
  // Type 3 = Logify
  return 3;
}
int _sumType(){
  // Type 1 = PairwiseSum
  return 1;
}

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
  // Get calllables on load (So it isn't called in threads).
  _RxODE_solveLinB = (RxODE_solveLinB) R_GetCCallable("RxODE","RxODE_solveLinB");
  _RxODE_assign_ptr= (RxODE_assign_ptr) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
  _rxRmModelLibType_=(_rxRmModelLibType) R_GetCCallable("RxODE","rxRmModelLib");
  _rxGetModelLibType_=(_rxGetModelLibType) R_GetCCallable("RxODE","rxGetModelLib");
  _RxODE_ode_solver_old_c=(RxODE_ode_solver_old_c) R_GetCCallable("RxODE","rxSolveOldC");
  _Rx_pow_di=(RxODE_fn2i) R_GetCCallable("RxODE","RxODE_pow_di");
  _Rx_pow=(RxODE_fn2) R_GetCCallable("RxODE","RxODE_pow");
  _sign_exp=(RxODE_fn2) R_GetCCallable("RxODE","RxODE_sign_exp");
  _as_zero_=(RxODE_fn) R_GetCCallable("RxODE","RxODE_as_zero");
  _safe_log_=(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_log");
  _safe_zero_ = (RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_zero");
  _abs_log_=(RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log");
  _abs_log1p_=(RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log1p");
  _factorial_=(RxODE_fn) R_GetCCallable("RxODE","RxODE_factorial");
  _transit4P_=(RxODE_transit4P) R_GetCCallable("RxODE","RxODE_transit4P");
  _transit3P_= (RxODE_transit3P) R_GetCCallable("RxODE","RxODE_transit3P");
  _podo_=(RxODE_val) R_GetCCallable("RxODE","RxODE_podoP");
  _tlast_=(RxODE_val) R_GetCCallable("RxODE","RxODE_tlastP");
  _dadt_counter_inc_=(RxODE_inc) R_GetCCallable("RxODE","RxODE_dadt_counter_incP");
  _jac_counter_inc_=(RxODE_inc) R_GetCCallable("RxODE","RxODE_jac_counter_incP");
  _dadt_counter_val_= (RxODE_cnt) R_GetCCallable("RxODE","RxODE_dadt_counter_valP");
  _jac_counter_val_=(RxODE_cnt) R_GetCCallable("RxODE","RxODE_jac_counter_valP");
  _update_par_ptr_=(RxODE_update_par_ptr) R_GetCCallable("RxODE","RxODE_update_par_ptrP");
  _par_ptr_=(RxODE_vec) R_GetCCallable("RxODE","RxODE_par_ptrP");
  _InfusionRate_=(RxODE_vec) R_GetCCallable("RxODE","RxODE_InfusionRateP");
  _rx_asgn_=(_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxAssignPtr");
  _rxIsCurrentC_=(int(*)(SEXP))R_GetCCallable("RxODE","rxIsCurrentC");
  _sum_= (double(*)(double *, int, long double *, int, int)) R_GetCCallable("PreciseSums","PreciseSums_sum_r");
  _prod_=(double(*)(double*, double*, int, int)) R_GetCCallable("PreciseSums","PreciseSums_prod_r");
}

void __R_UNLOAD__ (DllInfo *info){
  // Free resources required for single subject solve.
  SEXP _mv = PROTECT(_rxGetModelLib(__MODEL_VARS_STR__));
  if (!isNull(_mv)){
    _rxRmModelLib(__MODEL_VARS_STR__);
  }
  UNPROTECT(1);
}
