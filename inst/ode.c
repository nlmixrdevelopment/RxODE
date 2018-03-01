#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#define JAC_Rprintf Rprintf
#define JAC0_Rprintf if ( (&_solveData->subjects[_cSub])->jac_counter == 0) Rprintf
#define ODE_Rprintf Rprintf
#define ODE0_Rprintf if ( (&_solveData->subjects[_cSub])->dadt_counter == 0) Rprintf
#define LHS_Rprintf Rprintf
#define R_pow(a, b) (((a) == 0 && (b) <= 0) ? R_pow(DOUBLE_EPS, b) : R_pow(a, b))
#define R_pow_di(a, b) (((a) == 0 && (b) <= 0) ? R_pow_di(DOUBLE_EPS, b) : R_pow_di(a, b))
#define Rx_pow(a, b) (((a) == 0 && (b) <= 0) ? R_pow(DOUBLE_EPS, b) : R_pow(a, b))
#define Rx_pow_di(a, b) (((a) == 0 && (b) <= 0) ? R_pow_di(DOUBLE_EPS, b) : R_pow_di(a, b))


// Types for par pointers.r
typedef double (*RxODE_fn) (double x);
typedef double (*RxODE_fn2) (double x, double y);
typedef double (*RxODE_fn2i) (double x, int i);
typedef int (*RxODE_fn0i) ();
typedef double (*RxODE_vec) (int val, rx_solve *rx, unsigned int id);
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


RxODE_fn2 sign_exp = NULL;

RxODE_fn _as_zero = NULL;
RxODE_fn _safe_log = NULL;
RxODE_fn safe_zero = NULL;
RxODE_fn abs_log = NULL;
RxODE_fn abs_log1p = NULL;
RxODE_fn factorial = NULL;

RxODE_fn0i _ptrid=NULL;

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

double _prod(double *input, double *p, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  return _prodPS(input, p, n, type);
}

double _sign(unsigned int n, ...){
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

double _max(unsigned int n, ...){
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

double _min(unsigned int n, ...){
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

/* Authors: Robert Gentleman and Ross Ihaka and The R Core Team */
/* Taken directly from https://github.com/wch/r-source/blob/922777f2a0363fd6fe07e926971547dd8315fc24/src/library/stats/src/approx.c*/
/* Changed as follows:
   - Different Name
   - Use RxODE structure
   - Make inline
*/
inline double rx_approxP(double v, double *x, double *y, int n,
                         rx_solving_options *Meth, rx_solving_options_ind *id){
  /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
  int i, j, ij;

  if(!n) return R_NaN;

  i = 0;
  j = n - 1;

  /* handle out-of-domain points */
  if(v < x[i]) return id->ylow;
  if(v > x[j]) return id->yhigh;

  /* find the correct interval by bisection */
  while(i < j - 1) { /* x[i] <= v <= x[j] */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    if(v < x[ij]) j = ij; else i = ij;
    /* still i < j */
  }
  /* provably have i == j-1 */

  /* interpolation */

  if(v == x[j]) return y[j];
  if(v == x[i]) return y[i];
  /* impossible: if(x[j] == x[i]) return y[i]; */

  if(Meth->kind == 1) /* linear */
    return y[i] + (y[j] - y[i]) * ((v - x[i])/(x[j] - x[i]));
  else /* 2 : constant */
    return (Meth->f1 != 0.0 ? y[i] * Meth->f1 : 0.0)
      + (Meth->f2 != 0.0 ? y[j] * Meth->f2 : 0.0);
}/* approx1() */

/* End approx from R */


inline void _update_par_ptr(double t, unsigned int id){
  rx_solving_options_ind *ind;
  ind = (&_solveData->subjects[id]);
  rx_solving_options *op = _solveData->op;
  if (op->neq > 0){
    // Update all covariate parameters
    int k;
    double *par_ptr = ind->par_ptr;
    double *all_times = ind->all_times;
    double *cov_ptr = ind->cov_ptr;
    int ncov = op->ncov;
    if (op->do_par_cov){
      for (k = 0; k < ncov; k++){
        if (op->par_cov[k]){
          // Use the same methodology as approxfun.
          // There is some rumor the C function may go away...
          ind->ylow = cov_ptr[ind->n_all_times*k];
          ind->yhigh = cov_ptr[ind->n_all_times*k+ind->n_all_times-1];
          par_ptr[op->par_cov[k]-1] = rx_approxP(t, all_times, cov_ptr+ind->n_all_times*k, ind->n_all_times, op, ind);
        }
        /* if (global_debug){ */
        /*   REprintf("par_ptr[%d] (cov %d/%d) = %f\\n",op->par_cov[k]-1, k,ncov,cov_ptr[op->par_cov[k]-1]); */
        /* } */
      }
    }
  }
}

inline double _transit4P(double t, unsigned int id, double n, double mtt, double bio){
  double ktr = (n+1)/mtt;
  double lktr = log(n+1)-log(mtt);
  return exp(log(bio*(_solveData->subjects[id].podo))+lktr+n*(lktr+log(t))-ktr*t-lgamma1p(n));
}

inline double _transit3P(double t, unsigned int id, double n, double mtt){
  double ktr = (n+1)/mtt;
  double lktr = log(n+1)-log(mtt);
  return exp(log(_solveData->subjects[id].podo)+lktr+n*(lktr+log(t))-ktr*t-lgamma1p(n));
}

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
  if (_ptrid() != __TIMEID__ ){ _assign_ptr(__MODEL_VARS__());}
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
  sign_exp = (RxODE_fn2) R_GetCCallable("RxODE","RxODE_sign_exp");
  _as_zero = (RxODE_fn) R_GetCCallable("RxODE","RxODE_as_zero");
  _safe_log=(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_log");
  safe_zero=(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_zero");
  abs_log = (RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log");
  abs_log1p=(RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log1p");
  factorial = (RxODE_fn) R_GetCCallable("RxODE","RxODE_factorial");
  _RxODE_rxAssignPtr=(_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxAssignPtr");
  _rxIsCurrentC = (_rxIsCurrentC_type)R_GetCCallable("RxODE","rxIsCurrentC");
  _sumPS  = (_rxSumType) R_GetCCallable("PreciseSums","PreciseSums_sum_r");
  _prodPS = (_rxProdType) R_GetCCallable("PreciseSums","PreciseSums_prod_r");
  _prodType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_prod_get");
  _sumType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_sum_get");
  _ptrid=(RxODE_fn0i)R_GetCCallable("RxODE", "RxODE_current_fn_pointer_id");
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
