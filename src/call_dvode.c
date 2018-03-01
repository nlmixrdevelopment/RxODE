#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "dop853.h"
#define max(a, b) ((a) > (b) ? (a) : (b))
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>
#include <PreciseSums.h>
#include "solve.h"

extern void calc_lhs(int cSub, double t, double *A, double *lhs);

extern double RxODE_as_zero(double x){
  if (fabs(x) < sqrt(DOUBLE_EPS)){
    return(0.0);
  } else {
    return(x);
  }
}

extern double RxODE_safe_log(double x){
  if (x <= 0){
    // Warning?
    return log(DOUBLE_EPS);
  } else {
    return log(x);
  }
}

extern double RxODE_safe_zero(double x){
  if (x == 0){
    // Warning?
    return DOUBLE_EPS;
  } else {
    return(x);
  }
}

extern double RxODE_pow(double x, double y){
  if (x == 0 && y <= 0){
    return R_pow(DOUBLE_EPS, y);
  } else {
    return R_pow(x, y);
  }
}
extern double RxODE_pow_di(double x, int i){
  if (x == 0 && i <= 0){
    return R_pow_di(DOUBLE_EPS, i);
  } else {
    return R_pow_di(x, i);
  }
}

extern double RxODE_sign_exp(double sgn, double x){
  if (sgn > 0.0){
    return(exp(x));
  } else if (sgn < 0.0){
    return(-exp(x));
  } else {
    return(0.0);
  }
}

extern double RxODE_abs_log(double x){
  if  (fabs(x) <= sqrt(DOUBLE_EPS)){
    return log(sqrt(DOUBLE_EPS));
  } else if (x > 0.0){
    return log(x);
  } else if (x < 0.0){
    return log(-x);
  } else {
    return 0.0;
  }
}

extern double RxODE_abs_log1p(double x){
  if (x + 1.0 > 0.0){
    return(log1p(x));
  } else if (x + 1.0 > 0.0){
    return(log1p(-x));
  } else {
    return 0.0;
  }
}

extern double RxODE_factorial(double x){
  return exp(lgamma1p(x));
}

//--------------------------------------------------------------------------

// These are now allocated via R structures in Rcpp.
extern void RxODE_ode_free(){
}

void RxODE_ode_alloc(){
}

char __mv[1000];
extern void RxODE_assign_fn_pointers_(const char *mv){
  sprintf(__mv, "%s", mv);
}

void rxAssignPtrC(SEXP obj);
int RxODE_current_fn_pointer_id_ = 0;
extern int RxODE_current_fn_pointer_id(){
  return RxODE_current_fn_pointer_id_;
}
extern void RxODE_assign_fn_pointers(SEXP mv){
  RxODE_current_fn_pointer_id_ = INTEGER(VECTOR_ELT(mv, 15))[0];
  rxAssignPtrC(mv);
}

SEXP rxModelVarsC(char *ptr);

extern SEXP RxODE_get_mv(){
  /* if (!rxIsC(__mv,"rxModelVars")){ */
  /*   error("RxODE C functions were not setup correctly."); */
  /* } */
  return rxModelVarsC(__mv);
}

/* extern void rxode_assign_rx(rx_solve *rx); */

extern void rxode_assign_rx(rx_solve *rx);
extern rx_solve *rxSingle(SEXP object, const int stiff,const int transit_abs,
			  const double atol, const double rtol, const int maxsteps,
			  const double hmin, const double hini, const int maxordn,
			  const int maxords, const int cores, const int ncov,
			  int *par_cov, int do_par_cov, 
			  int is_locf,
			  // Other single solve option
			  double hmax, double *par,
			  double *amt, double *solve, double *lhs,
			  int *evid, int *rc, double *cov,
			  int nTimes, double *all_times);

rx_solving_options_ind *rxOptionsIniEnsure(int mx);
rx_solving_options *getRxOp(rx_solve *rx);
rx_solve *getRxSolve_();
int *global_BadDose(unsigned int mx);
double *global_InfusionRate(unsigned int mx);



extern double RxODE_sum(double *input, int len){
  return PreciseSums_sum(input, len);
}

extern double RxODE_sumV(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = PreciseSums_sum(p, n);
  Free(p);
  return s;
}

extern double RxODE_sumV_r(double *p, long double *pld, int m, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  return PreciseSums_sum_r(p, n, pld, m, type);
}

extern double RxODE_prod(double *input, int len){
  return PreciseSums_prod(input, len);
}

extern double RxODE_prodV(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = PreciseSums_prod(p, n);
  Free(p);
  return s;
}

extern double RxODE_prodV_r(double *input, double *p, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  return PreciseSums_prod_r(input, p, n, type);
}
