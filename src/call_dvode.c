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
#include "../inst/include/RxODE.h"

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
  int cur = INTEGER(VECTOR_ELT(mv, 13))[0];
  if (RxODE_current_fn_pointer_id_ != cur){
    rxAssignPtrC(mv);
    RxODE_current_fn_pointer_id_ = cur;
  } 
}

SEXP rxModelVarsC(char *ptr);

extern SEXP RxODE_get_mv(){
  return rxModelVarsC(__mv);
}

/* extern void rxode_assign_rx(rx_solve *rx); */

extern void rxode_assign_rx(rx_solve *rx);

void rxOptionsIniEnsure(int mx);
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
  for (unsigned int i = (unsigned int)n; i--;){
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
  for (unsigned int i = (unsigned int)n; i--;){
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
  for (unsigned int i = (unsigned int)n; i--;){
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
  for (unsigned int i = (unsigned int)n; i--;){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  return PreciseSums_prod_r(input, p, n, type);
}
