#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.

extern double RxODE_sum (double *input, int n);

extern double RxODE_pairwise_add_DOUBLE(double *a, int n);

extern double RxODE_safe_log(double x);

int RxODE_prod_type = 3;

extern void RxODE_prod_set(int i){
  RxODE_prod_type = i;
}

extern int RxODE_prod_get(){
  return RxODE_prod_type;
}


SEXP _rxSetProd(SEXP input){
  RxODE_prod_type = (int) INTEGER(input)[0];
  return R_NilValue;
}

extern double RxODE_prod_ld(double *input, int n){
  long double p = 1;
  for  (int i = 0; i < n; i++){
    if (input[i] == 0){
      return 0.0; 
    }
    p *= input[i];
  }
  return (double)p;
}

extern double RxODE_prod_d(double *input, int n){
  double p = 1;
  for  (int i = 0; i < n; i++){
    if (input[i] == 0){
      return 0.0; 
    }
    p *= input[i];
  }
  return p;
}

extern double RxODE_prod_logify(double *input, int n){
  double *p = Calloc(n,double);
  double s = 1.0;
  for (int i = 0; i < n; i++){
    if (input[i] == 0){
      Free(p);
      return 0.0;
    }
    s = sign(input[i])*s;
    p[i] = RxODE_safe_log(fabs(input[i]));
  }
  s = exp(RxODE_pairwise_add_DOUBLE(p, n))*s;
  Free(p);
  return s;
}

extern double RxODE_prod(double *input, int n){
  switch (RxODE_prod_type){
  case 1: // long double multiply, then convert back.
    return RxODE_prod_ld(input, n);
    break;
  case 2: // simple double multiply
    return RxODE_prod_d(input, n);
    break;
  case 3: // logify
    return RxODE_prod_logify(input, n);
  }
  return 0.0;
}



SEXP _rxProd(SEXP input){
  int len = length(input);
  double *dinput = REAL(input);
  SEXP rets = PROTECT(allocVector(REALSXP,1));
  REAL(rets)[0] = RxODE_prod(dinput, len);
  UNPROTECT(1);
  return rets;
}

extern double RxODE_prodV(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = RxODE_prod(p, n);
  Free(p);
  return s;
}


extern double RxODE_signV(int n, ...){
  va_list valist;
  va_start(valist, n);
  double s = 1;
  for (int i = 0; i < n; i++){
    s = sign(va_arg(valist, double))*s;
    if (s == 0){
      break;
    }
  }
  va_end(valist);
  return s;
}
