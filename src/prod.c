#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.

extern long double RxODE_pairwise_add_ld2(long double *a, unsigned int n);

unsigned int RxODE_prod_type = 1;
extern void RxODE_prod_set(unsigned int i){
  RxODE_prod_type = i;
}

extern unsigned int RxODE_prod_get(){
  return RxODE_prod_type;
}


SEXP _rxSetProd(SEXP input){
  RxODE_prod_type = (unsigned int) INTEGER(input)[0];
  return R_NilValue;
}

extern double RxODE_prod_ld(double *input, unsigned int n){
  long double p = 1;
  for  (unsigned int i = 0; i < n; i++){
    if (input[i] == 0){
      return 0.0; 
    }
    p *= input[i];
  }
  return (double)p;
}

extern double RxODE_prod_ldl(double *input, unsigned int n){
  long double s = 1.0;
  long double *p = Calloc(n,long double);
  p[0] = (long double)input[0];
  for  (unsigned int i = 0; i < n; i++){
    if (input[i] == 0){
      Free(p);
      return 0.0;
    }
    s =(long double)(sign(input[i]));
    p[i] = logl(fabsl((long double)input[i]));
  }
  s = (double)(s*expl(RxODE_pairwise_add_ld2(p, n)));
  Free(p);
  return (double)s;
}

extern double RxODE_prod_d(double *input, unsigned int n){
  double p = 1;
  for  (unsigned int i = 0; i < n; i++){
    if (input[i] == 0){
      return 0.0; 
    }
    p *= input[i];
  }
  return p;
}

extern double RxODE_prod(double *input, unsigned int n){
  switch (RxODE_prod_type){
  case 1: // long double multiply, then convert back.
    return RxODE_prod_ld(input, n);
    break;
  case 2: // simple double multiply
    return RxODE_prod_d(input, n);
    break;
  case 3: // long double exp(log add)
    return RxODE_prod_ldl(input, n);
    break;
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

extern double RxODE_prodV(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = RxODE_prod(p, n);
  Free(p);
  return s;
}


extern double RxODE_signV(unsigned int n, ...){
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
