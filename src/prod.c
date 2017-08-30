#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.

extern double RxODE_sum (double *input, unsigned int n);

extern double RxODE_prod(double *input, unsigned int n){
  double s = 1;
  double *p = Calloc(n, double);
  unsigned int i;
  for  (i = 0; i < n; i++){
    s = sign(input[i])*s;
    if (s == 0){
      // Return 0; 0*
      Free(p);
      return 0; 
    }
    p[i] = log(fabs(input[i]));
  }
  s = s*exp(RxODE_sum(p, n));
  Free(p);
  return s;
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
