#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.

// These are more precise sum algorithms.
// They are adapted from pseudo-code on wikipedia
// https://en.wikipedia.org/wiki/Kahan_summation_algorithm

//Not used: https://sourceforge.net/p/gmat/git/ci/264a12acad195e6a2467cfdc68abdcee801f73fc/tree/prototype/OptimalControl/ThirdParty/Intlab_V6/accsumdot/FastAccSum.m

extern double RxODE_DoubleSum(double *input, int n){
  double sum = input[0];
  for (int i = 1; i < n; i++){
    sum += input[i];
  }
  return sum;
}

extern double RxODE_KahanSum(double *input, int len){
  volatile double sum = 0.0;
  volatile double y;
  volatile double t, c = 0.0; // A running compensation for lost low-order bits.
  int i;
  for (i = 0; i < len; i++){
    y = input[i] - c; 
    t = sum + y;
    c = (t - sum) - y;       // (t - sum) cancels the high-order part of y; subtracting y recovers negative (low part of y)
    sum = t;                 // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
  }
  return sum;
}

SEXP _rxKahanSum(SEXP input){
  int len = length(input);
  double *dinput = REAL(input);
  SEXP rets = PROTECT(allocVector(REALSXP,1));
  REAL(rets)[0] = RxODE_KahanSum(dinput, len);
  UNPROTECT(1);
  return rets;
}

extern double RxODE_NeumaierSum(double *input, int len){
  double sum = input[0];
  volatile double t,  c = 0.0; // A running compensation for lost low-order bits.
  int i;
  for (i = 1; i < len; i++){
    t = sum + input[i];
    if (fabs(sum) >= fabs(input[i])){
      c += (sum - t) + input[i]; // If sum is bigger, low-order digits of input[i] are lost.
    } else {
      c += (input[i] - t) + sum; // Else low-order digits of sum are lost
    }
    sum = t;
  }
  return sum + c; // Correction only applied once in the very end
}

SEXP _rxNeumaierSum(SEXP input){
  int len = length(input);
  double *dinput = REAL(input);
  SEXP rets = PROTECT(allocVector(REALSXP,1));
  REAL(rets)[0] = RxODE_NeumaierSum(dinput, len);
  UNPROTECT(1);
  return rets;
}

#define NUM_PARTIALS  32  /* initial partials array size, on stack */

extern double RxODE_Python_fsum (double *iterable, unsigned int iterable_len){
  // See http://code.activestate.com/recipes/393090-binary-floating-point-summation-accurate-to-full-p/
  // Also https://github.com/python/cpython/blob/a0ce375e10b50f7606cb86b072fed7d8cd574fe7/Modules/mathmodule.c
  // Mostly the same as python's math.fsum
  double x, y, t;
  double xsave, special_sum = 0.0, inf_sum = 0.0, sum = 0.0;
  volatile double hi, yr, lo;
  unsigned int ix, i, j, n = 0, m = NUM_PARTIALS;
  double *p = Calloc(NUM_PARTIALS, double);
  // for x in input
  for (ix = 0; ix < iterable_len; ix++){
    x = iterable[ix];
    xsave = x;
    for (i = j = 0; j < n; j++) {
      y = p[j];
      if (fabs(x) < fabs(y)) {
        t = x; x = y; y = t;
      }
      hi = x + y;
      yr = hi - x;
      lo = y - yr;
      if (lo != 0.0)
        p[i++] = lo;
      x = hi;
    }
    
    n = i; 
    if (x != 0.0) {
      if (!R_FINITE(x)) {
	/* a nonfinite x could arise either as
	   a result of intermediate overflow, or
	   as a result of a nan or inf in the
	   summands */
	if (R_FINITE(xsave) || ISNAN(xsave)) {
	  Free(p);
	  error("intermediate overflow in fsum");
	} else {
	  inf_sum += xsave;
        }
	special_sum += xsave;
	/* reset partials */
	n = 0;
      } else {
	if (n >= m){
          //&& _fsum_realloc(&p, n, ps, &m)
          // Doubles the size of array.
          m += m;
          p = Realloc(p, m, double);
	}
	p[n++] = x;
      }
    }
  }
  if (special_sum != 0.0) {
    if (ISNAN(inf_sum)){
      Free(p);
      error("-inf + inf in fsum");
    }
    sum = special_sum;
    Free(p);
    return sum;
  }

  hi = 0.0;
  j = n;
  if (n > 0) {
    hi = p[--n];
    /* sum_exact(ps, hi) from the top, stop when the sum becomes
       inexact. */
    while (n > 0) {
      x = hi;
      y = p[--n];
      if (fabs(y) >= fabs(x)){
	Rprintf("Partial Sums:\n");
	for (i = 0; i < j; i++){
	  Rprintf("p[%d] = %f\n",i,p[i]);
	}
	Rprintf("Assertion Error:\n");
	Rprintf("fabs(y) >= fabs(x) or %f >= %f\n",fabs(y),fabs(x));
	Free(p);
	error("Error in parital sums.");
      }
      hi = x + y;
      yr = hi - x;
      lo = y - yr;
      if (lo != 0.0)
        break;
    }
    /* Make half-even rounding work across multiple partials.
       Needed so that sum([1e-16, 1, 1e16]) will round-up the last
       digit to two instead of down to zero (the 1e-16 makes the 1
       slightly closer to two).  With a potential 1 ULP rounding
       error fixed-up, math.fsum() can guarantee commutativity. */
    if (n > 0 && ((lo < 0.0 && p[n-1] < 0.0) ||
                  (lo > 0.0 && p[n-1] > 0.0))) {
      y  = lo * 2.0;
      x  = hi + y;
      yr = x - hi;
      if (y == yr)
        hi = x;
    }
  }
  sum = hi;
  Free(p);
  return sum;
}

SEXP _rxPythonSum(SEXP input){
  int len = length(input);
  double *dinput = REAL(input);
  SEXP rets = PROTECT(allocVector(REALSXP,1));
  REAL(rets)[0] = RxODE_Python_fsum(dinput, len);
  UNPROTECT(1);
  return rets;
}

unsigned int RxODE_sum_type = 1;
extern double RxODE_sum (double *input, unsigned int n){
  /* switch (RxODE_sum_type){ */
  /* case 0: */
  /*   return RxODE_DoubleSum(input, n); */
  /* case 1: */
    return RxODE_Python_fsum(input, n);
  /*   break; */
  /* } */
  /* error("Unknown sum type."); */
  /* return 0; */
}

SEXP _rxSum(SEXP input){
  int len = length(input);
  double *dinput = REAL(input);
  SEXP rets = PROTECT(allocVector(REALSXP,1));
  REAL(rets)[0] = RxODE_sum(dinput, len);
  UNPROTECT(1);
  return rets;
}

extern double RxODE_sumV(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = RxODE_sum(p, n);
  Free(p);
  return s;
}

