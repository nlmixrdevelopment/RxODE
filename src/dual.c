#include <R.h>
#include <Rversion.h>
#include <Rinternals.h>
#include "../inst/include/RxODE.h"

typedef struct dualN {
  double f;
  int n;
  double grad[4];
} dualN;

dualN sqrtD(dualN x) {
  dualN ret;
  ret.n = x.n;
  ret.f = sqrt(x.f);
  double gr = 0.5 / ret.f;
  for (int i = x.n; i--;) {
    ret.grad[i] = x.grad[i]*gr;
  }
  return ret;
}

dualN expD(dualN x) {
  dualN ret;
  ret.n = x.n;
  ret.f = exp(x.f);
  for (int i = x.n; i--;) {
    ret.grad[i] = x.grad[i]*ret.f;
  }
  return ret;
}


dualN add2(dualN x, dualN y) {
  dualN ret;
  ret.n = x.n;
  ret.f = x.f + y.f;
  for (int i = ret.n; i--;){
    ret.grad[i] = x.grad[i] + y.grad[i];
  }
  return ret;
}

dualN subtr2(dualN x, dualN y) {
  dualN ret;
  ret.n = x.n;
  ret.f = x.f - y.f;
  for (int i = ret.n; i--;){
    ret.grad[i] = x.grad[i] - y.grad[i];
  }
  return ret;
}

dualN prod2(dualN e1, dualN e2) {
  dualN ret;
  ret.n = e1.n;
  ret.f = e1.f * e2.f;
  for (int i = ret.n; i--;) {
    ret.grad[i] = e1.grad[i] * e2.f + e1.f * e2.grad[i];
  }
  return ret;
}

dualN prodd2(double e1, dualN e2) {
  dualN ret;
  ret.n = e2.n;
  ret.f = e2.f * e1;
  for (int i = ret.n; i--;) {
    ret.grad[i] = e2.grad[i]*e1;
  }
  return ret;
}

dualN iniD(double val, int which, int n){
  dualN ret;
  ret.n = n;
  ret.f = val;
  for (int i = ret.n; i--;) {
    ret.grad[i] = 0.0;
  }
  if (which >= 0) {
    ret.grad[which] = 1.0;
  }
  return ret;
}

void DtwoCmtKaRate(double *A, double *Alast,
		   double *t, double *b1, double *b2,
		   double *r1, double *r2,
		   double *ka,  double *k20, 
		   double *k23, double *k32) {
  dualN kad  = iniD(*ka,  0, 4);
  dualN k23d = iniD(*k23, 1, 4);
  dualN k32d = iniD(*k32, 2, 4);
  dualN k20d = iniD(*k20, 3, 4);
  //double E2 =  (*k20)+ (*k23);
  dualN E2   = add2(k20d, k32d);
  //double s = (*k23)+(*k32)+(*k20);
  dualN s    = add2(add2(k23d, k32d), k20d);
  //double beta  = 0.5*(s - sqrt(s*s - 4*(*k32)*(*k20)));
  dualN betad = prodd2(0.5, subtr2(s, sqrtD(subtr2(prod2(s,s), prodd2(4.0, prod2(k32d, k23d))))));
}
