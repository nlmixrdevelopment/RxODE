#pragma once
#ifndef __RxODE_H__
#define __RxODE_H__
#define isDose(evid) ((evid) == 3 || (evid) >= 100)
#define isObs(evid) ((evid) == 0 || (evid) == 2 || ((evid) >= 9 && (evid) <= 99))
#if defined(__cplusplus)
#include "RxODE_RcppExports.h"
extern "C" {
#endif
#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

typedef void (*t_dydt)(int *neq, double t, double *A, double *DADT);
typedef void (*t_calc_jac)(int *neq, double t, double *A, double *JAC, unsigned int __NROWPD__);
typedef void (*t_calc_lhs)(int cSub, double t, double *A, double *lhs);
typedef void (*t_update_inis)(int cSub, double *);
typedef void (*t_dydt_lsoda_dum)(int *neq, double *t, double *A, double *DADT);
typedef void (*t_jdum_lsoda)(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd);
typedef int (*t_dydt_liblsoda)(double t, double *y, double *ydot, void *data);
typedef void (*t_ode_current)();
typedef double (*t_F)(int _cSub,  int _cmt, double _amt, double t);
typedef double (*t_LAG)(int _cSub,  int _cmt, double t);
typedef double (*t_RATE)(int _cSub,  int _cmt, double _amt, double t);
typedef double (*t_DUR)(int _cSub,  int _cmt, double _amt, double t);

typedef void (*t_calc_mtime)(int cSub, double *mtime);

typedef struct sbuf {
  char *s;        /* curr print buffer */
  int sN;
  int o;                        /* offset of print buffer */
} sbuf;

typedef struct {
  // These options should not change based on an individual solve
  int badSolve;
  double ATOL;          //absolute error
  double RTOL;          //relative error
  double H0;
  double HMIN;
  int mxstep;
  int MXORDN;
  int MXORDS;
  //
  int do_transit_abs;
  int nlhs;
  int neq;
  int stiff;
  int ncov;
  char modNamePtr[1000];
  int *par_cov;
  double *inits;
  double *scale;
  int do_par_cov;
  // approx fun options
  double f1;
  double f2;
  int kind;
  int is_locf;
  int cores;
  int extraCmt;
  double hmax2; // Determined by diff
  double *rtol2;
  double *atol2;
  double *ssRtol;
  double *ssAtol;
  int nDisplayProgress;
  int ncoresRV;
  int isChol;
  int *svar;
  int nsvar;
  int abort;
  int minSS;
  int maxSS;
  int linLog;
  int strictSS;
  double infSSstep;
  int mxhnil;
  double hmxi;
} rx_solving_options;


typedef struct {
  int *slvr_counter;
  int *dadt_counter;
  int *jac_counter;
  double *InfusionRate;
  int *BadDose;
  int nBadDose;
  double HMAX; // Determined by diff
  double tlast;
  double podo;
  double *par_ptr;
  double *dose;
  double *ii;
  double *solve;
  double *mtime;
  double *solveSave;
  double *solveLast;
  double *solveLast2;
  double *lhs;
  int  *evid;
  int *rc;
  double *cov_ptr;
  int n_all_times;
  int nevid2;
  int ixds;
  int ndoses;
  double *all_times;
  int *ix;
  double *dv;
  int *idose;
  int *on;
  int idosen;
  int id;
  int sim;
  int idx;
  double ylow;
  double yhigh;
  double lambda;
  double yj;
  // Saved info
  int wh;
  int wh100;
  int cmt;
  int whI;
  int wh0;
  int doSS;
  int allCovWarn;
  int wrongSSDur;
  int timeReset;
  int _newind;
  int err;
} rx_solving_options_ind;

typedef struct {
  rx_solving_options_ind *subjects;
  rx_solving_options *op;
  int nsub;
  int nsim;
  int nall;
  int nevid9;
  int nobs;
  int nobs2;
  int nr;
  int add_cov;
  int matrix;
  int needSort;
  int nMtime;
  double stateTrim;
  int *stateIgnore;
  int nCov0;
  int *cov0;
  int nKeep0;
  int nKeepF;
  int istateReset;
} rx_solve;
  
typedef void (*t_set_solve)(rx_solve *);
typedef rx_solve *(*t_get_solve)();

typedef double *(*t_get_theta)(double *theta);
typedef void *(*t_assignFuns)();

rx_solve *getRxSolve_();
rx_solve *getRxSolve2_();
rx_solve *getRxSolve(SEXP ptr);

void par_solve(rx_solve *rx);

rx_solving_options *getRxOp(rx_solve *rx);

SEXP RxODE_df(int doDose, int doTBS);
SEXP RxODE_par_df();

void rxOptionsIniEnsure(int mx);

void rxUpdateFuns(SEXP trans);

#define _eps sqrt(DOUBLE_EPS)
static double _powerDi(double x, double lambda, int yj)  __attribute__((unused));
static double _powerDi(double x, double lambda, int yj){
  double x0=x, ret, l2;
  switch(yj){
  case 3:
    return exp(x);
  case 2: 
    return x;
  case 0:
    if (lambda == 1.0) return (x+1.0);
    if (lambda == 0) return exp(x);
    // (x^lambda-1)/lambda=y
    // (lambda*y+1)^(1/lambda)
    x0 = x*lambda+1.0;
    if (x0 <= _eps) return _eps;
    ret = pow(x0, 1.0/lambda);
    if (ISNA(ret)) {
      // Warning?
      return _eps;
    }
    return ret;
  case 1:
    if (lambda == 1.0) return x;
    if (x >= 0){
      // log(x+1)= y; exp(y)-1=x
      if (lambda == 0) return expm1(x);
      // ((x+1)^lambda-1)/lambda=y
      // (y*lambda+1)^(1/y)-1=y
      return pow(x*lambda+1.0, 1.0/lambda)-1.0;
    } else {
      // (-(1-x)^(2-lambda)-1)/(2-lambda)
      if (lambda ==  2.0) return -expm1(-x);
      // (-(1-x)^(2-lambda)-1)/(2-lambda) = y
      l2 = (2.0 - lambda);
      return 1.0 - pow(1.0 - l2*x, 1.0/l2);
    }
  }
  return NA_REAL;
}

static double _powerD(double x, double lambda, int yj)  __attribute__((unused));
static double _powerD(double x, double lambda, int yj){
  double x0=x, l2;
  switch (yj){
  case 3:
    if (x <= _eps) x0= _eps;
    return log(x0);
  case 2:
    return x;
  case 0:
    if (lambda == 1.0) return x-1.0;
    if (x <= _eps) x0= _eps;
    if (lambda ==  0.0) return log(x0);
    return (pow(x0, lambda) - 1.0)/lambda;
  case 1:
    if (lambda == 1.0) return x;
    if (x >= 0){
      if (lambda == 0) return log1p(x);
      return (pow(x + 1.0, lambda) - 1.0)/lambda;
    } else {
      if (lambda == 2.0) return -log1p(-x);
      l2 = 2.0 - lambda;
      return (1.0 - pow(1.0 - x, l2))/l2;
    }
  }
  return NA_REAL;
}

static double _powerDD(double x, double lambda, int yj)  __attribute__((unused));
static double _powerDD(double x, double lambda, int yj){
  double x0 = x;
  switch(yj){
  case 3:
    if (x <= _eps) return x0 = _eps;
    return 1/x0;
  case 2:
    return 1.0;
  case 0:
    if (lambda == 1.0) return 1.0;
    if (x <= _eps) return x0 = _eps;
    if (lambda == 0.0) return 1/x0;
    // pow(x,lambda)/lambda - 1/lambda
    return pow(x0, lambda-1);
  case 1:
    if (lambda ==  1.0) return 1.0;
    if (x >= 0){
      if (lambda == 0.0) return 1.0/(x + 1.0);
      return pow(x + 1.0, lambda-1.0);
    } else {
      if (lambda == 2.0) return -1/(1.0 - x);
      return pow(1.0 - x, 1.0-lambda);
    }
  }
  return NA_REAL;
}

static double _powerDDD(double x, double lambda, int yj) __attribute__((unused));
static double _powerDDD(double x, double lambda, int yj){
  double x0 = x;
  switch(yj){
  case 3:
    if (x <= _eps) x0 = _eps;
    return -1/(x0*x0);
  case 2: 
    return 0;
  case 0:
    if (lambda == 1.0) return 0;
    if (x <= _eps) return x0 = _eps;
    if (lambda == 0.0) return -1/(x0*x0);
    // pow(x,lambda)/lambda - 1/lambda
    return (lambda-1)*pow(x0, lambda-2);
  case 1:
    if (lambda == 1.0) return 0;
    if (x >= 0){
      if (lambda ==  0.0) return -1/((x + 1.0)*(x + 1.0));
      return (lambda-1.0)*pow(x + 1.0, lambda-2.0);
    } else {
      if (lambda == 2.0) return -1/((1.0 - x)*(1.0 - x));
      return -(1.0-lambda)*pow(1.0 - x, -lambda);
    }
  }
  return NA_REAL;
}

static double _powerL(double x, double lambda, int yj) __attribute__((unused));
static double _powerL(double x, double lambda, int yj){
  double x0 = x;
  switch(yj){
  case 3:
    if (x <= _eps) x0 = _eps;
    return -log(x0);
  case 2:
    return 0;
  case 0:
    if (lambda == 1.0) return 0;
    if (x <= _eps) x0 = _eps;
    return (lambda - 1.0)*log(x0);
  case 1:
    if (x >= 0) return (lambda - 1.0)*log1p(x);
    return (1.0-lambda)*log1p(-x);
  }
  return NA_REAL;
  // d = 0.0 for cox box
  // d = 1.0 fo  Yeo- Johnson
  // logLik approximation
  // y^(lambda)/lambda - 1/lambda
  // dh/dy = y^(lambda-1)
  // log(dh/dy) = (lambda-1)*log(y) + log(lambda) 
  //
  // (x + 1.0)^(lambda)/lambda - 1/lambda
  // dh/dy = (x+1.0)^(lambda-1)
  // log(dh/dy) = (lambda-1)*log(x+1.0)
  
  // For negative values yj becomes
  // (-x+1)^(2-lambda)/(2-lambda) - 1/(2-lambda)
  // dh/dy = (-x+1)^(1-lambda)
  // log(dh/dy) = (1-lambda)*log(-x+1)
}

static double _powerDL(double x, double lambda, int yj) __attribute__((unused));
static double _powerDL(double x, double lambda, int yj){
  // d(logLik/dlambda)
  double x0 = x;
  switch (yj){
  case 3:
    if (x <= _eps) x0 = _eps;
    return log(x0);
  case 2:
    return 0;
  case 0:
    if (lambda == 1.0) return 0;
    if (x <= _eps) x0 = _eps;
    return log(x0);
  case 1:
    if (lambda == 1.0) return 0;
    if (x >= 0) return log1p(x);
    return -log1p(x);
  }
  return NA_REAL;
  // d = 0.0 for cox box
  // d = 1.0 fo  Yeo- Johnson
  // logLik approximation
  // y^(lambda)/lambda - 1/lambda
  // dh/dy = y^(lambda-1)
  // log(dh/dy) = (lambda-1)*log(y) + log(lambda) 
  //
  // (x + 1.0)^(lambda)/lambda - 1/lambda
  // dh/dy = (x+1.0)^(lambda-1)
  // log(dh/dy) = (lambda-1)*log(x+1.0)
  
  // For negative values yj becomes
  // (-x+1)^(2-lambda)/(2-lambda) - 1/(2-lambda)
  // dh/dy = (-x+1)^(1-lambda)
  // log(dh/dy) = (1-lambda)*log(-x+1)

}


#endif
#if defined(__cplusplus)
}
#endif
