#pragma once
#ifndef __RxODE_H__
#define __RxODE_H__
#define isDose(evid) ((evid) == 3 || (evid) >= 100)
#define isObs(evid) ((evid) == 0 || (evid) == 2 || ((evid) >= 10 && (evid) <= 99))
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
  int nDisplayProgress;
  int ncoresRV;
  int isChol;
  int *svar;
  int abort;
  int minSS;
  int maxSS;
  double atolSS;
  double rtolSS;
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
  double *lhs;
  int  *evid;
  int *rc;
  double *cov_ptr;
  int n_all_times;
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
} rx_solving_options_ind;

typedef struct {
  rx_solving_options_ind *subjects;
  rx_solving_options *op;
  int nsub;
  int nsim;
  int nall;
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
} rx_solve;
  
typedef void (*t_set_solve)(rx_solve *);
typedef rx_solve *(*t_get_solve)();

typedef double *(*t_get_theta)(double *theta);


rx_solve *getRxSolve_();
rx_solve *getRxSolve(SEXP ptr);

void par_solve(rx_solve *rx);

rx_solving_options *getRxOp(rx_solve *rx);

SEXP RxODE_df(int doDose, int doTBS);
SEXP RxODE_par_df();

rx_solving_options_ind *rxOptionsIniEnsure(int mx);

void rxUpdateFuns(SEXP trans);

#endif
#if defined(__cplusplus)
}
#endif
