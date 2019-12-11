#ifndef __RxODE_model_H__
#define __RxODE_model_H__
#ifdef _isRxODE_
#include "RxODE.h"
#else
#include <RxODE.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define _evid (&_solveData->subjects[_cSub])->evid[(&_solveData->subjects[_cSub])->ix[(&_solveData->subjects[_cSub])->idx]]
#define amt (isDose(_evid) ?  (&_solveData->subjects[_cSub])->dose[(&_solveData->subjects[_cSub])->ixds] : NA_REAL)
#define JAC_Rprintf Rprintf
#define _idx (&_solveData->subjects[_cSub])->idx
#define JAC0_Rprintf if ( (&_solveData->subjects[_cSub])->jac_counter == 0) Rprintf
#define ODE_Rprintf Rprintf
#define ODE0_Rprintf if ( (&_solveData->subjects[_cSub])->dadt_counter == 0) Rprintf
#define LHS_Rprintf Rprintf
#define _safe_log(a) (&_solveData->safeZero ? (((a) <= 0) ? log(DOUBLE_EPS) : log(a)) : log(a))
#define safe_zero(a) (&_solveData->safeZero ? ((a) == 0 ? DOUBLE_EPS : (a)) : (a))
#define _as_zero(a) (fabs(a) < sqrt(DOUBLE_EPS) ? 0.0 : a)
#define factorial(a) exp(lgamma1p(a))
#define sign_exp(sgn, x)(((sgn) > 0.0) ? exp(x) : (((sgn) < 0.0) ? -exp(x) : 0.0))
#define R_pow(a, b) (((a) == 0 && (b) <= 0) ? R_pow(DOUBLE_EPS, b) : R_pow(a, b))
#define R_pow_di(a, b) (((a) == 0 && (b) <= 0) ? R_pow_di(DOUBLE_EPS, b) : R_pow_di(a, b))
#define Rx_pow(a, b) (((a) == 0 && (b) <= 0) ? R_pow(DOUBLE_EPS, b) : R_pow(a, b))
#define Rx_pow_di(a, b) (((a) == 0 && (b) <= 0) ? R_pow_di(DOUBLE_EPS, b) : R_pow_di(a, b))
#define abs_log1p(x) (((x) + 1.0 > 0.0) ? log1p(x) : (((x) + 1.0 > 0.0) ? log1p(-x) : 0.0))
#define abs_log(x) ((fabs(x) <= sqrt(DOUBLE_EPS)) ? log(sqrt(DOUBLE_EPS)) : (((x) > 0.0) ? log(x) ? (((x) == 0) ? 0.0 : log(-x))))
#define _IR (_solveData->subjects[_cSub].InfusionRate)
#define _ON (_solveData->subjects[_cSub].on)
#define _PP (_solveData->subjects[_cSub].par_ptr)
#define _PL (_solveData->subjects[_cSub].lhs)
#define _SR (INTEGER(stateRmS))
#define NEWIND (_solveData->subjects[_cSub]._newind)
#define newind (_solveData->subjects[_cSub]._newind)
#define rx_lambda_ _solveData->subjects[_cSub].lambda
#define rx_yj_ _solveData->subjects[_cSub].yj
#define rxTBS(x, lm, yj)  _powerD(x,  lm, (int)(yj))
#define rxTBSi(x, lm, yj) _powerDi(x,  lm, (int)(yj))
#define rxTBSd(x, lm, yj) _powerDD(x, lm, (int)(yj))
#define rxTBSd2(x, lm, yj) _powerDDD(x, lm, (int)(yj))
#undef rnorm
#define rnorm(x,y) rxnorm(x,y)
#define rxnorm1(x) rxnorm(x, 1.0)
#define rnorm1(x) rxnorm(x, 1.0)
#define rnormV(x,y) rxnormV(x,y)
#define rxnormV1(x) rxnormV(x, 1.0)
#define rnormV1(x) rxnormV(x, 1.0)
#define rxcauchy1(x) rxcauchy(x, 1.0)
#undef rchisq
#define rchisq(x) rxchisq(x)
#undef rexp
#define rexp(x) rxexp(x)
#undef rgamma
#define rgamma(x,y) rxgamma(x,y)
#define rgamma1(x) rxgamma(x,1.0)
#define rxgamma1(x) rxgamma(x,1.0)
#undef rgeom
#define rgeom(x) rxgeom(x)
#undef rpois
#define rpois(x) rxpois(x)
#undef runif
#define runif(x,y) rxunif(x,y)
#define runif1(x) rxunif(x,1.0)
#define rxunif1(x) rxunif(x,1.0)
#undef rweibull
#define rweibull(x,y) rxweibull(x,y)
#define rxweibull1(x) rxweibull(x,1.0)
#define rweibull1(x) rxweibull(x,1.0)

// Types for par pointers.r
typedef double (*RxODE_fn) (double x);
typedef int (*RxODE_ifn) (double x);
typedef double (*RxODE_fn2) (double x, double y);
typedef double (*RxODE_fn3i) (double x, double y, int i);
typedef double (*RxODE_fn2i) (double x, int i);
typedef int (*RxODE_fn0i) ();
typedef double (*RxODE_vec) (int val, rx_solve *rx, unsigned int id);
typedef double (*RxODE_val) (rx_solve *rx, unsigned int id);
typedef void (*RxODE_assign_ptr)(SEXP);
typedef void (*RxODE_ode_solver_old_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc);

typedef void (*_rxRmModelLibType)(const char *inp);
typedef SEXP (*_rxGetModelLibType)(const char *s);
typedef  SEXP (*_rx_asgn) (SEXP objectSEXP);
typedef int(*_rxIsCurrentC_type)(SEXP);
typedef double(*_rxSumType)(double *, int, double *, int, int);

double _sum(double *input, double *pld, int m, int type, int n, ...);

typedef double(*_rxProdType)(double*, double*, int, int);

double _prod(double *input, double *p, int type, int n, ...);

double _sign(unsigned int n, ...);

double _max(unsigned int n, ...);

double _min(unsigned int n, ...);

double _transit4P(double t, unsigned int id, double n, double mtt, double bio);

double _transit3P(double t, unsigned int id, double n, double mtt);

typedef double (*linCmtA_p) (rx_solve *rx, unsigned int id, double t, int linCmt,
			     int ncmt, int trans, double d_ka,
			     double p1, double v1,
			     double p2, double p3,
			     double p4, double p5,
			     double d_tlag, double d_tlag2, double d_F, double d_F2,
			     // Rate and dur can only apply to central compartment even w/ oral dosing
			     // Therefore, only 1 model rate is possible with RxODE
			     double d_rate, double d_dur);

typedef double (*linCmtB_p) (rx_solve *rx, unsigned int id, double t, int linCmt,
	       int i_cmt, int trans, int val,
	       double dd_p1, double dd_v1,
	       double dd_p2, double dd_p3,
	       double dd_p4, double dd_p5,
	       double dd_ka,
	       double dd_tlag, double dd_tlag2,
	       double dd_F, double dd_F2,
	       double dd_rate, double dd_dur);


typedef void (*_update_par_ptr_p)(double t, unsigned int id, rx_solve *rx, int idx);

typedef double (*_getParCov_p)(unsigned int id, rx_solve *rx, int parNo, int idx);

typedef rx_solve *(*_getRxSolve_t)();

_getRxSolve_t _getRxSolve_;

void _assignFuns();


extern RxODE_assign_ptr _assign_ptr;
extern _rxRmModelLibType _rxRmModelLib;
extern _rxGetModelLibType _rxGetModelLib;
extern RxODE_ode_solver_old_c _old_c;
extern RxODE_fn0i _ptr;
extern _rxIsCurrentC_type _rxIsCurrentC;
extern _rxSumType _sumPS;
extern _rxProdType _prodPS;
extern RxODE_fn0i _prodType;
extern RxODE_fn0i _sumType;
extern rx_solve *_solveData;

#ifdef _isRxODE_
double linCmtA(rx_solve *rx, unsigned int id, double t, int linCmt,
	       int ncmt, int trans, double d_ka,
	       double p1, double v1,
	       double p2, double p3,
	       double p4, double p5,
	       double d_tlag, double d_tlag2, double d_F, double d_F2,
	       // Rate and dur can only apply to central compartment even w/ oral dosing
	       // Therefore, only 1 model rate is possible with RxODE
	       double d_rate, double d_dur);

double linCmtB(rx_solve *rx, unsigned int id, double t, int linCmt,
	       int i_cmt, int trans, int val,
	       double dd_p1, double dd_v1,
	       double dd_p2, double dd_p3,
	       double dd_p4, double dd_p5,
	       double dd_ka,
	       double dd_tlag, double dd_tlag2,
	       double dd_F, double dd_F2,
	       double dd_rate, double dd_dur);
void _update_par_ptr(double t, unsigned int id, rx_solve *rx, int idx);
SEXP _RxODE_rxAssignPtr(SEXP);
#else
extern linCmtA_p linCmtA;
extern linCmtB_p linCmtB;
extern _update_par_ptr_p _update_par_ptr;
extern _getParCov_p _getParCov;
extern _rx_asgn _RxODE_rxAssignPtr;
typedef int (*RxODE_rxbinom) (int n, double prob);
extern RxODE_rxbinom rxbinom;
extern RxODE_fn2 rxcauchy;
extern RxODE_fn rxchisq;
extern RxODE_fn rxexp;
extern RxODE_fn2 rxf;
extern RxODE_ifn rxgeom;
extern RxODE_fn2 rxnorm;
extern RxODE_fn2 rxnormV;
extern RxODE_fn2 rxgamma;
extern RxODE_fn2 rxbeta;
extern RxODE_ifn rxpois;
extern RxODE_fn rxt_;
extern RxODE_fn2 rxunif;
extern RxODE_fn2 rxweibull;
#endif

#endif// __RxODE_model_H__
