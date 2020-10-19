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
#define _as_zero(a) (&_solveData->safeZero && fabs(a) < sqrt(DOUBLE_EPS) ? 0.0 : a)
#define _as_dbleps(a) (&_solveData->safeZero && fabs(a) < sqrt(DOUBLE_EPS) ? ((a) < 0 ? -sqrt(DOUBLE_EPS)  : sqrt(DOUBLE_EPS)) : a)
#define factorial(a) exp(lgamma1p(a))
#define sign_exp(sgn, x)(((sgn) > 0.0) ? exp(x) : (((sgn) < 0.0) ? -exp(x) : 0.0))
#define Rx_pow(a, b) R_pow(a, b)
#define Rx_pow_di(a, b) R_pow_di(a, b)
#define abs_log1p(x) (((x) + 1.0 > 0.0) ? log1p(x) : (((x) + 1.0 > 0.0) ? log1p(-x) : 0.0))
#define abs_log(x) ((&_solveData->safeZero && fabs(x) <= sqrt(DOUBLE_EPS)) ? log(sqrt(DOUBLE_EPS)) : (((x) > 0.0) ? log(x) ? (((x) == 0) ? 0.0 : log(-x))))
#define _IR (_solveData->subjects[_cSub].InfusionRate)
#define _ON (_solveData->subjects[_cSub].on)
#define _PP (_solveData->subjects[_cSub].par_ptr)
#define _PL (_solveData->subjects[_cSub].lhs)
#define _SR (INTEGER(stateRmS))
#define NEWIND (_solveData->subjects[_cSub]._newind)
#define newind (_solveData->subjects[_cSub]._newind)
#define rx_lambda_ _solveData->subjects[_cSub].lambda
#define rx_yj_ _solveData->subjects[_cSub].yj
#define rx_hi_ _solveData->subjects[_cSub].logitHi
#define rx_low_ _solveData->subjects[_cSub].logitLow
#define rxTBS(x, lm, yj, hi, low)  _powerD(x,  lm, (int)(yj), hi, low)
#define rxTBSi(x, lm, yj, hi, low) _powerDi(x,  lm, (int)(yj), hi, low)
#define rxTBSd(x, lm, yj, hi, low) _powerDD(x, lm, (int)(yj), hi, low)
#define rxTBSd2(x, lm, yj, hi, low) _powerDDD(x, lm, (int)(yj), hi, low)
#undef rbeta
#define rbeta(ind, x, y) rxbeta(ind, x, y)
#undef rnorm
#define rnorm(ind,x,y) rxnorm(ind, x,y)
#define rxnorm1(x) rxnorm(&_solveData->subjects[_cSub], x, 1.0)
#define rnorm1(x) rxnorm(&_solveData->subjects[_cSub],x, 1.0)
#define rnormV(ind, x,y) rxnormV(ind,x,y)
#define rxnormV1(x) rxnormV(&_solveData->subjects[_cSub],x, 1.0)
#define rnormV1(x) rxnormV(&_solveData->subjects[_cSub],x, 1.0)
#undef rcauchy 
#define rcauchy(ind, x, y) rxcauchy(ind,x,y)
#define rxcauchy1(x) rxcauchy(&_solveData->subjects[_cSub],x, 1.0)
#undef rchisq
#define rchisq(ind, x) rxchisq(ind, x)
#undef rexp
#define rexp(ind, x) rxexp(ind, x)
#undef rgamma
#define rgamma(ind, x,y) rxgamma(ind, x,y)
#define rgamma1(x) rxgamma(&_solveData->subjects[_cSub], x,1.0)
#define rxgamma1(x) rxgamma(&_solveData->subjects[_cSub], x,1.0)
#undef rgeom
#define rgeom(ind,x) rxgeom(ind,x)
#undef rpois
#define rpois(ind,x) rxpois(ind,x)
#undef runif
#define runif(ind,x,y) rxunif(ind,x,y)
#define runif1(x) rxunif(&_solveData->subjects[_cSub],x,1.0)
#define rxunif1(x) rxunif(&_solveData->subjects[_cSub],x,1.0)
#undef rweibull
#define rweibull(ind,x,y) rxweibull(ind,x,y)
#define rxweibull1(x) rxweibull(&_solveData->subjects[_cSub],x,1.0)
#define rweibull1(x) rxweibull(&_solveData->subjects[_cSub],x,1.0)
#define _pnorm1(x) pnorm(x, 0.0, 1.0, 1, 0)
#define _pnorm2(x, mu) pnorm(x, mu, 1.0, 1, 0)
#define _pnorm3(x, mu, sd) pnorm(x, mu, sd, 1, 0)
#define _qnorm1(x) qnorm(x, 0.0, 1.0, 1, 0)
#define _qnorm2(x, mu) qnorm(x, mu, 1.0, 1, 0)
#define _qnorm3(x, mu, sd) qnorm(x, mu, sd, 1, 0)
#define probitt(x) qnorm(x, 0.0, 1.0, 1, 0)
#define _logit1(x) logit(x, 0.0, 1.0)
#define _logit2(x, y) logit(x, y, 1.0)
#define _expit1(x) expit(x, 0.0, 1.0)
#define _expit2(x, y) expit(x, y, 1.0)
#define _tad0() (t-_solveData->subjects[_cSub].tlast)
#define _tad1(x) (t-_solveData->subjects[_cSub].tlastS[_swapCmt(x)])
#define _tafd0()  (t-_solveData->subjects[_cSub].tfirst)
#define _tafd1(x) (t-_solveData->subjects[_cSub].tfirstS[_swapCmt(x)])
#define _tlast0() _solveData->subjects[_cSub].tlast
#define _tlast1(x) _solveData->subjects[_cSub].tlastS[_swapCmt(x)]
#define _tfirst0()  _solveData->subjects[_cSub].tfirst
#define _tfirst1(x) _solveData->subjects[_cSub].tfirstS[_swapCmt(x)]
#undef rf
#define rf(ind, x, y) rxf(ind, x, y)
// int compareFactorVal(int val, const char *valStr, const char *cmpValue)
// equality_str2 : identifier_r ('!=' | '==' ) string;
#define _cmp2(val, valStr, type, cmpStr) (type ? _compareFactorVal(val, valStr, cmpStr) : !_compareFactorVal(val, valStr, cmpStr))
// equality_str1 : string ('!=' | '==' ) identifier_r; //type=1 is equal, type=0 not equal
#define _cmp1(cmpStr, type, val, valStr) (type ? _compareFactorVal(val, valStr, cmpStr) : !_compareFactorVal(val, valStr, cmpStr))

// Types for par pointers.r
typedef int (*RxODE_compareFactorVal_fn)(int val, const char *factor, const char *value);
typedef double (*RxODE_fn) (double x);
typedef int (*RxODE_ifn) (double x);
typedef double (*RxODE_fn2) (double x, double y);
typedef double (*RxODE_fn3) (double x, double y, double z);
typedef double (*RxODE_fn3i) (double x, double y, int i);
typedef double (*RxODE_fn2i) (double x, int i);
typedef int (*RxODE_fn0i) ();
typedef double (*RxODEi_fn) (rx_solving_options_ind* ind, double x);
typedef int (*RxODEi_ifn) (rx_solving_options_ind* ind, double x);
typedef double (*RxODEi_fn2) (rx_solving_options_ind* ind, double x, double y);
typedef double (*RxODEi_fn3i) (rx_solving_options_ind* ind, double x, double y, int i);
typedef double (*RxODEi_fn2i) (rx_solving_options_ind* ind, double x, int i);
typedef int (*RxODEi_fn0i) (rx_solving_options_ind* ind);
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
			     double d_rate, double d_dur,
			     double d_rate2, double d_dur2);

typedef double (*linCmtB_p) (rx_solve *rx, unsigned int id, double t, int linCmt,
			     int i_cmt, int trans, int val,
			     double dd_p1, double dd_v1,
			     double dd_p2, double dd_p3,
			     double dd_p4, double dd_p5,
			     double dd_ka,
			     double dd_tlag, double dd_tlag2,
			     double dd_F, double dd_F2,
			     double dd_rate, double dd_dur,
			     double dd_rate2, double dd_dur2);


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
	       double d_rate, double d_dur, double d_rate2, double d_dur2);

double linCmtB(rx_solve *rx, unsigned int id, double t, int linCmt,
	       int i_cmt, int trans, int val,
	       double dd_p1, double dd_v1,
	       double dd_p2, double dd_p3,
	       double dd_p4, double dd_p5,
	       double dd_ka,
	       double dd_tlag, double dd_tlag2,
	       double dd_F, double dd_F2,
	       double dd_rate, double dd_dur, double dd_rate2, double dd_dur2);
void _update_par_ptr(double t, unsigned int id, rx_solve *rx, int idx);
SEXP _RxODE_rxAssignPtr(SEXP);
#else
extern linCmtA_p linCmtA;
extern linCmtB_p linCmtB;
extern _update_par_ptr_p _update_par_ptr;
extern _getParCov_p _getParCov;
extern _rx_asgn _RxODE_rxAssignPtr;
typedef int (*RxODEi_rxbinom) (rx_solving_options_ind* ind, int n, double prob);
extern RxODEi_rxbinom rxbinom;
extern RxODEi_fn2 rxcauchy;
extern RxODEi_fn rxchisq;
extern RxODEi_fn rxexp;
extern RxODEi_fn2 rxf;
extern RxODEi_ifn rxgeom;
extern RxODEi_fn2 rxnorm;
extern RxODEi_fn2 rxnormV;
extern RxODEi_fn2 rxgamma;
extern RxODEi_fn2 rxbeta;
extern RxODEi_ifn rxpois;
extern RxODEi_fn rxt_;
extern RxODEi_fn2 rxunif;
extern RxODEi_fn2 rxweibull;
extern RxODE_compareFactorVal_fn _compareFactorVal;
#endif

#endif// __RxODE_model_H__
