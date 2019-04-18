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

#define JAC_Rprintf Rprintf
#define _idx (&_solveData->subjects[_cSub])->idx
#define JAC0_Rprintf if ( (&_solveData->subjects[_cSub])->jac_counter == 0) Rprintf
#define ODE_Rprintf Rprintf
#define ODE0_Rprintf if ( (&_solveData->subjects[_cSub])->dadt_counter == 0) Rprintf
#define LHS_Rprintf Rprintf
#define _safe_log(a) (((a) <= 0) ? log(DOUBLE_EPS) : log(a))
#define safe_zero(a) ((a) == 0 ? DOUBLE_EPS : (a))
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
#define _SR (INTEGER(stateRmS))
#define NEWIND (_solveData->subjects[_cSub]._newind)
#define newind (_solveData->subjects[_cSub]._newind)
#define rx_lambda_ _solveData->subjects[_cSub].lambda
#define rx_yj_ _solveData->subjects[_cSub].yj
#define rxTBS(x, lm, yj)  _powerD(x,  lm, (int)(yj))
#define rxTBSi(x, lm, yj) _powerDi(x,  lm, (int)(yj))
#define rxTBSd(x, lm, yj) _powerDD(x, lm, (int)(yj))
#define rxTBSd2(x, lm, yj) _powerDDD(x, lm, (int)(yj))

// Types for par pointers.r
typedef double (*RxODE_fn) (double x);
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

typedef double (*solveLinB_p) (rx_solve *rx, unsigned int id, double t, int linCmt,
			       double d_A, double d_A2, double d_alpha,
			       double d_B, double d_B2, double d_beta,
			       double d_C, double d_C2, double d_gamma,
			       double d_ka, double d_tlag, double d_tlag2, double d_F, double d_F2,
			       double d_rate, double d_dur);


typedef void (*_update_par_ptr_p)(double t, unsigned int id, rx_solve *rx, int idx);

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
double solveLinB(rx_solve *rx, unsigned int id, double t, int linCmt,
		 double d_A, double d_A2, double d_alpha,
		 double d_B, double d_B2, double d_beta,
		 double d_C, double d_C2, double d_gamma,
		 double d_ka, double d_tlag, double d_tlag2,
		 double d_F, double d_F2,
		 double d_rate, double d_dur);
void _update_par_ptr(double t, unsigned int id, rx_solve *rx, int idx);
SEXP _RxODE_rxAssignPtr(SEXP);
#else
extern solveLinB_p solveLinB;
extern _update_par_ptr_p _update_par_ptr;
extern _rx_asgn _RxODE_rxAssignPtr;
#endif

#endif// __RxODE_model_H__
