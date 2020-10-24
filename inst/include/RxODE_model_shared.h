#ifndef __RxODE_model_shared_H__
#define __RxODE_model_shared_H__
#include <RxODE.h>

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
#define normcdf(x) phi(x)
#define _getIndSim(id, val) (_solveData->subjects[_cSub].isIni == 1 ? \
			     (_solveData->subjects[_cSub].simIni[id] = (val)) : _solveData->subjects[_cSub].simIni[id])
#undef rbeta
#define rbeta(ind, x, y) rxbeta(ind, x, y)
#define ribeta(ind, id, x, y) _getIndSim(id, rxbeta(ind, x, y))
#undef rnorm
#define rnorm(ind,x,y) rxnorm(ind, x,y)
#define rxnorm1(x) rxnorm(&_solveData->subjects[_cSub], x, 1.0)
#define rnorm1(x) rxnorm(&_solveData->subjects[_cSub],x, 1.0)
#define rxnormV1(x) rxnorm(&_solveData->subjects[_cSub], x, 1.0)
#define rinorm1(id, x) rinorm(&_solveData->subjects[_cSub], id, x, 1.0)
#define rinormV1(id, x) rinorm(&_solveData->subjects[_cSub], id, x, 1.0)

// FIXME: need to use same scheme here
#define rnormV(ind, x,y) rxnormV(ind,x,y)
#define rnormV1(ind, id, x) rxnormV(ind, id, x, 1.0)
  
#undef rcauchy 
#define rcauchy(ind, x, y) rxcauchy(ind,x,y)
#define rxcauchy1(x) rxcauchy(&_solveData->subjects[_cSub],x, 1.0)
#define ricauchy1(id, x) ricauchy(&_solveData->subjects[_cSub], id, x, 1.0)
#undef rchisq
#define rchisq(ind, x) rxchisq(ind, x)
#undef rexp
#define rexp(ind, x) rxexp(ind, x)
#undef rgamma
#define rgamma(ind, x,y) rxgamma(ind, x,y)
#define rgamma1(x) rxgamma(&_solveData->subjects[_cSub], x,1.0)
#define rxgamma1(x) rxgamma(&_solveData->subjects[_cSub], x,1.0)
#define rigamma1(id, x) rigamma(&_solveData->subjects[_cSub], id, x,1.0)
#undef rgeom
#define rgeom(ind,x) rxgeom(ind,x)
#undef rpois
#define rpois(ind,x) rxpois(ind,x)
#undef runif
#define runif(ind,x,y) rxunif(ind,x,y)
#define runif1(x) rxunif(&_solveData->subjects[_cSub],x,1.0)
#define rxunif1(x) rxunif(&_solveData->subjects[_cSub],x,1.0)
#define riunif1(id, x) riunif(&_solveData->subjects[_cSub],id, x,1.0)
#undef rweibull
#define rweibull(ind,x,y) rxweibull(ind,x,y)
#define rxweibull1(x) rxweibull(&_solveData->subjects[_cSub], x, 1.0)
#define riweibull1(id, x) riweibull(&_solveData->subjects[_cSub], id, x, 1.0)
#define rweibull1(x) rxweibull(&_solveData->subjects[_cSub], x, 1.0)
#define _pnorm1(x) pnorm(x, 0.0, 1.0, 1, 0)
#define _pnorm2(x, mu) pnorm(x, mu, 1.0, 1, 0)
#define _pnorm3(x, mu, sd) pnorm(x, mu, sd, 1, 0)
#define _qnorm1(x) qnorm(x, 0.0, 1.0, 1, 0)
#define _qnorm2(x, mu) qnorm(x, mu, 1.0, 1, 0)
#define _qnorm3(x, mu, sd) qnorm(x, mu, sd, 1, 0)
#define norminv(x) qnorm(x, 0.0, 1.0, 1, 0)
#define probitt(x) qnorm(x, 0.0, 1.0, 1, 0)
#define _logit1(x) logit(x, 0.0, 1.0)
#define _logit2(x, y) logit(x, y, 1.0)
#define _expit1(x) expit(x, 0.0, 1.0)
#define _expit2(x, y) expit(x, y, 1.0)
#define _invLogit1(x) expit(x, 0.0, 1.0)
#define _invLogit2(x, y) expit(x, y, 1.0)
#define _logitInv1(x) expit(x, 0.0, 1.0)
#define _logitInv2(x, y) expit(x, y, 1.0)
#define _tad0() (t-_solveData->subjects[_cSub].tlast)
#define _tad1(x) (t-_solveData->subjects[_cSub].tlastS[x])
#define _tafd0()  (t-_solveData->subjects[_cSub].tfirst)
#define _tafd1(x) (t-_solveData->subjects[_cSub].tfirstS[x])
#define _tlast0() _solveData->subjects[_cSub].tlast
#define _tlast1(x) _solveData->subjects[_cSub].tlastS[x]
#define _tfirst0()  _solveData->subjects[_cSub].tfirst
#define _tfirst1(x) _solveData->subjects[_cSub].tfirstS[x]
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

typedef int (*RxODEi2_fn0i) (rx_solving_options_ind* ind, int id);
typedef double (*RxODEi2_fn) (rx_solving_options_ind* ind, int id, double x);
typedef int (*RxODEi2_ifn) (rx_solving_options_ind* ind, int id, double x);
typedef double (*RxODEi2_fn2) (rx_solving_options_ind* ind, int id, double x, double y);
typedef double (*RxODEi2_fn3i) (rx_solving_options_ind* ind, int id, double x, double y, int i);
typedef double (*RxODEi2_fn2i) (rx_solving_options_ind* ind, int id, double x, int i);

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

rx_solve *_solveData = NULL;
RxODE_assign_ptr _assign_ptr = NULL;
_rxRmModelLibType _rxRmModelLib = NULL;
_rxGetModelLibType _rxGetModelLib = NULL;
RxODE_ode_solver_old_c _old_c = NULL;
RxODE_fn0i _ptrid=NULL;
_rxIsCurrentC_type _rxIsCurrentC=NULL;
_rxSumType _sumPS = NULL;
_rxProdType _prodPS = NULL;

RxODE_fn0i _prodType = NULL;
RxODE_fn0i _sumType = NULL;


double _sum(double *input, double *pld, int m, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  double ret = _sumPS(input, n, pld, m, type);
  if (type == 2 && m < 0){
    for (int i = -m; i--;){
      pld[i] = 0.0;
    }
  }
  return ret;
}


double _prod(double *input, double *p, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  return _prodPS(input, p, n, type);
}

double _sign(unsigned int n, ...){
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

double _max(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double mx = NA_REAL;
  double tmp = 0;
  if (n >= 1){
    mx = va_arg(valist, double);
    for (unsigned int i = 1; i < n; i++){
      tmp = va_arg(valist, double);
      if (tmp>mx) mx=tmp;
    }
    va_end(valist);
  }
  return mx;
}

double _min(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double mn = NA_REAL;
  double tmp = 0;
  if (n >= 1){
    mn = va_arg(valist, double);
    for (unsigned int i = 1; i < n; i++){
      tmp = va_arg(valist, double);
      if (tmp<mn) mn=tmp;
    }
    va_end(valist);
  }
  return mn;
}

double _transit4P(double t, unsigned int id, double n, double mtt, double bio){
  double ktr = (n+1)/mtt;
  double lktr = log(n+1)-log(mtt);
  double tc = (t-(_solveData->subjects[id].tlast));
  return exp(log(bio*(_solveData->subjects[id].podo))+lktr+n*(lktr+log(tc))-ktr*(tc)-lgamma1p(n));
}

double _transit3P(double t, unsigned int id, double n, double mtt){
  double ktr = (n+1)/mtt;
  double lktr = log(n+1)-log(mtt);
  double tc = (t-(_solveData->subjects[id].tlast));
  return exp(log(_solveData->subjects[id].podo)+lktr+n*(lktr+log(tc))-ktr*(tc)-lgamma1p(n));
}

_update_par_ptr_p _update_par_ptr=NULL;
_getParCov_p _getParCov=NULL;
linCmtA_p linCmtA;
linCmtA_p linCmtC;
linCmtB_p linCmtB;
_rx_asgn _RxODE_rxAssignPtr =NULL;
typedef int (*RxODEi_rxbinom) (rx_solving_options_ind* ind, int n, double prob);
typedef int (*RxODEi2_ribinom) (rx_solving_options_ind* ind, int id, int n, double prob);

RxODE_fn phi;
RxODE_fn3 logit;
RxODE_fn3 expit;
RxODE_fn2 gammap;
RxODE_fn2 gammaq;
RxODE_fn2 lowergamma;
RxODE_fn2 uppergamma;
RxODE_fn2 gammapInv;
RxODE_fn2 gammapDer;
RxODE_fn2 gammapInva;
RxODE_fn2 gammaqInv;
RxODE_fn2 gammaqInva;

RxODEi_fn2 rxnorm;
RxODEi_fn2 rxnormV;
RxODEi_rxbinom rxbinom;
RxODEi_fn2 rxcauchy;
RxODEi_fn rxchisq;
RxODEi_fn rxexp;
RxODEi_fn2 rxf;
RxODEi_ifn rxgeom;
RxODEi_fn2 rxgamma;
RxODEi_fn2 rxbeta;
RxODEi_ifn rxpois;
RxODEi_fn rxt_;
RxODEi_fn2 rxunif;
RxODEi_fn2 rxweibull;

RxODEi2_fn2 rinorm;
RxODEi2_fn2 rinormV;
RxODEi2_ribinom ribinom;
RxODEi2_fn2 ricauchy;
RxODEi2_fn richisq;
RxODEi2_fn riexp;
RxODEi2_fn2 rif;
RxODEi2_ifn rigeom;
RxODEi2_fn2 rigamma;
RxODEi2_fn2 ribeta;
RxODEi2_ifn ripois;
RxODEi2_fn rit_;
RxODEi2_fn2 riunif;
RxODEi2_fn2 riweibull;

RxODE_compareFactorVal_fn _compareFactorVal;
void _assignFuns(){
  if (_assign_ptr == NULL){
    _getRxSolve_ = (_getRxSolve_t) R_GetCCallable("RxODE","getRxSolve_");
    _assign_ptr=(RxODE_assign_ptr) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
    _rxRmModelLib=(_rxRmModelLibType) R_GetCCallable("RxODE","rxRmModelLib");
    _rxGetModelLib=(_rxGetModelLibType) R_GetCCallable("RxODE","rxGetModelLib");
    _RxODE_rxAssignPtr=(_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxAssignPtr");
    _rxIsCurrentC = (_rxIsCurrentC_type)R_GetCCallable("RxODE","rxIsCurrentC");
    _sumPS  = (_rxSumType) R_GetCCallable("PreciseSums","PreciseSums_sum_r");
    _prodPS = (_rxProdType) R_GetCCallable("PreciseSums","PreciseSums_prod_r");
    _prodType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_prod_get");
    _sumType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_sum_get");
    _ptrid=(RxODE_fn0i)R_GetCCallable("RxODE", "RxODE_current_fn_pointer_id");
    linCmtA=(linCmtA_p)R_GetCCallable("RxODE", "linCmtA");
    linCmtB=(linCmtB_p)R_GetCCallable("RxODE", "linCmtB");
    linCmtC=(linCmtA_p)R_GetCCallable("RxODE", "linCmtC");
    
    rxnorm = (RxODEi_fn2)R_GetCCallable("RxODE", "rxnorm");
    rxnormV = (RxODEi_fn2)R_GetCCallable("RxODE", "rxnormV");
    rxbinom = (RxODEi_rxbinom)R_GetCCallable("RxODE","rxbinom") ;
    rxcauchy = (RxODEi_fn2)R_GetCCallable("RxODE","rxcauchy") ;
    rxchisq = (RxODEi_fn)R_GetCCallable("RxODE","rxchisq") ;
    rxexp = (RxODEi_fn)R_GetCCallable("RxODE","rxexp");
    rxf = (RxODEi_fn2)R_GetCCallable("RxODE","rxf") ;
    rxgeom = (RxODEi_ifn)R_GetCCallable("RxODE","rxgeom") ;
    rxgamma = (RxODEi_fn2)R_GetCCallable("RxODE","rxgamma") ;
    rxbeta = (RxODEi_fn2)R_GetCCallable("RxODE","rxbeta") ;
    rxpois = (RxODEi_ifn)R_GetCCallable("RxODE","rxpois") ;
    rxt_ = (RxODEi_fn)R_GetCCallable("RxODE","rxt_") ;
    rxunif = (RxODEi_fn2)R_GetCCallable("RxODE","rxunif") ;
    rxweibull = (RxODEi_fn2)R_GetCCallable("RxODE","rxweibull");

    rinorm = (RxODEi2_fn2)R_GetCCallable("RxODE", "rinorm");
    rinormV = (RxODEi2_fn2)R_GetCCallable("RxODE", "rinormV");
    ribinom = (RxODEi2_ribinom)R_GetCCallable("RxODE","ribinom") ;
    ricauchy = (RxODEi2_fn2)R_GetCCallable("RxODE","ricauchy") ;
    richisq = (RxODEi2_fn)R_GetCCallable("RxODE","richisq") ;
    riexp = (RxODEi2_fn)R_GetCCallable("RxODE","riexp");
    rif = (RxODEi2_fn2)R_GetCCallable("RxODE","rif") ;
    rigeom = (RxODEi2_ifn)R_GetCCallable("RxODE","rigeom") ;
    rigamma = (RxODEi2_fn2)R_GetCCallable("RxODE","rigamma") ;
    ribeta = (RxODEi2_fn2)R_GetCCallable("RxODE","ribeta") ;
    ripois = (RxODEi2_ifn)R_GetCCallable("RxODE","ripois") ;
    rit_ = (RxODEi2_fn)R_GetCCallable("RxODE","rit_") ;
    riunif = (RxODEi2_fn2)R_GetCCallable("RxODE","riunif") ;
    riweibull = (RxODEi2_fn2)R_GetCCallable("RxODE","riweibull");
    
    phi = (RxODE_fn)R_GetCCallable("RxODE","phi");
    gammap = (RxODE_fn2) R_GetCCallable("RxODE","gammap");
    gammaq = (RxODE_fn2) R_GetCCallable("RxODE","gammaq");
    gammapInv = (RxODE_fn2) R_GetCCallable("RxODE","gammapInv");
    gammapInva = (RxODE_fn2) R_GetCCallable("RxODE","gammapInva");
    gammaqInv = (RxODE_fn2) R_GetCCallable("RxODE","gammaqInv");
    gammaqInva = (RxODE_fn2) R_GetCCallable("RxODE","gammaqInva");
    uppergamma = (RxODE_fn2) R_GetCCallable("RxODE","uppergamma");
    lowergamma = (RxODE_fn2) R_GetCCallable("RxODE","lowergamma");
    gammapDer  = (RxODE_fn2) R_GetCCallable("RxODE","gammapDer");
    logit = (RxODE_fn3) R_GetCCallable("RxODE", "logit");
    expit = (RxODE_fn3) R_GetCCallable("RxODE", "expit");
    _compareFactorVal=(RxODE_compareFactorVal_fn) R_GetCCallable("RxODE", "compareFactorVal");
    _update_par_ptr = (_update_par_ptr_p) R_GetCCallable("RxODE","_update_par_ptr");
    _getParCov = (_getParCov_p) R_GetCCallable("RxODE","_getParCov");
    _solveData = _getRxSolve_();
  }
}

#endif // __RxODE_model_shared_H__
