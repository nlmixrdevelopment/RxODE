#ifndef __RxODE_model_shared_H__
#define __RxODE_model_shared_H__

#ifdef _isRxODE_
#include "RxODE_model.h"
#else
#include <RxODE_model.h>
#endif

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
  return _sumPS(input, n, pld, m, type);
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


#ifdef _isRxODE_
double solveLinB(rx_solve *rx, unsigned int id, double t, int linCmt,
		 double d_A, double d_A2, double d_alpha,
		 double d_B, double d_B2, double d_beta,
		 double d_C, double d_C2, double d_gamma,
		 double d_ka, double d_tlag, double d_tlag2,
		 double d_F, double d_F2,
		 double d_rate, double d_dur);
void _update_par_ptr(double t, unsigned int id, rx_solve *rx, int idx);
void _assignFuns(){
  if (_assign_ptr == NULL){
    _getRxSolve_ = (_getRxSolve_t) R_GetCCallable("RxODE","getRxSolve_");
    _assign_ptr=(RxODE_assign_ptr) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
    _rxRmModelLib=(_rxRmModelLibType) R_GetCCallable("RxODE","rxRmModelLib");
    _rxGetModelLib=(_rxGetModelLibType) R_GetCCallable("RxODE","rxGetModelLib");
    /* _RxODE_rxAssignPtr=(_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxAssignPtr"); */
    _rxIsCurrentC = (_rxIsCurrentC_type)R_GetCCallable("RxODE","rxIsCurrentC");
    _sumPS  = (_rxSumType) R_GetCCallable("PreciseSums","PreciseSums_sum_r");
    _prodPS = (_rxProdType) R_GetCCallable("PreciseSums","PreciseSums_prod_r");
    _prodType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_prod_get");
    _sumType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_sum_get");
    _ptrid=(RxODE_fn0i)R_GetCCallable("RxODE", "RxODE_current_fn_pointer_id");
    /* solveLinB=(solveLinB_p)R_GetCCallable("RxODE", "solveLinB"); */
    /* _update_par_ptr = (_update_par_ptr_p) R_GetCCallable("RxODE","_update_par_ptr"); */
    _solveData = _getRxSolve_();
  }
}
SEXP _RxODE_rxAssignPtr(SEXP);
#else
_update_par_ptr_p _update_par_ptr=NULL;
solveLinB_p solveLinB;
_rx_asgn _RxODE_rxAssignPtr =NULL;
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
    solveLinB=(solveLinB_p)R_GetCCallable("RxODE", "solveLinB");
    _update_par_ptr = (_update_par_ptr_p) R_GetCCallable("RxODE","_update_par_ptr");
    _solveData = _getRxSolve_();
  }
}
#endif // compiling in RxODE

#endif // __RxODE_model_shared_H__
