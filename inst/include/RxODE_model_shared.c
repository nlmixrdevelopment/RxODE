_getRxSolve_t _getRxSolve_;
_simfun simeps;
_simfun simeta;

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

_update_par_ptr_p _update_par_ptr=NULL;
_getParCov_p _getParCov=NULL;
linCmtA_p linCmtA;
linCmtA_p linCmtC;
linCmtB_p linCmtB;
_rx_asgn _RxODE_rxAssignPtr =NULL;
_rx_asgn _rxQr =NULL;

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

double _prod(double *input, double *p, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  return _prodPS(input, p, n, type);
}

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

double _sign(unsigned int n, ...) {
  va_list valist;
  va_start(valist, n);
  double s = 1;
  for (unsigned int i = 0; i < n; i++) {
    s = sign(va_arg(valist, double))*s;
    if (s == 0){
      break;
    }
  }
  va_end(valist);
  return s;
}

double _max(unsigned int n, ...) {
  va_list valist;
  va_start(valist, n);
  double mx = NA_REAL;
  double tmp = 0;
  if (n >= 1){
    mx = va_arg(valist, double);
    for (unsigned int i = 1; i < n; i++) {
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

void _assignFuns0() {
  _getRxSolve_ = (_getRxSolve_t) R_GetCCallable("RxODE","getRxSolve_");
  _assign_ptr=(RxODE_assign_ptr) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
  _rxRmModelLib=(_rxRmModelLibType) R_GetCCallable("RxODE","rxRmModelLib");
  _rxGetModelLib=(_rxGetModelLibType) R_GetCCallable("RxODE","rxGetModelLib");
  _RxODE_rxAssignPtr=(_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxAssignPtr");
  _rxQr=(_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxQr");
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
  simeta =(_simfun) R_GetCCallable("RxODE", "simeta");
  simeps =(_simfun) R_GetCCallable("RxODE", "simeps");
  _compareFactorVal=(RxODE_compareFactorVal_fn) R_GetCCallable("RxODE", "compareFactorVal");
  _update_par_ptr = (_update_par_ptr_p) R_GetCCallable("RxODE","_update_par_ptr");
  _getParCov = (_getParCov_p) R_GetCCallable("RxODE","_getParCov");
  _solveData = _getRxSolve_();
}

void _assignFuns() {
  if (_assign_ptr == NULL){
    _assignFuns0();
  }
}
