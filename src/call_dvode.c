#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "dop853.h"
#define max(a, b) ((a) > (b) ? (a) : (b))
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>
#include <PreciseSums.h>
#include "solve.h"

extern double RxODE_as_zero(double x){
  if (fabs(x) < sqrt(DOUBLE_EPS)){
    return(0.0);
  } else {
    return(x);
  }
}

extern double RxODE_safe_log(double x){
  if (x <= 0){
    // Warning?
    return log(DOUBLE_EPS);
  } else {
    return log(x);
  }
}

extern double RxODE_safe_zero(double x){
  if (x == 0){
    // Warning?
    return DOUBLE_EPS;
  } else {
    return(x);
  }
}

extern double RxODE_pow(double x, double y){
  if (x == 0 && y <= 0){
    return R_pow(DOUBLE_EPS, y);
  } else {
    return R_pow(x, y);
  }
}
extern double RxODE_pow_di(double x, int i){
  if (x == 0 && i <= 0){
    return R_pow_di(DOUBLE_EPS, i);
  } else {
    return R_pow_di(x, i);
  }
}

extern double RxODE_sign_exp(double sgn, double x){
  if (sgn > 0.0){
    return(exp(x));
  } else if (sgn < 0.0){
    return(-exp(x));
  } else {
    return(0.0);
  }
}

extern double RxODE_abs_log(double x){
  if  (fabs(x) <= sqrt(DOUBLE_EPS)){
    return log(sqrt(DOUBLE_EPS));
  } else if (x > 0.0){
    return log(x);
  } else if (x < 0.0){
    return log(-x);
  } else {
    return 0.0;
  }
}

extern double RxODE_abs_log1p(double x){
  if (x + 1.0 > 0.0){
    return(log1p(x));
  } else if (x + 1.0 > 0.0){
    return(log1p(-x));
  } else {
    return 0.0;
  }
}

extern double RxODE_factorial(double x){
  return exp(lgamma1p(x));
}

//--------------------------------------------------------------------------

// These are now allocated via R structures in Rcpp.
extern void RxODE_ode_free(){
}

void RxODE_ode_alloc(){
}

SEXP __mv;
extern void RxODE_assign_fn_pointers(SEXP mv){
  __mv = mv;
}

extern SEXP RxODE_get_mv(){
  return __mv;
}

extern void rxSolveOldC(int *neqa,
			double *theta,  //order:
			double *timep,
			int *evidp,
			int *ntime,
			double *initsp,
			double *dosep,
			double *retp,
			double *atol,
			double *rtol,
			int *stiffa,
			int *transit_abs,
			int *nlhsa,
			double *lhsp,
			int *rc){
  rx_solve *rx = (rx_solve*)(R_ExternalPtrAddr(VECTOR_ELT(VECTOR_ELT(__mv, 12), 11)));
  rx_solving_options *op = (rx_solving_options*)R_ExternalPtrAddr(rx->op);
  rx_solving_options_ind *inds = rx->subjects;
  rx_solving_options_ind *ind = &inds[0];
  ind->par_ptr = theta;
  ind->n_all_times  = *ntime;
  int neq = op->neq, i = 0;
  double *InfusionRate =ind->InfusionRate,
    *scale = op->scale;
  int *BadDose = ind->BadDose;
  // A bit paranoid -- make sure these are sane...
  for (i = 0; i < neq; i++){
    InfusionRate[i] = 0.0;
    scale[i] = 1.0;
    BadDose[i] = 0;
  }
  // Instead of having the correct length for idose, use idose length = length of ntime
  // Saves an additional for loop at the cost of a little memory.
  int *idose;
  idose = Calloc(*ntime,int);
  ind->idose = idose;
  ind->ndoses=0;
  for (i = 0; i < ind->n_all_times; i++){
    if (evidp[i]){
      ind->ndoses++;
      ind->idose[ind->ndoses-1] = i;
    }
  }
  op->do_par_cov = 0;
  // cov_ptr
  op->ncov              = 0;
  op->is_locf           = 0;
  // Solver Options
  op->ATOL = *atol;
  op->RTOL = *rtol;
  // Assign to default LSODA behvior, or 0
  op->HMIN           = 0;
  ind->HMAX          = 0;
  op->H0             = 0;
  op->MXORDN         = 0;
  op->MXORDS         = 0;
  op->mxstep         = 5000; // Not LSODA default but RxODE default
  // Counters
  ind->slvr_counter   = 0;
  ind->dadt_counter   = 0;
  ind->jac_counter    = 0;

  op->nlhs           = *nlhsa;
  op->neq            = *neqa;
  op->stiff          = *stiffa;
  
  ind->nBadDose = 0;
  op->do_transit_abs = *transit_abs;

  ind->all_times = timep;
  ind->par_ptr = theta;
  op->inits   = initsp;
  ind->dose    = dosep;
  ind->solve   = retp;
  ind->lhs     = lhsp;
  ind->evid    = evidp;
  ind->rc = rc;
  t_set_solve set_solve = (t_set_solve)(op->set_solve);
  SEXP sd = R_NilValue;
  set_solve(rx);
  par_solve(rx, sd, 0); // Solve without the option of updating residuals.
  t_calc_lhs calc_lhs = (t_calc_lhs)(op->calc_lhs);
  if (*nlhsa) {
    for (i=0; i<*ntime; i++){
      // 0 = first subject; Calc lhs changed...
      calc_lhs(0, timep[i], retp+i*(*neqa), lhsp+i*(*nlhsa));
    }
  }
  Free(idose);
}


void RxODE_ode_solve_env(SEXP sexp_rho){
  Rprintf("1\n");
  int pro = 0;
  SEXP sexp_theta = PROTECT(findVar(installChar(mkChar("params")),sexp_rho));pro++;
  SEXP sexp_inits = PROTECT(findVar(installChar(mkChar("inits")),sexp_rho)); pro++;
  SEXP sexp_lhs   = PROTECT(findVar(installChar(mkChar("lhs_vars")),sexp_rho)); pro++;
  // Events
  SEXP sexp_time = PROTECT(findVar(installChar(mkChar("time")),sexp_rho)); pro++;
  SEXP sexp_evid = PROTECT(findVar(installChar(mkChar("evid")),sexp_rho)); pro++;
  int *evidp = INTEGER(sexp_evid);
  SEXP sexp_dose = PROTECT(findVar(installChar(mkChar("amt")),sexp_rho)); pro++;
  // Covariates
  SEXP sexp_pcov = PROTECT(findVar(installChar(mkChar("pcov")),sexp_rho)); pro++;
  SEXP sexp_cov = PROTECT(findVar(installChar(mkChar("cov")),sexp_rho)); pro++;
  SEXP sexp_locf = PROTECT(findVar(installChar(mkChar("isLocf")),sexp_rho)); pro++;
  // Solver Options
  SEXP sexp_atol = PROTECT(findVar(installChar(mkChar("atol")),sexp_rho)); pro++;
  SEXP sexp_rtol = PROTECT(findVar(installChar(mkChar("rtol")),sexp_rho)); pro++;
  SEXP sexp_hmin = PROTECT(findVar(installChar(mkChar("hmin")),sexp_rho)); pro++;
  SEXP sexp_hmax = PROTECT(findVar(installChar(mkChar("hmax")),sexp_rho)); pro++;
  SEXP sexp_h0 = PROTECT(findVar(installChar(mkChar("hini")),sexp_rho)); pro++;
  SEXP sexp_mxordn = PROTECT(findVar(installChar(mkChar("maxordn")),sexp_rho)); pro++;
  SEXP sexp_mxords = PROTECT(findVar(installChar(mkChar("maxords")),sexp_rho)); pro++;
  SEXP sexp_mx = PROTECT(findVar(installChar(mkChar("maxsteps")),sexp_rho)); pro++;
  SEXP sexp_stiff = PROTECT(findVar(installChar(mkChar("stiff")),sexp_rho)); pro++;
  SEXP sexp_transit_abs = PROTECT(findVar(installChar(mkChar("transit_abs")),sexp_rho)); pro++;
  SEXP sexp_rc = PROTECT(findVar(installChar(mkChar("rc")),sexp_rho)); pro++;
  int *rce    = INTEGER(sexp_rc);
  Rprintf("2\n");
  rx_solve *rx = (rx_solve*)(R_ExternalPtrAddr(VECTOR_ELT(VECTOR_ELT(__mv, 12), 11)));
  rx_solving_options *op = (rx_solving_options*)R_ExternalPtrAddr(rx->op);
  rx_solving_options_ind *inds = rx->subjects;
  rx_solving_options_ind *ind = &inds[0];
  ind->par_ptr = REAL(sexp_theta);
  int neq = op->neq, i = 0;
  Rprintf("3\n");
  double *InfusionRate =ind->InfusionRate,
    *scale = op->scale;
  int *BadDose = ind->BadDose;
  // A bit paranoid -- make sure these are sane...
  for (i = 0; i < neq; i++){
    InfusionRate[i] = 0.0;
    scale[i] = 1.0;
    BadDose[i] = 0;
  }
  Rprintf("4\n");
  // Instead of having the correct length for idose, use idose length = length of ntime
  // Saves an additional for loop at the cost of a little memory.
  int *idose;
  idose = Calloc(length(sexp_time),int);
  ind->idose = idose;
  ind->ndoses=0;
  for (i = 0; i < ind->n_all_times; i++){
    if (evidp[i]){
      ind->ndoses++;
      ind->idose[ind->ndoses-1] = i;
    }
  }
  Rprintf("5\n");
  ind->all_times     = REAL(sexp_time);
  ind->n_all_times   = length(sexp_time);
  ind->evid          = INTEGER(sexp_evid);
  ind->dose          = REAL(sexp_dose);
  // Covariates
  Rprintf("6\n");
  op->par_cov       = INTEGER(sexp_pcov);
  op->do_par_cov    = 1;
  ind->cov_ptr       = REAL(sexp_cov);
  op->ncov          = length(sexp_pcov);
  op->is_locf       = INTEGER(sexp_locf)[0];
  // Solver options
  Rprintf("7\n");
  op->ATOL           = REAL(sexp_atol)[0];
  op->RTOL           = REAL(sexp_rtol)[0];
  op->HMIN           = REAL(sexp_hmin)[0];
  ind->HMAX          = REAL(sexp_hmax)[0];
  op->H0             = REAL(sexp_h0)[0];
  op->MXORDN         = INTEGER(sexp_mxordn)[0];
  op->MXORDS         = INTEGER(sexp_mxords)[0];
  op->mxstep         = INTEGER(sexp_mx)[0];
  op->do_transit_abs = INTEGER(sexp_transit_abs)[0];
  op->stiff          = INTEGER(sexp_stiff)[0];
  ind->slvr_counter   = 0;
  ind->dadt_counter   = 0;
  ind->jac_counter    = 0;
  // LOCF
  Rprintf("8\n");
  if (op->is_locf == 1){
    op->f2 = 0.0; //= f=0 
    op->f1 = 1.0; // = 1-f = 1;
    op->kind = 0;
  } else if (op->is_locf == 2) {
    // NOCB
    op->f2 = 1.0; //= f=1
    op->f1 = 0.0;
    op->kind = 0;
  } else if (op->is_locf == 3){
    op->f2 = 0.5; //= f=0.5
    op->f1 = 0.5;
    op->kind = 0;
  } else {
    // Linear
    op->f2 = 1.0; //= f=0
    op->f1 = 0.0;
    op->kind = 1;
  }
  op->nlhs          = length(sexp_lhs);
  op->neq           = length(sexp_inits);
  op->inits   = REAL(sexp_inits);
  ind->rc = rce;
  
   // Let R handle deallocating the solve and lhs expressions; Should disappear with evironment
  Rprintf("9\n");
  SEXP sexp_solve = PROTECT(allocVector(REALSXP,ind->n_all_times*op->neq)); pro++;
  SEXP sexp_lhsV = PROTECT(allocVector(REALSXP,ind->n_all_times*op->nlhs)); pro++;
  ind->solve = REAL(sexp_solve);
  ind->lhs     = REAL(sexp_lhsV);
  Rprintf("10\n");
  defineVar(install(".solve"), sexp_solve, sexp_rho);
  defineVar(install(".lhs"), sexp_lhsV, sexp_rho);
  Rprintf("11\n");
  t_set_solve set_solve = (t_set_solve)(op->set_solve);
  SEXP sd = R_NilValue;
  Rprintf("12\n");
  set_solve(rx);
  par_solve(rx, sd, 0); // Solve without the option of updating residuals.
  Rprintf("13\n");
  Free(idose);
  UNPROTECT(pro);
}

extern double RxODE_sum(double *input, int len){
  return PreciseSums_sum(input, len);
}

extern double RxODE_sumV(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = PreciseSums_sum(p, n);
  Free(p);
  return s;
}

extern double RxODE_prod(double *input, int len){
  return PreciseSums_prod(input, len);
}

extern double RxODE_prodV(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = PreciseSums_prod(p, n);
  Free(p);
  return s;
}
