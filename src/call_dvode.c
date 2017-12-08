#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "dop853.h"
#define NCMT 100
#define max(a, b) ((a) > (b) ? (a) : (b))
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>
#include <PreciseSums.h>

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
void RxODE_ode_solver_old_c(int *neqa,
                            double *theta,  //order:
                            double *time,
                            int *evidp,
                            int *ntime,
                            double *initsp,
                            double *dosep,
                            double *ret,
                            double *atol,
                            double *rtol,
                            int *stiffa,
                            int *transit_abs,
                            int *nlhsa,
                            double *lhsp,
                            int *rc){
  /* if (*neqa > NCMT){ */
  /*   error("RxODE does not support %d compartments (Currently only %d compartments)", neq, NCMT); */
  /* } */
  /* int i; */
  /* for (i=0; i< *neqa; i++) InfusionRate[i] = 0.0; */
  /* ndoses = -1; */
  /* all_times         = time; */
  /* n_all_times       = *ntime; */
  /* //RxODE_ode_dosing_ = RxODE_ode_get_dosing(); */
  /* // par_cov */
  /* do_par_cov        = 0; */
  /* // cov_ptr */
  /* ncov              = 0; */
  /* is_locf           = 0; */
  /* // Solver Options */
  /* ATOL = *atol; */
  /* RTOL = *rtol; */
  /* // Assign to default LSODA behvior, or 0 */
  /* HMIN           = 0; */
  /* HMAX           = 0; */
  /* H0             = 0; */
  /* MXORDN         = 0; */
  /* MXORDS         = 0; */
  /* mxstep         = 5000; // Not LSODA default but RxODE default */
  /* // Counters */
  /* slvr_counter   = 0; */
  /* dadt_counter   = 0; */
  /* jac_counter    = 0; */

  /* nlhs           = *nlhsa; */
  /* neq            = *neqa; */
  /* stiff          = *stiffa; */
  
  

  /* nBadDose = 0; */
  /* do_transit_abs = *transit_abs; */
  
  /* par_ptr = theta; */
  /* inits   = initsp; */
  /* dose    = dosep; */
  /* solve   = ret; */
  /* lhs     = lhsp; */
  /* evid    = evidp; */

  /* // Assign global time information */
  /* // Call solver */
  
  
  /* /\* Rprintf("Call Solver; par_ptr[0] = %f; evid[0]=%d; inits[0]=%d\n",par_ptr[0],evid[0],inits[0]); *\/ */
  /* RxODE_ode_solver_c(*neqa, *stiffa, evidp, initsp, dosep, ret, rc); */
  /* /\* Rprintf("Update LHS\n"); *\/ */
  /* // Update LHS */
  /* if (*nlhsa) { */
  /*   for (i=0; i<*ntime; i++){ */
  /*     calc_lhs(0, time[i], ret+i*(*neqa), lhsp+i*(*nlhsa)); */
  /*   } */
  /* } */
}


extern void RxODE_ode_free(){
}

void RxODE_ode_alloc(){
}

extern void RxODE_assign_fn_pointers(void (*fun_dydt)(int*, double, double *, double *),
				     void (*fun_calc_lhs)(int, double, double *, double *),
				     void (*fun_calc_jac)(int, double, double *, double *, unsigned int),
				     void (*fun_update_inis)(int, double *),
				     int fun_jt,
				     int fun_mf,
				     int fun_debug){
  // Assign functions pointers
  /* dydt     = fun_dydt; */
  /* calc_jac = fun_calc_jac; */
  /* calc_lhs = fun_calc_lhs; */
  /* update_inis = fun_update_inis; */
  /* // Assign solver options */
  /* global_jt     = fun_jt; */
  /* global_mf     = fun_mf; */
  /* global_debug  = fun_debug; */
}

void RxODE_ode_solve_env(SEXP sexp_rho){
  error("now.");
  /* int pro = 0; */
  /* SEXP sexp_theta = PROTECT(findVar(installChar(mkChar("params")),sexp_rho));pro++; */
  /* SEXP sexp_inits = PROTECT(findVar(installChar(mkChar("inits")),sexp_rho)); pro++; */
  /* SEXP sexp_lhs   = PROTECT(findVar(installChar(mkChar("lhs_vars")),sexp_rho)); pro++; */
  /* // Events */
  /* SEXP sexp_time = PROTECT(findVar(installChar(mkChar("time")),sexp_rho)); pro++; */
  /* SEXP sexp_evid = PROTECT(findVar(installChar(mkChar("evid")),sexp_rho)); pro++; */
  /* SEXP sexp_dose = PROTECT(findVar(installChar(mkChar("amt")),sexp_rho)); pro++; */
  /* // Covariates */
  /* SEXP sexp_pcov = PROTECT(findVar(installChar(mkChar("pcov")),sexp_rho)); pro++; */
  /* SEXP sexp_cov = PROTECT(findVar(installChar(mkChar("cov")),sexp_rho)); pro++; */
  /* SEXP sexp_locf = PROTECT(findVar(installChar(mkChar("isLocf")),sexp_rho)); pro++; */
  /* // Solver Options */
  /* SEXP sexp_atol = PROTECT(findVar(installChar(mkChar("atol")),sexp_rho)); pro++; */
  /* SEXP sexp_rtol = PROTECT(findVar(installChar(mkChar("rtol")),sexp_rho)); pro++; */
  /* SEXP sexp_hmin = PROTECT(findVar(installChar(mkChar("hmin")),sexp_rho)); pro++; */
  /* SEXP sexp_hmax = PROTECT(findVar(installChar(mkChar("hmax")),sexp_rho)); pro++; */
  /* SEXP sexp_h0 = PROTECT(findVar(installChar(mkChar("hini")),sexp_rho)); pro++; */
  /* SEXP sexp_mxordn = PROTECT(findVar(installChar(mkChar("maxordn")),sexp_rho)); pro++; */
  /* SEXP sexp_mxords = PROTECT(findVar(installChar(mkChar("maxords")),sexp_rho)); pro++; */
  /* SEXP sexp_mx = PROTECT(findVar(installChar(mkChar("maxsteps")),sexp_rho)); pro++; */
  /* SEXP sexp_stiff = PROTECT(findVar(installChar(mkChar("stiff")),sexp_rho)); pro++; */
  /* SEXP sexp_transit_abs = PROTECT(findVar(installChar(mkChar("transit_abs")),sexp_rho)); pro++; */
  /* SEXP sexp_rc = PROTECT(findVar(installChar(mkChar("rc")),sexp_rho)); pro++; */
  /* int *rce    = INTEGER(sexp_rc); */

  /* par_ptr       = REAL(sexp_theta); */
  /* inits         = REAL(sexp_inits); */

  /* RxODE_ode_setup(sexp_inits, sexp_lhs, sexp_time, sexp_evid, sexp_dose, sexp_pcov, sexp_cov, */
  /* 		  sexp_locf, sexp_atol, sexp_rtol, sexp_hmin, sexp_hmax, sexp_h0, sexp_mxordn, */
  /* 		  sexp_mxords, sexp_mx, sexp_stiff, sexp_transit_abs); */
  /* RxODE_ode_alloc(); */
  /* RxODE_ode_solver_c(neq, stiff, evid, inits, dose, solve, rc); */
  /* // Send rc to environment */
  /* rce[0] = rc[0]; */
  /* UNPROTECT(pro); */
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
