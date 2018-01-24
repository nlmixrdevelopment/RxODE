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
extern void calc_lhs(int cSub, double t, double *A, double *lhs);

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

void rxAddModelLib(SEXP);

SEXP __mv;
extern void RxODE_assign_fn_pointers_(SEXP mv, int addit){
  __mv = mv;
  if (addit){
    rxAddModelLib(mv);
  }
}

extern void RxODE_assign_fn_pointers(SEXP mv){
  RxODE_assign_fn_pointers_(mv, 1);
}


int rxIsC(SEXP obj, const char *cls);
extern SEXP RxODE_get_mv(){
  /* if (!rxIsC(__mv,"rxModelVars")){ */
  /*   error("RxODE C functions were not setup correctly."); */
  /* } */
  return __mv;
}

/* extern void rxode_assign_rx(rx_solve *rx); */

extern void rxode_assign_rx(rx_solve *rx);
extern rx_solve *rxSingle(SEXP object, const int stiff,const int transit_abs,
			  const double atol, const double rtol, const int maxsteps,
			  const double hmin, const double hini, const int maxordn,
			  const int maxords, const int cores, const int ncov,
			  int *par_cov, int do_par_cov, 
			  int is_locf,
			  // Other single solve option
			  double hmax, double *par,
			  double *amt, double *solve, double *lhs,
			  int *evid, int *rc, double *cov,
			  int nTimes, double *all_times);

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
  SEXP mv = RxODE_get_mv();
  rx_solve *rx;
  int par_cov[0];
  double cov[0];
  int i =0;
  rx = rxSingle(mv, *stiffa,*transit_abs, *atol, *rtol, 5000,//maxsteps
                0, 0, 12, 5, 1, 0, par_cov, 0, 0, 0, theta,
                dosep, retp, lhsp, evidp, rc, cov, *ntime,timep);
  rxode_assign_rx(rx);
  SEXP sd = R_NilValue;
  par_solve(rx, sd, 0); // Solve without the option of updating residuals.
  if (*nlhsa) {
    for (i=0; i<*ntime; i++){
      // 0 = first subject; Calc lhs changed...
      calc_lhs(0, timep[i], retp+i*(*neqa), lhsp+i*(*nlhsa));
    }
  }
}

#define aexists(a, env) if (Rf_findVarInFrame(env, Rf_install(a)) == R_UnboundValue){ error("need '%s' in environment for solving.",a);}
void RxODE_ode_solve_env(SEXP sexp_rho){
  if(!isEnvironment(sexp_rho)){
    error("Calling RxODE_ode_solve_env without an environment...");
  }
  int pro = 0, i = 0;
  SEXP mv = RxODE_get_mv();
  RxODE_assign_fn_pointers(mv);
  aexists("params",sexp_rho);
  SEXP sexp_theta = PROTECT(findVar(installChar(mkChar("params")),sexp_rho));pro++;
  aexists("inits",sexp_rho);
  SEXP sexp_inits = PROTECT(findVar(installChar(mkChar("inits")),sexp_rho)); pro++;
  aexists("lhs_vars",sexp_rho);
  SEXP sexp_lhs   = PROTECT(findVar(installChar(mkChar("lhs_vars")),sexp_rho)); pro++;
  // Events
  aexists("time",sexp_rho);
  SEXP sexp_time = PROTECT(findVar(installChar(mkChar("time")),sexp_rho)); pro++;
  aexists("evid",sexp_rho);
  SEXP sexp_evid = PROTECT(findVar(installChar(mkChar("evid")),sexp_rho)); pro++;
  aexists("amt",sexp_rho);
  SEXP sexp_dose = PROTECT(findVar(installChar(mkChar("amt")),sexp_rho)); pro++;
  // Covariates
  aexists("pcov",sexp_rho);
  SEXP sexp_pcov = PROTECT(findVar(installChar(mkChar("pcov")),sexp_rho)); pro++;
  aexists("covs",sexp_rho);
  SEXP sexp_cov = PROTECT(findVar(installChar(mkChar("covs")),sexp_rho)); pro++;
  aexists("isLocf",sexp_rho);
  SEXP sexp_locf = PROTECT(findVar(installChar(mkChar("isLocf")),sexp_rho)); pro++;
  // Solver Options
  aexists("atol",sexp_rho);
  SEXP sexp_atol = PROTECT(findVar(installChar(mkChar("atol")),sexp_rho)); pro++;
  aexists("rtol",sexp_rho);
  SEXP sexp_rtol = PROTECT(findVar(installChar(mkChar("rtol")),sexp_rho)); pro++;
  aexists("hmin",sexp_rho);
  SEXP sexp_hmin = PROTECT(findVar(installChar(mkChar("hmin")),sexp_rho)); pro++;
  aexists("hmax",sexp_rho);
  SEXP sexp_hmax = PROTECT(findVar(installChar(mkChar("hmax")),sexp_rho)); pro++;
  aexists("hini",sexp_rho);
  SEXP sexp_h0 = PROTECT(findVar(installChar(mkChar("hini")),sexp_rho)); pro++;
  aexists("maxordn",sexp_rho);
  SEXP sexp_mxordn = PROTECT(findVar(installChar(mkChar("maxordn")),sexp_rho)); pro++;
  aexists("maxords",sexp_rho);
  SEXP sexp_mxords = PROTECT(findVar(installChar(mkChar("maxords")),sexp_rho)); pro++;
  aexists("maxsteps",sexp_rho);
  SEXP sexp_mx = PROTECT(findVar(installChar(mkChar("maxsteps")),sexp_rho)); pro++;
  aexists("stiff",sexp_rho);
  SEXP sexp_stiff = PROTECT(findVar(installChar(mkChar("stiff")),sexp_rho)); pro++;
  aexists("transit_abs",sexp_rho);
  SEXP sexp_transit_abs = PROTECT(findVar(installChar(mkChar("transit_abs")),sexp_rho)); pro++;
  aexists("rc",sexp_rho);
  SEXP sexp_rc = PROTECT(findVar(installChar(mkChar("rc")),sexp_rho)); pro++;
  int *rce    = INTEGER(sexp_rc);  
  // Let R handle deallocating the solve and lhs expressions; Should disappear with evironment
  SEXP sexp_solve = PROTECT(allocVector(REALSXP,length(sexp_time)*length(sexp_inits))); pro++;
  double *solve = REAL(sexp_solve);
  for (i = 0; i < length(sexp_time)*length(sexp_inits); i++){
    solve[i] = 0;
  }
  SEXP sexp_lhsV = PROTECT(allocVector(REALSXP,length(sexp_time)*length(sexp_lhs))); pro++;
  double *lhs = REAL(sexp_lhsV);
  for (i = 0; i < length(sexp_time)*length(sexp_lhs); i++){
    lhs[i] = 0;
  }
  int stiff = INTEGER(sexp_stiff)[0];
  int transit_abs = INTEGER(sexp_transit_abs)[0];
  double atol = REAL(sexp_atol)[0];
  double rtol= REAL(sexp_rtol)[0];
  int mx = INTEGER(sexp_mx)[0];
  double hmin = REAL(sexp_hmin)[0];
  double h0 = REAL(sexp_h0)[0];
  int mxordn= INTEGER(sexp_mxordn)[0];
  int mxords = INTEGER(sexp_mxords)[0];
  int *pcov = INTEGER(sexp_pcov);
  int locf = INTEGER(sexp_locf)[0];
  double hmax = REAL(sexp_hmax)[0];
  double *theta = REAL(sexp_theta);
  double *dose = REAL(sexp_dose);
  int *evid = INTEGER(sexp_evid);
  double *cov = REAL(sexp_cov);
  double *time = REAL(sexp_time);
  rx_solve *rx = rxSingle(mv, stiff, transit_abs, atol, rtol, mx, hmin, h0,  mxordn,
		          mxords, 1, length(sexp_pcov), pcov, 1,  locf,
			  hmax, theta, dose, solve, lhs, evid, rce, cov,
			  length(sexp_time), time);
  rxode_assign_rx(rx);
  SEXP sd = R_NilValue;
  par_solve(rx, sd, 0); // Solve without the option of updating residuals.
  defineVar(install(".lhs"), sexp_lhsV, sexp_rho);
  defineVar(install(".solve"), sexp_solve, sexp_rho);
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
