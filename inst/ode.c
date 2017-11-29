#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#define JAC_Rprintf Rprintf
#define JAC0_Rprintf if (_jac_counter_val() == 0) Rprintf
#define ODE_Rprintf Rprintf
#define ODE0_Rprintf if (_dadt_counter_val() == 0) Rprintf
#define LHS_Rprintf Rprintf
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#define R_pow Rx_pow
#define R_pow_di Rx_pow_di

// Types for par pointers.r
typedef void (*RxODE_update_par_ptr)(double t);
typedef double (*RxODE_transit3)(double t, double n, double mtt);
typedef double (*RxODE_fn) (double x);
typedef double (*RxODE_fn2) (double x, double y);
typedef double (*RxODE_fn2i) (double x, int i);
typedef double (*RxODE_transit4)(double t, double n, double mtt, double bio);
typedef double (*RxODE_vec) (int val);
typedef long (*RxODE_cnt) ();
typedef void (*RxODE_inc) ();
typedef double (*RxODE_val) ();
typedef SEXP (*RxODE_ode_solver) (SEXP sexp_theta, SEXP sexp_inits, SEXP sexp_lhs, SEXP sexp_time, SEXP sexp_evid,SEXP sexp_dose, SEXP sexp_pcov, SEXP sexp_cov, SEXP sexp_locf, SEXP sexp_atol, SEXP sexp_rtol, SEXP sexp_hmin, SEXP sexp_hmax, SEXP sexp_h0, SEXP sexp_mxordn, SEXP sexp_mxords, SEXP sexp_mx,SEXP sexp_stiff, SEXP sexp_transit_abs, SEXP sexp_object, SEXP sexp_extra_args, SEXP sexp_matrix, SEXP sexp_add_cov);
typedef void (*RxODE_assign_fn_pointers)(void (*fun_dydt)(unsigned int, double, double *, double *),void (*fun_calc_lhs)(double, double *, double *),void (*fun_calc_jac)(unsigned int, double, double *, double *, unsigned int),void (*fun_update_inis)(SEXP _ini_sexp),int fun_jt,int fun_mf, int fun_debug);
typedef SEXP (*RxODE_assign_fn_xpointers)(void (*fun_dydt)(unsigned int, double, double *, double *),void (*fun_calc_lhs)(double, double *, double *),void (*fun_calc_jac)(unsigned int, double, double *, double *, unsigned int),void (*fun_update_inis)(SEXP _ini_sexp),int fun_jt,int fun_mf, int fun_debug);

typedef void (*RxODE_ode_solver_old_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc);
typedef void (*RxODE_ode_solver_0_6_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc,double hmin, double hmax,double h0,int mxordn,int mxords,int mxstep);
typedef double (*RxODE_solveLinB)(double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);
typedef double (*RxODE_sum_prod)(double *input, int n);
// Give par pointers
RxODE_vec _par_ptr, _InfusionRate;
RxODE_update_par_ptr _update_par_ptr;
RxODE_cnt _dadt_counter_val, _jac_counter_val;
RxODE_inc _dadt_counter_inc, _jac_counter_inc;
RxODE_val podo, tlast;
RxODE_transit4 _transit4;
RxODE_transit3 _transit3;
RxODE_fn _safe_log, safe_zero, factorial, _as_zero, abs_log, abs_log1p;
RxODE_fn2 sign_exp, Rx_pow;
RxODE_fn2i Rx_pow_di;
RxODE_assign_fn_pointers _assign_fn_pointers;
RxODE_assign_fn_xpointers _assign_fn_xpointers;
RxODE_ode_solver_old_c _old_c;
RxODE_ode_solver_0_6_c _c_0_6;
RxODE_solveLinB solveLinB;
RxODE_sum_prod _sum1, _prod1;

double _sum(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = _sum1(p, n);
  Free(p);
  return s;
}

double _prod(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = _prod1(p, n);
  Free(p);
  return s;
}

extern double _sign(unsigned int n, ...){
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

extern void __ODE_SOLVER_PTR__();
extern SEXP __ODE_SOLVER_XPTR__();

extern void __ODE_SOLVER__(
                    int *neq,
                    double *theta,      //order:
                    double *time,
                    int *evid,
                    int *ntime,
                    double *inits,
                    double *dose,
                    double *ret,
                    double *atol,
                    double *rtol,
                    int *stiff,
                    int *transit_abs,
                    int *nlhs,
                    double *lhs,
                    int *rc
                    ){
  // Backward compatible ode solver for 0.5* C interface
  __ODE_SOLVER_PTR__();
  _old_c(neq, theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc);
}

void __ODE_SOLVER_0_6__(int *neq,
                        double *theta,  //order:
                        double *time,
                        int *evid,
                        int *ntime,
                        double *inits,
                        double *dose,
                        double *ret,
                        double *atol,
                        double *rtol,
                        int *stiff,
                        int *transit_abs,
                        int *nlhs,
                        double *lhs,
                        int *rc,
                        double hmin,
                        double hmax,
                        double h0,
                        int mxordn,
                        int mxords,
                        int mxstep) {
  // Backward compatible ode solver for 0.5* C interface
  __ODE_SOLVER_PTR__();
  _c_0_6(neq, theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc,
	hmin, hmax, h0, mxordn, mxords, mxstep);
}

extern void __ODE_SOLVER_PTR__  (){
  _assign_fn_pointers(__DYDT__ , __CALC_LHS__ , __CALC_JAC__, __INIS__, __JT__ , __MF__,
#ifdef __DEBUG__
                      1
#else
                      0
#endif
                      );
}

extern SEXP __ODE_SOLVER_XPTR__  (){
  return _assign_fn_xpointers(__DYDT__ , __CALC_LHS__ , __CALC_JAC__, __INIS__, __JT__ , __MF__,
#ifdef __DEBUG__
                      1
#else
                      0
#endif
                      );
}

extern SEXP __ODE_SOLVER_SEXP__ (// Parameters
                                 SEXP sexp_theta,
                                 SEXP sexp_inits,
                                 SEXP sexp_lhs,
				 // Events
				 SEXP sexp_time,
				 SEXP sexp_evid,
				 SEXP sexp_dose,
				 // Covariates
				 SEXP sexp_pcov,
				 SEXP sexp_cov,
				 SEXP sexp_locf,
				 // Solver Options
				 SEXP sexp_atol,
				 SEXP sexp_rtol,
				 SEXP sexp_hmin,
				 SEXP sexp_hmax,
				 SEXP sexp_h0,
				 SEXP sexp_mxordn,
				 SEXP sexp_mxords,
				 SEXP sexp_mx,
				 SEXP sexp_stiff,
				 SEXP sexp_transit_abs,
				 // Object Creation
				 SEXP sexp_object,
				 SEXP sexp_extra_args,
				 SEXP sexp_matrix,
				 SEXP sexp_add_cov){
  __ODE_SOLVER_PTR__();
  RxODE_ode_solver ode_solver = (RxODE_ode_solver) R_GetCCallable("RxODE","RxODE_ode_solver");
  return ode_solver(sexp_theta,sexp_inits,sexp_lhs,sexp_time,sexp_evid,sexp_dose,sexp_pcov,sexp_cov,sexp_locf,sexp_atol,
		    sexp_rtol,sexp_hmin,sexp_hmax,sexp_h0,sexp_mxordn,sexp_mxords,sexp_mx,sexp_stiff,sexp_transit_abs,
		    sexp_object,sexp_extra_args,sexp_matrix,sexp_add_cov);
}

static R_NativePrimitiveArgType __ODE_SOLVER__rx_t[] = {
  //*neq, *theta, *time,  *evid, *ntime, *inits,   *dose,   *ret,     *atol,  *rtol,   *stiff, *transit_abs, *nlhs, *lhs, *rc
  INTSXP,REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP
};

static R_NativePrimitiveArgType  __ODE_SOLVER_0_6__rx_t[] = {
  //*neq, *theta, *time,  *evid, *ntime, *inits,   *dose,   *ret,     *atol,  *rtol,   *stiff, *transit_abs, *nlhs, *lhs, *rc
  INTSXP,REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

//Initilize the dll to match RxODE's calls
void __R_INIT__ (DllInfo *info){
  // Get the RxODE calling interfaces
  _InfusionRate   = (RxODE_vec) R_GetCCallable("RxODE","RxODE_InfusionRate");
  _update_par_ptr = (RxODE_update_par_ptr) R_GetCCallable("RxODE","RxODE_update_par_ptr");
  _par_ptr = (RxODE_vec) R_GetCCallable("RxODE","RxODE_par_ptr");
  _dadt_counter_val = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_dadt_counter_val");
  _jac_counter_val  = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_jac_counter_val");
  _dadt_counter_inc = (RxODE_inc) R_GetCCallable("RxODE","RxODE_dadt_counter_inc");
  _jac_counter_inc  = (RxODE_inc) R_GetCCallable("RxODE","RxODE_jac_counter_inc");
  podo  = (RxODE_val) R_GetCCallable("RxODE","RxODE_podo");
  tlast = (RxODE_val) R_GetCCallable("RxODE","RxODE_tlast");
  factorial=(RxODE_fn) R_GetCCallable("RxODE","RxODE_factorial");
  _transit3 = (RxODE_transit3) R_GetCCallable("RxODE","RxODE_transit3");
  _transit4 = (RxODE_transit4) R_GetCCallable("RxODE","RxODE_transit4");
  _safe_log =(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_log");
  safe_zero =(RxODE_fn) R_GetCCallable("RxODE","RxODE_safe_zero");
  _as_zero = (RxODE_fn) R_GetCCallable("RxODE","RxODE_as_zero");
  _assign_fn_pointers=(RxODE_assign_fn_pointers) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
  _assign_fn_xpointers=(RxODE_assign_fn_xpointers) R_GetCCallable("RxODE","RxODE_get_fn_pointers");
  _old_c = (RxODE_ode_solver_old_c) R_GetCCallable("RxODE","RxODE_ode_solver_old_c");
  _c_0_6 = (RxODE_ode_solver_0_6_c)R_GetCCallable("RxODE","RxODE_ode_solver_0_6_c");
  _sum1   = (RxODE_sum_prod)R_GetCCallable("RxODE","RxODE_sum");
  _prod1 = (RxODE_sum_prod) R_GetCCallable("RxODE","RxODE_prod");
  sign_exp = (RxODE_fn2) R_GetCCallable("RxODE","RxODE_sign_exp");
  Rx_pow = (RxODE_fn2) R_GetCCallable("RxODE","RxODE_pow");
  Rx_pow_di = (RxODE_fn2i) R_GetCCallable("RxODE","RxODE_pow_di");
  abs_log = (RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log");
  abs_log1p = (RxODE_fn) R_GetCCallable("RxODE","RxODE_abs_log1p");
  solveLinB = (RxODE_solveLinB) R_GetCCallable("RxODE","RxODE_solveLinB");
  // Register the outside functions
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_STR__,       (DL_FUNC) __ODE_SOLVER__);
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_SEXP_STR__,  (DL_FUNC) __ODE_SOLVER_SEXP__);
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_0_6_STR__,   (DL_FUNC) __ODE_SOLVER_0_6__);
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_PTR_STR__,   (DL_FUNC) __ODE_SOLVER_PTR__);

  /* R_registerRoutines(info, NULL, NULL, NULL, NULL); */
  /* R_useDynamicSymbols(info,TRUE); */

  static const R_CMethodDef cMethods[] = {
    {__ODE_SOLVER_STR__, (DL_FUNC) &__ODE_SOLVER__, 15, __ODE_SOLVER__rx_t},
    {__ODE_SOLVER_0_6_STR__, (DL_FUNC) &__ODE_SOLVER_0_6__, 21, __ODE_SOLVER_0_6__rx_t},
    {NULL, NULL, 0, NULL}
  };
  
  R_CallMethodDef callMethods[]  = {
    {__ODE_SOLVER_XPTR_STR__, (DL_FUNC) &__ODE_SOLVER_XPTR__, 0},
    {__ODE_SOLVER_PTR_STR__, (DL_FUNC) &__ODE_SOLVER_PTR__, 0},
    {__ODE_SOLVER_SEXP_STR__, (DL_FUNC) &__ODE_SOLVER_SEXP__, 23},
    {__MODEL_VARS_STR__, (DL_FUNC) &__MODEL_VARS__, 0},
    {"__INIS__", (DL_FUNC) &__INIS__, 1},
    {NULL, NULL, 0}
  };
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
  // Register the function pointers so if someone directly calls the
  // ode solvers directly, they use the last loaded RxODE model.
  __ODE_SOLVER_PTR__();
}

