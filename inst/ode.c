#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#define JAC_Rprintf Rprintf
#define JAC0_Rprintf if (jac_counter_val() == 0) Rprintf
#define ODE_Rprintf Rprintf
#define ODE0_Rprintf if (dadt_counter_val() == 0) Rprintf
#define LHS_Rprintf Rprintf
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
// Types for par pointers.r
typedef void (*RxODE_update_par_ptr)(double t);
typedef double (*RxODE_transit3)(double t, double n, double mtt);
typedef double (*RxODE_fn) (double x);
typedef double (*RxODE_transit4)(double t, double n, double mtt, double bio);
typedef double (*RxODE_vec) (int val);
typedef long (*RxODE_cnt) ();
typedef void (*RxODE_inc) ();
typedef double (*RxODE_val) ();
typedef SEXP (*RxODE_ode_solver) (SEXP sexp_theta, SEXP sexp_inits, SEXP sexp_lhs, SEXP sexp_time, SEXP sexp_evid,SEXP sexp_dose, SEXP sexp_pcov, SEXP sexp_cov, SEXP sexp_locf, SEXP sexp_atol, SEXP sexp_rtol, SEXP sexp_hmin, SEXP sexp_hmax, SEXP sexp_h0, SEXP sexp_mxordn, SEXP sexp_mxords, SEXP sexp_mx,SEXP sexp_stiff, SEXP sexp_transit_abs, SEXP sexp_object, SEXP sexp_extra_args,void (*fun_dydt)(unsigned int, double, double *, double *),void (*fun_calc_lhs)(double, double *, double *),void (*fun_calc_jac)(unsigned int, double, double *, double *, unsigned int),int fun_jt, int fun_mf, int fun_debug);
// Give par pointers
RxODE_vec par_ptr, InfusionRate;
RxODE_update_par_ptr update_par_ptr;
RxODE_cnt dadt_counter_val, jac_counter_val;
RxODE_inc dadt_counter_inc, jac_counter_inc;
RxODE_val podo, tlast;
RxODE_transit4 transit4;
RxODE_transit3 transit3;
RxODE_fn factorial;
//Initilize the dll to match RxODE's calls
void __R_INIT__ (DllInfo *info){
  InfusionRate   = (RxODE_vec) R_GetCCallable("RxODE","RxODE_InfusionRate");
  update_par_ptr = (RxODE_update_par_ptr) R_GetCCallable("RxODE","RxODE_update_par_ptr");
  par_ptr = (RxODE_vec) R_GetCCallable("RxODE","RxODE_par_ptr");
  dadt_counter_val = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_dadt_counter_val");
  jac_counter_val  = (RxODE_cnt) R_GetCCallable("RxODE","RxODE_jac_counter_val");
  dadt_counter_inc = (RxODE_inc) R_GetCCallable("RxODE","RxODE_dadt_counter_inc");
  jac_counter_inc  = (RxODE_inc) R_GetCCallable("RxODE","RxODE_jac_counter_inc");
  podo  = (RxODE_val) R_GetCCallable("RxODE","RxODE_podo");
  tlast = (RxODE_val) R_GetCCallable("RxODE","RxODE_tlast");
  transit3 = (RxODE_transit3) R_GetCCallable("RxODE","RxODE_transit3");
  transit4 = (RxODE_transit4) R_GetCCallable("RxODE","RxODE_transit4");
  factorial=(RxODE_fn) R_GetCCallable("RxODE","RxODE_factorial");
}


extern SEXP __ODE_SOLVER__ (// Parameters
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
			    SEXP sexp_extra_args){
  RxODE_ode_solver ode_solver=(RxODE_ode_solver) R_GetCCallable("RxODE","RxODE_ode_solver");
  ode_solver(sexp_theta,sexp_inits,sexp_lhs,sexp_time,sexp_evid,sexp_dose,sexp_pcov,sexp_cov,sexp_locf,sexp_atol,
	     sexp_rtol,sexp_hmin,sexp_hmax,sexp_h0,sexp_mxordn,sexp_mxords,sexp_mx,sexp_stiff,sexp_transit_abs,
	     sexp_object,sexp_extra_args, __DYDT__ , __CALC_LHS__ , __CALC_JAC__, __JT__ , __MF__,
#ifdef __DEBUG__
	     1
#else
	     0
#endif
	     );
}
