#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

static R_NativePrimitiveArgType RxODE_sign_exp_t[] = {
  REALSXP, REALSXP
};

static R_NativePrimitiveArgType RxODE_one_int_t[] = {
  INTSXP
};

static R_NativePrimitiveArgType RxODE_transit4_t[] = {
  REALSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType RxODE_transit3_t[] = {
  REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType RxODE_one_dbl_t[] = {
  REALSXP
};

SEXP RxODE_ode_solver(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP trans(SEXP orig_file, SEXP parse_file, SEXP c_file, SEXP extra_c, SEXP prefix, SEXP model_md5, SEXP parse_model,SEXP parse_model3);
SEXP _RxODE_linCmtEnv(SEXP rho);
SEXP _RxODE_rxInv(SEXP matrix);
SEXP _RxODE_removableDrive(SEXP letter);
SEXP _RxODE_rxCoutEcho(SEXP number);
SEXP _RxODE_RxODE_finalize_focei_omega(SEXP);
SEXP _RxODE_RxODE_finalize_log_det_OMGAinv_5(SEXP);
SEXP _rxCholInv(SEXP dms, SEXP theta, SEXP tn);
SEXP _RxODE_rxSymInvCholEnvCalculate(SEXP, SEXP, SEXP);
SEXP _RxODE_rxInvWishartVar(SEXP, SEXP);
SEXP _RxODE_rxSymInvChol(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxDataSetup(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxIs(SEXP,SEXP);
SEXP _RxODE_rxModelVars(SEXP);
SEXP _RxODE_rxState(SEXP, SEXP);
SEXP _RxODE_rxParams(SEXP);
SEXP _RxODE_rxDfdy(SEXP);
SEXP _RxODE_rxLhs(SEXP);
SEXP _RxODE_rxInits(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxUpdateResiduals(SEXP);
SEXP _RxODE_rxSetupIni(SEXP, SEXP);
SEXP _RxODE_rxDataParSetup(SEXP, SEXP, SEXP, SEXP, SEXP,
                           SEXP, SEXP, SEXP, SEXP, SEXP,
                           SEXP);
SEXP _RxODE_rxSolvingOptions(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxSolvingData(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxData(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,
		   SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,
		   SEXP,SEXP,SEXP);

double RxODE_solveLinB(double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);
static R_NativePrimitiveArgType RxODE_solveLinB_t[] = {
  //t,    linCmt,  diff1,  diff2,  A,       alpha,  B,       beta,     C,       gamma, double ka, double tlag)
  REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

SEXP _RxODE_rxToOmega(SEXP cholInv);

static R_NativePrimitiveArgType RxODE_Sum_t[] = {
  REALSXP, INTSXP
};

extern double RxODE_as_zero(double x);
extern double RxODE_safe_log(double x);
extern double RxODE_safe_zero(double x);
extern double RxODE_pow(double x, double y);
extern double RxODE_pow_di(double x, int i);
extern double RxODE_sign_exp(double sgn, double x);
extern double RxODE_abs_log(double x);
extern double RxODE_abs_log1p(double x);
extern double RxODE_factorial(double x);
extern double RxODE_transit4(double t, double n, double mtt, double bio);
extern double RxODE_transit3(double t, double n, double mtt);
extern double RxODE_sum(double *input, int len);
extern double RxODE_prod(double *input, int len);
extern void RxODE_ode_solve_env(SEXP sexp_rho);
extern int nEq ();
extern unsigned int nObs();
extern unsigned int nLhs ();
extern double RxODE_as_zero(double x);
extern double rxLhs(int i);
extern void rxCalcLhs(int i);
extern unsigned int nAllTimes ();
extern int rxEvid(int i);

extern SEXP RxODE_get_fn_pointers(void (*fun_dydt)(unsigned int, double, double *, double *),
                                  void (*fun_calc_lhs)(double, double *, double *),
                                  void (*fun_calc_jac)(unsigned int, double, double *, double *, unsigned int),
                                  void (*fun_update_inis)(SEXP _ini_sexp),
                                  void (*fun_dydt_lsoda_dum)(int *, double *, double *, double *),
                                  void (*fun_jdum_lsoda)(int *, double *, double *,int *, int *, double *, int *),
                                  int fun_jt,
                                  int fun_mf,
                                  int fun_debug);
extern void RxODE_ode_solver_old_c(int *neqa,
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
				   int *rc);

// Need to change to remove global variables
extern double RxODE_InfusionRate(int val);
extern double RxODE_par_ptr(int val);
extern long RxODE_jac_counter_val();
extern long RxODE_dadt_counter_val();
extern void RxODE_jac_counter_inc();
extern void RxODE_dadt_counter_inc();
extern double RxODE_podo();
extern double RxODE_tlast();
extern void update_par_ptr(double t);
extern void RxODE_ode_free();

// Remove these functions later...
extern void RxODE_assign_fn_pointers(void (*fun_dydt)(unsigned int, double, double *, double *),
                                     void (*fun_calc_lhs)(double, double *, double *),
                                     void (*fun_calc_jac)(unsigned int, double, double *, double *, unsigned int),
                                     void (*fun_update_inis)(SEXP _ini_sexp),
                                     int fun_jt,
                                     int fun_mf,
                                     int fun_debug);


void R_init_RxODE(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"RxODE_ode_solver", (DL_FUNC) &RxODE_ode_solver, 23},
    {"trans", (DL_FUNC) &trans, 8},
    {"_RxODE_rxInv", (DL_FUNC) &_RxODE_rxInv, 1},
    {"_RxODE_RxODE_finalize_focei_omega",(DL_FUNC) &_RxODE_RxODE_finalize_focei_omega, 1},
    {"_RxODE_RxODE_finalize_log_det_OMGAinv_5",(DL_FUNC) &_RxODE_RxODE_finalize_log_det_OMGAinv_5, 1},
    {"_RxODE_rxCoutEcho", (DL_FUNC) &_RxODE_rxCoutEcho, 1},
    {"_RxODE_removableDrive", (DL_FUNC) &_RxODE_removableDrive, 1},
    {"_rxCholInv", (DL_FUNC) &_rxCholInv, 3},
    {"_RxODE_rxToOmega", (DL_FUNC) &_RxODE_rxToOmega, 1},
    {"_RxODE_rxSymInvCholEnvCalculate", (DL_FUNC) &_RxODE_rxSymInvCholEnvCalculate, 3},
    {"_RxODE_rxInvWishartVar", (DL_FUNC) &_RxODE_rxInvWishartVar, 2},
    {"_RxODE_rxSymInvChol", (DL_FUNC) &_RxODE_rxSymInvChol, 4},
    {"_RxODE_rxDataSetup", (DL_FUNC) &_RxODE_rxDataSetup, 8},
    {"_RxODE_rxIs", (DL_FUNC) &_RxODE_rxIs, 2},
    {"_RxODE_rxModelVars", (DL_FUNC) &_RxODE_rxModelVars, 1},
    {"_RxODE_rxState", (DL_FUNC) &_RxODE_rxState, 2},
    {"_RxODE_rxParams", (DL_FUNC) &_RxODE_rxParams, 1},
    {"_RxODE_rxDfdy", (DL_FUNC) &_RxODE_rxDfdy, 1},
    {"_RxODE_rxLhs", (DL_FUNC) &_RxODE_rxLhs, 1},
    {"_RxODE_rxInits", (DL_FUNC) &_RxODE_rxInits, 6},
    {"_RxODE_rxUpdateResiduals", (DL_FUNC) &_RxODE_rxUpdateResiduals, 1},
    {"_RxODE_rxSetupIni", (DL_FUNC) &_RxODE_rxSetupIni, 2},
    {"_RxODE_rxDataParSetup", (DL_FUNC) &_RxODE_rxDataParSetup, 11},
    {"_RxODE_rxSolvingOptions",(DL_FUNC) &_RxODE_rxSolvingOptions, 12},
    {"_RxODE_rxSolvingData", (DL_FUNC) &_RxODE_rxSolvingData, 14},
    {"_RxODE_rxData", (DL_FUNC) &_RxODE_rxData, 23},
    {NULL, NULL, 0}
  };

  // C callables needed in FOCEi
  R_RegisterCCallable("RxODE","nEq",                 (DL_FUNC) nEq);
  R_RegisterCCallable("RxODE","nLhs",                (DL_FUNC) nLhs);
  R_RegisterCCallable("RxODE","rxLhs",               (DL_FUNC) rxLhs);
  R_RegisterCCallable("RxODE","nAllTimes",           (DL_FUNC) nAllTimes);
  R_RegisterCCallable("RxODE","rxEvid",              (DL_FUNC) rxEvid);
  R_RegisterCCallable("RxODE","rxCalcLhs",           (DL_FUNC) rxCalcLhs);
  R_RegisterCCallable("RxODE","nObs",                (DL_FUNC) nObs);

  R_RegisterCCallable("RxODE","RxODE_ode_solve_env", (DL_FUNC) RxODE_ode_solve_env);
  R_RegisterCCallable("RxODE","RxODE_ode_free",      (DL_FUNC) RxODE_ode_free);
  R_RegisterCCallable("RxODE","RxODE_safe_zero",     (DL_FUNC) RxODE_safe_zero);
  R_RegisterCCallable("RxODE","RxODE_safe_log",      (DL_FUNC) RxODE_safe_log);
  R_RegisterCCallable("RxODE","RxODE_sign_exp",      (DL_FUNC) RxODE_sign_exp);
  R_RegisterCCallable("RxODE","RxODE_abs_log",       (DL_FUNC) RxODE_abs_log);
  
  //Functions
  R_RegisterCCallable("RxODE","RxODE_ode_solver",       (DL_FUNC) RxODE_ode_solver);
  R_RegisterCCallable("RxODE","RxODE_assign_fn_pointers", (DL_FUNC) RxODE_assign_fn_pointers);
  R_RegisterCCallable("RxODE","RxODE_get_fn_pointers", (DL_FUNC) RxODE_get_fn_pointers);
  R_RegisterCCallable("RxODE","RxODE_ode_solver_old_c", (DL_FUNC) RxODE_ode_solver_old_c);
  
  //Infusion
  R_RegisterCCallable("RxODE","RxODE_InfusionRate",     (DL_FUNC) RxODE_InfusionRate);
  // Parameters
  R_RegisterCCallable("RxODE","RxODE_par_ptr",          (DL_FUNC) RxODE_par_ptr);
  R_RegisterCCallable("RxODE","RxODE_update_par_ptr",   (DL_FUNC) update_par_ptr);
  // Counters
  R_RegisterCCallable("RxODE","RxODE_dadt_counter_val", (DL_FUNC) RxODE_dadt_counter_val);
  R_RegisterCCallable("RxODE","RxODE_jac_counter_val",  (DL_FUNC) RxODE_jac_counter_val);
  R_RegisterCCallable("RxODE","RxODE_dadt_counter_inc", (DL_FUNC) RxODE_dadt_counter_inc);
  R_RegisterCCallable("RxODE","RxODE_jac_counter_inc",  (DL_FUNC) RxODE_jac_counter_inc);
  // podo or tlast
  R_RegisterCCallable("RxODE","RxODE_podo",             (DL_FUNC) RxODE_podo);
  R_RegisterCCallable("RxODE","RxODE_tlast",            (DL_FUNC) RxODE_tlast);
  // tranit compartment models
  R_RegisterCCallable("RxODE","RxODE_transit4",         (DL_FUNC) RxODE_transit4);
  R_RegisterCCallable("RxODE","RxODE_transit3",         (DL_FUNC) RxODE_transit3);
  R_RegisterCCallable("RxODE","RxODE_factorial",        (DL_FUNC) RxODE_factorial);
  R_RegisterCCallable("RxODE","RxODE_safe_log",         (DL_FUNC) RxODE_safe_log);
  R_RegisterCCallable("RxODE","RxODE_safe_zero",        (DL_FUNC) RxODE_safe_zero);
  R_RegisterCCallable("RxODE","RxODE_as_zero",          (DL_FUNC) RxODE_as_zero);
  R_RegisterCCallable("RxODE","RxODE_sign_exp",         (DL_FUNC) RxODE_sign_exp);
  R_RegisterCCallable("RxODE","RxODE_abs_log",          (DL_FUNC) RxODE_abs_log);
  R_RegisterCCallable("RxODE","RxODE_abs_log1p",        (DL_FUNC) RxODE_abs_log1p);
  R_RegisterCCallable("RxODE","RxODE_solveLinB",        (DL_FUNC) RxODE_solveLinB);

  R_RegisterCCallable("RxODE","RxODE_sum",              (DL_FUNC) RxODE_sum);
  R_RegisterCCallable("RxODE","RxODE_prod",             (DL_FUNC) RxODE_prod);

  R_RegisterCCallable("RxODE","RxODE_pow",              (DL_FUNC) RxODE_pow);
  R_RegisterCCallable("RxODE","RxODE_pow_di",           (DL_FUNC) RxODE_pow_di);


  static const R_CMethodDef cMethods[] = {
    {"RxODE_InfusionRate",      (DL_FUNC) &RxODE_InfusionRate, 1, RxODE_one_int_t},
    {"RxODE_par_ptr",           (DL_FUNC) &RxODE_par_ptr, 1, RxODE_one_int_t},
    {"RxODE_jac_counter_val",   (DL_FUNC) &RxODE_jac_counter_val, 0},
    {"RxODE_dadt_counter_val",  (DL_FUNC) &RxODE_dadt_counter_val, 0},
    {"RxODE_jac_counter_inc",   (DL_FUNC) &RxODE_jac_counter_inc, 0},
    {"RxODE_dadt_counter_inc",  (DL_FUNC) &RxODE_dadt_counter_inc, 0},
    {"RxODE_podo",              (DL_FUNC) &RxODE_podo, 0},
    {"RxODE_tlast",             (DL_FUNC) &RxODE_tlast, 0},
    {"RxODE_transit4",          (DL_FUNC) &RxODE_transit4, 4, RxODE_transit4_t},
    {"RxODE_transit3",          (DL_FUNC) &RxODE_transit3, 4, RxODE_transit3_t},
    {"RxODE_factorial",         (DL_FUNC) &RxODE_factorial, 1, RxODE_one_dbl_t},
    {"RxODE_safe_log",          (DL_FUNC) &RxODE_safe_log, 1, RxODE_one_dbl_t},
    {"RxODE_safe_zero",         (DL_FUNC) &RxODE_safe_zero, 1, RxODE_one_dbl_t},
    {"RxODE_as_zero",           (DL_FUNC) &RxODE_as_zero, 1, RxODE_one_dbl_t},
    {"RxODE_sign_exp",          (DL_FUNC) &RxODE_sign_exp, 2, RxODE_sign_exp_t},
    {"RxODE_abs_log",           (DL_FUNC) &RxODE_abs_log, 1, RxODE_one_dbl_t},
    {"RxODE_abs_log1p",         (DL_FUNC) &RxODE_abs_log1p, 1, RxODE_one_dbl_t},
    {"RxODE_solveLinB",         (DL_FUNC) &RxODE_solveLinB, 12, RxODE_solveLinB_t},
    {"RxODE_sum",               (DL_FUNC) &RxODE_sum, 2, RxODE_Sum_t},
    {"RxODE_prod",              (DL_FUNC) &RxODE_prod, 2, RxODE_Sum_t},
    {NULL, NULL, 0, NULL}
  };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

