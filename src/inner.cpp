#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
extern "C"{
#include "solve.h"
}

// These are focei inner options
typedef struct {
  int badSolve;
  double ATOL;          //absolute error
  double RTOL;          //relative error
  double H0;
  double HMIN;
  int mxstep;
  int MXORDN;
  int MXORDS;
  //
  int do_transit_abs;
  int nlhs;
  int neq;
  int stiff;
  int ncov;
  int *par_cov;
  double *inits;
  int do_par_cov;
  // approx fun options
  double f1;
  double f2;
  int kind;
  int is_locf;
  int cores;
  int extraCmt;
  double hmax2; // Determined by diff
  double *rtol2;
  double *atol2;
  int abort;
  // Integer of ETAs
  int *etas;
  unsigned int netas;
} focei_options;

focei_options op_focei;


typedef struct {
  int *slvr_counter;
  int *dadt_counter;
  int *jac_counter;
  double *InfusionRate;
  int *BadDose;
  int nBadDose;
  double HMAX; // Determined by diff
  double tlast;
  double podo;
  double *par_ptr; // Includes ETAs
  double *dose;
  double *solve;
  double *lhs;
  double *nzm; // Saved Hessian information
  int  *evid;
  int *rc;
  double *cov_ptr;
  int n_all_times;
  int ixds;
  int ndoses;
  double *all_times;
  int *idose;
  int idosen;
  int id;
  int sim;
  double ylow;
  double yhigh;
} focei_ind;

focei_ind *inds_focei = NULL;
int max_inds_focei = 0;

focei_ind *rxFoceiEnsure(int mx){
  if (mx >= max_inds_focei){
    Free(inds_focei);
    inds_focei =Calloc(mx+1024, focei_ind);
    max_inds_focei = mx+1024;
  }
  return inds_focei;
}


t_dydt inner_dydt = NULL;

t_calc_jac inner_calc_jac = NULL;

t_calc_lhs inner_calc_lhs = NULL;

t_update_inis inner_update_inis = NULL;

t_dydt_lsoda_dum inner_dydt_lsoda_dum = NULL;

t_dydt_liblsoda inner_dydt_liblsoda = NULL;

t_jdum_lsoda inner_jdum_lsoda = NULL;

t_set_solve inner_set_solve = NULL;

t_get_solve inner_get_solve = NULL;

int inner_global_jt = 2;
int inner_global_mf = 22;  
int inner_global_debug = 0;

void rxUpdateInnerFuns(SEXP trans){
  const char *lib, *s_dydt, *s_calc_jac, *s_calc_lhs, *s_inis, *s_dydt_lsoda_dum, *s_dydt_jdum_lsoda, 
    *s_ode_solver_solvedata, *s_ode_solver_get_solvedata, *s_dydt_liblsoda;
  lib = CHAR(STRING_ELT(trans, 0));
  s_dydt = CHAR(STRING_ELT(trans, 3));
  s_calc_jac = CHAR(STRING_ELT(trans, 4));
  s_calc_lhs = CHAR(STRING_ELT(trans, 5));
  s_inis = CHAR(STRING_ELT(trans, 8));
  s_dydt_lsoda_dum = CHAR(STRING_ELT(trans, 9));
  s_dydt_jdum_lsoda = CHAR(STRING_ELT(trans, 10));
  s_ode_solver_solvedata = CHAR(STRING_ELT(trans, 11));
  s_ode_solver_get_solvedata = CHAR(STRING_ELT(trans, 12));
  s_dydt_liblsoda = CHAR(STRING_ELT(trans, 13));
  inner_global_jt = 2;
  inner_global_mf = 22;  
  inner_global_debug = 0;
  if (strcmp(CHAR(STRING_ELT(trans, 1)),"fulluser") == 0){
    inner_global_jt = 1;
    inner_global_mf = 21;
  } else {
    inner_global_jt = 2;
    inner_global_mf = 22;
  }
  inner_calc_lhs =(t_calc_lhs) R_GetCCallable(lib, s_calc_lhs);
  inner_dydt =(t_dydt) R_GetCCallable(lib, s_dydt);
  inner_calc_jac =(t_calc_jac) R_GetCCallable(lib, s_calc_jac);
  inner_update_inis =(t_update_inis) R_GetCCallable(lib, s_inis);
  inner_dydt_lsoda_dum =(t_dydt_lsoda_dum) R_GetCCallable(lib, s_dydt_lsoda_dum);
  inner_jdum_lsoda =(t_jdum_lsoda) R_GetCCallable(lib, s_dydt_jdum_lsoda);
  inner_set_solve = (t_set_solve)R_GetCCallable(lib, s_ode_solver_solvedata);
  inner_get_solve = (t_get_solve)R_GetCCallable(lib, s_ode_solver_get_solvedata);
  inner_dydt_liblsoda = (t_dydt_liblsoda)R_GetCCallable(lib, s_dydt_liblsoda);
}

void rxClearInnerFuns(){
  inner_calc_lhs              = NULL;
  inner_dydt                  = NULL;
  inner_calc_jac              = NULL;
  inner_update_inis           = NULL;
  inner_dydt_lsoda_dum        = NULL;
  inner_jdum_lsoda            = NULL;
  inner_set_solve             = NULL;
  inner_get_solve             = NULL;
  inner_dydt_liblsoda         = NULL;
}

void innerEtaLik(unsigned int cid){
  focei_ind *ind = &(inds_focei[cid]);
  unsigned int nObs = ind->n_all_times - ind->ndoses;

  mat fpm = mat(nObs, op_focei.netas); // d(pred)/d(eta#)
    
  mat rp = mat(nObs, op_focei.netas);
    
  mat f(nObs, 1);
  mat err = mat(nObs,1);
  mat r = mat(nObs,1);

  mat B = mat(nObs,1);
  // Cant use list...
  // List c(neta);
  // List a(neta);
  mat c = mat(nObs, op_focei.netas);
  mat a = mat(nObs, op_focei.netas);
  double llik =0.0;
  mat lp = mat(op_focei.netas, 1);

}

// typedef void (*S2_fp) (int *, int *, double *, double *, double *, int *, float *, double *);
// typedef void (*fun_n1qn1) (S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
//                         int mode[], int niter[], int nsim[], int imp[], int lp[], double zm[], int izs[], float rzs[], double dzs[]);

// void n1qn1(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
// 	   int mode[], int niter[], int nsim[], int imp[], int lp[], double zm[], int izs[], float rzs[], double dzs[]){
//   static fun_n1qn1 fun=NULL;
//   if (fun == NULL) fun = (fun_n1qn1) R_GetCCallable("n1qn1","n1qn1_");
//   fun(simul, n, x,f, g, var, eps, mode,  niter, nsim, imp, zm, izs, rzs, dzs);
// }


