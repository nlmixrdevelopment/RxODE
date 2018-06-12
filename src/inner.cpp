#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
extern "C"{
#include "solve.h"
  void ind_solve(rx_solve *rx, unsigned int cid,
                 t_dydt_liblsoda dydt_lls,
                 t_dydt_lsoda_dum dydt_lsoda, t_jdum_lsoda jdum,
                 t_dydt c_dydt, t_update_inis u_inis);
}

// These are focei inner options
typedef struct {
  // Integer of ETAs
  int *etaTrans;
  unsigned int neta;
  unsigned int ntheta;
  double *theta;
  int *thetaTrans;
  double scaleTo;
  double epsilon;
  unsigned int max_iterations;
  unsigned int nzm;
} focei_options;

focei_options op_focei;


typedef struct {
  // F and varaibility
  mat f; // can change
  mat r; // can change
  //
  //mat dv; // doesn't change on update.
  mat err;
  //
  unsigned int neta;// = as<unsigned int>(e["neta"]);
  // Derivatives
  mat fpm;// = mat(nObs(), neta); or d(pred)/d(eta#)
  mat rp;// = mat(nObs(),neta);
    
  //
  mat B;// = mat(nObs(),1);
  mat c;// was List(neta); but should be (nObs(), neta)
  mat a;// was List(neta); but should be (nObs(), neta)
  
  // Likilihood gradient
  mat lp;// = mat(neta,1);
    
  double *eta;
  double *zm;
  unsigned int uzm;
} focei_ind;

focei_ind *inds_focei = NULL;
int max_inds_focei = 0;

rx_solve *getRxSolve_();

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

uvec lowerTri(mat H, bool diag = false){
  unsigned int d = H.n_rows;
  mat o(d, d, fill::ones);
  if (!diag){
    return find(trimatl(o,-1));
  } else {
    return find(trimatl(o));
  }
}

void updateZm(focei_ind *indF){
  if (!indF->uzm){
    // Udate the curvature to Hessian to restart n1qn1
    unsigned int n = op_focei.neta;
    mat L = eye(n,n);
    mat D = mat(n,n,fill::zeros);
    mat H = mat(n,n);
    unsigned int l_n = n*(n+1)/2;
    vec zmV(l_n);
    std::copy(&indF->zm[0], &indF->zm[0]+l_n, zmV.begin());
    H.elem(lowerTri(H,true)) = zmV;
    L.elem(lowerTri(H,false)) = H.elem(lowerTri(H,0));
    D.diag() = H.diag();
    H = L*D*L.t();
    // Hessian -> c.hess
    std::fill(&indF->zm[0], &indF->zm[0]+op_focei.nzm,0.0);
    vec hessV = H.elem(lowerTri(H,1));
    std::copy(hessV.begin(),hessV.end(),&indF->zm[0]);
    indF->uzm = 1;
  }
}
