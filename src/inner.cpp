#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

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
} focei_options;


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
  double *nzm;
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


typedef void (*S2_fp) (int *, int *, double *, double *, double *, int *, float *, double *);
typedef void (*fun_n1qn1) (S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
                        int mode[], int niter[], int nsim[], int imp[], int lp[], double zm[], int izs[], float rzs[], double dzs[]);

void n1qn1(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
	   int mode[], int niter[], int nsim[], int imp[], int lp[], double zm[], int izs[], float rzs[], double dzs[]){
  static fun_n1qn1 fun=NULL;
  if (fun == NULL) fun = (fun_n1qn1) R_GetCCallable("n1qn1","n1qn1_");
  fun(simul, n, x,f, g, var, eps, mode,  niter, nsim, imp, zm, izs, rzs, dzs);
}


