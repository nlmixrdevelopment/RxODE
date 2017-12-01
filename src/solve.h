typedef void (*t_dydt)(unsigned int neq, double t, double *A, double *DADT);
typedef void (*t_calc_jac)(unsigned int neq, double t, double *A, double *JAC, unsigned int __NROWPD__);
typedef void (*t_calc_lhs)(double t, double *A, double *lhs);
typedef void (*t_update_inis)(double *);
typedef void (*t_dydt_lsoda_dum)(int *neq, double *t, double *A, double *DADT);
typedef void (*t_jdum_lsoda)(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd);

typedef struct {
  // These options should not change based on an individual solve
  double ATOL;          //absolute error
  double RTOL;          //relative error
  double H0;
  double HMIN;
  int global_jt;
  int global_mf;
  int global_debug;
  int mxstep;
  int MXORDN;
  int MXORDS;
  //
  int do_transit_abs;
  int nlhs;
  int neq;
  int stiff;
  t_dydt dydt;
  t_calc_jac calc_jac;
  t_calc_lhs calc_lhs;
  t_update_inis update_inis;
  t_dydt_lsoda_dum dydt_lsoda_dum;
  t_jdum_lsoda jdum_lsoda;
  // approx fun options
  double f1;
  double f2;
  int kind;
  double ylow;
  double yhigh;
  int is_locf;
  int cores;
} rx_solving_options;


typedef struct {
  long slvr_counter;
  long dadt_counter;
  long jac_counter;
  double *InfusionRate;
  int *BadDose;
  int nBadDose;
  double HMAX; // Determined by diff
  double tlast;
  double podo;
  double *par_ptr;
  double *inits;
  double *dose;
  double *solve;
  double *lhs;
  int    *par_cov;
  int  *evid;
  int   do_par_cov;
  int *rc;
  double *cov_ptr;
  int ncov;
  int n_all_times;
  int ixds;
  int ndoses;
  double *all_times;
  int *idose;
  int idosen;
  int extraCmt;
  FILE *fp;
  int id;
  int sim;
} rx_solving_options_ind;

typedef struct {
  rx_solving_options_ind *subjects;
  int nsub;
  int nsim;
  SEXP op;
} rx_solve;
