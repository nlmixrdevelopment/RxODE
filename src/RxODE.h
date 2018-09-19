typedef void (*t_dydt)(int *neq, double t, double *A, double *DADT);
typedef void (*t_calc_jac)(int *neq, double t, double *A, double *JAC, unsigned int __NROWPD__);
typedef void (*t_calc_lhs)(int cSub, double t, double *A, double *lhs);
typedef void (*t_update_inis)(int cSub, double *);
typedef void (*t_dydt_lsoda_dum)(int *neq, double *t, double *A, double *DADT);
typedef void (*t_jdum_lsoda)(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd);
typedef int (*t_dydt_liblsoda)(double t, double *y, double *ydot, void *data);
typedef void (*t_ode_current)();

typedef struct {
  // These options should not change based on an individual solve
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
  char modNamePtr[1000];
  int *par_cov;
  double *inits;
  double *scale;
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
  int nDisplayProgress;
  int ncoresRV;
  int isChol;
  int *svar;
  int abort;
} rx_solving_options;


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
  double *par_ptr;
  double *dose;
  double *solve;
  double *lhs;
  int  *evid;
  int *rc;
  double *cov_ptr;
  int n_all_times;
  int ixds;
  int ndoses;
  double *all_times;
  double *dv;
  int *idose;
  int idosen;
  int id;
  int sim;
  int idx;
  double ylow;
  double yhigh;
  double lambda;
  double yj;
} rx_solving_options_ind;

typedef struct {
  rx_solving_options_ind *subjects;
  rx_solving_options *op;
  int nsub;
  int nsim;
  int nall;
  int nobs;
  int nr;
  int add_cov;
  int matrix;
  double stateTrim;
  int *stateIgnore;
} rx_solve;

typedef void (*t_set_solve)(rx_solve *);
typedef rx_solve *(*t_get_solve)();


rx_solve *getRxSolve_();
rx_solve *getRxSolve(SEXP ptr);

void par_solve(rx_solve *rx);

rx_solving_options *getRxOp(rx_solve *rx);

SEXP RxODE_df(int doDose, int doTBS);
SEXP RxODE_par_df();

rx_solving_options_ind *rxOptionsIniEnsure(int mx);

void rxUpdateFuns(SEXP trans);
