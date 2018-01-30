typedef void (*t_dydt)(int *neq, double t, double *A, double *DADT);
typedef void (*t_calc_jac)(int *neq, double t, double *A, double *JAC, unsigned int __NROWPD__);
typedef void (*t_calc_lhs)(int cSub, double t, double *A, double *lhs);
typedef void (*t_update_inis)(int cSub, double *);
typedef void (*t_dydt_lsoda_dum)(int *neq, double *t, double *A, double *DADT);
typedef void (*t_jdum_lsoda)(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd);
typedef int (*t_dydt_liblsoda)(double t, double *y, double *ydot, void *data);

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
  SEXP stateNames;
  SEXP lhsNames;
  SEXP paramNames;
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
  int *idose;
  int idosen;
  int id;
  int sim;
  double ylow;
  double yhigh;
} rx_solving_options_ind;

typedef struct {
  rx_solving_options_ind *subjects;
  int nsub;
  int nsim;
  int nobs;
  int add_cov;
  int matrix;
  int *stateIgnore;
  SEXP op;
} rx_solve;

typedef void (*t_set_solve)(rx_solve *);
typedef rx_solve *(*t_get_solve)();


rx_solve *getRxSolve_();
void rxSolveDataFree(SEXP ptr);
int rxUpdateResiduals_(SEXP md);
rx_solve *getRxSolve(SEXP ptr);

SEXP getSolvingOptionsPtr(double ATOL,          //absolute error
                          double RTOL,          //relative error
                          double H0,
                          double HMIN,
                          int mxstep,
                          int MXORDN,
                          int MXORDS,
                          // Approx options
                          int do_transit_abs,
                          int nlhs,
                          int neq,
                          int stiff,
                          double f1,
                          double f2,
                          int kind,
                          int is_locf,
                          int cores,
                          int ncov,
                          int *par_cov,
                          int do_par_cov,
                          double *inits,
			  double *scale,
                          SEXP stateNames,
                          SEXP lhsNames,
                          SEXP paramNames,
			  double hmax2,
                          double *atol2,
                          double *rtol2);
void getSolvingOptionsIndPtr(double *InfusionRate,
                             int *BadDose,
                             double HMAX, // Determined by diff
                             double *par_ptr,
                             double *dose,
                             int *idose,
                             double *solve,
                             double *lhs,
                             int *evid,
                             int *rc,
                             double *cov_ptr,
                             int n_all_times,
                             double *all_times,
                             int id,
                             int sim,
                             rx_solving_options_ind *o);
SEXP rxSolveData(rx_solving_options_ind *subjects,
                 int nsub,
                 int nsim,
                 int *stateIgnore,
                 int nobs,
                 int add_cov,
                 int matrix,
                 SEXP op);
void par_solve(rx_solve *rx);

rx_solving_options *getRxOp(rx_solve *rx);
SEXP RxODE_df(SEXP sd, int doDose, int ini_updateR);
SEXP RxODE_par_df(SEXP sd);

rx_solving_options_ind *rxOptionsIniEnsure(int mx);

void rxUpdateFuns(SEXP trans);
