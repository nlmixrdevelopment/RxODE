#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#define JAC_Rprintf Rprintf
#define JAC0_Rprintf if ( (&_solveData->subjects[_cSub])->jac_counter == 0) Rprintf
#define ODE_Rprintf Rprintf
#define ODE0_Rprintf if ( (&_solveData->subjects[_cSub])->dadt_counter == 0) Rprintf
#define LHS_Rprintf Rprintf
#define _safe_log(a) (((a) <= 0) ? log(DOUBLE_EPS) : log(a))
#define safe_zero(a) ((a) == 0 ? DOUBLE_EPS : (a))
#define _as_zero(a) (fabs(a) < sqrt(DOUBLE_EPS) ? 0.0 : a)
#define factorial(a) exp(lgamma1p(a))
#define sign_exp(sgn, x)(((sgn) > 0.0) ? exp(x) : (((sgn) < 0.0) ? -exp(x) : 0.0))
#define R_pow(a, b) (((a) == 0 && (b) <= 0) ? R_pow(DOUBLE_EPS, b) : R_pow(a, b))
#define R_pow_di(a, b) (((a) == 0 && (b) <= 0) ? R_pow_di(DOUBLE_EPS, b) : R_pow_di(a, b))
#define Rx_pow(a, b) (((a) == 0 && (b) <= 0) ? R_pow(DOUBLE_EPS, b) : R_pow(a, b))
#define Rx_pow_di(a, b) (((a) == 0 && (b) <= 0) ? R_pow_di(DOUBLE_EPS, b) : R_pow_di(a, b))
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
  double ylow;
  double yhigh;
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
  int *stateIgnore;
} rx_solve;

typedef void (*t_set_solve)(rx_solve *);
typedef rx_solve *(*t_get_solve)();


rx_solve *getRxSolve_();
rx_solve *getRxSolve(SEXP ptr);

void par_solve(rx_solve *rx);

rx_solving_options *getRxOp(rx_solve *rx);

SEXP RxODE_df(int doDose);
SEXP RxODE_par_df();

rx_solving_options_ind *rxOptionsIniEnsure(int mx);

void rxUpdateFuns(SEXP trans);
#define abs_log1p(x) (((x) + 1.0 > 0.0) ? log1p(x) : (((x) + 1.0 > 0.0) ? log1p(-x) : 0.0))
#define abs_log(x) ((fabs(x) <= sqrt(DOUBLE_EPS)) ? log(sqrt(DOUBLE_EPS)) : (((x) > 0.0) ? log(x) ? (((x) == 0) ? 0.0 : log(-x))))
#define _IR (_solveData->subjects[_cSub].InfusionRate)
#define _PP (_solveData->subjects[_cSub].par_ptr)
#define _SR (INTEGER(stateRmS))

// Types for par pointers.r
typedef double (*RxODE_fn) (double x);
typedef double (*RxODE_fn2) (double x, double y);
typedef double (*RxODE_fn2i) (double x, int i);
typedef int (*RxODE_fn0i) ();
typedef double (*RxODE_vec) (int val, rx_solve *rx, unsigned int id);
typedef double (*RxODE_val) (rx_solve *rx, unsigned int id);
typedef void (*RxODE_assign_ptr)(SEXP);
typedef void (*RxODE_ode_solver_old_c)(int *neq,double *theta,double *time,int *evid,int *ntime,double *inits,double *dose,double *ret,double *atol,double *rtol,int *stiff,int *transit_abs,int *nlhs,double *lhs,int *rc);

RxODE_assign_ptr _assign_ptr = NULL;

typedef void (*_rxRmModelLibType)(const char *inp);
_rxRmModelLibType _rxRmModelLib = NULL;

typedef SEXP (*_rxGetModelLibType)(const char *s);
_rxGetModelLibType _rxGetModelLib = NULL;

RxODE_ode_solver_old_c _old_c = NULL;

RxODE_fn0i _ptrid=NULL;

typedef  SEXP (*_rx_asgn) (SEXP objectSEXP);
_rx_asgn _RxODE_rxAssignPtr =NULL;

typedef int(*_rxIsCurrentC_type)(SEXP);
_rxIsCurrentC_type _rxIsCurrentC=NULL;

typedef double(*_rxSumType)(double *, int, double *, int, int);
_rxSumType _sumPS = NULL;

double _sum(double *input, double *pld, int m, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  return _sumPS(input, n, pld, m, type);
}

typedef double(*_rxProdType)(double*, double*, int, int);
_rxProdType _prodPS = NULL;

double _prod(double *input, double *p, int type, int n, ...){
  va_list valist;
  va_start(valist, n);
  for (unsigned int i = 0; i < n; i++){
    input[i] = va_arg(valist, double);
  }
  va_end(valist);
  return _prodPS(input, p, n, type);
}

double _sign(unsigned int n, ...){
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

double _max(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double mx = NA_REAL;
  double tmp = 0;
  if (n >= 1){
    mx = va_arg(valist, double);
    for (unsigned int i = 1; i < n; i++){
      tmp = va_arg(valist, double);
      if (tmp>mx) mx=tmp;
    }
    va_end(valist);
  }
  return mx;
}

double _min(unsigned int n, ...){
  va_list valist;
  va_start(valist, n);
  double mn = NA_REAL;
  double tmp = 0;
  if (n >= 1){
    mn = va_arg(valist, double);
    for (unsigned int i = 1; i < n; i++){
      tmp = va_arg(valist, double);
      if (tmp<mn) mn=tmp;
    }
    va_end(valist);
  }
  return mn;
}

rx_solve *_solveData = NULL;

/* Authors: Robert Gentleman and Ross Ihaka and The R Core Team */
/* Taken directly from https://github.com/wch/r-source/blob/922777f2a0363fd6fe07e926971547dd8315fc24/src/library/stats/src/approx.c*/
/* Changed as follows:
   - Different Name
   - Use RxODE structure
   - Make inline
*/
static inline double rx_approxP(double v, double *x, double *y, int n,
                         rx_solving_options *Meth, rx_solving_options_ind *id){
  /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
  int i, j, ij;

  if(!n) return R_NaN;

  i = 0;
  j = n - 1;

  /* handle out-of-domain points */
  if(v < x[i]) return id->ylow;
  if(v > x[j]) return id->yhigh;

  /* find the correct interval by bisection */
  while(i < j - 1) { /* x[i] <= v <= x[j] */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    if(v < x[ij]) j = ij; else i = ij;
    /* still i < j */
  }
  /* provably have i == j-1 */

  /* interpolation */

  if(v == x[j]) return y[j];
  if(v == x[i]) return y[i];
  /* impossible: if(x[j] == x[i]) return y[i]; */

  if(Meth->kind == 1) /* linear */
    return y[i] + (y[j] - y[i]) * ((v - x[i])/(x[j] - x[i]));
  else /* 2 : constant */
    return (Meth->f1 != 0.0 ? y[i] * Meth->f1 : 0.0)
      + (Meth->f2 != 0.0 ? y[j] * Meth->f2 : 0.0);
}/* approx1() */

/* End approx from R */


static inline void _update_par_ptr(double t, unsigned int id){
  rx_solving_options_ind *ind;
  ind = (&_solveData->subjects[id]);
  rx_solving_options *op = _solveData->op;
  if (op->neq > 0){
    // Update all covariate parameters
    int k;
    int ncov = op->ncov;
    if (op->do_par_cov){
      for (k = ncov; k--;){
        if (op->par_cov[k]){
	  double *par_ptr = ind->par_ptr;
          double *all_times = ind->all_times;
          double *cov_ptr = ind->cov_ptr;
          // Use the same methodology as approxfun.
          // There is some rumor the C function may go away...
          ind->ylow = cov_ptr[ind->n_all_times*k];
          ind->yhigh = cov_ptr[ind->n_all_times*k+ind->n_all_times-1];
          par_ptr[op->par_cov[k]-1] = rx_approxP(t, all_times, cov_ptr+ind->n_all_times*k, ind->n_all_times, op, ind);
        }
      }
    }
  }
}

static inline double _transit4P(double t, unsigned int id, double n, double mtt, double bio){
  double ktr = (n+1)/mtt;
  double lktr = log(n+1)-log(mtt);
  double tc = (t-(_solveData->subjects[id].tlast));
  return exp(log(bio*(_solveData->subjects[id].podo))+lktr+n*(lktr+log(tc))-ktr*(tc)-lgamma1p(n));
}

static inline double _transit3P(double t, unsigned int id, double n, double mtt){
  double ktr = (n+1)/mtt;
  double lktr = log(n+1)-log(mtt);
  double tc = (t-(_solveData->subjects[id].tlast));
  return exp(log(_solveData->subjects[id].podo)+lktr+n*(lktr+log(tc))-ktr*(tc)-lgamma1p(n));
}


// Linear compartment models/functions

static inline int _locateDoseIndex(const double obs_time,  rx_solving_options_ind *ind){
  // Uses bisection for slightly faster lookup of dose index.
  int i, j, ij;
  i = 0;
  j = ind->ndoses - 1;
  if (obs_time <= ind->all_times[ind->idose[i]]){
    while(i < ind->ndoses-2 && obs_time == ind->all_times[ind->idose[i+1]]){
      i++;
    }
    return i;
  }
  if (obs_time >= ind->all_times[ind->idose[j]]){
    return j;
  }
  while(i < j - 1) { /* x[i] <= obs_time <= x[j] */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    if(obs_time < ind->all_times[ind->idose[ij]])
      j = ij;
    else
      i = ij;
  }
  while(i < ind->ndoses-2 && obs_time == ind->all_times[ind->idose[i+1]]){
    i++;
  }
  return i;
}

static inline double solveLinB(rx_solve *rx, unsigned int id, double t, int linCmt, int diff1, int diff2, double d_A, double d_alpha, double d_B, double d_beta, double d_C, double d_gamma, double d_ka, double d_tlag){
  if (diff1 != 0 || diff2 != 0){
    error("Exact derivtives are no longer calculated.");
  }
  unsigned int ncmt = 1;
  double beta1=0, gamma1=0, alpha1=0;
  double alpha = d_alpha;
  double A = d_A;
  double beta = d_beta;
  double B = d_B;
  double gamma = d_gamma;
  double C = d_C;
  double ka = d_ka;
  double tlag = d_tlag;
  if (d_gamma > 0.){
    ncmt = 3;
    gamma1 = 1.0/gamma;
    beta1 = 1.0/beta;
    alpha1 = 1.0/alpha;
  } else if (d_beta > 0.){
    ncmt = 2;
    beta1 = 1.0/beta;
    alpha1 = 1.0/alpha;
  } else if (d_alpha > 0.){
    ncmt = 1;
    alpha1 = 1.0/alpha;
  } else {
    return 0.0;
    //error("You need to specify at least A(=%f) and alpha (=%f). (@t=%f, d1=%d, d2=%d)", d_A, d_alpha, t, diff1, diff2);
  }
  rx_solving_options *op = _solveData->op;
  double ATOL = op->ATOL;          //absolute error
  double RTOL = op->RTOL;          //relative error
  if (linCmt+1 > op->extraCmt){
    op->extraCmt = linCmt+1;
  }
  int oral, cmt;
  oral = (ka > 0) ? 1 : 0;
  double ret = 0,cur=0, tmp=0;
  unsigned int m = 0, l = 0, p = 0;
  int evid, evid100;
  double thisT = 0.0, tT = 0.0, res, t1, t2, tinf, dose = 0;
  double rate;
  rx_solving_options_ind *ind = &(_solveData->subjects[id]);
  if (ind->ndoses < 0){
    ind->ndoses=0;
    for (unsigned int i = 0; i < ind->n_all_times; i++){
      if (ind->evid[i]){
        ind->ndoses++;
        ind->idose[ind->ndoses-1] = i;
      }
    }
  }
  m = _locateDoseIndex(t, ind);
  int ndoses = ind->ndoses;
  for(l=m+1; l--;){// Optimized for loop as https://www.thegeekstuff.com/2015/01/c-cpp-code-optimization/
    cur=0;
    //superpostion
    evid = ind->evid[ind->idose[l]];
    dose = ind->dose[l];
    // Support 100+ compartments...
    evid100 = floor(evid/1e5);
    evid = evid- evid100*1e5;
    cmt = (evid%10000)/100 - 1 + 100*evid100;
    if (cmt != linCmt) continue;
    if (evid > 10000) {
      if (dose > 0){
        // During infusion
        tT = t - ind->all_times[ind->idose[l]] ;
        thisT = tT - tlag;
        p = l+1;
        while (p < ndoses && ind->dose[p] != -dose){
          p++;
        }
        if (ind->dose[p] != -dose){
          error("Could not find a error to the infusion.  Check the event table.");
        }
        tinf  = ind->all_times[ind->idose[p]] - ind->all_times[ind->idose[l]];
        rate  = dose;
        if (tT >= tinf) continue;
      } else {
        // After  infusion
        p = l-1;
        while (p > 0 && ind->dose[p] != -dose){
          p--;
        }
        if (ind->dose[p] != -dose){
          error("Could not find a start to the infusion.  Check the event table.");
        }
        tinf  = ind->all_times[ind->idose[l]] - ind->all_times[ind->idose[p]] - tlag;
        
        tT = t - ind->all_times[ind->idose[p]];
        thisT = tT -tlag;
        rate  = -dose;
      }
      t1 = ((thisT < tinf) ? thisT : tinf);        //during infusion
      t2 = ((thisT > tinf) ? thisT - tinf : 0.0);  // after infusion
      cur +=  rate*A*alpha1*(1.0-exp(-alpha*t1))*exp(-alpha*t2);
      if (ncmt >= 2){
        cur +=  rate*B*beta1*(1.0-exp(-beta*t1))*exp(-beta*t2);
        if (ncmt >= 3){
          cur +=  rate*C*gamma1*(1.0-exp(-gamma*t1))*exp(-gamma*t2);
        }
      }
    } else {
      tT = t - ind->all_times[ind->idose[l]];
      thisT = tT -tlag;
      if (thisT < 0) continue;
      res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
      cur +=  dose*A*(exp(-alpha*thisT)-res);
      if (ncmt >= 2){
        cur +=  dose*B*(exp(-beta*thisT)-res);
        if (ncmt >= 3){
          cur += dose*C*(exp(-gamma*thisT)-res);
        }
      }
    }
    // Since this starts with the most recent dose, and then goes
    // backward, you can use a tolerance calcuation to exit the loop
    // early.
    //
    // See  http://web.mit.edu/10.001/Web/Tips/Converge.htm
    //
    // | True value - Computed value | < RTOL*|True Value| + ATOL 
    // | (ret+cur) - ret| < RTOL*|ret+cur|+ATOL
    // | cur | < RTOL*|ret+cur|+ATOL
    //
    // For this calcuation all values should be > 0.  If they are less
    // than 0 then it is approximately zero.
    if (cur < 0) break;
    tmp = ret+cur;
    if (cur < RTOL*tmp+ATOL){ 
      ret=tmp;
      break;
    }
    ret = tmp;
  } //l
  return ret;
}

RxODE_fn0i _prodType = NULL;
RxODE_fn0i _sumType = NULL;

extern void m1_x64_ode_solver_solvedata (rx_solve *solve){
  _solveData = solve;
}

extern rx_solve *m1_x64_ode_solver_get_solvedata(){
  return _solveData;
}

SEXP m1_x64_model_vars();
extern void m1_x64_ode_solver(int *neq,
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
			   int *rc){
  // Backward compatible ode solver for 0.5* C interface
  //if (_ptrid() != 1529118994 ){ _assign_ptr(m1_x64_model_vars());}
  double *_theta = theta;
  _old_c(neq, _theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc);
}

static R_NativePrimitiveArgType m1_x64_ode_solverrx_t[] = {
  //*neq, *theta, *time,  *evid, *ntime, *inits,   *dose,   *ret,     *atol,  *rtol,   *stiff, *transit_abs, *nlhs, *lhs, *rc
  INTSXP,REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP
};


// prj-specific differential eqns
void m1_x64_dydt(int *_neq, double t, double *__zzStateVar__, double *__DDtStateVar__)
{
  int _cSub = _neq[1];
  double   C2,
  centr,
  V2,
  C3,
  peri,
  V3,
  depot,
  KA,
  CL,
  Q,
  eff,
  Kin,
  Kout,
  EC50;

  (void)t;
  (void)C2;
  (void)centr;
  (void)V2;
  (void)C3;
  (void)peri;
  (void)V3;
  (void)depot;
  (void)KA;
  (void)CL;
  (void)Q;
  (void)eff;
  (void)Kin;
  (void)Kout;
  (void)EC50;

  _update_par_ptr(t, _cSub);
  V2 = _PP[0];
  V3 = _PP[1];
  KA = _PP[2];
  CL = _PP[3];
  Q = _PP[4];
  Kin = _PP[5];
  Kout = _PP[6];
  EC50 = _PP[7];

  depot = __zzStateVar__[0];
  centr = __zzStateVar__[1];
  peri = __zzStateVar__[2];
  eff = __zzStateVar__[3];

  C2=centr/safe_zero(V2);
  C3=peri/safe_zero(V3);
  __DDtStateVar__[0]=_IR[0]-KA*depot;
  __DDtStateVar__[1]=_IR[1]+KA*depot-CL*C2-Q*C2+Q*C3;
  __DDtStateVar__[2]=_IR[2]+Q*C2-Q*C3;
  __DDtStateVar__[3]=_IR[3]+Kin-Kout*(1-C2/safe_zero((EC50+C2)))*eff;
  (&_solveData->subjects[_cSub])->dadt_counter[0]++;
}

// Jacobian derived vars
void m1_x64_calc_jac(int *_neq, double t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {
  int _cSub=_neq[1];
  (&_solveData->subjects[_cSub])->jac_counter[0]++;
}
// Functional based initial conditions.
void m1_x64_inis(int _cSub, double *__zzStateVar__){
  double t=0;
  double   C2,
  centr,
  V2,
  C3,
  peri,
  V3,
  depot,
  KA,
  CL,
  Q,
  eff,
  Kin,
  Kout,
  EC50;

  (void)t;
  (void)C2;
  (void)centr;
  (void)V2;
  (void)C3;
  (void)peri;
  (void)V3;
  (void)depot;
  (void)KA;
  (void)CL;
  (void)Q;
  (void)eff;
  (void)Kin;
  (void)Kout;
  (void)EC50;

  _update_par_ptr(0.0, _cSub);
  V2 = _PP[0];
  V3 = _PP[1];
  KA = _PP[2];
  CL = _PP[3];
  Q = _PP[4];
  Kin = _PP[5];
  Kout = _PP[6];
  EC50 = _PP[7];

  depot = __zzStateVar__[0];
  centr = __zzStateVar__[1];
  peri = __zzStateVar__[2];
  eff = __zzStateVar__[3];

  C2=centr/safe_zero(V2);
  C3=peri/safe_zero(V3);
  __zzStateVar__[0]=depot;
  __zzStateVar__[1]=centr;
  __zzStateVar__[2]=peri;
  __zzStateVar__[3]=eff;
}
// prj-specific derived vars
void m1_x64_calc_lhs(int _cSub, double t, double *__zzStateVar__, double *_lhs) {
  double   __DDtStateVar_0__,
  __DDtStateVar_1__,
  __DDtStateVar_2__,
  __DDtStateVar_3__,
  C2,
  centr,
  V2,
  C3,
  peri,
  V3,
  depot,
  KA,
  CL,
  Q,
  eff,
  Kin,
  Kout,
  EC50;

  (void)t;
  (void)__DDtStateVar_0__;
  (void)__DDtStateVar_1__;
  (void)__DDtStateVar_2__;
  (void)__DDtStateVar_3__;
  (void)C2;
  (void)centr;
  (void)V2;
  (void)C3;
  (void)peri;
  (void)V3;
  (void)depot;
  (void)KA;
  (void)CL;
  (void)Q;
  (void)eff;
  (void)Kin;
  (void)Kout;
  (void)EC50;

  _update_par_ptr(t, _cSub);
  V2 = _PP[0];
  V3 = _PP[1];
  KA = _PP[2];
  CL = _PP[3];
  Q = _PP[4];
  Kin = _PP[5];
  Kout = _PP[6];
  EC50 = _PP[7];

  depot = __zzStateVar__[0];
  centr = __zzStateVar__[1];
  peri = __zzStateVar__[2];
  eff = __zzStateVar__[3];

  C2=centr/safe_zero(V2);
  C3=peri/safe_zero(V3);
  __DDtStateVar_0__=_IR[0]-KA*depot;
  __DDtStateVar_1__=_IR[1]+KA*depot-CL*C2-Q*C2+Q*C3;
  __DDtStateVar_2__=_IR[2]+Q*C2-Q*C3;
  __DDtStateVar_3__=_IR[3]+Kin-Kout*(1-C2/safe_zero((EC50+C2)))*eff;

  _lhs[0]=C2;
  _lhs[1]=C3;
}
extern SEXP m1_x64_model_vars(){
  int pro=0;
  SEXP _mv = PROTECT(_rxGetModelLib("rx_5c2c6f8a65d272301b81504c87d75239_x64_model_vars"));pro++;
  if (!_rxIsCurrentC(_mv)){
    SEXP lst      = PROTECT(allocVector(VECSXP, 16));pro++;
    SEXP names    = PROTECT(allocVector(STRSXP, 16));pro++;
    SEXP params   = PROTECT(allocVector(STRSXP, 8));pro++;
    SEXP lhs      = PROTECT(allocVector(STRSXP, 2));pro++;
    SEXP state    = PROTECT(allocVector(STRSXP, 4));pro++;
    SEXP stateRmS = PROTECT(allocVector(INTSXP, 4));pro++;
    SEXP timeInt = PROTECT(allocVector(INTSXP, 1));pro++;
    INTEGER(timeInt)[0] = 1529118994;
    SEXP sens     = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP normState= PROTECT(allocVector(STRSXP, 4));pro++;
    SEXP fn_ini   = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP dfdy     = PROTECT(allocVector(STRSXP, 0));pro++;
    SEXP tran     = PROTECT(allocVector(STRSXP, 14));pro++;
    SEXP trann    = PROTECT(allocVector(STRSXP, 14));pro++;
    SEXP mmd5     = PROTECT(allocVector(STRSXP, 2));pro++;
    SEXP mmd5n    = PROTECT(allocVector(STRSXP, 2));pro++;
    SEXP model    = PROTECT(allocVector(STRSXP, 4));pro++;
    SEXP modeln   = PROTECT(allocVector(STRSXP, 4));pro++;
    SEXP solve    = PROTECT(allocVector(VECSXP, 4));pro++;
    SEXP solven   = PROTECT(allocVector(STRSXP, 4));pro++;
    SEXP initsr   = PROTECT(allocVector(REALSXP, 4));pro++;
    SEXP scaler   = PROTECT(allocVector(REALSXP, 4));pro++;
    SEXP infusionr= PROTECT(allocVector(REALSXP, 4));pro++;
    SEXP badDosei = PROTECT(allocVector(INTSXP, 4));pro++;
    SEXP version    = PROTECT(allocVector(STRSXP, 3));pro++;
    SEXP versionn   = PROTECT(allocVector(STRSXP, 3));pro++;
    SET_STRING_ELT(version,0,mkChar("0.7.2-2"));
    SET_STRING_ELT(version,1,mkChar("https://github.com/nlmixrdevelopment/RxODE"));
    SET_STRING_ELT(version,2,mkChar("fc4916a07950ad2bfbf61583bf729db2"));
    SET_STRING_ELT(versionn,0,mkChar("version"));
    SET_STRING_ELT(versionn,1,mkChar("repo"));
    SET_STRING_ELT(versionn,2,mkChar("md5"));
    SET_STRING_ELT(solven,0,mkChar("inits"));
    SET_VECTOR_ELT(solve,  0,initsr);
    SET_STRING_ELT(solven,1,mkChar("scale"));
    SET_VECTOR_ELT(solve,  1,scaler);
    SET_STRING_ELT(solven,2,mkChar("infusion"));
    SET_VECTOR_ELT(solve,  2,infusionr);
    SET_STRING_ELT(solven,3,mkChar("badDose"));
    SET_VECTOR_ELT(solve,  3,badDosei);
    setAttrib(solve, R_NamesSymbol, solven);
  SET_STRING_ELT(lhs,0,mkChar("C2"));
    SET_STRING_ELT(params,0,mkChar("V2"));
  SET_STRING_ELT(lhs,1,mkChar("C3"));
    SET_STRING_ELT(params,1,mkChar("V3"));
    SET_STRING_ELT(params,2,mkChar("KA"));
    SET_STRING_ELT(params,3,mkChar("CL"));
    SET_STRING_ELT(params,4,mkChar("Q"));
    SET_STRING_ELT(params,5,mkChar("Kin"));
    SET_STRING_ELT(params,6,mkChar("Kout"));
    SET_STRING_ELT(params,7,mkChar("EC50"));
    SET_STRING_ELT(state,0,mkChar("depot"));
    SET_STRING_ELT(normState,0,mkChar("depot"));
    _SR[0] = 0;
    SET_STRING_ELT(state,1,mkChar("centr"));
    SET_STRING_ELT(normState,1,mkChar("centr"));
    _SR[1] = 0;
    SET_STRING_ELT(state,2,mkChar("peri"));
    SET_STRING_ELT(normState,2,mkChar("peri"));
    _SR[2] = 0;
    SET_STRING_ELT(state,3,mkChar("eff"));
    SET_STRING_ELT(normState,3,mkChar("eff"));
    _SR[3] = 0;
    SET_STRING_ELT(modeln,0,mkChar("model"));
    SET_STRING_ELT(model,0,mkChar("\n   # A 4-compartment model, 3 PK and a PD (effect) compartment\n   # (notice state variable names 'depot', 'centr', 'peri', 'eff')\n\n   C2 = centr/V2;\n   C3 = peri/V3;\n   d/dt(depot) =-KA*depot;\n   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;\n   d/dt(peri)  =                    Q*C2 - Q*C3;\n   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;\n\n"));
    SET_STRING_ELT(modeln,1,mkChar("normModel"));
    SET_STRING_ELT(model,1,mkChar("C2=centr/V2;\nC3=peri/V3;\nd/dt(depot)=-KA*depot;\nd/dt(centr)=KA*depot-CL*C2-Q*C2+Q*C3;\nd/dt(peri)=Q*C2-Q*C3;\nd/dt(eff)=Kin-Kout*(1-C2/(EC50+C2))*eff;\n"));
    SET_STRING_ELT(modeln,2,mkChar("parseModel"));
    SET_STRING_ELT(model,2,mkChar("C2=centr/safe_zero(V2);\nC3=peri/safe_zero(V3);\n__DDtStateVar__[0] = _IR[0] -KA*depot;\n__DDtStateVar__[1] = _IR[1] + KA*depot-CL*C2-Q*C2+Q*C3;\n__DDtStateVar__[2] = _IR[2] + Q*C2-Q*C3;\n__DDtStateVar__[3] = _IR[3] + Kin-Kout*(1-C2/safe_zero((EC50+C2)))*eff;\n"));
    SET_STRING_ELT(modeln,3,mkChar("expandModel"));
    SET_STRING_ELT(model,3,mkChar("\n   # A 4-compartment model, 3 PK and a PD (effect) compartment\n   # (notice state variable names 'depot', 'centr', 'peri', 'eff')\n\n   C2 = centr/V2;\n   C3 = peri/V3;\n   d/dt(depot) =-KA*depot;\n   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;\n   d/dt(peri)  =                    Q*C2 - Q*C3;\n   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;\n\n"));
    SEXP ini    = PROTECT(allocVector(REALSXP,0));pro++;
    SEXP inin   = PROTECT(allocVector(STRSXP, 0));pro++;
    SET_STRING_ELT(names,0,mkChar("params"));
    SET_VECTOR_ELT(lst,  0,params);
    SET_STRING_ELT(names,1,mkChar("lhs"));
    SET_VECTOR_ELT(lst,  1,lhs);
    SET_STRING_ELT(names,2,mkChar("state"));
    SET_VECTOR_ELT(lst,  2,state);
    SET_STRING_ELT(names,3,mkChar("trans"));
    SET_VECTOR_ELT(lst,  3,tran);
    SET_STRING_ELT(names,5,mkChar("model"));
    SET_VECTOR_ELT(lst,  5,model);
    SET_STRING_ELT(names,4,mkChar("ini"));
    SET_VECTOR_ELT(lst,  4,ini);
    SET_STRING_ELT(names,6,mkChar("md5"));
    SET_VECTOR_ELT(lst,  6,mmd5);
    SET_STRING_ELT(names,7,mkChar("podo"));
    SET_VECTOR_ELT(lst,  7,ScalarLogical(0));
    SET_STRING_ELT(names,8,mkChar("dfdy"));
    SET_VECTOR_ELT(lst,  8,dfdy);
    SET_STRING_ELT(names,9,mkChar("sens"));
    SET_VECTOR_ELT(lst,  9,sens);
    SET_STRING_ELT(names,10,mkChar("fn.ini"));
    SET_VECTOR_ELT(lst,  10,fn_ini);
    SET_STRING_ELT(names,11,mkChar("state.ignore"));
    SET_VECTOR_ELT(lst,  11,stateRmS);
    SET_STRING_ELT(names,12,mkChar("solve"));
    SET_VECTOR_ELT(lst,  12,solve);
    SET_STRING_ELT(names,13,mkChar("version"));
    SET_VECTOR_ELT(lst,  13,version);
    SET_STRING_ELT(names,14,mkChar("normal.state"));
    SET_VECTOR_ELT(lst,  14,normState);
    SET_STRING_ELT(names,15,mkChar("timeId"));
    SET_VECTOR_ELT(lst,  15,timeInt);
    SET_STRING_ELT(mmd5n,0,mkChar("file_md5"));
    SET_STRING_ELT(mmd5,0,mkChar("4f09646f2b75c132ca3444f34bdc3a6e"));
    SET_STRING_ELT(mmd5n,1,mkChar("parsed_md5"));
    SET_STRING_ELT(mmd5,1,mkChar("1222155014066cc9c85935cfa10adeca"));
    SET_STRING_ELT(trann,0,mkChar("lib.name"));
    SET_STRING_ELT(tran, 0,mkChar("m1_x64"));
    SET_STRING_ELT(trann,1,mkChar("jac"));
    SET_STRING_ELT(tran,1,mkChar("fullint"));
    SET_STRING_ELT(trann,2,mkChar("prefix"));
    SET_STRING_ELT(tran, 2,mkChar("m1_x64_"));
    SET_STRING_ELT(trann,3,mkChar("dydt"));
    SET_STRING_ELT(tran, 3,mkChar("m1_x64_dydt"));
    SET_STRING_ELT(trann,4,mkChar("calc_jac"));
    SET_STRING_ELT(tran, 4,mkChar("m1_x64_calc_jac"));
    SET_STRING_ELT(trann,5,mkChar("calc_lhs"));
    SET_STRING_ELT(tran, 5,mkChar("m1_x64_calc_lhs"));
    SET_STRING_ELT(trann,6,mkChar("model_vars"));
    SET_STRING_ELT(tran, 6,mkChar("m1_x64_model_vars"));
    SET_STRING_ELT(trann,7,mkChar("ode_solver"));
    SET_STRING_ELT(tran, 7,mkChar("m1_x64_ode_solver"));
    SET_STRING_ELT(trann,8,mkChar("inis"));
    SET_STRING_ELT(tran, 8,mkChar("m1_x64_inis"));
    SET_STRING_ELT(trann,  9,mkChar("dydt_lsoda"));
    SET_STRING_ELT(tran,   9,mkChar("m1_x64_dydt_lsoda"));
    SET_STRING_ELT(trann,10,mkChar("calc_jac_lsoda"));
    SET_STRING_ELT(tran, 10,mkChar("m1_x64_calc_jac_lsoda"));
    SET_STRING_ELT(trann,11,mkChar("ode_solver_solvedata"));
    SET_STRING_ELT(tran, 11,mkChar("m1_x64_ode_solver_solvedata"));
    SET_STRING_ELT(trann,12,mkChar("ode_solver_get_solvedata"));
    SET_STRING_ELT(tran, 12,mkChar("m1_x64_ode_solver_get_solvedata"));
    SET_STRING_ELT(trann,13,mkChar("dydt_liblsoda"));
    SET_STRING_ELT(tran, 13,mkChar("m1_x64_dydt_liblsoda"));
    setAttrib(tran, R_NamesSymbol, trann);
    setAttrib(mmd5, R_NamesSymbol, mmd5n);
    setAttrib(model, R_NamesSymbol, modeln);
    setAttrib(ini, R_NamesSymbol, inin);
    setAttrib(version, R_NamesSymbol, versionn);
    setAttrib(lst, R_NamesSymbol, names);
    SEXP cls = PROTECT(allocVector(STRSXP, 1));pro++;
    SET_STRING_ELT(cls, 0, mkChar("rxModelVars"));
    classgets(lst, cls);
    _assign_ptr(lst);
    UNPROTECT(pro);
    return lst;
  } else {
    UNPROTECT(pro);
    return _mv;
  }
}

extern void m1_x64_dydt_lsoda(int *neq, double *t, double *A, double *DADT)
{
  m1_x64_dydt(neq, *t, A, DADT);
}

extern int m1_x64_dydt_liblsoda(double t, double *y, double *ydot, void *data)
{
  int *neq = (int*)(data);
  m1_x64_dydt(neq, t, y, ydot);
  return(0);
}

extern void m1_x64_calc_jac_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){
  // Update all covariate parameters
  m1_x64_calc_jac(neq, *t, A, JAC, *nrowpd);
}

//Initilize the dll to match RxODE's calls
void R_init_m1_x64 (DllInfo *info){
  // Get C callables on load; Otherwise it isn't thread safe
  _assign_ptr=(RxODE_assign_ptr) R_GetCCallable("RxODE","RxODE_assign_fn_pointers");
  _rxRmModelLib=(_rxRmModelLibType) R_GetCCallable("RxODE","rxRmModelLib");
  _rxGetModelLib=(_rxGetModelLibType) R_GetCCallable("RxODE","rxGetModelLib");
  _old_c = (RxODE_ode_solver_old_c) R_GetCCallable("RxODE","rxSolveOldC");
  _RxODE_rxAssignPtr=(_rx_asgn)R_GetCCallable("RxODE","_RxODE_rxAssignPtr");
  _rxIsCurrentC = (_rxIsCurrentC_type)R_GetCCallable("RxODE","rxIsCurrentC");
  _sumPS  = (_rxSumType) R_GetCCallable("PreciseSums","PreciseSums_sum_r");
  _prodPS = (_rxProdType) R_GetCCallable("PreciseSums","PreciseSums_prod_r");
  _prodType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_prod_get");
  _sumType=(RxODE_fn0i)R_GetCCallable("PreciseSums", "PreciseSums_sum_get");
  _ptrid=(RxODE_fn0i)R_GetCCallable("RxODE", "RxODE_current_fn_pointer_id");
  // Register the outside functions
  R_RegisterCCallable("m1_x64","m1_x64_ode_solver",       (DL_FUNC) m1_x64_ode_solver);
  R_RegisterCCallable("m1_x64","m1_x64_inis", (DL_FUNC) m1_x64_inis);
  R_RegisterCCallable("m1_x64","m1_x64_inis", (DL_FUNC) m1_x64_inis);
  R_RegisterCCallable("m1_x64","m1_x64_dydt", (DL_FUNC) m1_x64_dydt);
  R_RegisterCCallable("m1_x64","m1_x64_calc_lhs", (DL_FUNC) m1_x64_calc_lhs);
  R_RegisterCCallable("m1_x64","m1_x64_calc_jac", (DL_FUNC) m1_x64_calc_jac);
  R_RegisterCCallable("m1_x64","m1_x64_dydt_lsoda", (DL_FUNC) m1_x64_dydt_lsoda);
  R_RegisterCCallable("m1_x64","m1_x64_calc_jac_lsoda", (DL_FUNC) m1_x64_calc_jac_lsoda);
  R_RegisterCCallable("m1_x64","m1_x64_ode_solver_solvedata", (DL_FUNC) m1_x64_ode_solver_solvedata);
  R_RegisterCCallable("m1_x64","m1_x64_ode_solver_get_solvedata", (DL_FUNC) m1_x64_ode_solver_get_solvedata);
  R_RegisterCCallable("m1_x64","m1_x64_dydt_liblsoda", (DL_FUNC) m1_x64_dydt_liblsoda);
  
  static const R_CMethodDef cMethods[] = {
    {"m1_x64_ode_solver", (DL_FUNC) &m1_x64_ode_solver, 15, m1_x64_ode_solverrx_t},
    {NULL, NULL, 0, NULL}
  };
  
  R_CallMethodDef callMethods[]  = {
    {"m1_x64_model_vars", (DL_FUNC) &m1_x64_model_vars, 0},
    {NULL, NULL, 0}
  };
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
}

void R_unload_m1_x64 (DllInfo *info){
  // Free resources required for single subject solve.
  SEXP _mv = PROTECT(_rxGetModelLib("m1_x64_model_vars"));
  if (!isNull(_mv)){
    _rxRmModelLib("m1_x64_model_vars");
  }
  UNPROTECT(1);
}
