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

extern void __ODE_SOLVER_SOLVEDATA__ (rx_solve *solve){
  _solveData = solve;
}

extern rx_solve *__ODE_SOLVER_GET_SOLVEDATA__(){
  return _solveData;
}

SEXP __MODEL_VARS__();
extern void __ODE_SOLVER__(int *neq,
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
  //if (_ptrid() != __TIMEID__ ){ _assign_ptr(__MODEL_VARS__());}
  __FIX_INIS__
  _old_c(neq, _theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc);
}

static R_NativePrimitiveArgType __ODE_SOLVER__rx_t[] = {
  //*neq, *theta, *time,  *evid, *ntime, *inits,   *dose,   *ret,     *atol,  *rtol,   *stiff, *transit_abs, *nlhs, *lhs, *rc
  INTSXP,REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP
};

// CODE HERE

extern void __DYDT_LSODA__(int *neq, double *t, double *A, double *DADT)
{
  __DYDT__(neq, *t, A, DADT);
}

extern int __DYDT_LIBLSODA__(double t, double *y, double *ydot, void *data)
{
  int *neq = (int*)(data);
  __DYDT__(neq, t, y, ydot);
  return(0);
}

extern void __CALC_JAC_LSODA__(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){
  // Update all covariate parameters
  __CALC_JAC__(neq, *t, A, JAC, *nrowpd);
}

//Initilize the dll to match RxODE's calls
void __R_INIT__ (DllInfo *info){
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
  R_RegisterCCallable(__LIB_STR__,__ODE_SOLVER_STR__,       (DL_FUNC) __ODE_SOLVER__);
  R_RegisterCCallable(__LIB_STR__,"__INIS__", (DL_FUNC) __INIS__);
  R_RegisterCCallable(__LIB_STR__,"__INIS__", (DL_FUNC) __INIS__);
  R_RegisterCCallable(__LIB_STR__,"__DYDT__", (DL_FUNC) __DYDT__);
  R_RegisterCCallable(__LIB_STR__,"__CALC_LHS__", (DL_FUNC) __CALC_LHS__);
  R_RegisterCCallable(__LIB_STR__,"__CALC_JAC__", (DL_FUNC) __CALC_JAC__);
  R_RegisterCCallable(__LIB_STR__,"__DYDT_LSODA__", (DL_FUNC) __DYDT_LSODA__);
  R_RegisterCCallable(__LIB_STR__,"__CALC_JAC_LSODA__", (DL_FUNC) __CALC_JAC_LSODA__);
  R_RegisterCCallable(__LIB_STR__,"__ODE_SOLVER_SOLVEDATA__", (DL_FUNC) __ODE_SOLVER_SOLVEDATA__);
  R_RegisterCCallable(__LIB_STR__,"__ODE_SOLVER_GET_SOLVEDATA__", (DL_FUNC) __ODE_SOLVER_GET_SOLVEDATA__);
  R_RegisterCCallable(__LIB_STR__,"__DYDT_LIBLSODA__", (DL_FUNC) __DYDT_LIBLSODA__);
  
  static const R_CMethodDef cMethods[] = {
    {__ODE_SOLVER_STR__, (DL_FUNC) &__ODE_SOLVER__, 15, __ODE_SOLVER__rx_t},
    {NULL, NULL, 0, NULL}
  };
  
  R_CallMethodDef callMethods[]  = {
    {__MODEL_VARS_STR__, (DL_FUNC) &__MODEL_VARS__, 0},
    {NULL, NULL, 0}
  };
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info,FALSE);
}

void __R_UNLOAD__ (DllInfo *info){
  // Free resources required for single subject solve.
  SEXP _mv = PROTECT(_rxGetModelLib(__MODEL_VARS_STR__));
  if (!isNull(_mv)){
    _rxRmModelLib(__MODEL_VARS_STR__);
  }
  UNPROTECT(1);
}
