#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include "../inst/include/RxODE.h"
#define safe_zero(a) ((a) == 0 ? DOUBLE_EPS : (a))
#define _as_zero(a) (fabs(a) < sqrt(DOUBLE_EPS) ? 0.0 : a)
#define _as_dbleps(a) (fabs(a) < sqrt(DOUBLE_EPS) ? ((a) < 0 ? -sqrt(DOUBLE_EPS)  : sqrt(DOUBLE_EPS)) : a)


#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

// From https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
double log1mex(double a){
  if (a < M_LN2) return log(-expm1(-a));
  return(log1p(-exp(-a)));
}

void getWh(int evid, int *wh, int *cmt, int *wh100, int *whI, int *wh0);

// Linear compartment models/functions

extern int _locateDoseIndex(const double obs_time,  rx_solving_options_ind *ind){
  // Uses bisection for slightly faster lookup of dose index.
  int i, j, ij, idose;
  i = 0;
  j = ind->ndoses - 1;
  idose = ind->idose[i];
  if (obs_time < ind->all_times[idose]){
    return i;
  }
  idose = ind->idose[j];
  if (obs_time > ind->all_times[idose]){
    return j;
  }
  while(i < j - 1) { /* x[i] <= obs_time <= x[j] */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    idose = ind->idose[ij];
    if(obs_time < ind->all_times[idose])
      j = ij;
    else
      i = ij;
  }
  /* if (i == 0) return 0; */
  while(i != 0 && obs_time == ind->all_times[ind->idose[i]]){
    i--;
  }
  if (i == 0){
    while(i < ind->ndoses-2 && fabs(obs_time  - ind->all_times[ind->idose[i+1]])<= sqrt(DOUBLE_EPS)){
      i++;
    }
  }
  return i;
}

static inline double _getDur(int l, rx_solving_options_ind *ind, int backward, unsigned int *p){
  double dose = ind->dose[l];
  if (backward){
    p[0] = l-1;
    while (p[0] > 0 && ind->dose[p[0]] != -dose){
      p[0]--;
    }
    if (ind->dose[p[0]] != -dose){
      error(_("could not find a start to the infusion"));
    }
    return ind->all_times[ind->idose[l]] - ind->all_times[ind->idose[p[0]]];
  } else {
    p[0] = l+1;
    while (p[0] < ind->ndoses && ind->dose[p[0]] != -dose){
      p[0]++;
    }
    if (ind->dose[p[0]] != -dose){
      error(_("could not find an end to the infusion"));
    }
    return ind->all_times[ind->idose[p[0]]] - ind->all_times[ind->idose[l]];
  }
}


extern double getTime(int idx, rx_solving_options_ind *ind);

static inline int _locateTimeIndex(double obs_time,  rx_solving_options_ind *ind){
  // Uses bisection for slightly faster lookup of dose index.
  int i, j, ij;
  i = 0;
  j = ind->n_all_times - 1;
  if (obs_time < getTime(ind->ix[i], ind)){
    return i;
  }
  if (obs_time > getTime(ind->ix[j], ind)){
    return j;
  }
  while(i < j - 1) { /* x[i] <= obs_time <= x[j] */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    if(obs_time < getTime(ind->ix[ij], ind))
      j = ij;
    else
      i = ij;
  }
  /* if (i == 0) return 0; */
  while(i != 0 && obs_time == getTime(ind->ix[i], ind)){
    i--;
  }
  if (i == 0){
    while(i < ind->ndoses-2 && fabs(obs_time  - getTime(ind->ix[i+1], ind))<= sqrt(DOUBLE_EPS)){
      i++;
    }
  }
  return i;
}

/* Authors: Robert Gentleman and Ross Ihaka and The R Core Team */
/* Taken directly from https://github.com/wch/r-source/blob/922777f2a0363fd6fe07e926971547dd8315fc24/src/library/stats/src/approx.c*/
/* Changed as follows:
   - Different Name
   - Use RxODE structure
   - Use getTime to allow model-based changes to dose timing
   - Use getValue to ignore NA values for time-varying covariates
*/
static inline double getValue(int idx, double *y, rx_solving_options_ind *ind){
  int i = idx;
  double ret = y[ind->ix[idx]];
  if (ISNA(ret)){
    // Go backward.
    while (ISNA(ret) && i != 0){
      i--; ret = y[ind->ix[i]];
    }
    if (ISNA(ret)){
      // Still not found go forward.
      i = idx;
      while (ISNA(ret) && i != ind->n_all_times){
	i++; ret = y[ind->ix[i]];
      }
      if (ISNA(ret)){
	// All Covariates values for a single individual are NA.
	ind->allCovWarn=1;
      }
    }
  }
  return ret;
}
#define T(i) getTime(id->ix[i], id)
#define V(i) getValue(i, y, id)
double rx_approxP(double v, double *y, int n,
		  rx_solving_options *Meth, rx_solving_options_ind *id){
  /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
  int i, j, ij;

  if(!n) return R_NaN;

  i = 0;
  j = n - 1;

  /* handle out-of-domain points */
  if(v < T(i)) return id->ylow;
  if(v > T(j)) return id->yhigh;

  /* find the correct interval by bisection */
  while(i < j - 1) { /* T(i) <= v <= T(j) */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    if(v < T(ij)) j = ij; else i = ij;
    /* still i < j */
  }
  /* provably have i == j-1 */

  /* interpolation */

  if(v == T(j)) return V(j);
  if(v == T(i)) return V(i);
  /* impossible: if(T(j) == T(i)) return V(i); */

  if(Meth->kind == 1){ /* linear */
    return V(i) + (V(j) - V(i)) * ((v - T(i))/(T(j) - T(i)));
  } else { /* 2 : constant */
    return (Meth->f1 != 0.0 ? V(i) * Meth->f1 : 0.0)
      + (Meth->f2 != 0.0 ? V(j) * Meth->f2 : 0.0);
  }
}/* approx1() */

#undef T

/* End approx from R */

// getParCov first(parNo, idx=0) last(parNo, idx=ind->n_all_times-1)
double _getParCov(unsigned int id, rx_solve *rx, int parNo, int idx0){
  rx_solving_options_ind *ind;
  ind = &(rx->subjects[id]);
  rx_solving_options *op = rx->op;
  int idx=0;
  if (idx0 == NA_INTEGER){
    idx=0;
    if (ind->evid[ind->ix[idx]] == 9) idx++;
  } else if (idx0 >= ind->n_all_times) {
    return NA_REAL;
  } else {
    idx=idx0;
  }
  if (idx < 0 || idx > ind->n_all_times) return NA_REAL;
  if (op->do_par_cov){
    for (int k = op->ncov; k--;){
      if (op->par_cov[k] == parNo+1){
	double *y = ind->cov_ptr + ind->n_all_times*k;
	return y[ind->ix[idx]];
      }
    }
  }
  return ind->par_ptr[parNo];
}

void _update_par_ptr(double t, unsigned int id, rx_solve *rx, int idx){
  if (rx == NULL) error(_("solve data is not loaded"));
  if (ISNA(t)){
    rx_solving_options_ind *ind;
    ind = &(rx->subjects[id]);
    rx_solving_options *op = rx->op;
    // Update all covariate parameters
    int k;
    int ncov = op->ncov;
    if (op->do_par_cov){
      for (k = ncov; k--;){
	if (op->par_cov[k]){
	  double *y = ind->cov_ptr + ind->n_all_times*k;
	  ind->par_ptr[op->par_cov[k]-1] = getValue(idx, y, ind);
	  if (idx == 0){
	    ind->cacheME=0;
	  } else if (getValue(idx, y, ind) != getValue(idx-1, y, ind)) {
	    ind->cacheME=0;
	  }
	}
      }
    }
  } else {
    rx_solving_options_ind *ind;
    ind = &(rx->subjects[id]);
    rx_solving_options *op = rx->op;
    // Update all covariate parameters
    int k;
    int ncov = op->ncov;
    if (op->do_par_cov){
      for (k = ncov; k--;){
	if (op->par_cov[k]){
	  double *par_ptr = ind->par_ptr;
	  double *all_times = ind->all_times;
	  double *y = ind->cov_ptr + ind->n_all_times*k;
	  if (idx > 0 && idx < ind->n_all_times && t == all_times[idx]){
	    par_ptr[op->par_cov[k]-1] = getValue(idx, y, ind);
	    if (idx == 0){
	      ind->cacheME=0;
	    } else if (getValue(idx, y, ind) != getValue(idx-1, y, ind)) {
	      ind->cacheME=0;
	    }
	  } else {
	    // Use the same methodology as approxfun.
	    ind->ylow = getValue(0, y, ind);/* cov_ptr[ind->n_all_times*k]; */
	    ind->yhigh = getValue(ind->n_all_times-1, y, ind);/* cov_ptr[ind->n_all_times*k+ind->n_all_times-1]; */
	    par_ptr[op->par_cov[k]-1] = rx_approxP(t, y, ind->n_all_times, op, ind);
	    // Don't need to reset ME because solver doesn't use the
	    // times in-between.
	  }
	}
      }
    }
  }
}

void doSort(rx_solving_options_ind *ind);
void calcMtime(int solveid, double *mtime);
static inline void setLinCmt(rx_solving_options_ind *ind,
			     int linCmt, double lag, double lag2, double f, double f2,
			     double rate, double dur, double rate2, double dur2){
  ind->linCmt=linCmt;
  ind->lag = lag;
  ind->lag2 = lag2;
  ind->f = f;
  ind->f2 = f2;
  ind->rate = rate;
  ind->dur = dur;
  ind->rate2 = rate2;
  ind->dur2 = dur2;
}
// Advan-style linCmt solutions


/*
rxOptExpr(rxNorm(RxODE({
  A1=r1/ka-((r1+(-b1-A1last)*ka)*exp(-ka*t))/ka;
  A2=((r1+(-b1-A1last)*ka)*exp(-ka*t))/(ka-k20)-(((ka-k20)*r2+ka*r1+(-b2-b1-A2last-A1last)*k20*ka+(b2+A2last)*k20^2)*exp(-k20*t))/(k20*ka-k20^2)+(r2+r1)/k20}))) = 
*/

////////////////////////////////////////////////////////////////////////////////
// 1-3 oral absorption with rates
// From Richard Upton with RxODE Expression optimization (and some manual edits)
////////////////////////////////////////////////////////////////////////////////
static inline void oneCmtKaRate(double *A1, double *A2,
				double *A1last, double *A2last,
				double *t,
				double *b1, double *b2,
				double *r1, double *r2,
				double *ka, double *k20) {
  double eKa = exp(-(*ka)*(*t));
  double e20 = exp(-(*k20)*(*t));
  *A1 = (*r1)/(*ka)-(((*r1)-(*A1last)*(*ka))*eKa)/(*ka) + (*b1);
  double A21 = (((*r1)-(*A1last)*(*ka))*eKa)/((*ka)-(*k20));
  double A22=((((*ka)-(*k20))*(*r2)+(*ka)*(*r1)+(-(*A2last)-(*A1last))*(*k20)*(*ka)+(*A2last)*(*k20)*(*k20))*e20)/((*k20)*(*ka)-(*k20)*(*k20));
  double A23 =((*r2)+(*r1))/(*k20);
  *A2 = A21 - A22 + A23 + (*b2);
}


/*
Two compartment with rates in each
 */
static inline void twoCmtKaRate(double *A1, double *A2, double *A3,
				double *A1last, double *A2last, double *A3last,
				double *t,
				double *b1, double *b2,
				double *r1, double *r2,
				double *ka,  double *k20, 
				double *k23, double *k32) {
  double rx_expr_14 =exp(-(*t)*(*ka));
  *A1=(*b1)+(*r1)/safe_zero((*ka))-rx_expr_14*((*r1)-(*ka)*(*A1last))/safe_zero((*ka));
  double rx_expr_0 =R_pow(_as_dbleps((*ka)),2);
  double rx_expr_1 =(*k32);
  double rx_expr_2 =(*r1)*(*ka);
  double rx_expr_3 =R_pow(_as_dbleps((*k32)),2);
  double rx_expr_4 =2*(*k32);
  double rx_expr_7 =(*k32)*(*k20);
  double rx_expr_8 =(*k20)+(*k23);
  double rx_expr_9 =0.5*(*k32);
  double rx_expr_10 =2*rx_expr_3;
  double rx_expr_11 =rx_expr_1*(*k20);
  double rx_expr_12 =rx_expr_4*(*k20);
  double rx_expr_15 =rx_expr_8+(*k32);
  double rx_expr_17 =rx_expr_10*(*k20);
  double rx_expr_18 =(*A1last)+(*A2last);
  double rx_expr_19 =(*A2last)+(*A3last);
  double rx_expr_20 =R_pow(_as_dbleps((rx_expr_15)),2);
  *A2=(*b2)+((*r1)*(*k32)+(*r2)*(*k32))/safe_zero((rx_expr_7))+rx_expr_14*(-rx_expr_0*(*A1last)+(*r1)*(-(*k32)+(*ka))+(*ka)*(*k32)*(*A1last))/safe_zero((rx_expr_11+(*ka)*(-2*(*k32)*(*k20)/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20)))-0.5*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20)))+rx_expr_0))+exp(-0.5*(*t)*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20)))*((*ka)*(0.25*(rx_expr_18)*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),2)+rx_expr_9*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))*(-(*A1last)-(*A2last)-(*A3last)))+(*r2)*(-0.5*(*k32)*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))+(*ka)*((*k32)-0.5*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20)))+0.25*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),2))-0.125*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),3)*(*A2last)+rx_expr_2*((*k32)-0.5*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20)))+0.25*(rx_expr_19)*(*k32)*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),2))/safe_zero(((*ka)*(-(*k32)*(*k20)+0.25*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),2))+rx_expr_9*(*k20)*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),1)-0.125*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),3)))-exp(-2*(*t)*(*k32)*(*k20)/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))*(-8*(*A2last)*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),3)+(*ka)*(4*(rx_expr_18)*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),2)+rx_expr_17*(-(*A1last)-(*A2last)-(*A3last))/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))+(*r2)*((*ka)*((*k32)-rx_expr_12/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))-rx_expr_17/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20)))+4*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),2))+rx_expr_2*((*k32)-rx_expr_12/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))+4*(rx_expr_19)*(*k32)*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),2))/safe_zero(((*ka)*(rx_expr_11-4*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),2))-2*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),2)+8*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),3)));
  double rx_expr_5 =(*r1)*(*k23);
  double rx_expr_6 =(*ka)*(*k23);
  double rx_expr_13 =rx_expr_2*(*k23);
  *A3=(rx_expr_5+(*r2)*(*k23))/safe_zero((rx_expr_7))-rx_expr_14*(rx_expr_5-rx_expr_6*(*A1last))/safe_zero((rx_expr_11+(*ka)*(-2*(*k32)*(*k20)/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20)))-0.5*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20)))+rx_expr_0))+exp(-0.5*(*t)*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20)))*((*ka)*(0.25*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),2)*(*A3last)+0.5*(-(*A1last)-(*A2last))*(*k23)*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))-0.5*(rx_expr_8)*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))*(*A3last))+(*r2)*(-0.5*(*k23)*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))+rx_expr_6)-0.125*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),3)*(*A3last)+0.25*(*k23)*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),2)*(*A2last)+rx_expr_13+0.25*(rx_expr_8)*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),2)*(*A3last))/safe_zero(((*ka)*(-(*k32)*(*k20)+0.25*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),2))+rx_expr_9*(*k20)*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),1)-0.125*R_pow(_as_dbleps((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))),3)))-exp(-2*(*t)*(*k32)*(*k20)/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))*(-8*(*A3last)*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),3)+(*ka)*(4*(*A3last)*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),2)+2*(-(*A1last)-(*A2last))*(*k23)*(*k32)*(*k20)/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20)))-2*(rx_expr_8)*(*k32)*(*k20)*(*A3last)/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))+(*r2)*(rx_expr_6-2*(*k23)*(*k32)*(*k20)/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))+4*(*k23)*(*A2last)*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),2)+rx_expr_13+4*(rx_expr_8)*(*A3last)*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),2))/safe_zero(((*ka)*(rx_expr_11-4*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),2))-2*(rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),2)+8*R_pow(_as_dbleps((rx_expr_7/safe_zero((rx_expr_15-sqrt(-4*(*k32)*(*k20)+rx_expr_20))))),3)));
}


static inline void threeCmtKaRate(double *A1, double *A2, double *A3, double *A4,
				  double *A1last, double *A2last, double *A3last, double *A4last,
				  double *t,
				  double *b1, double *b2,
				  double *r1, double *r2,
				  double *ka, double *k20,
				  double *k23, double *k32,
				  double *k24, double *k42) {
  double E2=(*k20)+(*k23)+(*k24);
  double j=(*k23)+(*k20)+(*k32)+(*k42)+(*k24);
  double rx_expr_2 =(*k20)*(*k32);
  double k=(*k23)*(*k42)+rx_expr_2+(*k20)*(*k42)+(*k32)*(*k42)+(*k24)*(*k32);
  double l=rx_expr_2*(*k42);
  double m=(3*k-j*j)/3.0;
  double n=(2*j*j*j-9*j*k+27.0*l)/27.0;
  double Q=(n*n/4+m*m*m)/27.0;
  double alpha=sqrt(-Q);
  double beta=-n/2.0;
  double rho=sqrt(beta*beta+alpha*alpha);
  double theta=atan2(alpha,beta);
  double rx_expr_0 =j/3.0;
  double rx_expr_1 =0.3333333333333333333333;
  double rx_expr_3 =theta/3.0;
  double rx_expr_4 =1.732050807568877193177;
  double rx_expr_5 =R_pow(_as_dbleps(rho),(rx_expr_1));
  double rx_expr_6 =cos(rx_expr_3);
  double rx_expr_7 =sin(rx_expr_3);
  double rx_expr_8 =rx_expr_4*rx_expr_7;
  double lam1=rx_expr_0+rx_expr_5*(rx_expr_6+rx_expr_8);
  double lam2=rx_expr_0+rx_expr_5*(rx_expr_6-rx_expr_8);
  double lam3=rx_expr_0-(2*rx_expr_5*rx_expr_6);
  
  double _30 =exp(-(*t)*(*ka));
  *A1=(*b1)+(*r1)/safe_zero((*ka))-_30*((*r1)-(*ka)*(*A1last))/safe_zero((*ka));
  double _0 =(*ka)*(*ka);
  double _1 =_0*(*ka);
  double _2 =(*ka)*(*k42);
  double _3 =lam1*lam1;
  double _4 =_3*lam1;
  double _5 =_4*lam1;
  double _6 =lam2*lam2;
  double _7 =_6*lam2;
  double _8 =_7*lam2;
  double _9 =lam3*lam3;
  double _10 =_9*lam3;
  double _11 =_10*lam3;
  double _12 =(*r1)*(*k42);
  double _13 =(*r2)*(*k42);
  double _17 =(*k42)*(*k32);
  double _18 =(*k32)+(*k42);
  double _19 =(*ka)+lam1;
  double _20 =(*ka)*lam2;
  double _21 =(*ka)*lam1;
  double _22 =_0*lam1;
  double _23 =(*ka)*_3;
  double _24 =(*ka)*_6;
  double _25 =lam2*lam1;
  double _31 =_2*(*k32);
  double _40 =exp(-(*t)*lam1);
  double _41 =(_18)*(*ka);
  double _42 =exp(-(*t)*lam2);
  double _43 =_20*lam1;
  double _44 =exp(-(*t)*lam3);
  double _45 =(*A2last)+(*A4last);
  double _46 =(*A1last)+(*A2last);
  double _47 =(*A2last)+(*A3last);
  double _48 =_24*lam1;
  double _49 =lam2*(_19);
  double _50 =_25*lam3;
  double _52 =_23-_4;
  double _53 =_7*(_19);
  double _54 =_6*(-(*ka)-lam1);
  double _55 =_43*lam3;
  double _56 =(_46)*(*ka);
  double _57 =(*k42)*(_45);
  double _58 =(_47)*(*k32);
  double _64 =_17+_41;
  double _65 =_46+(*A3last);
  double _70 =_21+_49;
  double _71 =(-(*ka)-lam1-lam2)*_10;
  double _72 =lam2*(_52);
  double _73 =(*k32)*(-(*A1last)-(*A2last)-(*A3last));
  double _74 =(*k42)*(-(*A1last)-(*A2last)-(*A4last));
  double _75 =_65+(*A4last);
  double _80 =_54+_43;
  double _81 =_17*(-(*A2last)-(*A3last)-(*A4last));
  double _82 =_9*(_70);
  double _85 =_80+_7;
  double _86 =_57+_56;
  double _88 =_31*(_75);
  double _93 =lam3*(_85);
  double _99 =_82+_71;
  double _100 =_73+_74;
  double _101 =_86+_58;
  double _102 =(*ka)*(_100);
  double _106 =_99-_55;
  double _108 =_106+_11;
  double _110 =_102+_81;
  *A2=(*b2)-_30*(_0*((*k32)*(*A1last)+(*k42)*(*A1last))-_1*(*A1last)+(*r1)*(_17+(-(*k32)-(*k42))*(*ka)+_0)-_31*(*A1last))/safe_zero((_22+lam2*(-(*ka)*lam1+_0)+lam3*(-(*ka)*lam1+lam2*(-(*ka)+lam1)+_0)-_1))-_40*(_3*(_110)+_4*(_101)-_5*(*A2last)+(*r1)*(-(*ka)*_3-_31+_41*lam1)+(*r2)*(lam1*(_64)+(-(*k32)-(*k42)-(*ka))*_3-_31+_4)+_88*lam1)/safe_zero((-(*ka)*_4+_72+lam3*(_23+lam2*(-(*ka)*lam1+_3)-_4)+_5))+_42*(_6*(_110)+_7*(_101)-_8*(*A2last)+(*r1)*(-(*ka)*_6-_31+_41*lam2)+(*r2)*(lam2*(_64)+(-(*k32)-(*k42)-(*ka))*_6-_31+_7)+_88*lam2)/safe_zero((_53+_93-_48-_8))-_44*(_9*(_110)+_10*(_101)-_11*(*A2last)+(*r1)*(-(*ka)*_9-_31+_41*lam3)+(*r2)*(lam3*(_64)+(-(*k32)-(*k42)-(*ka))*_9-_31+_10)+_88*lam3)/safe_zero((_108))+(_12*(*k32)+_13*(*k32))/safe_zero((_50));
  double _14 =(*ka)*(*k23);
  double _26 =E2*(*A3last);
  double _27 =(*ka)*(*A3last);
  double _32 =_2*(*k23);
  double _33 =(*k24)*(*A3last);
  double _34 =(*k42)*(*A3last);
  double _35 =(*k23)*(*A2last);
  double _59 =(*k23)*(_45);
  double _60 =(_46)*(*k23);
  double _66 =_26+_34;
  double _76 =(*k23)*(-(*A1last)-(*A2last)-(*A4last));
  double _78 =_26+_59;
  double _89 =_78-_33;
  double _90 =_66+_60;
  double _94 =(*ka)*(_90);
  double _97 =(*k42)*(_89);
  *A3=_30*((*r1)*(-(*k42)*(*k23)+_14)-_0*(*k23)*(*A1last)+_32*(*A1last))/safe_zero((_22+lam2*(-(*ka)*lam1+_0)+lam3*(-(*ka)*lam1+lam2*(-(*ka)+lam1)+_0)-_1))+_40*(_3*(_97+_94)+_5*(*A3last)+(*r1)*(-(*ka)*(*k23)*lam1+_32)+(*r2)*((*k23)*_3+lam1*(-(*k42)*(*k23)-_14)+_32)+(-E2*(*A3last)-_35-_34-_27)*_4+_2*lam1*(-E2*(*A3last)+_76+_33))/safe_zero((-(*ka)*_4+_72+lam3*(_23+lam2*(-(*ka)*lam1+_3)-_4)+_5))-_42*(_6*(_97+_94)+_8*(*A3last)+(*r1)*(-(*ka)*(*k23)*lam2+_32)+(*r2)*((*k23)*_6+lam2*(-(*k42)*(*k23)-_14)+_32)+(-E2*(*A3last)-_35-_34-_27)*_7+_2*lam2*(-E2*(*A3last)+_76+_33))/safe_zero((_53+_93-_48-_8))+_44*(_9*(_97+_94)+_11*(*A3last)+(*r1)*(-(*ka)*(*k23)*lam3+_32)+(*r2)*((*k23)*_9+lam3*(-(*k42)*(*k23)-_14)+_32)+(-E2*(*A3last)-_35-_34-_27)*_10+_2*lam3*(-E2*(*A3last)+_76+_33))/safe_zero((_108))+(_12*(*k23)+_13*(*k23))/safe_zero((_50));
  double _15 =(*ka)*(*k24);
  double _16 =(*ka)*(*k32);
  double _28 =E2*(*A4last);
  double _29 =(*ka)*(*A4last);
  double _36 =_15*(*k32);
  double _37 =(*k23)*(*A4last);
  double _38 =(*k32)*(*A4last);
  double _39 =(*k24)*(*A2last);
  double _61 =(_47)*(*k24);
  double _62 =(_46)*(*k24);
  double _67 =_28-_37;
  double _68 =_28+_38;
  double _77 =(*k24)*(-(*A1last)-(*A2last)-(*A3last));
  double _91 =_67+_61;
  double _92 =_68+_62;
  double _95 =(*ka)*(_92);
  double _98 =(*k32)*(_91);
  *A4=_30*((*r1)*(-(*k24)*(*k32)+_15)-_0*(*k24)*(*A1last)+_36*(*A1last))/safe_zero((_22+lam2*(-(*ka)*lam1+_0)+lam3*(-(*ka)*lam1+lam2*(-(*ka)+lam1)+_0)-_1))+_40*(_3*(_98+_95)+_4*(-E2*(*A4last)-_39-_38-_29)+_5*(*A4last)+(*r1)*(_36-_15*lam1)+(*r2)*((*k24)*_3+lam1*(-(*k24)*(*k32)-_15)+_36)+_16*lam1*(-E2*(*A4last)+_37+_77))/safe_zero((-(*ka)*_4+_72+lam3*(_23+lam2*(-(*ka)*lam1+_3)-_4)+_5))-_42*(_6*(_98+_95)+_7*(-E2*(*A4last)-_39-_38-_29)+_8*(*A4last)+(*r1)*(_36-_15*lam2)+(*r2)*((*k24)*_6+lam2*(-(*k24)*(*k32)-_15)+_36)+_16*lam2*(-E2*(*A4last)+_37+_77))/safe_zero((_53+_93-_48-_8))+_44*(_9*(_98+_95)+_10*(-E2*(*A4last)-_39-_38-_29)+_11*(*A4last)+(*r1)*(_36-_15*lam3)+(*r2)*((*k24)*_9+lam3*(-(*k24)*(*k32)-_15)+_36)+_16*lam3*(-E2*(*A4last)+_37+_77))/safe_zero((_108))+((*r1)*(*k24)*(*k32)+(*r2)*(*k24)*(*k32))/safe_zero((_50));
}


////////////////////////////////////////////////////////////////////////////////
// 1-3 oral absorption with rates.
// Adapted from Richard Upton's maxima notebooks and supplementary material from
// Abuhelwa2015
////////////////////////////////////////////////////////////////////////////////
static inline void oneCmtKa(double *A1, double *A2,
			    double *A1last, double *A2last,
			    double *t,
			    double *b1, double *b2,
			    double *ka, double *k20) {
  double rx_expr_0=exp(-(*t)*(*ka));
  *A1=(*A1last)*rx_expr_0+(*b1);
  double rx_expr_1=exp(-(*t)*(*k20));
  *A2=(*A1last)*(*ka)/((*ka)-(*k20))*(rx_expr_1-rx_expr_0)+(*A2last)*rx_expr_1+(*b2);
}

static inline void twoCmtKa(double *A1, double *A2, double *A3,
			    double *A1last, double *A2last, double *A3last,
			    double *t,
			    double *b1, double *b2,
			    double *ka, double *kel,
			    double *k12, double *k21)  {
  double rxe2=exp(-(*t)*(*ka));
  *A1=(*b1)+rxe2*(*A1last);
  double rxe0=(*k12)+(*k21);
  double rxe1=(*k12)+(*kel);
  double rxe3=(*k21)*(*A2last);
  double rxe4=(*k21)*(*A3last);
  double rxe6=rxe0+(*kel);
  double rxe7=(rxe1)*(*k21);
  double rxe8=rxe6*rxe6;
  double rxe10=sqrt(-4*(-(*k12)*(*k21)+rxe7)+rxe8);
  *A2=(*b2)+(-exp(-0.5*(*t)*(rxe6-rxe10))*(-0.5*(*A2last)*(rxe6-rxe10)+rxe3+rxe4)+exp(-0.5*(*t)*(rxe6+rxe10))*(-0.5*(*A2last)*(rxe6+rxe10)+rxe3+rxe4))/(0.5*(rxe6-rxe10)-0.5*(rxe6+rxe10))+(*ka)*(rxe2*((*k21)-(*ka))/((-(*ka)+0.5*(rxe6-rxe10))*(-(*ka)+0.5*(rxe6+rxe10)))+exp(-0.5*(*t)*(rxe6-rxe10))*((*k21)-0.5*(rxe6-rxe10))/((-0.5*(rxe6-rxe10)+0.5*(rxe6+rxe10))*((*ka)-0.5*(rxe6-rxe10)))+exp(-0.5*(*t)*(rxe6+rxe10))*((*k21)-0.5*(rxe6+rxe10))/((0.5*(rxe6-rxe10)-0.5*(rxe6+rxe10))*((*ka)-0.5*(rxe6+rxe10))))*(*A1last);
  double rxe5=(*k12)*(*A2last);
  double rxe9=(rxe1)*(*A3last);
  *A3=(-exp(-0.5*(*t)*(rxe6-rxe10))*(-0.5*(*A3last)*(rxe6-rxe10)+rxe5+rxe9)+exp(-0.5*(*t)*(rxe6+rxe10))*(-0.5*(*A3last)*(rxe6+rxe10)+rxe5+rxe9))/(0.5*(rxe6-rxe10)-0.5*(rxe6+rxe10))+(*ka)*(*k12)*(*A1last)*(rxe2/((-(*ka)+0.5*(rxe6-rxe10))*(-(*ka)+0.5*(rxe6+rxe10)))+exp(-0.5*(*t)*(rxe6-rxe10))/((-0.5*(rxe6-rxe10)+0.5*(rxe6+rxe10))*((*ka)-0.5*(rxe6-rxe10)))+exp(-0.5*(*t)*(rxe6+rxe10))/((0.5*(rxe6-rxe10)-0.5*(rxe6+rxe10))*((*ka)-0.5*(rxe6+rxe10))));
}


static inline void threeCmtKa(double *A1, double *A2, double *A3, double *A4,
			      double *A1last, double *A2last, double *A3last, double *A4last,
			      double *t,
			      double *b1, double *b2,
			      double *KA, double *k20,
			      double *k23, double *k32,
			      double *k24, double *k42) {
  double E2 = (*k20)+(*k23)+(*k24);
  double E3 = (*k32);
  double E4 = (*k42);

  //calculate hybrid rate constants
  double a = E2+E3+E4;
  double b = E2*E3+E4*(E2+E3)-(*k23)*(*k32)-(*k24)*(*k42);
  double c = E2*E3*E4-E4*(*k23)*(*k32)-E3*(*k24)*(*k42);

  double a2 = a*a;
  double m = 0.333333333333333*(3.0*b - a2);
  double n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
  double Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

  double alpha = sqrt(-Q);
  double beta = -0.5*n;
  double gamma = sqrt(_as_zero(beta*beta+alpha*alpha));
  double theta = atan2(alpha,beta);
  double theta3 = 0.333333333333333*theta;
  double ctheta3 = cos(theta3);
  double stheta3 = 1.7320508075688771932*sin(theta3);
  double gamma3 = R_pow(gamma,0.333333333333333);
  
  double lambda1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
  double lambda2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
  double lambda3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);

  double B = (*A3last)*(*k32)+(*A4last)*(*k42);
  double C = E4*(*A3last)*(*k32)+E3*(*A4last)*(*k42);
  double I = (*A2last)*(*k23)*E4-(*A3last)*(*k24)*(*k42)+(*A4last)*(*k23)*(*k42);
  double J = (*A2last)*(*k24)*E3+(*A3last)*(*k24)*(*k32)-(*A4last)*(*k23)*(*k32);

  double eL1 = exp(-(*t)*lambda1);
  double eL2 = exp(-(*t)*lambda2);
  double eL3 = exp(-(*t)*lambda3);
  double eKA = exp(-(*t)*(*KA));

  double l12 = (lambda1-lambda2);
  double l13 = (lambda1-lambda3);
  double l21 = (lambda2-lambda1);
  double l23 = (lambda2-lambda3);
  double l31 = (lambda3-lambda1);
  double l32 = (lambda3-lambda2);

  double e2l1 = (E2-lambda1);
  double e2l2 = (E2-lambda2);
  double e3l1 = (E3-lambda1);
  double e3l2 = (E3-lambda2);
  double e3l3 = (E3-lambda3);
  double e4l1 = (E4-lambda1);
  double e4l2 = (E4-lambda2);
  double e4l3 = (E4-lambda3);
  
  double A2term1 = (*A2last)*(eL1*e3l1*e4l1/(l21*l31)+eL2*e3l2*e4l2/(l12*l32)+eL3*e3l3*e4l3/(l13*l23));
  
  double A2term2 = eL1*(C-B*lambda1)/(l12*l13)+eL2*(B*lambda2-C)/(l12*l23)+eL3*(B*lambda3-C)/(l13*l32);
  
  double A2term3 = (*A1last)*(*KA)*(eL1*e3l1*e4l1/(l21*l31*((*KA)-lambda1))+eL2*e3l2*e4l2/(l12*l32*((*KA)-lambda2))+eL3*e3l3*e4l3/(l13*l23*((*KA)-lambda3))+eKA*(E3-(*KA))*(E4-(*KA))/((lambda1-(*KA))*(lambda2-(*KA))*(lambda3-(*KA))));
  
  *A2 = A2term1+A2term2+A2term3 + (*b2);

  double A3term1 = (*A3last)*(eL1*e2l1*e4l1/(l21*l31)+eL2*e2l2*e4l2/(l12*l32)+eL3*(E2-lambda3)*e4l3/(l13*l23));
  
  double A3term2 = eL1*(I-(*A2last)*(*k23)*lambda1)/(l12*l13)+eL2*((*A2last)*(*k23)*lambda2-I)/(l12*l23)+eL3*((*A2last)*(*k23)*lambda3-I)/(l13*l32);
  
  double A3term3 = (*A1last)*(*KA)*(*k23)*(eL1*e4l1/(l21*l31*((*KA)-lambda1))+eL2*e4l2/(l12*l32*((*KA)-lambda2))+eL3*e4l3/(l13*l23*((*KA)-lambda3))+eKA*(E4-(*KA))/((lambda1-(*KA))*(lambda2-(*KA))*(lambda3-(*KA))));
  
  *A3 = A3term1+A3term2+A3term3;// Amount in the first-peripheral compartment

  double A4term1 = (*A4last)*(eL1*e2l1*e3l1/(l21*l31)+eL2*e2l2*e3l2/(l12*l32)+eL3*(E2-lambda3)*e3l3/(l13*l23));
  
  double A4term2 = eL1*(J-(*A2last)*(*k24)*lambda1)/(l12*l13)+eL2*((*A2last)*(*k24)*lambda2-J)/(l12*l23)+eL3*((*A2last)*(*k24)*lambda3-J)/(l13*l32);
  
  double A4term3 = (*A1last)*(*KA)*(*k24)*(eL1*e3l1/(l21*l31*((*KA)-lambda1))+eL2*e3l2/(l12*l32*((*KA)-lambda2))+eL3*e3l3/(l13*l23*((*KA)-lambda3))+eKA*(E3-(*KA))/((lambda1-(*KA))*(lambda2-(*KA))*(lambda3-(*KA))));
  *A4 = A4term1+A4term2+A4term3;

  *A1 = (*A1last)*eKA + (*b1);
}

////////////////////////////////////////////////////////////////////////////////
// 1-3 compartment bolus during infusion
////////////////////////////////////////////////////////////////////////////////
static inline void oneCmtRateSS4(double *A1, double *A1last, 
				 double *t,
				 double *b1, double *r1,
				 double *k10) {
  *A1 = (*r1)/(*k10);
}
static inline void oneCmtRate(double *A1, double *A1last, 
			      double *t,
			      double *b1, double *r1,
			      double *k10) {
  double eT = exp(-(*k10)*(*t));
  *A1 = (*r1)/(*k10)*(1-eT)+(*A1last)*eT + (*b1);
}

static inline void twoCmtRate(double *A1, double *A2, 
			      double *A1last, double *A2last,
			      double *t,
			      double *b1, double *Doserate,
			      double *k10, double *k12, double *k21) {
  double E1 = (*k10)+(*k12);
  double E2 = (*k21);

  //#calculate hybrid rate constants
  double s = E1+E2;
  double sqr = sqrt(s*s-4*(E1*E2-(*k12)*(*k21)));
  double lambda1 = 0.5*(s+sqr);
  double lambda2 = 0.5*(s-sqr);

  double eT1 = exp(-(*t)*lambda1);
  double eT2 = exp(-(*t)*lambda2);

  double l12 = (lambda1-lambda2);
  double l21 = (lambda2-lambda1);

  double c10 = ((*A1last)*E2+(*Doserate)+(*A2last)*(*k21));
  double c11 = (c10-(*A1last)*lambda1)/l21;
  double c12 = (c10-(*A1last)*lambda2)/l21;
  double A1term1 = c11*eT1 - c12*eT2;
  double A1term2 = (*Doserate)*E2*(1/(lambda1*lambda2)+eT1/(lambda1*l12)-eT2/(lambda2*l12));
  *A1 = A1term1+A1term2 + (*b1);//Amount in the central compartment

  double c20 = ((*A2last)*E1+(*A1last)*(*k12));
  double c21 = (c20-(*A2last)*lambda1)/l21;
  double c22 = (c20-(*A2last)*lambda2)/l21;
  double A2term1 = c21*eT1-c22*eT2;
  double A2term2 = (*Doserate)*(*k12)*(1/(lambda1*lambda2)+eT1/(lambda1*l12)-eT2/(lambda2*(lambda1-lambda2)));
  *A2 = A2term1+A2term2;//Amount in the peripheral compartment
}

static inline void threeCmtRate(double *A1, double *A2, double *A3,
				double *A1last, double *A2last, double *A3last,
				double *t, double *b1, double *Doserate,
				double *k10, double *k12, double *k21,
				double *k13, double *k31) {
  double E1 = (*k10)+(*k12)+(*k13);
  double E2 = (*k21);
  double E3 = (*k31);

  //#calculate hybrid rate constants
  double a = E1+E2+E3;
  double b = E1*E2+E3*(E1+E2)-(*k12)*(*k21)-(*k13)*(*k31);
  double c = E1*E2*E3-E3*(*k12)*(*k21)-E2*(*k13)*(*k31);

  double a2 = a*a;
  double m = 0.333333333333333*(3.0*b - a2);
  double n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
  double Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

  double alpha = sqrt(-Q);
  double beta = -0.5*n;
  double gamma = sqrt(_as_zero(beta*beta+alpha*alpha));
  double theta = atan2(alpha,beta);
  double theta3 = 0.333333333333333*theta;
  double ctheta3 = cos(theta3);
  double stheta3 = 1.7320508075688771932*sin(theta3);
  double gamma3 = R_pow(gamma,0.333333333333333);
  
  double lambda1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
  double lambda2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
  double lambda3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);

  double B = (*A2last)*(*k21)+(*A3last)*(*k31);
  double C = E3*(*A2last)*(*k21)+E2*(*A3last)*(*k31);
  double I = (*A1last)*(*k12)*E3-(*A2last)*(*k13)*(*k31)+(*A3last)*(*k12)*(*k31);
  double J = (*A1last)*(*k13)*E2+(*A2last)*(*k13)*(*k21)-(*A3last)*(*k12)*(*k21);

  double eL1 = exp(-(*t)*lambda1);
  double eL2 = exp(-(*t)*lambda2);
  double eL3 = exp(-(*t)*lambda3);

  double l12 = (lambda1-lambda2);
  double l13 = (lambda1-lambda3);
  double l21 = (lambda2-lambda1);
  double l23 = (lambda2-lambda3);
  double l31 = (lambda3-lambda1);
  double l32 = (lambda3-lambda2);
  
  double e1l1 = (E1-lambda1);
  double e1l2 = (E1-lambda2);
  double e1l3 = (E1-lambda3);
  double e2l1 = (E2-lambda1);
  double e2l2 = (E2-lambda2);
  double e2l3 = (E2-lambda3);
  double e3l1 = (E3-lambda1);
  double e3l2 = (E3-lambda2);
  double e3l3 = (E3-lambda3);

  double A1term1 = (*A1last)*(eL1*e2l1*e3l1/(l21*l31)+eL2*e2l2*e3l2/(l12*l32)+eL3*e2l3*e3l3/(l13*l23));
  double A1term2 = eL1*(C-B*lambda1)/(l12*l13)+eL2*(B*lambda2-C)/(l12*l23)+eL3*(B*lambda3-C)/(l13*l32);
  double A1term3 = (*Doserate)*((E2*E3)/(lambda1*lambda2*lambda3)-eL1*e2l1*e3l1/(lambda1*l21*l31)-eL2*e2l2*e3l2/(lambda2*l12*l32)-eL3*e2l3*e3l3/(lambda3*l13*l23));

  *A1 = A1term1+A1term2+A1term3 + (*b1);//Amount in the central compartment

  double A2term1 = (*A2last)*(eL1*e1l1*e3l1/(l21*l31)+eL2*e1l2*e3l2/(l12*l32)+eL3*e1l3*e3l3/(l13*l23));
  double A2term2 = eL1*(I-(*A1last)*(*k12)*lambda1)/(l12*l13)+eL2*((*A1last)*(*k12)*lambda2-I)/(l12*l23)+eL3*((*A1last)*(*k12)*lambda3-I)/(l13*l32);
  double A2term3 = (*Doserate)*(*k12)*(E3/(lambda1*lambda2*lambda3)-eL1*e3l1/(lambda1*l21*l31)-eL2*e3l2/(lambda2*l12*l32)-eL3*e3l3/(lambda3*l13*l23));

  *A2 = A2term1+A2term2+A2term3;// Amount in the first-peripheral compartment

  double A3term1 = (*A3last)*(eL1*e1l1*e2l1/(l21*l31)+eL2*e1l2*e2l2/(l12*l32)+eL3*e1l3*e2l3/(l13*l23));
  double A3term2 = eL1*(J-(*A1last)*(*k13)*lambda1)/(l12*l13)+eL2*((*A1last)*(*k13)*lambda2-J)/(l12*l23)+eL3*((*A1last)*(*k13)*lambda3-J)/(l13*l32);
  double A3term3 = (*Doserate)*(*k13)*(E2/(lambda1*lambda2*lambda3)-eL1*e2l1/(lambda1*l21*l31)-eL2*e2l2/(lambda2*l12*l32)-eL3*e2l3/(lambda3*l13*l23));

  *A3 = A3term1+A3term2+A3term3;//Amount in the second-peripheral compartment
}

////////////////////////////////////////////////////////////////////////////////
// 1-3 compartment bolus only
//
static inline void oneCmtBolusSS(double *A1, double *A1last, 
				 double *tau,
				 double *b1, double *k10) {
  double eT = exp(-(*k10)*(*tau));
  *A1 = (*b1)*(eT/(1-eT));
}
static inline void oneCmtBolus(double *A1, double *A1last, 
			       double *t,
			       double *b1, double *k10) {
  *A1 = fabs((*A1last)*exp(-(*k10)*(*t)) + (*b1));
}



static inline void twoCmtBolus(double *A1, double *A2,
			       double *A1last, double *A2last,
			       double *t, double *b1, double *k10,
			       double *k12, double *k21){
  double E1 = (*k10)+(*k12);
  double E2 = (*k21);

  double s = (*k12)+(*k21)+(*k10);
  double sqr = sqrt(s*s-4*(*k21)*(*k10));
  //calculate hybrid rate constants
  double lambda1 = 0.5*(s+sqr);
  double lambda2 = 0.5*(s-sqr);

  double eT1= exp(-(*t)*lambda1);
  double eT2= exp(-(*t)*lambda2);

  double A1term = ((((*A1last)*E2+(*A2last)*(*k21))-(*A1last)*lambda1)*eT1-(((*A1last)*E2+(*A2last)*(*k21))-(*A1last)*lambda2)*eT2)/(lambda2-lambda1);
  
  *A1 = fabs(A1term + (*b1)); //Amount in the central compartment

  double A2term = ((((*A2last)*E1+(*A1last)*(*k12))-(*A2last)*lambda1)*eT1-(((*A2last)*E1+(*A1last)*(*k12))-(*A2last)*lambda2)*eT2)/(lambda2-lambda1);
    *A2 = fabs(A2term);//            #Amount in the peripheral compartment
}

static inline void threeCmtBolus(double *A1, double *A2, double *A3,
				 double *A1last, double *A2last, double *A3last,
				 double *t, double *b1, double *k10,
				 double *k12, double *k21,
				 double *k13, double *k31){
  double E1 = (*k10)+(*k12)+(*k13);
  double E2 = (*k21);
  double E3 = (*k31);

  //calculate hybrid rate constants
  double a = E1+E2+E3;
  double b = E1*E2+E3*(E1+E2)-(*k12)*(*k21)-(*k13)*(*k31);
  double c = E1*E2*E3-E3*(*k12)*(*k21)-E2*(*k13)*(*k31);

  double a2 = a*a;
  double m = 0.333333333333333*(3.0*b - a2);
  double n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
  double Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

  double alpha = sqrt(-Q);
  double beta = -0.5*n;
  double gamma = sqrt(_as_zero(beta*beta+alpha*alpha));
  double theta = atan2(alpha,beta);

  double theta3 = 0.333333333333333*theta;
  double ctheta3 = cos(theta3);
  double stheta3 = 1.7320508075688771932*sin(theta3);
  double gamma3 = R_pow(gamma,0.333333333333333);
  double lambda1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
  double lambda2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
  double lambda3 = 0.333333333333333*a -(2*gamma3*ctheta3);

  double B = (*A2last)*(*k21)+(*A3last)*(*k31);
  double C = E3*(*A2last)*(*k21)+E2*(*A3last)*(*k31);
  double I = (*A1last)*(*k12)*E3-(*A2last)*(*k13)*(*k31)+(*A3last)*(*k12)*(*k31);
  double J = (*A1last)*(*k13)*E2+(*A2last)*(*k13)*(*k21)-(*A3last)*(*k12)*(*k21);

  double eL1 = exp(-(*t)*lambda1);
  double eL2 = exp(-(*t)*lambda2);
  double eL3 = exp(-(*t)*lambda3);

  double l12 = (lambda1-lambda2);
  double l13 = (lambda1-lambda3);
  double l21 = (lambda2-lambda1);
  double l23 = (lambda2-lambda3);
  double l31 = (lambda3-lambda1);
  double l32 = (lambda3-lambda2);
  
  double e1l1 = (E1-lambda1);
  double e1l2 = (E1-lambda2);
  double e1l3 = (E1-lambda3);
  double e2l1 = (E2-lambda1);
  double e2l2 = (E2-lambda2);
  double e2l3 = (E2-lambda3);
  double e3l1 = (E3-lambda1);
  double e3l2 = (E3-lambda2);
  double e3l3 = (E3-lambda3);

  double A1term1 = (*A1last)*(eL1*e2l1*e3l1/(l21*l31)+eL2*e2l2*e3l2/(l12*l32)+eL3*e2l3*e3l3/(l13*l23));
  double A1term2 = eL1*(C-B*lambda1)/(l12*l13)+eL2*(B*lambda2-C)/(l12*l23)+eL3*(B*lambda3-C)/(l13*l32);

  *A1 = fabs((*b1)+(A1term1+A1term2));

  double A2term1 = (*A2last)*(eL1*e1l1*e3l1/(l21*l31)+eL2*e1l2*e3l2/(l12*l32)+eL3*e1l3*e3l3/(l13*l23));
  double A2term2 = eL1*(I-(*A1last)*(*k12)*lambda1)/(l12*l13)+eL2*((*A1last)*(*k12)*lambda2-I)/(l12*l23)+eL3*((*A1last)*(*k12)*lambda3-I)/(l13*l32);

  *A2 = fabs(A2term1+A2term2);

  double A3term1 = (*A3last)*(eL1*e1l1*e2l1/(l21*l31)+eL2*e1l2*e2l2/(l12*l32)+eL3*e1l3*e2l3/(l13*l23));
  double A3term2 = eL1*(J-(*A1last)*(*k13)*lambda1)/(l12*l13)+eL2*((*A1last)*(*k13)*lambda2-J)/(l12*l23)+eL3*((*A1last)*(*k13)*lambda3-J)/(l13*l32);
  *A3 = fabs(A3term1+A3term2);
}

static inline void doAdvan(double *A,// Amounts
			   double *Alast, // Last amounts
			   double tlast, // Time of last amounts
			   double ct, // Time of the dose
			   int ncmt, // Number of compartments
			   int oral0, // Indicator of if this is an oral system
			   double *b1, // Amount of the dose in compartment #1
			   double *b2, // Amount of the dose in compartment #2
			   double *r1, // Rate in Compartment #1
			   double *r2, // Rate in Compartment #2
			   double *ka, // ka (for oral doses)
			   double *kel,  //double rx_v,
			   double *k12, double *k21,
			   double *k13, double *k31){
  double t = ct - tlast;
  if ((*r1) > DOUBLE_EPS  || (*r2) > DOUBLE_EPS){
    if (oral0){
      switch (ncmt){
      case 1: {
	oneCmtKaRate(&A[0], &A[1], &Alast[0], &Alast[1],
		     &t, b1, b2, r1, r2, ka, kel);
      } break;
      case 2: {
	twoCmtKaRate(&A[0], &A[1], &A[2],
		     &Alast[0], &Alast[1], &Alast[2],
		     &t, b1, b2, r1, r2,
		     ka,  kel, k12, k21);
      } break;
      case 3: {
	threeCmtKaRate(&A[0], &A[1], &A[2], &A[3],
		       &Alast[0], &Alast[1], &Alast[2], &Alast[3],
		       &t, b1, b2, r1, r2,
		       ka, kel, k12,  k21, k13,  k31);
      } break;
      }
    } else {
      switch (ncmt){
      case 1: {
	oneCmtRate(&A[0], &Alast[0], 
		   &t, b1, r1, kel);
	return;
      } break;
      case 2: {
	twoCmtRate(&A[0], &A[1], 
		   &Alast[0], &Alast[1],
		   &t, b1, r1,
		   kel, k12, k21);
	return;
      } break;
      case 3: {
	threeCmtRate(&A[0], &A[1], &A[2],
		     &Alast[0], &Alast[1], &Alast[2],
		     &t, b1, r1, kel,
		     k12, k21, k13, k31);
	return;
      } break;

      }
    }
  } else {
    // Bolus doses only
    if (oral0){
      switch (ncmt){
      case 1: {
	oneCmtKa(&A[0], &A[1], &Alast[0], &Alast[1],
		 &t, b1, b2, ka, kel);
      } break;
      case 2: {
	twoCmtKa(&A[0], &A[1], &A[2],
		 &Alast[0], &Alast[1], &Alast[2],
		 &t, b1, b2, ka,  kel, k12, k21);
      } break;
      case 3: {
	threeCmtKa(&A[0], &A[1], &A[2], &A[3],
		   &Alast[0], &Alast[1], &Alast[2], &Alast[3],
		   &t, b1, b2, ka, kel,
		   k12,  k21, k13,  k31);
      } break;
      }
    } else {
      // Bolus
      switch (ncmt){
      case 1: {
	oneCmtBolus(&A[0], &Alast[0],
		    &t, b1, kel);
	return;
      } break;
      case 2: {
	twoCmtBolus(&A[0], &A[1],
		    &Alast[0], &Alast[1], &t, b1,
		    kel, k12, k21);
	return;
      } break;
      case 3: {
	threeCmtBolus(&A[0], &A[1], &A[2],
		      &Alast[0], &Alast[1], &Alast[2],
		      &t, b1, kel, k12, k21, k13, k31);
      } break;
      }
    }
  }
}

double linCmtA(rx_solve *rx, unsigned int id, double t, int linCmt,
	       int i_cmt, int trans, 
	       double p1, double v1,
	       double p2, double p3,
	       double p4, double p5,
	       double d_ka, double d_tlag, double d_tlag2,
	       double d_F, double d_F2,
	       // Rate and dur can apply to depot and
	       // Therefore, only 1 model rate is possible with RxODE
	       double d_rate1, double d_dur1,
	       double d_rate2, double d_dur2) {
  rx_solving_options_ind *ind = &(rx->subjects[id]);
  rx_solving_options *op = rx->op;
  int oral0;
  oral0 = (d_ka > 0) ? 1 : 0;
  int idx = ind->idx;
  double *all_times = ind->all_times;
  double *A;
  double A0[4];
  double *Alast;
  double Alast0[4] = {0, 0, 0, 0};
  double tlast;
  double it = getTime(ind->ix[idx], ind);
  double curTime;
  if (t != it) {
    // Try to get another idx by bisection
    /* REprintf("it pre: %f", it); */
    idx = _locateTimeIndex(t, ind);
    it = getTime(ind->ix[idx], ind);
    /* REprintf("it post: %f", it); */
  }
  /* REprintf("idx: %d; solved: %d; t: %f fabs: %f\n", idx, ind->solved[idx], t, fabs(t-it)); */
  if (ind->solved[idx] && fabs(t-it) < sqrt(DOUBLE_EPS)){
    // Pull from last solved value (cached)
    A = ind->linCmtAdvan+(op->nlin)*idx;
    if (trans == 10) {
      return(A[oral0]*v1);
    } else {
      return(A[oral0]/v1);
    }
  }
  int setSolved=0;
  /* REprintf("sanity #2; %f; N: %d; idx: %d\n", ind->linCmtAdvan[0], */
  /* 	   op->nlin, idx); */
  if (t == it){
    A = ind->linCmtAdvan+(op->nlin)*idx;
    if (idx == 0) {
      Alast = Alast0;
      tlast = getTime(ind->ix[0], ind);
    } else {
      if (!ind->solved[idx-1]){
	ind->idx = idx-1;
	linCmtA(rx, id, all_times[idx-1], linCmt,
		i_cmt, trans, p1, v1, p2, p3, p4, p5,
		d_ka, d_tlag, d_tlag2, d_F, d_F2,
		// Rate and dur can apply to depot and
		// Therefore, only 1 model rate is possible with RxODE
		d_rate1, d_dur1, d_rate2, d_dur2);
	ind->idx = idx;
      }
      tlast = getTime(ind->ix[idx-1], ind);
      Alast = ind->linCmtAdvan+(op->nlin)*(idx-1);
    }
    setSolved=1;
  } else if (t > it){
    if (!ind->solved[idx])
      linCmtA(rx, id, all_times[idx], linCmt,
		  i_cmt, trans, p1, v1, p2, p3, p4, p5,
		  d_ka, d_tlag, d_tlag2, d_F, d_F2,
		  // Rate and dur can apply to depot and
		  // Therefore, only 1 model rate is possible with RxODE
		  d_rate1, d_dur1, d_rate2, d_dur2);
    tlast = getTime(ind->ix[idx], ind);
    Alast = ind->linCmtAdvan+(op->nlin)*idx;
    /* REprintf("dont assign #1\n"); */
    A = A0;
  } else {
    // Current 
    if (idx == 0){
      tlast = getTime(ind->ix[0], ind);
      Alast=Alast0;
    } else {
      if (!ind->solved[idx-1]){
	ind->idx = idx-1;
	linCmtA(rx, id, all_times[idx], linCmt,
		i_cmt, trans, p1, v1, p2, p3, p4, p5,
		d_ka, d_tlag, d_tlag2, d_F, d_F2,
		// Rate and dur can apply to depot and
		// Therefore, only 1 model rate is possible with RxODE
		d_rate1, d_dur1, d_rate2, d_dur2);
	ind->idx = idx;
      }
      tlast = getTime(ind->ix[idx-1], ind);
      Alast = ind->linCmtAdvan+(op->nlin)*(idx-1);
    }
    /* REprintf("dont assign #2; %f %f, idx: %d\n", t, it, idx); */
    A = A0;
  }
  curTime = getTime(ind->ix[idx], ind);
  unsigned int ncmt = 1;
  double rx_k=0, rx_v=0;
  double rx_k12=0;
  double rx_k21=0;
  double rx_k13=0;
  double rx_k31=0;
  double btemp, ctemp, dtemp;
  if (p5 > 0.){
    ncmt = 3;
    switch (trans){
    case 1: // cl v q vp
      rx_k = p1/v1; // k = CL/V
      rx_v = v1;
      rx_k12 = p2/v1; // k12 = Q/V
      rx_k21 = p2/p3; // k21 = Q/Vp
      rx_k13 = p4/v1; // k31 = Q2/V
      rx_k31 = p4/p5; // k31 = Q2/Vp2
      break;
    case 2: // k=p1 v=v1 k12=p2 k21=p3 k13=p4 k31=p5
      rx_k = p1;
      rx_v = v1;
      rx_k12 = p2;
      rx_k21 = p3;
      rx_k13 = p4;
      rx_k31 = p5;
      break;
    case 11:
#undef beta
#define A (1/v1)
#define B p3
#define C p5
#define alpha p1
#define beta p2
#define gamma p4
      ncmt=3;
      rx_v=1/(A+B+C);
      btemp = -(alpha*C + alpha*B + gamma*A + gamma*B + beta*A + beta*C)*rx_v;
      ctemp = (alpha*beta*C + alpha*gamma*B + beta*gamma*A)*rx_v;
      dtemp = sqrt(btemp*btemp-4*ctemp);
      rx_k21 = 0.5*(-btemp+dtemp);
      rx_k31 = 0.5*(-btemp-dtemp);
      rx_k   = alpha*beta*gamma/rx_k21/rx_k31;
      rx_k12 = ((beta*gamma + alpha*beta + alpha*gamma) -
		rx_k21*(alpha+beta+gamma) - rx_k * rx_k31 + rx_k21*rx_k21)/(rx_k31 - rx_k21);
      rx_k13 = alpha + beta + gamma - (rx_k + rx_k12 + rx_k21 + rx_k31);
      break;
    case 10:
#undef A
#define A v1
      ncmt=3;
      rx_v=1/(A+B+C);
      btemp = -(alpha*C + alpha*B + gamma*A + gamma*B + beta*A + beta*C)*rx_v;
      ctemp = (alpha*beta*C + alpha*gamma*B + beta*gamma*A)*rx_v;
      dtemp = sqrt(_as_zero(btemp*btemp-4*ctemp));
      rx_k21 = 0.5*(-btemp+dtemp);
      rx_k31 = 0.5*(-btemp-dtemp);
      rx_k   = alpha*beta*gamma/rx_k21/rx_k31;
      rx_k12 = ((beta*gamma + alpha*beta + alpha*gamma) -
		rx_k21*(alpha+beta+gamma) - rx_k * rx_k31 + rx_k21*rx_k21)/(rx_k31 - rx_k21);
      rx_k13 = alpha + beta + gamma - (rx_k + rx_k12 + rx_k21 + rx_k31);
#undef A
#undef B
#undef alpha
#undef beta
#define beta Rf_beta
      break;
    default:
      REprintf(_("invalid trans (3 cmt trans %d)\n"), trans);
      return NA_REAL;
    }
  } else if (p3 > 0.){
    ncmt = 2;
    switch (trans){
    case 1: // cl=p1 v=v1 q=p2 vp=p3
      rx_k = p1/v1; // k = CL/V
      rx_v = v1;
      rx_k12 = p2/v1; // k12 = Q/V
      rx_k21 = p2/p3; // k21 = Q/Vp
      break;
    case 2: // k=p1, v1=v k12=p2 k21=p3
      rx_k = p1;
      rx_v = v1;
      rx_k12 = p2;
      rx_k21 = p3;
      break;
    case 3: // cl=p1 v=v1 q=p2 vss=p3
      rx_k = p1/v1; // k = CL/V
      rx_v = v1;
      rx_k12 = p2/v1; // k12 = Q/V
      rx_k21 = p2/(p3-v1); // k21 = Q/(Vss-V)
      break;
    case 4: // alpha=p1 beta=p2 k21=p3
      rx_v = v1;
      rx_k21 = p3;
      rx_k = p1*p2/rx_k21; // p1 = alpha p2 = beta
      rx_k12 = p1 + p2 - rx_k21 - rx_k;
      break;
    case 5: // alpha=p1 beta=p2 aob=p3
      rx_v=v1;
      rx_k21 = (p3*p2+p1)/(p3+1);
      rx_k = (p1*p2)/rx_k21;
      rx_k12 = p1+p2 - rx_k21 - rx_k;
      break;
    case 11: // A2 V, alpha=p1, beta=p2, k21
#undef beta
#define A (1/v1)
#define B p3
#define alpha p1
#define beta p2
      ncmt=2;
      rx_v   = 1/(A+B);
      rx_k21 = (A*beta + B*alpha)*rx_v;
      rx_k   = alpha*beta/rx_k21;
      rx_k12 = alpha+beta-rx_k21-rx_k;
      break;
    case 10: // A=v1, alpha=p1, beta=p2, B=p3
      // Convert to A (right now A=v1 or A=1/v1)
#undef A
#define A v1
      ncmt=2;
      rx_v   = 1/(A + B);
      rx_k21 = (A*beta + B*alpha)*rx_v;
      rx_k   = alpha*beta/rx_k21;
      rx_k12 = alpha + beta - rx_k21 - rx_k;
      /* REprintf("A: %f, B: %f, alpha: %f, beta: %f\n", A, B, alpha, beta); */
      /* REprintf("V: %f, k10: %f, k12: %f, k21: %f\n", rx_v, rx_k, rx_k12, rx_k21); */
#undef A
#undef B
#undef alpha
#undef beta
#define beta Rf_beta
      break;
    default:
      REprintf(_("invalid trans (2 cmt trans %d)\n"), trans);
      return NA_REAL;
    }
  } else if (p1 > 0.){
    ncmt = 1;
    switch(trans){
    case 1: // cl v
      rx_k = p1/v1; // k = CL/V
      rx_v = v1;
      break;
    case 2: // k V
      rx_k = p1;
      rx_v = v1;
      break;
    case 11: // alpha V
      rx_k = p1;
      rx_v = v1;
      break;
    case 10: // alpha A
      rx_k = p1;
      rx_v = 1/v1;
      break;
    default:
      REprintf(_("invalid trans (1 cmt trans %d)\n"), trans);
      return NA_REAL;
    }
  } else {
    return NA_REAL;
  }
  int evid, wh, cmt, wh100, whI, wh0;
  evid = ind->evid[ind->ix[idx]];
  double *rate = ind->linCmtRate;
  double b1=0, b2=0, r1 = 0, r2 = 0;
  if (op->nlinR == 2){
    r1 = rate[0];
    r2 = rate[1];
  } else {
    r1 = rate[0];
  }
  int doMultiply = 0, doReplace=0;
  double amt;
  double rateAdjust= 0;
  int doRate=0;
  if (!setSolved || isObs(evid)){
    // Only apply doses when you need to set the solved system.
    // When it is an observation of course you do not need to apply the doses
  } else if (evid == 3){
    // Reset event
    Alast=Alast0;
  } else {
    getWh(evid, &wh, &cmt, &wh100, &whI, &wh0);
    /* REprintf("evid: %d; cmt: %d; %d\n", evid, cmt, linCmt); */
    int cmtOff = cmt-linCmt;
    if ((oral0 && cmtOff > 1) ||
	(!oral0 && cmtOff != 0)) {
    } else {
      amt = ind->dose[ind->ixds];
      if (!ISNA(amt) && (amt > 0) && (wh0 == 10 || wh0 == 20)) {
	// dosing to cmt
	// Steady state doses; wh0 == 20 is equivalent to SS=2 in NONMEM
	double tau = ind->ii[ind->ixds];
	// For now advance based solving to steady state (just like ODE)
	double aSave[4];
	if (wh0 == 20) {
	  doAdvan(A, Alast, tlast, // Time of last amounts
		  curTime, ncmt, oral0, &b1, &b2, &r1, &r2,
		  &d_ka, &rx_k, &rx_k12, &rx_k21,
		  &rx_k13, &rx_k31);
	  for (int i = ncmt + oral0; i--;){
	    aSave[i] = A[i];
	  }
	}
	double aLast0[4] = {0, 0, 0, 0};
	double aLast1[4] = {0, 0, 0, 0};
	double solveLast[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	// Reset all the rates
	if (op->nlinR == 2){
	  rate[0]=0.0;
	  rate[1]=0.0;
	  r1=0; r2=0;
	} else {
	  rate[0] = 0.0;
	  r1=0; r2=0;
	}
	double *aCur;
	Alast = aLast0;
	aCur  = aLast1;
	tlast = 0;
	curTime = tau;
	double tinf, r0;
	int doInf=0;
	switch (whI){
	case 0: {
	  // Get bolus dose
	  if (cmtOff == 0){
	    b1 = amt*d_F;
	  } else {
	    b2 = amt*d_F2;
	  }
	  doAdvan(aCur, Alast, tlast, // Time of last amounts
		  curTime, ncmt, oral0, &b1, &b2, &r1, &r2,
		  &d_ka, &rx_k, &rx_k12, &rx_k21,
		  &rx_k13, &rx_k31);
	  int canBreak=1;
	  for (int j = 0; j < op->maxSS; j++){
	    /* aTmp = Alast; */
	    canBreak=1;
	    Alast = aCur;
	    doAdvan(aCur, Alast, tlast, // Time of last amounts
		    curTime, ncmt, oral0, &b1, &b2, &r1, &r2,
		    &d_ka, &rx_k, &rx_k12, &rx_k21,
		    &rx_k13, &rx_k31);
	    if (j <= op->minSS -1) {
	      canBreak = 0;
	    } else {
	      for (int k = ncmt + oral0; k--;){
		if (op->RTOL*fabs(aCur[k]) + op->ATOL <= fabs(aCur[k]-solveLast[k])){
		  canBreak=0;
		}
		solveLast[k] = aCur[k];
	      }
	    }
	    for (int i = ncmt + oral0; i--;){
	      solveLast[i] = aCur[i];
	    }
	    if (canBreak){
	      break;
	    }
	  }
	  for (int i = ncmt + oral0; i--;){
	    A[i] = aCur[i];
	  }
	} break;
	case 8: // Duration is modeled
	case 9: { // Rate is modeled
	  if (whI == 9) {
	    if (cmtOff == 0)  {
	      // Infusion to central compartment with oral dosing
	      r0 = d_rate1;
	      tinf = amt/r0*d_F;
	    } else {
	      // Infusion to central compartment or depot
	      r0 = d_rate2;
	      tinf = amt/r0*d_F2;
	    }
	  } else {
	    // duration is modeled
	    if (cmtOff == 0) {
	      // With oral dosing infusion to central compartment
	      tinf = d_dur1;
	      r0 = amt/tinf*d_F;
	    } else {
	      // Infusion to compartment #1 or depot
	      tinf = d_dur2;
	      r0 = amt/tinf*d_F2;
	    }
	  }
	  doInf=1;
	} break;
	case 1:
	case 2: {
	  if (ISNA(amt)){
	  } else if (amt > 0) {
	    doInf=1;
	    unsigned int p;
	    r0 = amt;
	    tinf = _getDur(ind->ixds, ind, 1, &p);
	    if (whI == 1){
	      // Duration changes
	      if (cmtOff == 0){
		tinf *= d_F;
	      } else {
		tinf *= d_F2;
	      }
	    } else {
	      // Rate changes
	      if (cmtOff == 0){
		r0 *= d_F;
	      } else {
		r0 *= d_F2;
	      }
	    }
	  }
	}
	}
	if (doInf){
	  // Infusion steady state
	  if (tinf >= tau){
	    ind->wrongSSDur=1;
	    for (int i = ncmt + oral0; i--;){
	      A[i] += NA_REAL;
	    }
	    return NA_REAL;
	    /* error(_("tau: %f; tinf: %f; rate: %f, amt: %f, t: %f"), */
	    /* 	  tau, tinf, r0, amt, t); */
	  }
	  double r1o = 0.0;
	  if (cmtOff == 0) {
	    r1 = r0;
	  } else {
	    r2 = r0;
	  }
	  double extraT = tau-tinf;
	  //
	  Alast = aCur;
	  curTime = tinf;
	  doAdvan(aCur, Alast, tlast, // Time of last amounts
		  curTime, ncmt, oral0, &b1, &b2, &r1, &r2,
		  &d_ka, &rx_k, &rx_k12, &rx_k21,
		  &rx_k13, &rx_k31);
	  Alast = aCur;
	  curTime = extraT;
	  doAdvan(aCur, Alast, tlast, // Time of last amounts
		  curTime, ncmt, oral0, &b1, &b2, &r1o, &r1o,
		  &d_ka, &rx_k, &rx_k12, &rx_k21,
		  &rx_k13, &rx_k31);
	  int canBreak=1;
	  for (int j = 0; j < op->maxSS; j++) {
	    /* aTmp = Alast; */
	    canBreak=1;
	    Alast = aCur;
	    curTime = tinf;
	    doAdvan(aCur, Alast, tlast, // Time of last amounts
		    curTime, ncmt, oral0, &b1, &b2, &r1, &r2,
		    &d_ka, &rx_k, &rx_k12, &rx_k21,
		    &rx_k13, &rx_k31);
	    if (j <= op->minSS -1) {
	      canBreak = 0;
	    } else {
	      for (int k = ncmt + oral0; k--;){
		if (op->RTOL*fabs(aCur[k]) + op->ATOL <= fabs(aCur[k]-solveLast[k])){
		  canBreak=0;
		}
		solveLast[k] = aCur[k];
	      }
	    }
	    Alast = aCur;
	    curTime = extraT;
	    doAdvan(aCur, Alast, tlast, // Time of last amounts
		    curTime, ncmt, oral0, &b1, &b2, &r1o, &r1o,
		    &d_ka, &rx_k, &rx_k12, &rx_k21,
		    &rx_k13, &rx_k31);
	    if (j <= op->minSS -1) {
	      canBreak = 0;
	    } else {
	      for (int k = ncmt + oral0; k--;){
		if (op->RTOL*fabs(aCur[k]) + op->ATOL <= fabs(aCur[k]-solveLast[4+k])){
		  canBreak=0;
		}
		solveLast[4+k] = aCur[k];
	      }
	    }
	    if (canBreak){
	      break;
	    }
	  }
	  for (int i = ncmt + oral0; i--;){
	    A[i] = aCur[i];
	  }
	  // Turn back on the rate.
	  if (op->nlinR == 2){
	    if (cmtOff == 0){
	      rate[0]=r0;
	      rate[1]=0.0;
	    } else {
	      rate[0]=0.0;
	      rate[1]=r0;
	    }
	  } else {
	    rate[0] = r0;
	  }
	}
	// Now calculate steady state
	if (wh0 == 20) {
	  for (int i = ncmt + oral0; i--;){
	    A[i] += aSave[i];
	  }
	}
	if (setSolved) ind->solved[idx] = 1;
	return A[oral0]/rx_v;
      }
      if (wh0 == 30){
	// Turning off a compartment; Not supported put everything to NaN
	for (int i = ncmt + oral0; i--;){
	  A[i] = R_NaN;
	  if (setSolved) ind->solved[idx] = 1;
	  return R_NaN;
	}
      }
      if (wh0 == 40){
	// Steady state infusion
	// Now advance to steady state dosing
	// These are easy to solve
      }
      // dosing to cmt
      amt = ind->dose[ind->ixds];
      switch (whI){
      case 0: { // Bolus dose
	// base dose
	if (cmtOff == 0){
	  b1 = amt*d_F;
	} else {
	  b2 = amt*d_F2;
	}
      } break;
      case 4: {
	doReplace = cmtOff+1;
      } break;
      case 5: {
	doMultiply= cmtOff+1;
      } break;
      case 9: // modeled rate.
      case 8: { // modeled duration. 
	// Rate already calculated and saved in the next dose record
	rateAdjust = -ind->dose[ind->ixds+1];
	doRate = cmtOff+1;
	//InfusionRate[cmt] -= dose[ind->ixds+1];
      } break;
      case 1: // Begin infusion
      case 7: // End modeled rate
      case 6: { // end modeled duration
	// If cmt is off, don't remove rate....
	// Probably should throw an error if the infusion rate is on still.
	rateAdjust = amt; // Amt is negative when turning off
	doRate = cmtOff+1;
      }	break;
      case 2: {
	// In this case bio-availability changes the rate, but the duration remains constant.
	// rate = amt/dur
	if (cmtOff == 0){
	  rateAdjust = amt*d_F;
	} else {
	  rateAdjust = amt*d_F2;
	}
	doRate = cmtOff+1;
      } break;
      }
    }
  }
  /* REprintf("evid: %d; wh: %d; cmt: %d; wh100: %d; whI: %d; wh0: %d; %f\n", */
  /* 	   evid, wh, cmt, wh100, whI, wh0, A[oral0]); */
  /* REprintf("curTime: t:%f, it: %f curTime:%f, tlast: %f, b1: %f ", t, it, curTime, tlast, b1); */
  doAdvan(A, Alast, tlast, // Time of last amounts
	  curTime, ncmt, oral0, &b1, &b2, &r1, &r2,
	  &d_ka, &rx_k, &rx_k12, &rx_k21,
	  &rx_k13, &rx_k31);
  if (doReplace){
    A[doReplace-1] = amt;
  } else if (doMultiply){
    A[doMultiply-1] *= amt;
  } else if (doRate){
    rate[doRate-1] += rateAdjust; 
  } 
  if (setSolved) ind->solved[idx] = 1;
  /* REprintf("%f,%f,%f\n", A[oral0], rx_v, A[oral0]/rx_v); */
  return A[oral0]/rx_v;
}
