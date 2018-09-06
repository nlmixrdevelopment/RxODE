#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include "RxODE.h"

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

double solveLinB(rx_solve *rx, unsigned int id, double t, int linCmt, int diff1, int diff2, double d_A, double d_alpha, double d_B, double d_beta, double d_C, double d_gamma, double d_ka, double d_tlag){
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
  rx_solving_options *op = rx->op;
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
  rx_solving_options_ind *ind = &(rx->subjects[id]);
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


/* Authors: Robert Gentleman and Ross Ihaka and The R Core Team */
/* Taken directly from https://github.com/wch/r-source/blob/922777f2a0363fd6fe07e926971547dd8315fc24/src/library/stats/src/approx.c*/
/* Changed as follows:
   - Different Name
   - Use RxODE structure
   - Make inline
*/
double rx_approxP(double v, double *x, double *y, int n,
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


void _update_par_ptr(double t, unsigned int id, rx_solve *rx, int idx){
  rx_solving_options_ind *ind;
  ind = &(rx->subjects[id]);
  rx_solving_options *op = rx->op;
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
	  if (idx > 0 && idx < ind->n_all_times && t == all_times[idx]){
	    par_ptr[op->par_cov[k]-1] = cov_ptr[ind->n_all_times*k+idx];
	  } else {
            // Use the same methodology as approxfun.
            ind->ylow = cov_ptr[ind->n_all_times*k];
            ind->yhigh = cov_ptr[ind->n_all_times*k+ind->n_all_times-1];
            par_ptr[op->par_cov[k]-1] = rx_approxP(t, all_times, cov_ptr+ind->n_all_times*k, ind->n_all_times, op, ind);
	  }
        }
      }
    }
  }
}
