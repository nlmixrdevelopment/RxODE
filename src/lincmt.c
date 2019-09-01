#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include "../inst/include/RxODE.h"

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
  /* if (idose < 0 || idose >= ind->n_all_times) error("idose corrupted #1."); */
  if (obs_time < ind->all_times[idose]){
    return i;
  }
  idose = ind->idose[j];
  /* if (idose < 0 || idose >= ind->n_all_times) error("idose corrupted #2"); */
  if (obs_time > ind->all_times[idose]){
    return j;
  }
  while(i < j - 1) { /* x[i] <= obs_time <= x[j] */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    idose = ind->idose[ij];
    /* if (idose < 0 || idose >= ind->n_all_times) error("idose corrupted #3"); */
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

static inline double _getDur(int l, rx_solving_options_ind *ind, int backward){
  double dose = ind->dose[l];
  if (backward){
    int p = l-1;
    while (p > 0 && ind->dose[p] != -dose){
      p--;
    }
    if (ind->dose[p] != -dose){
      error("Could not find a start to the infusion.  Check the event table.");	
    }
    return ind->all_times[ind->idose[l]] - ind->all_times[ind->idose[p]];
  } else {
    int p = l+1;
    while (p < ind->ndoses && ind->dose[p] != -dose){
      p++;
    }
    if (ind->dose[p] != -dose){
      error("Could not find an end to the infusion.  Check the event table.");
    }
    return ind->all_times[ind->idose[p]] - ind->all_times[ind->idose[l]];
  }
}


extern double getTime(int idx, rx_solving_options_ind *ind);

static inline double linCmtAA(rx_solve *rx, unsigned int id, double t, int linCmt,
			      int i_cmt, int trans,
			      double p1, double v1,
			      double p2, double p3,
			      double p4, double p5,
			      double d_ka, double d_tlag, double d_tlag2, double d_F, double d_F2,
			      // Rate and dur can only apply to central compartment even w/ oral dosing
			      // Therefore, only 1 model rate is possible with RxODE
			      double d_rate, double d_dur){
  double rx_k=0;
  double rx_v=0;

  double rx_k12=0;
  double rx_k21=0;
  
  double rx_k13=0;
  double rx_k31=0;
  double d_alpha = 0;
  double d_A = 0;
  double d_A2 = 0;
  double d_beta = 0;
  double d_B = 0;
  double d_B2 = 0;
  double d_gamma = 0;
  double d_C = 0;
  double d_C2 = 0;
  if (trans >= 10){
    // Direct translation
    if (trans == 11){
      d_A = 1/v1;
    } else {
      d_A = v1;
    }
    d_alpha = p1;
    d_beta = p2;
    d_B = p3;
    d_gamma = p4;
    d_C = p5;
    if (d_ka > 0){
      if (d_gamma > 0){
	d_A2 = d_A;
	d_B2 = d_B;
	d_C2 = d_C;
	d_A = d_ka / (d_ka - d_alpha) * d_A;
	d_B = d_ka / (d_ka - d_beta) * d_B;
	d_C = d_ka / (d_ka - d_gamma) * d_C;
      } else if (d_beta > 0){
	d_A2 = d_A;
	d_B2 = d_B;
	d_A = d_ka / (d_ka - d_alpha) * d_A;
	d_B = d_ka / (d_ka - d_beta) * d_B;
      } else {
	d_A2 = d_A;
	d_A = d_ka / (d_ka - d_alpha) * d_A;
      }
    }
  } else {
    if (i_cmt == 1){
      switch(trans){
      case 1: // cl v
	rx_k = p1/v1; // k = CL/V
	rx_v = v1;
	break;
      case 2: // k V
	rx_k = p1;
	rx_v = v1;
	break;
      default:
	error("invalid trans (1 cmt trans %d).", trans);
      }
      if (d_ka > 0){
	d_alpha = rx_k;
	d_A = d_ka / (d_ka - d_alpha) / rx_v;
	d_A2 = 1.0 / rx_v;
      } else {
	d_alpha = rx_k;
	d_A = 1.0 / rx_v;
      }
    } else if (i_cmt == 2){
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
      default:
	error("invalid trans (2 cmt trans %d).", trans);
      }
      double rx_tmp = rx_k12+rx_k21+rx_k;
      d_beta = 0.5 * (rx_tmp - sqrt(rx_tmp * rx_tmp - 4.0 * rx_k21 * rx_k));
      d_alpha = rx_k21 * rx_k / d_beta;
      if (d_ka > 0){
	d_A = d_ka / (d_ka - d_alpha) * (d_alpha - rx_k21) / (d_alpha - d_beta) / rx_v;
	d_B = d_ka / (d_ka - d_beta) * (d_beta - rx_k21) / (d_beta - d_alpha) / rx_v;
	d_A2 = (d_alpha - rx_k21) / (d_alpha - d_beta) / rx_v;
	d_B2 = (d_beta - rx_k21) / (d_beta - d_alpha) / rx_v;
      } else {
	d_A = (d_alpha - rx_k21) / (d_alpha - d_beta) / rx_v;
	d_B = (d_beta - rx_k21) / (d_beta - d_alpha) / rx_v;
	d_A2 = 0;
	d_B2 = 0;
      }
    } else if (i_cmt == 3){
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
      default:
	error("invalid trans (3 cmt trans %d).", trans);
      }
      double rx_a0 = rx_k * rx_k21 * rx_k31;
      double rx_a1 = rx_k * rx_k31 + rx_k21 * rx_k31 + rx_k21 * rx_k13 + rx_k * rx_k21 + rx_k31 * rx_k12;
      double rx_a2 = rx_k + rx_k12 + rx_k13 + rx_k21 + rx_k31;
      double rx_p = rx_a1 - rx_a2 * rx_a2 / 3.0;
      double rx_q = 2.0 * rx_a2 * rx_a2 * rx_a2 / 27.0 - rx_a1 * rx_a2 /3.0 + rx_a0;
      double rx_r1 = sqrt(-rx_p * rx_p * rx_p / 27.0);
      double rx_r2 = 2 * pow(rx_r1,1.0/3.0);
      double rx_theta = acos(-rx_q / (2.0 * rx_r1)) / 3.0;
      d_alpha = -(cos(rx_theta) * rx_r2 - rx_a2 / 3.0);
      d_beta = -(cos(rx_theta + M_2PI/3.0) * rx_r2 - rx_a2 / 3.0);
      d_gamma = -(cos(rx_theta + 4.0 / 3.0 * M_PI) * rx_r2 - rx_a2 / 3.0);
      d_A = (rx_k21 - d_alpha) * (rx_k31 - d_alpha) / (d_alpha - d_beta) / (d_alpha - d_gamma) / rx_v;
      d_B = (rx_k21 - d_beta) * (rx_k31 - d_beta) / (d_beta - d_alpha) / (d_beta - d_gamma) / rx_v;
      d_C = (rx_k21 - d_gamma) * (rx_k31 - d_gamma) / (d_gamma - d_alpha) / (d_gamma - d_beta) / rx_v;
      if (d_ka > 0){
	d_A2 = d_A;
	d_B2 = d_B;
	d_C2 = d_C;
	d_A = d_ka / (d_ka - d_alpha) * d_A;
	d_B = d_ka / (d_ka - d_beta) * d_B;
	d_C = d_ka / (d_ka - d_gamma) * d_C;
      }
    } else {
      error("Only 1-3 compartment linCmt() are supported");
    }
  }
  unsigned int ncmt = 1;
  double beta1=0, gamma1=0, alpha1=0;
  double alpha = d_alpha;
  double beta = d_beta;
  double gamma = d_gamma;
  double ka = d_ka;
  double tlag = d_tlag;
  double F = d_F;
  double A = d_A;
  double B = d_B;
  double C = d_C;
  rx_solving_options *op = rx->op;
  if (d_gamma > 0.){
    ncmt = 3;
    if (op->linLog){
      gamma1 = -log(gamma); // 1/gamma log(1)-log(gamma) = -log(gamma)
      beta1 = -log(beta);
      alpha1 = -log(alpha);
    } else {
      gamma1 = 1.0/gamma; // 1/gamma log(1)-log(gamma) = -log(gamma)
      beta1 = 1.0/beta;
      alpha1 = 1.0/alpha;
    }
    
  } else if (d_beta > 0.){
    ncmt = 2;
    if (op->linLog){
      beta1 = -log(beta);
      alpha1 = -log(alpha);
    } else {
      beta1 = 1.0/beta;
      alpha1 = 1.0/alpha;
    }
  } else if (d_alpha > 0.){
    ncmt = 1;
    if (op->linLog){
      alpha1 = -log(alpha);
    } else {
      alpha1 = 1.0/alpha;
    }
  } else {
    return 0.0;
    //error("You need to specify at least A(=%f) and alpha (=%f). (@t=%f, d1=%d, d2=%d)", d_A, d_alpha, t, diff1, diff2);
  }
  double ATOL = op->ATOL; //absolute error
  double RTOL = op->RTOL; //relative error
  int oral0, oral, cmt;
  oral0 = (ka > 0) ? 1 : 0;
  double ret = 0,cur=0;
  unsigned int m=0, l = 0, p = 0;
  int evid, wh, wh100, whI, wh0;
  double thisT = 0.0, tT = 0.0, res, t1, t2, tinf = 0, dose = 0, tau = 0, expr1;
  double logRate=0, rate=0, tmp;
  rx_solving_options_ind *ind = &(rx->subjects[id]);
  // don't need to adjust based on tlag t is the most conservative.
  // When tadr - tlag < 0 ignore the dose.
  m = _locateDoseIndex(t, ind);
  for(l=m+1; l--;){// Optimized for loop as https://www.thegeekstuff.com/2015/01/c-cpp-code-optimization/
    cur=0;
    //super-position
    evid = ind->evid[ind->idose[l]];
    if (evid == 3) return ret; // Was a reset event.
    getWh(evid, &wh, &cmt, &wh100, &whI, &wh0);
    dose = ind->dose[l];
    oral = oral0;
    A = d_A;
    B = d_B;
    C = d_C;
    tlag = d_tlag;
    F = d_F;
    if (cmt != linCmt){
      if (!oral){
	continue;
      } else if (cmt != linCmt+1) continue;
      oral = 0; // Switch to "central" compartment
      A = d_A2;
      B = d_B2;
      C = d_C2;
      tlag = d_tlag2;
      F = d_F2;
    }
    switch(whI){
    case 7:
      continue;
    case 6:
      continue;
    case 8: // Duration is modeled
    case 9: // Rate is modeled
      tT = t - ind->all_times[ind->idose[l]] ;
      thisT = tT - tlag;
      tau = ind->ii[l];
      if (whI == 9){
	if (op->linLog){
	  tinf  = exp(log(dose)-log(d_rate));
	  logRate  = log(d_rate);
	} else {
	  tinf  = dose/d_rate;
	  rate  = d_rate;
	}
      } else {
	if (op->linLog){
	  tinf  = d_dur;
	  logRate  = log(dose)-log(d_dur);
	} else {
	  tinf  = d_dur;
	  rate  = dose/d_dur;
	}
	
      }
      dose=NA_REAL;
    case 2:
    case 1:
      if (oral) error("Infusions to depot are not possible with the linear solved system");
      if (wh0 == 30){
	error("You cannot turn off a compartment with a solved system.");
      }
      // Steady state
      if (ISNA(dose)){
      } else if (dose > 0){
	// During infusion
	tT = t - ind->all_times[ind->idose[l]] ;
	thisT = tT - tlag;
	tinf = _getDur(ind->ixds, ind, 0);
	tau = ind->ii[l];
	if (op->linLog){
	  logRate = log(dose);
	} else {
	  rate  = dose;
	}
	if (tT >= tinf) continue;
      } else {
	// After  infusion
	tinf = _getDur(ind->ixds, ind, 1);
	tau = ind->ii[p];
	tT = t - ind->all_times[ind->idose[p]];
	thisT = tT -tlag;
	if (op->linLog){
	  logRate  = log(-dose);
	} else {
	  rate  = -dose;
	}
      }
      if (thisT < 0) continue;
      if (F <= 0) error("Bioavailability cannot be negative or zero.");
      if (whI == 1){ // Duration changes
	tinf *=F;
      } else { // Rate Changes
	if (op->linLog){
	  logRate += log(F);
	} else {
	  rate *= F;
	}
      }
      if (wh0 == 10 || wh0 == 20){
	if (tinf >= tau){
	  error("Infusion time greater then inter-dose interval, ss cannot be calculated.");
	} 
	if (thisT < tinf){ // during infusion
	  if (op->linLog){
	    expr1= logRate+alpha1;
	    cur += A*exp(expr1+log1mex(alpha*thisT))+
	      A*exp(expr1-alpha*tau+log1mex(alpha*tinf)-alpha*(thisT-tinf)-log1mex(alpha*tau));
	    if (ncmt >= 2){
	      expr1= logRate+beta1;
	      cur += B*exp(expr1+log1mex(beta*thisT))+
		B*exp(expr1-beta*tau+log1mex(beta*tinf)-beta*(thisT-tinf)-log1mex(beta*tau));
	      if (ncmt >= 3){
		expr1= logRate+gamma1;
		cur += C*exp(expr1+log1mex(gamma*thisT))+
		  C*exp(expr1-gamma*tau+log1mex(gamma*tinf)-gamma*(thisT-tinf)-log1mex(gamma*tau));
	      }
	    }
	  } else {
	    cur += rate*A*alpha1*((1-exp(-alpha*thisT))+
				  exp(-alpha*tau)*(1-exp(-alpha*tinf))*
				  exp(-alpha*(thisT-tinf))/(1-exp(-alpha*tau)));
	    if (ncmt >= 2){
	      cur += rate*B*beta1*((1-exp(-beta*thisT))+
				   exp(-beta*tau)*(1-exp(-beta*tinf))*
				   exp(-beta*(thisT-tinf))/
				   (1-exp(-beta*tau)));
	      if (ncmt >= 3){
		cur += rate*C*gamma1*((1-exp(-gamma*thisT))+
				      exp(-gamma*tau)*(1-exp(-gamma*tinf))*exp(-gamma*(thisT-tinf))/
				      (1-exp(-gamma*tau)));
	      }
	    }
	  }
	  if (wh0 == 10) return (ret+cur);
	} else { // after infusion
	  if (op->linLog){
	    expr1=logRate+alpha1;
	    cur += A*exp(expr1+log1mex(alpha*tinf)-alpha*(thisT-tinf)-log1mex(alpha*tau));
	    if (ncmt >= 2){
	      expr1=logRate+beta1;
	      cur += B*exp(expr1+log1mex(beta*tinf)-beta*(thisT-tinf)-log1mex(beta*tau));
	      if (ncmt >= 3){
		expr1=logRate+gamma1;
		cur += C*exp(expr1+log1mex(gamma*tinf)-gamma*(thisT-tinf)-log1mex(gamma*tau));
	      }
	    }
	  } else {
	    cur += rate*A*alpha1*((1-exp(-alpha*tinf))*exp(-alpha*(thisT-tinf))/(1-exp(-alpha*tau)));
	    if (ncmt >= 2){
	      cur += rate*B*beta1*((1-exp(-beta*tinf))*exp(-beta*(thisT-tinf))/(1-exp(-beta*tau)));
	      if (ncmt >= 3){
		cur += rate*C*gamma1*((1-exp(-gamma*tinf))*exp(-gamma*(thisT-tinf))/(1-exp(-gamma*tau)));
	      }
	    }
	  }
	  if (wh0 == 10) return (ret+cur);
	}
      } else {
	t1  = ((thisT < tinf) ? thisT : tinf);        //during infusion
	t2  = ((thisT > tinf) ? thisT - tinf : 0.0);  // after infusion
	if (op->linLog){
	  cur +=  A*exp(logRate+alpha1+log1mex(alpha*t1)-alpha*t2);
	  if (ncmt >= 2){
	    cur +=  B*exp(logRate+beta1+log1mex(beta*t1)-beta*t2);
	    if (ncmt >= 3){
	      cur +=  C*exp(logRate+gamma1+log1mex(gamma*t1)-gamma*t2);
	    }
	  }
	} else {
	  cur +=  rate*A*alpha1*(1.0-exp(-alpha*t1))*exp(-alpha*t2);
	  if (ncmt >= 2){
	    cur +=  rate*B*beta1*(1.0-exp(-beta*t1))*exp(-beta*t2);
	    if (ncmt >= 3){
	      cur +=  rate*C*gamma1*(1.0-exp(-gamma*t1))*exp(-gamma*t2);
	    }
	  }
	}
      }
      break;
    case 0:
      if (wh0 == 10 || wh0 == 20){
	// steady state
	tT = t - ind->all_times[ind->idose[l]];
	thisT = tT -tlag;
	if (thisT < 0) continue;
	tau = ind->ii[l];
	if (op->linLog){
	  res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	  cur += dose*F*A*(exp(-alpha*thisT)/(1-exp(-alpha*tau))-res);
	  if (ncmt >= 2){
	    cur +=  dose*F*B*(exp(-beta*thisT)/(1-exp(-beta*tau))-res);
	    if (ncmt >= 3){
	      cur += dose*F*C*(exp(-gamma*thisT)/(1-exp(-gamma*tau))-res);
	    }
	  }
	} else {
	  expr1 = log(dose)+log(F);
	  res = ((oral == 1) ? exp(expr1-ka*thisT-log1mex(ka*tau)) : 0.0);
	  cur += A*(exp(expr1-alpha*thisT-log1mex(alpha*tau))-res);
	  if (ncmt >= 2){
	    cur +=  B*(exp(expr1-beta*thisT-log1mex(beta*tau))-res);
	    if (ncmt >= 3){
	      cur += C*(exp(expr1-gamma*thisT-log1mex(gamma*tau))-res);
	    }
	  }
	}
	// ss=1 is equivalent to a reset + ss dose
	if (wh0 == 10) return(ret+cur);
      } else if (wh0 == 30) {
	error("You cannot turn off a compartment with a solved system.");
      } else {
	tT = t - ind->all_times[ind->idose[l]];
	thisT = tT -tlag;
	if (thisT < 0) continue;
	if (op->linLog){
	  expr1 = log(dose)+log(F);
	  res = ((oral == 1) ? exp(expr1-ka*thisT) : 0.0);
	  cur +=  A*(exp(expr1-alpha*thisT)-res);
	  if (ncmt >= 2){
	    cur +=  B*(exp(expr1-beta*thisT)-res);
	    if (ncmt >= 3){
	      cur += C*(exp(expr1-gamma*thisT)-res);
	    }
	  }
	} else {
	  res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	  cur +=  dose*F*A*(exp(-alpha*thisT)-res);
	  if (ncmt >= 2){
	    cur +=  dose*F*B*(exp(-beta*thisT)-res);
	    if (ncmt >= 3){
	      cur += dose*F*C*(exp(-gamma*thisT)-res);
	    }
	  }
	}
	
      }
      break;
    default:
      error("Invalid evid in linear solved system.");
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
    tmp = fabs(ret+cur);
    if (fabs(cur) < RTOL*tmp+ATOL){
      ret=ret+cur;
      break;
    }
    ret = ret+cur;
  } //l
  return ret;
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


void _update_par_ptr(double t, unsigned int id, rx_solve *rx, int idx){
  if (rx == NULL) error("solve data is not loaded.");
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

// This updateds Alast and tlast by adding a bolus to the ""
static inline void realizeBolus(double *Alast, // Last amounts
				double tlast, // Time of last amounts
				int ncmt, // Number of compartments
				int oral0, // Indicator of if this is an oral system
				int cmtOff, // Compartment offest for dose
				double amt, // Amount of the dose
				double Doserate, // Rate
				double ct, // Time of the dose
				double KA, // ka (for oral doses)
				double k10,  //double rx_v,
				double k12, double k21,
				double k13, double k31){
  double t = ct - tlast;
  if (Doserate > 0){
    if (oral0){
      // FIXME bolus with infusion in oral one-compartment model
      error("Mixed oral and iv infusions are not supported with advan compartments");
    } else {
      // Bolus with infusion
      if (ncmt == 1){
	Alast[0] = Doserate/k10*(1-exp(-t*k10))+Alast[0]*exp(-t*k10)+amt;
      } else if (ncmt == 2){
	double k20 = 0;
	double E1 = k10+k12;
	double E2 = k21+k20;
	//calculate hybrid rate constants
	double e1e2 = E1+E2;
	double sqr = sqrt(e1e2*e1e2-4*(E1*E2-k12*k21));
	double lambda1 = 0.5*(e1e2+sqr);
	double lambda2 = 0.5*(e1e2-sqr);
	double A1term1 = (((Alast[0]*E2+Doserate+Alast[1]*k21)-Alast[0]*lambda1)*exp(-t*lambda1)-((Alast[0]*E2+Doserate+Alast[1]*k21)-Alast[0]*lambda2)*exp(-t*lambda2))/(lambda2-lambda1);
	double A1term2 = Doserate*E2*(1/(lambda1*lambda2)+exp(-t*lambda1)/(lambda1*(lambda1-lambda2))-exp(-t*lambda2)/(lambda2*(lambda1-lambda2)));

	double A2term1 = (((Alast[1]*E1+Alast[0]*k12)-Alast[1]*lambda1)*exp(-t*lambda1)-((Alast[1]*E1+Alast[0]*k12)-Alast[1]*lambda2)*exp(-t*lambda2))/(lambda2-lambda1);
	double A2term2 = Doserate*k12*(1/(lambda1*lambda2)+exp(-t*lambda1)/(lambda1*(lambda1-lambda2))-exp(-t*lambda2)/(lambda2*(lambda1-lambda2)));

	Alast[0] = amt+A1term1+A1term2;
	Alast[1] = A2term1+A2term2;
      } else if (ncmt == 3){
	double k20 = 0;
	double k30 = 0;
	double E1 = k10+k12+k13;
	double E2 = k21+k20;
	double E3 = k31+k30;

	//calculate hybrid rate constants
	double a = E1+E2+E3;
	double b = E1*E2+E3*(E1+E2)-k12*k21-k13*k31;
	double c = E1*E2*E3-E3*k12*k21-E2*k13*k31;

	double m = (3.*b - a*a)/3.;
	double n = (2.*a*a*a - 9.*a*b + 27.*c)/27.;
	double Q = (n*n)/4. + (m*m*m)/27.;

	double alpha = sqrt(-1*Q);
	double beta = -1*n/2.;
	double gamma = sqrt(beta*beta+alpha*alpha);
	double theta = atan2(alpha,beta);

	double g13 = pow(gamma, 1./3.);
	double c3 = cos(theta/3.);
	double s3 = sqrt(3.)*sin(theta/3.);
	double lambda1 = a/3 + g13*(c3 + s3);
	double lambda2 = a/3 + g13*(c3 - s3);
	double lambda3 = a/3 -(2*g13*c3);

	double B = Alast[1]*k21+Alast[2]*k31;
	double C = E3*Alast[1]*k21+E2*Alast[2]*k31;
	double I = Alast[0]*k12*E3-Alast[1]*k13*k31+Alast[2]*k12*k31;
	double J = Alast[0]*k13*E2+Alast[1]*k13*k21-Alast[2]*k12*k21;

	double A1term1 = Alast[0]*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	double A1term2 = exp(-t*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2));
	double A1term3 = Doserate*((E2*E3)/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));


	double A2term1 = Alast[1]*(exp(-t*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	double A2term2 = exp(-t*lambda1)*(I-Alast[0]*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(Alast[0]*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(Alast[0]*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2));
	double A2term3 = Doserate*k12*(E3/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));


	double A3term1 = Alast[2]*(exp(-t*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	double A3term2 = exp(-t*lambda1)*(J-Alast[0]*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(Alast[0]*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(Alast[0]*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2));
	double A3term3 = Doserate*k13*(E2/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E2-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)));

	Alast[0] = A1term1+A1term2+A1term3+amt;
	Alast[1] = A2term1+A2term2+A2term3;
	Alast[2] = A3term1+A3term2+A3term3;
      }
    }
  } else {
    if (oral0){
      //
      if (ncmt == 1){
	double A2last = Alast[0]*KA/(KA-k10)*(exp(-t*k10)-exp(-t*KA))+Alast[1]*exp(-t*k10);
	double A1last = Alast[0]*exp(-t*KA);
	if (cmtOff == 0){
	  Alast[1] = A2last; //Amount in the central compartment
	  Alast[0] = A1last + amt;  //Amount in the absorption compartment
	} else {
	  Alast[1] = A2last + amt; //Amount in the central compartment
	  Alast[0] = A1last;  //Amount in the absorption compartment
	}
      } else if (ncmt == 2) {
	double k30 = 0;
	double k23 = k12;
	double k32 = k21;
	double k20 = k10;
	double E2 = k20+k23;
	double E3 = k32+k30;

	double e2e3 = E2+E3;
	double sqr = sqrt(e2e3*e2e3-4*(E2*E3-k23*k32));
	//calculate hybrid rate constants
	double lambda1 = 0.5*(e2e3+sqr);
	double lambda2 = 0.5*(e2e3-sqr);
	
	double A2term1 = (((Alast[1]*E3+Alast[2]*k32)-Alast[1]*lambda1)*exp(-t*lambda1)-((Alast[1]*E3+Alast[2]*k32)-Alast[1]*lambda2)*exp(-t*lambda2))/(lambda2-lambda1);
	double A2term2 = Alast[0]*KA*(exp(-t*KA)*(E3-KA)/((lambda1-KA)*(lambda2-KA))+exp(-t*lambda1)*(E3-lambda1)/((lambda2-lambda1)*(KA-lambda1))+exp(-t*lambda2)*(E3-lambda2)/((lambda1-lambda2)*(KA-lambda2)));

	double A3term1 = (((Alast[2]*E2+Alast[1]*k23)-Alast[2]*lambda1)*exp(-t*lambda1)-((Alast[2]*E2+Alast[1]*k23)-Alast[2]*lambda2)*exp(-t*lambda2))/(lambda2-lambda1);
	double A3term2 = Alast[0]*KA*k23*(exp(-t*KA)/((lambda1-KA)*(lambda2-KA))+exp(-t*lambda1)/((lambda2-lambda1)*(KA-lambda1))+exp(-t*lambda2)/((lambda1-lambda2)*(KA-lambda2)));

	double Alast0 = Alast[0]*exp(-t*KA);
	if (cmtOff == 0){
	  Alast[1] = A2term1+A2term2;  //Amount in the central compartment
	  Alast[2] = A3term1+A3term2;  //Amount in the peripheral compartment
	  Alast[0] = Alast0 + amt;
	} else {
	  Alast[1] = A2term1+A2term2 + amt;  //Amount in the central compartment
	  Alast[2] = A3term1+A3term2;  //Amount in the peripheral compartment
	  Alast[0] = Alast0;
	}
      } else if (ncmt == 3){
	double k30 = 0;
	double k40 = 0;
#define k23 k12
#define k32 k21
#define k24 k13
#define k42 k31
#define k20 k10
	double E2 = k20+k23+k24;
	double E3 = k32+k30;
	double E4 = k42+k40;

	//#calculate hybrid rate constants
	double a = E2+E3+E4;
	double b = E2*E3+E4*(E2+E3)-k23*k32-k24*k42;
	double c = E2*E3*E4-E4*k23*k32-E3*k24*k42;

	double m = (3.*b - a*a)/3.;
	double n = (2.*a*a*a - 9.*a*b + 27.*c)/27.;
	double Q = (n*n)/4. + (m*m*m)/27.;

	double alpha = sqrt(-1*Q);
	double beta = -1*n/2;
	double gamma = sqrt(beta*beta+alpha*alpha);
	double theta = atan2(alpha,beta);

	double g13 = pow(gamma, 1./3.);
	double c3 = cos(theta/3.);
	double s3 = sqrt(3.)*sin(theta/3.);
	double lambda1 = a/3. + g13*(c3 + s3);
	double lambda2 = a/3. + g13*(c3 - s3);
	double lambda3 = a/3. -(2.*g13*c3);
	
	double B = Alast[2]*k32+Alast[3]*k42;
	double C = E4*Alast[2]*k32+E3*Alast[3]*k42;
	double I = Alast[1]*k23*E4-Alast[2]*k24*k42+Alast[3]*k23*k42;
	double J = Alast[1]*k24*E3+Alast[2]*k24*k32-Alast[3]*k23*k32;

	double A2term1 = Alast[1]*(exp(-t*lambda1)*(E3-lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E3-lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E3-lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	double A2term2 = exp(-t*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2));
	double A2term3 = Alast[0]*KA*(exp(-t*lambda1)*(E3-lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(KA-lambda1))+exp(-t*lambda2)*(E3-lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2)*(KA-lambda2))+exp(-t*lambda3)*(E3-lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)*(KA-lambda3))+exp(-t*KA)*(E3-KA)*(E4-KA)/((lambda1-KA)*(lambda2-KA)*(lambda3-KA)));

	double A3term1 = Alast[2]*(exp(-t*lambda1)*(E2-lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	double A3term2 = exp(-t*lambda1)*(I-Alast[1]*k23*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(Alast[1]*k23*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(Alast[1]*k23*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2));
	double A3term3 = Alast[0]*KA*k23*(exp(-t*lambda1)*(E4-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(KA-lambda1))+exp(-t*lambda2)*(E4-lambda2)/((lambda1-lambda2)*(lambda3-lambda2)*(KA-lambda2))+exp(-t*lambda3)*(E4-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)*(KA-lambda3))+exp(-t*KA)*(E4-KA)/((lambda1-KA)*(lambda2-KA)*(lambda3-KA)));

	double A4term1 = Alast[3]*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	double A4term2 = exp(-t*lambda1)*(J-Alast[1]*k24*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(Alast[1]*k24*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(Alast[1]*k24*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2));
	double A4term3 = Alast[0]*KA*k24*(exp(-t*lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1)*(KA-lambda1))+exp(-t*lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2)*(KA-lambda2))+exp(-t*lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)*(KA-lambda3))+exp(-t*KA)*(E3-KA)/((lambda1-KA)*(lambda2-KA)*(lambda3-KA)));

	double Alast1 = Alast[0]*exp(-t*KA);
	if (cmtOff == 0){
	  Alast[1] = A2term1+A2term2+A2term3; //Amount in the central compartment
	  Alast[2] = A3term1+A3term2+A3term3; //Amount in the first-peripheral compartment
	  Alast[3] = A4term1+A4term2+A4term3;  //Amount in the second-peripheral compartment
	  Alast[0] = Alast1 + amt;
	} else {
	  Alast[1] = A2term1+A2term2+A2term3+amt; //Amount in the central compartment
	  Alast[2] = A3term1+A3term2+A3term3; //Amount in the first-peripheral compartment
	  Alast[3] = A4term1+A4term2+A4term3;  //Amount in the second-peripheral compartment
	  Alast[0] = Alast1;
	}
#undef k23 
#undef k32 
#undef k24
#undef k42 
#undef k20
      }
    } else {
      // Bolus dose without an infusion
      if (ncmt == 1){
	Alast[0] = amt + Alast[0]*exp(-t*k10);
      } else if (ncmt == 2){
	double k20 = 0;
	double E1 = k10+k12;
	double E2 = k21+k20;
	//calculate hybrid rate constants
	double tmp = (k12+k21+k10);
	double sqr = sqrt(tmp*tmp-4*k21*k10);
	double lambda1 = 0.5*(tmp+sqr);
	double lambda2 = 0.5*(tmp-sqr);
	double A1term = (((Alast[0]*E2+Alast[1]*k21)-Alast[0]*lambda1)*exp(-t*lambda1)-((Alast[0]*E2+Alast[1]*k21)-Alast[0]*lambda2)*exp(-t*lambda2))/(lambda2-lambda1);
	double A2term = (((Alast[1]*E1+Alast[0]*k12)-Alast[1]*lambda1)*exp(-t*lambda1)-((Alast[1]*E1+Alast[0]*k12)-Alast[1]*lambda2)*exp(-t*lambda2))/(lambda2-lambda1);
	Alast[0] = A1term+amt; // Amount in the central compartment
	Alast[1] = A2term;
      } else if (ncmt == 3){
	double k20 = 0;
	double k30 = 0;
	double E1 = k10+k12+k13;
	double E2 = k21+k20;
	double E3 = k31+k30;

	//calculate hybrid rate constants
	double a = E1+E2+E3;
	double b = E1*E2+E3*(E1+E2)-k12*k21-k13*k31;
	double c = E1*E2*E3-E3*k12*k21-E2*k13*k31;
	
	double m = (3.*b - a*a)/3.;
	double n = (2.*a*a*a - 9.*a*b + 27.*c)/27.;
	double Q = (n*n)/4. + (m*m*m)/27.;

	double alpha = sqrt(-1*Q);
	double beta = -1*n/2.;
	double gamma = sqrt(beta*beta+alpha*alpha);
	double theta = atan2(alpha,beta);

	double g13 = pow(gamma, 1./3.);
	double c3 = cos(theta/3.);
	double s3 = sqrt(3)*sin(theta/3.);
	double lambda1 = a/3 + g13*(c3 + s3);
	double lambda2 = a/3 + g13*(c3 - s3);
	double lambda3 = a/3 -(2*g13*c3);

	double B = Alast[1]*k21+Alast[2]*k31;
	double C = E3*Alast[1]*k21+E2*Alast[2]*k31;
	double I = Alast[0]*k12*E3-Alast[1]*k13*k31+Alast[2]*k12*k31;
	double J = Alast[0]*k13*E2+Alast[1]*k13*k21-Alast[2]*k12*k21;

	double A1term1 = Alast[0]*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	double A1term2 = exp(-t*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2));


	double A2term1 = Alast[1]*(exp(-t*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	double A2term2 = exp(-t*lambda1)*(I-Alast[0]*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(Alast[0]*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(Alast[0]*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2));

	  
	double A3term1 = Alast[2]*(exp(-t*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)));
	double A3term2 = exp(-t*lambda1)*(J-Alast[0]*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(Alast[0]*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(Alast[0]*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2));

	Alast[0] = amt + A1term1+A1term2;
	Alast[1] = A2term1+A2term2;
	Alast[2] = A3term1+A3term2;
      }
    }
  }
}

void doSort(rx_solving_options_ind *ind);
void calcMtime(int solveid, double *mtime);
void setLinCmt(int linCmt, double lag, double lag2, double f, double f2,
	       double rate, double dur);

void updateDoseIdx(rx_solving_options_ind *ind);
// Advan-style linCmt solutions
double linCmtAB(rx_solve *rx, unsigned int id, double t, int linCmt,
		int i_cmt, int trans, 
		double p1, double v1,
		double p2, double p3,
		double p4, double p5,
		double d_ka, double d_tlag, double d_tlag2, double d_F, double d_F2,
		// Rate and dur can only apply to central compartment even w/ oral dosing
		// Therefore, only 1 model rate is possible with RxODE
		double d_rate, double d_dur){
  rx_solving_options *op = rx->op;
  int oral0;
  oral0 = (d_ka > 0) ? 1 : 0;
  unsigned int ncmt = 1;
  rx_solving_options_ind *ind = &(rx->subjects[id]);
  // don't need to adjust based on tlag t is the most conservative.
  // When tadr - tlag < 0 ignore the dose.
  double rx_k=0, rx_v=0;
  double rx_k12=0;
  double rx_k21=0;
  double rx_k13=0;
  double rx_k31=0;
  /* k10 <- d$CL[i]/d$V1[i] */
  /* k12 <- d$Q[i]/d$V1[i] */
  /* k21 <- d$Q[i]/d$V2[i] */
  /* k20 <- 0 */
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
      // FIXME -- add warning
      return linCmtAA(rx, id, t, linCmt, i_cmt, trans, p1, v1,
		      p2, p3, p4, p5, d_ka, d_tlag, d_tlag2,  d_F,  d_F2, d_rate, d_dur);
      /* REprintf("V, alpha, beta, k21 are not supported with ADVAN routines"); */
      /* return NA_REAL; */
      break;
    case 10:
      // FIXME -- add warning
      return linCmtAA(rx, id, t, linCmt, i_cmt, trans, p1, v1,
		      p2, p3, p4, p5, d_ka, d_tlag, d_tlag2,  d_F,  d_F2, d_rate, d_dur);
      break;
    default:
      REprintf("invalid trans (3 cmt trans %d).\n", trans);
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
    case 11: // A2 V, alpha, beta, k21
      // FIXME -- add warning
      return linCmtAA(rx, id, t, linCmt, i_cmt, trans, p1, v1,
		      p2, p3, p4, p5, d_ka, d_tlag, d_tlag2,  d_F,  d_F2, d_rate, d_dur);
      /* REprintf("V, alpha, beta, k21 are not supported with ADVAN routines"); */
      /* return NA_REAL; */
      break;
    case 10: // A, alpha, B, beta
      // FIXME -- add warning
      return linCmtAA(rx, id, t, linCmt, i_cmt, trans, p1, v1,
		      p2, p3, p4, p5, d_ka, d_tlag, d_tlag2,  d_F,  d_F2, d_rate, d_dur);
      break;
    default:
      REprintf("invalid trans (2 cmt trans %d).\n", trans);
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
      REprintf("invalid trans (1 cmt trans %d).\n", trans);
      return NA_REAL;
    }
  } else {
    return 0.0;
  }
  if (ind->linCmtAdvanSetup == 0){
    ind->linCmtAdvan = Calloc((ncmt+oral0+2)*ind->n_all_times, double);
    for (int ii = ind->n_all_times; ii--; ) ind->linCmtAdvan[ii] = NA_REAL;
    ind->linCmtAdvanSetup=1;
  } else if (t < ind->all_times[0] && ind->linCmtAdvanSetup!=2){
    ind->linCmtAdvan[0] = NA_REAL;
    ind->linCmtAdvanSetup=3;
  } else if (t == ind->all_times[0] && ind->linCmtAdvanSetup!=2){
    for (int ii = ind->n_all_times; ii--; ) ind->linCmtAdvan[ii] = NA_REAL;
    ind->linCmtAdvanSetup=2;
  } else {
    ind->linCmtAdvanSetup=3;
  }
  double F[2] = {d_F, d_F2};
  double Alast[6] = {0};
  double Asave[6] = {0};
  double tlast = 0;
  double xout=0;
  int cmtOff=0;
#define doDose realizeBolus(Alast, tlast, ncmt, oral0, cmtOff, ind->dose[ind->ixds]*F[cmtOff], Alast[ncmt+oral0], xout, d_ka, rx_k, rx_k12, rx_k21, rx_k13, rx_k31)
#define doObs realizeBolus(Alast, tlast, ncmt, oral0, cmtOff, 0.0, Alast[ncmt+oral0], xout, d_ka, rx_k, rx_k12, rx_k21, rx_k13, rx_k31)
  int needSort = 1;
  for (int i = 0; i < ind->n_all_times; i++){
    if (ISNA(ind->linCmtAdvan[i])){
      if (needSort){
	setLinCmt(linCmt, d_tlag, d_tlag2, d_F, d_F2, d_rate, d_dur);
	if (rx->nMtime) calcMtime(id, ind->mtime);
	// FIXME do partial sort when i > 0
	if (rx->needSort) doSort(ind);
	needSort = 0;
      }
      ind->idx = i;
      xout = getTime(ind->ix[i], ind);
      // Value hasn't been calculated yet.
      int evid, wh, cmt, wh100, whI, wh0;
      evid = ind->evid[ind->ix[i]];
      if (isObs(evid)){
	cmtOff=0;
	doObs;
	tlast = xout;
      } else {
	if (evid == 3){ // Reset
	  for (int ii = ncmt + oral0+1; ii--; ){ 
	    Alast[ii] = 0;
	  }
	  tlast = xout;
	  continue;
	}
	getWh(evid, &wh, &cmt, &wh100, &whI, &wh0);
	cmtOff = linCmt-cmt;
	if (cmtOff < 0) {ind->ixds++;continue;} // Dose to comparment outside of system.
	if (oral0 && cmtOff > 1) {ind->ixds++;continue;}
	if (!oral0 && cmtOff != 0) {ind->ixds++;continue;}
	if (ind->ix[ind->idx] != ind->idose[ind->ixds]){
	  int foundIt=0;
	  // FIXME use bisection
	  for (int j = 0; j < ind->ndoses; j++){
	    if (ind->idose[j] == ind->ix[ind->idx]){
	      ind->ixds = j;
	      foundIt=1;
	      break;
	    }
	  }
	  if (foundIt==0) error("Corrupted event table.");
	}
	xout = ind->all_times[ind->idose[ind->ixds]];
	if (t < xout) break;
	if (wh0 == 30){
	  // Reset before dosing
	  error("Cannot turn off a compartment with a linear solved system");
	}
	switch(whI){
	case 9: // modeled rate.
	case 8: // modeled duration.
	  // Rate already calculated and saved in the next dose record
	  doObs;
	  tlast=xout;
	  Alast[ncmt+oral0]  -= ind->dose[ind->ixds+1];
	  /* if (wh0 == 20 && AMT(id, cmt, dose[ind->ixds], xout) != dose[ind->ixds]){ */
	  /*   if (!(ind->err & 1048576)){ */
	  /*     ind->err += 1048576; */
	  /*   } */
	  /*   return 0; */
	  /* } */
	  /* error("SS=2 & Modeled F does not work"); */
	  break;
	case 7: // End modeled rate
	case 6: // end modeled duration
	  doObs;
	  tlast=xout;
	  Alast[ncmt+oral0] += ind->dose[ind->ixds];
	  break;
	case 2:
	  // In this case bio-availability changes the rate, but the duration remains constant.
	  // rate = amt/dur
	  doObs;
	  tlast=xout;
	  Alast[ncmt+oral0] += ind->dose[ind->ixds]*F[cmtOff];
	  break;
	case 1:
	  doObs;
	  tlast=xout;
	  Alast[ncmt+oral0] += ind->dose[ind->ixds];
	}
	double rate  = 0;
	if (whI != 0) rate = -ind->dose[ind->ixds+1];

	if ((whI == 0 || rate > 0) && (wh0 == 10 || wh0 == 20)) {
	  // steady state event SS=2
	  // steady state event SS=1
	  double tau = ind->ii[ind->ixds];
	  double tmp = xout;
	  double lastSum = 0.0, curSum = 0.0, tinf = 0;
	  // Turn off RATE and reset compartments
	  if (wh0 == 20){
	    xout=tmp;
	    doObs;
	    for (int ii = ncmt + oral0+2; ii--; ){
	      Asave[ii] = Alast[ii];
	      Alast[ii] = 0;
	    }
	  } else {
	    for (int ii = ncmt + oral0+2; ii--; ){
	      Alast[ii] = 0;
	    }
	  }
	  if (whI == 1 || whI == 2){
	    if (rate <= 0) continue;
	    tinf = _getDur(ind->ixds, ind, 0);
	    if (whI == 2)rate*=F[cmtOff];
	    else tinf*=F[cmtOff];
	  }
	  for (int j = 0; j < op->maxSS; j++){
	    xout = 0;
	    tlast = 0;
	    if (whI == 0){
	      doDose;
	      xout = tau;
	      doObs;
	    } else {
	      doObs;
	      Alast[ncmt+oral0] += rate;
	      xout  = tinf;
	      doObs;
	      tlast = tinf;
	      Alast[ncmt+oral0] -= rate;
	      xout = tau;
	      doObs;
	    }
	    if (j == op->minSS -1){
	      lastSum =0.0;
	      for (int k = ncmt + oral0 ; k--;) lastSum += Alast[k];
	    } else if (j >= op->minSS){
	      curSum = 0.0;
	      for (int k = ncmt + oral0; k--;) curSum += Alast[k];
	      if (fabs(curSum-lastSum) < op->rtolSS*fabs(curSum) + op->atolSS){
		j = op->maxSS+1;
	      }
	      lastSum=curSum;
	    }
	  }
	  xout = 0;
	  tlast =0;
	  if (whI == 0){
	    doDose;
	  } else {
	    Alast[ncmt+oral0] += rate;
	  }
	  xout = tmp;
	  tlast = tmp;
	  if (wh0 == 20){
	    for (int ii = ncmt + oral0; ii--; ){
	      Alast[ii] += Asave[ii];
	    }
	  }
	} else {
	  if (whI == 0){
	    doDose;
	    tlast = xout;
	  }
	}
      }
      if (xout == t){
	for (int ii = ncmt + oral0+1; ii--; ){
	  ind->linCmtAdvan[ii*ind->n_all_times+i] = Alast[ii];
	}
	ind->linCmtAdvan[(ncmt+oral0+1)*ind->n_all_times+i] = tlast;
	/* setLinCmt(-100, 0, 0, 1, 1, 0, 0); */
	/* return Alast[oral0]/rx_v; */
      }
      if (i+1 < ind->n_all_times){
	xout = getTime(ind->ix[i+1], ind);
	if (xout > t){
	  xout = t;
	  cmtOff=0;
	  doObs;
	  tlast = xout;
	  setLinCmt(-100, 0, 0, 1, 1, 0, 0);
	  return Alast[oral0]/rx_v;
	}
      } else if (i >= ind->n_all_times){
	xout = t;
	cmtOff=0;
	doObs;
	tlast = xout;
	setLinCmt(-100, 0, 0, 1, 1, 0, 0);
	return Alast[oral0]/rx_v;
      }
    } else {
      // Set last values
      for (int ii = ncmt + oral0+1; ii--; ){
	Alast[ii] = ind->linCmtAdvan[ii*ind->n_all_times+i];
      }
      tlast = ind->linCmtAdvan[(ncmt+oral0+1)*ind->n_all_times+i];
    }
  }
  if (t > tlast){
    xout = t;
    cmtOff=0;
    doObs;
  }
  setLinCmt(-100, 0, 0, 1, 1, 0, 0);
  return (Alast[oral0]/rx_v);
}

double linCmtA(rx_solve *rx, unsigned int id, double t, int linCmt,
	       int i_cmt, int trans, 
	       double p1, double v1,
	       double p2, double p3,
	       double p4, double p5,
	       double d_ka, double d_tlag, double d_tlag2, double d_F, double d_F2,
	       // Rate and dur can only apply to central compartment even w/ oral dosing
	       // Therefore, only 1 model rate is possible with RxODE
	       double d_rate, double d_dur){
  rx_solving_options *op = rx->op;
  if (op->advanLinCmt){
    return linCmtAB(rx, id, t, linCmt, i_cmt, trans, p1, v1, p2, p3, p4, p5, d_ka, d_tlag, d_tlag2, d_F, d_F2,
  		    d_rate, d_dur);
  } else {
    return linCmtAA(rx, id, t, linCmt, i_cmt, trans, p1, v1, p2, p3, p4, p5, d_ka, d_tlag, d_tlag2, d_F, d_F2,
		  d_rate, d_dur);
  }
}
#undef V
