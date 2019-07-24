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

static inline int _locateDoseIndex(const double obs_time,  rx_solving_options_ind *ind){
  // Uses bisection for slightly faster lookup of dose index.
  int i, j, ij, idose;
  i = 0;
  j = ind->ndoses - 1;
  /* Rprintf("_locateDoseIndex ind->ndoses: %d; ind->id: %d\n",ind->ndoses, ind->id); */
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
extern double getTime(int idx, rx_solving_options_ind *ind);

double solveLinB(rx_solve *rx, unsigned int id, double t, int linCmt,
		 double d_A, double d_A2, double d_alpha,
		 double d_B, double d_B2, double d_beta,
		 double d_C, double d_C2, double d_gamma,
		 double d_ka,
		 // Lag and F can apply to 2 compartments, depot and/or central
		 double d_tlag, double d_tlag2, double d_F, double d_F2,
		 // Rate and dur can only apply to central compartment even w/ oral dosing
		 // Therefore, only 1 model rate is possible with RxODE
		 double d_rate, double d_dur){
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
  int ndoses = ind->ndoses;
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
	p = l+1;
	while (p < ndoses && ind->dose[p] != -dose){
	  p++;
	}
	if (ind->dose[p] != -dose){
	  error("Could not find an end to the infusion.  Check the event table.");
	}
	tinf  = ind->all_times[ind->idose[p]] - ind->all_times[ind->idose[l]];
	tau = ind->ii[l];
	if (op->linLog){
	  logRate = log(dose);
	} else {
	  rate  = dose;
	}
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
	tinf  = ind->all_times[ind->idose[l]] - ind->all_times[ind->idose[p]];
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
	tinf = tinf*F;
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
double linCmtA(rx_solve *rx, unsigned int id, double t, int linCmt,
	       int cmt, int trans, 
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
    if (cmt == 1){
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
    } else if (cmt == 2){
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
      double tmp = rx_k12+rx_k21+rx_k;
      d_beta = 0.5 * (tmp - sqrt(tmp * tmp - 4.0 * rx_k21 * rx_k));
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
    } else if (cmt == 3){
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
  return solveLinB(rx, id, t, linCmt, d_A, d_A2, d_alpha, d_B, d_B2, d_beta,
		   d_C, d_C2, d_gamma, d_ka, d_tlag, d_tlag2, d_F, d_F2,
		   d_rate, d_dur);
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
    /* Rprintf("NA->%f for id=%d\n", ret, ind->id); */
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
    if (op->neq > 0){
      // Update all covariate parameters
      int k;
      int ncov = op->ncov;
      if (op->do_par_cov){
	for (k = ncov; k--;){
	  if (op->par_cov[k]){
	    double *y = ind->cov_ptr + ind->n_all_times*k;
	    ind->par_ptr[op->par_cov[k]-1] = getValue(idx, y, ind);
	  }
	}
      }
    }
  } else {
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
	    double *y = ind->cov_ptr + ind->n_all_times*k;
	    if (idx > 0 && idx < ind->n_all_times && t == all_times[idx]){
	      par_ptr[op->par_cov[k]-1] = getValue(idx, y, ind);
	    } else {
	      // Use the same methodology as approxfun.
	      ind->ylow = getValue(0, y, ind);/* cov_ptr[ind->n_all_times*k]; */
	      ind->yhigh = getValue(ind->n_all_times-1, y, ind);/* cov_ptr[ind->n_all_times*k+ind->n_all_times-1]; */
	      par_ptr[op->par_cov[k]-1] = rx_approxP(t, y, ind->n_all_times, op, ind);
	    }
	  }
	}
      }
    }
  }
}
#undef V
