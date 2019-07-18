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
		 double d_A, double d_A2, double d_alpha, //
		 double d_B, double d_B2, double d_beta,
		 double d_C, double d_C2, double d_gamma,
		 double d_ka,
		 // Lag and F can apply to 2 compartments, depot and/or central
		 double d_tlag, double d_tlag2, double d_F, double d_F2,
		 // Rate and dur can only apply to central compartment even w/ oral dosing
		 // Therefore, only 1 model rate is possible with RxODE
		 double d_rate, double d_dur, int dPar){
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
  int in2=0;
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
    in2=0;
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
      in2=1;
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
      if (wh0 == 10 || wh0 == 20){ // steady state
	if (tinf >= tau){
	  error("Infusion time greater then inter-dose interval, ss cannot be calculated.");
	} 
	if (thisT < tinf){ // during infusion
	  if (dPar==0){
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
	    break;
	  } else if (dPar == 1){
	    // dA
	    if (in2==0){
	      cur += rate*alpha1*((1-exp(-alpha*thisT))+
				  exp(-alpha*tau)*(1-exp(-alpha*tinf))*
				  exp(-alpha*(thisT-tinf))/(1-exp(-alpha*tau)));
	    }
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 2){
	    // dA2
	    if (in2!=0){
	      cur += rate*alpha1*((1-exp(-alpha*thisT))+
				  exp(-alpha*tau)*(1-exp(-alpha*tinf))*
				  exp(-alpha*(thisT-tinf))/(1-exp(-alpha*tau)));
	    }
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 3){
	    // dalpha
	    cur += -rate*A*(1 + exp(-alpha*(thisT - tinf) - tau*alpha)*(1 - exp(-tinf*alpha))/(1 - exp(-tau*alpha)) - exp(-thisT*alpha))/(alpha*alpha) + rate*(exp(-thisT*alpha)*thisT + exp(-alpha*(thisT - tinf) - tau*alpha - tinf*alpha)*tinf/(1 - exp(-tau*alpha)) - exp(-alpha*(thisT - tinf) - 2*tau*alpha)*tau*(1 - exp(-tinf*alpha))/R_pow_di(1 - exp(-tau*alpha),2) + exp(-alpha*(thisT - tinf) - tau*alpha)*(-tau - thisT + tinf)*(1 - exp(-tinf*alpha))/(1 - exp(-tau*alpha)))/alpha;
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 4 && ncmt >= 2) {
	    // dB1
	    if (in2 == 0 && ncmt >= 2){
	      cur += rate*beta1*((1-exp(-beta*thisT))+
				 exp(-beta*tau)*(1-exp(-beta*tinf))*
				 exp(-beta*(thisT-tinf))/
				 (1-exp(-beta*tau)));
	    }
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 5 && ncmt >= 2){
	    // dB2
	    if (in2 != 0){
	      cur += rate*beta1*((1-exp(-beta*thisT))+
				 exp(-beta*tau)*(1-exp(-beta*tinf))*
				 exp(-beta*(thisT-tinf))/
				 (1-exp(-beta*tau)));
	    }
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 6 && ncmt >= 2){
	    // dBeta
	    cur += -B*rate*(1 + exp(-beta*(thisT - tinf) - tau*beta)*(1 - exp(-tinf*beta))/(1 - exp(-tau*beta)) - exp(-beta*thisT))/(beta*beta) + B*rate*(exp(-beta*thisT)*thisT + exp(-beta*(thisT - tinf) - tau*beta - tinf*beta)*tinf/(1 - exp(-tau*beta)) - exp(-beta*(thisT - tinf) - 2*tau*beta)*tau*(1 - exp(-tinf*beta))/R_pow_di(1 - exp(-tau*beta), 2) + exp(-beta*(thisT - tinf) - tau*beta)*(1 - exp(-tinf*beta))*(-tau - thisT + tinf)/(1 - exp(-tau*beta)))/beta;
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 7 && ncmt >= 3){
	    // dC1
	    if (in2 == 0){
	      cur += rate*gamma1*((1-exp(-gamma*thisT))+
				    exp(-gamma*tau)*(1-exp(-gamma*tinf))*exp(-gamma*(thisT-tinf))/
				    (1-exp(-gamma*tau)));
	    }
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 8 && ncmt >= 3){
	    // dC2
	    if (in2 != 0){
	      cur += rate*gamma1*((1-exp(-gamma*thisT))+
				    exp(-gamma*tau)*(1-exp(-gamma*tinf))*exp(-gamma*(thisT-tinf))/
				    (1-exp(-gamma*tau)));
	    }
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 9 && ncmt >= 3){
	    // dGamma
	    cur += -C*rate*(1 + exp(-gamma*(thisT - tinf) - tau*gamma)*(1 - exp(-tinf*gamma))/(1 - exp(-tau*gamma)) - exp(-thisT*gamma))/(gamma*gamma) + C*rate*(exp(-thisT*gamma)*thisT + exp(-gamma*(thisT - tinf) - tau*gamma - tinf*gamma)*tinf/(1 - exp(-tau*gamma)) - exp(-gamma*(thisT - tinf) - 2*tau*gamma)*tau*(1 - exp(-tinf*gamma))/R_pow_di(1 - exp(-tau*gamma),2) + exp(-gamma*(thisT - tinf) - tau*gamma)*(-tau - thisT + tinf)*(1 - exp(-tinf*gamma))/(1 - exp(-tau*gamma)))/gamma;
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 10){
	    // dKa
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 11){
	    // dTlag1
	    if (in2==0){
	      cur += -A*rate*alpha1*(exp(-thisT*alpha)*alpha - exp(-alpha*(thisT - tinf) - tau*alpha)*alpha*(1 - exp(-tinf*alpha))/(1 - exp(-tau*alpha)));
	      if (ncmt >= 2){
		cur += -B*rate*beta1*(exp(-beta*thisT)*beta - exp(-beta*(thisT - tinf) - tau*beta)*beta*(1 - exp(-tinf*beta))/(1 - exp(-tau*beta)));
		if (ncmt >= 3){
		  cur += -C*rate*gamma1*(exp(-thisT*gamma)*gamma - exp(-gamma*(thisT - tinf) - tau*gamma)*gamma*(1 - exp(-tinf*gamma))/(1 - exp(-tau*gamma)));
		}
	      }
	    }
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 12){
	    // dTlag2
	    if (in2!=0){
	      cur += -A*rate*alpha1*(exp(-thisT*alpha)*alpha - exp(-alpha*(thisT - tinf) - tau*alpha)*alpha*(1 - exp(-tinf*alpha))/(1 - exp(-tau*alpha)));
	      if (ncmt >= 2){
		cur += -B*rate*beta1*(exp(-beta*thisT)*beta - exp(-beta*(thisT - tinf) - tau*beta)*beta*(1 - exp(-tinf*beta))/(1 - exp(-tau*beta)));
		if (ncmt >= 3){
		  cur += -C*rate*gamma1*(exp(-thisT*gamma)*gamma - exp(-gamma*(thisT - tinf) - tau*gamma)*gamma*(1 - exp(-tinf*gamma))/(1 - exp(-tau*gamma)));
		}
	      }
	    }
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 13){
	    // dF1
	    if (in2==0){
	      if (whI == 1){
		// Duration where tinfM=tinf*F
		cur += A*rate*alpha1*(exp(-alpha*(thisT - tinf) - tau*alpha - tinf*alpha)*alpha/(1 - exp(-tau*alpha)) + exp(-alpha*(thisT - tinf) - tau*alpha)*alpha*(1 - exp(-tinf*alpha))/(1 - exp(-tau*alpha)))*tinf/F;
		if (ncmt >= 2){
		  cur += B*rate*beta1*(exp(-beta*(thisT - tinf) - tau*beta - tinf*beta)*beta/(1 - exp(-tau*beta)) + exp(-beta*(thisT - tinf) - tau*beta)*beta*(1 - exp(-tinf*beta))/(1 - exp(-tau*beta)))*tinf/F;
		  if (ncmt >= 3){
		    cur += C*rate*(exp(-gamma*(thisT - tinf) - tau*gamma - tinf*gamma)*gamma/(1 - exp(-tau*gamma)) + exp(-gamma*(thisT - tinf) - tau*gamma)*gamma*(1 - exp(-tinf*gamma))/(1 - exp(-tau*gamma)))*gamma1*tinf/F;
		  }
		}
	      } else {
		// Rate where rate=F*rate0
		cur += A*alpha1*((1-exp(-alpha*thisT))+
				 exp(-alpha*tau)*(1-exp(-alpha*tinf))*
				 exp(-alpha*(thisT-tinf))/(1-exp(-alpha*tau)))*rate/F;
		if (ncmt >= 2){
		  cur += B*beta1*((1-exp(-beta*thisT))+
				  exp(-beta*tau)*(1-exp(-beta*tinf))*
				  exp(-beta*(thisT-tinf))/
				  (1-exp(-beta*tau)))*rate/F;
		  if (ncmt >= 3){
		    cur += C*gamma1*((1-exp(-gamma*thisT))+
				     exp(-gamma*tau)*(1-exp(-gamma*tinf))*exp(-gamma*(thisT-tinf))/
				     (1-exp(-gamma*tau)))*rate/F;
		  }
		}
	      }
	    }
	  } else if (dPar == 14){
	    // dF2
	    if (in2!=0){
	      if (whI == 1){
		// Duration where tinfM=tinf*F
		cur += A*rate*alpha1*(exp(-alpha*(thisT - tinf) - tau*alpha - tinf*alpha)*alpha/(1 - exp(-tau*alpha)) + exp(-alpha*(thisT - tinf) - tau*alpha)*alpha*(1 - exp(-tinf*alpha))/(1 - exp(-tau*alpha)))*tinf/F;
		if (ncmt >= 2){
		  cur += B*rate*beta1*(exp(-beta*(thisT - tinf) - tau*beta - tinf*beta)*beta/(1 - exp(-tau*beta)) + exp(-beta*(thisT - tinf) - tau*beta)*beta*(1 - exp(-tinf*beta))/(1 - exp(-tau*beta)))*tinf/F;
		  if (ncmt >= 3){
		    cur += C*rate*(exp(-gamma*(thisT - tinf) - tau*gamma - tinf*gamma)*gamma/(1 - exp(-tau*gamma)) + exp(-gamma*(thisT - tinf) - tau*gamma)*gamma*(1 - exp(-tinf*gamma))/(1 - exp(-tau*gamma)))*gamma1*tinf/F;
		  }
		}
	      } else {
		// Rate where rate=F*rate0
		cur += A*alpha1*((1-exp(-alpha*thisT))+
				 exp(-alpha*tau)*(1-exp(-alpha*tinf))*
				 exp(-alpha*(thisT-tinf))/(1-exp(-alpha*tau)))*rate/F;
		if (ncmt >= 2){
		  cur += B*beta1*((1-exp(-beta*thisT))+
				  exp(-beta*tau)*(1-exp(-beta*tinf))*
				  exp(-beta*(thisT-tinf))/
				  (1-exp(-beta*tau)))*rate/F;
		  if (ncmt >= 3){
		    cur += C*gamma1*((1-exp(-gamma*thisT))+
				     exp(-gamma*tau)*(1-exp(-gamma*tinf))*exp(-gamma*(thisT-tinf))/
				     (1-exp(-gamma*tau)))*rate/F;
		  }
		}
	      }
	    }
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 15){
	    // dRate -- Fixed infusion not modeled
	    if (wh0 == 10) return (ret+cur);
	  } else if (dPar == 16){
	    // dDuration -- Fixed infusion not modeled
	    if (wh0 == 10) return (ret+cur);
	  }
	} else { // after infusion
	  if (dPar==0){
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
	  } else if (dPar == 1){
	    // dA
	    if (in2==0){
	      cur += rate*alpha1*((1-exp(-alpha*tinf))*exp(-alpha*(thisT-tinf))/(1-exp(-alpha*tau)));
	    }
	  } else if (dPar == 2){
	    // dA2
	    if (in2 != 0){
	      cur += rate*alpha1*((1-exp(-alpha*tinf))*exp(-alpha*(thisT-tinf))/(1-exp(-alpha*tau)));
	    }
	  } else if (dPar == 3){
	    // dalpha
	    cur += -exp(-alpha*(thisT - tinf))*rate*(1 - exp(-tinf*alpha))/(alpha*alpha*(1 - exp(-tau*alpha))) + exp(-alpha*(thisT - tinf) - tinf*alpha)*tinf*rate/(alpha*(1 - exp(-tau*alpha))) - exp(-alpha*(thisT - tinf))*rate*(thisT - tinf)*(1 - exp(-tinf*alpha))/(alpha*(1 - exp(-tau*alpha))) - exp(-alpha*(thisT - tinf) - tau*alpha)*tau*rate*(1 - exp(-tinf*alpha))/(alpha*R_pow_di(1 - exp(-tau*alpha),2));
	  } else if (dPar == 4 && ncmt >= 2){
	    // dB1
	    if (in2==0){
	      cur += rate*beta1*((1-exp(-beta*tinf))*exp(-beta*(thisT-tinf))/(1-exp(-beta*tau)));
	    }
	  } else if (dPar == 5 && ncmt >= 2){
	    // dB2
	    if (in2!=0){
	      cur += rate*beta1*((1-exp(-beta*tinf))*exp(-beta*(thisT-tinf))/(1-exp(-beta*tau)));
	    }
	  } else if (dPar == 6 && ncmt >= 2){
	    // dBeta
	    cur += -B*exp(-beta*(thisT - tinf))*rate*(1 - exp(-tinf*beta))/(beta*beta*(1 - exp(-tau*beta))) + B*exp(-beta*(thisT - tinf) - tinf*beta)*tinf*rate/(beta*(1 - exp(-tau*beta))) - B*exp(-beta*(thisT - tinf))*rate*(1 - exp(-tinf*beta))*(thisT - tinf)/(beta*(1 - exp(-tau*beta))) - B*exp(-beta*(thisT - tinf) - tau*beta)*tau*rate*(1 - exp(-tinf*beta))/(beta*R_pow_di(1 - exp(-tau*beta),2));
	  } else if (dPar == 7 && ncmt >= 3){
	    // dC1
	    if (in2==0){
	      cur += rate*gamma1*((1-exp(-gamma*tinf))*exp(-gamma*(thisT-tinf))/(1-exp(-gamma*tau)));
	    }
	  } else if (dPar == 8 && ncmt >= 3){
	    // dC2
	    if (in2!=0){
	      cur += rate*gamma1*((1-exp(-gamma*tinf))*exp(-gamma*(thisT-tinf))/(1-exp(-gamma*tau)));
	    }
	  } else if (dPar == 9 && ncmt >= 3){
	    // dGamma
	    cur += -C*exp(-gamma*(thisT - tinf))*rate*(1 - exp(-tinf*gamma))/(gamma*gamma*(1 - exp(-tau*gamma))) + C*exp(-gamma*(thisT - tinf) - tinf*gamma)*tinf*rate/(gamma*(1 - exp(-tau*gamma))) - C*exp(-gamma*(thisT - tinf))*rate*(thisT - tinf)*(1 - exp(-tinf*gamma))/(gamma*(1 - exp(-tau*gamma))) - C*exp(-gamma*(thisT - tinf) - tau*gamma)*tau*rate*(1 - exp(-tinf*gamma))/(gamma*R_pow_di(1 - exp(-tau*gamma),2));
	  } else if (dPar == 10){
	    // dKa
	  } else if (dPar == 11){
	    // dTlag1
	    if (in2 == 0){
	      cur += A*exp(-alpha*(thisT - tinf))*rate*alpha*alpha1*(1 - exp(-tinf*alpha))/(1 - exp(-tau*alpha));
	      if (ncmt >= 2){
		cur += B*exp(-beta*(thisT - tinf))*rate*beta*beta1*(1 - exp(-tinf*beta))/(1 - exp(-tau*beta));
		if (ncmt >= 3){
		  cur += C*exp(-gamma*(thisT - tinf))*rate*gamma*gamma1*(1 - exp(-tinf*gamma))/(1 - exp(-tau*gamma));
		}
	      }
	    }
	  } else if (dPar == 12){
	    // dTlag2
	    if (in2 != 0){
	      cur += A*exp(-alpha*(thisT - tinf))*rate*alpha*alpha1*(1 - exp(-tinf*alpha))/(1 - exp(-tau*alpha));
	      if (ncmt >= 2){
		cur += B*exp(-beta*(thisT - tinf))*rate*beta*beta1*(1 - exp(-tinf*beta))/(1 - exp(-tau*beta));
		if (ncmt >= 3){
		  cur += C*exp(-gamma*(thisT - tinf))*rate*gamma*gamma1*(1 - exp(-tinf*gamma))/(1 - exp(-tau*gamma));
		}
	      }
	    }
	  } else if (dPar == 13){
	    // dF1
	    if (in2==0){
	      if (whI == 1){
		// Duration where tinfM=tinf*F
		cur += A*exp(-alpha*(thisT - tinf) - tinf*alpha)*rate*alpha*alpha1/(1 - exp(-tau*alpha)) + A*exp(-alpha*(thisT - tinf))*rate*alpha*alpha1*(1 - exp(-tinf*alpha))/(1 - exp(-tau*alpha))*tinf/F;
		if (ncmt >= 2){
		  cur += B*exp(-beta*(thisT - tinf) - tinf*beta)*rate*beta*beta1/(1 - exp(-tau*beta)) + B*exp(-beta*(thisT - tinf))*rate*beta*beta1*(1 - exp(-tinf*beta))/(1 - exp(-tau*beta))*tinf/F;
		  if (ncmt >= 3){
		    cur += C*exp(-gamma*(thisT - tinf) - tinf*gamma)*rate*gamma*gamma1/(1 - exp(-tau*gamma)) + C*exp(-gamma*(thisT - tinf))*rate*gamma*gamma1*(1 - exp(-tinf*gamma))/(1 - exp(-tau*gamma))*tinf/F;
		  }
		}
	      } else {
		// Rate where rate=F*rate0
		cur += rate*A*alpha1*((1-exp(-alpha*tinf))*exp(-alpha*(thisT-tinf))/(1-exp(-alpha*tau)))/F;
		if (ncmt >= 2){
		  cur += rate*B*beta1*((1-exp(-beta*tinf))*exp(-beta*(thisT-tinf))/(1-exp(-beta*tau)))/F;
		  if (ncmt >= 3){
		    cur += rate*C*gamma1*((1-exp(-gamma*tinf))*exp(-gamma*(thisT-tinf))/(1-exp(-gamma*tau)))/F;
		  }
		}
	      }
	    }
	  } else if (dPar == 14){
	    // dF2
	    if (in2!=0){
	      if (whI == 1){
		// Duration where tinfM=tinf*F
		cur += A*exp(-alpha*(thisT - tinf) - tinf*alpha)*rate*alpha*alpha1/(1 - exp(-tau*alpha)) + A*exp(-alpha*(thisT - tinf))*rate*alpha*alpha1*(1 - exp(-tinf*alpha))/(1 - exp(-tau*alpha))*tinf/F;
		if (ncmt >= 2){
		  cur += B*exp(-beta*(thisT - tinf) - tinf*beta)*rate*beta*beta1/(1 - exp(-tau*beta)) + B*exp(-beta*(thisT - tinf))*rate*beta*beta1*(1 - exp(-tinf*beta))/(1 - exp(-tau*beta))*tinf/F;
		  if (ncmt >= 3){
		    cur += C*exp(-gamma*(thisT - tinf) - tinf*gamma)*rate*gamma*gamma1/(1 - exp(-tau*gamma)) + C*exp(-gamma*(thisT - tinf))*rate*gamma*gamma1*(1 - exp(-tinf*gamma))/(1 - exp(-tau*gamma))*tinf/F;
		  }
		}
	      } else {
		// Rate where rate=F*rate0
		cur += rate*A*alpha1*((1-exp(-alpha*tinf))*exp(-alpha*(thisT-tinf))/(1-exp(-alpha*tau)))/F;
		if (ncmt >= 2){
		  cur += rate*B*beta1*((1-exp(-beta*tinf))*exp(-beta*(thisT-tinf))/(1-exp(-beta*tau)))/F;
		  if (ncmt >= 3){
		    cur += rate*C*gamma1*((1-exp(-gamma*tinf))*exp(-gamma*(thisT-tinf))/(1-exp(-gamma*tau)))/F;
		  }
		}
	      }
	    }
	  } else if (dPar == 15){
	    // dRate -- Fixed infusion not modeled
	  } else if (dPar == 16){
	    // dDuration -- Fixed infusion not modeled
	  }
	}
	if (wh0 == 10) return (ret+cur);
      } else {
	// Non steady state
	t1  = ((thisT < tinf) ? thisT : tinf);        //during infusion
	t2  = ((thisT > tinf) ? thisT - tinf : 0.0);  // after infusion	  
	if (dPar==0){
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
	} else if (dPar == 2){
	  // dA
	  if (in2 == 0){
	    cur += rate*alpha1*(1.0-exp(-alpha*t1))*exp(-alpha*t2);
	  }
	} else if (dPar == 3){
	  // dA2
	  if (in2 != 0){
	    cur += rate*alpha1*(1.0-exp(-alpha*t1))*exp(-alpha*t2);
	  }
	} else if (dPar == 4){
	  // dAlpha
	  cur += -A*exp(-t2*alpha)*rate*(1.0 - exp(-t1*alpha))/(alpha*alpha) + A*exp(-t1*alpha - t2*alpha)*t1*rate/alpha - A*exp(-t2*alpha)*t2*rate*(1.0 - exp(-t1*alpha))/alpha;
	} else if (dPar == 5 && ncmt >= 2){
	  // dB1
	  if (in2 == 0){
	    cur +=  rate*beta1*(1.0-exp(-beta*t1))*exp(-beta*t2);
	  }
	} else if (dPar == 5 && ncmt >= 2){
	  // dB2
	  if (in2 != 0){
	    cur +=  rate*beta1*(1.0-exp(-beta*t1))*exp(-beta*t2);
	  }
	} else if (dPar == 6 && ncmt >= 2){
	  // dBeta
	  cur += -B*exp(-t2*beta)*rate*(1.0 - exp(-t1*beta))/(beta*beta) + B*exp(-t1*beta - t2*beta)*t1*rate/beta - B*exp(-t2*beta)*t2*rate*(1.0 - exp(-t1*beta))/beta;
	} else if (dPar == 7 && ncmt >= 3){
	  //dC1
	  if (in2==0){
	    cur += rate*gamma1*(1.0-exp(-gamma*t1))*exp(-gamma*t2);
	  }
	} else if (dPar == 8 && ncmt >= 3){
	  //dC2
	  if (in2 != 0){
	    cur += rate*gamma1*(1.0-exp(-gamma*t1))*exp(-gamma*t2);
	  }
	} else if (dPar == 9 && ncmt >= 3){
	  //dGamma
	  cur += -C*exp(-t2*gamma)*rate*(1.0 - exp(-t1*gamma))/(gamma*gamma) + C*exp(-t1*gamma - t2*gamma)*t1*rate/gamma - C*exp(-t2*gamma)*t2*rate*(1.0 - exp(-t1*gamma))/gamma;
	} else if (dPar == 10){
	  // dKa
	} else if (dPar == 11){
	  // dTlag1
	  if (in2 == 0){
	    if (thisT < tinf){
	      /* t1 = thisT; */
	      /* t2 = 0; */
	      cur += -A*exp(-thisT*alpha)*rate*alpha*alpha1;
	      if (ncmt >= 2){
		cur +=  -B*exp(-beta*thisT)*rate*beta*beta1;
		if (ncmt >= 3){
		  cur +=  -C*exp(-thisT*gamma)*rate*gamma*gamma1;
		}
	      }
	    } else {
	      /* t1 = tinf; */
	      /* t2 = ThisT-tinf; */
	      cur +=  A*exp(-alpha*(thisT - tinf))*rate*alpha*alpha1*(1.0 - exp(-tinf*alpha));
	      if (ncmt >= 2){
		cur +=  B*exp(-beta*(thisT - tinf))*rate*beta*beta1*(1.0 - exp(-tinf*beta));
		if (ncmt >= 3){
		  cur +=  C*exp(-gamma*(thisT - tinf))*rate*gamma*gamma1*(1.0 - exp(-tinf*gamma));
		}
	      }
	    }
	  }
	} else if (dPar == 12){
	  // dTlag2
	  if (in2 != 0){
	    if (thisT < tinf){
	      /* t1 = thisT; */
	      /* t2 = 0; */
	      cur += -A*exp(-thisT*alpha)*rate*alpha*alpha1;
	      if (ncmt >= 2){
		cur +=  -B*exp(-beta*thisT)*rate*beta*beta1;
		if (ncmt >= 3){
		  cur +=  -C*exp(-thisT*gamma)*rate*gamma*gamma1;
		}
	      }
	    } else {
	      /* t1 = tinf; */
	      /* t2 = ThisT-tinf; */
	      cur +=  A*exp(-alpha*(thisT - tinf))*rate*alpha*alpha1*(1.0 - exp(-tinf*alpha));
	      if (ncmt >= 2){
		cur +=  B*exp(-beta*(thisT - tinf))*rate*beta*beta1*(1.0 - exp(-tinf*beta));
		if (ncmt >= 3){
		  cur +=  C*exp(-gamma*(thisT - tinf))*rate*gamma*gamma1*(1.0 - exp(-tinf*gamma));
		}
	      }
	    }
	  }
	} else if (dPar == 13){
	  // dF1
	  if (in2!=0){
	    if (whI == 1){
	      // Duration where tinfM=tinf*F
	      if (thisT < tinf){
		/* t1 = thisT; */
		/* t2 = 0; */
	      } else {
		/* t1 = tinf; */
		/* t2 = ThisT-tinf; */
		cur +=  A*exp(-alpha*(thisT - tinf) - tinf*alpha)*rate*alpha*alpha1 + A*exp(-alpha*(thisT - tinf))*rate*alpha*alpha1*(1.0 - exp(-tinf*alpha))*tinf/F;
		if (ncmt >= 2){
		  cur +=  B*exp(-beta*(thisT - tinf) - tinf*beta)*rate*beta*beta1 + B*exp(-beta*(thisT - tinf))*rate*beta*beta1*(1.0 - exp(-tinf*beta))*tinf/F;
		  if (ncmt >= 3){
		    cur +=  C*exp(-gamma*(thisT - tinf) - tinf*gamma)*rate*gamma*gamma1 + C*exp(-gamma*(thisT - tinf))*rate*gamma*gamma1*(1.0 - exp(-tinf*gamma))*tinf/F;
		  }
		}
	      }
	    } else {
	      // Rate where rate=F*rate0
	      cur +=  rate*A*alpha1*(1.0-exp(-alpha*t1))*exp(-alpha*t2)/F;
	      if (ncmt >= 2){
		cur +=  rate*B*beta1*(1.0-exp(-beta*t1))*exp(-beta*t2)/F;
		if (ncmt >= 3){
		  cur +=  rate*C*gamma1*(1.0-exp(-gamma*t1))*exp(-gamma*t2)/F;
		}
	      }
	    }
	  }
	} else if (dPar == 14){
	  // dF2
	  if (in2 == 0){
	    if (whI == 1){
	      // Duration where tinfM=tinf*F
	      if (thisT < tinf){
		/* t1 = thisT; */
		/* t2 = 0; */
	      } else {
		/* t1 = tinf; */
		/* t2 = ThisT-tinf; */
		cur +=  A*exp(-alpha*(thisT - tinf) - tinf*alpha)*rate*alpha*alpha1 + A*exp(-alpha*(thisT - tinf))*rate*alpha*alpha1*(1.0 - exp(-tinf*alpha))*tinf/F;
		if (ncmt >= 2){
		  cur +=  B*exp(-beta*(thisT - tinf) - tinf*beta)*rate*beta*beta1 + B*exp(-beta*(thisT - tinf))*rate*beta*beta1*(1.0 - exp(-tinf*beta))*tinf/F;
		  if (ncmt >= 3){
		    cur +=  C*exp(-gamma*(thisT - tinf) - tinf*gamma)*rate*gamma*gamma1 + C*exp(-gamma*(thisT - tinf))*rate*gamma*gamma1*(1.0 - exp(-tinf*gamma))*tinf/F;
		  }
		}
	      }
	    } else {
	      // Rate where rate=F*rate0
	      cur +=  rate*A*alpha1*(1.0-exp(-alpha*t1))*exp(-alpha*t2)/F;
	      if (ncmt >= 2){
		cur +=  rate*B*beta1*(1.0-exp(-beta*t1))*exp(-beta*t2)/F;
		if (ncmt >= 3){
		  cur +=  rate*C*gamma1*(1.0-exp(-gamma*t1))*exp(-gamma*t2)/F;
		}
	      }
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
	if (dPar==0){
	  if (op->linLog){
	    expr1 = log(dose)+log(F);
	    res = ((oral == 1) ? exp(expr1-ka*thisT-log1mex(ka*tau)) : 0.0);
	    cur += A*(exp(expr1-alpha*thisT-log1mex(alpha*tau))-res);
	    if (ncmt >= 2){
	      cur +=  B*(exp(expr1-beta*thisT-log1mex(beta*tau))-res);
	      if (ncmt >= 3){
		cur += C*(exp(expr1-gamma*thisT-log1mex(gamma*tau))-res);
	      }
	    }
	  } else {
	    res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	    cur += dose*F*A*(exp(-alpha*thisT)/(1-exp(-alpha*tau))-res);
	    if (ncmt >= 2){
	      cur +=  dose*F*B*(exp(-beta*thisT)/(1-exp(-beta*tau))-res);
	      if (ncmt >= 3){
		cur += dose*F*C*(exp(-gamma*thisT)/(1-exp(-gamma*tau))-res);
	      }
	    }
	  }
	  // ss=1 is equivalent to a reset + ss dose
	  if (wh0 == 10) return(ret+cur);
	} else if (dPar == 1){
	  // dA1
	  if (in2 == 0){
	    res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	    cur += dose*F*(exp(-alpha*thisT)/(1-exp(-alpha*tau))-res);
	  }
	} else if (dPar == 2) {
	  // dA2
	  if (in2 != 0){
	    res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	    cur += dose*F*(exp(-alpha*thisT)/(1-exp(-alpha*tau))-res);
	  }
	} else if (dPar == 3){
	  // dAlpha
	  res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	  cur += A*F*dose*(-exp(-thisT*alpha)*thisT/(1-exp(-tau*alpha))-exp(-tau*alpha-thisT*alpha)*tau/R_pow_di(1-exp(-tau*alpha),2));
	} else if (dPar == 4 && ncmt >= 2){
	  // dB1
	  if (in2 == 0){
	    res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	    cur += dose*F*(exp(-alpha*thisT)/(1-exp(-alpha*tau))-res);
	  }
	} else if (dPar == 5 && ncmt >= 2){
	  // dB2
	  if (in2 != 0){
	    res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	    cur += dose*F*(exp(-alpha*thisT)/(1-exp(-alpha*tau))-res);
	  }
	} else if (dPar == 6 && ncmt >= 2){
	  // dBeta
	  res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	  cur +=  B*F*dose*(-exp(-beta*thisT)*thisT/(1-exp(-tau*beta))-exp(-beta*thisT-tau*beta)*tau/R_pow_di(1-exp(-tau*beta),2));
	} else if (dPar == 7 && ncmt >= 3){
	  // dC1
	  if (in2 == 0){
	    res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	    cur += dose*F*C*(exp(-gamma*thisT)/(1-exp(-gamma*tau))-res);
	  }
	} else if (dPar == 8 && ncmt >= 3){
	  if (in2 != 0){
	    res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	    cur += dose*F*C*(exp(-gamma*thisT)/(1-exp(-gamma*tau))-res);
	  }
	} else if (dPar == 9 && ncmt >= 3){
	  res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	  cur += C*F*dose*(-exp(-thisT*gamma)*thisT/(1-exp(-tau*gamma))-exp(-tau*gamma-thisT*gamma)*tau/R_pow_di(1-exp(-tau*gamma),2));
	} else if (dPar == 10){
	  /* res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0); */
	  if (oral == 1){
	    cur += A*F*dose*(exp(-ka*thisT)*thisT/(1-exp(-ka*tau))+exp(-ka*tau-ka*thisT)*tau/R_pow_di(1-exp(-ka*tau),2));
	    if (ncmt >= 2){
	      cur +=  B*F*dose*(exp(-ka*thisT)*thisT/(1-exp(-ka*tau))+exp(-ka*tau-ka*thisT)*tau/R_pow_di(1-exp(-ka*tau),2));
	      if (ncmt >= 3){
		cur += C*F*dose*(exp(-ka*thisT)*thisT/(1-exp(-ka*tau))+exp(-ka*tau-ka*thisT)*tau/R_pow_di(1-exp(-ka*tau),2));
	      }
	    }
	  }
	} else if (dPar == 11){
	  // dTlag1
	  // thisT = tT -tlag;
	  if (in2==0){
	    res = ((oral == 1) ? -exp(-ka*thisT)*ka/(1-exp(-ka*tau)) : 0.0);
	    cur += A*F*dose*(-exp(-thisT*alpha)*alpha/(1-exp(-tau*alpha))-res);
	    if (ncmt >= 2){
	      cur += -B*F*dose*(-exp(-beta*thisT)*beta/(1-exp(-tau*beta))-res);
	      if (ncmt >= 3){
		cur += -C*F*dose*(-exp(-thisT*gamma)*gamma/(1-exp(-tau*gamma))-res);
	      }
	    }
	  }
	} else if (dPar == 12){
	  // dTlag2
	  if (in2!=0){
	    res = ((oral == 1) ? -exp(-ka*thisT)*ka/(1-exp(-ka*tau)) : 0.0);
	    cur += A*F*dose*(-exp(-thisT*alpha)*alpha/(1-exp(-tau*alpha))-res);
	    if (ncmt >= 2){
	      cur += -B*F*dose*(-exp(-beta*thisT)*beta/(1-exp(-tau*beta))-res);
	      if (ncmt >= 3){
		cur += -C*F*dose*(-exp(-thisT*gamma)*gamma/(1-exp(-tau*gamma))-res);
	      }
	    }
	  }
	} else if (dPar == 13){
	  // dF1
	  if (in2==0){
	    res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	    cur += dose*A*(exp(-alpha*thisT)/(1-exp(-alpha*tau))-res);
	    if (ncmt >= 2){
	      cur +=  dose*B*(exp(-beta*thisT)/(1-exp(-beta*tau))-res);
	      if (ncmt >= 3){
		cur += dose*C*(exp(-gamma*thisT)/(1-exp(-gamma*tau))-res);
	      }
	    }
	  }
	} else if (dPar == 14){
	  // dF2
	  if (in2!=0){
	    res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	    cur += dose*A*(exp(-alpha*thisT)/(1-exp(-alpha*tau))-res);
	    if (ncmt >= 2){
	      cur +=  dose*B*(exp(-beta*thisT)/(1-exp(-beta*tau))-res);
	      if (ncmt >= 3){
		cur += dose*C*(exp(-gamma*thisT)/(1-exp(-gamma*tau))-res);
	      }
	    }
	  }
	}
      } else if (wh0 == 30) {
	error("You cannot turn off a compartment with a solved system.");
      } else {
	tT = t - ind->all_times[ind->idose[l]];
	thisT = tT -tlag;
	if (dPar == 0){
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
	} else if (dPar == 1){
	  // dA1
	  if (in2==0){
	    res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	    cur +=  dose*F*(exp(-alpha*thisT)-res);
	  }
	} else if (dPar == 2){
	  // dA2
	  if (in2!=0){
	    res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	    cur +=  dose*F*(exp(-alpha*thisT)-res);
	  }
	} else if (dPar == 3){
	  // dAlpha
	  res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	  cur +=  -F*exp(-thisT*alpha)*dose*thisT;
	} else if (dPar == 4 && ncmt >= 2){
	  // dB1
	  if (in2==0){
	    res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	    cur +=  dose*F*(exp(-beta*thisT)-res);
	  }
	} else if (dPar == 5 && ncmt >= 2){
	  // dB2
	  if (in2!=0){
	    res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	    cur +=  dose*F*(exp(-beta*thisT)-res);
	  }
	} else if (dPar == 6 && ncmt >= 2){
	  // dBeta
	  res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	  cur += -B*F*exp(-beta*thisT)*dose*thisT;
	} else if (dPar == 7 && ncmt >= 3){
	  // dC1
	  res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	  if (in2==0){
	    cur += dose*F*(exp(-gamma*thisT)-res);
	  }
	} else if (dPar == 8 && ncmt >= 3) {
	  // dC2
	  res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	  if (in2!=0){
	    cur += dose*F*(exp(-gamma*thisT)-res);
	  }
	} else if (dPar == 9 && ncmt >= 3){
	  // dGamma
	  res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	  cur += -C*F*exp(-thisT*gamma)*dose*thisT;
	} else if (dPar == 10){
	  //dKa
	  if (oral == 1){
	    cur +=  A*F*exp(-ka*thisT)*dose*thisT;
	    if (ncmt >= 2){
	      cur +=  B*F*exp(-ka*thisT)*dose*thisT;
	      if (ncmt >= 3){
		cur += C*F*exp(-ka*thisT)*dose*thisT;
	      }
	    }
	  }
	} else if (dPar == 11){
	  //dTlag1
	  if (in2==0){
	    res = ((oral == 1) ? -exp(-ka*thisT)*ka : 0.0);
	    cur +=  -A*F*dose*(-exp(-thisT*alpha)*alpha - res);
	    if (ncmt >= 2){
	      cur +=  -B*F*dose*(-exp(-beta*thisT)*beta - res);
	      if (ncmt >= 3){
		cur += -C*F*dose*(-exp(-thisT*gamma)*gamma - res);
	      }
	    }
	  }
	} else if (dPar == 12){
	  //dTlag2
	  if (in2!=1){
	    res = ((oral == 1) ? -exp(-ka*thisT)*ka : 0.0);
	    cur +=  -A*F*dose*(-exp(-thisT*alpha)*alpha - res);
	    if (ncmt >= 2){
	      cur +=  -B*F*dose*(-exp(-beta*thisT)*beta - res);
	      if (ncmt >= 3){
		cur += -C*F*dose*(-exp(-thisT*gamma)*gamma - res);
	      }
	    }
	  }
	} else if (dPar == 13){
	  //dF1
	  if (in2==0){
	    res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	    cur +=  dose*A*(exp(-alpha*thisT)-res);
	    if (ncmt >= 2){
	      cur +=  dose*B*(exp(-beta*thisT)-res);
	      if (ncmt >= 3){
		cur += dose*C*(exp(-gamma*thisT)-res);
	      }
	    }
	  }
	} else if (dPar == 14){
	  //dF2
	  if (in2!=0){
	    res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	    cur +=  dose*A*(exp(-alpha*thisT)-res);
	    if (ncmt >= 2){
	      cur +=  dose*B*(exp(-beta*thisT)-res);
	      if (ncmt >= 3){
		cur += dose*C*(exp(-gamma*thisT)-res);
	      }
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
