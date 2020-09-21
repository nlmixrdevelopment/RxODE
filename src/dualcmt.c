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

#include "dual.h"
#include "lincmtB1g.h"
#include "lincmtB2g.h"
#include "lincmtB3g.h"


int handle_evidL(int evid, double *yp, double xout, int id, rx_solving_options_ind *ind);
double getTime(int idx, rx_solving_options_ind *ind);
double _getDur(int l, rx_solving_options_ind *ind, int backward, unsigned int *p);


typedef struct parTr {
  int trans;
  int ncmt;
  int oral0;
  int linCmt;
  rx_solving_options_ind *ind;
  double d_tlag;
  double d_F;
  double d_rate1;
  double d_dur1;
  // Oral parameters
  double d_tlag2;
  double d_F2;
  double d_rate2;
  double d_dur2;
  // Input parameters 
  dualN p1;
  dualN v1;
  dualN p2;
  dualN p3;
  dualN p4;
  dualN p5;
  dualN ka;
  // Output parameters
  dualN rx_k;
  dualN rx_v;
  dualN rx_k12;
  dualN rx_k21;
  dualN rx_k13;
  dualN rx_k31;
} parTr;

#undef beta
static inline parTr parTrans(int *trans, 
			     double *p1, double *v1,
			     double *p2, double *p3,
			     double *p4, double *p5,
			     double *ka, rx_solving_options_ind *ind,
			     int linCmt,
			     double d_tlag, double d_F, double d_rate1, double d_dur1,
			     // Oral parameters
			     double d_tlag2, double d_F2,  double d_rate2, double d_dur2){
  parTr ret;
  ret.ind = ind;
  ret.linCmt = linCmt;
  ret.d_tlag = d_tlag;
  ret.d_F = d_F;
  ret.d_rate1 = d_rate1;
  ret.d_dur1 = d_dur1;
  // Oral parameters
  ret.d_tlag2 = d_tlag2;
  ret.d_F2 = d_F2;
  ret.d_rate2 = d_rate2;
  ret.d_dur2 = d_dur2;
  ret.ka = iniD(*ka, dKa);
  ret.p1 = iniD(*p1, dP1);
  ret.p2 = iniD(*p2, dP2);
  ret.p3 = iniD(*p3, dP3);
  ret.p4 = iniD(*p4, dP4);
  ret.p5 = iniD(*p5, dP5);
  ret.v1 = iniD(*v1, dV1);
  ret.oral0 = (*ka > 0) ? 1 : 0;
  dualN alpha, beta, gamma, A, B, C, btemp, ctemp, dtemp;
  if ((*p5) > 0.) {
    ret.ncmt = 3;
    ret.trans = *trans;
    switch (*trans) {
    case 1: // cl v q vp
      ret.rx_k = div2(ret.p1,ret.v1); // k = CL/V
      ret.rx_v = ret.v1;
      ret.rx_k12 = div2(ret.p2,ret.v1); // k12 = Q/V
      ret.rx_k21 = div2(ret.p2,ret.p3); // k21 = Q/Vp
      ret.rx_k13 = div2(ret.p4,ret.v1); // k31 = Q2/V
      ret.rx_k31 = div2(ret.p4,ret.p5); // k31 = Q2/Vp2
      break;
    case 2: // k=(*p1) v=(*v1) k12=(*p2) k21=(*p3) k13=(*p4) k31=(*p5)
      ret.rx_k = ret.p1;
      ret.rx_v = ret.v1;
      ret.rx_k12 = ret.p2;
      ret.rx_k21 = ret.p3;
      ret.rx_k13 = ret.p4;
      ret.rx_k31 = ret.p5;
      break;
    case 10:
    case 11:
      if (*trans == 11) A = divd2(1,ret.v1);
      else A = ret.v1;
      B  = ret.p3;
      C =  ret.p5;
      alpha = ret.p1;
      beta = ret.p2;
      gamma = ret.p4;
      ret.rx_v= divd2(1,(add2(add2(A,B),C)));
      btemp = prod2(negD((add2(add2(add2(add2(add2(prod2(alpha,C),prod2(alpha,B)),prod2(gamma,A)),prod2(gamma,B)),prod2(beta,A)),prod2(beta,C)))),ret.rx_v);
      ctemp = prod2((add2(add2(prod2(prod2(alpha,beta),C),prod2(prod2(alpha,gamma),B)),prod2(prod2(beta,gamma),A))),ret.rx_v);
      dtemp = sqrtD(subtr2(prod2(btemp,btemp),prodd2(4,ctemp)));
      ret.rx_k21 = prodd2(0.5,(add2(negD(btemp),dtemp)));
      ret.rx_k31 = prodd2(0.5,(subtr2(negD(btemp),dtemp)));
      ret.rx_k   = div2(div2(prod2(prod2(alpha,beta),gamma),ret.rx_k21),ret.rx_k31);
      ret.rx_k12 = div2((add2(subtr2(subtr2((add2(add2(prod2(beta,gamma),prod2(alpha,beta)),prod2(alpha,gamma))),
					    prod2(ret.rx_k21,(add2(add2(alpha,beta),gamma)))),
				     prod2(ret.rx_k,ret.rx_k31)),prod2(ret.rx_k21,ret.rx_k21))),
			(subtr2(ret.rx_k31,ret.rx_k21)));
      ret.rx_k13 = subtr2(add2(add2(alpha,beta),gamma),(add2(add2(add2(ret.rx_k,ret.rx_k12),ret.rx_k21),ret.rx_k31)));
      break;
    default:
      /* REprintf(_("invalid trans (3 cmt trans %d)\n"), *trans); */
      Rf_errorcall(R_NilValue, _("invalid translation"));
    }
  } else if ((*p3) > 0.) {
    ret.ncmt = 2;
    ret.trans = *trans;
    switch (*trans){
    case 1: // cl=(*p1) v=(*v1) q=(*p2) vp=(*p3)
      ret.rx_k = div2(ret.p1,ret.v1); // k = CL/V
      ret.rx_v = ret.v1;
      ret.rx_k12 = div2(ret.p2,ret.v1); // k12 = Q/V
      ret.rx_k21 = div2(ret.p2,ret.p3); // k21 = Q/Vp
      break;
    case 2: // k=(*p1), (*v1)=v k12=(*p2) k21=(*p3)
      ret.rx_k = ret.p1;
      ret.rx_v = ret.v1;
      ret.rx_k12 = ret.p2;
      ret.rx_k21 = ret.p3;
      break;
    case 3: // cl=(*p1) v=(*v1) q=(*p2) vss=(*p3)
      ret.rx_k = div2(ret.p1,ret.v1); // k = CL/V
      ret.rx_v = ret.v1;
      ret.rx_k12 = div2(ret.p2, ret.v1); // k12 = Q/V
      ret.rx_k21 = div2(ret.p2,(subtr2(ret.p3,ret.v1))); // k21 = Q/(Vss-V)
      break;
    case 4: // alpha=(*p1) beta=(*p2) k21=(*p3)
      ret.rx_v = ret.v1;
      ret.rx_k21 = ret.p3;
      ret.rx_k = div2(prod2(ret.p1,ret.p2),ret.rx_k21); // (*p1) = alpha (*p2) = beta
      ret.rx_k12 = subtr2(subtr2(add2(ret.p1,ret.p2),ret.rx_k21),ret.rx_k);
      break;
    case 5: // alpha=(*p1) beta=(*p2) aob=(*p3)
      ret.rx_v= ret.v1;
      ret.rx_k21 = div2((add2(prod2(ret.p3,ret.p2),ret.p1)),(add2d(ret.p3,1.0)));
      ret.rx_k = div2((prod2(ret.p1,ret.p2)),ret.rx_k21);
      ret.rx_k12 = subtr2(subtr2(add2(ret.p1,ret.p2),ret.rx_k21),ret.rx_k);
      break;
    case 10:
    case 11: // A2 V, alpha=(*p1), beta=(*p2), k21
      if (*trans == 11) A  = divd2(1,ret.v1);
      else A  = ret.v1;
      B = ret.p3;
      alpha = ret.p1;
      beta =ret.p2;
      ret.rx_v   = divd2(1,(add2(A,B)));
      ret.rx_k21 = prod2((add2(prod2(A,beta),prod2(B,alpha))),ret.rx_v);
      ret.rx_k   = div2(prod2(alpha,beta),ret.rx_k21);
      ret.rx_k12 = subtr2(subtr2(add2(alpha,beta),ret.rx_k21),ret.rx_k);
      break;
    default:
      /* REprintf(_("invalid trans (2 cmt trans %d)\n"), trans); */
      Rf_errorcall(R_NilValue, _("invalid translation"));
    }
  } else if ((*p1) > 0.) {
    ret.ncmt = 1;
    ret.trans = *trans;
    switch(*trans){
    case 1: // cl v
      ret.rx_k = div2(ret.p1,ret.v1); // k = CL/V
      ret.rx_v = ret.v1;
      break;
    case 2: // k V
      ret.rx_k = ret.p1;
      ret.rx_v = ret.v1;
      break;
    case 11: // alpha V
      ret.rx_k = ret.p1;
      ret.rx_v = ret.v1;
      break;
    case 10: // alpha A
      ret.rx_k = ret.p1;
      ret.rx_v = divd2(1,ret.v1);
      break;
    default:
      Rf_errorcall(R_NilValue, _("invalid translation"));
    }
  } else {
    Rf_errorcall(R_NilValue, _("invalid translation"));
  }
  return ret;
}

static inline dualN ssRateTauG(double *A,
			       double *tinf,
			       double *tau,
			       parTr *tr,
			       double *r1,
			       double *r2){
  if (tr->oral0){
    if ((*r1) > 0 ){
      switch (tr->ncmt){
      case 1: {
	return oneCmtKaRateSStr1G(A, tinf, tau, r1, tr->ka, tr->rx_k);
      } break;
      case 2: {
	return twoCmtKaRateSStr1G(A, tinf, tau, r1, tr->ka, tr->rx_k, tr->rx_k12, tr->rx_k21);
      } break;
      case 3: {
	return threeCmtKaRateSStr1G(A, tinf, tau, r1, tr->ka, tr->rx_k,
				    tr->rx_k12, tr->rx_k21, tr->rx_k13, tr->rx_k31);
      } break;
      }
    } else {
      switch (tr->ncmt){
      case 1: {
	return  oneCmtKaRateSStr2G(A, tinf, tau, r2, tr->ka, tr->rx_k);
      } break;
      case 2: {
	return twoCmtKaRateSStr2G(A, tinf, tau, r2, tr->ka, tr->rx_k, tr->rx_k12, tr->rx_k21);
      } break;
      case 3: {
	return threeCmtKaRateSStr2G(A, tinf, tau, r2, tr->ka, tr->rx_k, tr->rx_k12,  tr->rx_k21, tr->rx_k13, tr->rx_k31);
      } break;
      }
    }
  } else {
    switch (tr->ncmt){
    case 1: {
      return oneCmtRateSSG(A, tinf, tau, r1, tr->rx_k);
    } break;
    case 2: {
      return twoCmtRateSSG(A, tinf, tau, r1, tr->rx_k, tr->rx_k12, tr->rx_k21);
    } break;
    case 3: {
      return threeCmtRateSSG(A, tinf, tau, r1, tr->rx_k, tr->rx_k12,  tr->rx_k21, tr->rx_k13,  tr->rx_k31);
    } break;
    }
  }
  Rf_errorcall(R_NilValue, _("invalid model"));
  dualN ret;
  return ret;
}


static inline dualN ssTauG(double *A,
			   double *tau,
			   parTr *tr,
			   double *b1,
			   double *b2){
  if (tr->oral0){
    if ((*b1) > 0 ){
      switch (tr->ncmt){
      case 1: {
	return oneCmtKaSSb1G(A, tau, b1, tr->ka, tr->rx_k);
      } break;
      case 2: {
	return twoCmtKaSSb1G(A, tau, b1, tr->ka, tr->rx_k, tr->rx_k12, tr->rx_k21);
      } break;
      case 3: {
	return threeCmtKaSSb1G(A, tau, b1, tr->ka, tr->rx_k, tr->rx_k12,  tr->rx_k21, tr->rx_k13, tr->rx_k31);
      } break;
      }
    } else {
      switch (tr->ncmt){
      case 1: {
	return oneCmtKaSSb2G(A, tau, b2, tr->ka, tr->rx_k);
      } break;
      case 2: {
	return twoCmtKaSSb2G(A, tau, b2, tr->ka, tr->rx_k, tr->rx_k12, tr->rx_k21);
      } break;
      case 3: {
	return threeCmtKaSSb2G(A, tau, b2, tr->ka, tr->rx_k, tr->rx_k12,  tr->rx_k21, tr->rx_k13, tr->rx_k31);
      } break;
      }
    }
  } else {
    switch (tr->ncmt){
    case 1: {
      return oneCmtBolusSSG(A, tau, b1, tr->rx_k);
    } break;
    case 2: {
      return twoCmtBolusSSG(A, tau, b1, tr->rx_k, tr->rx_k12, tr->rx_k21);
    } break;
    case 3: {
      return threeCmtBolusSSG(A, tau, b1, tr->rx_k, tr->rx_k12,  tr->rx_k21, tr->rx_k13,  tr->rx_k31);
    } break;
    }
  }
}


static inline dualN handleSSLG(double *A,// Amounts
			       double *Alast, // Last amounts
			       double tlast, // Time of last amounts
			       double ct, // Time of the dose
			       parTr *tr,
			       double *b1, // Amount of the dose in compartment #1
			       double *b2, // Amount of the dose in compartment #2
			       double *r1, // Rate in Compartment #1
			       double *r2, // Rate in Compartment #2
			       int *nSave,
			       double *aSave,
			       dualN ret0
			       ) {
  rx_solving_options_ind *ind = tr->ind;
  // handle_evid has been called, so ind->wh0 and like have already been called
  double *rate = ind->linCmtRate;
  // note ind->ixds has already advanced
  double amt = ind->dose[ind->ixds-1];
  dualN ret = ret0;
  switch(ind->wh0) {
  case 40: { // Steady state constant infusion
    // Already advanced ind->ixds
    rate[0] =*r1  = *r2 = 0;
    if (tr->oral0){
      rate[1] = 0;
    }
    int cmtOff = ind->cmt-tr->linCmt;
    if (cmtOff == 0){
      // Infusion to central compartment with oral dosing
      *r1 = amt;
    } else {
      // Infusion to central compartment or depot
      *r2 = amt;
    }
    /* ssRateD(A, ncmt, oral0, r1, r2, */
    /* 	    ka, kel, k12, k21, k13, k31); */
    //
  } break;
  case 20: // Steady state + last observed event
  case 10: { // Steady state
    double tau = ind->ii[ind->ixds-1];
    rate[0] =*r1  = *r2 = 0;
    if (tr->oral0){
      rate[1] = 0;
    }
    int cmtOff = ind->cmt-tr->linCmt;
    //*tlast=0; // is this necessary
    switch (ind->whI) {
    case 0: { // bolus dose
      if (cmtOff == 0){
	*b1 = amt*(tr->d_F);
	*b2 = 0;
      } else {
	*b1 = 0;
	*b2 = amt*(tr->d_F2);
      }
      ret = ssTauG(A, &tau, tr, b1, b2);
    } break;
    case 8: // Duration is modeled
    case 9: { // Rate is modeled
      double tinf;
      if (ind->whI == 9) {
	if (cmtOff == 0){
	  // Infusion to central compartment with oral dosing
	  *r1 = tr->d_rate1;
	  tinf = amt*(*r1)/(tr->d_F);
	  rate[0] = *r1;
	} else {
	  // Infusion to central compartment or depot
	  *r2 = tr->d_rate2;
	  tinf = amt*(*r2)/(tr->d_F2);
	  rate[1] = *r2;
	}
      } else {
	// duration is modeled
	if (cmtOff == 0) {
	  // With oral dosing infusion to central compartment
	  tinf = tr->d_dur1;
	  *r1 = amt/tinf*(tr->d_F);
	  rate[0] = *r1;
	} else {
	  // Infusion to compartment #1 or depot
	  tinf = tr->d_dur2;
	  *r2 = amt/tinf*(tr->d_F2);
	  rate[1] = *r2;
	}
      }
      if (tinf >= tau){
	ind->wrongSSDur=1;
	for (int i = tr->ncmt + tr->oral0; i--;){
	  A[i] += R_NaN;
	}
      } else {
	ret = ssRateTauG(A, &tinf, &tau, tr, r1, r2);
      }
    } break;
    case 1: // Infusion rate is fixed
    case 2: { // Infusion, duration is fixed
      double tinf;
      if (ISNA(amt)){
      } else if (amt > 0) {
	unsigned int p;
	tinf = _getDur(ind->ixds-1, ind, 1, &p);
	if (ind->whI == 1){
	  // Duration changes with F
	  if (cmtOff == 0){
	    tinf *= (tr->d_F);
	    *r1 = amt;
	    rate[0] = *r1;
	  } else {
	    tinf *= (tr->d_F2);
	    *r2 = amt;
	    rate[1] = *r2;
	  }
	} else {
	  // rate changes with F; tinf remains constant
	  if (cmtOff == 0){
	    *r1 = amt*(tr->d_F);
	    rate[0] = *r1;
	  } else {
	    *r2 = amt*(tr->d_F2);
	    rate[1] = *r2;
	  }
	}
	if (tinf >= tau){
	  ind->wrongSSDur=1;
	  for (int i = tr->ncmt + tr->oral0; i--;){
	    A[i] += R_NaN;
	  }
	} else {
	  ret = ssRateTauG(A, &tinf, &tau, tr, r1, r2);
	}
      }
    } break;
    }
    // Add back saved values
    if (ind->wh0 == 20) {
      if (amt > 0) {
	for (int i = *nSave; i--;){
	  A[i] += aSave[i];
	}
      }
    }
  } break;
  }
  return ret;
}

static inline dualN doAdvanG(double *A,// Amounts
			     double *Alast, // Last amounts
			     double tlast, // Time of last amounts
			     double ct, // Time of the dose
			     parTr *tr,
			     double *b1, // Amount of the dose in compartment #1
			     double *b2, // Amount of the dose in compartment #2
			     double *r1, // Rate in Compartment #1
			     double *r2 // Rate in Compartment #2
			     ) {
  double t = ct - tlast;
  if (tr->oral0) {
    switch (tr->ncmt) {
    case 1: {
      return oneCmtKaRateG(A, Alast, &t, b1, b2, r1, r2, tr->ka, tr->rx_k);
    } break;
    case 2: {
      return twoCmtKaRateG(A, Alast, &t, b1, b2, r1, r2,
			   tr->ka,  tr->rx_k, tr->rx_k12, tr->rx_k21);
    } break;
    case 3: {
      return threeCmtKaRateG(A, Alast, &t, b1, b2, r1, r2,
  		      tr->ka,  tr->rx_k, tr->rx_k12, tr->rx_k21, tr->rx_k13, tr->rx_k31);
    } break;
    }
  } else {
    switch (tr->ncmt){
    case 1: {
      return oneCmtRateG(A, Alast, &t, b1, r1, tr->rx_k);
    } break;
    case 2: {
      return twoCmtRateG(A, Alast, &t, b1, r1,
  		  tr->rx_k, tr->rx_k12, tr->rx_k21);
    } break;
    case 3: {
      return threeCmtRateG(A, Alast, &t, b1, r1,
  		    tr->rx_k, tr->rx_k12, tr->rx_k21, tr->rx_k13, tr->rx_k31);
    }
    }
  }
  Rf_errorcall(R_NilValue, _("invalid model"));
  dualN ret;
  return ret;
}

double linCmtG(rx_solve *rx, unsigned int id, double t, int linCmt,
	       int i_cmt, int trans, int val,
	       double p1, double v1,
	       double p2, double p3,
	       double p4, double p5,
	       double d_tlag, double d_F, double d_rate1, double d_dur1,
	       // Oral parameters
	       double d_ka, double d_tlag2, double d_F2,  double d_rate2, double d_dur2) {
  rx_solving_options_ind *ind = &(rx->subjects[id]);
  int evid;
  /* evid = ind->evid[ind->ix[ind->idx]]; */
  /* if (evid) REprintf("evid0[%d:%d]: %d; curTime: %f\n", id, ind->idx, evid, t); */
  int idx = ind->idx;
  double Alast0[31] = {0.0, 0.0, 0.0, 0.0, 0.0, //5
		       0.0, 0.0, 0.0, 0.0, 0.0, //10
		       0.0, 0.0, 0.0, 0.0, 0.0, //15
		       0.0, 0.0, 0.0, 0.0, 0.0, //20
		       0.0, 0.0, 0.0, 0.0, 0.0, //25
		       0.0, 0.0, 0.0, 0.0, 0.0, //30
		       0.0}; //31
  rx_solving_options *op = rx->op;
  double *A;
  double *Alast;
  /* A = Alast0; Alast=Alast0; */
  double tlast;
  double curTime= getTime(ind->ix[idx], ind); // t0
  while (t < curTime) {
    idx--;
    if (idx < 0) return 0.0;
    curTime = getTime(ind->ix[idx], ind);
  }
  int sameTime = isSameTime(t, curTime);
  
  parTr tr =  parTrans(&trans, &p1, &v1, &p2, &p3, &p4, &p5,
		       &d_ka, ind, linCmt,
		       d_tlag, d_F, d_rate1, d_dur1, d_tlag2, d_F2,  d_rate2, d_dur2);
  if (tr.ncmt != i_cmt) {
    Rf_errorcall(R_NilValue, _("cmt mismatch"));
  }
  dualN ret;
  double *rate = ind->linCmtRate;
  double b1=0, b2=0, r1 = 0, r2 = 0;
  A = getAdvan(idx);
  if (idx == 0) {
    Alast = Alast0;
    tlast = getTime(ind->ix[0], ind);
  } else {
    tlast = getTime(ind->ix[idx-1], ind);
    Alast = getAdvan(idx-1);
  }
  curTime = getTime(ind->ix[idx], ind);
  evid = ind->evid[ind->ix[idx]];
  if (op->nlinR == 2){
    r1 = rate[0];
    r2 = rate[1];
  } else {
    r1 = rate[0];
  }
  ind->wh0 = 0;
  if (evid == 3){
    // Reset event
    Alast=Alast0;
  } else {
    ret = doAdvanG(A, Alast, tlast, // Time of last amounts
	     curTime, &tr, &b1, &b2, &r1, &r2);
    double aSave[31];
    int nSave = 0;
    switch(tr.ncmt) {
    case 1:
      nSave = (tr.oral0 ? 8 : 3);
      break;
    case 2:
      nSave = (tr.oral0 ? 23 : 10);
      break;
    case 3:
      nSave = (tr.oral0 ? 32 : 21);
    }
    for (int i = nSave; i--;){
      aSave[i] = A[i];
    }
    if (handle_evidL(evid, A, curTime, id, ind)){
      ret = handleSSLG(A, Alast, tlast, curTime, &tr,
		       &b1, &b2, &r1, &r2, &nSave, aSave, ret);

    }
  }
  if (!sameTime){
    // Compute the advan solution of a t outside of the mesh.
    Alast = A;
    /* Ac = Alast0; */
    tlast = curTime;
    curTime = t;
    b1 = b2 = 0;
    ret = doAdvanG(A, Alast, tlast, // Time of last amounts
		   curTime, &tr,  &b1, &b2, &r1, &r2);
  }
  ret = div2(ret,tr.rx_v);
  // Now the derivatives are in ret.grad[0-6]
  // and ret.f = cp
  return ret.f;
  return R_NaN;
}