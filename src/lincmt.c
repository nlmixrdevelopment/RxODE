#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include "solve.h"

/* extern SEXP RxODE_ode_dosing(); */
// extern "C" double getLinDeriv(int ncmt, int diff1, int diff2, double rate, double tinf, double Dose, double ka, double tlag, double T, double tT, mat par);

extern double RxODE_prodV_r(double *input, double *p, int type, int n, ...);
extern double RxODE_sumV_r(double *p, long double *pld, int m, int type, int n, ...);
extern double RxODE_safe_zero(double);

extern double rxDosingTimeP(int i, rx_solve *rx, unsigned int id);
extern unsigned int nDosesP(rx_solve *rx, unsigned int id);
extern void setExtraCmtP(int xtra, rx_solve *rx);
extern int rxDosingEvidP(int i, rx_solve *rx, unsigned int id);
extern double rxDoseP(int i, rx_solve *rx, unsigned int id);
extern double RxODE_solveLinB(rx_solve *rx, unsigned int id, double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);

// Type 1 = PairwiseSum
#define sum(...) RxODE_sumV_r(ps, pld, -6, 1, __VA_ARGS__)
// Type 3 = Logify
// Type 2 = product
#define prod(...) RxODE_prodV_r(ps, pi, 2, __VA_ARGS__)
#define safe_zero RxODE_safe_zero


int locateDoseIndex(const double obs_time, rx_solve *rx, unsigned int id){
  // Uses bisection for slightly faster lookup of dose index.
  int i, j, ij;
  i = 0;
  j = nDosesP(rx, id) - 1;
  if (obs_time <= rxDosingTimeP(i, rx, id)){
    return i;
  }
  if (obs_time >= rxDosingTimeP(j, rx, id)){
    return j;
  }
  while(i < j - 1) { /* x[i] <= obs_time <= x[j] */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    if(obs_time < rxDosingTimeP(ij, rx, id))
      j = ij;
    else
      i = ij;
  }
  return i;
}//subscript of dose

double RxODE_solveLinB(rx_solve *rx, unsigned int id, double t, int linCmt, int diff1, int diff2, double d_A, double d_alpha, double d_B, double d_beta, double d_C, double d_gamma, double d_ka, double d_tlag){
  if (diff1 != 0 || diff2 != 0){
    error("Exact derivtives are no longer calculated.");
  }
  double ps[6], pi[6];
  long double pld[6];
  unsigned int ncmt = 1;
  if (d_gamma > 0.){
    ncmt = 3;
  } else if (d_beta > 0.){
    ncmt = 2;
  } else if (d_alpha > 0.){
    ncmt = 1;
  } else {
    return 0.0;
    //error("You need to specify at least A(=%f) and alpha (=%f). (@t=%f, d1=%d, d2=%d)", d_A, d_alpha, t, diff1, diff2);
  }
  setExtraCmtP(linCmt+1, rx);
  
  double alpha = d_alpha;
  double A = d_A;
  double beta = d_beta;
  double B = d_B;
  double gamma = d_gamma;
  double C = d_C;
  double ka = d_ka;
  double tlag = d_tlag;
  
  int oral, cmt;
  oral = (ka > 0) ? 1 : 0;
  double ret = 0;
  unsigned int m = 0, l = 0, p = 0;
  int evid, evid100;
  double thisT = 0.0, tT = 0.0, res, t1, t2, tinf, dose = 0;
  double rate;
  m = locateDoseIndex(t, rx, id);
  
  for(l=0; l <= m; l++){
    //superpostion
    evid = rxDosingEvidP(l, rx, id);
    dose = rxDoseP(l, rx, id);
    // Support 100+ compartments...
    evid100 = floor(evid/1e5);
    evid = evid- evid100*1e5;
    cmt = (evid%10000)/100 - 1 + 100*evid100;
    if (cmt != linCmt) continue;
    if (evid > 10000) {
      if (dose > 0){
	// During infusion
	tT = sum(2, t, - rxDosingTimeP(l, rx, id));
	thisT = sum(2, tT, - tlag);
	p = l+1;
	while (p < nDosesP(rx, id) && rxDoseP(p, rx, id) != -dose){
	  p++;
	}
	if (rxDoseP(p,rx,id) != -dose){
	  error("Could not find a error to the infusion.  Check the event table.");
	}
	tinf  = sum(2,rxDosingTimeP(p,rx,id),-rxDosingTimeP(l,rx,id));
	rate  = dose;
	if (tT >= tinf) continue;
      } else {
	// After  infusion
	p = l-1;
	while (p > 0 && rxDoseP(p,rx,id) != -dose){
	  p--;
	}
	if (rxDoseP(p,rx, id) != -dose){
	  error("Could not find a start to the infusion.  Check the event table.");
	}
	tinf  = sum(3, rxDosingTimeP(l,rx,id), - rxDosingTimeP(p,rx,id), - tlag);
	
	tT = sum(2, t, - rxDosingTimeP(p, rx, id));
        thisT = sum(2, tT, -tlag);

	rate  = -dose;
      }
      t1 = ((thisT < tinf) ? thisT : tinf);        //during infusion
      t2 = ((thisT > tinf) ? sum(2, thisT, - tinf) : 0.0);  // after infusion
      ret +=  prod(5, rate,A,1.0/safe_zero(alpha),sum(2, 1.0,-exp(prod(2, -alpha, t1))),exp(prod(2, -alpha,t2)));
      if (ncmt >= 2){
	ret +=  prod(5, rate,B , 1.0/safe_zero(beta),sum(2, 1.0,-exp(prod(2, -beta, t1))),exp(prod(2, -beta,t2)));
	if (ncmt >= 3){
	  ret +=  prod(5, rate,C , 1.0/safe_zero(gamma),sum(2, 1.0,-exp(prod(2, -gamma, t1))),exp(prod(2, -gamma,t2)));
	}
      }
    } else {
      tT = sum(2, t, -rxDosingTimeP(l, rx, id));
      thisT = sum(2, tT, -tlag);
      if (thisT < 0) continue;
      res = ((oral == 1) ? exp(prod(2, -ka , thisT)) : 0.0);
      ret +=  prod(3, dose, A, sum(2, exp(prod(2, -alpha, thisT)),-res));
      if (ncmt >= 2){
	ret +=  prod(3, dose, B, sum(2, exp(prod(2, -beta, thisT)),-res));
        if (ncmt >= 3){
	  ret +=  prod(3, dose, C, sum(2, exp(prod(2, -gamma, thisT)),-res));
        }
      }
    }
  } //l
  return ret;
}

#undef sum
#undef prod
