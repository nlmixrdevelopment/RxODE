#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>

/* extern SEXP RxODE_ode_dosing(); */
// extern "C" double getLinDeriv(int ncmt, int diff1, int diff2, double rate, double tinf, double Dose, double ka, double tlag, double T, double tT, mat par);
extern double rxSolveLinBdInf(int diff1, int diff2, int dA, int dAlpha, double rate, double tT, double t1, double t2, double tinf, double A, double alpha, double tlag);
extern double rxSolveLinBDiff(int diff1, int diff2, int dA, int dAlpha, double dose, double tT, double A, double alpha, double ka, double tlag);
extern double RxODE_solveLinB(double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);
extern double RxODE_sum (double *input, unsigned int n);
extern unsigned int nDoses();
extern double rxDosingTime(int i);
extern int rxDosingEvid(int i);
extern double RxODE_prodV(unsigned int n, ...);
extern double RxODE_sumV(unsigned int n, ...);
extern double RxODE_safe_zero(double);
extern void setExtraCmt(int xtra);

#define sum RxODE_sumV
#define prod RxODE_prodV
#define safe_zero RxODE_safe_zero
double rxDose(int i);

int locateDoseIndex(const double obs_time){
  // Uses bisection for slightly faster lookup of dose index.
  int i, j, ij;
  i = 0;
  j = nDoses() - 1;
  if (obs_time <= rxDosingTime(i)){
    return i;
  }
  if (obs_time >= rxDosingTime(j)){
    return j;
  }
  while(i < j - 1) { /* x[i] <= obs_time <= x[j] */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    if(obs_time < rxDosingTime(ij))
      j = ij;
    else
      i = ij;
  }
  return i;
}//subscript of dose

#define NUM_PARTIALS  32  /* initial partials array size, on stack */

double RxODE_solveLinB(double t, int linCmt, int diff1, int diff2, double d_A, double d_alpha, double d_B, double d_beta, double d_C, double d_gamma, double d_ka, double d_tlag){
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
  setExtraCmt(linCmt+1);
  
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
  int evid;
  double thisT = 0.0, tT = 0.0, res, t1, t2, tinf, dose = 0;
  double rate;
  m = locateDoseIndex(t);
  unsigned int nsum = 0, np = NUM_PARTIALS;
  double *sumar = Calloc(NUM_PARTIALS, double);
  
  for(l=0; l <= m; l++){
    if (nsum >= NUM_PARTIALS){
      np += np;
      sumar = Realloc(sumar, np, double);
    }
    //superpostion
    evid = rxDosingEvid(l);
    dose = rxDose(l);
    cmt = (evid%10000)/100 - 1;
    if (cmt != linCmt) continue;
    if (evid > 10000) {
      if (dose > 0){
	// During infusion
	tT = sum(2, t, - rxDosingTime(l));
	thisT = sum(2, tT, - tlag);
	p = l+1;
	while (p < nDoses() && rxDose(p) != -dose){
	  p++;
	}
	if (rxDose(p) != -dose){
	  error("Could not find a error to the infusion.  Check the event table.");
	}
	tinf  = sum(2,rxDosingTime(p),-rxDosingTime(l));
	rate  = dose;
	if (tT >= tinf) continue;
      } else {
	// After  infusion
	p = l-1;
	while (p >= 0 && rxDose(p) != -dose){
	  p--;
	}
	if (rxDose(p) != -dose){
	  error("Could not find a start to the infusion.  Check the event table.");
	}
	tinf  = sum(3, rxDosingTime(l), - rxDosingTime(p), - tlag);
	
	tT = sum(2, t, - rxDosingTime(p));
        thisT = sum(2, tT, -tlag);

	rate  = -dose;
      }
      t1 = ((thisT < tinf) ? thisT : tinf);        //during infusion
      t2 = ((thisT > tinf) ? sum(2, thisT, - tinf) : 0.0);  // after infusion
      if (diff1 == 0 && diff2 == 0){
	sumar[nsum++] = prod(5, rate,A,1.0/safe_zero(alpha),sum(2, 1.0,-exp(prod(2, -alpha, t1))),exp(prod(2, -alpha,t2)));
      } else {
	sumar[nsum++] = rxSolveLinBdInf(diff1, diff2, 1, 2, rate, tT, t1, t2, tinf, A, alpha, tlag);
      }
      if (ncmt >= 2){
	if (diff1 == 0 && diff2 == 0){
	  sumar[nsum++] = prod(5, rate,B , 1.0/safe_zero(beta),sum(2, 1.0,-exp(prod(2, -beta, t1))),exp(prod(2, -beta,t2)));
	} else {
	  sumar[nsum++] = rxSolveLinBdInf(diff1, diff2, 3, 4, rate, tT, t1, t2, tinf, B, beta, tlag);
	}
	if (ncmt >= 3){
	  if (diff1 == 0 && diff2 == 0){
	    sumar[nsum++] = prod(5, rate,C , 1.0/safe_zero(gamma),sum(2, 1.0,-exp(prod(2, -gamma, t1))),exp(prod(2, -gamma,t2)));
          } else {
	    sumar[nsum++] = rxSolveLinBdInf(diff1, diff2, 5, 6, rate, tT, t1, t2, tinf, C, gamma, tlag);
	  }
	}
      }
    } else {
      tT = sum(2, t, -rxDosingTime(l));
      thisT = sum(2, tT, -tlag);
      if (thisT < 0) continue;
      res = ((oral == 1) ? exp(prod(2, -ka , thisT)) : 0.0);
      if (diff1 == 0 && diff2 == 0){
	sumar[nsum++] = prod(3, dose, A, sum(2, exp(prod(2, -alpha, thisT)),-res));
      } else {
	sumar[nsum++] = rxSolveLinBDiff(diff1, diff2, 1, 2, dose, tT, A, alpha, ka, tlag);
      }
      if (ncmt >= 2){
	if (diff1 == 0 && diff2 == 0){
          sumar[nsum++] = prod(3, dose, B, sum(2, exp(prod(2, -beta, thisT)),-res));
	} else {
	  sumar[nsum++] = rxSolveLinBDiff(diff1, diff2, 3, 4, dose, tT, B, beta, ka, tlag);
	}
        if (ncmt >= 3){
	  if (diff1 == 0 && diff2 == 0){
            sumar[nsum++] = prod(3, dose, C, sum(2, exp(prod(2, -gamma, thisT)),-res));
	  } else {
	    sumar[nsum++] = rxSolveLinBDiff(diff1, diff2, 5, 6, dose, tT, C, gamma, ka, tlag);
          }
        }
      }
    }
  } //l
  for (m = 0; m < nsum;m++){
    if (ISNAN(sumar[m])){
      error("NaN produced at A(=%f) and alpha (=%f). (@t=%f, d1=%d, d2=%d)", d_A, d_alpha, t, diff1, diff2);
    }
  }
  ret = RxODE_sum(sumar, nsum);
  Free(sumar);
  return ret;
}
