// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include "RxODE_types.h"

using namespace Rcpp;
using namespace R;
using namespace arma;

extern "C" SEXP RxODE_ode_dosing();
extern "C" double getLinDeriv(int ncmt, int diff1, int diff2, double rate, double tinf, double Dose, double ka, double tlag, double T, double tT, mat par);
extern "C" double solveLinB(double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);


int locateDoseIndex(mat dosing, const double obs_time){
  // Uses bisection for slightly faster lookup of dose index.
  int i, j, ij;
  i = 0;
  j = dosing.n_rows - 1;
  if (obs_time <= dosing(i, 0)){
    return i;
  }
  if (obs_time >= dosing(j, 0)){
    return j;
  }
  while(i < j - 1) { /* x[i] <= obs_time <= x[j] */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    if(obs_time < dosing(ij, 0))
      j = ij;
    else
      i = ij;
  }
  return i;
}//subscript of dose

double solveLinB(double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag){
  unsigned int ncmt = 1;
  if (C > 0 && gamma > 0){
    ncmt = 3;
  } else if (B > 0 && beta > 0){
    ncmt = 2;
  } else if (A > 0 && alpha > 0){
    ncmt = 1;
  } else {
    stop("You need to specify at least A(=%f) and alpha (=%f). (@t=%f, d1=%d, d2=%d)", A, alpha, t, diff1, diff2);
  }
  mat par = mat(ncmt, 2);
  par(0,0) = alpha;
  par(0,1) = A;
  if (ncmt >= 2){
    par(1,0) = beta;
    par(1,1) = B;
    if (ncmt >= 3){
      par(2,0) = gamma;
      par(2,1) = C;
    }
  }
  int oral, cmt;
  oral = (ka > 0) ? 1 : 0;
  mat dosing  = as<mat>(RxODE_ode_dosing()); // et$get.dosing()
  double ret = 0;
  unsigned int m = 0, l = 0, i = 0, p = 0;
  double thisT = 0.0, tT = 0.0, res, t1, t2, tinf, rate;
  m = locateDoseIndex(dosing, t);
  for(l=0; l <= m; l++){
    //superpostion
    cmt = (((int)(dosing(l, 1)))%10000)/100 - 1;
    if (cmt != linCmt) continue;
    if (dosing(l, 1) > 10000) {
      if (dosing(l, 2) > 0){
	// During infusion
	tT = t - dosing(l, 0);
	thisT = tT - tlag;
	p = l+1;
	while (p < dosing.n_rows && dosing(p, 2) != -dosing(l, 2)){
	  p++;
	}
	if (dosing(p, 2) != -dosing(l, 2)){
	  stop("Could not find a stop to the infusion.  Check the event table.");
	}
	tinf  = dosing(p, 0)-dosing(l, 0);
	rate  = dosing(l, 2);
      } else {
	// After  infusion
	tT = t - dosing(l - 1, 0);
	thisT = tT- tlag;
	
	p = l-1;
	while (p >= 0 && dosing(p, 2) != -dosing(l, 2)){
	  p--;
	}
	if (dosing(p, 2) != -dosing(l, 2)){
	  stop("Could not find a start to the infusion.  Check the event table.");
	}
	tinf  = dosing(l, 0) - dosing(p, 0) - tlag;
	rate  = -dosing(l, 2);
      }
      t1 = ((thisT < tinf) ? thisT : tinf);        //during infusion
      t2 = ((thisT > tinf) ? thisT - tinf : 0.0);  // after infusion
      for (i = 0; i < ncmt; i++){
	if (diff1 == 0 && diff2 == 0){
	  ret += rate*par(i,1) / par(i,0) * (1 - exp(-par(i,0) * t1)) * exp(-par(i,0) * t2);
	} else {
	  ret += getLinDeriv(ncmt, diff1, diff2, rate, tinf, 0.0, 0.0, 0.0, thisT, 0.0, par);
	}
      }
    } else {
      tT = t - dosing(l, 0);
      thisT = tT - tlag;
      if (thisT < 0) continue;
      res = ((oral == 1) ? exp(-ka * thisT) : 0.0);
      for (i = 0; i < ncmt; i++){
	if (diff1 == 0 && diff2 ==0){
	  ret += dosing(l, 2) * par(i,1) *(exp(-par(i,0) * thisT) - res);
	} else {
	  ret += getLinDeriv(ncmt, diff1, diff2, 0.0, 0.0, dosing(l, 2), ka, tlag, thisT, tT, par);
	}
      } //i
    }
  } //l
  return ret;
}
