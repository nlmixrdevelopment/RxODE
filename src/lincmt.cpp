// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include "RxODE_types.h"

using namespace Rcpp;
using namespace R;
using namespace arma;


// [[Rcpp::export]]
void getMacroConstants(SEXP rho){
  // Get the parameterizations based on what is present in the environment rho.
  Environment e = as<Environment>(rho);
  int parameterization = 1;
  int ncmt = 1;
  int oral = 0;
  double k = 0, volume = 0, alpha = 0, A = 0, ka = 0,
    k12 = 0, k21 = 0, beta = 0, B = 0, k13 = 0, k31 = 0,
    a0 = 0, a1 = 0, a2 = 0, p = 0, q = 0, r1 = 0, r2 = 0,
    theta = 0, gamma = 0, C = 0, aob = 0;
  if (e.exists("KA")){
    oral = 1;
    ka   = e["KA"];
  }
  if (e.exists("CL")){
    volume  = as<double>(e["V"]);
    parameterization = 1;
    k  = as<double>(e["CL"]) / volume;  // ke = CL/V
    if (e.exists("V2")){
      ncmt = 2;
      k12     = as<double>(e["Q"]) / volume;               // k12 = Q/V
      k21     = as<double>(e["Q"]) / as<double>(e["V2"]);  // k21 = Q/V2
    }
    if (e.exists("VSS")){
      ncmt = 2;
      k12     = as<double>(e["Q"]) / volume;               // k12 = Q/V
      k21     = as<double>(e["Q"]) / (as<double>(e["VSS"])-volume);  // k21 = Q/V2; V2=VSS-V
      parameterization = 3;
    }
    if (e.exists("V3")){
      ncmt = 3;
      //parameterization: CL V Q V2 Q2 V3
      k13     = as<double>(e["Q2"]) / volume;              // k12 = Q2/V
      k31     = as<double>(e["Q2"]) / as<double>(e["V3"]); // k21 = Q2/V3
    }
  } else if (e.exists("K")){
    volume  = as<double>(e["V"]);
    parameterization = 2;
    k  = as<double>(e["K"]);
    if (e.exists("K12")){
      ncmt =2;
      k12 = as<double>(e["K12"]);
      k21 = as<double>(e["K21"]);
    }
    if (e.exists("K13")){
      ncmt = 3;
      k13 = as<double>(e["K13"]);
      k31 = as<double>(e["K31"]);
    }
  } else if (e.exists("AOB")){
    ncmt = 2;
    parameterization = 4;
    volume  = as<double>(e["V"]);
    aob = as<double>(e["AOB"]);
    alpha = as<double>(e["ALPHA"]);
    beta = as<double>(e["BETA"]);
    k21 = (aob*beta+alpha)/(aob+1);
    k   = (alpha*beta)/k21;
    k12 = alpha+beta-k21-k;
  } else if (e.exists("ALPHA") && e.exists("BETA") && e.exists("K21")){
    ncmt = 2;
    parameterization = 5;
    volume  = as<double>(e["V"]);
    k21 = as<double>(e["K21"]);
    alpha = as<double>(e["ALPHA"]);
    beta = as<double>(e["BETA"]);
    k    = (alpha*beta)/k21;
    k12 = alpha+beta-k21-k;
  } else {
    stop("Cannot figure out the parameterization for the linear compartments.");
  }
  mat g = mat(ncmt, 2);
  if (ncmt == 1){
    alpha = k;
    A = 1.0/volume;
    g(0, 0) = alpha;
    if (oral == 1){
      g(0, 1) = ka / (ka - alpha) * A;
    } else {
      g(0, 1) = A;
    }
    if (parameterization == 1){
      k  = as<double>(e["CL"]) / volume;  // ke = CL/V or CL=ke*V
    }
    e["alpha"] = alpha;
    e["A"] = A;
  } else if (ncmt == 2) {
    beta  = 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
    alpha = k21 * k / beta;
    
    A     = (alpha - k21) / (alpha - beta) / volume;
    B     = (beta - k21) / (beta - alpha) / volume;
    
    g(0,0) = alpha;
    g(1,0) = beta;
    if (oral==1) {
      g(0,1) = ka / (ka - alpha) * A;
      g(1,1) = ka / (ka - beta) * B;
    } else { // IV
      g(0,1) = A;
      g(1,1) = B;
    }
    e["alpha"] = alpha;
    e["A"] = A;
    e["beta"] = beta;
    e["B"] = B;
  } else if (ncmt == 3){
    a0      = k * k21 * k31;
    a1      = k * k31 + k21 * k31 + k21 * k13 + k * k21 + k31 * k12;
    a2      = k + k12 + k13 + k21 + k31;

    p       = a1 - a2 * a2 / 3.0;
    q       = 2.0 * a2 * a2 * a2 / 27.0 - a1 * a2 /3.0 + a0;

    r1      = sqrt(-p * p * p / 27.0);
    r2      = 2 * pow(r1 , 1.0 / 3.0);

    theta   = acos(-q / (2.0 * r1)) / 3.0;

    alpha   = -(cos(theta) * r2 - a2 / 3.0);
    beta    = -(cos(theta + 2.0 / 3.0 * M_PI) * r2 - a2 / 3.0);
    gamma   = -(cos(theta + 4.0 / 3.0 * M_PI) * r2 - a2 / 3.0);

    A       = (k21 - alpha) * (k31 - alpha) / (alpha - beta) / (alpha - gamma) / volume;
    B       = (k21 - beta) * (k31 - beta) / (beta - alpha) / (beta - gamma) / volume;
    C       = (k21 - gamma) * (k31 - gamma) / (gamma - alpha) / (gamma - beta) / volume;

    g(0,0) = alpha;
    g(1,0) = beta;
    g(2,0) = gamma;

    if (oral==1) {
      g(0,1) = ka / (ka - alpha) * A;
      g(1,1) = ka / (ka - beta) * B;
      g(2,1) = ka / (ka - gamma) * C;
    } else {
      g(0,1) = A;
      g(1,1) = B;
      g(2,1) = C;
    }

    e["alpha"] = alpha;
    e["A"] = A;
    e["beta"] = beta;
    e["B"] = B;
    e["gamma"] = gamma;
    e["C"] = C;
  }
  e["parameterization"] = parameterization;
  e["ncmt"] = ncmt;
  e["oral"] = oral;
  e["g"] = g;
}

// Unlinke nlmixr's solved equtions, this uses an event table.
void getLinDerivs(SEXP rho);

// [[Rcpp::export]]
void linCmtEnv(SEXP rho){
  getMacroConstants(rho);
  getLinDerivs(rho);
  Environment e = as<Environment>(rho);
  int linCmt = as<int>(e["cmt"]);
  int parameterization = as<int>(e["parameterization"]);
  unsigned int ncmt = as<unsigned int>(e["ncmt"]);
  int oral = as<int>(e["oral"]);
  mat g = as<mat>(e["g"]);
  mat et = as<mat>(e["events"]); // et2$get.EventTable()
  mat dosing = as<mat>(e["dosing"]); // et$get.dosing()
  mat sampling = as<mat>(e["sampling"]); // et$get.sampling()
  mat par = as<mat>(e["g"]);
  int cmt = -1;
  double tlag = 0.0;
  if (e.exists("tlag")){
    tlag = as<double>(e["tlag"]);
  }
  double ka = 0;
  if (oral == 1){
    ka = as<double>(e["KA"]);
  }
  mat ret(sampling.n_rows,
	  2*ncmt+oral+((tlag == 0.0) ? 0 : 1) + 1,
	  fill::zeros);
  unsigned int b = 0, m = 0, l = 0, i = 0, j = 0;
  double thisT = 0.0, tT = 0.0, sum = 0.0, res, t1, t2, tinf, rate;
  mat cur;
  
  for (b = 0; b < sampling.n_rows; b++){
    // Get the next dose, index if necessary.
    while (sampling(b,0) > dosing(m,0) &&
	   m < dosing.n_rows){
      m++;
    }
    for(l=0; l < m; l++){
      //superpostion
      cmt = (((int)(dosing(l, 1)))%10000)/100 - 1;
      if (cmt != linCmt) continue;
      if (dosing(l, 1) > 10000) {
	sum = 0.0;
        if (dosing(l, 2) > 0){
          // During infusion
	  tT = sampling(b,0) - dosing(l, 0);
	  thisT = tT - tlag;
	  tinf  = dosing(l+1, 0)-dosing(l, 0);
	  rate  = dosing(l, 2);
	} else {
	  // After  infusion
	  tT = sampling(b,0) - dosing(l - 1, 0);
	  thisT = tT- tlag;
	  tinf  = dosing(l, 0) - dosing(l-1, 0) - tlag;
	  rate  = -dosing(l, 2);
	  if (dosing(l,2) != -dosing(l-1,2)){
	    stop("This can only handle constant on/off infusion rates...(currently)");
	  }
	}
	t1 = ((thisT < tinf) ? thisT : tinf);        //during infusion
	t2 = ((thisT > tinf) ? thisT - tinf : 0.0);  // after infusion
	for (i = 0; i < ncmt; i++)
	  sum += par(i,1) / par(i,0) * (1 - exp(-par(i,0) * t1)) * exp(-par(i,0) * t2);
	ret(b,0) += rate * sum;
      } else {
        thisT = sampling(b,0) - dosing(l, 0) - tlag;
	if (thisT < 0) continue;
	res = ((oral == 1) ? exp(-ka * thisT) : 0.0);
        for (i = 0; i < ncmt; i++){
	  ret(b, 0) += dosing(l, 2) * par(i,1) * (exp(-par(i,0) * thisT) - res);
	  // Now add the derivs
	  for (j = 0; j < ncmt*2; j++){
	    if (parameterization == 1){
	      switch(j){
	      case 0:
		cur = as<mat>(e["dV"]);
		break;
	      case 1:
                cur = as<mat>(e["dCL"]);
		break;
	      case 2:
                cur = as<mat>(e["dV2"]);
		break;
	      case 3:
                cur = as<mat>(e["dQ"]);
		break;
	      case 4:
                cur = as<mat>(e["dV3"]);
		break;
	      case 5:
                cur = as<mat>(e["dQ2"]);
		break;
	      }
	    } else if (parameterization == 2){
	      switch(j){
	      case 0:
                cur = as<mat>(e["dV"]);
                break;
	      case 1:
                cur = as<mat>(e["dK"]);
                break;
	      case 2:
                cur = as<mat>(e["dK12"]);
                break;
	      case 3:
                cur = as<mat>(e["dK21"]);
                break;
	      case 4:
                cur = as<mat>(e["dK13"]);
                break;
	      case 5:
                cur = as<mat>(e["dK31"]);
                break;
	      }
	    } else if (parameterization == 3){
	      switch(j){
              case 0:
                cur = as<mat>(e["dV"]);
                break;
              case 1:
                cur = as<mat>(e["dCL"]);
                break;
              case 2:
                cur = as<mat>(e["dVSS"]);
                break;
              case 3:
                cur = as<mat>(e["dQ"]);
                break;
              }
	    }
	    ret(b,j+1) = -dosing(l, 2)*(thisT*par(i,1)*cur(i, 0) + (res*exp(thisT*par(i, 0)) - 1)*cur(i, 1))*exp(-thisT*par(i, 0));
          } //j
	  if (oral == 1){ // dKA
	    cur = as<mat>(e["dKA"]);
	    // rxSymPy("diff(Dose*p1(V)*(exp(-p0(V)*thisT)-res), V)")
            // [1] "-Dose*thisT*p1(V)*exp(-thisT*p0(V))*Derivative(p0(V), V) + Dose*(-res + exp(-thisT*p0(V)))*Derivative(p1(V), V)"
	    ret(b,j) = dosing(l, 2)*(-thisT*exp(-thisT*par(i, 0))*cur(i, 0) + thisT*exp(-ka*thisT))*par(i, 1) +
	      dosing(l, 2)*(exp(-thisT*par(i, 0)) - exp(-ka*thisT))*cur(i, 1);
	    j++;
	  }
	  // dTlag
	  if (tlag > 0){
	    if (oral == 1){
	      // rxSymPy("diff(Dose*p1*(exp(-p0*(tT-tlag))-res), tlag)")
              // [1] "Dose*p1*(-ka*exp(-ka*(tT - tlag)) + p0*exp(-p0*(tT - tlag)))"
	      ret(b,j) = dosing(l, 2)*par(i, 1)*((-ka*(tT-tlag)) + par(i, 0)*exp(-par(i, 0)*(tT - tlag)));
            } else {
	      // rxSymPyExec("res=0");rxSymPy("diff(Dose*p1*(exp(-p0*(tT-tlag))-res), tlag)")
              // [1] "Dose*p0*p1*exp(-p0*(tT - tlag))"
	      ret(b,j) = dosing(l, 2)*par(i, 1)*par(i, 0)*exp(-par(i, 0)*(tT - tlag));
	    }
	  }
        } //i
      }
    } //l
  } // b
  e["ret"] = ret;
}
