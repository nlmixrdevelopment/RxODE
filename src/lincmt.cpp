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
