// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include "RxODE_types.h"

using namespace Rcpp;
using namespace R;
using namespace arma;

extern "C" SEXP RxODE_ode_dosing();

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

int locateDoseIndex(mat dosing, const double obs_time){
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


// Unlinke nlmixr's solved equtions, this uses an event table.
void getLinDerivs(SEXP rho);

// [[Rcpp::export]]
void linCmtEnv(SEXP rho){
  getMacroConstants(rho);
  getLinDerivs(rho);
  Environment e = as<Environment>(rho);
  int linCmt = as<int>(e["cmt"])-1;
  int parameterization = as<int>(e["parameterization"]);
  unsigned int ncmt = as<unsigned int>(e["ncmt"]);
  int oral = as<int>(e["oral"]);
  mat g = as<mat>(e["g"]);
  mat dosing = as<mat>(e["dosing"]); // et$get.dosing()
  mat sampling = as<mat>(e["t"]); // et$get.sampling()[,]
  mat par = as<mat>(e["g"]);
  int cmt = -1;
  double tlag = 0.0;
  if (e.exists("TLAG")){
    tlag = as<double>(e["TLAG"]);
  }
  double ka = 0;
  if (oral == 1){
    ka = as<double>(e["KA"]);
  }
  mat ret(sampling.n_rows,
	  2*ncmt+oral+((tlag == 0.0) ? 0 : 1) + 2,
	  fill::zeros);
  StringVector cname(2*ncmt+oral+((tlag == 0.0) ? 0 : 1) + 2);
  cname(0) = "t";
  cname(1) = "f";
  unsigned int b = 0, m = 0, l = 0, i = 0, j = 0, p = 0;
  double thisT = 0.0, tT = 0.0, res, t1, t2, tinf, rate;
  mat cur;
  for (b = 0; b < sampling.n_rows; b++){
    m = locateDoseIndex(dosing, sampling(b,0));
    ret(b,0)=sampling(b,0);
    for(l=0; l <= m; l++){
      //superpostion
      cmt = (((int)(dosing(l, 1)))%10000)/100 - 1;
      if (cmt != linCmt) continue;
      if (dosing(l, 1) > 10000) {
        if (dosing(l, 2) > 0){
          // During infusion
	  tT = sampling(b,0) - dosing(l, 0);
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
	  tT = sampling(b,0) - dosing(l - 1, 0);
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
	  ret(b,1) += rate*par(i,1) / par(i,0) * (1 - exp(-par(i,0) * t1)) * exp(-par(i,0) * t2);
	  // FIXME derivs
	  for (j = 0; j < ncmt*2; j++){
            if (parameterization == 1){
              switch(j){
              case 0:
                cur = as<mat>(e["dV"]);
		cname(2+j) = "dV";
                break;
              case 1:
                cur = as<mat>(e["dCL"]);
		cname(2+j) = "dCL";
                break;
              case 2:
                cur = as<mat>(e["dV2"]);
		cname(2+j) = "dV2";
                break;
              case 3:
                cur = as<mat>(e["dQ"]);
		cname(2+j) = "dQ";
                break;
              case 4:
                cur = as<mat>(e["dV3"]);
		cname(2+j) = "dV3";
                break;
              case 5:
                cur = as<mat>(e["dQ2"]);
		cname(2+j) = "dQ2";
                break;
              }
            } else if (parameterization == 2){
              switch(j){
              case 0:
                cur = as<mat>(e["dV"]);
		cname(2+j) = "dV";
                break;
              case 1:
                cur = as<mat>(e["dK"]);
		cname(2+j) = "dK";
                break;
              case 2:
                cur = as<mat>(e["dK12"]);
		cname(2+j) = "dK12";
                break;
              case 3:
                cur = as<mat>(e["dK21"]);
		cname(2+j) = "dK21";
                break;
              case 4:
                cur = as<mat>(e["dK13"]);
		cname(2+j) = "dK13";
                break;
              case 5:
                cur = as<mat>(e["dK31"]);
		cname(2+j) = "dK31";
                break;
              }
            } else if (parameterization == 3){
              switch(j){
              case 0:
                cur = as<mat>(e["dV"]);
		cname(2+j) = "dV";
                break;
              case 1:
                cur = as<mat>(e["dCL"]);
		cname(2+j) = "dCL";
                break;
              case 2:
                cur = as<mat>(e["dVSS"]);
		cname(2+j) = "dVSS";
                break;
              case 3:
                cur = as<mat>(e["dQ"]);
		cname(2+j) = "dQ";
                break;
              }
            }
	    // > rxSymPy("simplify(diff(rate*p1(V)/p0(V)*(1-exp(-p0(V)*t1))*exp(p0(V)*t2),V))")
            // [1] "rate*(-(exp(t1*p0(V)) - 1)*p1(V)*Derivative(p0(V), V) + (t1*p1(V)*Derivative(p0(V), V) + t2*(exp(t1*p0(V)) - 1)*p1(V)*Derivative(p0(V), V) + (exp(t1*p0(V)) - 1)*Derivative(p1(V), V))*p0(V))*exp(-2*t1*p0(V))*exp((t1 + t2)*p0(V))/p0(V)**2"
            ret(b,j+2) = rate*(-(exp(t1*par(i, 0)) - 1)*par(i, 1)*cur(i, 0) + (t1*par(i, 1)*cur(i, 0) + t2*(exp(t1*par(i, 0)) - 1)*par(i, 1)*cur(i, 0) + (exp(t1*par(i, 0)) - 1)*cur(i, 1))*par(i, 0))*exp(-2*t1*par(i, 0))*exp((t1 + t2)*par(i, 0))/pow(par(i, 0), 2);
          }
	  if (oral == 1){
	    ret(b, j+2) = 0;
	    cname(j+2) = "dKA";
	    j++;
	  }
	  // FIXME tlag
        }
      } else {
        thisT = sampling(b,0) - dosing(l, 0) - tlag;
	if (thisT < 0) continue;
	res = ((oral == 1) ? exp(-ka * thisT) : 0.0);
        for (i = 0; i < ncmt; i++){
	  ret(b, 1) += dosing(l, 2) * par(i,1) *
	    (exp(-par(i,0) * thisT) - res);
	  // Now add the derivs
	  for (j = 0; j < ncmt*2; j++){
	    if (parameterization == 1){
              switch(j){
              case 0:
                cur = as<mat>(e["dV"]);
                cname(2+j) = "dV";
                break;
              case 1:
                cur = as<mat>(e["dCL"]);
                cname(2+j) = "dCL";
                break;
              case 2:
                cur = as<mat>(e["dV2"]);
                cname(2+j) = "dV2";
                break;
              case 3:
                cur = as<mat>(e["dQ"]);
                cname(2+j) = "dQ";
                break;
              case 4:
                cur = as<mat>(e["dV3"]);
                cname(2+j) = "dV3";
                break;
              case 5:
                cur = as<mat>(e["dQ2"]);
                cname(2+j) = "dQ2";
                break;
              }
            } else if (parameterization == 2){
              switch(j){
              case 0:
                cur = as<mat>(e["dV"]);
                cname(2+j) = "dV";
                break;
              case 1:
                cur = as<mat>(e["dK"]);
                cname(2+j) = "dK";
                break;
              case 2:
                cur = as<mat>(e["dK12"]);
                cname(2+j) = "dK12";
                break;
              case 3:
                cur = as<mat>(e["dK21"]);
                cname(2+j) = "dK21";
                break;
              case 4:
                cur = as<mat>(e["dK13"]);
                cname(2+j) = "dK13";
                break;
              case 5:
                cur = as<mat>(e["dK31"]);
                cname(2+j) = "dK31";
                break;
              }
            } else if (parameterization == 3){
              switch(j){
              case 0:
                cur = as<mat>(e["dV"]);
                cname(2+j) = "dV";
                break;
              case 1:
                cur = as<mat>(e["dCL"]);
                cname(2+j) = "dCL";
                break;
              case 2:
                cur = as<mat>(e["dVSS"]);
                cname(2+j) = "dVSS";
                break;
              case 3:
                cur = as<mat>(e["dQ"]);
                cname(2+j) = "dQ";
                break;
              }
            }
	    ret(b,j+2) = -dosing(l, 2)*(thisT*par(i,1)*cur(i, 0) + (res*exp(thisT*par(i, 0)) - 1)*cur(i, 1))*exp(-thisT*par(i, 0));
          } //j
	  if (oral == 1){ // dKA
	    cur = as<mat>(e["dKA"]);
	    // rxSymPy("diff(Dose*p1(V)*(exp(-p0(V)*thisT)-res), V)")
            // [1] "-Dose*thisT*p1(V)*exp(-thisT*p0(V))*Derivative(p0(V), V) + Dose*(-res + exp(-thisT*p0(V)))*Derivative(p1(V), V)"
	    ret(b,j+2) = dosing(l, 2)*(-thisT*exp(-thisT*par(i, 0))*cur(i, 0) + thisT*exp(-ka*thisT))*par(i, 1) +
	      dosing(l, 2)*(exp(-thisT*par(i, 0)) - exp(-ka*thisT))*cur(i, 1);
	    cname(2+j) = "dKA";
	    j++;
	  }
	  // dTlag
	  if (tlag > 0){
	    cname(2+j) = "dTLAG";
	    if (oral == 1){
	      // rxSymPy("diff(Dose*p1*(exp(-p0*(tT-tlag))-res), tlag)")
              // [1] "Dose*p1*(-ka*exp(-ka*(tT - tlag)) + p0*exp(-p0*(tT - tlag)))"
	      ret(b,j+2) = dosing(l, 2)*par(i, 1)*((-ka*(tT-tlag)) + par(i, 0)*exp(-par(i, 0)*(tT - tlag)));
            } else {
	      // rxSymPyExec("res=0");rxSymPy("diff(Dose*p1*(exp(-p0*(tT-tlag))-res), tlag)")
              // [1] "Dose*p0*p1*exp(-p0*(tT - tlag))"
	      ret(b,j+2) = dosing(l, 2)*par(i, 1)*par(i, 0)*exp(-par(i, 0)*(tT - tlag));
	    }
	  }
        } //i
      }
    } //l
  } // b
  NumericMatrix retn = wrap(ret);
  colnames(retn) = cname;
  e["ret"] = retn;
}

extern "C" double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);

Environment solvedEnv;

void clearSolved(){
  const char* strs[] = { "t", "cmt", "V", "CL",
			 "V2", "Q", "V3", "Q2",
			 "K","K12","K21", "K13",
			 "K31", "VSS",
			 "KA","TLAG"};
  for (unsigned int i=0; i<sizeof(strs)/sizeof(const char*); i++) {
    const char *name = strs[i];
    if (solvedEnv.exists(name))
      solvedEnv.remove(name);
  }
}

double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8){
  double ret;
  Environment solvedEnv;
  clearSolved();
  mat tmat = mat(1,1);
  tmat(0,0) = t;
  solvedEnv["t"]=tmat;
  solvedEnv["cmt"]=cmt;
  if (parameterization == 1){
    if (p1 > 0 && p2 > 0){
      solvedEnv["V"] = p1;
      solvedEnv["CL"] = p2;
    } else {
      stop("Clearance and Volume must be above 0.");
    }
    if (p3 > 0 && p4 > 0){
      solvedEnv["V2"] = p3;
      solvedEnv["Q"] = p4;
    } else if (!(p3 <= 0 && p4 <= 0)){
      stop("Both V2 and Q have to be above 0.");
    }
    if (p5 > 0 && p6 > 0){
      solvedEnv["V3"] = p5;
      solvedEnv["Q2"] = p6;
    } else if (!(p5 <= 0 && p6 <= 0)){
      stop("Both V2 and Q have to be above 0.");
    }
  } else if (parameterization == 2){
    if (p1 > 0 && p2 > 0){
      solvedEnv["V"] = p1;
      solvedEnv["K"] = p2;
    } else {
      stop("K and Volume must be above 0.");
    }
    if (p3 > 0 && p4 > 0){
      solvedEnv["K12"] = p3;
      solvedEnv["K21"] = p4;
    } else if (!(p3 <= 0 && p4 <= 0)){
      stop("Both K12 and K21 have to be above 0.");
    }

    if (p5 > 0 && p6 > 0){
      solvedEnv["K13"] = p5;
      solvedEnv["K31"] = p6;
    } else if (!(p5 <= 0 && p6 <= 0)){
      stop("Both K13 and K31 have to be above 0.");
    }
  } else if (parameterization == 3) {
    if (p1 > 0 && p2 > 0){
      solvedEnv["V"] = p1;
      solvedEnv["CL"] = p2;
    } else {
      stop("Clearance and Volume must be above 0.");
    }
    if (p3 > 0 && p4 > 0){
      solvedEnv["VSS"] = p3;
      solvedEnv["Q"] = p4;
    } else if (!(p3 <= 0 && p4 <= 0)){
      stop("Both V2 and Q have to be above 0.");
    }
    if (p5 > 0 && p6 > 0){
      stop("This parametrizaiton is not yet supported for 3 cmt model");
      // solvedEnv["V3"] = p5;
      // solvedEnv["Q2"] = p6;
    } else if (!(p5 <= 0 && p6 <= 0)){
      stop("Both V2 and Q have to be above 0.");
    }
  } else {
    stop("Other parameterzations not yet supported.");
  }
  if (p7 > 0){
    solvedEnv["KA"] = p7;
  }
  if (p8 > 0){
    solvedEnv["TLAG"] = p8;
  }
  solvedEnv["dosing"] = RxODE_ode_dosing();
  linCmtEnv(as<SEXP>(solvedEnv));
  mat retm = as<mat>(solvedEnv["ret"]);
  if (col >= retm.n_cols){
    stop("Requested variable outside of bounds.");
  } else {
    ret = retm(0, col);
  }
  return ret;
}
