#define ARMA_DONT_USE_OPENMP // Known to cause speed problems
#include <iostream>
#include <RcppArmadillo.h>
#include <algorithm>
#include "../inst/include/RxODE.h"
#define ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_USE_OPENMP // Known to cause speed problems
using namespace Rcpp;

std::string symengineRes(std::string val){
  if (val == "e" ||
      val == "E" ||
      val == "EulerGamma" ||
      val == "Catalan" ||
      val == "GoldenRatio" ||
      val == "I"){
    return "rx_SymPy_Res_" + val;
  }
  return val;
}

// Create R source for creating a Inductive linearization matrix
// Assume .states=the states in the model
// Assume .env= symengine environment
//[[Rcpp::export]]
std::string rxIndLin_(CharacterVector states){
  std::string ret = "matrix(c(";
  std::string n = "c(";
  for (int i = 0; i < states.size(); i++){
    ret += ".rxIndLinLine(.env$rx__d_dt_"+symengineRes(as<std::string>(states[i]))+ "__" + ",.states),";
    n += "\"" + states[i] +"\",";
  }
  ret += "NULL)," + std::to_string(states.size()) + "," + std::to_string(states.size()+1) +
    ",TRUE,list(" + n +"NULL)," + n + "\"_rxF\")))";
  return ret;
}

bool expm_assign=false;
SEXP expm_s;

//' Inductive linearization solver
//'
//' @param cSub = Current subject number
//' @param neq - Number of equations
//' @param tp - Prior time point/time zeor
//' @param yp - Prior state;  vector size = neq; Final state is updated here
//' @param tf - Final Time
//' @param InfusionRate = Rates of each comparment;  vector size = neq
//' @param on Indicator for if the compartment is "on"
//' @param rtol - rtol based on cmt#; vector size = neq
//' @param atol - atol based on cmt#
//' @param maxsteps Maximum number of steps
//' @param doIndLin Integer to say if inductive linearization is needed.
//' @param locf Do LOCF interpolation for covariates
//' @param delta The delta added to the matrix to make it invertible
//' @param ME the RxODE matrix exponential function
//' @param IndF The RxODE Inductive Linearization function F
//' 
//' @return Returns a status for solving
//' 
//'   1 = Successful solve
//' 
//'   -1 = Maximum number of iterations reached when doing
//'        inductive linearization
//' 
//'   -2 = Maximum number of iterations reached when trying to
//'        make the matrix invertable.
extern "C" int indLin(int cSub, int neq, double tp, double *yp_, double tf,
		      double *InfusionRate_, int *on_, double *rtol, double *atol,
		      int maxsteps, int doIndLin, int locf,
		      int perterbMatrix,
		      double delta, double *rwork,
		      t_ME ME, t_IndF IndF){
  std::ostream nullstream(0);
  arma::set_cerr_stream(nullstream);
  double *ptr = &rwork[0];
  arma::mat m0(ptr, neq, neq, false, false);
  ptr += neq*neq;
  double tcov = tf;
  if (locf) tcov = tp;
  // LOCF=tp; If NOCB tp should be tf
  ME(cSub, tcov, m0.memptr()); // Calculate the initial A matrix based on current time/parameters
  if (!doIndLin){
    // Total possible enhanced matrix is (neq+neq)x(neq+neq)
    // Total possible initial value is (neq+neq)
    // expAt is (neq+neq)x(neq+neq)
    // Total possible output is (neq+neq)
    // =4*neq + 8*neq^2
    // These are simple linear with no f
    // Hence there is no need for matrix inversion
    const arma::vec InfusionRate(InfusionRate_, neq, false, false);
    const arma::vec yp(yp_, neq, false, false);
    // arma::mat inMat;
    // arma::mat mexp;
    // arma::mat ypout;
    unsigned int i, nInf=0;
    arma::vec ypExtra(ptr, neq, false, false);
    ptr += neq;
    arma::mat m0extra(ptr, neq, neq, false, false);
    ptr += neq*neq;
    m0extra.zeros();
    // arma::mat mout;
    for (i = 0; i < (unsigned int)neq; i++){
      if (InfusionRate[i] != 0.0){
	nInf++;
	m0extra[neq*(nInf-1)+i]=1;
	ypExtra[i] = InfusionRate[i];
      }
    }
    if (nInf == 0){
      arma::mat expAT(ptr, neq, neq, false, false);
      ptr += neq*neq;
      expAT = arma::expmat(m0*(tf-tp));
      arma::vec meSol(ptr, neq, false, false);
      ptr += neq;
      meSol = expAT*yp;
      std::copy(meSol.begin(), meSol.end(), yp_);
      return 1;
      // mout = m0;
      // ypout=yp;
    } else {
      arma::mat mout(ptr, neq+nInf, neq+nInf, false, false);
      ptr += (neq+nInf)*(neq+nInf);
      mout.zeros();
      for (int j = neq; j--;){
	std::copy(m0.colptr(j), m0.colptr(j)+neq, mout.colptr(j));
      }
      for (int j = nInf; j--;){
	std::copy(m0extra.colptr(j),m0extra.colptr(j)+neq, mout.colptr(neq+j));
      }
      arma::vec ypout(ptr, neq+nInf, false, false);
      ptr += (neq+nInf);
      arma::mat expAT(ptr, neq+nInf, neq+nInf, false, false);
      ptr += (neq+nInf)*(neq+nInf);
      expAT = arma::expmat(mout*(tf-tp));
      arma::vec meSol(ptr, neq+nInf, false, false);
      std::copy(meSol.begin(), meSol.begin()+neq, yp_);
      return 1;
    }
  } else {
    // In this case the inital matrix should not be expanded. The
    // infusions are put into the F function
    //
    // invA = neq*neq
    // m0 = neq*neq
    // E = neq*neq
    // expAT = neq*neq;
    // facM = neq*neq
    
    // extraF = neq
    // expATy0 = neq
    // f = neq
    // neq*3+ 5*neq^2
    arma::mat invA(ptr, neq, neq, false, false);
    ptr += neq*neq;
    bool canInvert = inv(invA, m0);
    arma::vec extraF(ptr, neq, false, false);
    extraF.zeros();
    ptr += neq;
    arma::mat E(ptr, neq, neq, false, false);
    ptr += neq*neq;
    E.eye();
    int invCount=0;
    while (!canInvert){
      // Add to the diagonal until you can invert
      //
      // At the same remove from the extraF so you can still use the
      // inductive linearization approach.
      m0 = m0 + delta*E;
      extraF = extraF-delta;
      canInvert = inv(invA, m0);
      invCount++;
      if (invCount > maxsteps){
	std::fill_n(&yp_[0], neq, NA_REAL);
	return -2;
      }
    }
    arma::mat expAT(ptr, neq, neq, false, false);
    ptr += neq*neq;
    expAT = arma::expmat(m0*(tf-tp));
    const arma::vec y0(yp_, neq);
    arma::mat expATy0(ptr, neq, neq, false, false);
    ptr+= neq;
    expATy0 = expAT*y0;
    arma::mat facM(ptr, neq, neq, false, false);
    ptr += neq*neq;
    facM = (expAT-E)*invA;
    arma::vec f(ptr, neq, false, false);
    ptr+= neq;
    double *fptr = f.memptr();
    // For LOCF tp for NOCB tf
    IndF(cSub, tcov, tf, fptr, yp_, InfusionRate_, extraF.memptr());
    arma::vec yLast(ptr, neq, false, false);
    ptr += neq;
    yLast = expATy0+facM*f;
    IndF(cSub, tcov, tf, fptr, yLast.memptr(), InfusionRate_, extraF.memptr());
    arma::vec yCur(ptr, neq, false, false);
    ptr += neq;
    yCur = expATy0+facM*f;
    bool converge=false;
    for (int i = 0; i < maxsteps; ++i){
      converge=true;
      for (int j=neq;j--;){
	if (fabs(yCur[j]-yLast[j]) >= rtol[j]*fabs(yCur[j])+atol[j]){
	  converge = false;
	  break;
	}
      }
      if (converge){
	break;
      }
      if (i % perterbMatrix == 0){
	// Reset matricies
	m0 = m0 + delta*E;
	extraF = extraF-delta;
	canInvert = inv(invA, m0);
	while (!canInvert){
	  m0 = m0 + delta*E;
	  extraF = extraF-delta;
	  canInvert = inv(invA, m0);
	  invCount++;
	  if (invCount > maxsteps){
	    std::fill_n(&yp_[0], neq, NA_REAL);
	    return -2;
	  }
	}
	expAT=arma::expmat(m0*(tf-tp));
	expATy0 = expAT*y0;
	facM = (expAT-E)*invA;
      }
      yLast = yCur;
      IndF(cSub, tcov, tf, fptr, yLast.memptr(), InfusionRate_, extraF.memptr());
      yCur = expATy0+facM*f;
    }
    if (!converge){
      std::copy(yCur.begin(), yCur.end(), &yp_[0]);
      // std::fill_n(&yp_[0], neq, NA_REAL);
      return 1;
    } else {
      std::copy(yCur.begin(), yCur.end(), &yp_[0]);
      return 1;
    }
  }
  return 1;
}
