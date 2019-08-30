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

// expm::expm method="PadeRBS" is faster than rexpokit::expm BUT
// expm::expm method="Higham08" is "better" and has a C interface.  It
// is slower, though.

//' Enhance the Matrix for expm
//' This adds the columns indicating infusions if needed.
//' It should also take away columns depending on if a compartment is "on"
void enhanceMatrix(const arma::mat& m0, const arma::vec& InfusionRate,
		   const arma::vec& yp, const int *on,
		   arma::mat &mout, arma::mat &ypout){
  // FIXME on is not used right now.
  unsigned int nrow = m0.n_rows;
  if (nrow != m0.n_cols)
    stop("m0 needs to be a square matrix");
  if (yp.n_elem != nrow)
    stop("yp needs to be the same dimension as m0.");
  if (InfusionRate.n_elem != nrow)
    stop("InfusionRate needs to be the same dimension as m0.");
  
}

//' Enhance the Matrix for expm
//' This addst the columns indicating infusions if needed.
//' @param m0 Initial matrix
//' @param InfusionRate is a vector of infusion rates
//' @param yp is the last known state concentrations
//' @param on is the vector of on/off compartment states
//'
//' This is mostly for testing
//' @noRd
//[[Rcpp::export]]
List rxExpmMat(const arma::mat& m0, const arma::vec& InfusionRate,
	       const arma::vec& yp, const arma::ivec &on){
  arma::mat mout;
  arma::vec ypout;
  enhanceMatrix(m0, InfusionRate, yp, on.memptr(), mout, ypout);
  List ret(2);
  ret[0] = wrap(mout);
  ret[1] = wrap(ypout);
  return ret;
}

//' Armadillo interface to R package expm
//'
//' @param inMat is the in matrix
//' @param t is the time for the calculation
//'
//' @inheritParams expm::expm
//'
//' @return expm(t*X)
//' @noRd
//[[Rcpp::export]]
arma::mat rxExpm(const arma::mat& inMat, double t = 1,
		 std::string method="PadeRBS"){
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment expmNS = loadNamespace("expm");
  Function expm = expmNS["expm"];
  arma::mat out0 =t*inMat;
  return as<arma::mat>(expm(out0,_["method"]=method));
}

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
		      int maxsteps, int doIndLin, int locf, double delta,
		      t_ME ME, t_IndF IndF){
  std::ostream nullstream(0);
  arma::set_cerr_stream(nullstream);
  arma::mat m0(neq, neq);
  double tcov = tf;
  if (locf) tcov = tp;
  // LOCF=tp; If NOCB tp should be tf
  ME(cSub, tcov, m0.memptr()); // Calculate the initial A matrix based on current time/parameters
  if (!doIndLin){
    // These are simple linear with no f
    // Hence there is no need for matrix inversion
    const arma::vec InfusionRate(InfusionRate_, neq);
    const arma::vec yp(yp_, neq);
    arma::mat inMat;
    arma::mat mexp;
    arma::mat ypout;
    unsigned int i, nInf=0;
    arma::vec ypExtra(neq);
    arma::mat m0extra(neq, neq, arma::fill::zeros);
    arma::mat mout;
    for (i = 0; i < (unsigned int)neq; i++){
      if (InfusionRate[i] != 0.0){
	nInf++;
	m0extra[neq*(nInf-1)+i]=1;
	ypExtra[i] = InfusionRate[i];
      }
    }
    if (nInf == 0){
      mout = m0;
      ypout=yp;
    } else {
      mout = join_cols(join_rows(m0, m0extra.cols(0,nInf-1)),
		       arma::mat(nInf, nInf+neq, arma::fill::zeros));
      ypout= join_cols(yp, ypExtra.head(nInf));
    }
    arma::mat expAT = rxExpm(mout, tf-tp);
    arma::vec meSol = expAT*ypout;
    std::copy(meSol.begin(), meSol.end()-nInf, yp_);
    return 1;
  } else {
    // In this case the inital matrix should not be expanded. The
    // infusions are put into the F function
    arma::mat invA;
    bool canInvert = inv(invA, m0);
    arma::vec extraF(neq, arma::fill::zeros);
    arma::mat E(neq, neq, arma::fill::eye);
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
    arma::mat expAT = rxExpm(m0, tf-tp);
    const arma::vec y0(yp_, neq);
    arma::vec expATy0 = expAT*y0;

    arma::mat facM = (expAT-E)*invA;
    arma::vec f(neq);
    double *fptr = f.memptr();
    // For LOCF tp for NOCB tf
    IndF(cSub, tcov, tf, fptr, yp_, InfusionRate_, extraF.memptr());
    arma::vec yLast = expATy0+facM*f;
    IndF(cSub, tcov, tf, fptr, yLast.memptr(), InfusionRate_, extraF.memptr());
    arma::vec yCur = expATy0+facM*f;
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
      if (i % 100 == 0){
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
	expAT=rxExpm(m0, tf-tp);
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
