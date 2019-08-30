#include <RcppArmadillo.h>
#include <algorithm>
#include "../inst/include/RxODE.h"
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
		   const arma::vec& yp, const arma::ivec& on,
		   arma::mat &mout, arma::mat &ypout){
  // FIXME on is not used right now.
  unsigned int nrow = m0.n_rows;
  if (nrow != m0.n_cols)
    stop("m0 needs to be a square matrix");
  if (yp.n_elem != nrow)
    stop("yp needs to be the same dimension as m0.");
  if (InfusionRate.n_elem != nrow)
    stop("InfusionRate needs to be the same dimension as m0.");
  unsigned int i, nInf=0;
  arma::vec ypExtra(nrow);
  arma::mat m0extra(nrow, nrow);
  for (i = 0; i < nrow; i++){
    if (InfusionRate[i] != 0.0){
      nInf++;
      m0extra.resize(nrow*nInf,0);
      m0extra[nrow*(nInf-1)+i]=1;
      ypExtra[i] = InfusionRate[i];
    }
  }
  if (nInf == 0){
    mout = m0;
    ypout=yp;
    return;
  }
  mout = join_cols(join_rows(m0, m0extra.cols(0,nInf)),
		   arma::mat(nInf, nInf+nrow, arma::fill::zeros));
  ypout= join_cols(yp, ypExtra.head(nInf));
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
  enhanceMatrix(m0, InfusionRate, yp, on, mout, ypout);
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
		 std::string method="Higham08.b"){
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
//' @param rtol - rtol based on cmt#; vector size = neq
//' @param atol - atol based on cmt#
//' 
//' @return Returns a status for solving
extern "C" int indLin(int cSub, int neq, double tp, double *yp_, double tf,
		      double *InfusionRate_, int *on_, double *rtol, double *atol,
		      int maxsteps, int doIndLin, int locf, t_ME ME, 
		      t_IndF IndF){
  arma::mat m0(neq, neq);
  double tcov = tf;
  if (locf) tcov = tp;
  // For now this is LOCF; If NOCB tp should be tf
  ME(cSub, tcov, m0.memptr()); // Calculate the initial A matrix based on current time/parameters
  if (doIndLin){
    // These are simple linear with no f
    // Hence there is no need for the 
    const arma::vec InfusionRate(InfusionRate_, neq, false, true);
    const arma::vec yp(yp_, neq);
    const arma::ivec on(on_, neq, false, true);
    arma::mat inMat;
    arma::mat mexp;
    arma::mat ypout;
    enhanceMatrix(m0, InfusionRate, yp, on, inMat, ypout);
    arma::mat expAT = rxExpm(inMat, tf-tp);
    arma::vec meSol = expAT*yp;
    std::copy(meSol.begin(), meSol.end(), yp_);
    return 0;
  } else {
    // In this case the inital matrix should not be expanded. The
    // infusions are put into the F function
    arma::mat invA;
    bool canInvert = inv(invA, m0);
    arma::vec extraF(neq, arma::fill::zeros);
    arma::mat E(invA.n_rows, invA.n_rows, arma::fill::eye);
    int invCount=0;
    while (!canInvert){
      // Add to the diagonal until you can invert
      //
      // At the same remove from the extraF so you can use the
      // inductive linearization approach while adjusting the A matrix
      // to be non-singular.
      //
      // FIXME max number of tries
      // FIXME change 0.1
      m0 = m0 + 0.1*E;
      extraF = extraF-0.1;
      canInvert = inv(invA, m0);
      invCount++;
      if (invCount > maxsteps){
	std::fill_n(&yp_[0], neq, NA_REAL);
	return 2;
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
      yLast = yCur;
      IndF(cSub, tcov, tf, fptr, yLast.memptr(), InfusionRate_, extraF.memptr());
      yCur = expATy0+facM*f;
    }
    if (!converge){
      std::fill_n(&yp_[0], neq, NA_REAL);
      return 1;
    } else {
      std::copy(yCur.begin(), yCur.end(), &yp_[0]);
      return 0;
    }
  }
  return 0;
}
