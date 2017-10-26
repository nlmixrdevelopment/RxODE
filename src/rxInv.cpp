// [[Rcpp::depends(RcppArmadillo)]]
#include <stdarg.h>
#include <RcppArmadillo.h>
#include <R.h>
using namespace Rcpp;
using namespace R;
using namespace arma;
extern "C" SEXP _rxCholInv(SEXP dms, SEXP theta, SEXP tn);

//' Invert matrix using Rcpp Armadilo.  
//'
//' @param matrix matrix to be inverted.
//' @return inverse or pseudo inverse of matrix.
//' @export
// [[Rcpp::export]]
NumericVector rxInv(SEXP matrix){
  mat smatrix= as<mat>(matrix);
  mat imat;
  try{
    imat = inv(smatrix);
  } catch(...){
    Rprintf("Warning: matrix seems singular; Using pseudo-inverse\n");
    imat = pinv(smatrix);
  }
  NumericVector ret;
  ret = wrap(imat);
  return(ret);
}

arma::mat rxToCholOmega(arma::mat cholMat){
  // Only the cholesky is needed for the liklihood caclation
  return inv(trimatu(cholMat));
}

// [[Rcpp::export]]
arma::mat rxToOmega(arma::mat cholMat){
  // The Omega is need.
  // U^-1*trans(U^1) = Omega
  arma::mat U1 = inv(trimatu(cholMat));
  return U1*trans(U1);
}


// [[Rcpp::export]]
void rxSymInvCholEnvCalculate(Environment e, std::string what, Function invFn){
  if (!e.exists(what)){
    List invObj;
    if (e.exists("invobj")){
      invObj = as<List>(e["invobj"]);
    } else {
      stop("error in rxSymInvCholEnvCalculate environment");
    }
    NumericVector theta;
    if (e.exists("theta")){
       theta = as<NumericVector>(e["theta"]);
    } else {
      stop("theta for omega calculations not setup yet.");
    }
    int ntheta = theta.size(), i=0;
    if (what == "chol.omegaInv"){
      e["chol.omegaInv"]=as<NumericMatrix>(invFn(invObj, theta, "cholOmegaInv"));
    } else if (what == "omegaInv"){
      e["omegaInv"]= as<NumericMatrix>(invFn(invObj, theta, "omegaInv",-1));
    } else if (what == "d.omegaInv"){
      List ret(ntheta);
      for (i =0; i < ntheta; i++){
	ret[i] = as<NumericMatrix>(invFn(invObj, theta, "d(omegaInv)",i+1));
      }
      e["d.omegaInv"] = ret;
    } else if (what == "d.D.omegaInv"){
      List ret(ntheta);
      for (i =0; i < ntheta; i++){
        ret[i] = as<NumericVector>(invFn(invObj, theta, "d(D)",i+1));
      }
      e["d.D.omegaInv"] = ret;
    } else if (what == "chol.omega"){
      rxSymInvCholEnvCalculate(e, "chol.omegaInv", invFn);
      arma::mat ret = rxToCholOmega(as<arma::mat>(e["chol.omegaInv"]));
      e["chol.omega"] = ret; 
    } else if (what == "omega"){
      rxSymInvCholEnvCalculate(e, "chol.omega", invFn);
      arma::mat U1 = as<mat>(e["chol.omega"]);
      arma::mat omega = U1*trans(U1);
      e["omega"] = omega;
    } else if (what == "log.det.OMGAinv.5"){
      // Note this does NOT include the 2 pi bit
      rxSymInvCholEnvCalculate(e,"chol.omegaInv",invFn);
      arma::mat c = as<arma::mat>(e["chol.omegaInv"]);
      arma::vec diag = c.diag();
      arma::vec ldiag = log(diag);
      NumericVector ret = as<NumericVector>(wrap(sum(ldiag)));
      e["log.det.OMGAinv.5"] = ret;
    }
  }
}

//' Calculate Wishart Variance based on Omega matrix
//'
//' @param Omega is the square positive definite matrix
//' @param nu is the degrees of freedom of the Wishart Distribution.
//' @return Variance matrix
//' @export
//[[Rcpp::export]]
arma::mat rxInvWishartVar(arma::mat Omega, double nu){
  if (Omega.n_rows != Omega.n_cols){
    stop("Matrix is required to be positive definite (square) matrix.");
  }
  double p = (double)(Omega.n_rows);
  double vp3 = nu - p - 3;
  if (vp3 <= 0){
    stop("The degrees of freedom are too small.");
  }
  double vpm1 = nu-p-1;
  double vpp1 = nu-p+1;
  double vp   = nu-p;
  double tmp;
  arma::mat var(Omega.n_rows, Omega.n_cols);
  for (int i =0; i < p; i++){
    for (int j=i; j<p; j++){
      tmp =Omega(i, j);
      if (j==i){
	var(i,i)=2*tmp*tmp/(vpm1*vpm1*vp3);
      } else {
	var(i,j) = (vpp1*tmp*tmp+vpm1*tmp*Omega(j,i))/(vp*vpm1*vpm1*vp3);
	var(j,i) = var(i,j);
      }
    }
  }
  return var;
}

// [[Rcpp::export]]
void RxODE_finalize_focei_omega(SEXP rho){
  Environment e = as<Environment>(rho);
  List dOmega = as<List>(e["dOmega"]);
  mat omegaInv = as<mat>(e["omegaInv"]);
  mat c;
  vec diag;
  int ntheta = dOmega.length();
  int neta = omegaInv.n_rows;
  mat cEta = zeros(neta,1);
  NumericVector trInv(ntheta);
  List prod1(ntheta);
  List prod2(ntheta);
  int i,j;
  for (i = 0; i < ntheta; i++){
    c = omegaInv * as<mat>(dOmega[i]);
    diag = c.diag();
    trInv[i] = 0.5*sum(diag);
    c = c * omegaInv;
    prod1[i] = c;
    List prodI(neta);
    for (j = 0; j < neta; j++){
      cEta(j,0) = 1;
      prodI[j] = c*cEta;
      cEta(j,0) = 0;
    }
    prod2[i] = prodI;
  }
  // = tr(omegaInv*dOmega) = tr(omegaInv*omega*d(Omega^-1)*omega) = tr(d(Omega^-1)*omega)
  // This requires an expensive inverse....
  //e["tr.omegaInv.dOmega.0.5"] = trInv;
  // omegaInv.dOmega.omegaInv = d(Omega^-1)
  e["omegaInv.dOmega.omegaInv"] = prod1;
  // omegaInv.dOmega.omegaInv.dEta = d(Omega^-1)*dEta
  e["omegaInv.dOmega.omegaInv.dEta"] = prod2;
}
// [[Rcpp::export]]
NumericVector RxODE_finalize_log_det_OMGAinv_5(SEXP rho){
  // log(det(omegaInv^1/2)) = 1/2*log(det(omegaInv))
  Environment e = as<Environment>(rho);
  NumericVector reset;
  mat c;
  try{
    c = chol(as<mat>(e["omegaInv"]));
  } catch(...){
    e["reset"] = 0;
    c = as<mat>(e["omegaInv"]);
    Function nearPD = as<Function>(e["nearPD"]);
    c = as<mat>(nearPD(c, rho));
    reset = as<NumericVector>(e["reset"]);
    if (reset[0] != 1){
      Rprintf("Warning: The Omega^-1 is non-positive definite, correcting with nearPD\n");
    }
    c = chol(-as<mat>(e["omegaInv"]));
  }
  vec diag = c.diag();
  vec ldiag = log(diag);
  NumericVector ret = as<NumericVector>(wrap(sum(ldiag)));
  e["log.det.OMGAinv.5"] = ret;
  return ret;
}

