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
//' Get Omega^-1 and derivatives
//'
//' @param invObjOrMatrix Object for inverse-type calculations.  If this is a matrix,
//'     setup the object for inversion by \code{\link{rxSymInvCholCreate}} with the default arguments and return
//'     a reactive s3 object.  Otherwise, use the inversion object to calculate the requested derivative/inverse.
//' @param theta Thetas to be used for calculation.  If missing (\code{NULL}), a
//'     special s3 class is created and returned to access Omega^1
//'     objects as needed and cache them based on the theta that is
//'     used.
//' @param type The type of object.  Currently the following types are
//'     supported:
//' \itemize{
//' \item \code{cholOmegaInv} gives the
//'     Cholesky decomposition of the Omega Inverse matrix.
//' \item \code{omegaInv} gives the Omega Inverse matrix.
//' \item \code{d(omegaInv)} gives the d(Omega^-1) withe respect to the
//'     theta parameter specified in \code{thetaNumber}.
//' \item \code{d(D)} gives the d(diagonal(Omega^-1)) with respect to
//'     the theta parameter specified in the \code{thetaNumber}
//'     parameter
//' }
//' @param thetaNumber For types \code{d(omegaInv)} and \code{d(D)},
//'     the theta number that the derivative is taken against.  This
//'     must be positive from 1 to the number of thetas defining the
//'     Omega matrix.
//' @return Matrix based on parameters or environment with all the
//'     matrixes calculated in variables omega, omegaInv, dOmega,
//'     dOmegaInv.
//' @author Matthew L. Fidler
//' @export
// [[Rcpp::export]]
RObject rxSymInvChol(RObject invObjOrMatrix, Nullable<NumericVector> theta = R_NilValue, std::string type = "cholOmegaInv", int thetaNumber = 0){
  if (invObjOrMatrix.isObject()){
    List invObj  = as<List>(invObjOrMatrix);
    if (theta.isNull()){
      // Missing theta
      Environment base = R_BaseEnv;
      Function newEnv = as<Function>(base["new.env"]);
      Environment e = newEnv(_["parent"] = R_EmptyEnv);
      e["invobj"] = invObj;
      List ret = Rcpp::List::create(Rcpp::Named("env")=e);
      ret.attr("class") = "rxSymInvCholEnv";
      return ret;
    } else {
      NumericVector par(theta);
      int tn = thetaNumber;
      if (type == "cholOmegaInv"){
        tn = 0;
      } else if (type == "omegaInv"){
        tn = -1;
      } else if (type == "d(omegaInv)"){
        if (tn <= 0){
          stop("Theta number must be positive for d(omegaInv).");
        }
      } else if (type == "d(D)"){
        if (tn <= 0){
          stop("Theta number must be positive for d(D).");
        }
        tn = -2 - tn;
      } else if (type == "ntheta"){
        tn = -2;
      }
      try {
        Function fn = as<Function>(invObj["fn"]);
        return fn(par, tn);
      } catch (...) {
        stop("Unspported invobj type.");
      }
    }
  } else  {
    Environment rxode("package:RxODE");
    Function rxSymInvCholCreate = as<Function>(rxode["rxSymInvCholCreate"]);
    return rxSymInvChol(rxSymInvCholCreate(invObjOrMatrix), R_NilValue, "cholOmegaInv", 0);
  }
  return R_NilValue;
}

// [[Rcpp::export]]
RObject rxSymInvCholEnvCalculate(List obj, std::string what, Nullable<NumericVector> theta = R_NilValue){
  Environment e = as<Environment>(obj["env"]);
  if (theta.isNull()){
    if (e.exists(what)){
      return e[what];
    } else if (what == "theta"){
      return R_NilValue;
    } else {
      List invObj;
      if (e.exists("invobj")){
        invObj = as<List>(e["invobj"]);
      } else {
        stop("Error in rxSymInvCholEnvCalculate environment.");
      }
      if (what == "ntheta"){
        e["ntheta"] = rxSymInvChol(invObj,NumericVector::create(1),"ntheta",0);
        return(e["ntheta"]);
      }
      NumericVector theta;
      if (e.exists("theta")){
        theta = as<NumericVector>(e["theta"]);
      } else {
        stop("theta for omega calculations not setup yet.");
      }
      int ntheta = theta.size(), i=0;
      if (what == "chol.omegaInv"){
        e["chol.omegaInv"]=as<NumericMatrix>(rxSymInvChol(invObj, theta, "cholOmegaInv"));
      } else if (what == "omegaInv"){
        e["omegaInv"]= as<NumericMatrix>(rxSymInvChol(invObj, theta, "omegaInv",-1));
      } else if (what == "d.omegaInv"){
        List ret(ntheta);
        for (i =0; i < ntheta; i++){
          ret[i] = as<NumericMatrix>(rxSymInvChol(invObj, theta, "d(omegaInv)",i+1));
        }
        e["d.omegaInv"] = ret;
      } else if (what == "d.D.omegaInv"){
        List ret(ntheta);
        for (i =0; i < ntheta; i++){
          ret[i] = as<NumericVector>(rxSymInvChol(invObj, theta, "d(D)",i+1));
        }
        e["d.D.omegaInv"] = ret;
      } else if (what == "chol.omega1"){
        rxSymInvCholEnvCalculate(obj, "chol.omegaInv", R_NilValue);
        arma::mat ret = rxToCholOmega(as<arma::mat>(e["chol.omegaInv"]));
        e["chol.omega1"] = ret; 
      } else if (what == "omega"){
        rxSymInvCholEnvCalculate(obj, "chol.omega1", R_NilValue);
        arma::mat U1 = as<mat>(e["chol.omega1"]);
        arma::mat omega = U1*trans(U1);
        e["omega"] = omega;
      } else if (what == "chol.omega"){
	rxSymInvCholEnvCalculate(obj, "omega", R_NilValue);
        arma::mat omega = as<mat>(e["omega"]);
        e["chol.omega"] = chol(omega);
      } else if (what == "log.det.OMGAinv.5"){
        // Note this does NOT include the 2 pi bit
        rxSymInvCholEnvCalculate(obj,"chol.omegaInv", R_NilValue);
        arma::mat c = as<arma::mat>(e["chol.omegaInv"]);
        arma::vec diag = c.diag();
        arma::vec ldiag = log(diag);
        NumericVector ret = as<NumericVector>(wrap(sum(ldiag)));
        e["log.det.OMGAinv.5"] = ret;
      }
      return e[what];
    }
  } else {
    if (what == "theta"){
      NumericVector par(theta);
      int ntheta = as<int>(rxSymInvCholEnvCalculate(obj, "ntheta", R_NilValue));
      if (par.size() == ntheta){
        // Clear cache with the exception of
        CharacterVector sym = e.ls(TRUE);
	// Clear the cache
        for (int i = 0; i < sym.size(); i++){
          if (sym[i] != "invobj" && sym[i] != "ntheta") {
            e.remove(as<std::string>(sym[i]));
          }
        }
        e["theta"] = par;
        return(obj);
      } else {
	stop("theta has to have %d elements.", ntheta);
      }
    } else {
      stop("Can only assign 'theta' in this environment.");
    }
  }
  return R_NilValue;
}

// [[Rcpp::export]]
void RxODE_finalize_focei_omega(RObject rho){
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
  e["tr.omegaInv.dOmega.0.5"] = trInv;
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

