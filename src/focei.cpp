// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include "RxODE_types.h"

using namespace Rcpp;
using namespace R;
using namespace arma;

extern "C" SEXP RxODE_ode_solver_focei_eta (SEXP sexp_eta, SEXP sexp_rho);
extern "C" SEXP RxODE_ode_solver_focei_hessian(SEXP sexp_rho);

// [[Rcpp::export]]
NumericVector RxODE_focei_eta_lik(SEXP sexp_eta, SEXP sexp_rho){
  SEXP solve_env = RxODE_ode_solver_focei_eta(sexp_eta, sexp_rho);
  Environment e = as<Environment>(solve_env);
  mat omegaInv = as<mat>(e["omegaInv"]);
  vec eta = as<vec>(e["eta"]);
  vec aret = -(as<vec>(e["llik"])-0.5*(eta.t() * omegaInv * eta));
  e["llik2"] = wrap(aret);
  NumericVector ret = as<NumericVector>(wrap(aret));
  return ret;
}

// [[Rcpp::export]]
NumericVector RxODE_focei_eta_lp(SEXP sexp_eta, SEXP sexp_rho){
  SEXP solve_env = RxODE_ode_solver_focei_eta(sexp_eta, sexp_rho);
  Environment e = as<Environment>(solve_env);
  mat omegaInv = as<mat>(e["omegaInv"]);
  vec eta = as<vec>(e["eta"]);
  vec aret = -(as<vec>(e["lp"])- omegaInv * eta);
  e.assign("ep2",aret);
  NumericVector ret = as<NumericVector>(wrap(aret));
  return ret;
}

// [[Rcpp::export]]
XPtr<rxFn2> RxODE_focei_eta(std::string fstr){
  if (fstr == "lik")
    return(XPtr<rxFn2>(new rxFn2(&RxODE_focei_eta_lik)));
  else if (fstr == "lp")
    return(XPtr<rxFn2>(new rxFn2(&RxODE_focei_eta_lp)));
  else 
    return XPtr<rxFn2>(R_NilValue); // runtime error as NULL no XPtr
}

// [[Rcpp::export]]
NumericVector RxODE_focei_finalize_llik(SEXP rho){
  RxODE_ode_solver_focei_hessian(rho);
  Environment e = as<Environment>(rho);
  // Calculate -1/2 log(det(-H)) by chol.
  mat c = chol(-as<mat>(e["H"]));
  vec diag = c.diag();
  vec ldiag = log(diag);
  NumericVector ret = -as<NumericVector>(e["llik2"]);
  // log(det(omegaInv^1/2)) = 1/2*log(det(omegaInv))
  ret += as<NumericVector>(e["log.det.OMGAinv.5"]);
  ret += -as<NumericVector>(wrap(sum(ldiag)));
  ret.attr("fitted") = as<NumericVector>(e["f"]);
  ret.attr("posthoc") = as<NumericVector>(e["eta"]);
  e["ret"] = ret;
  return ret;
}

// [[Rcpp::export]]
NumericVector RxODE_finalize_log_det_OMGAinv_5(SEXP rho){
  // log(det(omegaInv^1/2)) = 1/2*log(det(omegaInv))
  Environment e = as<Environment>(rho);
  mat c = chol(as<mat>(e["omegaInv"]));
  vec diag = c.diag();
  vec ldiag = log(diag);
  NumericVector ret = as<NumericVector>(wrap(sum(ldiag)));
  e["log.det.OMGAinv.5"] = ret;
  return ret;
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
  e["tr.omegaInv.dOmega.0.5"] = trInv;
  e["omegaInv.dOmega.omegaInv"] = prod1;
  e["omegaInv.dOmega.omegaInv.dEta"] = prod2;
}

//' Calculate d(eta)/d(omega)
//'
//' @param eta the eta to caluclate the differential for.
//'
//' @param rho environment where omegaInv.dOmega.omegaInv and tr.omegaInv.dOmega.0.5
//' are calculated.  This is done with the rxSymInv function.
//'
//' @return Nothing.  Add omega.28 and omega.47 to the environment rho.
//' 
//' @keywords internal
//' @export
// [[Rcpp::export]]
void rxDetaDomega(SEXP rho){
  // Used in  Eq #28
  Environment e = as<Environment>(rho);
  List dOmega = as<List>(e["omegaInv.dOmega.omegaInv"]);
  NumericVector  omegaInv = as<NumericVector>(e["tr.omegaInv.dOmega.0.5"]);
  List dOmega2 = as<List>(e["omegaInv.dOmega.omegaInv.dEta"]);
  mat eta = as<mat>(e["eta.mat"]);
  mat c,c2;
  vec ret;
  int ntheta = dOmega.length();
  NumericVector dEta(ntheta);
  mat o0 = as<mat>(dOmega[0]);
  int neta = o0.n_rows;
  int i,j;
  mat dEta47 = mat(neta,ntheta);
  for (i = 0; i < ntheta; i++){
    c = 0.5*(eta.t() * as<mat>(dOmega[i]) * eta);
    ret = c.diag();
    dEta[i] = sum(ret)-omegaInv[i];
    List dOmega2T = dOmega2[i];
    for (j = 0; j < neta; j++){
      c2 = eta.t() * as<mat>(dOmega2T[j]);
      ret = c2.diag();
      dEta47(j, i) = sum(ret);
    }
  }
  e["omega.28"] = dEta;
  e["omega.47"] = dEta47;
}

// [[Rcpp::export]]
void rxDetaDtheta(SEXP rho, Function f){
  int i,h,n;
  Environment e = as<Environment>(rho);
  mat H2 = as<mat>(e["H2"]);
  mat lDnDt = as<mat>(e["l.dEta.dTheta"]);
  mat iH2 = inv(H2);
  mat DnDt = -iH2 * lDnDt;
  e["dEta.dTheta"] = DnDt;
  // Now  (dErr/dTheta)*
  mat dErrdTheta= as<mat>(e["dErr.dTheta"]);
  mat dErr = as<mat>(e["dErr"]);
  int ntheta = dErrdTheta.n_cols;
  // matrix(tmp2$dErr.dTheta[,1]) + tmp2$dErr %*% matrix(tmp2$dEta.dTheta[,1])
  mat dErrdTheta_ = mat(dErrdTheta.n_rows,0);
  for (i = 0; i < ntheta; i++){
    mat cur = dErrdTheta.col(i)+dErr * DnDt.col(i);
    dErrdTheta_ = join_rows(cur,dErrdTheta_);
  }
  e["dErr.dTheta."] = dErrdTheta_;
  // Now (dR/dTheta)*
  mat dRdTheta= as<mat>(e["dR.dTheta"]);
  mat dR = as<mat>(e["dR"]);
  // matrix(tmp2$dErr.dTheta[,1]) + tmp2$dErr %*% matrix(tmp2$dEta.dTheta[,1])
  mat dRdTheta_ = mat(dRdTheta.n_rows,0);
  mat cur;
  for (i = 0; i < ntheta; i++){
    cur = dRdTheta.col(i)+dR * DnDt.col(i);
    dRdTheta_ = join_rows(cur,dRdTheta_);
  }
  e["dR.dTheta."] = dRdTheta_;
  // Now #37
  // tmp2$dErr.dEta.dTheta[[theta]][,eta]-sum(over eta1,matrix(rowSums(tmp2$dErr2[[eta]][,eta1]*tmp2$dEta.dTheta[eta1,theta])))
  // tmp2$dErr.dEta.dTheta[[2]][,2]-matrix(rowSums(tmp2$dErr2[[1]]*tmp2$dEta.dTheta[2,2]))
  int neta = DnDt.n_rows;
  List dErrdEtadTheta_(neta);
  List dErrdEtadTheta = as<List>(e["dErr.dEta.dTheta"]);
  List dErr2 = as<List>(e["dErr2"]);
  mat mat1, mat2;
  for (i = 0; i < neta; i++){
    cur = mat(dErrdTheta.n_rows,0);
    for (h=0; h < ntheta; h++){
      mat1 = as<mat>(dErrdEtadTheta[h]);
      mat1 = -mat1.col(i);
      mat2 = as<mat>(dErr2[i]);
      for (n = 0; n < neta; n++){
	mat1 = mat1 - mat2.col(n)*DnDt(n,h);
      }
      cur = join_rows(cur,mat1);
    }
    dErrdEtadTheta_[i]=cur;
  }
  e["dErr.dEta.dTheta."] = dErrdEtadTheta_;
  // And #37 equavialent for dR.
  List dRdEtadTheta_(neta);
  List dRdEtadTheta = as<List>(e["dR.dEta.dTheta"]);
  List dR2 = as<List>(e["dR2"]);
  for (i = 0; i < neta; i++){
    cur = mat(dRdTheta.n_rows,0);
    for (h=0; h < ntheta; h++){
      mat1 = as<mat>(dRdEtadTheta[h]);
      mat1 = mat1.col(i);
      mat2 = as<mat>(dR2[i]);
      for (n = 0; n < neta; n++){
        mat1 = mat1 + mat2.col(n)*DnDt(n,h);
      }
      cur = join_rows(cur,mat1);
    }
    dRdEtadTheta_[i]=cur;
  }
  e["dR.dEta.dTheta."] = dRdEtadTheta_;
  // Now dc*/dTheta #32
  // as.matrix(tmp2$dR.dTheta.[,theta]) * tmp2$dR[,eta]/(tmp2$R*tmp2$R)+tmp2$dR.dEta.dTheta.[[eta]][,theta]/tmp2$R
  List DcDh(neta);
  mat R = as<mat>(e["R"]);
  for (n = 0; n < neta; n++){
    cur = mat(dRdTheta.n_rows,0);
    for (h = 0; h < ntheta; h++){
      mat1 = as<mat>(dRdEtadTheta_[n]);
      mat2 = -dRdTheta_.col(h)  % dR.col(n)/(R % R)+mat1.col(h)/R;
      cur = join_rows(cur,mat2);
    }
    DcDh[n]=cur;
  }
  e["dc.dTheta"] = DcDh;
  // Now dB*/dTheta #31
  mat DbDh = mat(dRdTheta.n_rows,0);
  for (h = 0; h < ntheta; h++){
    mat1 = -2*dRdTheta_.col(h)/(R % R);
    DbDh = join_rows(DbDh,mat1);
  }
  e["dB.dTheta"] = DbDh;
  // Now da*/dTheta #30
  // matrix(tmp2$dErr.dEta.dTheta.[[eta]][,theta]) + matrix(tmp2$dErr.dTheta.[,theta])*matrix(tmp2$dR[,eta])/tmp2$R + matrix(tmp2$err)*matrix(tmp2$dR.dTheta.[,theta])*matrix(tmp2$dR[,eta])/(tmp2$R*tmp2$R) - tmp2$err*matrix(tmp2$tmp2$dR.dEta.dTheta.[[eta]][,theta])/tmp2$R
  // matrix(tmp2$dErr.dEta.dTheta.[[1]][,1]) + matrix(tmp2$dErr.dTheta.[,1]) * matrix(tmp2$dR[,1])/tmp2$R + matrix(tmp2$err)*matrix(tmp2$dR.dTheta.[,1])*matrix(tmp2$dR[,1])/(tmp2$R*tmp2$R) - tmp2$err*matrix(tmp2$dR.dEta.dTheta.[[1]][,1])/tmp2$R
  List DaDh(neta);
  mat mat3;
  mat err =as<mat>(e["err"]);
  NumericVector do_nonmem_v = as<NumericVector>(e["nonmem"]);
  int do_nonmem = (int)(do_nonmem_v[0]);
  for (n = 0; n < neta; n++){
    cur = mat(dRdTheta.n_rows,0);
    for (h = 0; h < ntheta; h++){
      mat1 = as<mat>(dRdEtadTheta_[n]);
      mat2 = as<mat>(dErrdEtadTheta_[n]);
      if (do_nonmem){
	mat3 = mat2.col(h);
      } else {
	mat3 = mat2.col(h) + dErrdTheta_.col(h) % dR.col(n) /R + err % dRdTheta_.col(h) % dR.col(n)/(R % R) - err % mat1.col(h)/R;
      }
      cur = join_rows(cur,mat3);
    }
    DaDh[n]=cur;
  }
  e["da.dTheta"] = DaDh;
  // Now calculate dH/dTheta (Eq 29)
  List DhDh(ntheta);
  int k, l;
  mat al, ak, dal, dak, cl, ck, dcl, dck;
  List a = as<List>(e["a"]);
  List c = as<List>(e["c"]);
  mat B = as<mat>(e["B"]);
  int ptheta = as<int>(e["ntheta"]);
  List dOmegainv = as<List>(e["dOmegaInv"]);
  for (h = 0; h < ntheta; h++){
    mat1 = mat(neta, neta);
    for (k = 0; k < neta; k++){
      for (l = 0; l <= k; l++){
	al = as<mat>(a[l]);
	ak = as<mat>(a[k]);
	dal = as<mat>(DaDh[l]);
	dak = as<mat>(DaDh[k]);
	cl = as<mat>(c[l]);
        ck = as<mat>(c[k]);
        dcl = as<mat>(DcDh[l]);
        dck = as<mat>(DcDh[k]);
	if (do_nonmem){
	  mat1(k,l) = -0.5*sum(dal.col(h) % B % al + al % DbDh.col(h) % ak + al % B % dak.col(h) + dcl.col(h) % ck + cl % dck.col(h));
	} else {
	  mat1(k,l) = -0.5*sum(dal.col(h) % B % al + al % DbDh.col(h) % ak + al % B % dak.col(h) - dcl.col(h) % ck - cl % dck.col(h));
	}
	if (h >= ptheta){
	  // Put in dOmega^-1/dTheta term.
	  mat2 = as<mat>(dOmegainv[h-ptheta]);
	  mat1(k,l) = mat1(k,l)-mat2(k,l);
	}
	mat1(l,k) = mat1(k,l);
      }
    }
    DhDh[h]=mat1;
  }
  e["dH.dTheta"] = DhDh;
  // Now caluclate dl(eta)/dTheta (Eq 28) and add to the overall dl/dTheta
  RxODE_ode_solver_focei_hessian(e); // Calculate Hessian
  mat H = as<mat>(e["H"]);
  mat Hinv = inv(H);
  e["Hinv"] = Hinv;
  NumericVector dEta = as<NumericVector>(e["omega.28"]);
  NumericVector dLdTheta(ntheta);
  for (h = 0; h < ntheta; h++){
    mat1 = -0.5*sum(2 * err % dErrdTheta.col(h) / R - err % err % dRdTheta.col(h) / (R % R) + dRdTheta.col(h) / R);
    if (h >= ptheta){
      mat1 = mat1+ dEta[h-ptheta];
    }
    // Now add -1/2*tr(Hinv*DhDh[h])
    mat2 = as<mat>(DhDh[h]);
    mat3 = Hinv * mat2;
    dLdTheta[h] = mat1(0,0)-0.5*sum(mat3.diag());
  }
  e["l.dTheta"] = dLdTheta;
  mat omegaInv = as<mat>(e["omegaInv"]);
  mat eta = as<mat>(e["eta.mat"]);
  mat aret = -(as<vec>(e["llik"])-0.5*(eta.t() * omegaInv * eta));
  e["llik2"] = aret(0,0);
  mat cH = chol(-as<mat>(e["H"]));
  vec diag = cH.diag();
  vec ldiag = log(diag);
  NumericVector ret = aret(0,0);
  // log(det(omegaInv^1/2)) = 1/2*log(det(omegaInv))
  ret = as<NumericVector>(e["log.det.OMGAinv.5"])-ret;
  ret += -as<NumericVector>(wrap(sum(ldiag)));
  ret.attr("fitted") = as<NumericVector>(e["f"]);
  ret.attr("posthoc") = as<NumericVector>(wrap(e["eta.mat"]));
  ret.attr("grad") = as<NumericVector>(wrap(dLdTheta));
  ret.attr("dEta.dTheta") = as<NumericVector>(wrap(DnDt));
  e["ret"] = ret;
}

