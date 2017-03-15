// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include "RxODE_types.h"

using namespace Rcpp;
using namespace R;
using namespace arma;

extern "C" SEXP RxODE_ode_solver_focei_eta (SEXP sexp_eta, SEXP sexp_rho);
extern "C" int nEq();
extern "C" unsigned int nLhs();
extern "C" unsigned int nAllTimes();
extern "C" int rxEvid(int i);
extern "C" double rxLhs(int i);
extern "C" void rxCalcLhs(int i);
extern "C" void RxODE_ode_solve_env(SEXP sexp_rho);
extern "C" unsigned int nObs();
extern "C" void RxODE_ode_solve_env(SEXP sexp_rho);
extern "C" void RxODE_ode_free();
extern "C" double RxODE_safe_zero(double x);
extern "C" double RxODE_safe_log(double x);

// [[Rcpp::export]]
void rxInner(SEXP etanews, SEXP rho){
  Environment e = as<Environment>(rho);
  if (!(e.exists("neta") && e.exists("ntheta") && e.exists("dOmega") &&
	e.exists("DV") && e.exists("nonmem") && e.exists("eta") &&
	e.exists("eta.mat") && e.exists("eta.trans") &&
	e.exists("params")
	)){
    stop("Environment not setup correctly for rxInner.");
  }
  NumericVector par_ptr = as<NumericVector>(e["params"]);
  IntegerVector eta_i = as<IntegerVector>(e["eta.trans"]);
  NumericVector etanew = as<NumericVector>(etanews);
  NumericVector eta = as<NumericVector>(e["eta"]);
  mat etam = as<NumericVector>(e["eta.mat"]);
  unsigned int recalc = 0;
  unsigned int i = 0, j = 0, k = 0;
  if (!e.exists("llik")){
    recalc = 1;
  } else if (eta.size() != etanew.size()){
    stop("Inconsistent eta size for rxInner.");
  } else {
    for (i = 0; i < (unsigned int)(eta.size()); i++){
      if (eta[i] != etanew[i]){
	recalc = 1;
	break;
      }
    }
  }
  if (recalc){
    for (i = 0; i < (unsigned int)(etanew.size()); i++){
      /* Rprintf("\tpar[%d] from %f to %f\n",eta_i[i],par_ptr[eta_i[i]], eta[i]); */
      par_ptr[eta_i[i]] = etanew[i];
      eta[i] = etanew[i];
      etam(i,0) = eta[i];
    }
    e["params"] = par_ptr;
    e["eta"] = eta;
    e["eta.mat"] = etam;

    RxODE_ode_solve_env(rho);
    
    unsigned int neta = as<unsigned int>(e["neta"]);
    unsigned int ntheta = as<unsigned int>(e["ntheta"]);
    List dOmega = as<List>(e["dOmega"]);
    NumericVector DV = as<NumericVector>(e["DV"]);
    int do_nonmem = as<int>(e["nonmem"]);
  
    unsigned int nomega = (unsigned int)(dOmega.size());
  
    mat fpm = mat(nObs(), neta);
    mat fpt = mat(nObs(),ntheta+nomega);
    NumericVector fp2(neta);
  
    mat rp = mat(nObs(),neta);
    NumericVector rp2(neta);
    mat rpt = mat(nObs(),ntheta+nomega);

    NumericVector f(nObs());
    mat err = mat(nObs(),1);
    mat r = mat(nObs(),1);

    mat B = mat(nObs(),1);
    List c(neta);
    List a(neta);
  
    List fpte(ntheta+nomega);
    List rpte(ntheta+nomega);

    NumericVector llik(1);
    mat lp = mat(neta,1);

    mat lDnDt = mat(neta,ntheta+nomega);
    mat lDn = mat(neta,neta);

    for (i = 0; i < neta; i++){
      a[i] = mat(nObs(),1);
      c[i] = mat(nObs(),1);
      lp[i] = 0;
    }

    llik[0]=0;

    /* // Now create the pred vector and d(pred)/d(eta) matrix. */
    /* // Assuming rxLhs(0) = pred and rxLhs(1:n) = d(pred)/d(eta#) */
    // Solve
    mat cur, cuR;
    for (i = 0; i < nAllTimes(); i++){
      if (!rxEvid(i)){
	rxCalcLhs(i);
	f[k] = rxLhs(0); // Pred
	err(k, 0) = DV[k] - f[k];
	// d(pred)/d(eta#)
	for (j = 1; j < neta+1; j++){
	  fpm(k, j-1) = rxLhs(j);
	  cur = as<mat>(a[j-1]);
	  if (do_nonmem){
	    cur(k,0) =  rxLhs(j);
	  } else {
	    cur(k,0) = rxLhs(j)-err(k, 0)/RxODE_safe_zero(rxLhs(neta+1))*rxLhs(j+neta+1);
	  }
	  a[j-1]=cur;
	}
	if (rxLhs(j) < 0){
	  for (j = 0; j < nLhs(); j++){
	    Rprintf("rxLhs(%d) = %f\n", j, rxLhs(j));
	  }
	  Rprintf("\n");
	  // temp = getAttrib(sexp_theta, R_NamesSymbol);
	  // for (j = 0; j < length(sexp_theta); j++){
	  //   Rprintf("params[%d] = %f\n", j, par_ptr[j]);
	  // }
	  RxODE_ode_free();
	  stop("A covariance term is zero or negative and should remain positive");
	}
	r(k, 0)=rxLhs(j); // R always has to be positive.
	/* logR[k]=log(rxLhs(j)); */
	/* Rinv[k]=1/rxLhs(j); */
	B(k, 0)=2/rxLhs(j);
	for (j=neta+2; j < nLhs(); j++){
	  /* Rprintf("j: %d; Adj: %d; k: %d\n",j, j-neta-2,k); */
	  rp(k,j-neta-2) = rxLhs(j);
	  cur = as<mat>(c[j-neta-2]);
	  cur(k,0) = rxLhs(j)/RxODE_safe_zero(r(k, 0));
	  c[j-neta-2] = cur;
	}
	for (j = 0; j < neta; j++){
	  // .5*apply(eps*fp*B + .5*eps^2*B*c - c, 2, sum) - OMGAinv %*% ETA
	  cur = as<mat>(c[j]);
	  lp[j] += 0.5 * err(k, 0)* fpm(k, j) * B(k, 0)  +
	    0.25 * err(k, 0) * err(k, 0) * B(k, 0) * cur(k,0) -
	    0.5 * cur(k,0);
	}
	llik[0] += -0.5*(err(k, 0)*err(k, 0)/RxODE_safe_zero(r(k, 0))+RxODE_safe_log(r(k, 0)));
	k++;
      }
    }
    // Free
    RxODE_ode_free();

    mat omegaInv = as<mat>(e["omegaInv"]);

    NumericVector llik2(1);
    mat llikm = mat(1,1);
    llikm(0, 0)=  llik[0];
    llikm = -(llikm - 0.5*(etam.t() * omegaInv * etam));
    llik2[0] = llikm(0, 0);

    mat ep2 = -(lp- omegaInv * etam);
  
    // Assign in env
    e["err"] = err;
    e["f"] = f;
    e["dErr"] = fpm;
    e["dR"] = rp;
    e["c"] = c;
    e["R"] = r;
    e["B"] = B;
    e["a"] = a;
    e["llik"] = llik;
    e["lp"] = lp;
    e["llik2"] = wrap(llik2);
    e["ep2"] = wrap(ep2);
  }
}

// [[Rcpp::export]]
void rxHessian(SEXP rho){
  Environment e = as<Environment>(rho);
  int do_nonmem = as<int>(e["nonmem"]);
  int neta = as<int>(e["neta"]);
  mat omegaInv = as<mat>(e["omegaInv"]);
  mat B = as<mat>(e["B"]);
  List c = as<List>(e["c"]);
  List a = as<List>(e["a"]);
  NumericVector f = as<NumericVector>(e["f"]);
  mat H(neta, neta);
  int k, l;
  mat al, ak, cl, ck, tmp(1,1);
  for (k = 0; k < neta; k++){
    for (l = 0; l <= k; l++){
      al = as<mat>(a[l]);
      ak = as<mat>(a[k]);
      cl = as<mat>(c[l]);
      ck = as<mat>(c[k]);
      if (do_nonmem){
	tmp  = -0.5*sum(al % B % ak + cl % ck)-omegaInv(k,l);
      } else {
	tmp  = -0.5*sum(al % B % ak - cl % ck)-omegaInv(k,l);
      }
      H(k,l) = tmp(0,0);
      // Fill out the mirror compenent.
      if (l != k){
	H(l,k)=H(k,l);
      }
    }
  }
  e["H"] = H;
}


// [[Rcpp::export]]
NumericVector RxODE_focei_eta_lik(SEXP sexp_eta, SEXP sexp_rho){
  rxInner(sexp_eta, sexp_rho);
  Environment e = as<Environment>(sexp_rho);
  NumericVector ret = as<NumericVector>(wrap(e["llik2"]));
  return ret;
}

// [[Rcpp::export]]
NumericVector RxODE_focei_eta_lp(SEXP sexp_eta, SEXP sexp_rho){
  rxInner(sexp_eta, sexp_rho);
  Environment e = as<Environment>(sexp_rho);
  NumericVector ret = as<NumericVector>(wrap(e["ep2"]));
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
  rxHessian(rho);
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
      dEta47(j, i) = c2(0, 0);
    }
  }
  e["omega.28"] = dEta;
  e["omega.47"] = dEta47;
}

// [[Rcpp::export]]
void rxOuter_ (SEXP rho){
  //Outer problem gradient for lbfgs
  Environment e = as<Environment>(rho);
  unsigned int i, j, k=0, h, n, i0 = 0,e1,e2;
  
  mat omegaInv = as<mat>(e["omegaInv"]);
  
  unsigned int neta = as<unsigned int>(e["neta"]);
  unsigned int ntheta = as<unsigned int>(e["ntheta"]);
  List dOmega = as<List>(e["dOmega"]);
  NumericVector DV = as<NumericVector>(e["DV"]);
  int do_nonmem = as<int>(e["nonmem"]);
  
  unsigned int nomega = (unsigned int)(dOmega.size());

  RxODE_ode_solve_env(rho);
  
  mat fpm = mat(nObs(), neta);
  mat fpt = mat(nObs(),ntheta+nomega);
  List fp2(neta);
  
  mat rp = mat(nObs(),neta);
  mat rpt = mat(nObs(),ntheta+nomega);
  List rp2(neta);

  NumericVector f(nObs());
  mat err = mat(nObs(),1);
  mat r = mat(nObs(),1);

  mat B = mat(nObs(),1);
  List c(neta);
  List a(neta);
  
  List fpte(ntheta+nomega);
  List rpte(ntheta+nomega);

  NumericVector llik(1);
  mat lp = mat(neta,1);

  mat lDnDt = mat(neta,ntheta+nomega);
  mat lDn = mat(neta,neta);

  for (i = 0; i < neta; i++){
    a[i] = mat(nObs(),1);
    c[i] = mat(nObs(),1);
    fp2[i] = mat(nObs(),neta);
    rp2[i] = mat(nObs(),neta);
    lp(i,0) = 0;
    for (j = 0; j < neta; j++){
      lDn(i, j) = -omegaInv(i, j);
    }
  }
  for (j = 0; j < ntheta+nomega; j++){
    fpte[j] = mat(nObs(),neta);
    rpte[j] = mat(nObs(),neta);
    for (i = 0; i < neta; i++){
      lDnDt(i,j) = 0;
    }
  }

  llik[0]=0;
  
  // Now create the pred vector and d(pred)/d(eta) matrix.
  // Assuming rxLhs(0) = pred and rxLhs(1:n) = d(pred)/d(eta#)
  mat cur, cuR;
  for (i = 0; i < nAllTimes(); i++){
    if (!rxEvid(i)){
      rxCalcLhs(i);
      f[k] = rxLhs(0); // Pred
      err(k, 0) = DV[k] - f[k];
      // d(pred)/d(eta#)
      // Rprintf("d(pred)/d(eta#)\n");
      for (j = 1; j < neta+1; j++){
        fpm(k,j-1) = rxLhs(j);
	cur = as<mat>(a[j-1]);
	cur(k,0) =  rxLhs(j);
        a[j-1]=cur;
      }
      /* // d(pred)/d(theta#) */
      // Rprintf("d(pred)/d(theta#)\n");
      i0 = 1+neta;
      for (j = i0; j < i0+ntheta; j++){
	fpt(k,j-i0) = rxLhs(j);
      }
      /* // d^2(pred)/d^2(eta#) */
      // Rprintf("d^2(pred)/d^2(eta#)\n");
      i0 += ntheta;
      e1=0; e2=0;
      for (j = i0; j < i0+(neta)*(neta+1)/2; j++){
        /* fp2[(nAllTimes()-ixds)*(j-i0)+k] = rxLhs(j); */
	cur = as<mat>(fp2[e1]);
	cur(k,e2) = rxLhs(j);
	fp2[e1] = cur;
        if (e1 == e2){
          e1=0;
          e2++;
        } else {
	  cur = as<mat>(fp2[e2]);
	  cur(k,e1) = rxLhs(j);
	  fp2[e2]=cur;
          e1++;
        }
      }
      /* // d^2(pred)/(d(eta#)d(theta#)) */
      // Rprintf("d^2(pred)/d2(eta#)d(theta#)\n");
      i0 += neta*(neta+1)/2;
      h = 0;
      for (j = i0; j < i0 + neta*ntheta; ){
        for (n = 0; n < neta; n++){
	  cur =as<mat>(fpte[h]);
	  cur(k,n) =rxLhs(j);
	  fpte[h] =cur;
          j++;
        }
        h++;
      }
      i0 += neta*ntheta;
      j=i0;
      // Now
      if (rxLhs(j) <= 0){
	Rprintf("R = rxLhs(%d) = %f\n\n", j, rxLhs(j));
	for (j = 0; j < nLhs(); j++){
          Rprintf("rxLhs(%d) = %f\n", j, rxLhs(j));
        }
        Rprintf("\n");
        RxODE_ode_free();
        stop("A covariance term is zero or negative and should remain positive.");
      }
      r(k, 0)=rxLhs(j); // R always has to be positive.
      /* logR[k]=log(rxLhs(j)); */
      /* Rinv[k]=1/rxLhs(j); */
      B(k, 0)=2/rxLhs(j);
      /* // d(R)/d(eta#) */
      // Rprintf("d(R)/d(eta#)\n");
      i0++;
      for (j=i0; j < i0+neta; j++){
        /* Rprintf("j: %d; Adj: %d; k: %d\n",j, j-neta-2,k); */
        rp(k, j-i0) = rxLhs(j);
	cur = as<mat>(c[j-i0]);
	cur(k,0) = rxLhs(j)/RxODE_safe_zero(r(k, 0));
	c[j-i0] = cur;
        if (!do_nonmem){
          // tmp1[["_sens_rx_pred__ETA_1_"]],ncol=1) - err/R*matrix(tmp1[["_sens_rx_r__ETA_1_"]]
	  cur =as<mat>(a[j-i0]);
	  cur(k, 0)+= -err(k, 0)/RxODE_safe_zero(r(k, 0))*rxLhs(j);
	  a[j-i0] = cur;
        }
      }
      i0 += neta;
      /* // d(R)/d(theta#) */
      // Rprintf("d(R)/d(theta#)\n");
      for (j=i0; j < i0+ntheta; j++){
        /* Rprintf("j: %d; Adj: %d; k: %d\n",j, j-neta-2,k); */
        rpt(k, j-i0) = rxLhs(j);
      }
      i0 += ntheta;
      /* // d^2(R)/d^2(eta) */
      // Rprintf("d(R)/d^2(eta#)\n");
      e1=0; e2=0;
      for (j=i0; j < i0+neta*(neta+1)/2; j++){
        /* Rprintf("j: %d; Adj: %d; k: %d\n",j, j-neta-2,k); */
	cur = as<mat>(rp2[e1]);
	cur(k,e2) = rxLhs(j);
	rp2[e1] = cur;
        if (e1 == e2){
          e1=0;
          e2++;
        } else {
	  cur = as<mat>(rp2[e2]);
	  cur(k,e1) = rxLhs(j);
	  rp2[e2] = cur;
          e1++;
        }
      }
      // d^2(R)/(d(eta#)d(theta#))
      // Rprintf("d^2(R)/d(eta#)d(theta#)\n");
      i0 += neta*(neta+1)/2;
      h = 0;
      for (j = i0; j < i0+ntheta*neta; ){
        for (n = 0; n < neta; n++){
	  cur = as<mat>(rpte[h]);
	  cur(k,n) = rxLhs(j);
	  rpte[h] = cur;
          j++;
        }
        h++;
      }
      // Rprintf("lp\n");
      i0 += ntheta*neta;
      for (j = 0; j < neta; j++){
        //.5*apply(eps*fp*B + .5*eps^2*B*c - c, 2, sum) - OMGAinv %*% ETA
	// eq 12
	cur = as<mat>(c[j]);
        lp(j, 0) += 0.5 * err(k, 0)* fpm(k, j)*B(k, 0)  +
          0.25 * err(k, 0) * err(k, 0) * B(k, 0) * cur(k,0) -
          0.5 * cur(k,0);
      }
      // Rprintf("47\n");
      for (h=0; h < ntheta; h++){
        for (n = 0; n < neta; n++){
          // Eq #47 Almquist 2015
	  cur = as<mat>(fpte[h]);
	  lDnDt(n, h) += -(fpt(k, h)*fpm(k, n)/r(k, 0)-
			   err(k, 0)*fpm(k, n)*rpt(k, h)/(r(k, 0)*r(k, 0))+
			   err(k, 0)*cur(k, n)/r(k, 0));
	  cur = as<mat>(rpte[h]);
	  lDnDt(n, h) += -(-0.5*err(k, 0)*err(k, 0)*cur(k, n)/(r(k, 0)*r(k, 0))+
			   err(k, 0)*err(k, 0)*rp(k, n)*rpt(k, h)/(r(k, 0)*r(k, 0)*r(k, 0))-
			   err(k, 0)*rp(k, n)*fpt(k, h)/(r(k, 0)*r(k, 0))+ // trace is not needed since R is a scalar, not a vector
			   0.5*rp(k, n)*rpt(k, h)/(r(k, 0)*r(k, 0))+
			   0.5*cur(k, n)/r(k, 0));
        }
      }
      // Rprintf("13\n");
      // Eq #13 Almquist 2015
      for (e1 = 0; e1 < neta; e1++){
        for (e2 = 0; e2 <= e1; e2++){
	  // fpm = d(err)/d(eta)
          // fpt = d(err)/d(theta)
          // fp2 = d^2(err)/d(eta)^2
          // fpte = d^2(err)/d(eta)d(theta)
	  cur = as<mat>(fp2[e1]);
	  cuR = as<mat>(rp2[e1]);
	  lDn(e1, e2) += -0.5*(2 * fpm(k, e1) * fpm(k, e2) / r(k, 0)-
			       2 * err(k, 0) *  rp(k, e2) * fpm(k,e1)/(r(k, 0) * r(k, 0))+
			       2 * err(k, 0) * cur(k,e2) / r(k, 0) -
			       err(k, 0) * err(k, 0) * cuR(k, e2) / (r(k, 0) * r(k, 0)) +
			       2 * err(k, 0) * err(k, 0) * rp(k, e1) * rp(k, e2) / (r(k, 0) * r(k, 0) * r(k, 0)) -
                               2 * err(k, 0) * rp(k, e1) * fpm(k, e2) / (r(k, 0) * r(k, 0)) -
                               rp(k, e1) * rp(k, e2) / (r(k, 0) * r(k, 0)) + cuR(k, e2) / r(k, 0));
          
        }
      }
      llik[0] += -0.5*(err(k, 0)*err(k, 0)/RxODE_safe_zero(r(k, 0))+RxODE_safe_log(r(k, 0)));
      k++;
    }
  }
  /* Finalize Eq #47 in Almquist 2015*/
  mat omega47 = as<mat>(e["omega.47"]);
  for (h=ntheta; h < ntheta+nomega; h++){
    for (n = 0; n < neta; n++){
      lDnDt(n, h) += -omega47(n,h-ntheta);
      // Finalize Eta2 and R2.
      for (i = 0; i < nObs(); i++){
	cur = as<mat>(fpte[h]);
	cur(i,n) = 0;
	fpte[h] = cur;
	cur = as<mat>(rpte[h]);
        cur(i,n) = 0;
        rpte[h] = cur;
      }
    }
    // Finalize dErr.dTheta to contain 0 for omega terms.
    for (i = 0; i < nObs(); i++){
      fpt(i, h) = 0;
      rpt(i, h) = 0;
    }
  }
  for (e1 = 0; e1 < neta; e1++){
    for (e2 = 0; e2 <= e1; e2++){
      lDn(e2, e1) = lDn(e1, e2);
    }
  }
  /* llik = -.5*sum(eps^2/(f^2*sig2) + log(f^2*sig2)) - .5*t(ETA) %*% OMGAinv %*% ETA */
  e["f"] = f;
  e["err"] = err;
  
  e["dErr"] = fpm;
  e["dErr2"] = fp2;
  e["dErr.dTheta"] = fpt;
  e["dErr.dEta.dTheta"] = fpte;

  e["R"] = r;
  
  e["dR"] = rp;
  e["dR.dTheta"] = rpt;
  e["dR2"] = rp2;
  e["dR.dEta.dTheta"] = rpte;

  e["a"] = a;
  e["B"] = B;
  e["c"] = c;

  e["llik"] = llik;
  e["lp"] = lp;

  e["l.dEta.dTheta"] = lDnDt;
  e["H2"] = lDn;
  RxODE_ode_free();
}

// [[Rcpp::export]]
void rxDetaDtheta(SEXP rho){
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
    dErrdTheta_ = join_rows(dErrdTheta_,cur);
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
    dRdTheta_ = join_rows(dRdTheta_,cur);
  }
  e["dR.dTheta."] = dRdTheta_;
  // Now #37
  // tmp2$dErr.dEta.dTheta[[theta]][,eta]-sum(over eta1,matrix(rowSums(tmp2$dErr2[[eta]][,eta1]*tmp2$dEta.dTheta[eta1,theta])))
  // tmp2$dErr.dEta.dTheta[[2]][,2]-matrix(rowSums(tmp2$dErr2[[1]]*tmp2$dEta.dTheta[2,2]))
  int neta = DnDt.n_rows;
  List dErrdEtadTheta_(neta);
  List dErrdEtadTheta = as<List>(e["dErr.dEta.dTheta"]); // dErr.dEta.dTheta
  List dErr2 = as<List>(e["dErr2"]);
  mat mat0, mat1, mat2;
  for (i = 0 ; i < neta; i++){
    cur = mat(dErrdTheta.n_rows,0);
    for(h = 0; h < ntheta; h++){
      cur = join_rows(cur, -(as<mat>(dErrdEtadTheta[h])).col(i) -
		      (as<mat>(dErr2[i])) * DnDt.col(h));
    }
    dErrdEtadTheta_[i]=cur;
  }
  e["dErr.dEta.dTheta."] = dErrdEtadTheta_;
  // And #37 equavialent for dR.
  List dRdEtadTheta_(neta);
  List dRdEtadTheta = as<List>(e["dR.dEta.dTheta"]);
  List dR2 = as<List>(e["dR2"]);
  for (i = 0 ; i < neta; i++){
    cur = mat(dRdTheta.n_rows,0);
    for(h = 0; h < ntheta; h++){
      cur = join_rows(cur, (as<mat>(dRdEtadTheta[h])).col(i) +
                      (as<mat>(dR2[i])) * DnDt.col(h));
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
      cur = join_rows(cur, mat2);
    }
    DcDh[n]=cur;
  }
  e["dc.dTheta"] = DcDh;
  // Now dB*/dTheta #31
  mat DbDh = mat(dRdTheta.n_rows,0);
  for (h = 0; h < ntheta; h++){
    mat1 = -2*dRdTheta_.col(h)/(R % R);
    DbDh = join_rows(DbDh, mat1);
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
	mat3 = mat2.col(h) - dErrdTheta_.col(h) % dR.col(n) /R + err % dRdTheta_.col(h) % dR.col(n)/(R % R) - err % mat1.col(h)/R;
      }
      cur = join_rows(cur, mat3);
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
  rxHessian(e); // Calculate Hessian
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

// [[Rcpp::export]]
NumericVector rxOuter(SEXP rho){
  rxDetaDomega(rho); // setup omega.28 and omega.47
  rxOuter_(rho);
  rxDetaDtheta(rho);
  Environment e = as<Environment>(rho);
  NumericVector ret = as<NumericVector>(e["ret"]);
  return ret;
}

