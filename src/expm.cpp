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

static inline arma::mat matrixExp(arma::mat& mat, int& type){
  switch(type){
  default:
    return (arma::expmat(mat));
  }
}

//' phiv3
//'
//' Approximates w = exp(t*A)*v + t*phi(t*A)*u using Krylov
//' subspace projection, where phi(z) = (exp(z)-1)/z and w is
//' the solution of the nonhomogeneous ODE w' = Aw + u, w(0) = v.
//' 
//' @param t Solve time
//' @param A Exponential matrix A
//' @param u U vector; Equivalent to F in indLin
//' @param v Initial condition vector
//' @param anorm absoulte normalization value for round-off error.
//' @param tol Tolerance
//' @param maxsteps Maximum number of steps
//' 
//' @author Roger B. Sidje, Wenping Wang, Matthew Fidler
//' 
arma::mat phiv3(double t, const arma::mat& A, const arma::vec& u,
		const arma::vec& v,
		double anorm, double tol, int m,
		int type){
  int n = A.n_rows, mxrej=100, mb=m, k1=3;
  double btol=1.0e-7, gamma=0.9, delta=1.2;
  double t_out=fabs(t), t_new=0, t_now=0, s_error=0;
  double eps=std::numeric_limits<double>::epsilon();
  double rndoff=anorm*eps; if(tol < eps) tol = sqrt(eps);
  int mh=m+3;
  arma::mat H0(mh, mh, arma::fill::zeros);
  arma::mat w = v;
  double xm=1.0/double(m);
  arma::mat F, H, Ht;
  arma::mat V(n, m+1);
  arma::mat p;
  arma::mat av;
  double beta;
  double err_loc, avnorm=0;
  double s;
  double h=0;
  int i, j, mx;
  
  for (int istep = 0; t_now < t_out; ++istep) {
    H = H0;
    V.col(0) = A*w + u;
    beta = norm(V.col(0));
    if (beta==0) break;
    V.col(0) /= beta;
    if (istep == 0) {
       double fact = (pow((m+1)/M_E, m+1))*sqrt(M_2PI*(m+1));
       t_new = (1.0/anorm)*pow((fact*tol)/(4*beta*anorm),xm);
       s = pow(10.0,floor(log10(t_new))-1.0);
       t_new = ceil(t_new/s)*s;
    }
    double t_step = fmin( t_out-t_now,t_new );
    for (j = 0; j < m; ++j) {
       p = A*V.col(j);
       for (i = 0; i <= j; ++i) {
	 H(i,j) = arma::dot(V.col(i),p.col(0));
	 p -= H(i,j)*V.col(i);
       }
       s = norm(p);
       if (s < btol) {
          k1 = 0;
          mb = j;
          t_step = t_out-t_now;
          break;
       }
       H(j+1,j) = s;
       V.col(j+1) = p/s;
    }
    H(0,mb) = 1;
    if (k1 != 0) {
       H(m,m+1) = 1; H(m+1,m+2) = 1;
       h = H(m,m-1); H(m,m-1) = 0;
       av=A*V.col(m);
       avnorm = norm(av);
    }
    for (int irej = 0; irej < mxrej; ++irej) {
      //cout << t_step << endl;
      mx = mb + std::max(1,k1);
      Ht = t_step*H.submat(0, 0, mx-1, mx-1);
      F = matrixExp(Ht, type);
      if (k1 == 0) {
	err_loc = btol;
	break;
      } else {
	F(m  ,m) = h*F(m-1,m+1);
	F(m+1,m) = h*F(m-1,m+2);
	double p1 = fabs( beta*F(m,  m));
	double p2 = fabs( beta*F(m+1,m) * avnorm);
	if (p1 > 10.0*p2) {
	  err_loc = p2;
	  xm = 1.0/double(m+1);
	} else if (p1 > p2) {
	  err_loc = (p1*p2)/(p1-p2);
	  xm = 1.0/double(m+1);
	} else {
	  err_loc = p1;
	  xm = 1.0/double(m);
	}
      }
      if (err_loc <= delta*t_step*tol)
	break;
      else {
	t_step = gamma * t_step * pow(t_step*tol/err_loc, xm);
	double s = pow(10.0, floor(log10(t_step))-1.0);
	t_step = ceil(t_step/s) * s;
	if (irej == mxrej){
	  // FIXME must reject the point and put NA in...
	  stop("The requested tolerance is too high.\n");
	}
      }
    }
    mx = mb + std::max( 0, k1-2 );
    if (mx > 0){
      w += V.submat(0, 0, n-1, mx-1)*(beta*F.submat(0, mb, mx-1, mb));
    }
    t_now += t_step;
    t_new = gamma * t_step * pow(t_step*tol/err_loc, xm);
    s = pow(10.0, floor(log10(t_new))-1.0);
    t_new = ceil(t_new/s) * s;

    err_loc = fmax(err_loc, rndoff);
    s_error += err_loc;
  }
  return w;
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
//' @param m Expokits m parameter
//' @param tol Expokit's tol parameter
//' @param cache
//'    0 = no Cache
//'    When doIndLin == 0, cache > 0 = nInf-1
//' @param ME the RxODE matrix exponential function
//' @param IndF The RxODE Inductive Linearization function F
//' 
//' @return Returns a status for solving
//' 
//'   1 = Successful solve
//' 
//'   -1 = Maximum number of iterations reached when doing
//'        inductive linearization
extern "C" int indLin(int cSub, int neq, double tp, double *yp_, double tf,
		      double *InfusionRate_, int *on_, double *rtol, double *atol,
		      int maxsteps, int doIndLin, int locf,
		      int phiM, double phiTol, double phiAnorm,
		      double *rwork, int *cache, int type,
		      t_ME ME, t_IndF  IndF){
  std::ostream nullstream(0);
  arma::set_cerr_stream(nullstream);
  double *ptr = &rwork[0];
  arma::mat m0(ptr, neq, neq, false, false);
  ptr += neq*neq;
  double tcov = tf;
  if (locf) tcov = tp;
  // LOCF=tp; If NOCB tp should be tf
  if (*cache == 0) ME(cSub, tcov, m0.memptr()); // Calculate the initial A matrix based on current time/parameters
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
    if (*cache == 0){
      m0extra.zeros();
      // arma::mat mout;
      for (i = 0; i < (unsigned int)neq; i++){
	if (InfusionRate[i] != 0.0){
	  nInf++;
	  m0extra[neq*(nInf-1)+i]=1;
	  ypExtra[i] = InfusionRate[i];
	}
      }
    } else {
      nInf = *cache-1;
    }
    if (nInf == 0){
      arma::mat expAT(ptr, neq, neq, false, false);
      ptr += neq*neq;
      expAT = arma::expmat(m0*(tf-tp));
      arma::vec meSol(ptr, neq, false, false);
      ptr += neq;
      meSol = expAT*yp;
      std::copy(meSol.begin(), meSol.end(), yp_);
      *cache = 1;
      return 1;
      // mout = m0;
      // ypout=yp;
    } else {
      arma::mat mout(ptr, neq+nInf, neq+nInf, false, false);
      ptr += (neq+nInf)*(neq+nInf);
      arma::vec ypout(ptr, neq+nInf, false, false);
      ptr += (neq+nInf);
      if (*cache == 0){
	mout.zeros();
	for (int j = neq; j--;){
	  std::copy(m0.colptr(j), m0.colptr(j)+neq, mout.colptr(j));
	}
	for (int j = nInf; j--;){
	  std::copy(m0extra.colptr(j),m0extra.colptr(j)+neq, mout.colptr(neq+j));
	}
	std::copy(yp.begin(),yp.end(),ypout.begin());
	std::copy(ypExtra.begin(),ypExtra.end(), ypout.begin()+neq);
      }
      arma::mat expAT(ptr, neq+nInf, neq+nInf, false, false);
      ptr += (neq+nInf)*(neq+nInf);
      // Unfortunately the tf-tp may change so we cann't cache this.
      expAT = mout*(tf-tp);
      expAT = matrixExp(expAT, type);
      arma::vec meSol(ptr, neq+nInf, false, false);
      meSol = expAT*ypout;
      std::copy(meSol.begin(), meSol.begin()+neq, yp_);
      return 1;
    }
  } else {
    // In this case the inital matrix should not be expanded. The
    // infusions are put into the F function
    const arma::vec InfusionRate(InfusionRate_, neq, false, false);

    const arma::vec yp(yp_, neq, false, false);
    
    arma::vec u(ptr, neq, false, false);
    ptr += neq;
    arma::vec w(ptr, neq, false, false);
    ptr += neq;
    arma::vec wLast(ptr, neq, false, false);
    ptr += neq;
    double *fptr = u.memptr();
    // For LOCF tp for NOCB tf
    IndF(cSub, tcov, tf, fptr, yp_, InfusionRate_);
    // IndF(cSub, tcov, tf, fptr, yp_, InfusionRate_, u.memptr())
    wLast = phiv3((tf-tp), m0, u, yp, phiAnorm, phiTol, phiM, type);
    IndF(cSub, tcov, tf, fptr, wLast.memptr(), InfusionRate_);
    w=phiv3((tf-tp), m0, u, yp, phiAnorm, phiTol, phiM, type);
    bool converge = false;
    for (int i = 0; i < maxsteps; ++i){
      converge=true;
      for (int j=neq;j--;){
	if (fabs(w[j]-wLast[j]) >= rtol[j]*fabs(w[j])+atol[j]){
	  converge = false;
	  break;
	}
      }
      if (converge){
	break;
      }
      wLast = w;
      IndF(cSub, tcov, tf, fptr, wLast.memptr(), InfusionRate_);
      w=phiv3((tf-tp), m0, u, yp, phiAnorm, phiTol, phiM, type);
    }
    if (!converge){
      std::copy(w.begin(), w.end(), &yp_[0]);
      // std::fill_n(&yp_[0], neq, NA_REAL);
      return 1;
    } else {
      std::copy(w.begin(), w.end(), &yp_[0]);
      return 1;
    }
  }
  return 1;
}
