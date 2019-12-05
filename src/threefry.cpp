#define min2( a , b )  ( (a) < (b) ? (a) : (b) )
#include <RcppArmadillo.h>
#include <threefry.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <R.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
SEXP rxRmvn_(NumericMatrix A_, arma::rowvec mu, arma::mat sigma,
	     int ncores=1, bool isChol=false){
  int n = A_.nrow();
  int d = mu.n_elem;
  arma::mat ch;
  if (isChol){
    ch=arma::trimatu(sigma);
  } else {
    ch=arma::trimatu(arma::chol(sigma));
  }
  if (n < 1) stop(_("n should be a positive integer"));
  if (ncores < 1) stop(_("'ncores' has to be greater than one"));
  if (d != (int)sigma.n_cols) stop("length(mu) != ncol(sigma)");
  if (d != (int)sigma.n_rows) stop("length(mu) != ncol(sigma)");
  if (d != (int)A_.ncol()) stop("length(mu) != ncol(A)");

  double seedD = runif(1, 1.0, std::numeric_limits<uint32_t>::max())[0];
  uint32_t seed = static_cast<uint32_t>(seedD);
  seed = min2(seed, std::numeric_limits<uint32_t>::max() - ncores - 1);
#ifdef _OPENMP
#pragma omp parallel num_threads(ncores) if(ncores > 1)
  {
#endif
    arma::mat A(A_.begin(), A_.nrow(), A_.ncol(), false, true);
    sitmo::threefry eng;
    eng.seed(seed);
       
    std::normal_distribution<> snorm(0.0, 1.0);

    double acc;
    arma::rowvec work(d);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int i = 0; i < n*d; ++i){
      A[i] = snorm(eng);
    }
    if (d == 1){
      double sd = ch(0, 0);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for (int i = 0; i < n; i++){
	A[i] = A[i]*sd+mu(0);
      }
    } else {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for(int ir = 0; ir < n; ++ir){
	for(int ic = d; ic--;){
	  acc = 0.0;
	  for (int ii = 0; ii <= ic; ++ii){
	    acc += A.at(ir,ii) * ch.at(ii,ic);
	  }
	  work.at(ic) = acc; 
	}
	work += mu;
	A(arma::span(ir), arma::span::all) = work;
      }
#ifdef _OPENMP
    }
#endif
  }
  return R_NilValue;
}


// Adapted from https://github.com/cran/TruncatedNormal/blob/7364c5bc3f7c84d00eb4a767807b103b4232b648/R/ntail.R
double ntail(double l, double u, sitmo::threefry& eng){
  // samples a column vector of length=length(l)=length(u)
  // from the standard multivariate normal distribution,
  // truncated over the region [l,u], where l>0 and
  // l and u are column vectors;
  // uses acceptance-rejection from Rayleigh distr;
  // method due to Marsaglia (1964);
  std::uniform_real_distribution<> unif(0.0, 1.0);
  double c=l*l/2.0;
  double f = expm1(c-u*u/2.0);
  double x =0.0;
  bool accept = true;
  while (accept){
    double tmp = unif(eng);
    x = c - log(1+unif(eng)*f); // sample using Rayleigh
    // Reject
    if (tmp*tmp*x <= c){ // accepted
      accept = false;
    }
  }
  // this Rayleigh transform can be delayed till the end
  return sqrt(2.0*x);
}

double tn(double l, double u, sitmo::threefry& eng, double tol = 2.05){
  // samples a column vector of length=length(l)=length(u)
  // from the standard multivariate normal distribution,
  // truncated over the region [l,u], where -a<l<u<a for some
  // 'a' and l and u are column vectors;
  // uses acceptance rejection and inverse-transform method;
  // tol=2.05 # controls switch between methods
  // # threshold can be tuned for maximum speed for each platform
  // case: abs(u-l)>tol, uses accept-reject from randn
  std::normal_distribution<> rnorm(0.0, 1.0);
  std::uniform_real_distribution<> runif(0.0, 1.0);
  double x=0;
  if (fabs(u - l) > tol){
    x = rnorm(eng);
    while (x < l || x > u){
	x = rnorm(eng);
    }
  } else {
    double pl=R::pnorm(l, 0.0, 1.0, true, false);
    double pu=R::pnorm(u, 0.0, 1.0, true, false);
    x = R::qnorm(pl+(pu-pl)*runif(eng), 0.0, 1.0, true, false);
  }
  return x;
}

double trandn(double l, double u, sitmo::threefry& eng, double a=0.4,
	      double tol = 2.05){
  // truncated normal generator
  // * efficient generator of a vector of length(l)=length(u)
  // from the standard multivariate normal distribution,
  // truncated over the region [l,u];
  // infinite values for 'u' and 'l' are accepted;
  // * Remark:
  // If you wish to simulate a random variable
  // 'Z' from the non-standard Gaussian N(m,s^2)
  // conditional on l<Z<u, then first simulate
  // X=trandn((l-m)/s,(u-m)/s) and set Z=m+s*X;
  // 
  // Reference:
  // Z. I. Botev (2015),
  // "The Normal Law Under Linear Restrictions:
  //  Simulation and Estimation via Minimax Tilting", submitted to JRSS(B)
  // a=.4; # treshold for switching between methods
  // threshold can be tuned for maximum speed for each Matlab version
  // three cases to consider:
  double x=0;
  if (l > a){
    // case 1: a<l<u
    x = ntail(l, u, eng);
  } else if (u < -a){
      // case 2: l<u<-a
    x = -ntail(-u, -l, eng);
  } else {
    x = tn(l, u, eng, tol);
  }
  return x;
}

double lnNpr(double a,double b) {
  // computes ln(P(a<Z<b))
  // where Z~N(0,1) very accurately for any 'a', 'b'
  if (a > 0){
    double pa = R::pnorm(a, 0.0, 1.0, false, true);
    double pb = R::pnorm(b, 0.0, 1.0, false, true);
    return pa + log1p(-exp(pb-pa));
  }
  // case a < b < 0
  if (b < 0){
    double pa = R::pnorm(a, 0.0, 1.0, true, true);
    double pb = R::pnorm(b, 0.0, 1.0, true, true);
    return pb+log1p(-exp(pa-pb));
  }
  // case a<0<b
  double pa = R::pnorm(a, 0.0, 1.0, true, false);
  double pb = R::pnorm(b, 0.0, 1.0, false, false);
  return log1p(-pa-pb);
}

typedef struct {
  arma::mat Z;
  arma::vec p;
  arma::vec u;
} rx_mvnrnd;

rx_mvnrnd mvnrnd(int n, arma::mat& L, arma::vec& l,
		 arma::vec& u, arma::vec mu,
		 sitmo::threefry& eng,
		 double a=0.4, double tol = 2.05){
  // generates the proposals from the exponentially tilted 
  // sequential importance sampling pdf;
  // output:    'logpr', log-likelihood of sample
  //              Z, random sample
  rx_mvnrnd ret;
  int d=l.n_elem; // Initialization
  mu[d-1]=0;
  arma::mat Z(d,n,arma::fill::zeros); // create array for variables
  arma::vec p(n, arma::fill::zeros);
  arma::vec uu(n, arma::fill::zeros);
  std::uniform_real_distribution<> unif(0.0, 1.0);
  for (int k = 0; k < d; ++k){
    //# compute matrix multiplication L*Z
    arma::vec col=trans(L(k,arma::span(0,k)) * Z.rows(0, k));
    //# compute limits of truncation
    arma::vec tl=l[k]-mu[k]-col;
    arma::vec tu=u[k]-mu[k]-col;
    //#simulate N(mu,1) conditional on [tl,tu]
    for (int j = n; j--;){
      Z(k,j) = mu[k] + trandn(tl[j], tu[j], eng, a, tol);
      // # update likelihood ratio
      p[j] += lnNpr(tl[j], tu[j]) + 0.5*mu[k]*mu[k] - mu[k]*Z(k,j);
      if (k == 0) uu[j] = -log(unif(eng));
    }
  }
  ret.Z = Z;
  ret.p = p;
  ret.u = uu;
  return ret;
}

//[[Rcpp::export]]
List rxMvnrnd(int n, arma::mat& L, arma::vec& l,
	      arma::vec& u, arma::vec mu,
	      double a=0.4, double tol = 2.05){
  double seedD = runif(1, 1.0, std::numeric_limits<uint32_t>::max())[0];
  uint32_t seed = static_cast<uint32_t>(seedD);
  sitmo::threefry eng;
  eng.seed(seed);
  rx_mvnrnd retI = mvnrnd(n, L, l, u, mu, eng,
			  a, tol);
  List ret(2);
  NumericVector po(retI.p.size());
  std::copy(retI.p.begin(), retI.p.end(), po.begin());
  ret[0] = po;
  ret[1] = wrap(retI.Z);
  ret.attr("names") = CharacterVector::create("logpr", "Z");
  return ret;
}



typedef struct {
  // return(list(L=L,l=l,u=u,perm=perm))
  arma::mat L;
  arma::vec l;
  arma::vec u;
  arma::uvec perm;
  
} rx_cholperms;

rx_cholperms cholperm(arma::mat Sig, arma::vec& cl, arma::vec& cu,
		      double eps=1e-10){
  // #  Computes permuted lower Cholesky factor L for Sig
  // #  by permuting integration limit vectors l and u.
  // #  Outputs perm, such that Sig(perm,perm)=L%*%t(L).
  // #
  // # Reference: 
  // #  Gibson G. J., Glasbey C. A., Elston D. A. (1994),
  // #  "Monte Carlo evaluation of multivariate normal integrals and
  // #  sensitivity to variate ordering", 
  // #  In: Advances in Numerical Methods and Applications, pages 120--126
  arma::vec l=cl;
  arma::vec u=cu;
  int d = l.n_elem;
  arma::uvec perm(d);
  std::iota(perm.begin(),perm.end(),0);
  arma::mat L(d,d, arma::fill::zeros);
  arma::vec z(d, arma::fill::zeros);
  for (int j = 0; j < d; ++j){
    arma::vec pr(d);
    std::fill_n(pr.begin(), d, R_PosInf);
    arma::span I = span(j,d-1);
    arma::vec D = Sig.diag();
    arma::vec sv(1);
    if (j > 1){
      arma::vec r(j);
      std::fill_n(r.begin(), j, 1.0);
      sv = D(I) - (L(I, span(0, j-1)) % L(I, span(0, j-1))) * r;
    } else if (j == 1){
      sv =D(I)-L(I,0) % L(I,0);
    } else {
      sv = D(I);
    }
    for(int kk=sv.size(); kk--;){
      if (sv[kk] < 0) sv[kk] = sqrt(eps);
      else sv[kk] = sqrt(sv[kk]);
    }
    arma::vec colsV;
    if (j >1){
      colsV=L(I,span(0, j-1)) * z(span(0,j-1));
    } else if (j==1){
      colsV=L(I,0)*z(0);
    }
    arma::vec tl;
    arma::vec tu;
    if (j == 0){
      tl = l(I)/sv;
      tu = u(I)/sv;
    } else {
      tl = (l(I) - colsV)/sv;
      tu = (u(I) - colsV)/sv;
    }
    double minPr=R_PosInf;
    int k=0;
    for (int kk = j; kk < d; ++kk){
      pr[kk] = lnNpr(tl[kk-j], tu[kk-j]);
      if (pr[kk] < minPr) {
	minPr = pr[kk];
	k = kk;
      }
    }
    // find smallest marginal dimension
    // flip dimensions k-->j
    arma::uvec jk(2); jk(0) = j; jk(1) = k;
    arma::uvec kj(2); kj(0) = k; kj(1) = j;
    // update rows and cols of Sig
    Sig.rows(jk)=Sig.rows(kj);
    Sig.cols(jk)=Sig.cols(kj);
    // Update only rows of L.
    L.rows(jk) = L.rows(kj);
    // keep track of permutation
    perm(jk) = perm(kj);
    // update integration limits
    l(jk) = l(kj); u(jk) = u(kj);
    //construct L sequentially via Cholesky computation
    if (j == 0){
      sv=Sig(j,j)-L(j,0) * L(j,0);
    } else {
      sv=Sig(j,j)-accu(L(j,span(0,j-1)) % L(j,span(0,j-1)));
    }
    if (sv(0) < -0.001){
      stop("'Sigma' is not positive semi-definite");
    }
    if (sv(0) < 0) sv(0) =eps;
    L(j,j)= sqrt(sv(0));
    if (j < d-1 ){
      if (j > 1){
	L(span(j+1,d-1), j) = (Sig(span(j+1,d-1),j)-
			       L(span(j+1,d-1),span(0,j-1)) *
			       trans(L(j,span(0,(j-1))))) / L(j,j);
      } else if (j == 1){
	L(span(j+1,d-1),j) = (Sig(span(j+1,d-1),j)-
					L(span(j+1,d-1),0) *
					L(j,0)) / L(j,j);
      } else if (j == 0){
	L(span(j+1,d-1),j)=Sig(span((j+1),d-1),j)/L(j,j);
      }
    }
    // find mean value, z(j), of truncated normal:
    tl = (l(j) - L(j, span(0, j)) * z(span(0, j))) / L(j,j);
    tu = (u(j) - L(j, span(0, j)) * z(span(0, j))) / L(j, j);
    double w = lnNpr(tl(0), tu(0));
    // 1/sqrt(2*pi) = M_1_SQRT_2PI
    z[j] = (exp(-0.5 * tl(0)*tl(0) - w) - exp(-0.5* tu(0) * tu(0) - w)) * M_1_SQRT_2PI;
    // Rprintf("tl: %f, tu: %f, w: %f, z: %f\n", tl(0), tu(0), w, z[j]);
  }
  rx_cholperms ret;
  ret.L = L;
  ret.l = l;
  ret.u = u;
  ret.perm=perm;
  return ret;
}
// Exported for testing
//[[Rcpp::export]]
List rxCholperm(arma::mat Sig, arma::vec l, arma::vec u,
		double eps=1e-10){
  rx_cholperms retI = cholperm(Sig, l, u, eps);
  List ret(4);
  NumericVector lI(retI.l.size());
  std::copy(retI.l.begin(),retI.l.end(), lI.begin());
  NumericVector uI(retI.u.size());
  std::copy(retI.u.begin(),retI.u.end(), uI.begin());
  IntegerVector permI(retI.perm.size());
  std::copy(retI.perm.begin(),retI.perm.end(), permI.begin());
  ret[0] = wrap(retI.L);
  ret[1] = lI;
  ret[2] = uI;
  ret[3] = permI;
  ret.attr("names") = CharacterVector::create("L", "l", "u", "perm");
  return ret;
}

typedef struct {
  arma::vec grad;
  arma::mat Jac;
} rx_gradpsi;

rx_gradpsi gradpsi(arma::vec y, arma::mat L, arma::vec l, arma::vec u,
		   int ncores=1){
  rx_gradpsi ret;
  //# implements grad_psi(x) to find optimal exponential twisting;
  //  # assume scaled 'L' with zero diagonal;
  int d = u.n_elem;
  arma::vec c(d,arma::fill::zeros);
  arma::vec x, mu;
  x = c; mu = x;
  if (d <= 1){
    stop(_("dimension wrong in 'gradpsi' (d=%d)"), d);
  }
  x(span(0, d-2)) = y(span(0, d-2));
  mu(span(0, d-2)) = y(span(d-1, 2*d-3));
  //
  // compute now ~l and ~u
  c(span(1, d-1)) = L(span(1,d-1), span(0, d-1)) * x;
  arma::vec lt = l - mu - c;
  arma::vec ut = u - mu - c;
  // compute gradients avoiding catastrophic cancellation
  arma::vec w(d);
  arma::vec pl(d);
  arma::vec pu(d);
  arma::vec P(d);
#ifdef _OPENMP
#pragma omp parallel num_threads(ncores) if(ncores > 1)
  {
#pragma omp for schedule(static)
#endif
  for (int j = 0; j < d; ++j){
    w[j] = lnNpr(lt[j], ut[j]);
    pl[j] = exp(-0.5*lt[j]*lt[j] - w[j])*M_1_SQRT_2PI;
    pu[j] = exp(-0.5*ut[j]*ut[j] - w[j])*M_1_SQRT_2PI;
    P[j] = pl[j] - pu[j];
  }
  arma::vec dfdx = -mu(span(0, d-2)) + trans(trans(P) * L(span(0, d-1), span(0, d-2)));
  arma::vec dfdm = mu - x + P;
  arma::vec grad(dfdx.size()+dfdm.size()-1);
  std::copy(dfdx.begin(),dfdx.end(), grad.begin());
  std::copy(dfdm.begin(),dfdm.end()-1, grad.begin()+dfdx.size());
  // here compute Jacobian matrix
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (int j = 0; j < d; ++j){
    if (!R_FINITE(lt[j])){
      lt[j] = 0;
    }
    if (!R_FINITE(ut[j])){
      ut[j] = 0;
    }
  }
  arma::vec dP = -(P % P) + lt % pl - ut % pu; // dPdm
  arma::mat dPm(d, d);
  for (int j = d; j--;){
    std::copy(dP.begin(),dP.end(),dPm.begin()+d*j);
  }
  arma::mat DL = dPm % L;
  arma::mat mx = -eye(d, d) + DL;
  arma::mat xx = trans(L) * DL;
  arma::mat Jac;
  if (d > 2){
    mx = mx(span(0, d-2), span(0, d-2));
    xx = xx(span(0, d-2), span(0, d-2));
    Jac = join_cols(join_rows(xx,trans(mx)),
		    join_rows(mx, diagmat(1.0+dP(span(0, d-2)))));
  } else {
    mx = mx(0, 0);
    xx = xx(0, 0);
    Jac  = mat(2,2);
    Jac(0,0) = xx(0,0);
    Jac(1,0) = mx(0,0);
    Jac(0,1) = mx(0,0);
    Jac(1,1) = 1+dP(0);
  }
  ret.grad = grad;
  ret.Jac = Jac;
#ifdef _OPENMP
  }
#endif
  return ret;
}


// Exported for testing
//[[Rcpp::export]]
List rxGradpsi(arma::vec y, arma::mat L, arma::vec l, arma::vec u){
  rx_gradpsi retI = gradpsi(y, L, l, u);
  List ret(2);
  NumericVector lgrad(retI.grad.size());
  std::copy(retI.grad.begin(),retI.grad.end(), lgrad.begin());
  ret[1] = wrap(retI.Jac);
  ret[0] = lgrad;
  ret.attr("names") = CharacterVector::create("grad", "Jac");
  return ret;
}


arma::vec nleq(arma::vec l, arma::vec u, arma::mat L, double tol= 1e-10, int maxiter=100){
  int d=l.n_elem;
  arma::vec x(2*d-2, arma::fill::zeros); // initial point for Newton iteration
  double err=R_PosInf;
  int iter=0;
  while (err > tol){
    rx_gradpsi f=gradpsi(x,L,l,u);
    // Jac=f$Jac
    // grad=f$grad
    arma::vec del=solve(f.Jac,-f.grad, solve_opts::likely_sympd); // Newton correction
    x=x+del;
    err=accu(f.grad%f.grad);
    iter++;
    if (iter> maxiter){
      stop(_("covariance matrix is ill-conditioned and method failed (truncated mvn)"));
    }
  }
  return x;
}

// Exported for testing
//[[Rcpp::export]]
NumericVector rxNleq(arma::vec l, arma::vec u, arma::mat L){
  arma::vec retA = nleq(l, u, L, 1e-10, 100);
  NumericVector ret(retA.n_elem);
  std::copy(retA.begin(),retA.end(),ret.begin());
  return ret;
}

double psy(arma::vec x,arma::mat L,arma::vec l, arma::mat u, arma::vec mu,
	   int ncores = 1){
  double p=0;
  // implements psi(x,mu); assume scaled 'L' without diagonal;
  int d=u.n_elem;
  x.resize(d);
  x(d-1)=0;
  mu.resize(d);
  mu(d-1)=0;
  // compute now ~l and ~u
  arma::vec c = L*x;
  l=l-mu-c;
  u=u-mu-c;
#ifdef _OPENMP
#pragma omp parallel num_threads(ncores) if(ncores > 1)
  {
#pragma omp for schedule(static)
#endif
  for (int j = 0; j < d; ++j){
    p+= lnNpr(l[j], u[j]) + 0.5*mu[j]*mu[j]-x[j]*mu[j];
  }
#ifdef _OPENMP
  }
#endif
  return p;
}
arma::mat mvrandn(arma::vec lin, arma::vec uin, arma::mat Sig, int n,
		  sitmo::threefry& eng, double a=0.4, double tol = 2.05,
		  double nlTol=1e-10, int nlMaxiter=100, int ncores=1){
  if (ncores < 1) stop(_("'ncores' has to be greater than one"));
  int d = lin.n_elem;
  if ((int)uin.n_elem != d) stop(_("'lower' and 'upper' must have the same number of elements."));
  if ((int)Sig.n_rows != d || (int)Sig.n_cols != d) stop(_("'sigma' must be a square matrix with the same dimension as 'upper' and 'lower'"));
  if (any(lin>uin)){
    stop(_("'lower' is bigger than 'upper' for at least one item"));
  }
  // Cholesky decomposition of matrix
  rx_cholperms out=cholperm(Sig,lin,uin);
  arma::mat Lfull=out.L;
  arma::vec l=out.l;
  arma::vec u=out.u;
  arma::vec D=Lfull.diag();
  arma::uvec perm=out.perm;
  if (any(D < 1e-10)){
    warning(_("truncated multivariate normal may fail as covariance matrix is singular"));
  }
  // rescale
  arma::mat L=Lfull.each_col()/D;
  u=u/D;
  l=l/D; 
  L=L-eye(d, d); // remove diagonal
  // find optimal tilting parameter via non-linear equation solver
  arma::vec xmu = nleq(l,u,L,nlTol, nlMaxiter); // nonlinear equation solver
  // assign saddlepoint x* and mu*
  arma::vec x = xmu(span(0, d-2));
  arma::vec mu = xmu(span(d-1, 2*d-3));
  // compute psi star
  double psistar=psy(x,L,l,u,mu, ncores);
  //  start acceptance rejection sampling
  int iter=0;
  // rv=c();
  int accepted=0;

  arma::mat ret(d, n);
  while (accepted < n){
    // rx_mvnrnd out=mvnrnd(n,L,l,u,mu);
    rx_mvnrnd out = mvnrnd(n, L, l, u, mu, eng, a, tol);
    arma::vec logpr = out.p;
    arma::mat curZ  = out.Z;
    // idx=-log(runif(n))>(psistar-logpr); # acceptance tests
    for (int i = n; i--;){
      if (out.u[i] > (psistar-logpr[i])){
	ret.col(accepted) =curZ.col(i);
	accepted++;
	if (accepted == n) break;
      }
    }
    iter++;
    if (iter == 1e3){
      warning(_("acceptance probability smaller than 0.001"));
    } else if (iter> 1e4){
      if (accepted == 0) {
	stop(_("could not sample from truncated normal"));
      } else if (accepted > 1) {
	warning(_("sample of size %d which is smaller than requested 'n' returned"), accepted);
	ret = ret.rows(0, accepted);
	break;
      }
    }
  }
  ret = trans(ret);
  ret = ret.cols(out.perm);
  return ret;
}

//[[Rcpp::export]]
arma::mat rxMvrandn_(NumericMatrix A_,
		     arma::rowvec mu, arma::mat sigma, arma::vec lower,
		     arma::vec upper, int ncores=1,
		     double a=0.4, double tol = 2.05,
		     double nlTol=1e-10, int nlMaxiter=100){
  int n = A_.nrow();
  int d = mu.n_elem;
  arma::mat ch;
  if (n < 1) stop(_("n should be a positive integer"));
  if (ncores < 1) stop(_("'ncores' has to be greater than one"));
  if (d != (int)sigma.n_cols) stop("length(mu) != ncol(sigma)");
  if (d != (int)sigma.n_rows) stop("length(mu) != ncol(sigma)");
  if (d != (int)A_.ncol()) stop("length(mu) != ncol(A)");

  double seedD = runif(1, 1.0, std::numeric_limits<uint32_t>::max())[0];
  uint32_t seed = static_cast<uint32_t>(seedD);
  seed = min2(seed, std::numeric_limits<uint32_t>::max() - ncores - 1);
  sitmo::threefry eng;
  arma::mat A(A_.begin(), A_.nrow(), A_.ncol(), false, true);

  arma::vec low = lower-trans(mu);
  arma::vec up = upper-trans(mu);
  if (d == 1){
    double sd = sqrt(sigma(0,0));
    double l=low(0)/sd;
    double u=up(0)/sd;
    for (int i = 0; i < n; ++i){
      A[i] = sd*trandn(l, u, eng, a, tol)+mu(0);
    }
  } else {
    arma::mat ret = mvrandn(low, up, sigma, n, eng, a, tol,
			    nlTol, nlMaxiter, ncores);
    ret.each_row() += mu;
    std::copy(ret.begin(), ret.end(), A.begin());
  }
  return A;
}


sitmo::threefry _eng;

extern "C" void seedEng(int ncores){
  double seedD = runif(1, 1.0, std::numeric_limits<uint32_t>::max())[0];
  uint32_t seed = static_cast<uint32_t>(seedD);
  seed = min2(seed, std::numeric_limits<uint32_t>::max() - ncores - 1);
  _eng.seed(seed);
}

extern "C" int rxbinom(int n, double prob){
  std::binomial_distribution<int> d(prob);
  return d(_eng);
}

extern "C" double rxcauchy(double location, double scale){
  std::cauchy_distribution<double> d(location, scale);
  return d(_eng);
}

extern "C" double rxchisq(double df){
  // Non central not supported in C++11
  std::chi_squared_distribution<double> d(df);
  return d(_eng);
}

extern "C" double rxexp(double rate){
  std::exponential_distribution<double> d(rate);
  return d(_eng);
}

extern "C" double rxf(double df1, double df2){
  std::fisher_f_distribution<double> d(df1, df2);
  return d(_eng);
}

extern "C" double rxgamma(double shape, double rate){
  std::gamma_distribution<double> d(shape, rate);
  return d(_eng);
}

extern "C" double rxbeta(double shape1, double shape2){
  // Efficient simulation when shape1 and shape2 are "large"
  // (p 658) Intro Prob Stats 8th Ed by Sheldon Ross
  double x = rxgamma(shape1,1.0);
  return x/(x+rxgamma(shape2, 1.0));
}

extern "C" int rxgeom(double prob){
  std::geometric_distribution<int> d(prob);
  return d(_eng);
}

extern "C" double rxlnorm(double meanlog, double sdlog){
  std::lognormal_distribution<double> d(meanlog, sdlog);
  return d(_eng);
}

// FIXME rnbinom

extern "C" double rxnorm(double mean, double sd){
  std::normal_distribution<double> d(mean, sd);
  return d(_eng);
}

extern "C" int rxpois(double lambda){
  std::poisson_distribution<int> d(lambda);
  return d(_eng);
}

extern "C" double rxt_(double df){
  std::student_t_distribution<double> d(df);
  return d(_eng);
}

extern "C" double rxunif(double low, double hi){
  std::uniform_real_distribution<double> d(low, hi);
  return d(_eng);
}

extern "C" double rxweibull(double shape, double scale){
  std::weibull_distribution<double> d(shape, scale);
  return d(_eng);
}
