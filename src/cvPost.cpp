#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
bool rxIs(const RObject &obj, std::string cls);

arma::mat rwish5(double nu, int p){
  // GetRNGstate();
  arma::mat Z(p,p, fill::zeros);
  double curp = nu;
  double tmp =sqrt(Rf_rchisq(curp--));
  Z(0,0) = (tmp < 1e-100) ? 1e-100 : tmp;
  int i, j;
  if (p > 1){
    for (i = 1; i < (int)p; i++){
      tmp = sqrt(Rf_rchisq(curp--));
      Z(i,i) = (tmp < 1e-100) ? 1e-100 : tmp;
      for (j = 0; j < i; j++){
        // row,col
        Z(j,i) = norm_rand();
      }
    }
  }
  // PutRNGstate();
  return Z;
}

NumericMatrix cvPost0(double nu, NumericMatrix omega, bool omegaIsChol = false,
                      bool returnChol = false){
  arma::mat S =as<arma::mat>(omega);
  int p = S.n_rows;
  if (p == 1){
    // GetRNGstate();
    NumericMatrix ret(1,1);
    if (omegaIsChol){
      ret[0] = nu*omega[0]*omega[0]/(Rf_rgamma(nu/2.0,2.0));
    } else {
      ret[0] = nu*omega[0]/(Rf_rgamma(nu/2.0,2.0));
    }
    if (returnChol) ret[0] = sqrt(ret[0]);
    // PutRNGstate();
    return ret;
  } else {
    arma::mat Z = rwish5(nu, p);
    arma::mat Z2 = arma::trans(arma::solve(trimatu(Z), eye(p, p)));
    arma::mat cv5;
    if (omegaIsChol){
      cv5 = S;
    } else {
      cv5 = arma::chol(S);
    }
    arma::mat mat1 = Z2 * cv5;
    mat1 = mat1.t() * mat1;
    mat1 = mat1 * nu;
    if (returnChol) mat1 = arma::chol(mat1);
    return wrap(mat1);
  }
}

//' Sample a covariance Matrix from the Posterior Inverse Wishart
//' distribution.
//'
//' Note this Inverse wishart rescaled to match the original scale of
//' the covariance matrix.
//'
//' If your covariance matrix is a 1x1 matrix, this uses an scaled
//' inverse chi-squared which is equivalent to the Inverse Wishart
//' distribution in the uni-directional case.
//'
//' @param nu Degrees of Freedom (Number of Observations) for 
//'        covariance matrix simulation.
//' 
//' @param omega Estimate of Covariance matrix.
//' 
//' @param n Number of Matrices to sample.  By default this is 1.
//' 
//' @param omegaIsChol is an indicator of if the omega matrix is in
//'   the Cholesky decomposition.
//' 
//' @param returnChol Return the Cholesky decomposition of the
//'   covariance matrix sample.
//'
//' @return a matrix (n=1) or a list of matrices  (n > 1)
//'
//' @author Matthew L.Fidler & Wenping Wang
//'
//' @examples
//' 
//' ## Sample a single covariance.
//' draw1 <- cvPost(3, matrix(c(1,.3,.3,1),2,2))
//'
//' ## Sample 3 covariances
//' set.seed(42)
//' draw3 <- cvPost(3, matrix(c(1,.3,.3,1),2,2), n=3)
//' 
//' ## Sample 3 covariances, but return the cholesky decomposition
//' set.seed(42)
//' draw3c <- cvPost(3, matrix(c(1,.3,.3,1),2,2), n=3, returnChol=TRUE)
//' @export
//[[Rcpp::export]]
RObject cvPost(double nu, RObject omega, int n = 1, bool omegaIsChol = false, bool returnChol = false){
  if (n == 1){
    if (rxIs(omega,"numeric.matrix") || rxIs(omega,"integer.matrix")){
      return as<RObject>(cvPost0(nu, as<NumericMatrix>(omega), omegaIsChol, returnChol));
    } else if (rxIs(omega, "numeric") || rxIs(omega, "integer")){
      NumericVector om1 = as<NumericVector>(omega);
      if (om1.size() % 2 == 0){
        int n1 = om1.size()/2;
        NumericMatrix om2(n1,n1);
        for (int i = 0; i < om1.size();i++){
          om2[i] = om1[i];
        }
        return as<RObject>(cvPost0(nu, om2, omegaIsChol, returnChol));
      }
    }
  } else {
    List ret(n);
    for (int i = 0; i < n; i++){
      ret[i] = cvPost(nu, omega, 1, omegaIsChol, returnChol);
    }
    return(as<RObject>(ret));
  }
  stop("omega needs to be a matrix or a numberic vector that can be converted to a matrix.");
  return R_NilValue;
}

//' Scaled Inverse Chi Squared distribution
//'
//' @param n Number of random samples
//' @param nu degrees of freedom of inverse chi square
//' @param scale  Scale of inverse chi squared distribution 
//'         (default is 1).
//' @return a vector of inverse chi squared deviates .
//' @examples
//' rinvchisq(3, 4, 1) ## Scale = 1, degrees of freedom = 4
//' rinvchisq(2, 4, 2) ## Scale = 2, degrees of freedom = 4
//' @export
//[[Rcpp::export]]
NumericVector rinvchisq(const int n = 1, const double &nu = 1.0, const double &scale = 1){
  NumericVector ret(n);
  // GetRNGstate();
  for (int i = 0; i < n; i++){
    ret[i] = nu*scale/(Rf_rgamma(nu/2.0,2.0));
  }
  // PutRNGstate();
  return ret;
}

// Adapted from banocc and ported to armadillo for speed.
// https://github.com/biobakery/banocc/blob/master/R/rlkj.R
void rgbeta(int d, double shape, double* out){
  if (std::isinf(shape)) {
    std::fill_n(out, d, 0.0);
  } else if (shape > 0){
    for (int j = d; j--;){
      out[j] = 2.0*Rf_rbeta(shape, shape) - 1.0;
    }
  } else if (shape == 0){
    for (int j = d; j--;){
      out[j] = 2.0*Rf_rbinom(1, 0.5) - 1.0;
    }
  } else {
    stop("'shape' must be non-negative");
  }
}
//' One correlation sample from the LKJ distribution
//'
//' @param d The dimension of the correlation matrix
//' @param eta The scaling parameter of the LKJ distribution.
//'   Must be > 1.  Also related to the degrees of freedom nu.
//'   eta = (nu-1)/2.
//' @param cholesky boolean; If \code{TRUE} return the cholesky
//'   decomposition.
//' @author Matthew Fidler (translated to RcppArmadillo) and Emma Schwager
//' @export
//[[Rcpp::export]]
arma::mat rLKJ1(int d, double eta = 1.0, bool cholesky = false){
  if (d < 2){
    stop("dimension, 'd' of correlation matrix must be > 1");
  }
  if (eta < 1){
    stop("'eta' must be >= 1");
  }
  double alpha = eta + ((double)(d) - 2.0)/2.0;
  arma::mat L(d,d,arma::fill::zeros);
  L(0,0) = 1.0;
  arma::vec partials(d-1);
  rgbeta(d-1, alpha, partials.memptr());
  std::copy(partials.begin(), partials.end(), L.memptr()+1);
  if (d == 2){
    L(1,1) = sqrt(1-L(1,0)*L(1,0));
    if (!cholesky){
      L = L * L.t();
    }
    return L;
  }
  arma::vec W = log(1-partials%partials);
  for (int i = 2; i <= d-1; i++){
    alpha -= 0.5;
    rgbeta(d-i, alpha, partials.memptr());
    // construct a vector pointing to the partials vector without
    // allocating new memory:
    arma::vec partials2 = arma::vec(partials.memptr(), d-i, false, true);
    L(i-1,i-1) = exp(0.5*W(i-2));
    L(arma::span(i,d-1),i-1) = partials2 % exp(0.5*W(arma::span(i-1,d-2)));
    W(arma::span(i-1,d-2)) = W(arma::span(i-1,d-2)) +log(1-partials2%partials2);
  }
  L(d-1,d-1) = exp(0.5*W(d-2));
  if (!cholesky){
    L = L * L.t();
  }
  return L;
}

//[[Rcpp::export]]
arma::mat rLKJcv1(arma::vec sd, double eta = 1.0){
  int d = sd.size();
  arma::mat r = rLKJ1(d, eta, false);
  arma::mat dSd = diagmat(sd);
  return dSd*r*dSd;
}

//[[Rcpp::export]]
arma::mat rLKJcvLsd1(arma::vec logSd, arma::vec logSdSD, double eta = 1.0){
  unsigned int d = logSd.size();
  if (d != logSdSD.size()){
    stop("log standard deviation size needs to be the same size as the log standard error of the estimate");
  }
  arma::vec sd(d);
  for (unsigned int j = d; j--;){
    sd[j] = exp(Rf_rnorm(logSd[j], logSdSD[j]));
  }
  return rLKJcv1(sd, eta);
}

//' One correlation sample from the Inverse Wishart distribution
//'
//' This correlation is constructed by transformation of the Inverse Wishart
//' random covariate to a correlation.
//'
//' @inheritParams rLKJ1
//' @param nu Degrees of freedom of the Wishart distribution
//' @inheritParams cvPost
//' @author Matthew Fidler
//' @export
//[[Rcpp::export]]
arma::mat invWR1d(int d, double nu, bool omegaIsChol = false){
  if (nu <= d - 1) stop("'nu' must be greater than 'd'-1");
  arma::mat I(d,d,arma::fill::eye);
  arma::mat invW = as<arma::mat>(cvPost0(nu, wrap(I),
					 omegaIsChol, false));
  arma::mat Dinv = diagmat(1/sqrt(invW.diag()));
  return Dinv * invW * Dinv;
}

arma::mat rinvWRcv1(arma::vec sd, double nu = 1.0){
  int d = sd.size();
  arma::mat r = invWR1d(nu, d, false);
  arma::mat dSd = diagmat(sd);
  return dSd*r*dSd;
}

//[[Rcpp::export]]
arma::mat rcvC1(arma::vec sdEst, double nu = 3.0,
		int diagXformType = 1, int rType = 1){
  // the sdEst should come from the multivariate normal distribution
  // with the appropriate transformation.
  unsigned int d = sdEst.size();
  // Nlmixr models variances as chol(omega^1)
  // var = diag(omega)
  // Assuming off-diagonals are zero
  // diag(omega^1) = (1/var)
  // chol(diag(omega^1)) = sqrt(1/var) = 1/sd
  // With diagonals this becomes
  // diagXform = c("sqrt", "log", "identity")
  //
  arma::vec sd(d);
  switch(diagXformType){
  case 1:
    // sqrt; In this case we estimate x^2
    // x^2 = 1/sd
    // sd = 1/x^2
    for (int j = d; j--;){
      sd[j] = 1/(sdEst[j]*sdEst[j]);
    }
    break;
  case 2:
    // log
    // In this case we estimate exp(x)
    // exp(x) = 1/sd
    // sd = 1/exp(x)
    for (int j = d; j--;){
      sd[j] = 1/exp(sdEst[j]);
    }
    break;
  case 3:
    // identity
    // In this case we estimate x
    // sd = 1/x
    for (int j = d; j--;){
      sd[j] = 1/sdEst[j];
    }
    break;
  case 4: // direct identity
    for (int j = d; j--;){
      sd[j] = sdEst[j];
    }
    break;
  case 5: // lognormal
    for (int j = d; j--;){
      sd[j] = exp(sdEst[j]);
    }
    break;
  }
  if (rType == 1){
    return rLKJcv1(sd, (nu-1.0)/2.0);
  } else {
    return rinvWRcv1(sd, nu);
  }
}
