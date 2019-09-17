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
    // Backsolve isn't available in armadillo
    arma::mat Z2 = arma::trans(arma::inv(trimatu(Z)));
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
