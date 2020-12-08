// [[Rcpp::interfaces(r, cpp)]]
//#undef NDEBUG
#define min2( a , b )  ( (a) < (b) ? (a) : (b) )
#include <RcppArmadillo.h>
#include "../inst/include/RxODE.h"
#include <vandercorput.h>
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

sitmo::vandercorput _engV;
void seedEngV(uint32_t seed){
  _engV.seed(seed);
}

extern "C" double rxnormV(rx_solving_options_ind* ind, double mean, double sd){
  if (!ind->inLhs) return 0;
  std::normal_distribution<double> d(mean, sd);
  return d(_engV);
}

extern "C" double rinormV(rx_solving_options_ind* ind, int id, double mean, double sd){
  if (ind->isIni == 1) {
    std::normal_distribution<double> d(mean, sd);
    ind->simIni[id] = d(_engV);
  }
  return ind->simIni[id];
}

//[[Rcpp::export]]
arma::mat rxrandnV(unsigned int nrow, unsigned int ncol){
  arma::mat ret(nrow, ncol);
  std::normal_distribution<double> d(0.0, 1.0);
  for (int j = nrow; j--;) {
    for (int i = ncol; i--;) {
      ret(j,i) = d(_engV);
    }
  }
  return ret;
}


//[[Rcpp::export]]
NumericVector rxnormV_(double mean, double sd, int n, int ncores){
  NumericVector ret(n);
  int n2 = ret.size();
  double seedD = runif(1, 1.0, std::numeric_limits<uint32_t>::max())[0];
  uint32_t seed = static_cast<uint32_t>(seedD);
  seed = min2(seed, std::numeric_limits<uint32_t>::max() - ncores - 1);
  sitmo::vandercorput eng;
  double *A  = ret.begin();
  std::normal_distribution<double> d(mean, sd);
  #ifdef _OPENMP
#pragma omp parallel num_threads(ncores) if(ncores > 1)
  {
    seed += omp_get_thread_num();
#endif
    eng.seed(seed);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int i = 0; i < n2; ++i){
      A[i] = d(eng);
    }
  #ifdef _OPENMP
  }
  #endif
  return ret;
}
