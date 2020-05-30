//#undef NDEBUG
#include <R.h>
#include <stan/math.hpp>
#include "../inst/include/RxODE.h"
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

namespace stan {
  namespace math {

    using std::exp;
    using stan::math::exp;
    using std::sqrt;
    using stan::math::sqrt;
    using std::pow;
    using stan::math::pow;
    using std::acos;
    using stan::math::acos;
    using std::cos;
    using stan::math::cos;

    template<class T>
    Eigen::Matrix<T, Eigen::Dynamic, 2>
    micros2macros(const Eigen::Matrix<T, Eigen::Dynamic, 1>& p,
		  const int& ncmt,
		  const int& trans){
      Eigen::Matrix<T, Eigen::Dynamic, 2> g(ncmt,3);
      T btemp, ctemp, dtemp;
#define p1    p[0]
#define v1    p[1]
#define p2    p[2]
#define p3    p[3]
#define p4    p[4]
#define p5    p[5]
#define v     g(0, 0)
#define k     g(0, 1)
#define k12   g(1, 0)
#define k21   g(1, 1)
#define k13   g(2, 0)
#define k31   g(2, 1)
      switch (ncmt) {
      case 3: { // 3 compartment model 
	switch (trans){
	case 1: // cl v q vp
	  k = p1/v1; // k = CL/V
	  v = v1;
	  k12 = p2/v1; // k12 = Q/V
	  k21 = p2/p3; // k21 = Q/Vp
	  k13 = p4/v1; // k31 = Q2/V
	  k31 = p4/p5; // k31 = Q2/Vp2
	  break;
	case 2: // k=(*p1) v=(*v1) k12=(*p2) k21=(*p3) k13=(*p4) k31=(*p5)
	  k = p1;
	  v = v1;
	  k12 = p2;
	  k21 = p3;
	  k13 = p4;
	  k31 = p5;
	  break;
	case 11:
#undef beta
#define A (1/v1)
#define B (p3)
#define C (p5)
#define alpha (p1)
#define beta (p2)
#define gamma (p4)
	  v=1/(A+B+C);
	  btemp = -(alpha*C + alpha*B + gamma*A + gamma*B + beta*A + beta*C)*v;
	  ctemp = (alpha*beta*C + alpha*gamma*B + beta*gamma*A)*v;
	  dtemp = sqrt(btemp*btemp-4*ctemp);
	  k21 = 0.5*(-btemp+dtemp);
	  k31 = 0.5*(-btemp-dtemp);
	  k   = alpha*beta*gamma/k21/k31;
	  k12 = ((beta*gamma + alpha*beta + alpha*gamma) -
		 k21*(alpha+beta+gamma) - k * k31 + k21*k21)/(k31 - k21);
	  k13 = alpha + beta + gamma - (k + k12 + k21 + k31);
	  break;
	case 10:
#undef A
#define A v1
	  v=1/(A+B+C);
	  btemp = -(alpha*C + alpha*B + gamma*A + gamma*B + beta*A + beta*C)*v;
	  ctemp = (alpha*beta*C + alpha*gamma*B + beta*gamma*A)*v;
	  dtemp = sqrt(_as_zero(btemp*btemp-4*ctemp));
	  k21 = 0.5*(-btemp+dtemp);
	  k31 = 0.5*(-btemp-dtemp);
	  k   = alpha*beta*gamma/k21/k31;
	  k12 = ((beta*gamma + alpha*beta + alpha*gamma) -
		 k21*(alpha+beta+gamma) - k * k31 + k21*k21)/(k31 - k21);
	  k13 = alpha + beta + gamma - (k + k12 + k21 + k31);
#undef A
#undef B
#undef C
#undef alpha
#undef beta
#undef gamma
#define beta Rf_beta
	  break;
	}
      } break;
      case 2:{ // 2 compartment model
	switch (trans){
	case 1: // cl=(*p1) v=(*v1) q=(*p2) vp=(*p3)
	  k = p1/v1; // k = CL/V
	  v = v1;
	  k12 = p2/v1; // k12 = Q/V
	  k21 = p2/p3; // k21 = Q/Vp
	  break;
	case 2: // k=(*p1), (*v1)=v k12=(*p2) k21=(*p3)
	  k = p1;
	  v = v1;
	  k12 = p2;
	  k21 = p3;
	  break;
	case 3: // cl=(*p1) v=(*v1) q=(*p2) vss=(*p3)
	  k = p1/v1; // k = CL/V
	  v = v1;
	  k12 = p2/v1; // k12 = Q/V
	  k21 = p2/(p3-v1); // k21 = Q/(Vss-V)
	  break;
	case 4: // alpha=(*p1) beta=(*p2) k21=(*p3)
	  v = v1;
	  k21 = p3;
	  k = p1*p2/k21; // (*p1) = alpha (*p2) = beta
	  k12 = p1 + p2 - k21 - k;
	  break;
	case 5: // alpha=(*p1) beta=(*p2) aob=(*p3)
	  v=v1;
	  k21 = (p3*p2+p1)/(p3+1);
	  k = (p1*p2)/k21;
	  k12 = p1+ p2 - k21 - k;
	  break;
	case 11: // A2 V, alpha=(*p1), beta=(*p2), k21
#undef beta
#define A (1/v1)
#define B (p3)
#define alpha (p1)
#define beta (p2)
	  v   = 1/(A+B);
	  k21 = (A*beta + B*alpha)*v;
	  k   = alpha*beta/k21;
	  k12 = alpha+beta-k21-k;
	  break;
	case 10: // A=(*v1), alpha=(*p1), beta=(*p2), B=(*p3)
	  // Convert to A (right now A=(*v1) or A=1/(*v1))
#undef A
#define A (v1)
	  v   = 1/(A + B);
	  k21 = (A*beta + B*alpha)*v;
	  k   = alpha*beta/k21;
	  k12 = alpha + beta - k21 - k;
#undef A
#undef B
#undef alpha
#undef beta
#define beta Rf_beta
	  break;
	default:
	  REprintf(_("invalid trans (2 cmt trans %d)\n"), trans);
	  return NA_REAL;
	}
      } break;
      case 1:{ // One compartment model
	switch(trans){
	case 1: // cl v
	  k = p1/v1; // k = CL/V
	  v = v1;
	  break;
	case 2: // k V
	  k = p1;
	  v = v1;
	  break;
	case 11: // alpha V
	  k = p1;
	  v = v1;
	  break;
	case 10: // alpha A
	  k = p1;
	  v = 1/v1;
	  break;
	default:
	  return 0;
	}
      } break;
      }
#undef p1
#undef v1
#undef p2
#undef p3
#undef p4
#undef p4
      
#undef k
#undef v
#undef k12
#undef k21
#undef k13
#undef k31
      return g;
    }
    
    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    genericCmtInterface(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
			const double t,
			const int oral0,
			const int trans,
			const int ncmt,
			const int linCmt,
			rx_solving_options_ind *ind){
      Eigen::Matrix<T, Eigen::Dynamic, 2> par(ncmt, 3);
      par = micros2macros(params, ncmt, trans);
    }
  }
}

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> MatrixPd;

extern "C" void sortIfNeeded(rx_solve *rx, rx_solving_options_ind *ind, unsigned int id,
			     int *linCmt,
			     double *d_tlag, double *d_tlag2,
			     double *d_F, double *d_F2,
			     double *d_rate1, double *d_dur1,
			     double *d_rate2, double *d_dur2);

extern "C" double getTime(int idx, rx_solving_options_ind *ind);

extern "C" int _locateTimeIndex(double obs_time,  rx_solving_options_ind *ind);

extern "C" double linCmtB(rx_solve *rx, unsigned int id, double t, int linCmt,
			  int i_cmt, int trans, int val,
			  double dd_p1, double dd_v1,
			  double dd_p2, double dd_p3,
			  double dd_p4, double dd_p5,
			  double dd_tlag, double dd_F, double dd_rate, double dd_dur,
			  // oral extra parameters
			  double dd_ka, double dd_tlag2, double dd_F2, double dd_rate2, double dd_dur2){
  
  // Get  Alast
  rx_solving_options_ind *ind = &(rx->subjects[id]);
  int evid, wh, cmt, wh100, whI, wh0;
  /* evid = ind->evid[ind->ix[ind->idx]]; */
  /* if (evid) REprintf("evid0[%d:%d]: %d; curTime: %f\n", id, ind->idx, evid, t); */
  int idx = ind->idx;
  sortIfNeeded(rx, ind, id, &linCmt, &dd_tlag, &dd_tlag2, &dd_F, &dd_F2,
	       &dd_rate, &dd_dur, &dd_rate2, &dd_dur2);
  rx_solving_options *op = rx->op;
  int oral0;
  oral0 = (dd_ka > 0) ? 1 : 0;
  double *A;
  double *Alast;
  double Alast0[4] = {0, 0, 0, 0};
  /* A = Alast0; Alast=Alast0; */
  double tlast;
  double it = getTime(ind->ix[idx], ind);
  double curTime=0.0;
  
  if (t != it) {
    // Try to get another idx by bisection
    /* REprintf("it pre: %f", it); */
    idx = _locateTimeIndex(t, ind);
    it = getTime(ind->ix[idx], ind);
    /* REprintf("it post: %f", it); */
  }
  /* REprintf("idx: %d; solved: %d; t: %f fabs: %f\n", idx, ind->solved[idx], t, fabs(t-it)); */
  int sameTime = fabs(t-it) < sqrt(DOUBLE_EPS);
  unsigned int ncmt = 1;
  if (dd_p4 > 0.){
    ncmt = 3;
  } else if (dd_p2 > 0.){
    ncmt = 2;
  }
  if (ind->solved[idx] && sameTime){
    // Pull from last solved value (cached)
    A = ind->linCmtAdvan+(op->nlin)*idx;
    // stop("fixme");
    if (val == 0){
      if (trans == 10) {
	return(A[oral0]*dd_v1);
      } else {
	return(A[oral0]/dd_v1);
      }    
    } else {
      return A[ncmt+oral0+val-1];
    }
  }
  MatrixPd params(2*ncmt + 4 + oral0*5, 1);
  params(0,0) = dd_p1;
  params(1,0) = dd_v1;
  if (ncmt >=2 ){
    params(2,0) = dd_p2;
    params(3,0) = dd_p3;
    if (ncmt >= 3){
      params(4,0) = dd_p4;
      params(5,0) = dd_p5;
    }
  }
  if (oral0) {
    params(2*ncmt,     0) = dd_ka;
    params(2*ncmt + 1, 0) = dd_tlag;
    params(2*ncmt + 2, 0) = dd_tlag2;
    params(2*ncmt + 3, 0) = dd_F;
    params(2*ncmt + 4, 0) = dd_F2;
    params(2*ncmt + 5, 0) = dd_rate;
    params(2*ncmt + 6, 0) = dd_dur;
    params(2*ncmt + 7, 0) = dd_rate2;
    params(2*ncmt + 8, 0) = dd_dur2;
  } else {
    params(2*ncmt,     0) = dd_tlag;
    params(2*ncmt + 1, 0) = dd_F;
    params(2*ncmt + 2, 0) = dd_rate;
    params(2*ncmt + 3, 0) = dd_dur;
  }
  
  
  return 0.0;
}
