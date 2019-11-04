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
extern "C" int _locateDoseIndex(const double obs_time,  rx_solving_options_ind *ind);
extern "C" void getWh(int evid, int *wh, int *cmt, int *wh100, int *whI, int *wh0);
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

    //micro to macro conversion
    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 2>
    micros2macros(const Eigen::Matrix<T, Eigen::Dynamic, 1>& p,
		  const int& ncmt,
		  const int& oral,
		  const int& trans,
		  T ka){
      Eigen::Matrix<T, Eigen::Dynamic, 2> g(ncmt,3);
#define A     g(0,0)
#define A2    g(0,1)
#define ALPHA g(0,2)
#define B     g(1,0)
#define B2    g(1,1)
#define BETA  g(1,2)
#define C     g(2,0)
#define C2    g(2,1)
#define GAMMA g(2,2)
#define P1    p[0]
#define V1    p[1]
#define P2    p[2]
#define P3    p[3]
#define P4    p[4]
#define P5    p[5]
      if (trans >= 10){
	ALPHA = P1; // alpha
	// A
	if (trans == 11){
	  A = 1/(V1+DOUBLE_EPS);
	} else {
	  A = V1;
	}
	A2 = 0;// A2
	if (ncmt == 2 || ncmt == 3){
	  BETA = P2; // beta
	  B = P3; // B
	  B2 = 0; // B2
	  if (ncmt == 3){
	    GAMMA = P4; //gamma
	    C = P5; //C
	    C2 = 0; // C2
	  }
	}
	if (oral == 1){
	  if (ncmt == 3){
	    A2 = A; // A2=A
	    B2 = B; // B2=B
	    C2 = C; // C2=C
	    T tmp = ka - ALPHA +DOUBLE_EPS;
	    A = ka / tmp * A;
	    tmp = ka - BETA + DOUBLE_EPS;
	    B = ka / tmp * B;
	    tmp = ka - GAMMA + DOUBLE_EPS;
	    C = ka / tmp * C;
	  } else if (ncmt == 2){
	    A2 = A;
	    B2 = B;
	    T tmp = ka - ALPHA + DOUBLE_EPS;
	    A = ka / tmp * A;
	    tmp = ka - BETA + DOUBLE_EPS;
	    B = ka / tmp * B;
	  } else {
	    A2 = A;
	    T tmp = ka - ALPHA + DOUBLE_EPS;
	    A = ka / tmp * A;
	  }
	}
      } else {
	T rx_k, rx_v, rx_k12, rx_k21, rx_k13, rx_k31;
	if (ncmt == 1){
	  switch(trans){
	  case 1: // cl v
	    rx_k = P1/(V1+DOUBLE_EPS); // k = CL/V
	    rx_v = V1;
	    break;
	  case 2: // k V
	    rx_k = P1;
	    rx_v = V1;
	    break;
	  default:
	    throw std::invalid_argument("invalid trans.");
	  }
	  if (ka > 0){
	    ALPHA = rx_k;
	    T tmp = ka - ALPHA + DOUBLE_EPS;
	    A = ka / tmp / (rx_v+DOUBLE_EPS);
	    A2 = 1.0 / (rx_v+DOUBLE_EPS);
	  } else {
	    ALPHA = rx_k;
	    A = 1.0 / (rx_v+DOUBLE_EPS);
	  }
	} else if (ncmt == 2){
	  switch (trans){
	  case 1: // cl=p1 v=v1 q=p2 vp=p3
	    rx_k = P1/(V1+DOUBLE_EPS); // k = CL/V
	    rx_v = V1;
	    rx_k12 = P2/(V1+DOUBLE_EPS); // k12 = Q/V
	    rx_k21 = P2/(P3+DOUBLE_EPS); // k21 = Q/Vp
	    break;
	  case 2: // k=p1, v1=v k12=p2 k21=p3
	    rx_k = P1;
	    rx_v = V1;
	    rx_k12 = P2;
	    rx_k21 = P3;
	    break;
	  case 3: // cl=p1 v=v1 q=p2 vss=p3
	    rx_k = P1/(V1+DOUBLE_EPS); // k = CL/V
	    rx_v = V1;
	    rx_k12 = P2/(V1+DOUBLE_EPS); // k12 = Q/V
	    rx_k21 = P2/(P3-V1+DOUBLE_EPS); // k21 = Q/(Vss-V)
	    break;
	  case 4: // alpha=p1 beta=p2 k21=p3
	    rx_v = V1;
	    rx_k21 = P3;
	    rx_k = P1*P2/(rx_k21+DOUBLE_EPS); // p1 = alpha p2 = beta
	    rx_k12 = P1 + P2 - rx_k21 - rx_k;
	    break;
	  case 5: // alpha=p1 beta=p2 aob=p3
	    rx_v=V1;
	    rx_k21 = (P3*P2+P1)/(P3+1+DOUBLE_EPS);
	    rx_k = (P1*P2)/(rx_k21+DOUBLE_EPS);
	    rx_k12 = P1+P2 - rx_k21 - rx_k;
	    break;
	  default:
	    throw std::invalid_argument("invalid trans (2cmt).");
	  }
	  T rx_tmp = rx_k12+rx_k21+rx_k;
	  BETA = 0.5 * (rx_tmp - pow(rx_tmp * rx_tmp - 4.0 * rx_k21 * rx_k, 0.5));
	  ALPHA = rx_k21 * rx_k / (BETA+DOUBLE_EPS);
	  T AmB = ALPHA - BETA + DOUBLE_EPS;
	  T BmA = -AmB;
	  if (oral == 1){
	    T kmA = ka - ALPHA+ DOUBLE_EPS;
	    T kmB = ka - BETA + DOUBLE_EPS;
	    T tmpV = rx_v+DOUBLE_EPS;
	    A = ka / (kmA) * (ALPHA - rx_k21) / (AmB) / tmpV;
	    B = ka / (kmB) * (BETA - rx_k21) / (BmA) /tmpV;
	    A2 = (ALPHA - rx_k21) / (AmB) /tmpV;
	    B2 = (BETA - rx_k21) / (BmA) /tmpV;
	  } else {
	    T tmpV = rx_v+DOUBLE_EPS;
	    A = (ALPHA - rx_k21) / (AmB) / tmpV;
	    B = (BETA - rx_k21) / (BmA) / tmpV;
	    A2 = 0;
	    B2 = 0;
	  }
	} else if (ncmt == 3){
	  switch (trans){
	  case 1: // cl v q vp
	    rx_k = P1/(V1+DOUBLE_EPS); // k = CL/V
	    rx_v = V1;
	    rx_k12 = P2/(V1+DOUBLE_EPS); // k12 = Q/V
	    rx_k21 = P2/(P3+DOUBLE_EPS); // k21 = Q/Vp
	    rx_k13 = P4/(V1+DOUBLE_EPS); // k31 = Q2/V
	    rx_k31 = P4/(P5+DOUBLE_EPS); // k31 = Q2/Vp2
	    break;
	  case 2: // k=p1 v=v1 k12=p2 k21=p3 k13=p4 k31=p5
	    rx_k = P1;
	    rx_v = V1;
	    rx_k12 = P2;
	    rx_k21 = P3;
	    rx_k13 = P4;
	    rx_k31 = P5;
	    break;
	  default:
	    throw std::invalid_argument("invalid trans (3cmt).");
	  }
	  T rx_a0, rx_a1, rx_a2, rx_p, rx_q, rx_r1, rx_r2, rx_theta;
	  rx_a0 = rx_k * rx_k21 * rx_k31;
	  rx_a1 = rx_k * rx_k31 + rx_k21 * rx_k31 + rx_k21 * rx_k13 + rx_k * rx_k21 + rx_k31 * rx_k12;
	  rx_a2 = rx_k + rx_k12 + rx_k13 + rx_k21 + rx_k31;
	  rx_p = rx_a1 - rx_a2 * rx_a2 / 3.0;
	  rx_q = 2.0 * rx_a2 * rx_a2 * rx_a2 / 27.0 - rx_a1 * rx_a2 /3.0 + rx_a0;
	  rx_r1 = pow(-rx_p * rx_p * rx_p / 27.0, 0.5);
	  rx_r2 = 2 * pow(rx_r1,1.0/3.0);
	  rx_theta = acos(-rx_q / (2.0 * rx_r1 + DOUBLE_EPS)) / 3.0;
	  ALPHA = -(cos(rx_theta) * rx_r2 - rx_a2 / 3.0);
	  BETA = -(cos(rx_theta + M_2PI/3.0) * rx_r2 - rx_a2 / 3.0);
	  GAMMA = -(cos(rx_theta + 4.0 / 3.0 * M_PI) * rx_r2 - rx_a2 / 3.0);
	  T AmB = ALPHA - BETA + DOUBLE_EPS;
	  T BmA = -AmB;
	  T AmG = ALPHA - GAMMA + DOUBLE_EPS;
	  T GmA = -AmG;
	  T BmG = BETA - GAMMA + DOUBLE_EPS;
	  T GmB = -BmG;
	  T tmpV = rx_v+DOUBLE_EPS;
	  A = (rx_k21 - ALPHA) * (rx_k31 - ALPHA) / (AmB) / (AmG) / tmpV;
	  B = (rx_k21 - BETA) * (rx_k31 - BETA) / (BmA) / (BmG) / tmpV;
	  C = (rx_k21 - GAMMA) * (rx_k31 - GAMMA) / (GmA) / (GmB) / tmpV;
	  if (oral == 1){
	    A2 = A;
	    B2 = B;
	    C2 = C;
	    T kmA = ka - ALPHA + DOUBLE_EPS;
	    T kmB = ka - BETA + DOUBLE_EPS;
	    T kmG = ka - GAMMA + DOUBLE_EPS;
	    A = ka / (kmA) * A;
	    B = ka / (kmB) * B;
	    C = ka / (kmG) * C;
	  }
	} else {
	  throw std::invalid_argument("1-3 compartments are supported");
	}
      }
#undef ALPHA
#undef A
#undef A2
#undef BETA
#undef B
#undef B2
#undef GAMMA
#undef C
#undef C2
#undef P1
#undef V1
#undef P2
#undef P3
#undef P4
#undef P5
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
      Eigen::Matrix<T, Eigen::Dynamic, 2> FF(1,2);
      Eigen::Matrix<T, Eigen::Dynamic, 2> tlag(1,2);
      unsigned int  l = 0, p = 0;
      int evid, wh, wh100, whI, wh0;
      T thisT = 0.0;
      double dose = 0;
      T res, tinf = 0;
      double tau = 0;
      T tT = 0.0;
      T rate=0;
      // After parameters (2*ncmt)
      // ka, tlag1, tlag2, f1, f2, rate1, rate2, dur1, dur2
      T ka	= params[2 * ncmt];
      tlag(0,0)	= params[2 * ncmt + 1];
      tlag(0,1)	= params[2 * ncmt + 2];
      FF(0,0)   = params[2 * ncmt + 3];
      FF(0,1)	= params[2 * ncmt + 4];
      T d_rate	= params[2 * ncmt + 5];
      T d_dur	= params[2 * ncmt + 6];
      Eigen::Matrix<T, Eigen::Dynamic, 2> par(ncmt, 3);
      par = micros2macros(params, ncmt, oral0, trans, ka);
      Eigen::Matrix<T, Eigen::Dynamic, 1> g(1);
      Eigen::Matrix<T, Eigen::Dynamic, 1> obs_time(ind->n_all_times);
      unsigned int shiftF = 0, oral=oral0;
      int cmt;
#define ret g(0, 0)
#define A par(jj,shiftF)
#define ALPHA par(jj,2)
#define F FF(0,shiftF)
#define TLAG tlag(0,shiftF)
      unsigned int m = _locateDoseIndex(t, ind);
      unsigned int ndoses = (unsigned int)(ind->ndoses);
      ret=0;
      for(l=m+1; l--;){// Optimized for loop as https://www.thegeekstuff.com/2015/01/c-cpp-code-optimization/
	//super-position
	evid = ind->evid[ind->idose[l]];
	if (evid == 3){
	  return g; 
	}
	getWh(evid, &wh, &cmt, &wh100, &whI, &wh0);
	dose = ind->dose[l];
	shiftF=0;
	oral=oral0;
	if (cmt != linCmt){
	  if (!oral){
	    continue;
	  } else if (cmt != linCmt+1) continue;
	  oral = 0; // Switch to "central" compartment
	  shiftF=1;
	}
	switch(whI){
	case 7:
	  continue;
	case 6:
	  continue;
	case 8: // Duration is modeled
	case 9: // Rate is modeled
	  tT = t - ind->all_times[ind->idose[l]] ;
	  thisT = tT - TLAG;
	  tau = ind->ii[l];
	  if (whI == 9){
	    tinf  = dose/d_rate;
	    rate  = d_rate;
	  } else {
	    tinf  = d_dur;
	    rate  = dose/d_dur;
	  }
	  dose=NA_REAL;
	case 2:
	case 1:
	  if (oral) {
	    REprintf(_("infusions to depot are not possible with the linear solved system\n"));
	    ret = NA_REAL;
	    return g;
	  }
	  if (wh0 == 30){
	    REprintf(_("you cannot turn off a compartment with a solved system\n"));
	    ret = NA_REAL;
	    return g;
	  }
	  // Steady state
	  if (ISNA(dose)){
	  } else if (dose > 0){
	    // During infusion
	    tT = t - ind->all_times[ind->idose[l]] ;
	    thisT = tT - TLAG;
	    p = l+1;
	    while (p < ndoses && ind->dose[p] != -dose){
	      p++;
	    }
	    if (ind->dose[p] != -dose){
	      REprintf(_("could not find an end to the infusion\n"));
	      ret = NA_REAL;
	      return g;
	    }
	    tinf  = ind->all_times[ind->idose[p]] - ind->all_times[ind->idose[l]];
	    tau = ind->ii[l];
	    rate  = dose;
	    if (tT >= tinf) continue;
	  } else {
	    // After  infusion
	    p = l-1;
	    while (p > 0 && ind->dose[p] != -dose){
	      p--;
	    }
	    if (ind->dose[p] != -dose){
	      REprintf(_("Could not find a start to the infusion\n"));
	      ret = NA_REAL;
	      return g;
	    }
	    tinf  = ind->all_times[ind->idose[l]] - ind->all_times[ind->idose[p]];
	    tau = ind->ii[p];
	    tT = t - ind->all_times[ind->idose[p]];
	    thisT = tT -TLAG;
	    rate  = -dose;
	  }
	  if (thisT < 0) continue;
	  if (F <= 0) {
	    REprintf(_("bioavailability cannot be negative or zero\n"));
	    ret = NA_REAL;
	    return g;
	  }
	  if (whI == 1){ // Duration changes
	    tinf = tinf*F;
	  } else { // Rate Changes
	    rate *= F;
	  }
	  if (wh0 == 10 || wh0 == 20){
	    if (tinf >= tau){
	      REprintf("infusion time greater then inter-dose interval, 'ss' cannot be calculated\n");
	      ret = NA_REAL;
	      return g;
	    }
	    if (thisT < tinf){ // during infusion
	      for (int jj = ncmt; jj--;){
		ret += rate*A*(1/(ALPHA+DOUBLE_EPS))*((1-exp(-ALPHA*thisT))+
					 exp(-ALPHA*tau)*(1-exp(-ALPHA*tinf))*
					 exp(-ALPHA*(thisT-tinf))/(1-exp(-ALPHA*tau)+DOUBLE_EPS));
	      }
	      if (wh0 == 10){
		return g;
	      } 
	    } else { // after infusion
	      for (int jj = ncmt; jj--;){
		ret += rate*A*(1/(ALPHA+DOUBLE_EPS))*((1-exp(-ALPHA*tinf))*exp(-ALPHA*(thisT-tinf))/(1-exp(-ALPHA*tau)+DOUBLE_EPS));
	      }
	      if (wh0 == 10){
		return g; // Was a reset event.
	      } 
	    }
	  } else {
	    //during infusion
#define t1 ((thisT < tinf) ? thisT : tinf)
	    // after infusion
#define t2 ((thisT > tinf) ? thisT - tinf : 0.0)
	    for (int jj = ncmt; jj--;){
	      ret +=  rate*A*(1/(ALPHA+DOUBLE_EPS))*(1.0-exp(-ALPHA*t1))*exp(-ALPHA*t2);
	    }
	  }
#undef t1
#undef t2
	  break;
	case 0:
	  if (wh0 == 10 || wh0 == 20){
	    // steady state
	    tT = t - ind->all_times[ind->idose[l]];
	    thisT = tT -TLAG;
	    if (thisT < 0) continue;
	    tau = ind->ii[l];
	    res = ((oral == 1) ? exp(-ka*thisT)/(1-exp(-ka*tau)) : 0.0);
	    for (int jj = ncmt; jj--;){
	      ret += dose*F*A*(exp(-ALPHA*thisT)/(1-exp(-ALPHA*tau))-res);
	    }
	    // ss=1 is equivalent to a reset + ss dose
	    if (wh0 == 10){
	      return g;  // Was a reset event.
	    }
	  } else if (wh0 == 30) {
	    REprintf(_("you cannot turn off a compartment with a solved system\n"));
	    ret = NA_REAL;
	    return g;
	  } else {
	    tT = t - ind->all_times[ind->idose[l]];
	    thisT = tT -TLAG;
	    if (thisT < 0) continue;
	    res = ((oral == 1) ? exp(-ka*thisT) : 0.0);
	    for (int jj = ncmt; jj--;)
	      ret +=  dose*F*A*(exp(-ALPHA*thisT)-res);
	  }
	  break;
	default:
	  REprintf(_("invalid evid in linear solved system\n"));
	  ret = NA_REAL;
	  return g;
	}
      }
#undef ret
#undef A
#undef F
#undef TLAG
#undef ALPHA
      return g;
    }
    struct linCmtFun {
      const double t_;
      const int ncmt_, linCmt_, oral_, trans_;
      rx_solving_options_ind *ind_;
      linCmtFun(const double t,
		const int ncmt,
		const int oral,
		const int trans,
		const int linCmt,
		rx_solving_options_ind *ind) :
	t_(t),
	ncmt_(ncmt),
	linCmt_(linCmt),
	oral_(oral),
	trans_(trans),
	ind_(ind)
      { }
      template <typename T>
	Eigen::Matrix<T, Eigen::Dynamic, 1> operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& params) const {
	return genericCmtInterface(params, t_, oral_, trans_, ncmt_, linCmt_, ind_);
      }
    };
  } // math
}
#undef beta1
#undef gamma1
#undef alpha1
#undef alpha
#undef beta
#undef gamma
#undef ka
#undef tlag

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> MatrixPd;

extern "C" double linCmtB(rx_solve *rx, unsigned int id, double t, int linCmt,
			  int i_cmt, int trans, int val,
			  double dd_p1, double dd_v1,
			  double dd_p2, double dd_p3,
			  double dd_p4, double dd_p5,
			  double dd_ka,
			  double dd_tlag, double dd_tlag2,
			  double dd_F, double dd_F2,
			  double dd_rate, double dd_dur){
  rx_solving_options_ind *ind;
  ind = &(rx->subjects[id]);
  if (ind->linCmtT == t){
    return ind->linCmtB[val];
  }
  if (i_cmt <= 0 || i_cmt >= 4) throw std::invalid_argument("invalid i_cmt");
  unsigned int ncmt = 1;
  if (dd_p4 > 0.){
    ncmt = 3;
  } else if (dd_p2 > 0.){
    ncmt = 2;
  }
  MatrixPd params(2*ncmt+7,1);
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
  int oral0 = (dd_ka > 0);
  params(2*ncmt,0) = dd_ka;
  params(2*ncmt+1,0) = dd_tlag;
  params(2*ncmt+2,0) = dd_tlag2;
  params(2*ncmt+3,0) = dd_F;
  params(2*ncmt+4,0) = dd_F2;
  params(2*ncmt+5,0) = dd_rate;
  params(2*ncmt+6,0) = dd_dur;

  stan::math::linCmtFun f(t, ncmt, oral0, trans, linCmt, ind);
  Eigen::VectorXd fx;
  Eigen::Matrix<double, -1, -1> J;
  // ncol = npar
  /// nrow = nobs
  stan::math::jacobian(f, params, fx, J);
  ind->linCmtT=t;
  ind->linCmtB[0] = fx[0];
  ind->linCmtB[1] = J(0,0);
  ind->linCmtB[2] = J(0,1);
  if (ncmt >= 2){
    ind->linCmtB[3] = J(0,2);
    ind->linCmtB[4] = J(0,3);
    if (ncmt == 3){
      ind->linCmtB[5] = J(0,4);
      ind->linCmtB[6] = J(0,5);
    } else {
      ind->linCmtB[5] = 0;
      ind->linCmtB[6] = 0;
    }
  } else {
    ind->linCmtB[3] = 0;
    ind->linCmtB[4] = 0;
    ind->linCmtB[5] = 0;
    ind->linCmtB[6] = 0;
  }
  ind->linCmtB[7]  = J(0,2*ncmt);
  ind->linCmtB[8]  = J(0,2*ncmt+1);
  ind->linCmtB[9]  = J(0,2*ncmt+2);
  ind->linCmtB[10] = J(0,2*ncmt+3);
  ind->linCmtB[11] = J(0,2*ncmt+4);
  ind->linCmtB[12] = J(0,2*ncmt+5);
  ind->linCmtB[13] = J(0,2*ncmt+6);
  // Rprintf("id:%d\t%f\t%f\t%f\t%f\t%f\t%f\n", ind->id, t, ind->linCmtB[0], ind->linCmtB[1], ind->linCmtB[2], ind->linCmtB[3], ind->linCmtB[4]);
  // Rprintf("\t%f\t%f\t%f\t%f\t%f\n", params(0,0), params(1,0), params(2,0), params(3,0), params(2*ncmt,0));
  return ind->linCmtB[val];
}
