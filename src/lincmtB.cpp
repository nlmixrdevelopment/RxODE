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
#define k23   g(1, 0)
      
#define k21   g(1, 1)
#define k32   g(1, 1)
      
#define k13   g(2, 0)
#define k24   g(2, 0)
      
#define k31   g(2, 1)
#define k42   g(2, 1)
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

#define v     g(0, 0)
#define k20   g(0, 1)
#define kel   g(0, 1)
#define k23   g(1, 0)
#define k32   g(1, 1)
#define k24   g(2, 0)
#define k42   g(2, 1)

#define k10   g(0, 1)
#define k12   g(1, 0)
#define k21   g(1, 1)
#define k13   g(2, 0)
#define k31   g(2, 1)
    
#define A1    A(0, 0)
#define A2    A(1, 0)
#define A3    A(2, 0)
#define A4    A(3, 0)
#define r1    rate(0, 0)
#define Doserate    rate(0, 0)
#define r2    rate(1, 0)
#define b1    bolus(0, 0)
#define b2    bolus(1, 0)
#define A1last Alast(0, 0)
#define A2last Alast(1, 0)
#define A3last Alast(2, 0)
#define A4last Alast(3, 0)
    // one compartment ka translations ncmt=1
#define tlag  params(2, 0)
#define F     params(3, 0)
#define rate1 params(4, 0)
#define dur1  params(5, 0)
#define ka    params(6, 0)
#define tlag2 params(7, 0)
#define f2    params(8, 0)
#define dur2  params(9, 0)    

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    oneCmtKaRateSSr1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		     Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		     Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(2, 1);
      A1 = r1/ka;
      A2 = r1/k20;
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    oneCmtKaRateSSr2(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		     Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		     Eigen::Matrix<T, Eigen::Dynamic, 2>& rate){
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(2, 1);
      A1 = 0;
      A2 = r2/k20;
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    oneCmtKaRateSStr1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		      Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		      Eigen::Matrix<T, Eigen::Dynamic, 2>& rate,
		      T tinf, T tau){
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(2, 1);
      T eKa = exp(-ka*(tau-tinf))/(1.0-exp(-tau*ka));
      T eiKa = exp(-ka*tinf);
      T eiK = exp(-k20*tinf);
      T eK = exp(-k20*(tau-tinf))/(1.0-exp(-tau*k20));
      A1=eKa*(r1/ka - eiKa*r1/ka);
      A2=eK*(r1/k20 + eiKa*r1/(-k20 + ka) - eiK*r1*ka/(ka*k20 - k20*k20)) +
	ka*(eK - eKa)*(r1/ka - eiKa*r1/ka)/(-k20 + ka);
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    oneCmtKaRateSStr2(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		      Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		      Eigen::Matrix<T, Eigen::Dynamic, 2>& rate,
		      T tinf, T tau){
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(2, 1);
      T eiK = exp(-k20*tinf);
      T eK = exp(-k20*(tau-tinf))/(1.0-exp(-k20*tau));
      A1=0.0;
      A2=eK*(r2/k20 - eiK*r2*(-k20 + ka)/(ka*k20 - k20*k20));
      return A;
    }
    
    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    oneCmtKaRate(T t, Eigen::Matrix<T, Eigen::Dynamic, 2>& Alast,
		 Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(2, 1);
      T eKa = exp(-ka*t);
      T e20 = exp(-k20*t);
      A1 = r1/ka-((r1-A1last*ka)*eKa)/ka + b1;
      T A21 = ((r1-A1last*ka)*eKa)/(ka-k20);
      T A22=(((ka-k20)*r2+ka*r1+(-A2last-A1last)*k20*ka+A2last*k20*k20)*e20)/(k20*ka-k20*k20);
      T A23 =(r2+r1)/k20;
      A2 = A21 - A22 + A23 + b2;
      return A;
    }
    // Undefine extras
#undef tlag
#undef F
#undef rate1
#undef dur1
#undef ka
#undef tlag2
#undef f2
#undef dur2
    // two compartment ka translations ncmt=1
#define tlag  params(4,  0)
#define F     params(5,  0)
#define rate1 params(6,  0)
#define dur1  params(7,  0)
#define ka    params(8,  0)
#define tlag2 params(9,  0)
#define f2    params(10, 0)
#define dur2  params(11, 0)

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    twoCmtKaRateSSr1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		     Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		     Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(3, 1);
      T s = k23+k32+k20;
      //#Calculate roots
      T beta  = 0.5*(s - sqrt(s*s - 4*k32*k20));
      T alpha = k32*k20/beta;
      A1=r1/ka;
      A2=r1*k32/(beta*alpha);
      A3=r1*k23/(beta*alpha);
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    twoCmtKaRateSSr2(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		     Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		     Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(3, 1);
      T s = k23+k32+k20;
      //#Calculate roots
      T beta  = 0.5*(s - sqrt(s*s - 4*k32*k20));
      T alpha = k32*k20/beta;
      A1=0;
      A2=r2*k32/(beta*alpha);
      A3=r2*k23/(beta*alpha);
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    twoCmtKaRateSStr1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		      Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		      Eigen::Matrix<T, Eigen::Dynamic, 2>& rate,
		      T tinf, T tau){
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(3, 1);
      T s = k23+k32+k20;
      //#Calculate roots
      T beta  = 0.5*(s - sqrt(s*s - 4*k32*k20));
      T alpha = k32*k20/beta;

      T eA = exp(-alpha*(tau-tinf))/(1.0-exp(-alpha*tau));
      T eB = exp(-beta*(tau-tinf))/(1.0-exp(-beta*tau));

      T eiA = exp(-alpha*tinf);
      T eiB = exp(-beta*tinf);

      T alpha2 = alpha*alpha;
      T alpha3 = alpha2*alpha;

      T beta2 = beta*beta;
      T beta3 = beta2*beta;

      T ka2 = ka*ka;

      T eKa = exp(-ka*(tau-tinf))/(1.0-exp(-ka*tau));
      T eiKa = exp(-ka*tinf);

      A1=eKa*(r1/ka - eiKa*r1/ka);
      A2=(eA*(-alpha*(r1*k32/(beta*alpha) + eiKa*r1*(-k32 + ka)/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*(-alpha + k32)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*(-beta + k32)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + k32*(r1*k23/(beta*alpha) - eiKa*r1*k23/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*k23/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*k23/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + k32*(r1*k32/(beta*alpha) + eiKa*r1*(-k32 + ka)/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*(-alpha + k32)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*(-beta + k32)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3))) - eB*(-beta*(r1*k32/(beta*alpha) + eiKa*r1*(-k32 + ka)/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*(-alpha + k32)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*(-beta + k32)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + k32*(r1*k23/(beta*alpha) - eiKa*r1*k23/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*k23/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*k23/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + k32*(r1*k32/(beta*alpha) + eiKa*r1*(-k32 + ka)/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*(-alpha + k32)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*(-beta + k32)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3))))/(-alpha + beta) + ka*(eA*(-alpha + k32)/((-alpha + beta)*(-alpha + ka)) + eB*(-beta + k32)/((-beta + ka)*(alpha - beta)) + eKa*(k32 - ka)/((beta - ka)*(alpha - ka)))*(r1/ka - eiKa*r1/ka);
      A3=(eA*(-alpha*(r1*k23/(beta*alpha) - eiKa*r1*k23/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*k23/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*k23/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + k23*(r1*k32/(beta*alpha) + eiKa*r1*(-k32 + ka)/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*(-alpha + k32)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*(-beta + k32)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + (k20 + k23)*(r1*k23/(beta*alpha) - eiKa*r1*k23/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*k23/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*k23/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3))) - eB*(-beta*(r1*k23/(beta*alpha) - eiKa*r1*k23/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*k23/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*k23/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + k23*(r1*k32/(beta*alpha) + eiKa*r1*(-k32 + ka)/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*(-alpha + k32)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*(-beta + k32)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + (k20 + k23)*(r1*k23/(beta*alpha) - eiKa*r1*k23/(beta*alpha + ka*(-alpha - beta) + ka2) - eiA*r1*ka*k23/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r1*ka*k23/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3))))/(-alpha + beta) + ka*k23*(eA/((-alpha + beta)*(-alpha + ka)) + eB/((-beta + ka)*(alpha - beta)) + eKa/((beta - ka)*(alpha - ka)))*(r1/ka - eiKa*r1/ka);
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    twoCmtKaRateSStr2(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		      Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		      Eigen::Matrix<T, Eigen::Dynamic, 2>& rate,
		      T tinf, T tau) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(3, 1);
      T E2 = k20+k23;
      T E3 = k32;
      T s = k23+k32+k20;
      //#Calculate roots
      T beta  = 0.5*(s - sqrt(s*s - 4*k32*k20));
      T alpha = k32*k20/beta;

      T eA = exp(-alpha*(tau-tinf))/(1.0-exp(-alpha*tau));
      T eB = exp(-beta*(tau-tinf))/(1.0-exp(-beta*tau));

      T eiA = exp(-alpha*tinf);
      T eiB = exp(-beta*tinf);

      T alpha2 = alpha*alpha;
      T alpha3 = alpha2*alpha;

      T beta2 = beta*beta;
      T beta3 = beta2*beta;
      A1=0.0;
      A2=(eA*(E3*(r2*k32/(beta*alpha) - eiA*r2*(-k32*alpha + ka*(-alpha + k32) + alpha2)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k32*beta + ka*(-beta + k32) + beta2)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) - alpha*(r2*k32/(beta*alpha) - eiA*r2*(-k32*alpha + ka*(-alpha + k32) + alpha2)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k32*beta + ka*(-beta + k32) + beta2)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + k32*(r2*k23/(beta*alpha) - eiA*r2*(-k23*alpha + ka*k23)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k23*beta + ka*k23)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3))) - eB*(E3*(r2*k32/(beta*alpha) - eiA*r2*(-k32*alpha + ka*(-alpha + k32) + alpha2)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k32*beta + ka*(-beta + k32) + beta2)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) - beta*(r2*k32/(beta*alpha) - eiA*r2*(-k32*alpha + ka*(-alpha + k32) + alpha2)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k32*beta + ka*(-beta + k32) + beta2)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + k32*(r2*k23/(beta*alpha) - eiA*r2*(-k23*alpha + ka*k23)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k23*beta + ka*k23)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3))))/(-alpha + beta);
      A3=(eA*(E2*(r2*k23/(beta*alpha) - eiA*r2*(-k23*alpha + ka*k23)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k23*beta + ka*k23)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) - alpha*(r2*k23/(beta*alpha) - eiA*r2*(-k23*alpha + ka*k23)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k23*beta + ka*k23)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + k23*(r2*k32/(beta*alpha) - eiA*r2*(-k32*alpha + ka*(-alpha + k32) + alpha2)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k32*beta + ka*(-beta + k32) + beta2)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3))) - eB*(E2*(r2*k23/(beta*alpha) - eiA*r2*(-k23*alpha + ka*k23)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k23*beta + ka*k23)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) - beta*(r2*k23/(beta*alpha) - eiA*r2*(-k23*alpha + ka*k23)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k23*beta + ka*k23)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3)) + k23*(r2*k32/(beta*alpha) - eiA*r2*(-k32*alpha + ka*(-alpha + k32) + alpha2)/(-beta*alpha2 + ka*(beta*alpha - alpha2) + alpha3) + eiB*r2*(-k32*beta + ka*(-beta + k32) + beta2)/(beta2*alpha + ka*(-beta*alpha + beta2) - beta3))))/(-alpha + beta);
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    twoCmtKaRate(T t, Eigen::Matrix<T, Eigen::Dynamic, 2>& Alast,
		 Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(3, 1);
      T E2 =  k20+ k23;
      T s = k23+k32+k20;
      //#Calculate roots
      T beta  = 0.5*(s - sqrt(s*s - 4*k32*k20));
      T alpha = k32*k20/beta;

      T eKa = exp(-ka*t);
      T eA = exp(-alpha*t);
      T eB = exp(-beta*t);

      T ka2 = ka*ka;
  
      T alpha2 = alpha*alpha;
      T alpha3 = alpha2*alpha;
  
      T beta2 = beta*beta;
      T beta3 = beta2*beta;
  
      A1 = b1+r1/ka-((r1-A1last*ka)*eKa)/ka;
      A2 = b2+(((ka-k32)*r1-A1last*ka2+A1last*k32*ka)*eKa)/(ka2+(-beta-alpha)*ka+alpha*beta)+((((k32-beta)*ka-beta*k32+beta2)*r2+(k32-beta)*ka*r1+((-A3last-A2last-A1last)*beta*k32+(A2last+A1last)*beta2)*ka+(A3last+A2last)*beta2*k32-A2last*beta3)*eB)/((beta2-alpha*beta)*ka-beta3+alpha*beta2)-((((k32-alpha)*ka-alpha*k32+alpha2)*r2+(k32-alpha)*ka*r1+((-A3last-A2last-A1last)*alpha*k32+(A2last+A1last)*alpha2)*ka+(A3last+A2last)*alpha2*k32-A2last*alpha3)*eA)/((alpha*beta-alpha2)*ka-alpha2*beta+alpha3)+(k32*r2+k32*r1)/(alpha*beta);
      A3 = -((k23*r1-A1last*k23*ka)*eKa)/(ka2+(-beta-alpha)*ka+alpha*beta)+(((k23*ka-beta*k23)*r2+k23*ka*r1+((-A2last-A1last)*beta*k23+A3last*beta2-A3last*E2*beta)*ka+A2last*beta2*k23-A3last*beta3+A3last*E2*beta2)*eB)/((beta2-alpha*beta)*ka-beta3+alpha*beta2)-(((k23*ka-alpha*k23)*r2+k23*ka*r1+((-A2last-A1last)*alpha*k23+A3last*alpha2-A3last*E2*alpha)*ka+A2last*alpha2*k23-A3last*alpha3+A3last*E2*alpha2)*eA)/((alpha*beta-alpha2)*ka-alpha2*beta+alpha3)+(k23*r2+k23*r1)/(alpha*beta);
      return A;
    }
    
        // undefine extras
#undef tlag
#undef F
#undef rate1
#undef dur1
#undef ka
#undef tlag2
#undef f2
#undef dur2
        // three compartment ka translations ncmt=1
#define tlag  params(6,  0)
#define F     params(7,  0)
#define rate1 params(8,  0)
#define dur1  params(9,  0)
#define ka    params(10, 0)
#define tlag2 params(11, 0)
#define f2    params(12, 0)
#define dur2  params(13, 0)

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    threeCmtKaRateSSr1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		       Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		       Eigen::Matrix<T, Eigen::Dynamic, 2>& rate){
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(4, 1);
      //##Calculate roots - see Upton, 2004
      T j = k23+k20+k32+k42+k24;
      T k = k23*k42+k20*k32+k20*k42+k32*k42+k24*k32;
      T l = k20*k32*k42;

      T m = 0.3333333333333333*(3.0*k - j*j);
      T n = 0.03703703703703703*(2.0*j*j*j - 9.0*j*k + 27.0*l);
      T Q = 0.25*n*n + 0.03703703703703703*m*m*m;

      T alpha = sqrt(-Q);
      T beta = -0.5*n;
      T rho=sqrt(beta*beta+alpha*alpha);
      T theta = atan2(alpha,beta);
      T ct3 = cos(0.3333333333333333*theta);
      T rho3 = pow(rho,0.3333333333333333);
      T st3 = 1.732050807568877193177*sin(0.3333333333333333*theta);
      T j3 = 0.3333333333333333*j;
      T lam1 = j3  + rho3*(ct3 + st3);
      T lam2 = j3 + rho3*(ct3 - st3);
      T lam3 = j3 -(2.0*rho3*ct3);
      T l123 = 1.0/(lam1*lam2*lam3);
      A1=r1/ka;
      A2=r1*k42*k32*l123;
      A3=r1*k42*k23*l123;
      A4=r1*k24*k32*l123;
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    threeCmtKaRateSSr2(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		       Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		       Eigen::Matrix<T, Eigen::Dynamic, 2>& rate){
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(4, 1);
      T j = k23+k20+k32+k42+k24;
      T k = k23*k42+k20*k32+k20*k42+k32*k42+k24*k32;
      T l = k20*k32*k42;

      T m = 0.3333333333333333*(3.0*k - j*j);
      T n = 0.03703703703703703*(2.0*j*j*j - 9.0*j*k + 27.0*l);
      T Q = 0.25*n*n + 0.03703703703703703*m*m*m;

      T alpha = sqrt(-Q);
      T beta = -0.5*n;
      T rho=sqrt(beta*beta+alpha*alpha);
      T theta = atan2(alpha,beta);
      T ct3 = cos(0.3333333333333333*theta);
      T rho3 = pow(rho,0.3333333333333333);
      T st3 = 1.732050807568877193177*sin(0.3333333333333333*theta);
      T j3 = 0.3333333333333333*j;
      T lam1 = j3  + rho3*(ct3 + st3);
      T lam2 = j3 + rho3*(ct3 - st3);
      T lam3 = j3 -(2.0*rho3*ct3);
      T l123 = 1.0/(lam1*lam2*lam3);
      A1=0;
      A2=r2*k42*k32*l123;
      A3=r2*k42*k23*l123;
      A4=r2*k24*k32*l123;
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    threeCmtKaRateSStr1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
			Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
			Eigen::Matrix<T, Eigen::Dynamic, 2>& rate,
			T tinf, T tau) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(4, 1);
      T E2 =  k20+ k23 + k24;
      T E3 = k32;
      T E4 = k42;
      //##Calculate roots - see Upton, 2004
      T j = k23+k20+k32+k42+k24;
      T k = k23*k42+k20*k32+k20*k42+k32*k42+k24*k32;
      T l = k20*k32*k42;

      T m = 0.3333333333333333*(3.0*k - j*j);
      T n = 0.03703703703703703*(2.0*j*j*j - 9.0*j*k + 27.0*l);
      T Q = 0.25*n*n + 0.03703703703703703*m*m*m;

      T alpha = sqrt(-Q);
      T beta = -0.5*n;
      T rho=sqrt(beta*beta+alpha*alpha);
      T theta = atan2(alpha,beta);
      T ct3 = cos(0.3333333333333333*theta);
      T rho3 = pow(rho,0.3333333333333333);
      T st3 = 1.732050807568877193177*sin(0.3333333333333333*theta);
      T j3 = 0.3333333333333333*j;
      T lam1 = j3  + rho3*(ct3 + st3);
      T lam2 = j3 + rho3*(ct3 - st3);
      T lam3 = j3 -(2.0*rho3*ct3);
  
      T eKa = exp(-ka*(tau-tinf))/(1.0-exp(-ka*tau));
      T eiKa = exp(-ka*tinf);

      T eL1 = exp(-lam1*(tau-tinf))/(1.0-exp(-lam1*tau));
      T eiL1 = exp(-lam1*tinf);

      T eL2 = exp(-lam2*(tau-tinf))/(1.0-exp(-lam2*tau));
      T eiL2 = exp(-lam2*tinf);

      T eL3 = exp(-lam3*(tau-tinf))/(1.0-exp(-lam3*tau));
      T eiL3 = exp(-lam3*tinf);

      T ka2 = ka*ka;
      T ka3 = ka2*ka;

      T lam12 = lam1*lam1;
      T lam13 = lam12*lam1;
      T lam14 = lam13*lam1;

      T lam22 = lam2*lam2;
      T lam23 = lam22*lam2;
      T lam24 = lam23*lam2;

      T lam32 = lam3*lam3;
      T lam33 = lam32*lam3;
      T lam34 = lam33*lam3;
      A1=eKa*(r1/ka - eiKa*r1/ka);
      A2=(eL1*(E4 - lam1)*(E3 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)) + eL2*(E4 - lam2)*(E3 - lam2)/((-lam2 + lam3)*(lam1 - lam2)) + eL3*(E3 - lam3)*(E4 - lam3)/((lam2 - lam3)*(lam1 - lam3)))*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) + ka*(eL1*(E4 - lam1)*(E3 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)*(ka - lam1)) + eL2*(E4 - lam2)*(E3 - lam2)/((-lam2 + lam3)*(lam1 - lam2)*(ka - lam2)) + eL3*(E3 - lam3)*(E4 - lam3)/((lam2 - lam3)*(lam1 - lam3)*(ka - lam3)) + eKa*(E3 - ka)*(E4 - ka)/((-ka + lam1)*(-ka + lam2)*(-ka + lam3)))*(r1/ka - eiKa*r1/ka) + eL1*(-lam1*(k32*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3)) + k42*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3))) + E3*k42*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3)) + E4*k32*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3)))/((lam1 - lam3)*(lam1 - lam2)) + eL3*(lam3*(k32*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3)) + k42*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3))) - (E3*k42*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3)) + E4*k32*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3))))/((-lam2 + lam3)*(lam1 - lam3)) + eL2*(lam2*(k32*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3)) + k42*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3))) - (E3*k42*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3)) + E4*k32*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3))))/((lam2 - lam3)*(lam1 - lam2));
      A3=(eL1*(E4 - lam1)*(E2 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)) + eL2*(E4 - lam2)*(E2 - lam2)/((-lam2 + lam3)*(lam1 - lam2)) + eL3*(E2 - lam3)*(E4 - lam3)/((lam2 - lam3)*(lam1 - lam3)))*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3)) + ka*k23*(eL1*(E4 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)*(ka - lam1)) + eL2*(E4 - lam2)/((-lam2 + lam3)*(lam1 - lam2)*(ka - lam2)) + eL3*(E4 - lam3)/((lam2 - lam3)*(lam1 - lam3)*(ka - lam3)) + eKa*(E4 - ka)/((-ka + lam1)*(-ka + lam2)*(-ka + lam3)))*(r1/ka - eiKa*r1/ka) + eL1*(E4*k23*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) - k23*lam1*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) + k42*k23*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3)) - k42*k24*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3)))/((lam1 - lam3)*(lam1 - lam2)) + eL3*(k23*lam3*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) - (E4*k23*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) + k42*k23*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3)) - k42*k24*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3))))/((-lam2 + lam3)*(lam1 - lam3)) + eL2*(k23*lam2*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) - (E4*k23*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) + k42*k23*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3)) - k42*k24*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3))))/((lam2 - lam3)*(lam1 - lam2));
      A4=(eL1*(E3 - lam1)*(E2 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)) + eL2*(E2 - lam2)*(E3 - lam2)/((-lam2 + lam3)*(lam1 - lam2)) + eL3*(E2 - lam3)*(E3 - lam3)/((lam2 - lam3)*(lam1 - lam3)))*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3)) + ka*k24*(eL1*(E3 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)*(ka - lam1)) + eL2*(E3 - lam2)/((-lam2 + lam3)*(lam1 - lam2)*(ka - lam2)) + eL3*(E3 - lam3)/((lam2 - lam3)*(lam1 - lam3)*(ka - lam3)) + eKa*(E3 - ka)/((-ka + lam1)*(-ka + lam2)*(-ka + lam3)))*(r1/ka - eiKa*r1/ka) + eL1*(E3*k24*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) - k23*k32*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3)) + k24*k32*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3)) - k24*lam1*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)))/((lam1 - lam3)*(lam1 - lam2)) + eL3*(k24*lam3*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) - (E3*k24*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) - k23*k32*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3)) + k24*k32*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3))))/((-lam2 + lam3)*(lam1 - lam3)) + eL2*(k24*lam2*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) - (E3*k24*(-eiKa*r1*(k42*k32 + (-k32 - k42)*ka + ka2)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) - eiL1*r1*(-ka*lam12 - ka*k42*k32 + (k32 + k42)*ka*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r1*(-ka*lam22 - ka*k42*k32 + (k32 + k42)*ka*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r1*(-ka*lam32 - ka*k42*k32 + (k32 + k42)*ka*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k32/(lam2*lam1*lam3)) - k23*k32*(eiKa*r1*(-k24*k32 + ka*k24)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(ka*k24*k32 - ka*k24*lam1)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(ka*k24*k32 - ka*k24*lam2)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(ka*k24*k32 - ka*k24*lam3)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k24*k32/(lam2*lam1*lam3)) + k24*k32*(eiKa*r1*(-k42*k23 + ka*k23)/(ka2*lam1 + lam2*(-ka*lam1 + ka2) + lam3*(-ka*lam1 + lam2*(-ka + lam1) + ka2) - ka3) + eiL1*r1*(-ka*k23*lam1 + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r1*(-ka*k23*lam2 + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r1*(-ka*k23*lam3 + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r1*k42*k23/(lam2*lam1*lam3))))/((lam2 - lam3)*(lam1 - lam2));
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    threeCmtKaRateSStr2(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
			Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
			Eigen::Matrix<T, Eigen::Dynamic, 2>& rate,
			T tinf, T tau) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(4, 1);
      T E2 =  k20+ k23 + k24;
      T E3 = k32;
      T E4 = k42;
      //##Calculate roots - see Upton, 2004
      T j = k23+k20+k32+k42+k24;
      T k = k23*k42+k20*k32+k20*k42+k32*k42+k24*k32;
      T l = k20*k32*k42;

      T m = 0.3333333333333333*(3.0*k - j*j);
      T n = 0.03703703703703703*(2.0*j*j*j - 9.0*j*k + 27.0*l);
      T Q = 0.25*n*n + 0.03703703703703703*m*m*m;

      T alpha = sqrt(-Q);
      T beta = -0.5*n;
      T rho=sqrt(beta*beta+alpha*alpha);
      T theta = atan2(alpha,beta);
      T ct3 = cos(0.3333333333333333*theta);
      T rho3 = pow(rho,0.3333333333333333);
      T st3 = 1.732050807568877193177*sin(0.3333333333333333*theta);
      T j3 = 0.3333333333333333*j;
      T lam1 = j3  + rho3*(ct3 + st3);
      T lam2 = j3 + rho3*(ct3 - st3);
      T lam3 = j3 -(2.0*rho3*ct3);
  
      /* T eKa = 1.0/(1.0-exp(-ka*tau)); */
      /* T eiKa = exp(-ka*tinf); */

      T eL1 = exp(-lam1*(tau-tinf))/(1.0-exp(-lam1*tau));
      T eiL1 = exp(-lam1*tinf);

      T eL2 = exp(-lam2*(tau-tinf))/(1.0-exp(-lam2*tau));
      T eiL2 = exp(-lam2*tinf);

      T eL3 = exp(-lam3*(tau-tinf))/(1.0-exp(-lam3*tau));
      T eiL3 = exp(-lam3*tinf);

      /* T ka2 = ka*ka; */
      /* T ka3 = ka2*ka; */

      T lam12 = lam1*lam1;
      T lam13 = lam12*lam1;
      T lam14 = lam13*lam1;

      T lam22 = lam2*lam2;
      T lam23 = lam22*lam2;
      T lam24 = lam23*lam2;

      T lam32 = lam3*lam3;
      T lam33 = lam32*lam3;
      T lam34 = lam33*lam3;
      A1=0.0;
      A2=(eL1*(E4 - lam1)*(E3 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)) + eL2*(E4 - lam2)*(E3 - lam2)/((-lam2 + lam3)*(lam1 - lam2)) + eL3*(E3 - lam3)*(E4 - lam3)/((lam2 - lam3)*(lam1 - lam3)))*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) + eL1*(-lam1*(k32*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3)) + k42*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3))) + E3*k42*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3)) + E4*k32*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3)))/((lam1 - lam3)*(lam1 - lam2)) + eL3*(lam3*(k32*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3)) + k42*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3))) - (E3*k42*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3)) + E4*k32*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3))))/((-lam2 + lam3)*(lam1 - lam3)) + eL2*(lam2*(k32*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3)) + k42*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3))) - (E3*k42*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3)) + E4*k32*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3))))/((lam2 - lam3)*(lam1 - lam2));
      A3=(eL1*(E4 - lam1)*(E2 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)) + eL2*(E4 - lam2)*(E2 - lam2)/((-lam2 + lam3)*(lam1 - lam2)) + eL3*(E2 - lam3)*(E4 - lam3)/((lam2 - lam3)*(lam1 - lam3)))*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3)) + eL1*(E4*k23*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) - k23*lam1*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) + k42*k23*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3)) - k42*k24*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3)))/((lam1 - lam3)*(lam1 - lam2)) + eL3*(k23*lam3*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) - (E4*k23*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) + k42*k23*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3)) - k42*k24*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3))))/((-lam2 + lam3)*(lam1 - lam3)) + eL2*(k23*lam2*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) - (E4*k23*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) + k42*k23*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3)) - k42*k24*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3))))/((lam2 - lam3)*(lam1 - lam2));
      A4=(eL1*(E3 - lam1)*(E2 - lam1)/((-lam1 + lam3)*(-lam1 + lam2)) + eL2*(E2 - lam2)*(E3 - lam2)/((-lam2 + lam3)*(lam1 - lam2)) + eL3*(E2 - lam3)*(E3 - lam3)/((lam2 - lam3)*(lam1 - lam3)))*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3)) + eL1*(E3*k24*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) - k23*k32*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3)) + k24*k32*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3)) - k24*lam1*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)))/((lam1 - lam3)*(lam1 - lam2)) + eL3*(k24*lam3*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) - (E3*k24*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) - k23*k32*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3)) + k24*k32*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3))))/((-lam2 + lam3)*(lam1 - lam3)) + eL2*(k24*lam2*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) - (E3*k24*(-eiL1*r2*(lam1*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam12 - ka*k42*k32 + lam13)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) + eiL2*r2*(lam2*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam22 - ka*k42*k32 + lam23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) - eiL3*r2*(lam3*(k42*k32 + (k32 + k42)*ka) + (-k32 - k42 - ka)*lam32 - ka*k42*k32 + lam33)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k32/(lam2*lam1*lam3)) - k23*k32*(eiL1*r2*(k24*lam12 + lam1*(-k24*k32 - ka*k24) + ka*k24*k32)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k24*lam22 + lam2*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k24*lam32 + lam3*(-k24*k32 - ka*k24) + ka*k24*k32)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k24*k32/(lam2*lam1*lam3)) + k24*k32*(eiL1*r2*(k23*lam12 + lam1*(-k42*k23 - ka*k23) + ka*k42*k23)/(-ka*lam13 + lam2*(ka*lam12 - lam13) + lam3*(ka*lam12 + lam2*(-ka*lam1 + lam12) - lam13) + lam14) - eiL2*r2*(k23*lam22 + lam2*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam23*(ka + lam1) + lam3*(lam22*(-ka - lam1) + ka*lam2*lam1 + lam23) - ka*lam22*lam1 - lam24) + eiL3*r2*(k23*lam32 + lam3*(-k42*k23 - ka*k23) + ka*k42*k23)/(lam32*(ka*lam1 + lam2*(ka + lam1)) + (-ka - lam1 - lam2)*lam33 - ka*lam2*lam1*lam3 + lam34) + r2*k42*k23/(lam2*lam1*lam3))))/((lam2 - lam3)*(lam1 - lam2));
      return A;
    }
    
    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    threeCmtKaRate(T t, Eigen::Matrix<T, Eigen::Dynamic, 2>& Alast,
		   Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(4, 1);
      T E2 =  k20+ k23 + k24;
      //##Calculate roots - see Upton, 2004
      T j = k23+k20+k32+k42+k24;
      T k = k23*k42+k20*k32+k20*k42+k32*k42+k24*k32;
      T l = k20*k32*k42;

      T m = 0.3333333333333333*(3.0*k - j*j);
      T n = 0.03703703703703703*(2.0*j*j*j - 9.0*j*k + 27.0*l);
      T Q = 0.25*n*n + 0.03703703703703703*m*m*m;

      T alpha = sqrt(-Q);
      T beta = -0.5*n;
      T rho=sqrt(beta*beta+alpha*alpha);
      T theta = atan2(alpha,beta);
      T ct3 = cos(0.3333333333333333*theta);
      T rho3 = pow(rho,0.3333333333333333);
      T st3 = 1.732050807568877193177*sin(0.3333333333333333*theta);
      T j3 = 0.3333333333333333*j;
      T lam1 = j3  + rho3*(ct3 + st3);
      T lam2 = j3 + rho3*(ct3 - st3);
      T lam3 = j3 -(2.0*rho3*ct3);
      T eKa = exp(-ka*t);
      A1 = b1+ r1/ka-((r1-A1last*ka)*eKa)/ka;
  
      T lam12 = lam1*lam1;
      T lam13 = lam12*lam1;
      T lam14 = lam13*lam1;
  
      T lam22 = lam2*lam2;
      T lam23 = lam22*lam2;
      T lam24 = lam23*lam2;
  
      T lam32 = lam3*lam3;
      T lam33 = lam32*lam3;
      T lam34 = lam33*lam3;

      T ka2 = ka*ka;
      T ka3 = ka2*ka;
  
      T a21 = (((lam33+(-ka-k42-k32)*lam32+((k42+k32)*ka+k32*k42)*lam3-k32*k42*ka)*r2+(-ka*lam32+(k42+k32)*ka*lam3-k32*k42*ka)*r1-A2last*lam34+((A2last+A1last)*ka+(A4last+A2last)*k42+(A3last+A2last)*k32)*lam33+(((-A4last-A2last-A1last)*k42+(-A3last-A2last-A1last)*k32)*ka+(-A4last-A3last-A2last)*k32*k42)*lam32+(A4last+A3last+A2last+A1last)*k32*k42*ka*lam3)*exp(-lam3*t))/(lam34+(-lam2-lam1-ka)*lam33+((lam1+ka)*lam2+ka*lam1)*lam32-ka*lam1*lam2*lam3);
      T a22 = (((lam23+(-ka-k42-k32)*lam22+((k42+k32)*ka+k32*k42)*lam2-k32*k42*ka)*r2+(-ka*lam22+(k42+k32)*ka*lam2-k32*k42*ka)*r1-A2last*lam24+((A2last+A1last)*ka+(A4last+A2last)*k42+(A3last+A2last)*k32)*lam23+(((-A4last-A2last-A1last)*k42+(-A3last-A2last-A1last)*k32)*ka+(-A4last-A3last-A2last)*k32*k42)*lam22+(A4last+A3last+A2last+A1last)*k32*k42*ka*lam2)*exp(-lam2*t))/((lam23+(-lam1-ka)*lam22+ka*lam1*lam2)*lam3-lam24+(lam1+ka)*lam23-ka*lam1*lam22);
      T a23 = (((lam13+(-ka-k42-k32)*lam12+((k42+k32)*ka+k32*k42)*lam1-k32*k42*ka)*r2+(-ka*lam12+(k42+k32)*ka*lam1-k32*k42*ka)*r1-A2last*lam14+((A2last+A1last)*ka+(A4last+A2last)*k42+(A3last+A2last)*k32)*lam13+(((-A4last-A2last-A1last)*k42+(-A3last-A2last-A1last)*k32)*ka+(-A4last-A3last-A2last)*k32*k42)*lam12+(A4last+A3last+A2last+A1last)*k32*k42*ka*lam1)*exp(-lam1*t))/(((lam12-ka*lam1)*lam2-lam13+ka*lam12)*lam3+(ka*lam12-lam13)*lam2+lam14-ka*lam13);
      T a24 = (((ka2+(-k42-k32)*ka+k32*k42)*r1-A1last*ka3+(A1last*k42+A1last*k32)*ka2-A1last*k32*k42*ka)*exp(-ka*t))/(((lam1-ka)*lam2-ka*lam1+ka2)*lam3+(ka2-ka*lam1)*lam2+ka2*lam1-ka3);
      T a25 = (k32*k42*r2+k32*k42*r1)/(lam1*lam2*lam3);
      A2 = b2-a21+a22-a23-a24+a25;
      T a31 = (((k23*lam32+(-k23*ka-k23*k42)*lam3+k23*k42*ka)*r2+(k23*k42*ka-k23*ka*lam3)*r1+A3last*lam34+(-A3last*ka-A3last*k42-A2last*k23-A3last*E2)*lam33+((A3last*k42+(A2last+A1last)*k23+A3last*E2)*ka+(-A3last*k24+(A4last+A2last)*k23+A3last*E2)*k42)*lam32+(A3last*k24+(-A4last-A2last-A1last)*k23-A3last*E2)*k42*ka*lam3)*exp(-lam3*t))/(lam34+(-lam2-lam1-ka)*lam33+((lam1+ka)*lam2+ka*lam1)*lam32-ka*lam1*lam2*lam3);
      T a32 = (((k23*lam22+(-k23*ka-k23*k42)*lam2+k23*k42*ka)*r2+(k23*k42*ka-k23*ka*lam2)*r1+A3last*lam24+(-A3last*ka-A3last*k42-A2last*k23-A3last*E2)*lam23+((A3last*k42+(A2last+A1last)*k23+A3last*E2)*ka+(-A3last*k24+(A4last+A2last)*k23+A3last*E2)*k42)*lam22+(A3last*k24+(-A4last-A2last-A1last)*k23-A3last*E2)*k42*ka*lam2)*exp(-lam2*t))/((lam23+(-lam1-ka)*lam22+ka*lam1*lam2)*lam3-lam24+(lam1+ka)*lam23-ka*lam1*lam22);
      T a33 = (((k23*lam12+(-k23*ka-k23*k42)*lam1+k23*k42*ka)*r2+(k23*k42*ka-k23*ka*lam1)*r1+A3last*lam14+(-A3last*ka-A3last*k42-A2last*k23-A3last*E2)*lam13+((A3last*k42+(A2last+A1last)*k23+A3last*E2)*ka+(-A3last*k24+(A4last+A2last)*k23+A3last*E2)*k42)*lam12+(A3last*k24+(-A4last-A2last-A1last)*k23-A3last*E2)*k42*ka*lam1)*exp(-lam1*t))/(((lam12-ka*lam1)*lam2-lam13+ka*lam12)*lam3+(ka*lam12-lam13)*lam2+lam14-ka*lam13);
      T a34 = (((k23*ka-k23*k42)*r1-A1last*k23*ka2+A1last*k23*k42*ka)*exp(-ka*t))/(((lam1-ka)*lam2-ka*lam1+ka2)*lam3+(ka2-ka*lam1)*lam2+ka2*lam1-ka3);
      T a35 = (k23*k42*r2+k23*k42*r1)/(lam1*lam2*lam3);
      A3=a31-a32+a33+a34+a35;
      T a41 = (((k24*lam32+(-k24*ka-k24*k32)*lam3+k24*k32*ka)*r2+(k24*k32*ka-k24*ka*lam3)*r1+A4last*lam34+(-A4last*ka-A4last*k32-A2last*k24-A4last*E2)*lam33+((A4last*k32+(A2last+A1last)*k24+A4last*E2)*ka+((A3last+A2last)*k24-A4last*k23+A4last*E2)*k32)*lam32+((-A3last-A2last-A1last)*k24+A4last*k23-A4last*E2)*k32*ka*lam3)*exp(-lam3*t))/(lam34+(-lam2-lam1-ka)*lam33+((lam1+ka)*lam2+ka*lam1)*lam32-ka*lam1*lam2*lam3);
      T a42 = (((k24*lam22+(-k24*ka-k24*k32)*lam2+k24*k32*ka)*r2+(k24*k32*ka-k24*ka*lam2)*r1+A4last*lam24+(-A4last*ka-A4last*k32-A2last*k24-A4last*E2)*lam23+((A4last*k32+(A2last+A1last)*k24+A4last*E2)*ka+((A3last+A2last)*k24-A4last*k23+A4last*E2)*k32)*lam22+((-A3last-A2last-A1last)*k24+A4last*k23-A4last*E2)*k32*ka*lam2)*exp(-lam2*t))/((lam23+(-lam1-ka)*lam22+ka*lam1*lam2)*lam3-lam24+(lam1+ka)*lam23-ka*lam1*lam22);
      T a43 = (((k24*lam12+(-k24*ka-k24*k32)*lam1+k24*k32*ka)*r2+(k24*k32*ka-k24*ka*lam1)*r1+A4last*lam14+(-A4last*ka-A4last*k32-A2last*k24-A4last*E2)*lam13+((A4last*k32+(A2last+A1last)*k24+A4last*E2)*ka+((A3last+A2last)*k24-A4last*k23+A4last*E2)*k32)*lam12+((-A3last-A2last-A1last)*k24+A4last*k23-A4last*E2)*k32*ka*lam1)*exp(-lam1*t))/(((lam12-ka*lam1)*lam2-lam13+ka*lam12)*lam3+(ka*lam12-lam13)*lam2+lam14-ka*lam13);
      T a44 = (((k24*ka-k24*k32)*r1-A1last*k24*ka2+A1last*k24*k32*ka)*exp(-ka*t))/(((lam1-ka)*lam2-ka*lam1+ka2)*lam3+(ka2-ka*lam1)*lam2+ka2*lam1-ka3);
      T a45 = (k24*k32*r2+k24*k32*r1)/(lam1*lam2*lam3);
      A4=a41-a42+a43+a44+a45;
      return A;
    }

        // undefine extras
#undef tlag
#undef F
#undef rate1
#undef dur1
#undef ka
#undef tlag2
#undef f2
#undef dur2


        // one compartment ka translations ncmt=1
#define tlag  params(2, 0)
#define F     params(3, 0)
#define rate1 params(4, 0)
#define dur1  params(5, 0)
#define ka    params(6, 0)
#define tlag2 params(7, 0)
#define f2    params(8, 0)
#define dur2  params(9, 0)    

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    oneCmtKaSSb1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
		 T tau) {
      
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(2, 1);
      T eKa = 1.0/(1.0-exp(-tau*ka));
      T eK =  1.0/(1.0-exp(-tau*k20));
      A1=eKa*b1;
      A2=ka*b1*(eK - eKa)/(-k20 + ka);
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    oneCmtKaSSb2(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
		 T tau) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(2, 1);
      T eK =  1.0/(1.0-exp(-tau*k20));
      A1=0.0;
      A2=eK*b2;
      return A;
    }
    
    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    oneCmtKa(T t, Eigen::Matrix<T, Eigen::Dynamic, 2>& Alast,
	     Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
	     Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
	     Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(2, 1);
      T rx_expr_0=exp(-t*ka);
      A1=A1last*rx_expr_0+b1;
      T rx_expr_1=exp(-t*k20);
      A2=A1last*ka/(ka-k20)*(rx_expr_1-rx_expr_0)+A2last*rx_expr_1+b2;
      return A;
    }

    // Undefine extras
#undef tlag
#undef F
#undef rate1
#undef dur1
#undef ka
#undef tlag2
#undef f2
#undef dur2
    // two compartment ka translations ncmt=1
#define tlag  params(4,  0)
#define F     params(5,  0)
#define rate1 params(6,  0)
#define dur1  params(7,  0)
#define ka    params(8,  0)
#define tlag2 params(9,  0)
#define f2    params(10, 0)
#define dur2  params(11, 0)

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    twoCmtKaSSb1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
		 T tau) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(3, 1);
      T E2 = k20+k23;
      T E3 = k32;
      T e2e3 = E2+E3;
      T s = sqrt(e2e3*e2e3-4*(E2*E3-k23*k32));

      //calculate hybrid rate constants
      T lambda1 = 0.5*(e2e3+s);
      T lambda2 = 0.5*(e2e3-s);
      T eKa=1.0/(1.0-exp(-tau*ka));
      T eL1=1.0/(1.0-exp(-tau*lambda1));
      T eL2=1.0/(1.0-exp(-tau*lambda2));
      A1=eKa*b1;
      A2=ka*b1*(eL1*(E3 - lambda1)/((-lambda1 + lambda2)*(ka - lambda1)) + eL2*(E3 - lambda2)/((lambda1 - lambda2)*(ka - lambda2)) + eKa*(E3 - ka)/((-ka + lambda2)*(-ka + lambda1)));
      A3=ka*b1*k23*(eL1/((-lambda1 + lambda2)*(ka - lambda1)) + eL2/((lambda1 - lambda2)*(ka - lambda2)) + eKa/((-ka + lambda2)*(-ka + lambda1)));
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    twoCmtKaSSb2(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
		 T tau) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(3, 1);
      T E2 = k20+k23;
      T E3 = k32;
      T e2e3 = E2+E3;
      T s = sqrt(e2e3*e2e3-4*(E2*E3-k23*k32));

      //calculate hybrid rate constants
      T lambda1 = 0.5*(e2e3+s);
      T lambda2 = 0.5*(e2e3-s);
      /* T eKa=1.0/(1.0+exp(-tau*ka)); */
      T eL1=1.0/(1.0-exp(-tau*lambda1));
      T eL2=1.0/(1.0-exp(-tau*lambda2));

      A1=0.0;
      A2=(eL1*(b2*E3 - b2*lambda1) - eL2*(b2*E3 - b2*lambda2))/(-lambda1 + lambda2);
      A3=(eL1*b2*k23 - eL2*b2*k23)/(-lambda1 + lambda2);
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    twoCmtKa(T t, Eigen::Matrix<T, Eigen::Dynamic, 2>& Alast,
	     Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
	     Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
	     Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(3, 1);
      T rxe2=exp(-t*ka);
      A1=b1+rxe2*A1last;
      T rxe0=k12+k21;
      T rxe1=k12+kel;
      T rxe3=k21*A2last;
      T rxe4=k21*A3last;
      T rxe6=rxe0+kel;
      T rxe7=(rxe1)*k21;
      T rxe8=rxe6*rxe6;
      T rxe10=sqrt(-4*(-k12*k21+rxe7)+rxe8);
      A2=b2+(-exp(-0.5*t*(rxe6-rxe10))*(-0.5*A2last*(rxe6-rxe10)+rxe3+rxe4)+exp(-0.5*t*(rxe6+rxe10))*(-0.5*A2last*(rxe6+rxe10)+rxe3+rxe4))/(0.5*(rxe6-rxe10)-0.5*(rxe6+rxe10))+ka*(rxe2*(k21-ka)/((-ka+0.5*(rxe6-rxe10))*(-ka+0.5*(rxe6+rxe10)))+exp(-0.5*t*(rxe6-rxe10))*(k21-0.5*(rxe6-rxe10))/((-0.5*(rxe6-rxe10)+0.5*(rxe6+rxe10))*(ka-0.5*(rxe6-rxe10)))+exp(-0.5*t*(rxe6+rxe10))*(k21-0.5*(rxe6+rxe10))/((0.5*(rxe6-rxe10)-0.5*(rxe6+rxe10))*(ka-0.5*(rxe6+rxe10))))*A1last;
      T rxe5=k12*A2last;
      T rxe9=(rxe1)*A3last;
      A3=(-exp(-0.5*t*(rxe6-rxe10))*(-0.5*A3last*(rxe6-rxe10)+rxe5+rxe9)+exp(-0.5*t*(rxe6+rxe10))*(-0.5*A3last*(rxe6+rxe10)+rxe5+rxe9))/(0.5*(rxe6-rxe10)-0.5*(rxe6+rxe10))+ka*k12*A1last*(rxe2/((-ka+0.5*(rxe6-rxe10))*(-ka+0.5*(rxe6+rxe10)))+exp(-0.5*t*(rxe6-rxe10))/((-0.5*(rxe6-rxe10)+0.5*(rxe6+rxe10))*(ka-0.5*(rxe6-rxe10)))+exp(-0.5*t*(rxe6+rxe10))/((0.5*(rxe6-rxe10)-0.5*(rxe6+rxe10))*(ka-0.5*(rxe6+rxe10))));
      return A;
    }
    
    // undefine extras
#undef tlag
#undef F
#undef rate1
#undef dur1
#undef ka
#undef tlag2
#undef f2
#undef dur2
    // three compartment ka translations ncmt=1
#define tlag  params(6,  0)
#define F     params(7,  0)
#define rate1 params(8,  0)
#define dur1  params(9,  0)
#define ka    params(10, 0)
#define tlag2 params(11, 0)
#define f2    params(12, 0)
#define dur2  params(13, 0)
    
    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    threeCmtKaSSb1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
		   T tau){
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(4, 1);
      T E2 = k20+k23+k24;
      T E3 = k32;
      T E4 = k42;

      //calculate hybrid rate constants
      T a = E2+E3+E4;
      T b = E2*E3+E4*(E2+E3)-k23*k32-k24*k42;
      T c = E2*E3*E4-E4*k23*k32-E3*k24*k42;

      T a2 = a*a;
      T m = 0.333333333333333*(3.0*b - a2);
      T n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
      T Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

      T alpha = sqrt(-Q);
      T beta = -0.5*n;
      T gamma = sqrt(_as_zero(beta*beta+alpha*alpha));
      T theta = atan2(alpha,beta);
      T theta3 = 0.333333333333333*theta;
      T ctheta3 = cos(theta3);
      T stheta3 = 1.7320508075688771932*sin(theta3);
      T gamma3 = R_pow(gamma,0.333333333333333);
  
      T lambda1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
      T lambda2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
      T lambda3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);

      T eKa = 1.0/(1.0-exp(-tau*ka));
      T eL1 = 1.0/(1.0-exp(-tau*lambda1));
      T eL2 = 1.0/(1.0-exp(-tau*lambda2));
      T eL3 = 1.0/(1.0-exp(-tau*lambda3));
  
      A1=eKa*b1;
      A2=ka*b1*(eL1*(E3 - lambda1)*(E4 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*(ka - lambda1)) + eL2*(E3 - lambda2)*(E4 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*(ka - lambda2)) + eL3*(E3 - lambda3)*(E4 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*(ka - lambda3)) + eKa*(E3 - ka)*(E4 - ka)/((-ka + lambda1)*(-ka + lambda3)*(-ka + lambda2)));
      A3=ka*b1*k23*(eL1*(E4 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*(ka - lambda1)) + eL2*(E4 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*(ka - lambda2)) + eL3*(E4 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*(ka - lambda3)) + eKa*(E4 - ka)/((-ka + lambda1)*(-ka + lambda3)*(-ka + lambda2)));
      A4=ka*b1*k24*(eL1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*(ka - lambda1)) + eL2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*(ka - lambda2)) + eL3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*(ka - lambda3)) + eKa*(E3 - ka)/((-ka + lambda1)*(-ka + lambda3)*(-ka + lambda2)));
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    threeCmtKaSSb2(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
		   T tau) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(4, 1);
      T E2 = k20+k23+k24;
      T E3 = k32;
      T E4 = k42;

      //calculate hybrid rate constants
      T a = E2+E3+E4;
      T b = E2*E3+E4*(E2+E3)-k23*k32-k24*k42;
      T c = E2*E3*E4-E4*k23*k32-E3*k24*k42;

      T a2 = a*a;
      T m = 0.333333333333333*(3.0*b - a2);
      T n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
      T Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

      T alpha = sqrt(-Q);
      T beta = -0.5*n;
      T gamma = sqrt(_as_zero(beta*beta+alpha*alpha));
      T theta = atan2(alpha,beta);
      T theta3 = 0.333333333333333*theta;
      T ctheta3 = cos(theta3);
      T stheta3 = 1.7320508075688771932*sin(theta3);
      T gamma3 = R_pow(gamma,0.333333333333333);
  
      T lambda1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
      T lambda2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
      T lambda3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);

      /* T eKa = 1.0/(1.0-exp(-tau*KA)); */
      T eL1 = 1.0/(1.0-exp(-tau*lambda1));
      T eL2 = 1.0/(1.0-exp(-tau*lambda2));
      T eL3 = 1.0/(1.0-exp(-tau*lambda3));
      *A1=0.0;
      *A2=b2*(eL1*(E3 - lambda1)*(E4 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)) + eL2*(E3 - lambda2)*(E4 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)) + eL3*(E3 - lambda3)*(E4 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)));
      *A3=eL2*(-b2*E4*k23 + b2*k23*lambda2)/((lambda1 - lambda2)*(lambda2 - lambda3)) + eL1*(b2*E4*k23 - b2*k23*lambda1)/((lambda1 - lambda3)*(lambda1 - lambda2)) + eL3*(-b2*E4*k23 + b2*k23*lambda3)/((lambda1 - lambda3)*(-lambda2 + lambda3));
      *A4=eL2*(-b2*E3*k24 + b2*k24*lambda2)/((lambda1 - lambda2)*(lambda2 - lambda3)) + eL1*(b2*E3*k24 - b2*k24*lambda1)/((lambda1 - lambda3)*(lambda1 - lambda2)) + eL3*(-b2*E3*k24 + b2*k24*lambda3)/((lambda1 - lambda3)*(-lambda2 + lambda3));
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    threeCmtKa(T t, Eigen::Matrix<T, Eigen::Dynamic, 2>& Alast,
	       Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
	       Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
	       Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(4, 1);
      T E2 = k20+k23+k24;
      T E3 = k32;
      T E4 = k42;

      //calculate hybrid rate constants
      T a = E2+E3+E4;
      T b = E2*E3+E4*(E2+E3)-k23*k32-k24*k42;
      T c = E2*E3*E4-E4*k23*k32-E3*k24*k42;

      T a2 = a*a;
      T m = 0.333333333333333*(3.0*b - a2);
      T n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
      T Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

      T alpha = sqrt(-Q);
      T beta = -0.5*n;
      T gamma = sqrt(_as_zero(beta*beta+alpha*alpha));
      T theta = atan2(alpha,beta);
      T theta3 = 0.333333333333333*theta;
      T ctheta3 = cos(theta3);
      T stheta3 = 1.7320508075688771932*sin(theta3);
      T gamma3 = R_pow(gamma,0.333333333333333);
  
      T lambda1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
      T lambda2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
      T lambda3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);

      T B = A3last*k32+A4last*k42;
      T C = E4*A3last*k32+E3*A4last*k42;
      T I = A2last*k23*E4-A3last*k24*k42+A4last*k23*k42;
      T J = A2last*k24*E3+A3last*k24*k32-A4last*k23*k32;

      T eL1 = exp(-t*lambda1);
      T eL2 = exp(-t*lambda2);
      T eL3 = exp(-t*lambda3);
      T eKA = exp(-t*ka);

      T l12 = (lambda1-lambda2);
      T l13 = (lambda1-lambda3);
      T l21 = (lambda2-lambda1);
      T l23 = (lambda2-lambda3);
      T l31 = (lambda3-lambda1);
      T l32 = (lambda3-lambda2);

      T e2l1 = (E2-lambda1);
      T e2l2 = (E2-lambda2);
      T e3l1 = (E3-lambda1);
      T e3l2 = (E3-lambda2);
      T e3l3 = (E3-lambda3);
      T e4l1 = (E4-lambda1);
      T e4l2 = (E4-lambda2);
      T e4l3 = (E4-lambda3);
  
      T A2term1 = A2last*(eL1*e3l1*e4l1/(l21*l31)+eL2*e3l2*e4l2/(l12*l32)+eL3*e3l3*e4l3/(l13*l23));
  
      T A2term2 = eL1*(C-B*lambda1)/(l12*l13)+eL2*(B*lambda2-C)/(l12*l23)+eL3*(B*lambda3-C)/(l13*l32);
  
      T A2term3 = A1last*ka*(eL1*e3l1*e4l1/(l21*l31*(ka-lambda1))+eL2*e3l2*e4l2/(l12*l32*(ka-lambda2))+eL3*e3l3*e4l3/(l13*l23*(ka-lambda3))+eKA*(E3-ka)*(E4-ka)/((lambda1-ka)*(lambda2-ka)*(lambda3-ka)));
  
      A2 = A2term1+A2term2+A2term3 + b2;

      T A3term1 = A3last*(eL1*e2l1*e4l1/(l21*l31)+eL2*e2l2*e4l2/(l12*l32)+eL3*(E2-lambda3)*e4l3/(l13*l23));
  
      T A3term2 = eL1*(I-A2last*k23*lambda1)/(l12*l13)+eL2*(A2last*k23*lambda2-I)/(l12*l23)+eL3*(A2last*k23*lambda3-I)/(l13*l32);
  
      T A3term3 = A1last*ka*k23*(eL1*e4l1/(l21*l31*(ka-lambda1))+eL2*e4l2/(l12*l32*(ka-lambda2))+eL3*e4l3/(l13*l23*(ka-lambda3))+eKA*(E4-ka)/((lambda1-ka)*(lambda2-ka)*(lambda3-ka)));
  
      A3 = A3term1+A3term2+A3term3;// Amount in the first-peripheral compartment

      T A4term1 = A4last*(eL1*e2l1*e3l1/(l21*l31)+eL2*e2l2*e3l2/(l12*l32)+eL3*(E2-lambda3)*e3l3/(l13*l23));
  
      T A4term2 = eL1*(J-A2last*k24*lambda1)/(l12*l13)+eL2*(A2last*k24*lambda2-J)/(l12*l23)+eL3*(A2last*k24*lambda3-J)/(l13*l32);
  
      T A4term3 = A1last*ka*k24*(eL1*e3l1/(l21*l31*(ka-lambda1))+eL2*e3l2/(l12*l32*(ka-lambda2))+eL3*e3l3/(l13*l23*(ka-lambda3))+eKA*(E3-ka)/((lambda1-ka)*(lambda2-ka)*(lambda3-ka)));
      A4 = A4term1+A4term2+A4term3;

      A1 = A1last*eKA + b1;
      return A;
    }

            // undefine extras
#undef tlag
#undef F
#undef rate1
#undef dur1
#undef ka
#undef tlag2
#undef f2
#undef dur2


        // one compartment ka translations ncmt=1
#define tlag  params(2, 0)
#define F     params(3, 0)
#define rate1 params(4, 0)
#define dur1  params(5, 0)
#define ka    params(6, 0)
#define tlag2 params(7, 0)
#define f2    params(8, 0)
#define dur2  params(9, 0)

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    oneCmtRateSSr1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(1, 1);
      A1 = r1/k10;
      return A;
    }
    
    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    oneCmtRateSS(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& rate,
		 T tinf, T tau) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(1, 1);
      T eiK = exp(-k10*tinf);
      T eK = exp(-k10*(tau-tinf))/(1.0-exp(-k10*tau));
      A1=r1*(1-eiK)*eK/(k10);
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    oneCmtRate(T t, Eigen::Matrix<T, Eigen::Dynamic, 2>& Alast,
	       Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
	       Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
	       Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
	       Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(1, 1);
      T eT = exp(-k10*t);
      A1 = r1/k10*(1-eT)+A1last*eT + b1;
      return A;
    }
        // Undefine extras
#undef tlag
#undef F
#undef rate1
#undef dur1
#undef ka
#undef tlag2
#undef f2
#undef dur2
    // two compartment ka translations ncmt=1
#define tlag  params(4,  0)
#define F     params(5,  0)
#define rate1 params(6,  0)
#define dur1  params(7,  0)
#define ka    params(8,  0)
#define tlag2 params(9,  0)
#define f2    params(10, 0)
#define dur2  params(11, 0)

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    twoCmtRateSSr1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(2, 1);
      T E1 = k10+k12;
      T E2 = k21;
      //#calculate hybrid rate constants
      T s = E1+E2;
      T sqr = sqrt(s*s-4*(E1*E2-k12*k21));
      T lambda1 = 0.5*(s+sqr);
      T lambda2 = 0.5*(s-sqr);
      T l12 = 1.0/(lambda1*lambda2);
      A1=r1*E2*l12;
      A2=r1*k12*l12;
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    twoCmtRateSS(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& rate,
		 T tinf, T tau) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(2, 1);
      T E1 = k10+k12;
      T E2 = k21;

      //#calculate hybrid rate constants
      T s = E1+E2;
      T sqr = sqrt(s*s-4*(E1*E2-k12*k21));
      T lambda1 = 0.5*(s+sqr);
      T lambda2 = 0.5*(s-sqr);

      T eTi1 = exp(-tinf*lambda1);
      T eTi2 = exp(-tinf*lambda2);
      T eT1 =exp(-lambda1*(tau-tinf))/(1.0-exp(-tau*lambda1));
      T eT2 =exp(-lambda2*(tau-tinf))/(1.0-exp(-tau*lambda2));
      A1=(eT1*(E2*((eTi1*r1 - eTi2*r1)/(-lambda1 + lambda2) + r1*E2*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2))) - lambda1*((eTi1*r1 - eTi2*r1)/(-lambda1 + lambda2) + r1*E2*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2))) + r1*k12*k21*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2))) - eT2*(E2*((eTi1*r1 - eTi2*r1)/(-lambda1 + lambda2) + r1*E2*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2))) - lambda2*((eTi1*r1 - eTi2*r1)/(-lambda1 + lambda2) + r1*E2*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2))) + r1*k12*k21*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2))))/(-lambda1 + lambda2);
      A2=(eT1*(k12*((eTi1*r1 - eTi2*r1)/(-lambda1 + lambda2) + r1*E2*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2))) + r1*E1*k12*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2)) - r1*k12*lambda1*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2))) - eT2*(k12*((eTi1*r1 - eTi2*r1)/(-lambda1 + lambda2) + r1*E2*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2))) + r1*E1*k12*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2)) - r1*k12*lambda2*(1.0/(lambda1*lambda2) + eTi1/((lambda1 - lambda2)*lambda1) - eTi2/((lambda1 - lambda2)*lambda2))))/(-lambda1 + lambda2);
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    twoCmtRate(T t, Eigen::Matrix<T, Eigen::Dynamic, 2>& Alast,
	       Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
	       Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
	       Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
	       Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(2, 1);
      T E1 = k10+k12;
      T E2 = k21;

      //#calculate hybrid rate constants
      T s = E1+E2;
      T sqr = sqrt(s*s-4*(E1*E2-k12*k21));
      T lambda1 = 0.5*(s+sqr);
      T lambda2 = 0.5*(s-sqr);

      T eT1 = exp(-t*lambda1);
      T eT2 = exp(-t*lambda2);

      T l12 = (lambda1-lambda2);
      T l21 = (lambda2-lambda1);

      T c10 = (A1last*E2+Doserate+A2last*k21);
      T c11 = (c10-A1last*lambda1)/l21;
      T c12 = (c10-A1last*lambda2)/l21;
      T A1term1 = c11*eT1 - c12*eT2;
      T A1term2 = Doserate*E2*(1/(lambda1*lambda2)+eT1/(lambda1*l12)-eT2/(lambda2*l12));
      A1 = A1term1+A1term2 + b1;//Amount in the central compartment

      T c20 = (A2last*E1+A1last*k12);
      T c21 = (c20-A2last*lambda1)/l21;
      T c22 = (c20-A2last*lambda2)/l21;
      T A2term1 = c21*eT1-c22*eT2;
      T A2term2 = Doserate*k12*(1/(lambda1*lambda2)+eT1/(lambda1*l12)-eT2/(lambda2*(lambda1-lambda2)));
      A2 = A2term1+A2term2;//Amount in the peripheral compartment
      return A;
    }

    // undefine extras
#undef tlag
#undef F
#undef rate1
#undef dur1
#undef ka
#undef tlag2
#undef f2
#undef dur2
    // three compartment ka translations ncmt=1
#define tlag  params(6,  0)
#define F     params(7,  0)
#define rate1 params(8,  0)
#define dur1  params(9,  0)
#define ka    params(10, 0)
#define tlag2 params(11, 0)
#define f2    params(12, 0)
#define dur2  params(13, 0)

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    threeCmtRateSSr1(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		     Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		     Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(3, 1);
      T E1 = k10+k12+k13;
      T E2 = k21;
      T E3 = k31;

      //#calculate hybrid rate constants
      T a = E1+E2+E3;
      T b = E1*E2+E3*(E1+E2)-k12*k21-k13*k31;
      T c = E1*E2*E3-E3*k12*k21-E2*k13*k31;

      T a2 = a*a;
      T m = 0.333333333333333*(3.0*b - a2);
      T n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
      T Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

      T alpha = sqrt(-Q);
      T beta = -0.5*n;
      T gamma = sqrt(_as_zero(beta*beta+alpha*alpha));
      T theta = atan2(alpha,beta);
      T theta3 = 0.333333333333333*theta;
      T ctheta3 = cos(theta3);
      T stheta3 = 1.7320508075688771932*sin(theta3);
      T gamma3 = R_pow(gamma,0.333333333333333);
  
      T lambda1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
      T lambda2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
      T lambda3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);
      T l123 = 1.0/(lambda1*lambda2*lambda3);
      A1=r1*E2*E3*l123;
      A2=r1*E3*k12*l123;
      A3=r1*E2*k13*l123;
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    threeCmtRateSS(Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		   Eigen::Matrix<T, Eigen::Dynamic, 2>& rate,
		   T tinf, T tau){
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(3, 1);
      T E1 = k10+k12+k13;
      T E2 = k21;
      T E3 = k31;

      //#calculate hybrid rate constants
      T a = E1+E2+E3;
      T b = E1*E2+E3*(E1+E2)-k12*k21-k13*k31;
      T c = E1*E2*E3-E3*k12*k21-E2*k13*k31;

      T a2 = a*a;
      T m = 0.333333333333333*(3.0*b - a2);
      T n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
      T Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

      T alpha = sqrt(-Q);
      T beta = -0.5*n;
      T gamma = sqrt(_as_zero(beta*beta+alpha*alpha));
      T theta = atan2(alpha,beta);
      T theta3 = 0.333333333333333*theta;
      T ctheta3 = cos(theta3);
      T stheta3 = 1.7320508075688771932*sin(theta3);
      T gamma3 = R_pow(gamma,0.333333333333333);
  
      T lambda1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
      T lambda2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
      T lambda3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);
      T eTi1 = exp(-tinf*lambda1);
      T eTi2 = exp(-tinf*lambda2);
      T eTi3 = exp(-tinf*lambda3);
      T eT1 = exp(-lambda1*(tau-tinf))/(1.0-exp(-tau*lambda1));
      T eT2 = exp(-lambda2*(tau-tinf))/(1.0-exp(-tau*lambda2));
      T eT3 = exp(-lambda3*(tau-tinf))/(1.0-exp(-tau*lambda3));
      A1=r1*(eT1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)) + eT2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)) + eT3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)))*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + eT2*(lambda2*(r1*k12*k21*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + r1*k13*k31*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))) - (r1*E2*k13*k31*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + r1*E3*k12*k21*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))))/((lambda1 - lambda2)*(lambda2 - lambda3)) + eT1*(-lambda1*(r1*k12*k21*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + r1*k13*k31*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))) + r1*E2*k13*k31*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + r1*E3*k12*k21*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)))/((lambda1 - lambda3)*(lambda1 - lambda2)) + eT3*(lambda3*(r1*k12*k21*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + r1*k13*k31*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))) - (r1*E2*k13*k31*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + r1*E3*k12*k21*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))))/((lambda1 - lambda3)*(-lambda2 + lambda3));
      A2=r1*k12*(eT1*(E3 - lambda1)*(E1 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)) + eT2*(E1 - lambda2)*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)) + eT3*(E1 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)))*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + eT2*(r1*k12*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))*lambda2 - (r1*E3*k12*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + r1*k13*k12*k31*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) - r1*k13*k12*k31*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))))/((lambda1 - lambda2)*(lambda2 - lambda3)) + eT1*(r1*E3*k12*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) - r1*k12*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))*lambda1 + r1*k13*k12*k31*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) - r1*k13*k12*k31*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)))/((lambda1 - lambda3)*(lambda1 - lambda2)) + eT3*(r1*k12*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))*lambda3 - (r1*E3*k12*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + r1*k13*k12*k31*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) - r1*k13*k12*k31*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))))/((lambda1 - lambda3)*(-lambda2 + lambda3));
      A3=r1*k13*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))*(eT1*(E2 - lambda1)*(E1 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)) + eT2*(E1 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)) + eT3*(E1 - lambda3)*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3))) + eT2*(r1*k13*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))*lambda2 - (r1*E2*k13*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) - r1*k13*k12*k21*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + r1*k13*k12*k21*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))))/((lambda1 - lambda2)*(lambda2 - lambda3)) + eT1*(r1*E2*k13*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) - r1*k13*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))*lambda1 - r1*k13*k12*k21*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + r1*k13*k12*k21*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)))/((lambda1 - lambda3)*(lambda1 - lambda2)) + eT3*(r1*k13*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))*lambda3 - (r1*E2*k13*(E2*E3/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) - r1*k13*k12*k21*(E2/(lambda1*lambda2*lambda3) - eTi1*(E2 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E2 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E2 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3)) + r1*k13*k12*k21*(E3/(lambda1*lambda2*lambda3) - eTi1*(E3 - lambda1)/((-lambda1 + lambda3)*(-lambda1 + lambda2)*lambda1) - eTi2*(E3 - lambda2)/((lambda1 - lambda2)*(-lambda2 + lambda3)*lambda2) - eTi3*(E3 - lambda3)/((lambda1 - lambda3)*(lambda2 - lambda3)*lambda3))))/((lambda1 - lambda3)*(-lambda2 + lambda3));
      return A;
    }

    template <class T>
    Eigen::Matrix<T, Eigen::Dynamic, 1>
    threeCmtRate(T t, Eigen::Matrix<T, Eigen::Dynamic, 2>& Alast,
		 Eigen::Matrix<T, Eigen::Dynamic, 1>& params,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& g,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& bolus,
		 Eigen::Matrix<T, Eigen::Dynamic, 2>& rate) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> A(3, 1);
      T E1 = k10+k12+k13;
      T E2 = k21;
      T E3 = k31;

      //#calculate hybrid rate constants
      T a = E1+E2+E3;
      T b = E1*E2+E3*(E1+E2)-k12*k21-k13*k31;
      T c = E1*E2*E3-E3*k12*k21-E2*k13*k31;

      T a2 = a*a;
      T m = 0.333333333333333*(3.0*b - a2);
      T n = 0.03703703703703703*(2.0*a2*a - 9.0*a*b + 27.0*c);
      T Q = 0.25*(n*n) + 0.03703703703703703*(m*m*m);

      T alpha = sqrt(-Q);
      T beta = -0.5*n;
      T gamma = sqrt(_as_zero(beta*beta+alpha*alpha));
      T theta = atan2(alpha,beta);
      T theta3 = 0.333333333333333*theta;
      T ctheta3 = cos(theta3);
      T stheta3 = 1.7320508075688771932*sin(theta3);
      T gamma3 = R_pow(gamma,0.333333333333333);
  
      T lambda1 = 0.333333333333333*a + gamma3*(ctheta3 + stheta3);
      T lambda2 = 0.333333333333333*a + gamma3*(ctheta3 -stheta3);
      T lambda3 = 0.333333333333333*a -(2.0*gamma3*ctheta3);

      T B = A2last*k21+A3last*k31;
      T C = E3*A2last*k21+E2*A3last*k31;
      T I = A1last*k12*E3-A2last*k13*k31+A3last*k12*k31;
      T J = A1last*k13*E2+A2last*k13*k21-A3last*k12*k21;

      T eL1 = exp(-t*lambda1);
      T eL2 = exp(-t*lambda2);
      T eL3 = exp(-t*lambda3);

      T l12 = (lambda1-lambda2);
      T l13 = (lambda1-lambda3);
      T l21 = (lambda2-lambda1);
      T l23 = (lambda2-lambda3);
      T l31 = (lambda3-lambda1);
      T l32 = (lambda3-lambda2);
  
      T e1l1 = (E1-lambda1);
      T e1l2 = (E1-lambda2);
      T e1l3 = (E1-lambda3);
      T e2l1 = (E2-lambda1);
      T e2l2 = (E2-lambda2);
      T e2l3 = (E2-lambda3);
      T e3l1 = (E3-lambda1);
      T e3l2 = (E3-lambda2);
      T e3l3 = (E3-lambda3);

      T A1term1 = A1last*(eL1*e2l1*e3l1/(l21*l31)+eL2*e2l2*e3l2/(l12*l32)+eL3*e2l3*e3l3/(l13*l23));
      T A1term2 = eL1*(C-B*lambda1)/(l12*l13)+eL2*(B*lambda2-C)/(l12*l23)+eL3*(B*lambda3-C)/(l13*l32);
      T A1term3 = Doserate*((E2*E3)/(lambda1*lambda2*lambda3)-eL1*e2l1*e3l1/(lambda1*l21*l31)-eL2*e2l2*e3l2/(lambda2*l12*l32)-eL3*e2l3*e3l3/(lambda3*l13*l23));

      A1 = A1term1+A1term2+A1term3 + b1;//Amount in the central compartment

      T A2term1 = A2last*(eL1*e1l1*e3l1/(l21*l31)+eL2*e1l2*e3l2/(l12*l32)+eL3*e1l3*e3l3/(l13*l23));
      T A2term2 = eL1*(I-A1last*k12*lambda1)/(l12*l13)+eL2*(A1last*k12*lambda2-I)/(l12*l23)+eL3*(A1last*k12*lambda3-I)/(l13*l32);
      T A2term3 = Doserate*k12*(E3/(lambda1*lambda2*lambda3)-eL1*e3l1/(lambda1*l21*l31)-eL2*e3l2/(lambda2*l12*l32)-eL3*e3l3/(lambda3*l13*l23));

      A2 = A2term1+A2term2+A2term3;// Amount in the first-peripheral compartment

      T A3term1 = A3last*(eL1*e1l1*e2l1/(l21*l31)+eL2*e1l2*e2l2/(l12*l32)+eL3*e1l3*e2l3/(l13*l23));
      T A3term2 = eL1*(J-A1last*k13*lambda1)/(l12*l13)+eL2*(A1last*k13*lambda2-J)/(l12*l23)+eL3*(A1last*k13*lambda3-J)/(l13*l32);
      T A3term3 = Doserate*k13*(E2/(lambda1*lambda2*lambda3)-eL1*e2l1/(lambda1*l21*l31)-eL2*e2l2/(lambda2*l12*l32)-eL3*e2l3/(lambda3*l13*l23));

      A3 = A3term1+A3term2+A3term3;//Amount in the second-peripheral compartment
      return A;
    }

#undef v
#undef k20
#undef kel
#undef k23
#undef k32
#undef k24
#undef k42
#undef k10
#undef k12
#undef k21
#undef k13
#undef k31
#undef A1
#undef A2
#undef A3
#undef A3
#undef r1
#undef Doserate
#undef r2
#undef b1
#undef b2
    // undefine extras
#undef tlag
#undef F
#undef rate1
#undef dur1
#undef ka
#undef tlag2
#undef f2
#undef dur2
    
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
  params(0, 0) = dd_p1;
  params(1, 0) = dd_v1;
  if (ncmt >=2 ){
    params(2,0) = dd_p2;
    params(3,0) = dd_p3;
    if (ncmt >= 3){
      params(4,0) = dd_p4;
      params(5,0) = dd_p5;
    }
  }
  params(2*ncmt,     0) = dd_tlag;
  params(2*ncmt + 1, 0) = dd_F;
  params(2*ncmt + 2, 0) = dd_rate;
  params(2*ncmt + 3, 0) = dd_dur;
  if (oral0) {
    params(2*ncmt + 4, 0) = dd_ka;
    params(2*ncmt + 5, 0) = dd_tlag2;
    params(2*ncmt + 6, 0) = dd_F2;
    params(2*ncmt + 7, 0) = dd_rate2;
    params(2*ncmt + 8, 0) = dd_dur2;
  }
  
  return 0.0;
}
