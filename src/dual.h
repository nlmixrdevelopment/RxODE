#ifndef dual_header
#define dual_header
#include <math.h>

#define dKa 0
#define dP1 1
#define dP2 2
#define dP3 3
#define dP4 4
#define dP5 5
#define dV1 6

void iniD(double val, int which, dualN *ret){
  ret->f = val;
  ret->grad[0] = 0.0;
  ret->grad[1] = 0.0;
  ret->grad[2] = 0.0;
  ret->grad[3] = 0.0;
  ret->grad[4] = 0.0;
  ret->grad[5] = 0.0;
  ret->grad[6] = 0.0;
  if (which >= 0) {
    ret->grad[which] = 1.0;
  }
}

#define iniD0() iniD(NAN, -1)

dualN sqrtD(dualN x) {
  dualN ret;
  ret.f = sqrt(x.f);
  double gr = 0.5 / ret.f;
  ret.grad[0] = x.grad[0]*gr;
  ret.grad[1] = x.grad[1]*gr;
  ret.grad[2] = x.grad[2]*gr;
  ret.grad[3] = x.grad[3]*gr;
  ret.grad[4] = x.grad[4]*gr;
  ret.grad[5] = x.grad[5]*gr;
  ret.grad[6] = x.grad[6]*gr;
  return ret;
}

dualN expD(dualN x) {
  dualN ret;
  ret.f = exp(x.f);
  ret.grad[0] = x.grad[0]*ret.f;
  ret.grad[1] = x.grad[1]*ret.f;
  ret.grad[2] = x.grad[2]*ret.f;
  ret.grad[3] = x.grad[3]*ret.f;
  ret.grad[4] = x.grad[4]*ret.f;
  ret.grad[5] = x.grad[5]*ret.f;
  ret.grad[6] = x.grad[6]*ret.f;
  return ret;
}

dualN add2(dualN x, dualN y) {
  dualN ret;
  ret.f = x.f + y.f;
  ret.grad[0] = x.grad[0] + y.grad[0];
  ret.grad[1] = x.grad[1] + y.grad[1];
  ret.grad[2] = x.grad[2] + y.grad[2];
  ret.grad[3] = x.grad[3] + y.grad[3];
  ret.grad[4] = x.grad[4] + y.grad[4];
  ret.grad[5] = x.grad[5] + y.grad[5];
  ret.grad[6] = x.grad[6] + y.grad[6];
  return ret;
}

dualN add2d(dualN x, double y) {
  dualN ret;
  ret.f = x.f + y;
  ret.grad[0] = x.grad[0];
  ret.grad[1] = x.grad[1];
  ret.grad[2] = x.grad[2];
  ret.grad[3] = x.grad[3];
  ret.grad[4] = x.grad[4];
  ret.grad[5] = x.grad[5];
  ret.grad[6] = x.grad[6];
  return ret;
}

#define addd2(x, y) add2d(y, x)

dualN subtr2(dualN x, dualN y) {
  dualN ret;
  ret.f = x.f - y.f;
  ret.grad[0] = x.grad[0] - y.grad[0];
  ret.grad[1] = x.grad[1] - y.grad[1];
  ret.grad[2] = x.grad[2] - y.grad[2];
  ret.grad[3] = x.grad[3] - y.grad[3];
  ret.grad[4] = x.grad[4] - y.grad[4];
  ret.grad[5] = x.grad[5] - y.grad[5];
  ret.grad[6] = x.grad[6] - y.grad[6];
  return ret;
}

dualN negD(dualN x) {
  dualN ret;
  ret.f = -x.f;
  ret.grad[0] = -x.grad[0];
  ret.grad[1] = -x.grad[1];
  ret.grad[2] = -x.grad[2];
  ret.grad[3] = -x.grad[3];
  ret.grad[4] = -x.grad[4];
  ret.grad[5] = -x.grad[5];
  ret.grad[6] = -x.grad[6];
  return ret;
}

dualN subtrd2(double x, dualN y) {
  dualN ret;
  ret.f = x - y.f;
  ret.grad[0] = -y.grad[0];
  ret.grad[1] = -y.grad[1];
  ret.grad[2] = -y.grad[2];
  ret.grad[3] = -y.grad[3];
  ret.grad[4] = -y.grad[4];
  ret.grad[5] = -y.grad[5];
  ret.grad[6] = -y.grad[6];
  return ret;
}

dualN subtr2d(dualN x, double y) {
  dualN ret;
  ret.f = x.f - y;
  ret.grad[0] = x.grad[0];
  ret.grad[1] = x.grad[1];
  ret.grad[2] = x.grad[2];
  ret.grad[3] = x.grad[3];
  ret.grad[4] = x.grad[4];
  ret.grad[5] = x.grad[5];
  ret.grad[6] = x.grad[6];
  return ret;
}

dualN prod2(dualN e1, dualN e2) {
  dualN ret;
  ret.f = e1.f * e2.f;
  ret.grad[0] = e1.grad[0] * e2.f + e1.f * e2.grad[0];
  ret.grad[1] = e1.grad[1] * e2.f + e1.f * e2.grad[1];
  ret.grad[2] = e1.grad[2] * e2.f + e1.f * e2.grad[2];
  ret.grad[3] = e1.grad[3] * e2.f + e1.f * e2.grad[3];
  ret.grad[4] = e1.grad[4] * e2.f + e1.f * e2.grad[4];
  ret.grad[5] = e1.grad[5] * e2.f + e1.f * e2.grad[5];
  ret.grad[6] = e1.grad[6] * e2.f + e1.f * e2.grad[6];
  return ret;
}

dualN prodd2(double e1, dualN e2) {
  dualN ret;
  ret.f = e2.f * e1;
  ret.grad[0] = e2.grad[0]*e1;
  ret.grad[1] = e2.grad[1]*e1;
  ret.grad[2] = e2.grad[2]*e1;
  ret.grad[3] = e2.grad[3]*e1;
  ret.grad[4] = e2.grad[4]*e1;
  ret.grad[5] = e2.grad[5]*e1;
  ret.grad[6] = e2.grad[6]*e1;
  return ret;
}

#define prod2d(x, y) prodd2(y, x)

dualN div2(dualN e1, dualN e2) {
  dualN ret;
  double invE2 = 1.0 / e2.f;
  ret.f = e1.f*invE2;
  ret.grad[0] = (e1.grad[0] - ret.f * e2.grad[0]) * invE2;
  ret.grad[1] = (e1.grad[1] - ret.f * e2.grad[1]) * invE2;
  ret.grad[2] = (e1.grad[2] - ret.f * e2.grad[2]) * invE2;
  ret.grad[3] = (e1.grad[3] - ret.f * e2.grad[3]) * invE2;
  ret.grad[4] = (e1.grad[4] - ret.f * e2.grad[4]) * invE2;
  ret.grad[5] = (e1.grad[5] - ret.f * e2.grad[5]) * invE2;
  ret.grad[6] = (e1.grad[6] - ret.f * e2.grad[6]) * invE2;
  return ret;
}

dualN divd2(double e1, dualN e2) {
  dualN ret;
  double invE2 = 1.0 / e2.f;
  ret.f = e1*invE2;
  ret.grad[0] = ret.f * invE2 * e2.grad[0];
  ret.grad[1] = ret.f * invE2 * e2.grad[1];
  ret.grad[2] = ret.f * invE2 * e2.grad[2];
  ret.grad[3] = ret.f * invE2 * e2.grad[3];
  ret.grad[4] = ret.f * invE2 * e2.grad[4];
  ret.grad[5] = ret.f * invE2 * e2.grad[5];
  ret.grad[6] = ret.f * invE2 * e2.grad[6];
  return ret;
}

dualN div2d(dualN e1, double e2) {
  dualN ret;
  double invE2 = 1.0 / e2;
  ret.f = e1.f*invE2;
  ret.grad[0] = invE2 * e1.grad[0];
  ret.grad[1] = invE2 * e1.grad[1];
  ret.grad[2] = invE2 * e1.grad[2];
  ret.grad[3] = invE2 * e1.grad[3];
  ret.grad[4] = invE2 * e1.grad[4];
  ret.grad[5] = invE2 * e1.grad[5];
  ret.grad[5] = invE2 * e1.grad[5];
  ret.grad[6] = invE2 * e1.grad[6];
  return ret;
}


dualN atan2D(dualN y, dualN x) {
  double derg = y.f * y.f;
  double derf = x.f * x.f;
  derg = 1.0 / (derg + derf);
  derf = x.f * derg;
  derg = y.f * derg;
  dualN ret;
  ret.f = atan2(y.f, x.f);
  ret.grad[0] = y.grad[0] * derf - derg * x.grad[0];
  ret.grad[1] = y.grad[1] * derf - derg * x.grad[1];
  ret.grad[2] = y.grad[2] * derf - derg * x.grad[2];
  ret.grad[3] = y.grad[3] * derf - derg * x.grad[3];
  ret.grad[4] = y.grad[4] * derf - derg * x.grad[4];
  ret.grad[5] = y.grad[5] * derf - derg * x.grad[5];
  ret.grad[6] = y.grad[6] * derf - derg * x.grad[6];
  return ret;
}

dualN cosD(dualN x){
  //dual(f = cos(x@f), grad = -x@grad * sin(x@f))
  dualN ret;
  ret.f = cos(x.f);
  double sf = sin(x.f);
  ret.grad[0] = x.grad[0] * sf;
  ret.grad[1] = x.grad[1] * sf;
  ret.grad[2] = x.grad[2] * sf;
  ret.grad[3] = x.grad[3] * sf;
  ret.grad[4] = x.grad[4] * sf;
  ret.grad[5] = x.grad[5] * sf;
  ret.grad[6] = x.grad[6] * sf;

  return ret;
}

dualN sinD(dualN x) {
  dualN ret;
  ret.f = sin(x.f);
  double cs = cos(x.f);
  ret.grad[0] = x.grad[0] * cs;
  ret.grad[1] = x.grad[1] * cs;
  ret.grad[2] = x.grad[2] * cs;
  ret.grad[3] = x.grad[3] * cs;
  ret.grad[4] = x.grad[4] * cs;
  ret.grad[5] = x.grad[5] * cs;
  ret.grad[6] = x.grad[6] * cs;
  return ret;
}

dualN pow2d(dualN base, double e) {
  dualN ret;
  ret.f = R_pow(base.f, e);
  double mult = e*pow(base.f, e - 1.0);
  ret.grad[0] = base.grad[0] * mult;
  ret.grad[1] = base.grad[1] * mult;
  ret.grad[2] = base.grad[2] * mult;
  ret.grad[3] = base.grad[3] * mult;
  ret.grad[4] = base.grad[4] * mult;
  ret.grad[5] = base.grad[5] * mult;
  ret.grad[6] = base.grad[6] * mult;
  return ret;
}

dualN pow2i(dualN base, int e) {
  dualN ret;
  ret.f = R_pow_di(base.f, e);
  double mult = (double)(e)*R_pow_di(base.f, e - 1);
  ret.grad[0] = base.grad[0] * mult;
  ret.grad[1] = base.grad[1] * mult;
  ret.grad[2] = base.grad[2] * mult;
  ret.grad[3] = base.grad[3] * mult;
  ret.grad[4] = base.grad[4] * mult;
  ret.grad[5] = base.grad[5] * mult;
  ret.grad[6] = base.grad[6] * mult;
  return ret;
}

#endif
