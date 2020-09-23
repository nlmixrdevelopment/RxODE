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

typedef struct dualN {
  double f;
  double grad[7];
} dualN;

void assignD(dualN *out, dualN *in) {
  out->f = in->f;
  out->grad[0] = in->grad[0];
  out->grad[1] = in->grad[1];
  out->grad[2] = in->grad[2];
  out->grad[3] = in->grad[3];
  out->grad[4] = in->grad[4];
  out->grad[5] = in->grad[5];
  out->grad[6] = in->grad[6];
}

dualN *iniD(double val, int which, dualN *ret){
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
  return ret;
}

dualN *sqrtD(dualN *x, dualN *ret) {
  ret->f = sqrt(x->f);
  double gr = 0.5 / ret->f;
  ret->grad[0] = x->grad[0]*gr;
  ret->grad[1] = x->grad[1]*gr;
  ret->grad[2] = x->grad[2]*gr;
  ret->grad[3] = x->grad[3]*gr;
  ret->grad[4] = x->grad[4]*gr;
  ret->grad[5] = x->grad[5]*gr;
  ret->grad[6] = x->grad[6]*gr;
  return ret;
}

dualN *sqrtDX(dualN *x) {
  x->f = sqrt(x->f);
  double gr = 0.5 / x->f;
  x->grad[0] = x->grad[0]*gr;
  x->grad[1] = x->grad[1]*gr;
  x->grad[2] = x->grad[2]*gr;
  x->grad[3] = x->grad[3]*gr;
  x->grad[4] = x->grad[4]*gr;
  x->grad[5] = x->grad[5]*gr;
  x->grad[6] = x->grad[6]*gr;
  return x;
}

dualN *expD(dualN *x, dualN *ret) {
  ret->f = exp(x->f);
  ret->grad[0] = x->grad[0]*ret->f;
  ret->grad[1] = x->grad[1]*ret->f;
  ret->grad[2] = x->grad[2]*ret->f;
  ret->grad[3] = x->grad[3]*ret->f;
  ret->grad[4] = x->grad[4]*ret->f;
  ret->grad[5] = x->grad[5]*ret->f;
  ret->grad[6] = x->grad[6]*ret->f;
  return ret;
}

dualN *expDX(dualN *x) {
  x->f = exp(x->f);
  x->grad[0] = x->grad[0]*x->f;
  x->grad[1] = x->grad[1]*x->f;
  x->grad[2] = x->grad[2]*x->f;
  x->grad[3] = x->grad[3]*x->f;
  x->grad[4] = x->grad[4]*x->f;
  x->grad[5] = x->grad[5]*x->f;
  x->grad[6] = x->grad[6]*x->f;
  return x;
}

dualN *add2(dualN *x, dualN *y, dualN *ret) {
  ret->f = x->f + y->f;
  ret->grad[0] = x->grad[0] + y->grad[0];
  ret->grad[1] = x->grad[1] + y->grad[1];
  ret->grad[2] = x->grad[2] + y->grad[2];
  ret->grad[3] = x->grad[3] + y->grad[3];
  ret->grad[4] = x->grad[4] + y->grad[4];
  ret->grad[5] = x->grad[5] + y->grad[5];
  ret->grad[6] = x->grad[6] + y->grad[6];
  return ret;
}

dualN *add2X(dualN *x, dualN *y) {
  x->f = x->f + y->f;
  x->grad[0] = x->grad[0] + y->grad[0];
  x->grad[1] = x->grad[1] + y->grad[1];
  x->grad[2] = x->grad[2] + y->grad[2];
  x->grad[3] = x->grad[3] + y->grad[3];
  x->grad[4] = x->grad[4] + y->grad[4];
  x->grad[5] = x->grad[5] + y->grad[5];
  x->grad[6] = x->grad[6] + y->grad[6];
  return x;
}

dualN *add2X2(dualN *x, dualN *y) {
  y->f = x->f + y->f;
  y->grad[0] = x->grad[0] + y->grad[0];
  y->grad[1] = x->grad[1] + y->grad[1];
  y->grad[2] = x->grad[2] + y->grad[2];
  y->grad[3] = x->grad[3] + y->grad[3];
  y->grad[4] = x->grad[4] + y->grad[4];
  y->grad[5] = x->grad[5] + y->grad[5];
  y->grad[6] = x->grad[6] + y->grad[6];
  return y;
}

dualN *add2d(dualN *x, double y, dualN *ret) {
  ret->f = x->f + y;
  ret->grad[0] = x->grad[0];
  ret->grad[1] = x->grad[1];
  ret->grad[2] = x->grad[2];
  ret->grad[3] = x->grad[3];
  ret->grad[4] = x->grad[4];
  ret->grad[5] = x->grad[5];
  ret->grad[6] = x->grad[6];
  return ret;
}

dualN *add2dX(dualN *x, double y) {
  x->f = x->f + y;
  x->grad[0] = x->grad[0];
  x->grad[1] = x->grad[1];
  x->grad[2] = x->grad[2];
  x->grad[3] = x->grad[3];
  x->grad[4] = x->grad[4];
  x->grad[5] = x->grad[5];
  x->grad[6] = x->grad[6];
  return x;
}

#define addd2(x, y) add2d(y, x, ret)
#define addd2X(x, y) add2dX(y, x)

dualN *subtr2(dualN *x, dualN *y, dualN *ret) {
  ret->f = x->f - y->f;
  ret->grad[0] = x->grad[0] - y->grad[0];
  ret->grad[1] = x->grad[1] - y->grad[1];
  ret->grad[2] = x->grad[2] - y->grad[2];
  ret->grad[3] = x->grad[3] - y->grad[3];
  ret->grad[4] = x->grad[4] - y->grad[4];
  ret->grad[5] = x->grad[5] - y->grad[5];
  ret->grad[6] = x->grad[6] - y->grad[6];
  return ret;
}

dualN *subtr2X(dualN *x, dualN *y) {
  x->f = x->f - y->f;
  x->grad[0] = x->grad[0] - y->grad[0];
  x->grad[1] = x->grad[1] - y->grad[1];
  x->grad[2] = x->grad[2] - y->grad[2];
  x->grad[3] = x->grad[3] - y->grad[3];
  x->grad[4] = x->grad[4] - y->grad[4];
  x->grad[5] = x->grad[5] - y->grad[5];
  x->grad[6] = x->grad[6] - y->grad[6];
  return x;
}

dualN *subtr2X2(dualN *x, dualN *y) {
  y->f = x->f - y->f;
  y->grad[0] = x->grad[0] - y->grad[0];
  y->grad[1] = x->grad[1] - y->grad[1];
  y->grad[2] = x->grad[2] - y->grad[2];
  y->grad[3] = x->grad[3] - y->grad[3];
  y->grad[4] = x->grad[4] - y->grad[4];
  y->grad[5] = x->grad[5] - y->grad[5];
  y->grad[6] = x->grad[6] - y->grad[6];
  return y;
}

dualN *negD(dualN *x, dualN *ret) {
  ret->f = -x->f;
  ret->grad[0] = -x->grad[0];
  ret->grad[1] = -x->grad[1];
  ret->grad[2] = -x->grad[2];
  ret->grad[3] = -x->grad[3];
  ret->grad[4] = -x->grad[4];
  ret->grad[5] = -x->grad[5];
  ret->grad[6] = -x->grad[6];
  return ret;
}

dualN *negDX(dualN *x) {
  x->f = -x->f;
  x->grad[0] = -x->grad[0];
  x->grad[1] = -x->grad[1];
  x->grad[2] = -x->grad[2];
  x->grad[3] = -x->grad[3];
  x->grad[4] = -x->grad[4];
  x->grad[5] = -x->grad[5];
  x->grad[6] = -x->grad[6];
  return x;
}

dualN *subtrd2(double x, dualN *y,   dualN *ret) {
  ret->f = x - y->f;
  ret->grad[0] = -y->grad[0];
  ret->grad[1] = -y->grad[1];
  ret->grad[2] = -y->grad[2];
  ret->grad[3] = -y->grad[3];
  ret->grad[4] = -y->grad[4];
  ret->grad[5] = -y->grad[5];
  ret->grad[6] = -y->grad[6];
  return ret;
}

dualN *subtrd2X(double x, dualN *y) {
  y->f = x - y->f;
  y->grad[0] = -y->grad[0];
  y->grad[1] = -y->grad[1];
  y->grad[2] = -y->grad[2];
  y->grad[3] = -y->grad[3];
  y->grad[4] = -y->grad[4];
  y->grad[5] = -y->grad[5];
  y->grad[6] = -y->grad[6];
  return y;
}

dualN *subtr2d(dualN *x, double y, dualN *ret) {
  ret->f = x->f - y;
  ret->grad[0] = x->grad[0];
  ret->grad[1] = x->grad[1];
  ret->grad[2] = x->grad[2];
  ret->grad[3] = x->grad[3];
  ret->grad[4] = x->grad[4];
  ret->grad[5] = x->grad[5];
  ret->grad[6] = x->grad[6];
  return ret;
}

dualN *subtr2dX(dualN *x, double y) {
  x->f = x->f - y;
  x->grad[0] = x->grad[0];
  x->grad[1] = x->grad[1];
  x->grad[2] = x->grad[2];
  x->grad[3] = x->grad[3];
  x->grad[4] = x->grad[4];
  x->grad[5] = x->grad[5];
  x->grad[6] = x->grad[6];
  return x;
}

dualN *prod2(dualN *e1, dualN *e2, dualN *ret) {
  ret->f = e1->f * e2->f;
  ret->grad[0] = e1->grad[0] * e2->f + e1->f * e2->grad[0];
  ret->grad[1] = e1->grad[1] * e2->f + e1->f * e2->grad[1];
  ret->grad[2] = e1->grad[2] * e2->f + e1->f * e2->grad[2];
  ret->grad[3] = e1->grad[3] * e2->f + e1->f * e2->grad[3];
  ret->grad[4] = e1->grad[4] * e2->f + e1->f * e2->grad[4];
  ret->grad[5] = e1->grad[5] * e2->f + e1->f * e2->grad[5];
  ret->grad[6] = e1->grad[6] * e2->f + e1->f * e2->grad[6];
  return ret;
}

dualN *prod2X(dualN *e1, dualN *e2) {
  e1->f = e1->f * e2->f;
  e1->grad[0] = e1->grad[0] * e2->f + e1->f * e2->grad[0];
  e1->grad[1] = e1->grad[1] * e2->f + e1->f * e2->grad[1];
  e1->grad[2] = e1->grad[2] * e2->f + e1->f * e2->grad[2];
  e1->grad[3] = e1->grad[3] * e2->f + e1->f * e2->grad[3];
  e1->grad[4] = e1->grad[4] * e2->f + e1->f * e2->grad[4];
  e1->grad[5] = e1->grad[5] * e2->f + e1->f * e2->grad[5];
  e1->grad[6] = e1->grad[6] * e2->f + e1->f * e2->grad[6];
  return e1;
}

dualN *prod2X2(dualN *e1, dualN *e2) {
  e2->f = e1->f * e2->f;
  e2->grad[0] = e1->grad[0] * e2->f + e1->f * e2->grad[0];
  e2->grad[1] = e1->grad[1] * e2->f + e1->f * e2->grad[1];
  e2->grad[2] = e1->grad[2] * e2->f + e1->f * e2->grad[2];
  e2->grad[3] = e1->grad[3] * e2->f + e1->f * e2->grad[3];
  e2->grad[4] = e1->grad[4] * e2->f + e1->f * e2->grad[4];
  e2->grad[5] = e1->grad[5] * e2->f + e1->f * e2->grad[5];
  e2->grad[6] = e1->grad[6] * e2->f + e1->f * e2->grad[6];
  return e2;
}

dualN *prodd2(double e1, dualN *e2, dualN *ret) {
  ret->f = e2->f * e1;
  ret->grad[0] = e2->grad[0]*e1;
  ret->grad[1] = e2->grad[1]*e1;
  ret->grad[2] = e2->grad[2]*e1;
  ret->grad[3] = e2->grad[3]*e1;
  ret->grad[4] = e2->grad[4]*e1;
  ret->grad[5] = e2->grad[5]*e1;
  ret->grad[6] = e2->grad[6]*e1;
  return ret;
}

dualN *prodd2X(double e1, dualN *e2) {
  e2->f = e2->f * e1;
  e2->grad[0] = e2->grad[0]*e1;
  e2->grad[1] = e2->grad[1]*e1;
  e2->grad[2] = e2->grad[2]*e1;
  e2->grad[3] = e2->grad[3]*e1;
  e2->grad[4] = e2->grad[4]*e1;
  e2->grad[5] = e2->grad[5]*e1;
  e2->grad[6] = e2->grad[6]*e1;
  return e2;
}

#define prod2d(x, y, dn) prodd2(y, x, dn)
#define prod2dX(x, y) prodd2X(y, x)

dualN *div2(dualN *e1, dualN *e2, dualN *ret) {
  double invE2 = 1.0 / e2->f;
  ret->f = e1->f*invE2;
  ret->grad[0] = (e1->grad[0] - ret->f * e2->grad[0]) * invE2;
  ret->grad[1] = (e1->grad[1] - ret->f * e2->grad[1]) * invE2;
  ret->grad[2] = (e1->grad[2] - ret->f * e2->grad[2]) * invE2;
  ret->grad[3] = (e1->grad[3] - ret->f * e2->grad[3]) * invE2;
  ret->grad[4] = (e1->grad[4] - ret->f * e2->grad[4]) * invE2;
  ret->grad[5] = (e1->grad[5] - ret->f * e2->grad[5]) * invE2;
  ret->grad[6] = (e1->grad[6] - ret->f * e2->grad[6]) * invE2;
  return ret;
}

dualN *div2X(dualN *e1, dualN *e2) {
  double invE2 = 1.0 / e2->f;
  e1->f = e1->f*invE2;
  e1->grad[0] = (e1->grad[0] - e1->f * e2->grad[0]) * invE2;
  e1->grad[1] = (e1->grad[1] - e1->f * e2->grad[1]) * invE2;
  e1->grad[2] = (e1->grad[2] - e1->f * e2->grad[2]) * invE2;
  e1->grad[3] = (e1->grad[3] - e1->f * e2->grad[3]) * invE2;
  e1->grad[4] = (e1->grad[4] - e1->f * e2->grad[4]) * invE2;
  e1->grad[5] = (e1->grad[5] - e1->f * e2->grad[5]) * invE2;
  e1->grad[6] = (e1->grad[6] - e1->f * e2->grad[6]) * invE2;
  return e1;
}

dualN *div2X2(dualN *e1, dualN *e2) {
  double invE2 = 1.0 / e2->f;
  e2->f = e1->f*invE2;
  e2->grad[0] = (e1->grad[0] - e2->f * e2->grad[0]) * invE2;
  e2->grad[1] = (e1->grad[1] - e2->f * e2->grad[1]) * invE2;
  e2->grad[2] = (e1->grad[2] - e2->f * e2->grad[2]) * invE2;
  e2->grad[3] = (e1->grad[3] - e2->f * e2->grad[3]) * invE2;
  e2->grad[4] = (e1->grad[4] - e2->f * e2->grad[4]) * invE2;
  e2->grad[5] = (e1->grad[5] - e2->f * e2->grad[5]) * invE2;
  e2->grad[6] = (e1->grad[6] - e2->f * e2->grad[6]) * invE2;
  return e2;
}

dualN *divd2(double e1, dualN *e2, dualN *ret) {
  double invE2 = 1.0 / e2->f;
  ret->f = e1*invE2;
  ret->grad[0] = ret->f * invE2 * e2->grad[0];
  ret->grad[1] = ret->f * invE2 * e2->grad[1];
  ret->grad[2] = ret->f * invE2 * e2->grad[2];
  ret->grad[3] = ret->f * invE2 * e2->grad[3];
  ret->grad[4] = ret->f * invE2 * e2->grad[4];
  ret->grad[5] = ret->f * invE2 * e2->grad[5];
  ret->grad[6] = ret->f * invE2 * e2->grad[6];
  return ret;
}

dualN *divd2X(double e1, dualN *e2) {
  double invE2 = 1.0 / e2->f;
  e2->f = e1*invE2;
  e2->grad[0] = e2->f * invE2 * e2->grad[0];
  e2->grad[1] = e2->f * invE2 * e2->grad[1];
  e2->grad[2] = e2->f * invE2 * e2->grad[2];
  e2->grad[3] = e2->f * invE2 * e2->grad[3];
  e2->grad[4] = e2->f * invE2 * e2->grad[4];
  e2->grad[5] = e2->f * invE2 * e2->grad[5];
  e2->grad[6] = e2->f * invE2 * e2->grad[6];
  return e2;
}

dualN *div2d(dualN *e1, double e2, dualN *ret) {
  double invE2 = 1.0 / e2;
  ret->f = e1->f*invE2;
  ret->grad[0] = invE2 * e1->grad[0];
  ret->grad[1] = invE2 * e1->grad[1];
  ret->grad[2] = invE2 * e1->grad[2];
  ret->grad[3] = invE2 * e1->grad[3];
  ret->grad[4] = invE2 * e1->grad[4];
  ret->grad[5] = invE2 * e1->grad[5];
  ret->grad[5] = invE2 * e1->grad[5];
  ret->grad[6] = invE2 * e1->grad[6];
  return ret;
}

dualN *div2dX(dualN *e1, double e2) {
  double invE2 = 1.0 / e2;
  e1->f = e1->f*invE2;
  e1->grad[0] = invE2 * e1->grad[0];
  e1->grad[1] = invE2 * e1->grad[1];
  e1->grad[2] = invE2 * e1->grad[2];
  e1->grad[3] = invE2 * e1->grad[3];
  e1->grad[4] = invE2 * e1->grad[4];
  e1->grad[5] = invE2 * e1->grad[5];
  e1->grad[5] = invE2 * e1->grad[5];
  e1->grad[6] = invE2 * e1->grad[6];
  return e1;
}


dualN *atan2D(dualN *y, dualN *x, dualN *ret) {
  double derg = y->f * y->f;
  double derf = x->f * x->f;
  derg = 1.0 / (derg + derf);
  derf = x->f * derg;
  derg = y->f * derg;
  ret->f = atan2(y->f, x->f);
  ret->grad[0] = y->grad[0] * derf - derg * x->grad[0];
  ret->grad[1] = y->grad[1] * derf - derg * x->grad[1];
  ret->grad[2] = y->grad[2] * derf - derg * x->grad[2];
  ret->grad[3] = y->grad[3] * derf - derg * x->grad[3];
  ret->grad[4] = y->grad[4] * derf - derg * x->grad[4];
  ret->grad[5] = y->grad[5] * derf - derg * x->grad[5];
  ret->grad[6] = y->grad[6] * derf - derg * x->grad[6];
  return ret;
}

dualN *atan2DX(dualN *y, dualN *x) {
  double derg = y->f * y->f;
  double derf = x->f * x->f;
  derg = 1.0 / (derg + derf);
  derf = x->f * derg;
  derg = y->f * derg;
  y->f = atan2(y->f, x->f);
  y->grad[0] = y->grad[0] * derf - derg * x->grad[0];
  y->grad[1] = y->grad[1] * derf - derg * x->grad[1];
  y->grad[2] = y->grad[2] * derf - derg * x->grad[2];
  y->grad[3] = y->grad[3] * derf - derg * x->grad[3];
  y->grad[4] = y->grad[4] * derf - derg * x->grad[4];
  y->grad[5] = y->grad[5] * derf - derg * x->grad[5];
  y->grad[6] = y->grad[6] * derf - derg * x->grad[6];
  return y;
}

dualN *atan2DX2(dualN *y, dualN *x) {
  double derg = y->f * y->f;
  double derf = x->f * x->f;
  derg = 1.0 / (derg + derf);
  derf = x->f * derg;
  derg = y->f * derg;
  x->f = atan2(y->f, x->f);
  x->grad[0] = y->grad[0] * derf - derg * x->grad[0];
  x->grad[1] = y->grad[1] * derf - derg * x->grad[1];
  x->grad[2] = y->grad[2] * derf - derg * x->grad[2];
  x->grad[3] = y->grad[3] * derf - derg * x->grad[3];
  x->grad[4] = y->grad[4] * derf - derg * x->grad[4];
  x->grad[5] = y->grad[5] * derf - derg * x->grad[5];
  x->grad[6] = y->grad[6] * derf - derg * x->grad[6];
  return x;
}


dualN *cosD(dualN *x, dualN *ret){
  //dual(f = cos(x@f), grad = -x@grad * sin(x@f))
  ret->f = cos(x->f);
  double sf = sin(x->f);
  ret->grad[0] = x->grad[0] * sf;
  ret->grad[1] = x->grad[1] * sf;
  ret->grad[2] = x->grad[2] * sf;
  ret->grad[3] = x->grad[3] * sf;
  ret->grad[4] = x->grad[4] * sf;
  ret->grad[5] = x->grad[5] * sf;
  ret->grad[6] = x->grad[6] * sf;
  return ret;
}

dualN *cosDX(dualN *x){
  double sf = sin(x->f);
  x->f = cos(x->f);
  x->grad[0] *= sf;
  x->grad[1] *= sf;
  x->grad[2] *= sf;
  x->grad[3] *= sf;
  x->grad[4] *= sf;
  x->grad[5] *= sf;
  x->grad[6] *= sf;
  return x;
}

dualN *sinD(dualN *x, dualN *ret) {
  ret->f = sin(x->f);
  double cs = cos(x->f);
  ret->grad[0] = x->grad[0] * cs;
  ret->grad[1] = x->grad[1] * cs;
  ret->grad[2] = x->grad[2] * cs;
  ret->grad[3] = x->grad[3] * cs;
  ret->grad[4] = x->grad[4] * cs;
  ret->grad[5] = x->grad[5] * cs;
  ret->grad[6] = x->grad[6] * cs;
  return ret;
}

dualN *sinDX(dualN *x) {
  double cs = cos(x->f);
  x->f = sin(x->f);
  x->grad[0] *= cs;
  x->grad[1] *= cs;
  x->grad[2] *= cs;
  x->grad[3] *= cs;
  x->grad[4] *= cs;
  x->grad[5] *= cs;
  x->grad[6] *= cs;
  return x;
}

dualN *pow2d(dualN *base, double e, dualN *ret) {
  ret->f = pow(base->f, e);
  double mult = e*pow(base->f, e - 1.0);
  ret->grad[0] = base->grad[0] * mult;
  ret->grad[1] = base->grad[1] * mult;
  ret->grad[2] = base->grad[2] * mult;
  ret->grad[3] = base->grad[3] * mult;
  ret->grad[4] = base->grad[4] * mult;
  ret->grad[5] = base->grad[5] * mult;
  ret->grad[6] = base->grad[6] * mult;
  return ret;
}

dualN *pow2dX(dualN *base, double e) {
  double f = pow(base->f, e);
  double mult=e*pow(base->f, e - 1.0);
  base->f = f;
  base->grad[0] *= mult;
  base->grad[1] *= mult;
  base->grad[2] *= mult;
  base->grad[3] *= mult;
  base->grad[4] *= mult;
  base->grad[5] *= mult;
  base->grad[6] *= mult;
  return base;
}

#endif
