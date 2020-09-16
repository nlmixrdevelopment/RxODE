#ifndef dual_header
#define dual_header
#include <math.h>

typedef struct dualN {
  double f;
  int n;
  double grad[4];
} dualN;

dualN sqrtD(dualN x) {
  dualN ret;
  ret.f = sqrt(x.f);
  double gr = 0.5 / ret.f;
  ret.grad[0] = x.grad[0]*gr;
  ret.grad[1] = x.grad[1]*gr;
  ret.grad[2] = x.grad[2]*gr;
  ret.grad[3] = x.grad[3]*gr;
  return ret;
}

dualN expD(dualN x) {
  dualN ret;
  ret.f = exp(x.f);
  ret.grad[0] = x.grad[0]*ret.f;
  ret.grad[1] = x.grad[1]*ret.f;
  ret.grad[2] = x.grad[2]*ret.f;
  ret.grad[3] = x.grad[3]*ret.f;
  return ret;
}

dualN add2(dualN x, dualN y) {
  dualN ret;
  ret.f = x.f + y.f;
  ret.grad[0] = x.grad[0] + y.grad[0];
  ret.grad[1] = x.grad[1] + y.grad[1];
  ret.grad[2] = x.grad[2] + y.grad[2];
  ret.grad[3] = x.grad[3] + y.grad[3];
  return ret;
}

dualN add2d(dualN x, double y) {
  dualN ret;
  ret.f = x.f + y;
  ret.grad[0] = x.grad[0];
  ret.grad[1] = x.grad[1];
  ret.grad[2] = x.grad[2];
  ret.grad[3] = x.grad[3];
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
  return ret;
}

dualN negD(dualN x) {
  dualN ret;
  ret.f = -x.f;
  ret.grad[0] = -x.grad[0];
  ret.grad[1] = -x.grad[1];
  ret.grad[2] = -x.grad[2];
  ret.grad[3] = -x.grad[3];
  return ret;
}

dualN subtrd2(double x, dualN y) {
  dualN ret;
  ret.f = x - y.f;
  ret.grad[0] = -y.grad[0];
  ret.grad[1] = -y.grad[1];
  ret.grad[2] = -y.grad[2];
  ret.grad[3] = -y.grad[3];
  return ret;
}

dualN subtr2d(dualN x, double y) {
  dualN ret;
  ret.f = x.f - y;
  ret.grad[0] = x.grad[0];
  ret.grad[1] = x.grad[1];
  ret.grad[2] = x.grad[2];
  ret.grad[3] = x.grad[3];
  return ret;
}

dualN prod2(dualN e1, dualN e2) {
  dualN ret;
  ret.f = e1.f * e2.f;
  ret.grad[0] = e1.grad[0] * e2.f + e1.f * e2.grad[0];
  ret.grad[1] = e1.grad[1] * e2.f + e1.f * e2.grad[1];
  ret.grad[2] = e1.grad[2] * e2.f + e1.f * e2.grad[2];
  ret.grad[3] = e1.grad[3] * e2.f + e1.f * e2.grad[3];
  return ret;
}

dualN prodd2(double e1, dualN e2) {
  dualN ret;
  ret.f = e2.f * e1;
  ret.grad[0] = e2.grad[0]*e1;
  ret.grad[1] = e2.grad[1]*e1;
  ret.grad[2] = e2.grad[2]*e1;
  ret.grad[3] = e2.grad[3]*e1;
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
  return ret;
}

dualN iniD(double val, int which){
  dualN ret;
  ret.f = val;
  ret.grad[0] = 0.0;
  ret.grad[1] = 0.0;
  ret.grad[2] = 0.0;
  ret.grad[3] = 0.0;
  if (which >= 0) {
    ret.grad[which] = 1.0;
  }
  return ret;
}

#endif
