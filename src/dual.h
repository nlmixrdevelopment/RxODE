#ifndef dual_header
#define dual_header

typedef struct dualN {
  double f;
  int n;
  double grad[4];
} dualN;

dualN sqrtD(dualN x) {
  dualN ret;
  ret.n = x.n;
  ret.f = sqrt(x.f);
  double gr = 0.5 / ret.f;
  for (int i = x.n; i--;) {
    ret.grad[i] = x.grad[i]*gr;
  }
  return ret;
}

dualN expD(dualN x) {
  dualN ret;
  ret.n = x.n;
  ret.f = exp(x.f);
  for (int i = x.n; i--;) {
    ret.grad[i] = x.grad[i]*ret.f;
  }
  return ret;
}


dualN add2(dualN x, dualN y) {
  dualN ret;
  ret.n = x.n;
  ret.f = x.f + y.f;
  for (int i = ret.n; i--;){
    ret.grad[i] = x.grad[i] + y.grad[i];
  }
  return ret;
}

dualN add2d(dualN x, double y) {
  dualN ret;
  ret.n = x.n;
  ret.f = x.f + y;
  for (int i = ret.n; i--;){
    ret.grad[i] = x.grad[i];
  }
  return ret;
}

dualN addd2(double y, dualN x) {
  dualN ret;
  ret.n = x.n;
  ret.f = x.f + y;
  for (int i = ret.n; i--;){
    ret.grad[i] = x.grad[i];
  }
  return ret;
}

dualN subtr2(dualN x, dualN y) {
  dualN ret;
  ret.n = x.n;
  ret.f = x.f - y.f;
  for (int i = ret.n; i--;){
    ret.grad[i] = x.grad[i] - y.grad[i];
  }
  return ret;
}

dualN negD(dualN x) {
  dualN ret;
  ret.n = x.n;
  ret.f = -x.f;
  for (int i = ret.n; i--;){
    ret.grad[i] = -x.grad[i];
  }
  return ret;
}

dualN subtrd2(double x, dualN y) {
  dualN ret;
  ret.n = y.n;
  ret.f = x - y.f;
  for (int i = ret.n; i--;){
    ret.grad[i] = -y.grad[i];
  }
  return ret;
}

dualN subtr2d(dualN x, double y) {
  dualN ret;
  ret.n = x.n;
  ret.f = x.f - y;
  for (int i = ret.n; i--;){
    ret.grad[i] = x.grad[i];
  }
  return ret;
}

dualN prod2(dualN e1, dualN e2) {
  dualN ret;
  ret.n = e1.n;
  ret.f = e1.f * e2.f;
  for (int i = ret.n; i--;) {
    ret.grad[i] = e1.grad[i] * e2.f + e1.f * e2.grad[i];
  }
  return ret;
}

dualN prodd2(double e1, dualN e2) {
  dualN ret;
  ret.n = e2.n;
  ret.f = e2.f * e1;
  for (int i = ret.n; i--;) {
    ret.grad[i] = e2.grad[i]*e1;
  }
  return ret;
}

dualN div2(dualN e1, dualN e2) {
  dualN ret;
  ret.n = e2.n;
  double invE2 = 1.0 / e2.f;
  ret.f = e1.f*invE2;
  for (int i = ret.n; i--;) {
    ret.grad[i] = (e1.grad[i] - ret.f * e2.grad[i]) * invE2;
  }
  return ret;
}

dualN divd2(double e1, dualN e2) {
  dualN ret;
  ret.n = e2.n;
  double invE2 = 1.0 / e2.f;
  ret.f = e1*invE2;
  for (int i = ret.n; i--;) {
    ret.grad[i] = ret.f * invE2 * ret.grad[i];
  }
  return ret;
}

dualN div2d(dualN e1, double e2) {
  dualN ret;
  ret.n = e1.n;
  double invE2 = 1.0 / e2;
  ret.f = e1.f*invE2;
  for (int i = ret.n; i--;) {
    ret.grad[i] = invE2 * ret.grad[i];
  }
  return ret;
}

dualN iniD(double val, int which, int n){
  dualN ret;
  ret.n = n;
  ret.f = val;
  for (int i = ret.n; i--;) {
    ret.grad[i] = 0.0;
  }
  if (which >= 0) {
    ret.grad[which] = 1.0;
  }
  return ret;
}

#endif
