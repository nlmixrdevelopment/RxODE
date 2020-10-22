#include "../inst/include/RxODE.h"

double powerDi(double x, double lambda, int yj){
  return _powerDi(x, lambda, yj);
}

double powerD(double x, double lambda, int yj){
  return _powerD(x, lambda, yj);
}

double powerDD(double x, double lambda, int yj){
  return _powerDD(x, lambda, yj);
}

double powerDDD(double x, double lambda, int yj){
  return _powerDDD(x, lambda, yj);
}

double powerL(double x, double lambda, int yj){
  return _powerL(x, lambda, yj);
}

double powerDL(double x, double lambda, int yj){
  return _powerDL(x, lambda, yj);
}
