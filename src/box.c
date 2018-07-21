#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double powerD(double x, double lambda, int yj){
  if (lambda == 1.0) return x;
  if (yj == 0){
    if (x < 0) return NA_REAL;
    if (lambda == 0.0) return log(x);
    return (pow(x, lambda) - 1.0)/lambda;
  } else {
    if (x >= 0){
      if (lambda == 0.0) return log(x + 1.0);
      return (pow(x + 1.0, lambda) - 1.0)/lambda;
    } else {
      if (lambda == 2.0) return -log(1.0 - x);
      double l2 = 2.0 - lambda;
      return (1.0 - pow(1.0 - x, l2))/l2;
    }
  }
}

double powerDD(double x, double lambda, int yj){
  if (lambda == 1) return 1;
  if (yj == 0){
    if (x < 0) return NA_REAL;
    if (lambda == 0.0) return 1/x;
    // pow(x,lambda)/lambda - 1/lambda
    return pow(x, lambda-1);
  } else {
    if (x >= 0){
      if (lambda == 0.0) return 1/(x + 1.0);
      return pow(x + 1.0, lambda-1.0);
    } else {
      if (lambda == 2.0) return -1/(1.0 - x);
      return pow(1.0 - x, 1.0-lambda);
    }
  }
}

double powerL(double x, double lambda, int yj){
  // logLik addition based on dTBS
  // yj is indicator for yeo-johson
  if (lambda == 1.0) return 0;
  if (yj == 0){
    if (x > 0) return (lambda - 1.0)*log(x);
    return NA_REAL;
  } else {
    if (x >= 0) return (lambda - 1.0)*log(x+1.0);
    return (1.0-lambda)*log(1.0-x);
  }
  // d = 0.0 for cox box
  // d = 1.0 fo  Yeo- Johnson
  // logLik approximation
  // y^(lambda)/lambda - 1/lambda
  // dh/dy = y^(lambda-1)
  // log(dh/dy) = (lambda-1)*log(y) + log(lambda) 
  //
  // (x + 1.0)^(lambda)/lambda - 1/lambda
  // dh/dy = (x+1.0)^(lambda-1)
  // log(dh/dy) = (lambda-1)*log(x+1.0)
  
  // For negative values yj becomes
  // (-x+1)^(2-lambda)/(2-lambda) - 1/(2-lambda)
  // dh/dy = (-x+1)^(1-lambda)
  // log(dh/dy) = (1-lambda)*log(-x+1)
}

double powerDL(double x, double lambda, int yj){
  // d(logLik/dlambda)
  if (!yj){
    if (x > 0) return log(x);
    return NA_REAL;
  } else {
    if (x >= 0) return log(x+1.0);
    return -log(1.0-x);
  }
  // d = 0.0 for cox box
  // d = 1.0 fo  Yeo- Johnson
  // logLik approximation
  // y^(lambda)/lambda - 1/lambda
  // dh/dy = y^(lambda-1)
  // log(dh/dy) = (lambda-1)*log(y) + log(lambda) 
  //
  // (x + 1.0)^(lambda)/lambda - 1/lambda
  // dh/dy = (x+1.0)^(lambda-1)
  // log(dh/dy) = (lambda-1)*log(x+1.0)
  
  // For negative values yj becomes
  // (-x+1)^(2-lambda)/(2-lambda) - 1/(2-lambda)
  // dh/dy = (-x+1)^(1-lambda)
  // log(dh/dy) = (1-lambda)*log(-x+1)

}
