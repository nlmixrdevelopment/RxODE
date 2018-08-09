#include <stdio.h>
#include <stdarg.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#define X0MIN 1e-5
#define eps sqrt(DOUBLE_EPS)

double powerDi(double x, double lambda, int yj){
  double x0=x, ret, l2;
  switch(yj){
  case 3:
    return exp(x);
  case 2: 
    return x;
  case 0:
    if (lambda == 1.0) return (x+1.0);
    if (lambda == 0) return exp(x);
    // (x^lambda-1)/lambda=y
    // (lambda*y+1)^(1/lambda)
    x0 = x*lambda+1.0;
    if (x0 <= eps) return eps;
    ret = pow(x0, 1.0/lambda);
    if (ISNA(ret)) {
      // Warning?
      return eps;
    }
    return ret;
  case 1:
    if (lambda == 1.0) return x;
    if (x >= 0){
      // log(x+1)= y; exp(y)-1=x
      if (lambda == 0) return expm1(x);
      // ((x+1)^lambda-1)/lambda=y
      // (y*lambda+1)^(1/y)-1=y
      return pow(x*lambda+1.0, 1.0/lambda)-1.0;
    } else {
      // (-(1-x)^(2-lambda)-1)/(2-lambda)
      if (lambda ==  2.0) return -expm1(-x);
      // (-(1-x)^(2-lambda)-1)/(2-lambda) = y
      l2 = (2.0 - lambda);
      return 1.0 - pow(1.0 - l2*x, 1.0/l2);
    }
  }
  return NA_REAL;
}

double powerD(double x, double lambda, int yj){
  double x0=x, l2;
  switch (yj){
  case 3:
    if (x <= eps) x0= eps;
    return log(x0);
  case 2:
    return x;
  case 0:
    if (lambda == 1.0) return x-1.0;
    if (x <= eps) x0= eps;
    if (lambda ==  0.0) return log(x0);
    return (pow(x0, lambda) - 1.0)/lambda;
  case 1:
    if (lambda == 1.0) return x;
    if (x >= 0){
      if (lambda == 0) return log1p(x);
      return (pow(x + 1.0, lambda) - 1.0)/lambda;
    } else {
      if (lambda == 2.0) return -log1p(-x);
      l2 = 2.0 - lambda;
      return (1.0 - pow(1.0 - x, l2))/l2;
    }
  }
  return NA_REAL;
}

double powerDD(double x, double lambda, int yj){
  double x0 = x;
  switch(yj){
  case 3:
    if (x <= eps) return x0 = eps;
    return 1/x0;
  case 2:
    return 1.0;
  case 0:
    if (lambda == 1.0) return 1.0;
    if (x <= eps) return x0 = eps;
    if (lambda == 0.0) return 1/x0;
    // pow(x,lambda)/lambda - 1/lambda
    return pow(x0, lambda-1);
  case 1:
    if (lambda ==  1.0) return 1.0;
    if (x >= 0){
      if (lambda == 0.0) return 1.0/(x + 1.0);
      return pow(x + 1.0, lambda-1.0);
    } else {
      if (lambda == 2.0) return -1/(1.0 - x);
      return pow(1.0 - x, 1.0-lambda);
    }
  }
  return NA_REAL;
}

double powerDDD(double x, double lambda, int yj){
  double x0 = x;
  switch(yj){
  case 3:
    if (x <= eps) x0 = eps;
    return -1/(x0*x0);
  case 2: 
    return 0;
  case 0:
    if (lambda == 1.0) return 0;
    if (x <= eps) return x0 = eps;
    if (lambda == 0.0) return -1/(x0*x0);
    // pow(x,lambda)/lambda - 1/lambda
    return (lambda-1)*pow(x0, lambda-2);
  case 1:
    if (lambda == 1.0) return 0;
    if (x >= 0){
      if (lambda ==  0.0) return -1/((x + 1.0)*(x + 1.0));
      return (lambda-1.0)*pow(x + 1.0, lambda-2.0);
    } else {
      if (lambda == 2.0) return -1/((1.0 - x)*(1.0 - x));
      return -(1.0-lambda)*pow(1.0 - x, -lambda);
    }
  }
  return NA_REAL;
}

double powerL(double x, double lambda, int yj){
  double x0 = x;
  switch(yj){
  case 3:
    if (x <= eps) x0 = eps;
    return -log(x0);
  case 2:
    return 0;
  case 0:
    if (lambda == 1.0) return 0;
    if (x <= eps) x0 = eps;
    return (lambda - 1.0)*log(x0);
  case 1:
    if (x >= 0) return (lambda - 1.0)*log1p(x);
    return (1.0-lambda)*log1p(-x);
  }
  return NA_REAL;
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
  double x0 = x;
  switch (yj){
  case 3:
    if (x <= eps) x0 = eps;
    return log(x0);
  case 2:
    return 0;
  case 0:
    if (lambda == 1.0) return 0;
    if (x <= eps) x0 = eps;
    return log(x0);
  case 1:
    if (lambda == 1.0) return 0;
    if (x >= 0) return log1p(x);
    return -log1p(x);
  }
  return NA_REAL;
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
