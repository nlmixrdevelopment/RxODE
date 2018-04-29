//Generated from refresh.R for 12 dimensions
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>
#include <Rmath.h>
SEXP _rxCholInv(SEXP dms, SEXP theta, SEXP tn){
int dm=INTEGER(dms)[0];
if (dm == 0){
  SEXP ret=  PROTECT(allocVector(INTSXP,1));
  INTEGER(ret)[0] = 12;
  UNPROTECT(1);
  return(ret);
}else if (dm == 1){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 1;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -3 || theta_n > 1){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 1){
    error("Requires vector with 1 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 1, 1));for (int i = 0; i < 1; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 1));for(int i = 0; i < 1; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 2){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 2;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -5 || theta_n > 3){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 3){
    error("Requires vector with 3 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 2, 2));for (int i = 0; i < 4; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[1];
      REAL(ret)[3] = exp(REAL(theta)[2]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[3] = R_pow_di(REAL(theta)[1], 2) + exp(2 * REAL(theta)[2]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[1] * exp(REAL(theta)[0]);
    }
    else if (theta_n == 2){
      REAL(ret)[1] = exp(REAL(theta)[0]);
      REAL(ret)[2] = exp(REAL(theta)[0]);
      REAL(ret)[3] = 2 * REAL(theta)[1];
    }
    else if (theta_n == 3){
      REAL(ret)[3] = 2 * exp(2 * REAL(theta)[2]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 2));for(int i = 0; i < 2; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 2 * exp(2 * REAL(theta)[2]);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 3){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 3;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -8 || theta_n > 6){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 6){
    error("Requires vector with 6 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 3, 3));for (int i = 0; i < 9; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[1];
      REAL(ret)[4] = exp(REAL(theta)[2]);
      REAL(ret)[6] = REAL(theta)[3];
      REAL(ret)[7] = REAL(theta)[4];
      REAL(ret)[8] = exp(REAL(theta)[5]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[4] = R_pow_di(REAL(theta)[1], 2) + exp(2 * REAL(theta)[2]);
      REAL(ret)[5] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[6] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[8] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + exp(2 * REAL(theta)[5]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[3] * exp(REAL(theta)[0]);
    }
    else if (theta_n == 2){
      REAL(ret)[1] = exp(REAL(theta)[0]);
      REAL(ret)[3] = exp(REAL(theta)[0]);
      REAL(ret)[4] = 2 * REAL(theta)[1];
      REAL(ret)[5] = REAL(theta)[3];
      REAL(ret)[7] = REAL(theta)[3];
    }
    else if (theta_n == 3){
      REAL(ret)[4] = 2 * exp(2 * REAL(theta)[2]);
      REAL(ret)[5] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[7] = REAL(theta)[4] * exp(REAL(theta)[2]);
    }
    else if (theta_n == 4){
      REAL(ret)[2] = exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[1];
      REAL(ret)[6] = exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[1];
      REAL(ret)[8] = 2 * REAL(theta)[3];
    }
    else if (theta_n == 5){
      REAL(ret)[5] = exp(REAL(theta)[2]);
      REAL(ret)[7] = exp(REAL(theta)[2]);
      REAL(ret)[8] = 2 * REAL(theta)[4];
    }
    else if (theta_n == 6){
      REAL(ret)[8] = 2 * exp(2 * REAL(theta)[5]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 3));for(int i = 0; i < 3; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 2 * exp(2 * REAL(theta)[2]);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 2 * exp(2 * REAL(theta)[5]);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 4){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 4;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -12 || theta_n > 10){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 10){
    error("Requires vector with 10 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 4, 4));for (int i = 0; i < 16; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[1];
      REAL(ret)[5] = exp(REAL(theta)[2]);
      REAL(ret)[8] = REAL(theta)[3];
      REAL(ret)[9] = REAL(theta)[4];
      REAL(ret)[10] = exp(REAL(theta)[5]);
      REAL(ret)[12] = REAL(theta)[6];
      REAL(ret)[13] = REAL(theta)[7];
      REAL(ret)[14] = REAL(theta)[8];
      REAL(ret)[15] = exp(REAL(theta)[9]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[5] = R_pow_di(REAL(theta)[1], 2) + exp(2 * REAL(theta)[2]);
      REAL(ret)[6] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[7] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[8] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[10] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + exp(2 * REAL(theta)[5]);
      REAL(ret)[11] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[12] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[15] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + exp(2 * REAL(theta)[9]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[12] = REAL(theta)[6] * exp(REAL(theta)[0]);
    }
    else if (theta_n == 2){
      REAL(ret)[1] = exp(REAL(theta)[0]);
      REAL(ret)[4] = exp(REAL(theta)[0]);
      REAL(ret)[5] = 2 * REAL(theta)[1];
      REAL(ret)[6] = REAL(theta)[3];
      REAL(ret)[7] = REAL(theta)[6];
      REAL(ret)[9] = REAL(theta)[3];
      REAL(ret)[13] = REAL(theta)[6];
    }
    else if (theta_n == 3){
      REAL(ret)[5] = 2 * exp(2 * REAL(theta)[2]);
      REAL(ret)[6] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[7] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[9] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[7] * exp(REAL(theta)[2]);
    }
    else if (theta_n == 4){
      REAL(ret)[2] = exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[1];
      REAL(ret)[8] = exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[1];
      REAL(ret)[10] = 2 * REAL(theta)[3];
      REAL(ret)[11] = REAL(theta)[6];
      REAL(ret)[14] = REAL(theta)[6];
    }
    else if (theta_n == 5){
      REAL(ret)[6] = exp(REAL(theta)[2]);
      REAL(ret)[9] = exp(REAL(theta)[2]);
      REAL(ret)[10] = 2 * REAL(theta)[4];
      REAL(ret)[11] = REAL(theta)[7];
      REAL(ret)[14] = REAL(theta)[7];
    }
    else if (theta_n == 6){
      REAL(ret)[10] = 2 * exp(2 * REAL(theta)[5]);
      REAL(ret)[11] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[14] = REAL(theta)[8] * exp(REAL(theta)[5]);
    }
    else if (theta_n == 7){
      REAL(ret)[3] = exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[1];
      REAL(ret)[11] = REAL(theta)[3];
      REAL(ret)[12] = exp(REAL(theta)[0]);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[14] = REAL(theta)[3];
      REAL(ret)[15] = 2 * REAL(theta)[6];
    }
    else if (theta_n == 8){
      REAL(ret)[7] = exp(REAL(theta)[2]);
      REAL(ret)[11] = REAL(theta)[4];
      REAL(ret)[13] = exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[4];
      REAL(ret)[15] = 2 * REAL(theta)[7];
    }
    else if (theta_n == 9){
      REAL(ret)[11] = exp(REAL(theta)[5]);
      REAL(ret)[14] = exp(REAL(theta)[5]);
      REAL(ret)[15] = 2 * REAL(theta)[8];
    }
    else if (theta_n == 10){
      REAL(ret)[15] = 2 * exp(2 * REAL(theta)[9]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 4));for(int i = 0; i < 4; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 2 * exp(2 * REAL(theta)[2]);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 2 * exp(2 * REAL(theta)[5]);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 2 * exp(2 * REAL(theta)[9]);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 5){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 5;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -17 || theta_n > 15){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 15){
    error("Requires vector with 15 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 5, 5));for (int i = 0; i < 25; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[1];
      REAL(ret)[6] = exp(REAL(theta)[2]);
      REAL(ret)[10] = REAL(theta)[3];
      REAL(ret)[11] = REAL(theta)[4];
      REAL(ret)[12] = exp(REAL(theta)[5]);
      REAL(ret)[15] = REAL(theta)[6];
      REAL(ret)[16] = REAL(theta)[7];
      REAL(ret)[17] = REAL(theta)[8];
      REAL(ret)[18] = exp(REAL(theta)[9]);
      REAL(ret)[20] = REAL(theta)[10];
      REAL(ret)[21] = REAL(theta)[11];
      REAL(ret)[22] = REAL(theta)[12];
      REAL(ret)[23] = REAL(theta)[13];
      REAL(ret)[24] = exp(REAL(theta)[14]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[6] = R_pow_di(REAL(theta)[1], 2) + exp(2 * REAL(theta)[2]);
      REAL(ret)[7] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[8] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[9] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[10] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[11] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[12] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + exp(2 * REAL(theta)[5]);
      REAL(ret)[13] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[14] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[15] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[16] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[18] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + exp(2 * REAL(theta)[9]);
      REAL(ret)[19] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[20] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[21] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[22] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[23] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[24] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + exp(2 * REAL(theta)[14]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[10] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[15] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[20] = REAL(theta)[10] * exp(REAL(theta)[0]);
    }
    else if (theta_n == 2){
      REAL(ret)[1] = exp(REAL(theta)[0]);
      REAL(ret)[5] = exp(REAL(theta)[0]);
      REAL(ret)[6] = 2 * REAL(theta)[1];
      REAL(ret)[7] = REAL(theta)[3];
      REAL(ret)[8] = REAL(theta)[6];
      REAL(ret)[9] = REAL(theta)[10];
      REAL(ret)[11] = REAL(theta)[3];
      REAL(ret)[16] = REAL(theta)[6];
      REAL(ret)[21] = REAL(theta)[10];
    }
    else if (theta_n == 3){
      REAL(ret)[6] = 2 * exp(2 * REAL(theta)[2]);
      REAL(ret)[7] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[8] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[9] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[11] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[21] = REAL(theta)[11] * exp(REAL(theta)[2]);
    }
    else if (theta_n == 4){
      REAL(ret)[2] = exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[1];
      REAL(ret)[10] = exp(REAL(theta)[0]);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[12] = 2 * REAL(theta)[3];
      REAL(ret)[13] = REAL(theta)[6];
      REAL(ret)[14] = REAL(theta)[10];
      REAL(ret)[17] = REAL(theta)[6];
      REAL(ret)[22] = REAL(theta)[10];
    }
    else if (theta_n == 5){
      REAL(ret)[7] = exp(REAL(theta)[2]);
      REAL(ret)[11] = exp(REAL(theta)[2]);
      REAL(ret)[12] = 2 * REAL(theta)[4];
      REAL(ret)[13] = REAL(theta)[7];
      REAL(ret)[14] = REAL(theta)[11];
      REAL(ret)[17] = REAL(theta)[7];
      REAL(ret)[22] = REAL(theta)[11];
    }
    else if (theta_n == 6){
      REAL(ret)[12] = 2 * exp(2 * REAL(theta)[5]);
      REAL(ret)[13] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[14] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[17] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[22] = REAL(theta)[12] * exp(REAL(theta)[5]);
    }
    else if (theta_n == 7){
      REAL(ret)[3] = exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[1];
      REAL(ret)[13] = REAL(theta)[3];
      REAL(ret)[15] = exp(REAL(theta)[0]);
      REAL(ret)[16] = REAL(theta)[1];
      REAL(ret)[17] = REAL(theta)[3];
      REAL(ret)[18] = 2 * REAL(theta)[6];
      REAL(ret)[19] = REAL(theta)[10];
      REAL(ret)[23] = REAL(theta)[10];
    }
    else if (theta_n == 8){
      REAL(ret)[8] = exp(REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[4];
      REAL(ret)[16] = exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[4];
      REAL(ret)[18] = 2 * REAL(theta)[7];
      REAL(ret)[19] = REAL(theta)[11];
      REAL(ret)[23] = REAL(theta)[11];
    }
    else if (theta_n == 9){
      REAL(ret)[13] = exp(REAL(theta)[5]);
      REAL(ret)[17] = exp(REAL(theta)[5]);
      REAL(ret)[18] = 2 * REAL(theta)[8];
      REAL(ret)[19] = REAL(theta)[12];
      REAL(ret)[23] = REAL(theta)[12];
    }
    else if (theta_n == 10){
      REAL(ret)[18] = 2 * exp(2 * REAL(theta)[9]);
      REAL(ret)[19] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[23] = REAL(theta)[13] * exp(REAL(theta)[9]);
    }
    else if (theta_n == 11){
      REAL(ret)[4] = exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[1];
      REAL(ret)[14] = REAL(theta)[3];
      REAL(ret)[19] = REAL(theta)[6];
      REAL(ret)[20] = exp(REAL(theta)[0]);
      REAL(ret)[21] = REAL(theta)[1];
      REAL(ret)[22] = REAL(theta)[3];
      REAL(ret)[23] = REAL(theta)[6];
      REAL(ret)[24] = 2 * REAL(theta)[10];
    }
    else if (theta_n == 12){
      REAL(ret)[9] = exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[4];
      REAL(ret)[19] = REAL(theta)[7];
      REAL(ret)[21] = exp(REAL(theta)[2]);
      REAL(ret)[22] = REAL(theta)[4];
      REAL(ret)[23] = REAL(theta)[7];
      REAL(ret)[24] = 2 * REAL(theta)[11];
    }
    else if (theta_n == 13){
      REAL(ret)[14] = exp(REAL(theta)[5]);
      REAL(ret)[19] = REAL(theta)[8];
      REAL(ret)[22] = exp(REAL(theta)[5]);
      REAL(ret)[23] = REAL(theta)[8];
      REAL(ret)[24] = 2 * REAL(theta)[12];
    }
    else if (theta_n == 14){
      REAL(ret)[19] = exp(REAL(theta)[9]);
      REAL(ret)[23] = exp(REAL(theta)[9]);
      REAL(ret)[24] = 2 * REAL(theta)[13];
    }
    else if (theta_n == 15){
      REAL(ret)[24] = 2 * exp(2 * REAL(theta)[14]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 5));for(int i = 0; i < 5; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 2 * exp(2 * REAL(theta)[2]);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 2 * exp(2 * REAL(theta)[5]);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 2 * exp(2 * REAL(theta)[9]);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 2 * exp(2 * REAL(theta)[14]);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 6){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 6;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -23 || theta_n > 21){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 21){
    error("Requires vector with 21 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 6, 6));for (int i = 0; i < 36; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[1];
      REAL(ret)[7] = exp(REAL(theta)[2]);
      REAL(ret)[12] = REAL(theta)[3];
      REAL(ret)[13] = REAL(theta)[4];
      REAL(ret)[14] = exp(REAL(theta)[5]);
      REAL(ret)[18] = REAL(theta)[6];
      REAL(ret)[19] = REAL(theta)[7];
      REAL(ret)[20] = REAL(theta)[8];
      REAL(ret)[21] = exp(REAL(theta)[9]);
      REAL(ret)[24] = REAL(theta)[10];
      REAL(ret)[25] = REAL(theta)[11];
      REAL(ret)[26] = REAL(theta)[12];
      REAL(ret)[27] = REAL(theta)[13];
      REAL(ret)[28] = exp(REAL(theta)[14]);
      REAL(ret)[30] = REAL(theta)[15];
      REAL(ret)[31] = REAL(theta)[16];
      REAL(ret)[32] = REAL(theta)[17];
      REAL(ret)[33] = REAL(theta)[18];
      REAL(ret)[34] = REAL(theta)[19];
      REAL(ret)[35] = exp(REAL(theta)[20]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[7] = R_pow_di(REAL(theta)[1], 2) + exp(2 * REAL(theta)[2]);
      REAL(ret)[8] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[9] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[10] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[11] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[12] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[14] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + exp(2 * REAL(theta)[5]);
      REAL(ret)[15] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[16] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[17] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[18] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[19] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[20] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[21] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + exp(2 * REAL(theta)[9]);
      REAL(ret)[22] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[23] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[24] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[25] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[26] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[27] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[28] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + exp(2 * REAL(theta)[14]);
      REAL(ret)[29] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[30] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[31] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[32] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[33] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[34] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[35] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + exp(2 * REAL(theta)[20]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[12] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[18] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[24] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[30] = REAL(theta)[15] * exp(REAL(theta)[0]);
    }
    else if (theta_n == 2){
      REAL(ret)[1] = exp(REAL(theta)[0]);
      REAL(ret)[6] = exp(REAL(theta)[0]);
      REAL(ret)[7] = 2 * REAL(theta)[1];
      REAL(ret)[8] = REAL(theta)[3];
      REAL(ret)[9] = REAL(theta)[6];
      REAL(ret)[10] = REAL(theta)[10];
      REAL(ret)[11] = REAL(theta)[15];
      REAL(ret)[13] = REAL(theta)[3];
      REAL(ret)[19] = REAL(theta)[6];
      REAL(ret)[25] = REAL(theta)[10];
      REAL(ret)[31] = REAL(theta)[15];
    }
    else if (theta_n == 3){
      REAL(ret)[7] = 2 * exp(2 * REAL(theta)[2]);
      REAL(ret)[8] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[9] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[10] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[11] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[19] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[25] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[31] = REAL(theta)[16] * exp(REAL(theta)[2]);
    }
    else if (theta_n == 4){
      REAL(ret)[2] = exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[1];
      REAL(ret)[12] = exp(REAL(theta)[0]);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[14] = 2 * REAL(theta)[3];
      REAL(ret)[15] = REAL(theta)[6];
      REAL(ret)[16] = REAL(theta)[10];
      REAL(ret)[17] = REAL(theta)[15];
      REAL(ret)[20] = REAL(theta)[6];
      REAL(ret)[26] = REAL(theta)[10];
      REAL(ret)[32] = REAL(theta)[15];
    }
    else if (theta_n == 5){
      REAL(ret)[8] = exp(REAL(theta)[2]);
      REAL(ret)[13] = exp(REAL(theta)[2]);
      REAL(ret)[14] = 2 * REAL(theta)[4];
      REAL(ret)[15] = REAL(theta)[7];
      REAL(ret)[16] = REAL(theta)[11];
      REAL(ret)[17] = REAL(theta)[16];
      REAL(ret)[20] = REAL(theta)[7];
      REAL(ret)[26] = REAL(theta)[11];
      REAL(ret)[32] = REAL(theta)[16];
    }
    else if (theta_n == 6){
      REAL(ret)[14] = 2 * exp(2 * REAL(theta)[5]);
      REAL(ret)[15] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[16] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[17] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[20] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[26] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[32] = REAL(theta)[17] * exp(REAL(theta)[5]);
    }
    else if (theta_n == 7){
      REAL(ret)[3] = exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[1];
      REAL(ret)[15] = REAL(theta)[3];
      REAL(ret)[18] = exp(REAL(theta)[0]);
      REAL(ret)[19] = REAL(theta)[1];
      REAL(ret)[20] = REAL(theta)[3];
      REAL(ret)[21] = 2 * REAL(theta)[6];
      REAL(ret)[22] = REAL(theta)[10];
      REAL(ret)[23] = REAL(theta)[15];
      REAL(ret)[27] = REAL(theta)[10];
      REAL(ret)[33] = REAL(theta)[15];
    }
    else if (theta_n == 8){
      REAL(ret)[9] = exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[4];
      REAL(ret)[19] = exp(REAL(theta)[2]);
      REAL(ret)[20] = REAL(theta)[4];
      REAL(ret)[21] = 2 * REAL(theta)[7];
      REAL(ret)[22] = REAL(theta)[11];
      REAL(ret)[23] = REAL(theta)[16];
      REAL(ret)[27] = REAL(theta)[11];
      REAL(ret)[33] = REAL(theta)[16];
    }
    else if (theta_n == 9){
      REAL(ret)[15] = exp(REAL(theta)[5]);
      REAL(ret)[20] = exp(REAL(theta)[5]);
      REAL(ret)[21] = 2 * REAL(theta)[8];
      REAL(ret)[22] = REAL(theta)[12];
      REAL(ret)[23] = REAL(theta)[17];
      REAL(ret)[27] = REAL(theta)[12];
      REAL(ret)[33] = REAL(theta)[17];
    }
    else if (theta_n == 10){
      REAL(ret)[21] = 2 * exp(2 * REAL(theta)[9]);
      REAL(ret)[22] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[23] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[27] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[33] = REAL(theta)[18] * exp(REAL(theta)[9]);
    }
    else if (theta_n == 11){
      REAL(ret)[4] = exp(REAL(theta)[0]);
      REAL(ret)[10] = REAL(theta)[1];
      REAL(ret)[16] = REAL(theta)[3];
      REAL(ret)[22] = REAL(theta)[6];
      REAL(ret)[24] = exp(REAL(theta)[0]);
      REAL(ret)[25] = REAL(theta)[1];
      REAL(ret)[26] = REAL(theta)[3];
      REAL(ret)[27] = REAL(theta)[6];
      REAL(ret)[28] = 2 * REAL(theta)[10];
      REAL(ret)[29] = REAL(theta)[15];
      REAL(ret)[34] = REAL(theta)[15];
    }
    else if (theta_n == 12){
      REAL(ret)[10] = exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[4];
      REAL(ret)[22] = REAL(theta)[7];
      REAL(ret)[25] = exp(REAL(theta)[2]);
      REAL(ret)[26] = REAL(theta)[4];
      REAL(ret)[27] = REAL(theta)[7];
      REAL(ret)[28] = 2 * REAL(theta)[11];
      REAL(ret)[29] = REAL(theta)[16];
      REAL(ret)[34] = REAL(theta)[16];
    }
    else if (theta_n == 13){
      REAL(ret)[16] = exp(REAL(theta)[5]);
      REAL(ret)[22] = REAL(theta)[8];
      REAL(ret)[26] = exp(REAL(theta)[5]);
      REAL(ret)[27] = REAL(theta)[8];
      REAL(ret)[28] = 2 * REAL(theta)[12];
      REAL(ret)[29] = REAL(theta)[17];
      REAL(ret)[34] = REAL(theta)[17];
    }
    else if (theta_n == 14){
      REAL(ret)[22] = exp(REAL(theta)[9]);
      REAL(ret)[27] = exp(REAL(theta)[9]);
      REAL(ret)[28] = 2 * REAL(theta)[13];
      REAL(ret)[29] = REAL(theta)[18];
      REAL(ret)[34] = REAL(theta)[18];
    }
    else if (theta_n == 15){
      REAL(ret)[28] = 2 * exp(2 * REAL(theta)[14]);
      REAL(ret)[29] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[34] = REAL(theta)[19] * exp(REAL(theta)[14]);
    }
    else if (theta_n == 16){
      REAL(ret)[5] = exp(REAL(theta)[0]);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[17] = REAL(theta)[3];
      REAL(ret)[23] = REAL(theta)[6];
      REAL(ret)[29] = REAL(theta)[10];
      REAL(ret)[30] = exp(REAL(theta)[0]);
      REAL(ret)[31] = REAL(theta)[1];
      REAL(ret)[32] = REAL(theta)[3];
      REAL(ret)[33] = REAL(theta)[6];
      REAL(ret)[34] = REAL(theta)[10];
      REAL(ret)[35] = 2 * REAL(theta)[15];
    }
    else if (theta_n == 17){
      REAL(ret)[11] = exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[4];
      REAL(ret)[23] = REAL(theta)[7];
      REAL(ret)[29] = REAL(theta)[11];
      REAL(ret)[31] = exp(REAL(theta)[2]);
      REAL(ret)[32] = REAL(theta)[4];
      REAL(ret)[33] = REAL(theta)[7];
      REAL(ret)[34] = REAL(theta)[11];
      REAL(ret)[35] = 2 * REAL(theta)[16];
    }
    else if (theta_n == 18){
      REAL(ret)[17] = exp(REAL(theta)[5]);
      REAL(ret)[23] = REAL(theta)[8];
      REAL(ret)[29] = REAL(theta)[12];
      REAL(ret)[32] = exp(REAL(theta)[5]);
      REAL(ret)[33] = REAL(theta)[8];
      REAL(ret)[34] = REAL(theta)[12];
      REAL(ret)[35] = 2 * REAL(theta)[17];
    }
    else if (theta_n == 19){
      REAL(ret)[23] = exp(REAL(theta)[9]);
      REAL(ret)[29] = REAL(theta)[13];
      REAL(ret)[33] = exp(REAL(theta)[9]);
      REAL(ret)[34] = REAL(theta)[13];
      REAL(ret)[35] = 2 * REAL(theta)[18];
    }
    else if (theta_n == 20){
      REAL(ret)[29] = exp(REAL(theta)[14]);
      REAL(ret)[34] = exp(REAL(theta)[14]);
      REAL(ret)[35] = 2 * REAL(theta)[19];
    }
    else if (theta_n == 21){
      REAL(ret)[35] = 2 * exp(2 * REAL(theta)[20]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 6));for(int i = 0; i < 6; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 2 * exp(2 * REAL(theta)[2]);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 2 * exp(2 * REAL(theta)[5]);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 2 * exp(2 * REAL(theta)[9]);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 2 * exp(2 * REAL(theta)[14]);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 2 * exp(2 * REAL(theta)[20]);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 7){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 7;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -30 || theta_n > 28){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 28){
    error("Requires vector with 28 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 7, 7));for (int i = 0; i < 49; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[1];
      REAL(ret)[8] = exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[3];
      REAL(ret)[15] = REAL(theta)[4];
      REAL(ret)[16] = exp(REAL(theta)[5]);
      REAL(ret)[21] = REAL(theta)[6];
      REAL(ret)[22] = REAL(theta)[7];
      REAL(ret)[23] = REAL(theta)[8];
      REAL(ret)[24] = exp(REAL(theta)[9]);
      REAL(ret)[28] = REAL(theta)[10];
      REAL(ret)[29] = REAL(theta)[11];
      REAL(ret)[30] = REAL(theta)[12];
      REAL(ret)[31] = REAL(theta)[13];
      REAL(ret)[32] = exp(REAL(theta)[14]);
      REAL(ret)[35] = REAL(theta)[15];
      REAL(ret)[36] = REAL(theta)[16];
      REAL(ret)[37] = REAL(theta)[17];
      REAL(ret)[38] = REAL(theta)[18];
      REAL(ret)[39] = REAL(theta)[19];
      REAL(ret)[40] = exp(REAL(theta)[20]);
      REAL(ret)[42] = REAL(theta)[21];
      REAL(ret)[43] = REAL(theta)[22];
      REAL(ret)[44] = REAL(theta)[23];
      REAL(ret)[45] = REAL(theta)[24];
      REAL(ret)[46] = REAL(theta)[25];
      REAL(ret)[47] = REAL(theta)[26];
      REAL(ret)[48] = exp(REAL(theta)[27]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[8] = R_pow_di(REAL(theta)[1], 2) + exp(2 * REAL(theta)[2]);
      REAL(ret)[9] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[10] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[11] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[12] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[16] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + exp(2 * REAL(theta)[5]);
      REAL(ret)[17] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[18] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[19] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[20] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[21] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[22] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[23] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[24] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + exp(2 * REAL(theta)[9]);
      REAL(ret)[25] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[26] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[27] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[28] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[29] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[30] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[31] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[32] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + exp(2 * REAL(theta)[14]);
      REAL(ret)[33] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[34] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[35] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[36] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[37] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[38] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[39] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[40] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + exp(2 * REAL(theta)[20]);
      REAL(ret)[41] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[42] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[43] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[44] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[45] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[46] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[47] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[48] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + exp(2 * REAL(theta)[27]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[14] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[21] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[28] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[35] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[42] = REAL(theta)[21] * exp(REAL(theta)[0]);
    }
    else if (theta_n == 2){
      REAL(ret)[1] = exp(REAL(theta)[0]);
      REAL(ret)[7] = exp(REAL(theta)[0]);
      REAL(ret)[8] = 2 * REAL(theta)[1];
      REAL(ret)[9] = REAL(theta)[3];
      REAL(ret)[10] = REAL(theta)[6];
      REAL(ret)[11] = REAL(theta)[10];
      REAL(ret)[12] = REAL(theta)[15];
      REAL(ret)[13] = REAL(theta)[21];
      REAL(ret)[15] = REAL(theta)[3];
      REAL(ret)[22] = REAL(theta)[6];
      REAL(ret)[29] = REAL(theta)[10];
      REAL(ret)[36] = REAL(theta)[15];
      REAL(ret)[43] = REAL(theta)[21];
    }
    else if (theta_n == 3){
      REAL(ret)[8] = 2 * exp(2 * REAL(theta)[2]);
      REAL(ret)[9] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[10] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[11] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[12] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[22] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[29] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[36] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[43] = REAL(theta)[22] * exp(REAL(theta)[2]);
    }
    else if (theta_n == 4){
      REAL(ret)[2] = exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[1];
      REAL(ret)[14] = exp(REAL(theta)[0]);
      REAL(ret)[15] = REAL(theta)[1];
      REAL(ret)[16] = 2 * REAL(theta)[3];
      REAL(ret)[17] = REAL(theta)[6];
      REAL(ret)[18] = REAL(theta)[10];
      REAL(ret)[19] = REAL(theta)[15];
      REAL(ret)[20] = REAL(theta)[21];
      REAL(ret)[23] = REAL(theta)[6];
      REAL(ret)[30] = REAL(theta)[10];
      REAL(ret)[37] = REAL(theta)[15];
      REAL(ret)[44] = REAL(theta)[21];
    }
    else if (theta_n == 5){
      REAL(ret)[9] = exp(REAL(theta)[2]);
      REAL(ret)[15] = exp(REAL(theta)[2]);
      REAL(ret)[16] = 2 * REAL(theta)[4];
      REAL(ret)[17] = REAL(theta)[7];
      REAL(ret)[18] = REAL(theta)[11];
      REAL(ret)[19] = REAL(theta)[16];
      REAL(ret)[20] = REAL(theta)[22];
      REAL(ret)[23] = REAL(theta)[7];
      REAL(ret)[30] = REAL(theta)[11];
      REAL(ret)[37] = REAL(theta)[16];
      REAL(ret)[44] = REAL(theta)[22];
    }
    else if (theta_n == 6){
      REAL(ret)[16] = 2 * exp(2 * REAL(theta)[5]);
      REAL(ret)[17] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[18] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[19] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[20] = REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[23] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[30] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[37] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[44] = REAL(theta)[23] * exp(REAL(theta)[5]);
    }
    else if (theta_n == 7){
      REAL(ret)[3] = exp(REAL(theta)[0]);
      REAL(ret)[10] = REAL(theta)[1];
      REAL(ret)[17] = REAL(theta)[3];
      REAL(ret)[21] = exp(REAL(theta)[0]);
      REAL(ret)[22] = REAL(theta)[1];
      REAL(ret)[23] = REAL(theta)[3];
      REAL(ret)[24] = 2 * REAL(theta)[6];
      REAL(ret)[25] = REAL(theta)[10];
      REAL(ret)[26] = REAL(theta)[15];
      REAL(ret)[27] = REAL(theta)[21];
      REAL(ret)[31] = REAL(theta)[10];
      REAL(ret)[38] = REAL(theta)[15];
      REAL(ret)[45] = REAL(theta)[21];
    }
    else if (theta_n == 8){
      REAL(ret)[10] = exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[4];
      REAL(ret)[22] = exp(REAL(theta)[2]);
      REAL(ret)[23] = REAL(theta)[4];
      REAL(ret)[24] = 2 * REAL(theta)[7];
      REAL(ret)[25] = REAL(theta)[11];
      REAL(ret)[26] = REAL(theta)[16];
      REAL(ret)[27] = REAL(theta)[22];
      REAL(ret)[31] = REAL(theta)[11];
      REAL(ret)[38] = REAL(theta)[16];
      REAL(ret)[45] = REAL(theta)[22];
    }
    else if (theta_n == 9){
      REAL(ret)[17] = exp(REAL(theta)[5]);
      REAL(ret)[23] = exp(REAL(theta)[5]);
      REAL(ret)[24] = 2 * REAL(theta)[8];
      REAL(ret)[25] = REAL(theta)[12];
      REAL(ret)[26] = REAL(theta)[17];
      REAL(ret)[27] = REAL(theta)[23];
      REAL(ret)[31] = REAL(theta)[12];
      REAL(ret)[38] = REAL(theta)[17];
      REAL(ret)[45] = REAL(theta)[23];
    }
    else if (theta_n == 10){
      REAL(ret)[24] = 2 * exp(2 * REAL(theta)[9]);
      REAL(ret)[25] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[26] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[27] = REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[31] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[38] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[45] = REAL(theta)[24] * exp(REAL(theta)[9]);
    }
    else if (theta_n == 11){
      REAL(ret)[4] = exp(REAL(theta)[0]);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[18] = REAL(theta)[3];
      REAL(ret)[25] = REAL(theta)[6];
      REAL(ret)[28] = exp(REAL(theta)[0]);
      REAL(ret)[29] = REAL(theta)[1];
      REAL(ret)[30] = REAL(theta)[3];
      REAL(ret)[31] = REAL(theta)[6];
      REAL(ret)[32] = 2 * REAL(theta)[10];
      REAL(ret)[33] = REAL(theta)[15];
      REAL(ret)[34] = REAL(theta)[21];
      REAL(ret)[39] = REAL(theta)[15];
      REAL(ret)[46] = REAL(theta)[21];
    }
    else if (theta_n == 12){
      REAL(ret)[11] = exp(REAL(theta)[2]);
      REAL(ret)[18] = REAL(theta)[4];
      REAL(ret)[25] = REAL(theta)[7];
      REAL(ret)[29] = exp(REAL(theta)[2]);
      REAL(ret)[30] = REAL(theta)[4];
      REAL(ret)[31] = REAL(theta)[7];
      REAL(ret)[32] = 2 * REAL(theta)[11];
      REAL(ret)[33] = REAL(theta)[16];
      REAL(ret)[34] = REAL(theta)[22];
      REAL(ret)[39] = REAL(theta)[16];
      REAL(ret)[46] = REAL(theta)[22];
    }
    else if (theta_n == 13){
      REAL(ret)[18] = exp(REAL(theta)[5]);
      REAL(ret)[25] = REAL(theta)[8];
      REAL(ret)[30] = exp(REAL(theta)[5]);
      REAL(ret)[31] = REAL(theta)[8];
      REAL(ret)[32] = 2 * REAL(theta)[12];
      REAL(ret)[33] = REAL(theta)[17];
      REAL(ret)[34] = REAL(theta)[23];
      REAL(ret)[39] = REAL(theta)[17];
      REAL(ret)[46] = REAL(theta)[23];
    }
    else if (theta_n == 14){
      REAL(ret)[25] = exp(REAL(theta)[9]);
      REAL(ret)[31] = exp(REAL(theta)[9]);
      REAL(ret)[32] = 2 * REAL(theta)[13];
      REAL(ret)[33] = REAL(theta)[18];
      REAL(ret)[34] = REAL(theta)[24];
      REAL(ret)[39] = REAL(theta)[18];
      REAL(ret)[46] = REAL(theta)[24];
    }
    else if (theta_n == 15){
      REAL(ret)[32] = 2 * exp(2 * REAL(theta)[14]);
      REAL(ret)[33] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[34] = REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[39] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[46] = REAL(theta)[25] * exp(REAL(theta)[14]);
    }
    else if (theta_n == 16){
      REAL(ret)[5] = exp(REAL(theta)[0]);
      REAL(ret)[12] = REAL(theta)[1];
      REAL(ret)[19] = REAL(theta)[3];
      REAL(ret)[26] = REAL(theta)[6];
      REAL(ret)[33] = REAL(theta)[10];
      REAL(ret)[35] = exp(REAL(theta)[0]);
      REAL(ret)[36] = REAL(theta)[1];
      REAL(ret)[37] = REAL(theta)[3];
      REAL(ret)[38] = REAL(theta)[6];
      REAL(ret)[39] = REAL(theta)[10];
      REAL(ret)[40] = 2 * REAL(theta)[15];
      REAL(ret)[41] = REAL(theta)[21];
      REAL(ret)[47] = REAL(theta)[21];
    }
    else if (theta_n == 17){
      REAL(ret)[12] = exp(REAL(theta)[2]);
      REAL(ret)[19] = REAL(theta)[4];
      REAL(ret)[26] = REAL(theta)[7];
      REAL(ret)[33] = REAL(theta)[11];
      REAL(ret)[36] = exp(REAL(theta)[2]);
      REAL(ret)[37] = REAL(theta)[4];
      REAL(ret)[38] = REAL(theta)[7];
      REAL(ret)[39] = REAL(theta)[11];
      REAL(ret)[40] = 2 * REAL(theta)[16];
      REAL(ret)[41] = REAL(theta)[22];
      REAL(ret)[47] = REAL(theta)[22];
    }
    else if (theta_n == 18){
      REAL(ret)[19] = exp(REAL(theta)[5]);
      REAL(ret)[26] = REAL(theta)[8];
      REAL(ret)[33] = REAL(theta)[12];
      REAL(ret)[37] = exp(REAL(theta)[5]);
      REAL(ret)[38] = REAL(theta)[8];
      REAL(ret)[39] = REAL(theta)[12];
      REAL(ret)[40] = 2 * REAL(theta)[17];
      REAL(ret)[41] = REAL(theta)[23];
      REAL(ret)[47] = REAL(theta)[23];
    }
    else if (theta_n == 19){
      REAL(ret)[26] = exp(REAL(theta)[9]);
      REAL(ret)[33] = REAL(theta)[13];
      REAL(ret)[38] = exp(REAL(theta)[9]);
      REAL(ret)[39] = REAL(theta)[13];
      REAL(ret)[40] = 2 * REAL(theta)[18];
      REAL(ret)[41] = REAL(theta)[24];
      REAL(ret)[47] = REAL(theta)[24];
    }
    else if (theta_n == 20){
      REAL(ret)[33] = exp(REAL(theta)[14]);
      REAL(ret)[39] = exp(REAL(theta)[14]);
      REAL(ret)[40] = 2 * REAL(theta)[19];
      REAL(ret)[41] = REAL(theta)[25];
      REAL(ret)[47] = REAL(theta)[25];
    }
    else if (theta_n == 21){
      REAL(ret)[40] = 2 * exp(2 * REAL(theta)[20]);
      REAL(ret)[41] = REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[47] = REAL(theta)[26] * exp(REAL(theta)[20]);
    }
    else if (theta_n == 22){
      REAL(ret)[6] = exp(REAL(theta)[0]);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[20] = REAL(theta)[3];
      REAL(ret)[27] = REAL(theta)[6];
      REAL(ret)[34] = REAL(theta)[10];
      REAL(ret)[41] = REAL(theta)[15];
      REAL(ret)[42] = exp(REAL(theta)[0]);
      REAL(ret)[43] = REAL(theta)[1];
      REAL(ret)[44] = REAL(theta)[3];
      REAL(ret)[45] = REAL(theta)[6];
      REAL(ret)[46] = REAL(theta)[10];
      REAL(ret)[47] = REAL(theta)[15];
      REAL(ret)[48] = 2 * REAL(theta)[21];
    }
    else if (theta_n == 23){
      REAL(ret)[13] = exp(REAL(theta)[2]);
      REAL(ret)[20] = REAL(theta)[4];
      REAL(ret)[27] = REAL(theta)[7];
      REAL(ret)[34] = REAL(theta)[11];
      REAL(ret)[41] = REAL(theta)[16];
      REAL(ret)[43] = exp(REAL(theta)[2]);
      REAL(ret)[44] = REAL(theta)[4];
      REAL(ret)[45] = REAL(theta)[7];
      REAL(ret)[46] = REAL(theta)[11];
      REAL(ret)[47] = REAL(theta)[16];
      REAL(ret)[48] = 2 * REAL(theta)[22];
    }
    else if (theta_n == 24){
      REAL(ret)[20] = exp(REAL(theta)[5]);
      REAL(ret)[27] = REAL(theta)[8];
      REAL(ret)[34] = REAL(theta)[12];
      REAL(ret)[41] = REAL(theta)[17];
      REAL(ret)[44] = exp(REAL(theta)[5]);
      REAL(ret)[45] = REAL(theta)[8];
      REAL(ret)[46] = REAL(theta)[12];
      REAL(ret)[47] = REAL(theta)[17];
      REAL(ret)[48] = 2 * REAL(theta)[23];
    }
    else if (theta_n == 25){
      REAL(ret)[27] = exp(REAL(theta)[9]);
      REAL(ret)[34] = REAL(theta)[13];
      REAL(ret)[41] = REAL(theta)[18];
      REAL(ret)[45] = exp(REAL(theta)[9]);
      REAL(ret)[46] = REAL(theta)[13];
      REAL(ret)[47] = REAL(theta)[18];
      REAL(ret)[48] = 2 * REAL(theta)[24];
    }
    else if (theta_n == 26){
      REAL(ret)[34] = exp(REAL(theta)[14]);
      REAL(ret)[41] = REAL(theta)[19];
      REAL(ret)[46] = exp(REAL(theta)[14]);
      REAL(ret)[47] = REAL(theta)[19];
      REAL(ret)[48] = 2 * REAL(theta)[25];
    }
    else if (theta_n == 27){
      REAL(ret)[41] = exp(REAL(theta)[20]);
      REAL(ret)[47] = exp(REAL(theta)[20]);
      REAL(ret)[48] = 2 * REAL(theta)[26];
    }
    else if (theta_n == 28){
      REAL(ret)[48] = 2 * exp(2 * REAL(theta)[27]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 7));for(int i = 0; i < 7; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 2 * exp(2 * REAL(theta)[2]);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 2 * exp(2 * REAL(theta)[5]);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 2 * exp(2 * REAL(theta)[9]);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 2 * exp(2 * REAL(theta)[14]);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 2 * exp(2 * REAL(theta)[20]);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 2 * exp(2 * REAL(theta)[27]);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 8){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 8;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -38 || theta_n > 36){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 36){
    error("Requires vector with 36 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 8, 8));for (int i = 0; i < 64; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[1];
      REAL(ret)[9] = exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[3];
      REAL(ret)[17] = REAL(theta)[4];
      REAL(ret)[18] = exp(REAL(theta)[5]);
      REAL(ret)[24] = REAL(theta)[6];
      REAL(ret)[25] = REAL(theta)[7];
      REAL(ret)[26] = REAL(theta)[8];
      REAL(ret)[27] = exp(REAL(theta)[9]);
      REAL(ret)[32] = REAL(theta)[10];
      REAL(ret)[33] = REAL(theta)[11];
      REAL(ret)[34] = REAL(theta)[12];
      REAL(ret)[35] = REAL(theta)[13];
      REAL(ret)[36] = exp(REAL(theta)[14]);
      REAL(ret)[40] = REAL(theta)[15];
      REAL(ret)[41] = REAL(theta)[16];
      REAL(ret)[42] = REAL(theta)[17];
      REAL(ret)[43] = REAL(theta)[18];
      REAL(ret)[44] = REAL(theta)[19];
      REAL(ret)[45] = exp(REAL(theta)[20]);
      REAL(ret)[48] = REAL(theta)[21];
      REAL(ret)[49] = REAL(theta)[22];
      REAL(ret)[50] = REAL(theta)[23];
      REAL(ret)[51] = REAL(theta)[24];
      REAL(ret)[52] = REAL(theta)[25];
      REAL(ret)[53] = REAL(theta)[26];
      REAL(ret)[54] = exp(REAL(theta)[27]);
      REAL(ret)[56] = REAL(theta)[28];
      REAL(ret)[57] = REAL(theta)[29];
      REAL(ret)[58] = REAL(theta)[30];
      REAL(ret)[59] = REAL(theta)[31];
      REAL(ret)[60] = REAL(theta)[32];
      REAL(ret)[61] = REAL(theta)[33];
      REAL(ret)[62] = REAL(theta)[34];
      REAL(ret)[63] = exp(REAL(theta)[35]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[9] = R_pow_di(REAL(theta)[1], 2) + exp(2 * REAL(theta)[2]);
      REAL(ret)[10] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[11] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[12] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[28] + REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[17] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[18] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + exp(2 * REAL(theta)[5]);
      REAL(ret)[19] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[20] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[21] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[22] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[23] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[24] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[25] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[26] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[27] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + exp(2 * REAL(theta)[9]);
      REAL(ret)[28] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[29] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[30] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[31] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[32] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[33] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[34] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[35] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[36] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + exp(2 * REAL(theta)[14]);
      REAL(ret)[37] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[38] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[39] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[40] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[41] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[42] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[43] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[44] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[45] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + exp(2 * REAL(theta)[20]);
      REAL(ret)[46] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[47] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[48] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[49] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[50] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[51] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[52] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[53] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[54] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + exp(2 * REAL(theta)[27]);
      REAL(ret)[55] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[56] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[57] = REAL(theta)[1] * REAL(theta)[28] + REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[58] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[59] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[60] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[61] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[62] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[63] = R_pow_di(REAL(theta)[28], 2) + R_pow_di(REAL(theta)[29], 2) + R_pow_di(REAL(theta)[30], 2) + R_pow_di(REAL(theta)[31], 2) + R_pow_di(REAL(theta)[32], 2) + R_pow_di(REAL(theta)[33], 2) + R_pow_di(REAL(theta)[34], 2) + exp(2 * REAL(theta)[35]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[16] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[24] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[32] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[40] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[48] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[56] = REAL(theta)[28] * exp(REAL(theta)[0]);
    }
    else if (theta_n == 2){
      REAL(ret)[1] = exp(REAL(theta)[0]);
      REAL(ret)[8] = exp(REAL(theta)[0]);
      REAL(ret)[9] = 2 * REAL(theta)[1];
      REAL(ret)[10] = REAL(theta)[3];
      REAL(ret)[11] = REAL(theta)[6];
      REAL(ret)[12] = REAL(theta)[10];
      REAL(ret)[13] = REAL(theta)[15];
      REAL(ret)[14] = REAL(theta)[21];
      REAL(ret)[15] = REAL(theta)[28];
      REAL(ret)[17] = REAL(theta)[3];
      REAL(ret)[25] = REAL(theta)[6];
      REAL(ret)[33] = REAL(theta)[10];
      REAL(ret)[41] = REAL(theta)[15];
      REAL(ret)[49] = REAL(theta)[21];
      REAL(ret)[57] = REAL(theta)[28];
    }
    else if (theta_n == 3){
      REAL(ret)[9] = 2 * exp(2 * REAL(theta)[2]);
      REAL(ret)[10] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[11] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[12] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[25] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[33] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[41] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[49] = REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[57] = REAL(theta)[29] * exp(REAL(theta)[2]);
    }
    else if (theta_n == 4){
      REAL(ret)[2] = exp(REAL(theta)[0]);
      REAL(ret)[10] = REAL(theta)[1];
      REAL(ret)[16] = exp(REAL(theta)[0]);
      REAL(ret)[17] = REAL(theta)[1];
      REAL(ret)[18] = 2 * REAL(theta)[3];
      REAL(ret)[19] = REAL(theta)[6];
      REAL(ret)[20] = REAL(theta)[10];
      REAL(ret)[21] = REAL(theta)[15];
      REAL(ret)[22] = REAL(theta)[21];
      REAL(ret)[23] = REAL(theta)[28];
      REAL(ret)[26] = REAL(theta)[6];
      REAL(ret)[34] = REAL(theta)[10];
      REAL(ret)[42] = REAL(theta)[15];
      REAL(ret)[50] = REAL(theta)[21];
      REAL(ret)[58] = REAL(theta)[28];
    }
    else if (theta_n == 5){
      REAL(ret)[10] = exp(REAL(theta)[2]);
      REAL(ret)[17] = exp(REAL(theta)[2]);
      REAL(ret)[18] = 2 * REAL(theta)[4];
      REAL(ret)[19] = REAL(theta)[7];
      REAL(ret)[20] = REAL(theta)[11];
      REAL(ret)[21] = REAL(theta)[16];
      REAL(ret)[22] = REAL(theta)[22];
      REAL(ret)[23] = REAL(theta)[29];
      REAL(ret)[26] = REAL(theta)[7];
      REAL(ret)[34] = REAL(theta)[11];
      REAL(ret)[42] = REAL(theta)[16];
      REAL(ret)[50] = REAL(theta)[22];
      REAL(ret)[58] = REAL(theta)[29];
    }
    else if (theta_n == 6){
      REAL(ret)[18] = 2 * exp(2 * REAL(theta)[5]);
      REAL(ret)[19] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[20] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[21] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[22] = REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[23] = REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[26] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[34] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[42] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[50] = REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[58] = REAL(theta)[30] * exp(REAL(theta)[5]);
    }
    else if (theta_n == 7){
      REAL(ret)[3] = exp(REAL(theta)[0]);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[19] = REAL(theta)[3];
      REAL(ret)[24] = exp(REAL(theta)[0]);
      REAL(ret)[25] = REAL(theta)[1];
      REAL(ret)[26] = REAL(theta)[3];
      REAL(ret)[27] = 2 * REAL(theta)[6];
      REAL(ret)[28] = REAL(theta)[10];
      REAL(ret)[29] = REAL(theta)[15];
      REAL(ret)[30] = REAL(theta)[21];
      REAL(ret)[31] = REAL(theta)[28];
      REAL(ret)[35] = REAL(theta)[10];
      REAL(ret)[43] = REAL(theta)[15];
      REAL(ret)[51] = REAL(theta)[21];
      REAL(ret)[59] = REAL(theta)[28];
    }
    else if (theta_n == 8){
      REAL(ret)[11] = exp(REAL(theta)[2]);
      REAL(ret)[19] = REAL(theta)[4];
      REAL(ret)[25] = exp(REAL(theta)[2]);
      REAL(ret)[26] = REAL(theta)[4];
      REAL(ret)[27] = 2 * REAL(theta)[7];
      REAL(ret)[28] = REAL(theta)[11];
      REAL(ret)[29] = REAL(theta)[16];
      REAL(ret)[30] = REAL(theta)[22];
      REAL(ret)[31] = REAL(theta)[29];
      REAL(ret)[35] = REAL(theta)[11];
      REAL(ret)[43] = REAL(theta)[16];
      REAL(ret)[51] = REAL(theta)[22];
      REAL(ret)[59] = REAL(theta)[29];
    }
    else if (theta_n == 9){
      REAL(ret)[19] = exp(REAL(theta)[5]);
      REAL(ret)[26] = exp(REAL(theta)[5]);
      REAL(ret)[27] = 2 * REAL(theta)[8];
      REAL(ret)[28] = REAL(theta)[12];
      REAL(ret)[29] = REAL(theta)[17];
      REAL(ret)[30] = REAL(theta)[23];
      REAL(ret)[31] = REAL(theta)[30];
      REAL(ret)[35] = REAL(theta)[12];
      REAL(ret)[43] = REAL(theta)[17];
      REAL(ret)[51] = REAL(theta)[23];
      REAL(ret)[59] = REAL(theta)[30];
    }
    else if (theta_n == 10){
      REAL(ret)[27] = 2 * exp(2 * REAL(theta)[9]);
      REAL(ret)[28] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[29] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[30] = REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[31] = REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[35] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[43] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[51] = REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[59] = REAL(theta)[31] * exp(REAL(theta)[9]);
    }
    else if (theta_n == 11){
      REAL(ret)[4] = exp(REAL(theta)[0]);
      REAL(ret)[12] = REAL(theta)[1];
      REAL(ret)[20] = REAL(theta)[3];
      REAL(ret)[28] = REAL(theta)[6];
      REAL(ret)[32] = exp(REAL(theta)[0]);
      REAL(ret)[33] = REAL(theta)[1];
      REAL(ret)[34] = REAL(theta)[3];
      REAL(ret)[35] = REAL(theta)[6];
      REAL(ret)[36] = 2 * REAL(theta)[10];
      REAL(ret)[37] = REAL(theta)[15];
      REAL(ret)[38] = REAL(theta)[21];
      REAL(ret)[39] = REAL(theta)[28];
      REAL(ret)[44] = REAL(theta)[15];
      REAL(ret)[52] = REAL(theta)[21];
      REAL(ret)[60] = REAL(theta)[28];
    }
    else if (theta_n == 12){
      REAL(ret)[12] = exp(REAL(theta)[2]);
      REAL(ret)[20] = REAL(theta)[4];
      REAL(ret)[28] = REAL(theta)[7];
      REAL(ret)[33] = exp(REAL(theta)[2]);
      REAL(ret)[34] = REAL(theta)[4];
      REAL(ret)[35] = REAL(theta)[7];
      REAL(ret)[36] = 2 * REAL(theta)[11];
      REAL(ret)[37] = REAL(theta)[16];
      REAL(ret)[38] = REAL(theta)[22];
      REAL(ret)[39] = REAL(theta)[29];
      REAL(ret)[44] = REAL(theta)[16];
      REAL(ret)[52] = REAL(theta)[22];
      REAL(ret)[60] = REAL(theta)[29];
    }
    else if (theta_n == 13){
      REAL(ret)[20] = exp(REAL(theta)[5]);
      REAL(ret)[28] = REAL(theta)[8];
      REAL(ret)[34] = exp(REAL(theta)[5]);
      REAL(ret)[35] = REAL(theta)[8];
      REAL(ret)[36] = 2 * REAL(theta)[12];
      REAL(ret)[37] = REAL(theta)[17];
      REAL(ret)[38] = REAL(theta)[23];
      REAL(ret)[39] = REAL(theta)[30];
      REAL(ret)[44] = REAL(theta)[17];
      REAL(ret)[52] = REAL(theta)[23];
      REAL(ret)[60] = REAL(theta)[30];
    }
    else if (theta_n == 14){
      REAL(ret)[28] = exp(REAL(theta)[9]);
      REAL(ret)[35] = exp(REAL(theta)[9]);
      REAL(ret)[36] = 2 * REAL(theta)[13];
      REAL(ret)[37] = REAL(theta)[18];
      REAL(ret)[38] = REAL(theta)[24];
      REAL(ret)[39] = REAL(theta)[31];
      REAL(ret)[44] = REAL(theta)[18];
      REAL(ret)[52] = REAL(theta)[24];
      REAL(ret)[60] = REAL(theta)[31];
    }
    else if (theta_n == 15){
      REAL(ret)[36] = 2 * exp(2 * REAL(theta)[14]);
      REAL(ret)[37] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[38] = REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[39] = REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[44] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[52] = REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[60] = REAL(theta)[32] * exp(REAL(theta)[14]);
    }
    else if (theta_n == 16){
      REAL(ret)[5] = exp(REAL(theta)[0]);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[21] = REAL(theta)[3];
      REAL(ret)[29] = REAL(theta)[6];
      REAL(ret)[37] = REAL(theta)[10];
      REAL(ret)[40] = exp(REAL(theta)[0]);
      REAL(ret)[41] = REAL(theta)[1];
      REAL(ret)[42] = REAL(theta)[3];
      REAL(ret)[43] = REAL(theta)[6];
      REAL(ret)[44] = REAL(theta)[10];
      REAL(ret)[45] = 2 * REAL(theta)[15];
      REAL(ret)[46] = REAL(theta)[21];
      REAL(ret)[47] = REAL(theta)[28];
      REAL(ret)[53] = REAL(theta)[21];
      REAL(ret)[61] = REAL(theta)[28];
    }
    else if (theta_n == 17){
      REAL(ret)[13] = exp(REAL(theta)[2]);
      REAL(ret)[21] = REAL(theta)[4];
      REAL(ret)[29] = REAL(theta)[7];
      REAL(ret)[37] = REAL(theta)[11];
      REAL(ret)[41] = exp(REAL(theta)[2]);
      REAL(ret)[42] = REAL(theta)[4];
      REAL(ret)[43] = REAL(theta)[7];
      REAL(ret)[44] = REAL(theta)[11];
      REAL(ret)[45] = 2 * REAL(theta)[16];
      REAL(ret)[46] = REAL(theta)[22];
      REAL(ret)[47] = REAL(theta)[29];
      REAL(ret)[53] = REAL(theta)[22];
      REAL(ret)[61] = REAL(theta)[29];
    }
    else if (theta_n == 18){
      REAL(ret)[21] = exp(REAL(theta)[5]);
      REAL(ret)[29] = REAL(theta)[8];
      REAL(ret)[37] = REAL(theta)[12];
      REAL(ret)[42] = exp(REAL(theta)[5]);
      REAL(ret)[43] = REAL(theta)[8];
      REAL(ret)[44] = REAL(theta)[12];
      REAL(ret)[45] = 2 * REAL(theta)[17];
      REAL(ret)[46] = REAL(theta)[23];
      REAL(ret)[47] = REAL(theta)[30];
      REAL(ret)[53] = REAL(theta)[23];
      REAL(ret)[61] = REAL(theta)[30];
    }
    else if (theta_n == 19){
      REAL(ret)[29] = exp(REAL(theta)[9]);
      REAL(ret)[37] = REAL(theta)[13];
      REAL(ret)[43] = exp(REAL(theta)[9]);
      REAL(ret)[44] = REAL(theta)[13];
      REAL(ret)[45] = 2 * REAL(theta)[18];
      REAL(ret)[46] = REAL(theta)[24];
      REAL(ret)[47] = REAL(theta)[31];
      REAL(ret)[53] = REAL(theta)[24];
      REAL(ret)[61] = REAL(theta)[31];
    }
    else if (theta_n == 20){
      REAL(ret)[37] = exp(REAL(theta)[14]);
      REAL(ret)[44] = exp(REAL(theta)[14]);
      REAL(ret)[45] = 2 * REAL(theta)[19];
      REAL(ret)[46] = REAL(theta)[25];
      REAL(ret)[47] = REAL(theta)[32];
      REAL(ret)[53] = REAL(theta)[25];
      REAL(ret)[61] = REAL(theta)[32];
    }
    else if (theta_n == 21){
      REAL(ret)[45] = 2 * exp(2 * REAL(theta)[20]);
      REAL(ret)[46] = REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[47] = REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[53] = REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[61] = REAL(theta)[33] * exp(REAL(theta)[20]);
    }
    else if (theta_n == 22){
      REAL(ret)[6] = exp(REAL(theta)[0]);
      REAL(ret)[14] = REAL(theta)[1];
      REAL(ret)[22] = REAL(theta)[3];
      REAL(ret)[30] = REAL(theta)[6];
      REAL(ret)[38] = REAL(theta)[10];
      REAL(ret)[46] = REAL(theta)[15];
      REAL(ret)[48] = exp(REAL(theta)[0]);
      REAL(ret)[49] = REAL(theta)[1];
      REAL(ret)[50] = REAL(theta)[3];
      REAL(ret)[51] = REAL(theta)[6];
      REAL(ret)[52] = REAL(theta)[10];
      REAL(ret)[53] = REAL(theta)[15];
      REAL(ret)[54] = 2 * REAL(theta)[21];
      REAL(ret)[55] = REAL(theta)[28];
      REAL(ret)[62] = REAL(theta)[28];
    }
    else if (theta_n == 23){
      REAL(ret)[14] = exp(REAL(theta)[2]);
      REAL(ret)[22] = REAL(theta)[4];
      REAL(ret)[30] = REAL(theta)[7];
      REAL(ret)[38] = REAL(theta)[11];
      REAL(ret)[46] = REAL(theta)[16];
      REAL(ret)[49] = exp(REAL(theta)[2]);
      REAL(ret)[50] = REAL(theta)[4];
      REAL(ret)[51] = REAL(theta)[7];
      REAL(ret)[52] = REAL(theta)[11];
      REAL(ret)[53] = REAL(theta)[16];
      REAL(ret)[54] = 2 * REAL(theta)[22];
      REAL(ret)[55] = REAL(theta)[29];
      REAL(ret)[62] = REAL(theta)[29];
    }
    else if (theta_n == 24){
      REAL(ret)[22] = exp(REAL(theta)[5]);
      REAL(ret)[30] = REAL(theta)[8];
      REAL(ret)[38] = REAL(theta)[12];
      REAL(ret)[46] = REAL(theta)[17];
      REAL(ret)[50] = exp(REAL(theta)[5]);
      REAL(ret)[51] = REAL(theta)[8];
      REAL(ret)[52] = REAL(theta)[12];
      REAL(ret)[53] = REAL(theta)[17];
      REAL(ret)[54] = 2 * REAL(theta)[23];
      REAL(ret)[55] = REAL(theta)[30];
      REAL(ret)[62] = REAL(theta)[30];
    }
    else if (theta_n == 25){
      REAL(ret)[30] = exp(REAL(theta)[9]);
      REAL(ret)[38] = REAL(theta)[13];
      REAL(ret)[46] = REAL(theta)[18];
      REAL(ret)[51] = exp(REAL(theta)[9]);
      REAL(ret)[52] = REAL(theta)[13];
      REAL(ret)[53] = REAL(theta)[18];
      REAL(ret)[54] = 2 * REAL(theta)[24];
      REAL(ret)[55] = REAL(theta)[31];
      REAL(ret)[62] = REAL(theta)[31];
    }
    else if (theta_n == 26){
      REAL(ret)[38] = exp(REAL(theta)[14]);
      REAL(ret)[46] = REAL(theta)[19];
      REAL(ret)[52] = exp(REAL(theta)[14]);
      REAL(ret)[53] = REAL(theta)[19];
      REAL(ret)[54] = 2 * REAL(theta)[25];
      REAL(ret)[55] = REAL(theta)[32];
      REAL(ret)[62] = REAL(theta)[32];
    }
    else if (theta_n == 27){
      REAL(ret)[46] = exp(REAL(theta)[20]);
      REAL(ret)[53] = exp(REAL(theta)[20]);
      REAL(ret)[54] = 2 * REAL(theta)[26];
      REAL(ret)[55] = REAL(theta)[33];
      REAL(ret)[62] = REAL(theta)[33];
    }
    else if (theta_n == 28){
      REAL(ret)[54] = 2 * exp(2 * REAL(theta)[27]);
      REAL(ret)[55] = REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[62] = REAL(theta)[34] * exp(REAL(theta)[27]);
    }
    else if (theta_n == 29){
      REAL(ret)[7] = exp(REAL(theta)[0]);
      REAL(ret)[15] = REAL(theta)[1];
      REAL(ret)[23] = REAL(theta)[3];
      REAL(ret)[31] = REAL(theta)[6];
      REAL(ret)[39] = REAL(theta)[10];
      REAL(ret)[47] = REAL(theta)[15];
      REAL(ret)[55] = REAL(theta)[21];
      REAL(ret)[56] = exp(REAL(theta)[0]);
      REAL(ret)[57] = REAL(theta)[1];
      REAL(ret)[58] = REAL(theta)[3];
      REAL(ret)[59] = REAL(theta)[6];
      REAL(ret)[60] = REAL(theta)[10];
      REAL(ret)[61] = REAL(theta)[15];
      REAL(ret)[62] = REAL(theta)[21];
      REAL(ret)[63] = 2 * REAL(theta)[28];
    }
    else if (theta_n == 30){
      REAL(ret)[15] = exp(REAL(theta)[2]);
      REAL(ret)[23] = REAL(theta)[4];
      REAL(ret)[31] = REAL(theta)[7];
      REAL(ret)[39] = REAL(theta)[11];
      REAL(ret)[47] = REAL(theta)[16];
      REAL(ret)[55] = REAL(theta)[22];
      REAL(ret)[57] = exp(REAL(theta)[2]);
      REAL(ret)[58] = REAL(theta)[4];
      REAL(ret)[59] = REAL(theta)[7];
      REAL(ret)[60] = REAL(theta)[11];
      REAL(ret)[61] = REAL(theta)[16];
      REAL(ret)[62] = REAL(theta)[22];
      REAL(ret)[63] = 2 * REAL(theta)[29];
    }
    else if (theta_n == 31){
      REAL(ret)[23] = exp(REAL(theta)[5]);
      REAL(ret)[31] = REAL(theta)[8];
      REAL(ret)[39] = REAL(theta)[12];
      REAL(ret)[47] = REAL(theta)[17];
      REAL(ret)[55] = REAL(theta)[23];
      REAL(ret)[58] = exp(REAL(theta)[5]);
      REAL(ret)[59] = REAL(theta)[8];
      REAL(ret)[60] = REAL(theta)[12];
      REAL(ret)[61] = REAL(theta)[17];
      REAL(ret)[62] = REAL(theta)[23];
      REAL(ret)[63] = 2 * REAL(theta)[30];
    }
    else if (theta_n == 32){
      REAL(ret)[31] = exp(REAL(theta)[9]);
      REAL(ret)[39] = REAL(theta)[13];
      REAL(ret)[47] = REAL(theta)[18];
      REAL(ret)[55] = REAL(theta)[24];
      REAL(ret)[59] = exp(REAL(theta)[9]);
      REAL(ret)[60] = REAL(theta)[13];
      REAL(ret)[61] = REAL(theta)[18];
      REAL(ret)[62] = REAL(theta)[24];
      REAL(ret)[63] = 2 * REAL(theta)[31];
    }
    else if (theta_n == 33){
      REAL(ret)[39] = exp(REAL(theta)[14]);
      REAL(ret)[47] = REAL(theta)[19];
      REAL(ret)[55] = REAL(theta)[25];
      REAL(ret)[60] = exp(REAL(theta)[14]);
      REAL(ret)[61] = REAL(theta)[19];
      REAL(ret)[62] = REAL(theta)[25];
      REAL(ret)[63] = 2 * REAL(theta)[32];
    }
    else if (theta_n == 34){
      REAL(ret)[47] = exp(REAL(theta)[20]);
      REAL(ret)[55] = REAL(theta)[26];
      REAL(ret)[61] = exp(REAL(theta)[20]);
      REAL(ret)[62] = REAL(theta)[26];
      REAL(ret)[63] = 2 * REAL(theta)[33];
    }
    else if (theta_n == 35){
      REAL(ret)[55] = exp(REAL(theta)[27]);
      REAL(ret)[62] = exp(REAL(theta)[27]);
      REAL(ret)[63] = 2 * REAL(theta)[34];
    }
    else if (theta_n == 36){
      REAL(ret)[63] = 2 * exp(2 * REAL(theta)[35]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 8));for(int i = 0; i < 8; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 2 * exp(2 * REAL(theta)[2]);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 2 * exp(2 * REAL(theta)[5]);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 2 * exp(2 * REAL(theta)[9]);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 2 * exp(2 * REAL(theta)[14]);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 2 * exp(2 * REAL(theta)[20]);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 2 * exp(2 * REAL(theta)[27]);
    }
    else if (theta_n == -38){
      REAL(ret)[7] = 2 * exp(2 * REAL(theta)[35]);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 9){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 9;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -47 || theta_n > 45){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 45){
    error("Requires vector with 45 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 9, 9));for (int i = 0; i < 81; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[1];
      REAL(ret)[10] = exp(REAL(theta)[2]);
      REAL(ret)[18] = REAL(theta)[3];
      REAL(ret)[19] = REAL(theta)[4];
      REAL(ret)[20] = exp(REAL(theta)[5]);
      REAL(ret)[27] = REAL(theta)[6];
      REAL(ret)[28] = REAL(theta)[7];
      REAL(ret)[29] = REAL(theta)[8];
      REAL(ret)[30] = exp(REAL(theta)[9]);
      REAL(ret)[36] = REAL(theta)[10];
      REAL(ret)[37] = REAL(theta)[11];
      REAL(ret)[38] = REAL(theta)[12];
      REAL(ret)[39] = REAL(theta)[13];
      REAL(ret)[40] = exp(REAL(theta)[14]);
      REAL(ret)[45] = REAL(theta)[15];
      REAL(ret)[46] = REAL(theta)[16];
      REAL(ret)[47] = REAL(theta)[17];
      REAL(ret)[48] = REAL(theta)[18];
      REAL(ret)[49] = REAL(theta)[19];
      REAL(ret)[50] = exp(REAL(theta)[20]);
      REAL(ret)[54] = REAL(theta)[21];
      REAL(ret)[55] = REAL(theta)[22];
      REAL(ret)[56] = REAL(theta)[23];
      REAL(ret)[57] = REAL(theta)[24];
      REAL(ret)[58] = REAL(theta)[25];
      REAL(ret)[59] = REAL(theta)[26];
      REAL(ret)[60] = exp(REAL(theta)[27]);
      REAL(ret)[63] = REAL(theta)[28];
      REAL(ret)[64] = REAL(theta)[29];
      REAL(ret)[65] = REAL(theta)[30];
      REAL(ret)[66] = REAL(theta)[31];
      REAL(ret)[67] = REAL(theta)[32];
      REAL(ret)[68] = REAL(theta)[33];
      REAL(ret)[69] = REAL(theta)[34];
      REAL(ret)[70] = exp(REAL(theta)[35]);
      REAL(ret)[72] = REAL(theta)[36];
      REAL(ret)[73] = REAL(theta)[37];
      REAL(ret)[74] = REAL(theta)[38];
      REAL(ret)[75] = REAL(theta)[39];
      REAL(ret)[76] = REAL(theta)[40];
      REAL(ret)[77] = REAL(theta)[41];
      REAL(ret)[78] = REAL(theta)[42];
      REAL(ret)[79] = REAL(theta)[43];
      REAL(ret)[80] = exp(REAL(theta)[44]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[10] = R_pow_di(REAL(theta)[1], 2) + exp(2 * REAL(theta)[2]);
      REAL(ret)[11] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[12] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[1] * REAL(theta)[28] + REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[1] * REAL(theta)[36] + REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[18] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[19] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[20] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + exp(2 * REAL(theta)[5]);
      REAL(ret)[21] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[22] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[23] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[24] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[25] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[26] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[27] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[28] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[29] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[30] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + exp(2 * REAL(theta)[9]);
      REAL(ret)[31] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[32] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[33] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[34] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[35] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[36] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[37] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[38] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[39] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[40] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + exp(2 * REAL(theta)[14]);
      REAL(ret)[41] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[42] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[43] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[44] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[45] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[46] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[47] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[48] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[49] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[50] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + exp(2 * REAL(theta)[20]);
      REAL(ret)[51] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[52] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[53] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[54] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[55] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[56] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[57] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[58] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[59] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[60] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + exp(2 * REAL(theta)[27]);
      REAL(ret)[61] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[62] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[63] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[64] = REAL(theta)[1] * REAL(theta)[28] + REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[65] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[66] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[67] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[68] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[69] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[70] = R_pow_di(REAL(theta)[28], 2) + R_pow_di(REAL(theta)[29], 2) + R_pow_di(REAL(theta)[30], 2) + R_pow_di(REAL(theta)[31], 2) + R_pow_di(REAL(theta)[32], 2) + R_pow_di(REAL(theta)[33], 2) + R_pow_di(REAL(theta)[34], 2) + exp(2 * REAL(theta)[35]);
      REAL(ret)[71] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[72] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[73] = REAL(theta)[1] * REAL(theta)[36] + REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[74] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[75] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[76] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[77] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[78] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[79] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[80] = R_pow_di(REAL(theta)[36], 2) + R_pow_di(REAL(theta)[37], 2) + R_pow_di(REAL(theta)[38], 2) + R_pow_di(REAL(theta)[39], 2) + R_pow_di(REAL(theta)[40], 2) + R_pow_di(REAL(theta)[41], 2) + R_pow_di(REAL(theta)[42], 2) + R_pow_di(REAL(theta)[43], 2) + exp(2 * REAL(theta)[44]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[18] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[27] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[36] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[45] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[54] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[63] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[72] = REAL(theta)[36] * exp(REAL(theta)[0]);
    }
    else if (theta_n == 2){
      REAL(ret)[1] = exp(REAL(theta)[0]);
      REAL(ret)[9] = exp(REAL(theta)[0]);
      REAL(ret)[10] = 2 * REAL(theta)[1];
      REAL(ret)[11] = REAL(theta)[3];
      REAL(ret)[12] = REAL(theta)[6];
      REAL(ret)[13] = REAL(theta)[10];
      REAL(ret)[14] = REAL(theta)[15];
      REAL(ret)[15] = REAL(theta)[21];
      REAL(ret)[16] = REAL(theta)[28];
      REAL(ret)[17] = REAL(theta)[36];
      REAL(ret)[19] = REAL(theta)[3];
      REAL(ret)[28] = REAL(theta)[6];
      REAL(ret)[37] = REAL(theta)[10];
      REAL(ret)[46] = REAL(theta)[15];
      REAL(ret)[55] = REAL(theta)[21];
      REAL(ret)[64] = REAL(theta)[28];
      REAL(ret)[73] = REAL(theta)[36];
    }
    else if (theta_n == 3){
      REAL(ret)[10] = 2 * exp(2 * REAL(theta)[2]);
      REAL(ret)[11] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[12] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[19] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[28] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[37] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[46] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[55] = REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[64] = REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[73] = REAL(theta)[37] * exp(REAL(theta)[2]);
    }
    else if (theta_n == 4){
      REAL(ret)[2] = exp(REAL(theta)[0]);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[18] = exp(REAL(theta)[0]);
      REAL(ret)[19] = REAL(theta)[1];
      REAL(ret)[20] = 2 * REAL(theta)[3];
      REAL(ret)[21] = REAL(theta)[6];
      REAL(ret)[22] = REAL(theta)[10];
      REAL(ret)[23] = REAL(theta)[15];
      REAL(ret)[24] = REAL(theta)[21];
      REAL(ret)[25] = REAL(theta)[28];
      REAL(ret)[26] = REAL(theta)[36];
      REAL(ret)[29] = REAL(theta)[6];
      REAL(ret)[38] = REAL(theta)[10];
      REAL(ret)[47] = REAL(theta)[15];
      REAL(ret)[56] = REAL(theta)[21];
      REAL(ret)[65] = REAL(theta)[28];
      REAL(ret)[74] = REAL(theta)[36];
    }
    else if (theta_n == 5){
      REAL(ret)[11] = exp(REAL(theta)[2]);
      REAL(ret)[19] = exp(REAL(theta)[2]);
      REAL(ret)[20] = 2 * REAL(theta)[4];
      REAL(ret)[21] = REAL(theta)[7];
      REAL(ret)[22] = REAL(theta)[11];
      REAL(ret)[23] = REAL(theta)[16];
      REAL(ret)[24] = REAL(theta)[22];
      REAL(ret)[25] = REAL(theta)[29];
      REAL(ret)[26] = REAL(theta)[37];
      REAL(ret)[29] = REAL(theta)[7];
      REAL(ret)[38] = REAL(theta)[11];
      REAL(ret)[47] = REAL(theta)[16];
      REAL(ret)[56] = REAL(theta)[22];
      REAL(ret)[65] = REAL(theta)[29];
      REAL(ret)[74] = REAL(theta)[37];
    }
    else if (theta_n == 6){
      REAL(ret)[20] = 2 * exp(2 * REAL(theta)[5]);
      REAL(ret)[21] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[22] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[23] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[24] = REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[25] = REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[26] = REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[29] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[38] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[47] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[56] = REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[65] = REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[74] = REAL(theta)[38] * exp(REAL(theta)[5]);
    }
    else if (theta_n == 7){
      REAL(ret)[3] = exp(REAL(theta)[0]);
      REAL(ret)[12] = REAL(theta)[1];
      REAL(ret)[21] = REAL(theta)[3];
      REAL(ret)[27] = exp(REAL(theta)[0]);
      REAL(ret)[28] = REAL(theta)[1];
      REAL(ret)[29] = REAL(theta)[3];
      REAL(ret)[30] = 2 * REAL(theta)[6];
      REAL(ret)[31] = REAL(theta)[10];
      REAL(ret)[32] = REAL(theta)[15];
      REAL(ret)[33] = REAL(theta)[21];
      REAL(ret)[34] = REAL(theta)[28];
      REAL(ret)[35] = REAL(theta)[36];
      REAL(ret)[39] = REAL(theta)[10];
      REAL(ret)[48] = REAL(theta)[15];
      REAL(ret)[57] = REAL(theta)[21];
      REAL(ret)[66] = REAL(theta)[28];
      REAL(ret)[75] = REAL(theta)[36];
    }
    else if (theta_n == 8){
      REAL(ret)[12] = exp(REAL(theta)[2]);
      REAL(ret)[21] = REAL(theta)[4];
      REAL(ret)[28] = exp(REAL(theta)[2]);
      REAL(ret)[29] = REAL(theta)[4];
      REAL(ret)[30] = 2 * REAL(theta)[7];
      REAL(ret)[31] = REAL(theta)[11];
      REAL(ret)[32] = REAL(theta)[16];
      REAL(ret)[33] = REAL(theta)[22];
      REAL(ret)[34] = REAL(theta)[29];
      REAL(ret)[35] = REAL(theta)[37];
      REAL(ret)[39] = REAL(theta)[11];
      REAL(ret)[48] = REAL(theta)[16];
      REAL(ret)[57] = REAL(theta)[22];
      REAL(ret)[66] = REAL(theta)[29];
      REAL(ret)[75] = REAL(theta)[37];
    }
    else if (theta_n == 9){
      REAL(ret)[21] = exp(REAL(theta)[5]);
      REAL(ret)[29] = exp(REAL(theta)[5]);
      REAL(ret)[30] = 2 * REAL(theta)[8];
      REAL(ret)[31] = REAL(theta)[12];
      REAL(ret)[32] = REAL(theta)[17];
      REAL(ret)[33] = REAL(theta)[23];
      REAL(ret)[34] = REAL(theta)[30];
      REAL(ret)[35] = REAL(theta)[38];
      REAL(ret)[39] = REAL(theta)[12];
      REAL(ret)[48] = REAL(theta)[17];
      REAL(ret)[57] = REAL(theta)[23];
      REAL(ret)[66] = REAL(theta)[30];
      REAL(ret)[75] = REAL(theta)[38];
    }
    else if (theta_n == 10){
      REAL(ret)[30] = 2 * exp(2 * REAL(theta)[9]);
      REAL(ret)[31] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[32] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[33] = REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[34] = REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[35] = REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[39] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[48] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[57] = REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[66] = REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[75] = REAL(theta)[39] * exp(REAL(theta)[9]);
    }
    else if (theta_n == 11){
      REAL(ret)[4] = exp(REAL(theta)[0]);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[22] = REAL(theta)[3];
      REAL(ret)[31] = REAL(theta)[6];
      REAL(ret)[36] = exp(REAL(theta)[0]);
      REAL(ret)[37] = REAL(theta)[1];
      REAL(ret)[38] = REAL(theta)[3];
      REAL(ret)[39] = REAL(theta)[6];
      REAL(ret)[40] = 2 * REAL(theta)[10];
      REAL(ret)[41] = REAL(theta)[15];
      REAL(ret)[42] = REAL(theta)[21];
      REAL(ret)[43] = REAL(theta)[28];
      REAL(ret)[44] = REAL(theta)[36];
      REAL(ret)[49] = REAL(theta)[15];
      REAL(ret)[58] = REAL(theta)[21];
      REAL(ret)[67] = REAL(theta)[28];
      REAL(ret)[76] = REAL(theta)[36];
    }
    else if (theta_n == 12){
      REAL(ret)[13] = exp(REAL(theta)[2]);
      REAL(ret)[22] = REAL(theta)[4];
      REAL(ret)[31] = REAL(theta)[7];
      REAL(ret)[37] = exp(REAL(theta)[2]);
      REAL(ret)[38] = REAL(theta)[4];
      REAL(ret)[39] = REAL(theta)[7];
      REAL(ret)[40] = 2 * REAL(theta)[11];
      REAL(ret)[41] = REAL(theta)[16];
      REAL(ret)[42] = REAL(theta)[22];
      REAL(ret)[43] = REAL(theta)[29];
      REAL(ret)[44] = REAL(theta)[37];
      REAL(ret)[49] = REAL(theta)[16];
      REAL(ret)[58] = REAL(theta)[22];
      REAL(ret)[67] = REAL(theta)[29];
      REAL(ret)[76] = REAL(theta)[37];
    }
    else if (theta_n == 13){
      REAL(ret)[22] = exp(REAL(theta)[5]);
      REAL(ret)[31] = REAL(theta)[8];
      REAL(ret)[38] = exp(REAL(theta)[5]);
      REAL(ret)[39] = REAL(theta)[8];
      REAL(ret)[40] = 2 * REAL(theta)[12];
      REAL(ret)[41] = REAL(theta)[17];
      REAL(ret)[42] = REAL(theta)[23];
      REAL(ret)[43] = REAL(theta)[30];
      REAL(ret)[44] = REAL(theta)[38];
      REAL(ret)[49] = REAL(theta)[17];
      REAL(ret)[58] = REAL(theta)[23];
      REAL(ret)[67] = REAL(theta)[30];
      REAL(ret)[76] = REAL(theta)[38];
    }
    else if (theta_n == 14){
      REAL(ret)[31] = exp(REAL(theta)[9]);
      REAL(ret)[39] = exp(REAL(theta)[9]);
      REAL(ret)[40] = 2 * REAL(theta)[13];
      REAL(ret)[41] = REAL(theta)[18];
      REAL(ret)[42] = REAL(theta)[24];
      REAL(ret)[43] = REAL(theta)[31];
      REAL(ret)[44] = REAL(theta)[39];
      REAL(ret)[49] = REAL(theta)[18];
      REAL(ret)[58] = REAL(theta)[24];
      REAL(ret)[67] = REAL(theta)[31];
      REAL(ret)[76] = REAL(theta)[39];
    }
    else if (theta_n == 15){
      REAL(ret)[40] = 2 * exp(2 * REAL(theta)[14]);
      REAL(ret)[41] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[42] = REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[43] = REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[44] = REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[49] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[58] = REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[67] = REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[76] = REAL(theta)[40] * exp(REAL(theta)[14]);
    }
    else if (theta_n == 16){
      REAL(ret)[5] = exp(REAL(theta)[0]);
      REAL(ret)[14] = REAL(theta)[1];
      REAL(ret)[23] = REAL(theta)[3];
      REAL(ret)[32] = REAL(theta)[6];
      REAL(ret)[41] = REAL(theta)[10];
      REAL(ret)[45] = exp(REAL(theta)[0]);
      REAL(ret)[46] = REAL(theta)[1];
      REAL(ret)[47] = REAL(theta)[3];
      REAL(ret)[48] = REAL(theta)[6];
      REAL(ret)[49] = REAL(theta)[10];
      REAL(ret)[50] = 2 * REAL(theta)[15];
      REAL(ret)[51] = REAL(theta)[21];
      REAL(ret)[52] = REAL(theta)[28];
      REAL(ret)[53] = REAL(theta)[36];
      REAL(ret)[59] = REAL(theta)[21];
      REAL(ret)[68] = REAL(theta)[28];
      REAL(ret)[77] = REAL(theta)[36];
    }
    else if (theta_n == 17){
      REAL(ret)[14] = exp(REAL(theta)[2]);
      REAL(ret)[23] = REAL(theta)[4];
      REAL(ret)[32] = REAL(theta)[7];
      REAL(ret)[41] = REAL(theta)[11];
      REAL(ret)[46] = exp(REAL(theta)[2]);
      REAL(ret)[47] = REAL(theta)[4];
      REAL(ret)[48] = REAL(theta)[7];
      REAL(ret)[49] = REAL(theta)[11];
      REAL(ret)[50] = 2 * REAL(theta)[16];
      REAL(ret)[51] = REAL(theta)[22];
      REAL(ret)[52] = REAL(theta)[29];
      REAL(ret)[53] = REAL(theta)[37];
      REAL(ret)[59] = REAL(theta)[22];
      REAL(ret)[68] = REAL(theta)[29];
      REAL(ret)[77] = REAL(theta)[37];
    }
    else if (theta_n == 18){
      REAL(ret)[23] = exp(REAL(theta)[5]);
      REAL(ret)[32] = REAL(theta)[8];
      REAL(ret)[41] = REAL(theta)[12];
      REAL(ret)[47] = exp(REAL(theta)[5]);
      REAL(ret)[48] = REAL(theta)[8];
      REAL(ret)[49] = REAL(theta)[12];
      REAL(ret)[50] = 2 * REAL(theta)[17];
      REAL(ret)[51] = REAL(theta)[23];
      REAL(ret)[52] = REAL(theta)[30];
      REAL(ret)[53] = REAL(theta)[38];
      REAL(ret)[59] = REAL(theta)[23];
      REAL(ret)[68] = REAL(theta)[30];
      REAL(ret)[77] = REAL(theta)[38];
    }
    else if (theta_n == 19){
      REAL(ret)[32] = exp(REAL(theta)[9]);
      REAL(ret)[41] = REAL(theta)[13];
      REAL(ret)[48] = exp(REAL(theta)[9]);
      REAL(ret)[49] = REAL(theta)[13];
      REAL(ret)[50] = 2 * REAL(theta)[18];
      REAL(ret)[51] = REAL(theta)[24];
      REAL(ret)[52] = REAL(theta)[31];
      REAL(ret)[53] = REAL(theta)[39];
      REAL(ret)[59] = REAL(theta)[24];
      REAL(ret)[68] = REAL(theta)[31];
      REAL(ret)[77] = REAL(theta)[39];
    }
    else if (theta_n == 20){
      REAL(ret)[41] = exp(REAL(theta)[14]);
      REAL(ret)[49] = exp(REAL(theta)[14]);
      REAL(ret)[50] = 2 * REAL(theta)[19];
      REAL(ret)[51] = REAL(theta)[25];
      REAL(ret)[52] = REAL(theta)[32];
      REAL(ret)[53] = REAL(theta)[40];
      REAL(ret)[59] = REAL(theta)[25];
      REAL(ret)[68] = REAL(theta)[32];
      REAL(ret)[77] = REAL(theta)[40];
    }
    else if (theta_n == 21){
      REAL(ret)[50] = 2 * exp(2 * REAL(theta)[20]);
      REAL(ret)[51] = REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[52] = REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[53] = REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[59] = REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[68] = REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[77] = REAL(theta)[41] * exp(REAL(theta)[20]);
    }
    else if (theta_n == 22){
      REAL(ret)[6] = exp(REAL(theta)[0]);
      REAL(ret)[15] = REAL(theta)[1];
      REAL(ret)[24] = REAL(theta)[3];
      REAL(ret)[33] = REAL(theta)[6];
      REAL(ret)[42] = REAL(theta)[10];
      REAL(ret)[51] = REAL(theta)[15];
      REAL(ret)[54] = exp(REAL(theta)[0]);
      REAL(ret)[55] = REAL(theta)[1];
      REAL(ret)[56] = REAL(theta)[3];
      REAL(ret)[57] = REAL(theta)[6];
      REAL(ret)[58] = REAL(theta)[10];
      REAL(ret)[59] = REAL(theta)[15];
      REAL(ret)[60] = 2 * REAL(theta)[21];
      REAL(ret)[61] = REAL(theta)[28];
      REAL(ret)[62] = REAL(theta)[36];
      REAL(ret)[69] = REAL(theta)[28];
      REAL(ret)[78] = REAL(theta)[36];
    }
    else if (theta_n == 23){
      REAL(ret)[15] = exp(REAL(theta)[2]);
      REAL(ret)[24] = REAL(theta)[4];
      REAL(ret)[33] = REAL(theta)[7];
      REAL(ret)[42] = REAL(theta)[11];
      REAL(ret)[51] = REAL(theta)[16];
      REAL(ret)[55] = exp(REAL(theta)[2]);
      REAL(ret)[56] = REAL(theta)[4];
      REAL(ret)[57] = REAL(theta)[7];
      REAL(ret)[58] = REAL(theta)[11];
      REAL(ret)[59] = REAL(theta)[16];
      REAL(ret)[60] = 2 * REAL(theta)[22];
      REAL(ret)[61] = REAL(theta)[29];
      REAL(ret)[62] = REAL(theta)[37];
      REAL(ret)[69] = REAL(theta)[29];
      REAL(ret)[78] = REAL(theta)[37];
    }
    else if (theta_n == 24){
      REAL(ret)[24] = exp(REAL(theta)[5]);
      REAL(ret)[33] = REAL(theta)[8];
      REAL(ret)[42] = REAL(theta)[12];
      REAL(ret)[51] = REAL(theta)[17];
      REAL(ret)[56] = exp(REAL(theta)[5]);
      REAL(ret)[57] = REAL(theta)[8];
      REAL(ret)[58] = REAL(theta)[12];
      REAL(ret)[59] = REAL(theta)[17];
      REAL(ret)[60] = 2 * REAL(theta)[23];
      REAL(ret)[61] = REAL(theta)[30];
      REAL(ret)[62] = REAL(theta)[38];
      REAL(ret)[69] = REAL(theta)[30];
      REAL(ret)[78] = REAL(theta)[38];
    }
    else if (theta_n == 25){
      REAL(ret)[33] = exp(REAL(theta)[9]);
      REAL(ret)[42] = REAL(theta)[13];
      REAL(ret)[51] = REAL(theta)[18];
      REAL(ret)[57] = exp(REAL(theta)[9]);
      REAL(ret)[58] = REAL(theta)[13];
      REAL(ret)[59] = REAL(theta)[18];
      REAL(ret)[60] = 2 * REAL(theta)[24];
      REAL(ret)[61] = REAL(theta)[31];
      REAL(ret)[62] = REAL(theta)[39];
      REAL(ret)[69] = REAL(theta)[31];
      REAL(ret)[78] = REAL(theta)[39];
    }
    else if (theta_n == 26){
      REAL(ret)[42] = exp(REAL(theta)[14]);
      REAL(ret)[51] = REAL(theta)[19];
      REAL(ret)[58] = exp(REAL(theta)[14]);
      REAL(ret)[59] = REAL(theta)[19];
      REAL(ret)[60] = 2 * REAL(theta)[25];
      REAL(ret)[61] = REAL(theta)[32];
      REAL(ret)[62] = REAL(theta)[40];
      REAL(ret)[69] = REAL(theta)[32];
      REAL(ret)[78] = REAL(theta)[40];
    }
    else if (theta_n == 27){
      REAL(ret)[51] = exp(REAL(theta)[20]);
      REAL(ret)[59] = exp(REAL(theta)[20]);
      REAL(ret)[60] = 2 * REAL(theta)[26];
      REAL(ret)[61] = REAL(theta)[33];
      REAL(ret)[62] = REAL(theta)[41];
      REAL(ret)[69] = REAL(theta)[33];
      REAL(ret)[78] = REAL(theta)[41];
    }
    else if (theta_n == 28){
      REAL(ret)[60] = 2 * exp(2 * REAL(theta)[27]);
      REAL(ret)[61] = REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[62] = REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[69] = REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[78] = REAL(theta)[42] * exp(REAL(theta)[27]);
    }
    else if (theta_n == 29){
      REAL(ret)[7] = exp(REAL(theta)[0]);
      REAL(ret)[16] = REAL(theta)[1];
      REAL(ret)[25] = REAL(theta)[3];
      REAL(ret)[34] = REAL(theta)[6];
      REAL(ret)[43] = REAL(theta)[10];
      REAL(ret)[52] = REAL(theta)[15];
      REAL(ret)[61] = REAL(theta)[21];
      REAL(ret)[63] = exp(REAL(theta)[0]);
      REAL(ret)[64] = REAL(theta)[1];
      REAL(ret)[65] = REAL(theta)[3];
      REAL(ret)[66] = REAL(theta)[6];
      REAL(ret)[67] = REAL(theta)[10];
      REAL(ret)[68] = REAL(theta)[15];
      REAL(ret)[69] = REAL(theta)[21];
      REAL(ret)[70] = 2 * REAL(theta)[28];
      REAL(ret)[71] = REAL(theta)[36];
      REAL(ret)[79] = REAL(theta)[36];
    }
    else if (theta_n == 30){
      REAL(ret)[16] = exp(REAL(theta)[2]);
      REAL(ret)[25] = REAL(theta)[4];
      REAL(ret)[34] = REAL(theta)[7];
      REAL(ret)[43] = REAL(theta)[11];
      REAL(ret)[52] = REAL(theta)[16];
      REAL(ret)[61] = REAL(theta)[22];
      REAL(ret)[64] = exp(REAL(theta)[2]);
      REAL(ret)[65] = REAL(theta)[4];
      REAL(ret)[66] = REAL(theta)[7];
      REAL(ret)[67] = REAL(theta)[11];
      REAL(ret)[68] = REAL(theta)[16];
      REAL(ret)[69] = REAL(theta)[22];
      REAL(ret)[70] = 2 * REAL(theta)[29];
      REAL(ret)[71] = REAL(theta)[37];
      REAL(ret)[79] = REAL(theta)[37];
    }
    else if (theta_n == 31){
      REAL(ret)[25] = exp(REAL(theta)[5]);
      REAL(ret)[34] = REAL(theta)[8];
      REAL(ret)[43] = REAL(theta)[12];
      REAL(ret)[52] = REAL(theta)[17];
      REAL(ret)[61] = REAL(theta)[23];
      REAL(ret)[65] = exp(REAL(theta)[5]);
      REAL(ret)[66] = REAL(theta)[8];
      REAL(ret)[67] = REAL(theta)[12];
      REAL(ret)[68] = REAL(theta)[17];
      REAL(ret)[69] = REAL(theta)[23];
      REAL(ret)[70] = 2 * REAL(theta)[30];
      REAL(ret)[71] = REAL(theta)[38];
      REAL(ret)[79] = REAL(theta)[38];
    }
    else if (theta_n == 32){
      REAL(ret)[34] = exp(REAL(theta)[9]);
      REAL(ret)[43] = REAL(theta)[13];
      REAL(ret)[52] = REAL(theta)[18];
      REAL(ret)[61] = REAL(theta)[24];
      REAL(ret)[66] = exp(REAL(theta)[9]);
      REAL(ret)[67] = REAL(theta)[13];
      REAL(ret)[68] = REAL(theta)[18];
      REAL(ret)[69] = REAL(theta)[24];
      REAL(ret)[70] = 2 * REAL(theta)[31];
      REAL(ret)[71] = REAL(theta)[39];
      REAL(ret)[79] = REAL(theta)[39];
    }
    else if (theta_n == 33){
      REAL(ret)[43] = exp(REAL(theta)[14]);
      REAL(ret)[52] = REAL(theta)[19];
      REAL(ret)[61] = REAL(theta)[25];
      REAL(ret)[67] = exp(REAL(theta)[14]);
      REAL(ret)[68] = REAL(theta)[19];
      REAL(ret)[69] = REAL(theta)[25];
      REAL(ret)[70] = 2 * REAL(theta)[32];
      REAL(ret)[71] = REAL(theta)[40];
      REAL(ret)[79] = REAL(theta)[40];
    }
    else if (theta_n == 34){
      REAL(ret)[52] = exp(REAL(theta)[20]);
      REAL(ret)[61] = REAL(theta)[26];
      REAL(ret)[68] = exp(REAL(theta)[20]);
      REAL(ret)[69] = REAL(theta)[26];
      REAL(ret)[70] = 2 * REAL(theta)[33];
      REAL(ret)[71] = REAL(theta)[41];
      REAL(ret)[79] = REAL(theta)[41];
    }
    else if (theta_n == 35){
      REAL(ret)[61] = exp(REAL(theta)[27]);
      REAL(ret)[69] = exp(REAL(theta)[27]);
      REAL(ret)[70] = 2 * REAL(theta)[34];
      REAL(ret)[71] = REAL(theta)[42];
      REAL(ret)[79] = REAL(theta)[42];
    }
    else if (theta_n == 36){
      REAL(ret)[70] = 2 * exp(2 * REAL(theta)[35]);
      REAL(ret)[71] = REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[79] = REAL(theta)[43] * exp(REAL(theta)[35]);
    }
    else if (theta_n == 37){
      REAL(ret)[8] = exp(REAL(theta)[0]);
      REAL(ret)[17] = REAL(theta)[1];
      REAL(ret)[26] = REAL(theta)[3];
      REAL(ret)[35] = REAL(theta)[6];
      REAL(ret)[44] = REAL(theta)[10];
      REAL(ret)[53] = REAL(theta)[15];
      REAL(ret)[62] = REAL(theta)[21];
      REAL(ret)[71] = REAL(theta)[28];
      REAL(ret)[72] = exp(REAL(theta)[0]);
      REAL(ret)[73] = REAL(theta)[1];
      REAL(ret)[74] = REAL(theta)[3];
      REAL(ret)[75] = REAL(theta)[6];
      REAL(ret)[76] = REAL(theta)[10];
      REAL(ret)[77] = REAL(theta)[15];
      REAL(ret)[78] = REAL(theta)[21];
      REAL(ret)[79] = REAL(theta)[28];
      REAL(ret)[80] = 2 * REAL(theta)[36];
    }
    else if (theta_n == 38){
      REAL(ret)[17] = exp(REAL(theta)[2]);
      REAL(ret)[26] = REAL(theta)[4];
      REAL(ret)[35] = REAL(theta)[7];
      REAL(ret)[44] = REAL(theta)[11];
      REAL(ret)[53] = REAL(theta)[16];
      REAL(ret)[62] = REAL(theta)[22];
      REAL(ret)[71] = REAL(theta)[29];
      REAL(ret)[73] = exp(REAL(theta)[2]);
      REAL(ret)[74] = REAL(theta)[4];
      REAL(ret)[75] = REAL(theta)[7];
      REAL(ret)[76] = REAL(theta)[11];
      REAL(ret)[77] = REAL(theta)[16];
      REAL(ret)[78] = REAL(theta)[22];
      REAL(ret)[79] = REAL(theta)[29];
      REAL(ret)[80] = 2 * REAL(theta)[37];
    }
    else if (theta_n == 39){
      REAL(ret)[26] = exp(REAL(theta)[5]);
      REAL(ret)[35] = REAL(theta)[8];
      REAL(ret)[44] = REAL(theta)[12];
      REAL(ret)[53] = REAL(theta)[17];
      REAL(ret)[62] = REAL(theta)[23];
      REAL(ret)[71] = REAL(theta)[30];
      REAL(ret)[74] = exp(REAL(theta)[5]);
      REAL(ret)[75] = REAL(theta)[8];
      REAL(ret)[76] = REAL(theta)[12];
      REAL(ret)[77] = REAL(theta)[17];
      REAL(ret)[78] = REAL(theta)[23];
      REAL(ret)[79] = REAL(theta)[30];
      REAL(ret)[80] = 2 * REAL(theta)[38];
    }
    else if (theta_n == 40){
      REAL(ret)[35] = exp(REAL(theta)[9]);
      REAL(ret)[44] = REAL(theta)[13];
      REAL(ret)[53] = REAL(theta)[18];
      REAL(ret)[62] = REAL(theta)[24];
      REAL(ret)[71] = REAL(theta)[31];
      REAL(ret)[75] = exp(REAL(theta)[9]);
      REAL(ret)[76] = REAL(theta)[13];
      REAL(ret)[77] = REAL(theta)[18];
      REAL(ret)[78] = REAL(theta)[24];
      REAL(ret)[79] = REAL(theta)[31];
      REAL(ret)[80] = 2 * REAL(theta)[39];
    }
    else if (theta_n == 41){
      REAL(ret)[44] = exp(REAL(theta)[14]);
      REAL(ret)[53] = REAL(theta)[19];
      REAL(ret)[62] = REAL(theta)[25];
      REAL(ret)[71] = REAL(theta)[32];
      REAL(ret)[76] = exp(REAL(theta)[14]);
      REAL(ret)[77] = REAL(theta)[19];
      REAL(ret)[78] = REAL(theta)[25];
      REAL(ret)[79] = REAL(theta)[32];
      REAL(ret)[80] = 2 * REAL(theta)[40];
    }
    else if (theta_n == 42){
      REAL(ret)[53] = exp(REAL(theta)[20]);
      REAL(ret)[62] = REAL(theta)[26];
      REAL(ret)[71] = REAL(theta)[33];
      REAL(ret)[77] = exp(REAL(theta)[20]);
      REAL(ret)[78] = REAL(theta)[26];
      REAL(ret)[79] = REAL(theta)[33];
      REAL(ret)[80] = 2 * REAL(theta)[41];
    }
    else if (theta_n == 43){
      REAL(ret)[62] = exp(REAL(theta)[27]);
      REAL(ret)[71] = REAL(theta)[34];
      REAL(ret)[78] = exp(REAL(theta)[27]);
      REAL(ret)[79] = REAL(theta)[34];
      REAL(ret)[80] = 2 * REAL(theta)[42];
    }
    else if (theta_n == 44){
      REAL(ret)[71] = exp(REAL(theta)[35]);
      REAL(ret)[79] = exp(REAL(theta)[35]);
      REAL(ret)[80] = 2 * REAL(theta)[43];
    }
    else if (theta_n == 45){
      REAL(ret)[80] = 2 * exp(2 * REAL(theta)[44]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 9));for(int i = 0; i < 9; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 2 * exp(2 * REAL(theta)[2]);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 2 * exp(2 * REAL(theta)[5]);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 2 * exp(2 * REAL(theta)[9]);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 2 * exp(2 * REAL(theta)[14]);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 2 * exp(2 * REAL(theta)[20]);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 2 * exp(2 * REAL(theta)[27]);
    }
    else if (theta_n == -38){
      REAL(ret)[7] = 2 * exp(2 * REAL(theta)[35]);
    }
    else if (theta_n == -47){
      REAL(ret)[8] = 2 * exp(2 * REAL(theta)[44]);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 10){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 10;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -57 || theta_n > 55){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 55){
    error("Requires vector with 55 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 10, 10));for (int i = 0; i < 100; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
      REAL(ret)[10] = REAL(theta)[1];
      REAL(ret)[11] = exp(REAL(theta)[2]);
      REAL(ret)[20] = REAL(theta)[3];
      REAL(ret)[21] = REAL(theta)[4];
      REAL(ret)[22] = exp(REAL(theta)[5]);
      REAL(ret)[30] = REAL(theta)[6];
      REAL(ret)[31] = REAL(theta)[7];
      REAL(ret)[32] = REAL(theta)[8];
      REAL(ret)[33] = exp(REAL(theta)[9]);
      REAL(ret)[40] = REAL(theta)[10];
      REAL(ret)[41] = REAL(theta)[11];
      REAL(ret)[42] = REAL(theta)[12];
      REAL(ret)[43] = REAL(theta)[13];
      REAL(ret)[44] = exp(REAL(theta)[14]);
      REAL(ret)[50] = REAL(theta)[15];
      REAL(ret)[51] = REAL(theta)[16];
      REAL(ret)[52] = REAL(theta)[17];
      REAL(ret)[53] = REAL(theta)[18];
      REAL(ret)[54] = REAL(theta)[19];
      REAL(ret)[55] = exp(REAL(theta)[20]);
      REAL(ret)[60] = REAL(theta)[21];
      REAL(ret)[61] = REAL(theta)[22];
      REAL(ret)[62] = REAL(theta)[23];
      REAL(ret)[63] = REAL(theta)[24];
      REAL(ret)[64] = REAL(theta)[25];
      REAL(ret)[65] = REAL(theta)[26];
      REAL(ret)[66] = exp(REAL(theta)[27]);
      REAL(ret)[70] = REAL(theta)[28];
      REAL(ret)[71] = REAL(theta)[29];
      REAL(ret)[72] = REAL(theta)[30];
      REAL(ret)[73] = REAL(theta)[31];
      REAL(ret)[74] = REAL(theta)[32];
      REAL(ret)[75] = REAL(theta)[33];
      REAL(ret)[76] = REAL(theta)[34];
      REAL(ret)[77] = exp(REAL(theta)[35]);
      REAL(ret)[80] = REAL(theta)[36];
      REAL(ret)[81] = REAL(theta)[37];
      REAL(ret)[82] = REAL(theta)[38];
      REAL(ret)[83] = REAL(theta)[39];
      REAL(ret)[84] = REAL(theta)[40];
      REAL(ret)[85] = REAL(theta)[41];
      REAL(ret)[86] = REAL(theta)[42];
      REAL(ret)[87] = REAL(theta)[43];
      REAL(ret)[88] = exp(REAL(theta)[44]);
      REAL(ret)[90] = REAL(theta)[45];
      REAL(ret)[91] = REAL(theta)[46];
      REAL(ret)[92] = REAL(theta)[47];
      REAL(ret)[93] = REAL(theta)[48];
      REAL(ret)[94] = REAL(theta)[49];
      REAL(ret)[95] = REAL(theta)[50];
      REAL(ret)[96] = REAL(theta)[51];
      REAL(ret)[97] = REAL(theta)[52];
      REAL(ret)[98] = REAL(theta)[53];
      REAL(ret)[99] = exp(REAL(theta)[54]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[45] * exp(REAL(theta)[0]);
      REAL(ret)[10] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[11] = R_pow_di(REAL(theta)[1], 2) + exp(2 * REAL(theta)[2]);
      REAL(ret)[12] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[1] * REAL(theta)[28] + REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[18] = REAL(theta)[1] * REAL(theta)[36] + REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[19] = REAL(theta)[1] * REAL(theta)[45] + REAL(theta)[46] * exp(REAL(theta)[2]);
      REAL(ret)[20] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[21] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[22] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + exp(2 * REAL(theta)[5]);
      REAL(ret)[23] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[24] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[25] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[26] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[27] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[28] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[29] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * exp(REAL(theta)[5]);
      REAL(ret)[30] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[31] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[32] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[33] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + exp(2 * REAL(theta)[9]);
      REAL(ret)[34] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[35] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[36] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[37] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[38] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[39] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * exp(REAL(theta)[9]);
      REAL(ret)[40] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[41] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[42] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[43] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[44] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + exp(2 * REAL(theta)[14]);
      REAL(ret)[45] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[46] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[47] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[48] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[49] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + REAL(theta)[49] * exp(REAL(theta)[14]);
      REAL(ret)[50] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[51] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[52] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[53] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[54] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[55] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + exp(2 * REAL(theta)[20]);
      REAL(ret)[56] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[57] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[58] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[59] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + REAL(theta)[50] * exp(REAL(theta)[20]);
      REAL(ret)[60] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[61] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[62] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[63] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[64] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[65] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[66] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + exp(2 * REAL(theta)[27]);
      REAL(ret)[67] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[68] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[69] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + REAL(theta)[51] * exp(REAL(theta)[27]);
      REAL(ret)[70] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[71] = REAL(theta)[1] * REAL(theta)[28] + REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[72] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[73] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[74] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[75] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[76] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[77] = R_pow_di(REAL(theta)[28], 2) + R_pow_di(REAL(theta)[29], 2) + R_pow_di(REAL(theta)[30], 2) + R_pow_di(REAL(theta)[31], 2) + R_pow_di(REAL(theta)[32], 2) + R_pow_di(REAL(theta)[33], 2) + R_pow_di(REAL(theta)[34], 2) + exp(2 * REAL(theta)[35]);
      REAL(ret)[78] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[79] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + REAL(theta)[52] * exp(REAL(theta)[35]);
      REAL(ret)[80] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[81] = REAL(theta)[1] * REAL(theta)[36] + REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[82] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[83] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[84] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[85] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[86] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[87] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[88] = R_pow_di(REAL(theta)[36], 2) + R_pow_di(REAL(theta)[37], 2) + R_pow_di(REAL(theta)[38], 2) + R_pow_di(REAL(theta)[39], 2) + R_pow_di(REAL(theta)[40], 2) + R_pow_di(REAL(theta)[41], 2) + R_pow_di(REAL(theta)[42], 2) + R_pow_di(REAL(theta)[43], 2) + exp(2 * REAL(theta)[44]);
      REAL(ret)[89] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + REAL(theta)[53] * exp(REAL(theta)[44]);
      REAL(ret)[90] = REAL(theta)[45] * exp(REAL(theta)[0]);
      REAL(ret)[91] = REAL(theta)[1] * REAL(theta)[45] + REAL(theta)[46] * exp(REAL(theta)[2]);
      REAL(ret)[92] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * exp(REAL(theta)[5]);
      REAL(ret)[93] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * exp(REAL(theta)[9]);
      REAL(ret)[94] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + REAL(theta)[49] * exp(REAL(theta)[14]);
      REAL(ret)[95] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + REAL(theta)[50] * exp(REAL(theta)[20]);
      REAL(ret)[96] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + REAL(theta)[51] * exp(REAL(theta)[27]);
      REAL(ret)[97] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + REAL(theta)[52] * exp(REAL(theta)[35]);
      REAL(ret)[98] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + REAL(theta)[53] * exp(REAL(theta)[44]);
      REAL(ret)[99] = R_pow_di(REAL(theta)[45], 2) + R_pow_di(REAL(theta)[46], 2) + R_pow_di(REAL(theta)[47], 2) + R_pow_di(REAL(theta)[48], 2) + R_pow_di(REAL(theta)[49], 2) + R_pow_di(REAL(theta)[50], 2) + R_pow_di(REAL(theta)[51], 2) + R_pow_di(REAL(theta)[52], 2) + R_pow_di(REAL(theta)[53], 2) + exp(2 * REAL(theta)[54]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[45] * exp(REAL(theta)[0]);
      REAL(ret)[10] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[20] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[30] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[40] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[50] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[60] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[70] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[80] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[90] = REAL(theta)[45] * exp(REAL(theta)[0]);
    }
    else if (theta_n == 2){
      REAL(ret)[1] = exp(REAL(theta)[0]);
      REAL(ret)[10] = exp(REAL(theta)[0]);
      REAL(ret)[11] = 2 * REAL(theta)[1];
      REAL(ret)[12] = REAL(theta)[3];
      REAL(ret)[13] = REAL(theta)[6];
      REAL(ret)[14] = REAL(theta)[10];
      REAL(ret)[15] = REAL(theta)[15];
      REAL(ret)[16] = REAL(theta)[21];
      REAL(ret)[17] = REAL(theta)[28];
      REAL(ret)[18] = REAL(theta)[36];
      REAL(ret)[19] = REAL(theta)[45];
      REAL(ret)[21] = REAL(theta)[3];
      REAL(ret)[31] = REAL(theta)[6];
      REAL(ret)[41] = REAL(theta)[10];
      REAL(ret)[51] = REAL(theta)[15];
      REAL(ret)[61] = REAL(theta)[21];
      REAL(ret)[71] = REAL(theta)[28];
      REAL(ret)[81] = REAL(theta)[36];
      REAL(ret)[91] = REAL(theta)[45];
    }
    else if (theta_n == 3){
      REAL(ret)[11] = 2 * exp(2 * REAL(theta)[2]);
      REAL(ret)[12] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[18] = REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[19] = REAL(theta)[46] * exp(REAL(theta)[2]);
      REAL(ret)[21] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[31] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[41] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[51] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[61] = REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[71] = REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[81] = REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[91] = REAL(theta)[46] * exp(REAL(theta)[2]);
    }
    else if (theta_n == 4){
      REAL(ret)[2] = exp(REAL(theta)[0]);
      REAL(ret)[12] = REAL(theta)[1];
      REAL(ret)[20] = exp(REAL(theta)[0]);
      REAL(ret)[21] = REAL(theta)[1];
      REAL(ret)[22] = 2 * REAL(theta)[3];
      REAL(ret)[23] = REAL(theta)[6];
      REAL(ret)[24] = REAL(theta)[10];
      REAL(ret)[25] = REAL(theta)[15];
      REAL(ret)[26] = REAL(theta)[21];
      REAL(ret)[27] = REAL(theta)[28];
      REAL(ret)[28] = REAL(theta)[36];
      REAL(ret)[29] = REAL(theta)[45];
      REAL(ret)[32] = REAL(theta)[6];
      REAL(ret)[42] = REAL(theta)[10];
      REAL(ret)[52] = REAL(theta)[15];
      REAL(ret)[62] = REAL(theta)[21];
      REAL(ret)[72] = REAL(theta)[28];
      REAL(ret)[82] = REAL(theta)[36];
      REAL(ret)[92] = REAL(theta)[45];
    }
    else if (theta_n == 5){
      REAL(ret)[12] = exp(REAL(theta)[2]);
      REAL(ret)[21] = exp(REAL(theta)[2]);
      REAL(ret)[22] = 2 * REAL(theta)[4];
      REAL(ret)[23] = REAL(theta)[7];
      REAL(ret)[24] = REAL(theta)[11];
      REAL(ret)[25] = REAL(theta)[16];
      REAL(ret)[26] = REAL(theta)[22];
      REAL(ret)[27] = REAL(theta)[29];
      REAL(ret)[28] = REAL(theta)[37];
      REAL(ret)[29] = REAL(theta)[46];
      REAL(ret)[32] = REAL(theta)[7];
      REAL(ret)[42] = REAL(theta)[11];
      REAL(ret)[52] = REAL(theta)[16];
      REAL(ret)[62] = REAL(theta)[22];
      REAL(ret)[72] = REAL(theta)[29];
      REAL(ret)[82] = REAL(theta)[37];
      REAL(ret)[92] = REAL(theta)[46];
    }
    else if (theta_n == 6){
      REAL(ret)[22] = 2 * exp(2 * REAL(theta)[5]);
      REAL(ret)[23] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[24] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[25] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[26] = REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[27] = REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[28] = REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[29] = REAL(theta)[47] * exp(REAL(theta)[5]);
      REAL(ret)[32] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[42] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[52] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[62] = REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[72] = REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[82] = REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[92] = REAL(theta)[47] * exp(REAL(theta)[5]);
    }
    else if (theta_n == 7){
      REAL(ret)[3] = exp(REAL(theta)[0]);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[23] = REAL(theta)[3];
      REAL(ret)[30] = exp(REAL(theta)[0]);
      REAL(ret)[31] = REAL(theta)[1];
      REAL(ret)[32] = REAL(theta)[3];
      REAL(ret)[33] = 2 * REAL(theta)[6];
      REAL(ret)[34] = REAL(theta)[10];
      REAL(ret)[35] = REAL(theta)[15];
      REAL(ret)[36] = REAL(theta)[21];
      REAL(ret)[37] = REAL(theta)[28];
      REAL(ret)[38] = REAL(theta)[36];
      REAL(ret)[39] = REAL(theta)[45];
      REAL(ret)[43] = REAL(theta)[10];
      REAL(ret)[53] = REAL(theta)[15];
      REAL(ret)[63] = REAL(theta)[21];
      REAL(ret)[73] = REAL(theta)[28];
      REAL(ret)[83] = REAL(theta)[36];
      REAL(ret)[93] = REAL(theta)[45];
    }
    else if (theta_n == 8){
      REAL(ret)[13] = exp(REAL(theta)[2]);
      REAL(ret)[23] = REAL(theta)[4];
      REAL(ret)[31] = exp(REAL(theta)[2]);
      REAL(ret)[32] = REAL(theta)[4];
      REAL(ret)[33] = 2 * REAL(theta)[7];
      REAL(ret)[34] = REAL(theta)[11];
      REAL(ret)[35] = REAL(theta)[16];
      REAL(ret)[36] = REAL(theta)[22];
      REAL(ret)[37] = REAL(theta)[29];
      REAL(ret)[38] = REAL(theta)[37];
      REAL(ret)[39] = REAL(theta)[46];
      REAL(ret)[43] = REAL(theta)[11];
      REAL(ret)[53] = REAL(theta)[16];
      REAL(ret)[63] = REAL(theta)[22];
      REAL(ret)[73] = REAL(theta)[29];
      REAL(ret)[83] = REAL(theta)[37];
      REAL(ret)[93] = REAL(theta)[46];
    }
    else if (theta_n == 9){
      REAL(ret)[23] = exp(REAL(theta)[5]);
      REAL(ret)[32] = exp(REAL(theta)[5]);
      REAL(ret)[33] = 2 * REAL(theta)[8];
      REAL(ret)[34] = REAL(theta)[12];
      REAL(ret)[35] = REAL(theta)[17];
      REAL(ret)[36] = REAL(theta)[23];
      REAL(ret)[37] = REAL(theta)[30];
      REAL(ret)[38] = REAL(theta)[38];
      REAL(ret)[39] = REAL(theta)[47];
      REAL(ret)[43] = REAL(theta)[12];
      REAL(ret)[53] = REAL(theta)[17];
      REAL(ret)[63] = REAL(theta)[23];
      REAL(ret)[73] = REAL(theta)[30];
      REAL(ret)[83] = REAL(theta)[38];
      REAL(ret)[93] = REAL(theta)[47];
    }
    else if (theta_n == 10){
      REAL(ret)[33] = 2 * exp(2 * REAL(theta)[9]);
      REAL(ret)[34] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[35] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[36] = REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[37] = REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[38] = REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[39] = REAL(theta)[48] * exp(REAL(theta)[9]);
      REAL(ret)[43] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[53] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[63] = REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[73] = REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[83] = REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[93] = REAL(theta)[48] * exp(REAL(theta)[9]);
    }
    else if (theta_n == 11){
      REAL(ret)[4] = exp(REAL(theta)[0]);
      REAL(ret)[14] = REAL(theta)[1];
      REAL(ret)[24] = REAL(theta)[3];
      REAL(ret)[34] = REAL(theta)[6];
      REAL(ret)[40] = exp(REAL(theta)[0]);
      REAL(ret)[41] = REAL(theta)[1];
      REAL(ret)[42] = REAL(theta)[3];
      REAL(ret)[43] = REAL(theta)[6];
      REAL(ret)[44] = 2 * REAL(theta)[10];
      REAL(ret)[45] = REAL(theta)[15];
      REAL(ret)[46] = REAL(theta)[21];
      REAL(ret)[47] = REAL(theta)[28];
      REAL(ret)[48] = REAL(theta)[36];
      REAL(ret)[49] = REAL(theta)[45];
      REAL(ret)[54] = REAL(theta)[15];
      REAL(ret)[64] = REAL(theta)[21];
      REAL(ret)[74] = REAL(theta)[28];
      REAL(ret)[84] = REAL(theta)[36];
      REAL(ret)[94] = REAL(theta)[45];
    }
    else if (theta_n == 12){
      REAL(ret)[14] = exp(REAL(theta)[2]);
      REAL(ret)[24] = REAL(theta)[4];
      REAL(ret)[34] = REAL(theta)[7];
      REAL(ret)[41] = exp(REAL(theta)[2]);
      REAL(ret)[42] = REAL(theta)[4];
      REAL(ret)[43] = REAL(theta)[7];
      REAL(ret)[44] = 2 * REAL(theta)[11];
      REAL(ret)[45] = REAL(theta)[16];
      REAL(ret)[46] = REAL(theta)[22];
      REAL(ret)[47] = REAL(theta)[29];
      REAL(ret)[48] = REAL(theta)[37];
      REAL(ret)[49] = REAL(theta)[46];
      REAL(ret)[54] = REAL(theta)[16];
      REAL(ret)[64] = REAL(theta)[22];
      REAL(ret)[74] = REAL(theta)[29];
      REAL(ret)[84] = REAL(theta)[37];
      REAL(ret)[94] = REAL(theta)[46];
    }
    else if (theta_n == 13){
      REAL(ret)[24] = exp(REAL(theta)[5]);
      REAL(ret)[34] = REAL(theta)[8];
      REAL(ret)[42] = exp(REAL(theta)[5]);
      REAL(ret)[43] = REAL(theta)[8];
      REAL(ret)[44] = 2 * REAL(theta)[12];
      REAL(ret)[45] = REAL(theta)[17];
      REAL(ret)[46] = REAL(theta)[23];
      REAL(ret)[47] = REAL(theta)[30];
      REAL(ret)[48] = REAL(theta)[38];
      REAL(ret)[49] = REAL(theta)[47];
      REAL(ret)[54] = REAL(theta)[17];
      REAL(ret)[64] = REAL(theta)[23];
      REAL(ret)[74] = REAL(theta)[30];
      REAL(ret)[84] = REAL(theta)[38];
      REAL(ret)[94] = REAL(theta)[47];
    }
    else if (theta_n == 14){
      REAL(ret)[34] = exp(REAL(theta)[9]);
      REAL(ret)[43] = exp(REAL(theta)[9]);
      REAL(ret)[44] = 2 * REAL(theta)[13];
      REAL(ret)[45] = REAL(theta)[18];
      REAL(ret)[46] = REAL(theta)[24];
      REAL(ret)[47] = REAL(theta)[31];
      REAL(ret)[48] = REAL(theta)[39];
      REAL(ret)[49] = REAL(theta)[48];
      REAL(ret)[54] = REAL(theta)[18];
      REAL(ret)[64] = REAL(theta)[24];
      REAL(ret)[74] = REAL(theta)[31];
      REAL(ret)[84] = REAL(theta)[39];
      REAL(ret)[94] = REAL(theta)[48];
    }
    else if (theta_n == 15){
      REAL(ret)[44] = 2 * exp(2 * REAL(theta)[14]);
      REAL(ret)[45] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[46] = REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[47] = REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[48] = REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[49] = REAL(theta)[49] * exp(REAL(theta)[14]);
      REAL(ret)[54] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[64] = REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[74] = REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[84] = REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[94] = REAL(theta)[49] * exp(REAL(theta)[14]);
    }
    else if (theta_n == 16){
      REAL(ret)[5] = exp(REAL(theta)[0]);
      REAL(ret)[15] = REAL(theta)[1];
      REAL(ret)[25] = REAL(theta)[3];
      REAL(ret)[35] = REAL(theta)[6];
      REAL(ret)[45] = REAL(theta)[10];
      REAL(ret)[50] = exp(REAL(theta)[0]);
      REAL(ret)[51] = REAL(theta)[1];
      REAL(ret)[52] = REAL(theta)[3];
      REAL(ret)[53] = REAL(theta)[6];
      REAL(ret)[54] = REAL(theta)[10];
      REAL(ret)[55] = 2 * REAL(theta)[15];
      REAL(ret)[56] = REAL(theta)[21];
      REAL(ret)[57] = REAL(theta)[28];
      REAL(ret)[58] = REAL(theta)[36];
      REAL(ret)[59] = REAL(theta)[45];
      REAL(ret)[65] = REAL(theta)[21];
      REAL(ret)[75] = REAL(theta)[28];
      REAL(ret)[85] = REAL(theta)[36];
      REAL(ret)[95] = REAL(theta)[45];
    }
    else if (theta_n == 17){
      REAL(ret)[15] = exp(REAL(theta)[2]);
      REAL(ret)[25] = REAL(theta)[4];
      REAL(ret)[35] = REAL(theta)[7];
      REAL(ret)[45] = REAL(theta)[11];
      REAL(ret)[51] = exp(REAL(theta)[2]);
      REAL(ret)[52] = REAL(theta)[4];
      REAL(ret)[53] = REAL(theta)[7];
      REAL(ret)[54] = REAL(theta)[11];
      REAL(ret)[55] = 2 * REAL(theta)[16];
      REAL(ret)[56] = REAL(theta)[22];
      REAL(ret)[57] = REAL(theta)[29];
      REAL(ret)[58] = REAL(theta)[37];
      REAL(ret)[59] = REAL(theta)[46];
      REAL(ret)[65] = REAL(theta)[22];
      REAL(ret)[75] = REAL(theta)[29];
      REAL(ret)[85] = REAL(theta)[37];
      REAL(ret)[95] = REAL(theta)[46];
    }
    else if (theta_n == 18){
      REAL(ret)[25] = exp(REAL(theta)[5]);
      REAL(ret)[35] = REAL(theta)[8];
      REAL(ret)[45] = REAL(theta)[12];
      REAL(ret)[52] = exp(REAL(theta)[5]);
      REAL(ret)[53] = REAL(theta)[8];
      REAL(ret)[54] = REAL(theta)[12];
      REAL(ret)[55] = 2 * REAL(theta)[17];
      REAL(ret)[56] = REAL(theta)[23];
      REAL(ret)[57] = REAL(theta)[30];
      REAL(ret)[58] = REAL(theta)[38];
      REAL(ret)[59] = REAL(theta)[47];
      REAL(ret)[65] = REAL(theta)[23];
      REAL(ret)[75] = REAL(theta)[30];
      REAL(ret)[85] = REAL(theta)[38];
      REAL(ret)[95] = REAL(theta)[47];
    }
    else if (theta_n == 19){
      REAL(ret)[35] = exp(REAL(theta)[9]);
      REAL(ret)[45] = REAL(theta)[13];
      REAL(ret)[53] = exp(REAL(theta)[9]);
      REAL(ret)[54] = REAL(theta)[13];
      REAL(ret)[55] = 2 * REAL(theta)[18];
      REAL(ret)[56] = REAL(theta)[24];
      REAL(ret)[57] = REAL(theta)[31];
      REAL(ret)[58] = REAL(theta)[39];
      REAL(ret)[59] = REAL(theta)[48];
      REAL(ret)[65] = REAL(theta)[24];
      REAL(ret)[75] = REAL(theta)[31];
      REAL(ret)[85] = REAL(theta)[39];
      REAL(ret)[95] = REAL(theta)[48];
    }
    else if (theta_n == 20){
      REAL(ret)[45] = exp(REAL(theta)[14]);
      REAL(ret)[54] = exp(REAL(theta)[14]);
      REAL(ret)[55] = 2 * REAL(theta)[19];
      REAL(ret)[56] = REAL(theta)[25];
      REAL(ret)[57] = REAL(theta)[32];
      REAL(ret)[58] = REAL(theta)[40];
      REAL(ret)[59] = REAL(theta)[49];
      REAL(ret)[65] = REAL(theta)[25];
      REAL(ret)[75] = REAL(theta)[32];
      REAL(ret)[85] = REAL(theta)[40];
      REAL(ret)[95] = REAL(theta)[49];
    }
    else if (theta_n == 21){
      REAL(ret)[55] = 2 * exp(2 * REAL(theta)[20]);
      REAL(ret)[56] = REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[57] = REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[58] = REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[59] = REAL(theta)[50] * exp(REAL(theta)[20]);
      REAL(ret)[65] = REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[75] = REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[85] = REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[95] = REAL(theta)[50] * exp(REAL(theta)[20]);
    }
    else if (theta_n == 22){
      REAL(ret)[6] = exp(REAL(theta)[0]);
      REAL(ret)[16] = REAL(theta)[1];
      REAL(ret)[26] = REAL(theta)[3];
      REAL(ret)[36] = REAL(theta)[6];
      REAL(ret)[46] = REAL(theta)[10];
      REAL(ret)[56] = REAL(theta)[15];
      REAL(ret)[60] = exp(REAL(theta)[0]);
      REAL(ret)[61] = REAL(theta)[1];
      REAL(ret)[62] = REAL(theta)[3];
      REAL(ret)[63] = REAL(theta)[6];
      REAL(ret)[64] = REAL(theta)[10];
      REAL(ret)[65] = REAL(theta)[15];
      REAL(ret)[66] = 2 * REAL(theta)[21];
      REAL(ret)[67] = REAL(theta)[28];
      REAL(ret)[68] = REAL(theta)[36];
      REAL(ret)[69] = REAL(theta)[45];
      REAL(ret)[76] = REAL(theta)[28];
      REAL(ret)[86] = REAL(theta)[36];
      REAL(ret)[96] = REAL(theta)[45];
    }
    else if (theta_n == 23){
      REAL(ret)[16] = exp(REAL(theta)[2]);
      REAL(ret)[26] = REAL(theta)[4];
      REAL(ret)[36] = REAL(theta)[7];
      REAL(ret)[46] = REAL(theta)[11];
      REAL(ret)[56] = REAL(theta)[16];
      REAL(ret)[61] = exp(REAL(theta)[2]);
      REAL(ret)[62] = REAL(theta)[4];
      REAL(ret)[63] = REAL(theta)[7];
      REAL(ret)[64] = REAL(theta)[11];
      REAL(ret)[65] = REAL(theta)[16];
      REAL(ret)[66] = 2 * REAL(theta)[22];
      REAL(ret)[67] = REAL(theta)[29];
      REAL(ret)[68] = REAL(theta)[37];
      REAL(ret)[69] = REAL(theta)[46];
      REAL(ret)[76] = REAL(theta)[29];
      REAL(ret)[86] = REAL(theta)[37];
      REAL(ret)[96] = REAL(theta)[46];
    }
    else if (theta_n == 24){
      REAL(ret)[26] = exp(REAL(theta)[5]);
      REAL(ret)[36] = REAL(theta)[8];
      REAL(ret)[46] = REAL(theta)[12];
      REAL(ret)[56] = REAL(theta)[17];
      REAL(ret)[62] = exp(REAL(theta)[5]);
      REAL(ret)[63] = REAL(theta)[8];
      REAL(ret)[64] = REAL(theta)[12];
      REAL(ret)[65] = REAL(theta)[17];
      REAL(ret)[66] = 2 * REAL(theta)[23];
      REAL(ret)[67] = REAL(theta)[30];
      REAL(ret)[68] = REAL(theta)[38];
      REAL(ret)[69] = REAL(theta)[47];
      REAL(ret)[76] = REAL(theta)[30];
      REAL(ret)[86] = REAL(theta)[38];
      REAL(ret)[96] = REAL(theta)[47];
    }
    else if (theta_n == 25){
      REAL(ret)[36] = exp(REAL(theta)[9]);
      REAL(ret)[46] = REAL(theta)[13];
      REAL(ret)[56] = REAL(theta)[18];
      REAL(ret)[63] = exp(REAL(theta)[9]);
      REAL(ret)[64] = REAL(theta)[13];
      REAL(ret)[65] = REAL(theta)[18];
      REAL(ret)[66] = 2 * REAL(theta)[24];
      REAL(ret)[67] = REAL(theta)[31];
      REAL(ret)[68] = REAL(theta)[39];
      REAL(ret)[69] = REAL(theta)[48];
      REAL(ret)[76] = REAL(theta)[31];
      REAL(ret)[86] = REAL(theta)[39];
      REAL(ret)[96] = REAL(theta)[48];
    }
    else if (theta_n == 26){
      REAL(ret)[46] = exp(REAL(theta)[14]);
      REAL(ret)[56] = REAL(theta)[19];
      REAL(ret)[64] = exp(REAL(theta)[14]);
      REAL(ret)[65] = REAL(theta)[19];
      REAL(ret)[66] = 2 * REAL(theta)[25];
      REAL(ret)[67] = REAL(theta)[32];
      REAL(ret)[68] = REAL(theta)[40];
      REAL(ret)[69] = REAL(theta)[49];
      REAL(ret)[76] = REAL(theta)[32];
      REAL(ret)[86] = REAL(theta)[40];
      REAL(ret)[96] = REAL(theta)[49];
    }
    else if (theta_n == 27){
      REAL(ret)[56] = exp(REAL(theta)[20]);
      REAL(ret)[65] = exp(REAL(theta)[20]);
      REAL(ret)[66] = 2 * REAL(theta)[26];
      REAL(ret)[67] = REAL(theta)[33];
      REAL(ret)[68] = REAL(theta)[41];
      REAL(ret)[69] = REAL(theta)[50];
      REAL(ret)[76] = REAL(theta)[33];
      REAL(ret)[86] = REAL(theta)[41];
      REAL(ret)[96] = REAL(theta)[50];
    }
    else if (theta_n == 28){
      REAL(ret)[66] = 2 * exp(2 * REAL(theta)[27]);
      REAL(ret)[67] = REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[68] = REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[69] = REAL(theta)[51] * exp(REAL(theta)[27]);
      REAL(ret)[76] = REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[86] = REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[96] = REAL(theta)[51] * exp(REAL(theta)[27]);
    }
    else if (theta_n == 29){
      REAL(ret)[7] = exp(REAL(theta)[0]);
      REAL(ret)[17] = REAL(theta)[1];
      REAL(ret)[27] = REAL(theta)[3];
      REAL(ret)[37] = REAL(theta)[6];
      REAL(ret)[47] = REAL(theta)[10];
      REAL(ret)[57] = REAL(theta)[15];
      REAL(ret)[67] = REAL(theta)[21];
      REAL(ret)[70] = exp(REAL(theta)[0]);
      REAL(ret)[71] = REAL(theta)[1];
      REAL(ret)[72] = REAL(theta)[3];
      REAL(ret)[73] = REAL(theta)[6];
      REAL(ret)[74] = REAL(theta)[10];
      REAL(ret)[75] = REAL(theta)[15];
      REAL(ret)[76] = REAL(theta)[21];
      REAL(ret)[77] = 2 * REAL(theta)[28];
      REAL(ret)[78] = REAL(theta)[36];
      REAL(ret)[79] = REAL(theta)[45];
      REAL(ret)[87] = REAL(theta)[36];
      REAL(ret)[97] = REAL(theta)[45];
    }
    else if (theta_n == 30){
      REAL(ret)[17] = exp(REAL(theta)[2]);
      REAL(ret)[27] = REAL(theta)[4];
      REAL(ret)[37] = REAL(theta)[7];
      REAL(ret)[47] = REAL(theta)[11];
      REAL(ret)[57] = REAL(theta)[16];
      REAL(ret)[67] = REAL(theta)[22];
      REAL(ret)[71] = exp(REAL(theta)[2]);
      REAL(ret)[72] = REAL(theta)[4];
      REAL(ret)[73] = REAL(theta)[7];
      REAL(ret)[74] = REAL(theta)[11];
      REAL(ret)[75] = REAL(theta)[16];
      REAL(ret)[76] = REAL(theta)[22];
      REAL(ret)[77] = 2 * REAL(theta)[29];
      REAL(ret)[78] = REAL(theta)[37];
      REAL(ret)[79] = REAL(theta)[46];
      REAL(ret)[87] = REAL(theta)[37];
      REAL(ret)[97] = REAL(theta)[46];
    }
    else if (theta_n == 31){
      REAL(ret)[27] = exp(REAL(theta)[5]);
      REAL(ret)[37] = REAL(theta)[8];
      REAL(ret)[47] = REAL(theta)[12];
      REAL(ret)[57] = REAL(theta)[17];
      REAL(ret)[67] = REAL(theta)[23];
      REAL(ret)[72] = exp(REAL(theta)[5]);
      REAL(ret)[73] = REAL(theta)[8];
      REAL(ret)[74] = REAL(theta)[12];
      REAL(ret)[75] = REAL(theta)[17];
      REAL(ret)[76] = REAL(theta)[23];
      REAL(ret)[77] = 2 * REAL(theta)[30];
      REAL(ret)[78] = REAL(theta)[38];
      REAL(ret)[79] = REAL(theta)[47];
      REAL(ret)[87] = REAL(theta)[38];
      REAL(ret)[97] = REAL(theta)[47];
    }
    else if (theta_n == 32){
      REAL(ret)[37] = exp(REAL(theta)[9]);
      REAL(ret)[47] = REAL(theta)[13];
      REAL(ret)[57] = REAL(theta)[18];
      REAL(ret)[67] = REAL(theta)[24];
      REAL(ret)[73] = exp(REAL(theta)[9]);
      REAL(ret)[74] = REAL(theta)[13];
      REAL(ret)[75] = REAL(theta)[18];
      REAL(ret)[76] = REAL(theta)[24];
      REAL(ret)[77] = 2 * REAL(theta)[31];
      REAL(ret)[78] = REAL(theta)[39];
      REAL(ret)[79] = REAL(theta)[48];
      REAL(ret)[87] = REAL(theta)[39];
      REAL(ret)[97] = REAL(theta)[48];
    }
    else if (theta_n == 33){
      REAL(ret)[47] = exp(REAL(theta)[14]);
      REAL(ret)[57] = REAL(theta)[19];
      REAL(ret)[67] = REAL(theta)[25];
      REAL(ret)[74] = exp(REAL(theta)[14]);
      REAL(ret)[75] = REAL(theta)[19];
      REAL(ret)[76] = REAL(theta)[25];
      REAL(ret)[77] = 2 * REAL(theta)[32];
      REAL(ret)[78] = REAL(theta)[40];
      REAL(ret)[79] = REAL(theta)[49];
      REAL(ret)[87] = REAL(theta)[40];
      REAL(ret)[97] = REAL(theta)[49];
    }
    else if (theta_n == 34){
      REAL(ret)[57] = exp(REAL(theta)[20]);
      REAL(ret)[67] = REAL(theta)[26];
      REAL(ret)[75] = exp(REAL(theta)[20]);
      REAL(ret)[76] = REAL(theta)[26];
      REAL(ret)[77] = 2 * REAL(theta)[33];
      REAL(ret)[78] = REAL(theta)[41];
      REAL(ret)[79] = REAL(theta)[50];
      REAL(ret)[87] = REAL(theta)[41];
      REAL(ret)[97] = REAL(theta)[50];
    }
    else if (theta_n == 35){
      REAL(ret)[67] = exp(REAL(theta)[27]);
      REAL(ret)[76] = exp(REAL(theta)[27]);
      REAL(ret)[77] = 2 * REAL(theta)[34];
      REAL(ret)[78] = REAL(theta)[42];
      REAL(ret)[79] = REAL(theta)[51];
      REAL(ret)[87] = REAL(theta)[42];
      REAL(ret)[97] = REAL(theta)[51];
    }
    else if (theta_n == 36){
      REAL(ret)[77] = 2 * exp(2 * REAL(theta)[35]);
      REAL(ret)[78] = REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[79] = REAL(theta)[52] * exp(REAL(theta)[35]);
      REAL(ret)[87] = REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[97] = REAL(theta)[52] * exp(REAL(theta)[35]);
    }
    else if (theta_n == 37){
      REAL(ret)[8] = exp(REAL(theta)[0]);
      REAL(ret)[18] = REAL(theta)[1];
      REAL(ret)[28] = REAL(theta)[3];
      REAL(ret)[38] = REAL(theta)[6];
      REAL(ret)[48] = REAL(theta)[10];
      REAL(ret)[58] = REAL(theta)[15];
      REAL(ret)[68] = REAL(theta)[21];
      REAL(ret)[78] = REAL(theta)[28];
      REAL(ret)[80] = exp(REAL(theta)[0]);
      REAL(ret)[81] = REAL(theta)[1];
      REAL(ret)[82] = REAL(theta)[3];
      REAL(ret)[83] = REAL(theta)[6];
      REAL(ret)[84] = REAL(theta)[10];
      REAL(ret)[85] = REAL(theta)[15];
      REAL(ret)[86] = REAL(theta)[21];
      REAL(ret)[87] = REAL(theta)[28];
      REAL(ret)[88] = 2 * REAL(theta)[36];
      REAL(ret)[89] = REAL(theta)[45];
      REAL(ret)[98] = REAL(theta)[45];
    }
    else if (theta_n == 38){
      REAL(ret)[18] = exp(REAL(theta)[2]);
      REAL(ret)[28] = REAL(theta)[4];
      REAL(ret)[38] = REAL(theta)[7];
      REAL(ret)[48] = REAL(theta)[11];
      REAL(ret)[58] = REAL(theta)[16];
      REAL(ret)[68] = REAL(theta)[22];
      REAL(ret)[78] = REAL(theta)[29];
      REAL(ret)[81] = exp(REAL(theta)[2]);
      REAL(ret)[82] = REAL(theta)[4];
      REAL(ret)[83] = REAL(theta)[7];
      REAL(ret)[84] = REAL(theta)[11];
      REAL(ret)[85] = REAL(theta)[16];
      REAL(ret)[86] = REAL(theta)[22];
      REAL(ret)[87] = REAL(theta)[29];
      REAL(ret)[88] = 2 * REAL(theta)[37];
      REAL(ret)[89] = REAL(theta)[46];
      REAL(ret)[98] = REAL(theta)[46];
    }
    else if (theta_n == 39){
      REAL(ret)[28] = exp(REAL(theta)[5]);
      REAL(ret)[38] = REAL(theta)[8];
      REAL(ret)[48] = REAL(theta)[12];
      REAL(ret)[58] = REAL(theta)[17];
      REAL(ret)[68] = REAL(theta)[23];
      REAL(ret)[78] = REAL(theta)[30];
      REAL(ret)[82] = exp(REAL(theta)[5]);
      REAL(ret)[83] = REAL(theta)[8];
      REAL(ret)[84] = REAL(theta)[12];
      REAL(ret)[85] = REAL(theta)[17];
      REAL(ret)[86] = REAL(theta)[23];
      REAL(ret)[87] = REAL(theta)[30];
      REAL(ret)[88] = 2 * REAL(theta)[38];
      REAL(ret)[89] = REAL(theta)[47];
      REAL(ret)[98] = REAL(theta)[47];
    }
    else if (theta_n == 40){
      REAL(ret)[38] = exp(REAL(theta)[9]);
      REAL(ret)[48] = REAL(theta)[13];
      REAL(ret)[58] = REAL(theta)[18];
      REAL(ret)[68] = REAL(theta)[24];
      REAL(ret)[78] = REAL(theta)[31];
      REAL(ret)[83] = exp(REAL(theta)[9]);
      REAL(ret)[84] = REAL(theta)[13];
      REAL(ret)[85] = REAL(theta)[18];
      REAL(ret)[86] = REAL(theta)[24];
      REAL(ret)[87] = REAL(theta)[31];
      REAL(ret)[88] = 2 * REAL(theta)[39];
      REAL(ret)[89] = REAL(theta)[48];
      REAL(ret)[98] = REAL(theta)[48];
    }
    else if (theta_n == 41){
      REAL(ret)[48] = exp(REAL(theta)[14]);
      REAL(ret)[58] = REAL(theta)[19];
      REAL(ret)[68] = REAL(theta)[25];
      REAL(ret)[78] = REAL(theta)[32];
      REAL(ret)[84] = exp(REAL(theta)[14]);
      REAL(ret)[85] = REAL(theta)[19];
      REAL(ret)[86] = REAL(theta)[25];
      REAL(ret)[87] = REAL(theta)[32];
      REAL(ret)[88] = 2 * REAL(theta)[40];
      REAL(ret)[89] = REAL(theta)[49];
      REAL(ret)[98] = REAL(theta)[49];
    }
    else if (theta_n == 42){
      REAL(ret)[58] = exp(REAL(theta)[20]);
      REAL(ret)[68] = REAL(theta)[26];
      REAL(ret)[78] = REAL(theta)[33];
      REAL(ret)[85] = exp(REAL(theta)[20]);
      REAL(ret)[86] = REAL(theta)[26];
      REAL(ret)[87] = REAL(theta)[33];
      REAL(ret)[88] = 2 * REAL(theta)[41];
      REAL(ret)[89] = REAL(theta)[50];
      REAL(ret)[98] = REAL(theta)[50];
    }
    else if (theta_n == 43){
      REAL(ret)[68] = exp(REAL(theta)[27]);
      REAL(ret)[78] = REAL(theta)[34];
      REAL(ret)[86] = exp(REAL(theta)[27]);
      REAL(ret)[87] = REAL(theta)[34];
      REAL(ret)[88] = 2 * REAL(theta)[42];
      REAL(ret)[89] = REAL(theta)[51];
      REAL(ret)[98] = REAL(theta)[51];
    }
    else if (theta_n == 44){
      REAL(ret)[78] = exp(REAL(theta)[35]);
      REAL(ret)[87] = exp(REAL(theta)[35]);
      REAL(ret)[88] = 2 * REAL(theta)[43];
      REAL(ret)[89] = REAL(theta)[52];
      REAL(ret)[98] = REAL(theta)[52];
    }
    else if (theta_n == 45){
      REAL(ret)[88] = 2 * exp(2 * REAL(theta)[44]);
      REAL(ret)[89] = REAL(theta)[53] * exp(REAL(theta)[44]);
      REAL(ret)[98] = REAL(theta)[53] * exp(REAL(theta)[44]);
    }
    else if (theta_n == 46){
      REAL(ret)[9] = exp(REAL(theta)[0]);
      REAL(ret)[19] = REAL(theta)[1];
      REAL(ret)[29] = REAL(theta)[3];
      REAL(ret)[39] = REAL(theta)[6];
      REAL(ret)[49] = REAL(theta)[10];
      REAL(ret)[59] = REAL(theta)[15];
      REAL(ret)[69] = REAL(theta)[21];
      REAL(ret)[79] = REAL(theta)[28];
      REAL(ret)[89] = REAL(theta)[36];
      REAL(ret)[90] = exp(REAL(theta)[0]);
      REAL(ret)[91] = REAL(theta)[1];
      REAL(ret)[92] = REAL(theta)[3];
      REAL(ret)[93] = REAL(theta)[6];
      REAL(ret)[94] = REAL(theta)[10];
      REAL(ret)[95] = REAL(theta)[15];
      REAL(ret)[96] = REAL(theta)[21];
      REAL(ret)[97] = REAL(theta)[28];
      REAL(ret)[98] = REAL(theta)[36];
      REAL(ret)[99] = 2 * REAL(theta)[45];
    }
    else if (theta_n == 47){
      REAL(ret)[19] = exp(REAL(theta)[2]);
      REAL(ret)[29] = REAL(theta)[4];
      REAL(ret)[39] = REAL(theta)[7];
      REAL(ret)[49] = REAL(theta)[11];
      REAL(ret)[59] = REAL(theta)[16];
      REAL(ret)[69] = REAL(theta)[22];
      REAL(ret)[79] = REAL(theta)[29];
      REAL(ret)[89] = REAL(theta)[37];
      REAL(ret)[91] = exp(REAL(theta)[2]);
      REAL(ret)[92] = REAL(theta)[4];
      REAL(ret)[93] = REAL(theta)[7];
      REAL(ret)[94] = REAL(theta)[11];
      REAL(ret)[95] = REAL(theta)[16];
      REAL(ret)[96] = REAL(theta)[22];
      REAL(ret)[97] = REAL(theta)[29];
      REAL(ret)[98] = REAL(theta)[37];
      REAL(ret)[99] = 2 * REAL(theta)[46];
    }
    else if (theta_n == 48){
      REAL(ret)[29] = exp(REAL(theta)[5]);
      REAL(ret)[39] = REAL(theta)[8];
      REAL(ret)[49] = REAL(theta)[12];
      REAL(ret)[59] = REAL(theta)[17];
      REAL(ret)[69] = REAL(theta)[23];
      REAL(ret)[79] = REAL(theta)[30];
      REAL(ret)[89] = REAL(theta)[38];
      REAL(ret)[92] = exp(REAL(theta)[5]);
      REAL(ret)[93] = REAL(theta)[8];
      REAL(ret)[94] = REAL(theta)[12];
      REAL(ret)[95] = REAL(theta)[17];
      REAL(ret)[96] = REAL(theta)[23];
      REAL(ret)[97] = REAL(theta)[30];
      REAL(ret)[98] = REAL(theta)[38];
      REAL(ret)[99] = 2 * REAL(theta)[47];
    }
    else if (theta_n == 49){
      REAL(ret)[39] = exp(REAL(theta)[9]);
      REAL(ret)[49] = REAL(theta)[13];
      REAL(ret)[59] = REAL(theta)[18];
      REAL(ret)[69] = REAL(theta)[24];
      REAL(ret)[79] = REAL(theta)[31];
      REAL(ret)[89] = REAL(theta)[39];
      REAL(ret)[93] = exp(REAL(theta)[9]);
      REAL(ret)[94] = REAL(theta)[13];
      REAL(ret)[95] = REAL(theta)[18];
      REAL(ret)[96] = REAL(theta)[24];
      REAL(ret)[97] = REAL(theta)[31];
      REAL(ret)[98] = REAL(theta)[39];
      REAL(ret)[99] = 2 * REAL(theta)[48];
    }
    else if (theta_n == 50){
      REAL(ret)[49] = exp(REAL(theta)[14]);
      REAL(ret)[59] = REAL(theta)[19];
      REAL(ret)[69] = REAL(theta)[25];
      REAL(ret)[79] = REAL(theta)[32];
      REAL(ret)[89] = REAL(theta)[40];
      REAL(ret)[94] = exp(REAL(theta)[14]);
      REAL(ret)[95] = REAL(theta)[19];
      REAL(ret)[96] = REAL(theta)[25];
      REAL(ret)[97] = REAL(theta)[32];
      REAL(ret)[98] = REAL(theta)[40];
      REAL(ret)[99] = 2 * REAL(theta)[49];
    }
    else if (theta_n == 51){
      REAL(ret)[59] = exp(REAL(theta)[20]);
      REAL(ret)[69] = REAL(theta)[26];
      REAL(ret)[79] = REAL(theta)[33];
      REAL(ret)[89] = REAL(theta)[41];
      REAL(ret)[95] = exp(REAL(theta)[20]);
      REAL(ret)[96] = REAL(theta)[26];
      REAL(ret)[97] = REAL(theta)[33];
      REAL(ret)[98] = REAL(theta)[41];
      REAL(ret)[99] = 2 * REAL(theta)[50];
    }
    else if (theta_n == 52){
      REAL(ret)[69] = exp(REAL(theta)[27]);
      REAL(ret)[79] = REAL(theta)[34];
      REAL(ret)[89] = REAL(theta)[42];
      REAL(ret)[96] = exp(REAL(theta)[27]);
      REAL(ret)[97] = REAL(theta)[34];
      REAL(ret)[98] = REAL(theta)[42];
      REAL(ret)[99] = 2 * REAL(theta)[51];
    }
    else if (theta_n == 53){
      REAL(ret)[79] = exp(REAL(theta)[35]);
      REAL(ret)[89] = REAL(theta)[43];
      REAL(ret)[97] = exp(REAL(theta)[35]);
      REAL(ret)[98] = REAL(theta)[43];
      REAL(ret)[99] = 2 * REAL(theta)[52];
    }
    else if (theta_n == 54){
      REAL(ret)[89] = exp(REAL(theta)[44]);
      REAL(ret)[98] = exp(REAL(theta)[44]);
      REAL(ret)[99] = 2 * REAL(theta)[53];
    }
    else if (theta_n == 55){
      REAL(ret)[99] = 2 * exp(2 * REAL(theta)[54]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 10));for(int i = 0; i < 10; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 2 * exp(2 * REAL(theta)[2]);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 2 * exp(2 * REAL(theta)[5]);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 2 * exp(2 * REAL(theta)[9]);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 2 * exp(2 * REAL(theta)[14]);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 2 * exp(2 * REAL(theta)[20]);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 2 * exp(2 * REAL(theta)[27]);
    }
    else if (theta_n == -38){
      REAL(ret)[7] = 2 * exp(2 * REAL(theta)[35]);
    }
    else if (theta_n == -47){
      REAL(ret)[8] = 2 * exp(2 * REAL(theta)[44]);
    }
    else if (theta_n == -57){
      REAL(ret)[9] = 2 * exp(2 * REAL(theta)[54]);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 11){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 11;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -68 || theta_n > 66){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 66){
    error("Requires vector with 66 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 11, 11));for (int i = 0; i < 121; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[12] = exp(REAL(theta)[2]);
      REAL(ret)[22] = REAL(theta)[3];
      REAL(ret)[23] = REAL(theta)[4];
      REAL(ret)[24] = exp(REAL(theta)[5]);
      REAL(ret)[33] = REAL(theta)[6];
      REAL(ret)[34] = REAL(theta)[7];
      REAL(ret)[35] = REAL(theta)[8];
      REAL(ret)[36] = exp(REAL(theta)[9]);
      REAL(ret)[44] = REAL(theta)[10];
      REAL(ret)[45] = REAL(theta)[11];
      REAL(ret)[46] = REAL(theta)[12];
      REAL(ret)[47] = REAL(theta)[13];
      REAL(ret)[48] = exp(REAL(theta)[14]);
      REAL(ret)[55] = REAL(theta)[15];
      REAL(ret)[56] = REAL(theta)[16];
      REAL(ret)[57] = REAL(theta)[17];
      REAL(ret)[58] = REAL(theta)[18];
      REAL(ret)[59] = REAL(theta)[19];
      REAL(ret)[60] = exp(REAL(theta)[20]);
      REAL(ret)[66] = REAL(theta)[21];
      REAL(ret)[67] = REAL(theta)[22];
      REAL(ret)[68] = REAL(theta)[23];
      REAL(ret)[69] = REAL(theta)[24];
      REAL(ret)[70] = REAL(theta)[25];
      REAL(ret)[71] = REAL(theta)[26];
      REAL(ret)[72] = exp(REAL(theta)[27]);
      REAL(ret)[77] = REAL(theta)[28];
      REAL(ret)[78] = REAL(theta)[29];
      REAL(ret)[79] = REAL(theta)[30];
      REAL(ret)[80] = REAL(theta)[31];
      REAL(ret)[81] = REAL(theta)[32];
      REAL(ret)[82] = REAL(theta)[33];
      REAL(ret)[83] = REAL(theta)[34];
      REAL(ret)[84] = exp(REAL(theta)[35]);
      REAL(ret)[88] = REAL(theta)[36];
      REAL(ret)[89] = REAL(theta)[37];
      REAL(ret)[90] = REAL(theta)[38];
      REAL(ret)[91] = REAL(theta)[39];
      REAL(ret)[92] = REAL(theta)[40];
      REAL(ret)[93] = REAL(theta)[41];
      REAL(ret)[94] = REAL(theta)[42];
      REAL(ret)[95] = REAL(theta)[43];
      REAL(ret)[96] = exp(REAL(theta)[44]);
      REAL(ret)[99] = REAL(theta)[45];
      REAL(ret)[100] = REAL(theta)[46];
      REAL(ret)[101] = REAL(theta)[47];
      REAL(ret)[102] = REAL(theta)[48];
      REAL(ret)[103] = REAL(theta)[49];
      REAL(ret)[104] = REAL(theta)[50];
      REAL(ret)[105] = REAL(theta)[51];
      REAL(ret)[106] = REAL(theta)[52];
      REAL(ret)[107] = REAL(theta)[53];
      REAL(ret)[108] = exp(REAL(theta)[54]);
      REAL(ret)[110] = REAL(theta)[55];
      REAL(ret)[111] = REAL(theta)[56];
      REAL(ret)[112] = REAL(theta)[57];
      REAL(ret)[113] = REAL(theta)[58];
      REAL(ret)[114] = REAL(theta)[59];
      REAL(ret)[115] = REAL(theta)[60];
      REAL(ret)[116] = REAL(theta)[61];
      REAL(ret)[117] = REAL(theta)[62];
      REAL(ret)[118] = REAL(theta)[63];
      REAL(ret)[119] = REAL(theta)[64];
      REAL(ret)[120] = exp(REAL(theta)[65]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[45] * exp(REAL(theta)[0]);
      REAL(ret)[10] = REAL(theta)[55] * exp(REAL(theta)[0]);
      REAL(ret)[11] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[12] = R_pow_di(REAL(theta)[1], 2) + exp(2 * REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[18] = REAL(theta)[1] * REAL(theta)[28] + REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[19] = REAL(theta)[1] * REAL(theta)[36] + REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[20] = REAL(theta)[1] * REAL(theta)[45] + REAL(theta)[46] * exp(REAL(theta)[2]);
      REAL(ret)[21] = REAL(theta)[1] * REAL(theta)[55] + REAL(theta)[56] * exp(REAL(theta)[2]);
      REAL(ret)[22] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[23] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[24] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + exp(2 * REAL(theta)[5]);
      REAL(ret)[25] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[26] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[27] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[28] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[29] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[30] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[31] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * exp(REAL(theta)[5]);
      REAL(ret)[32] = REAL(theta)[3] * REAL(theta)[55] + REAL(theta)[4] * REAL(theta)[56] + REAL(theta)[57] * exp(REAL(theta)[5]);
      REAL(ret)[33] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[34] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[35] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[36] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + exp(2 * REAL(theta)[9]);
      REAL(ret)[37] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[38] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[39] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[40] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[41] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[42] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * exp(REAL(theta)[9]);
      REAL(ret)[43] = REAL(theta)[55] * REAL(theta)[6] + REAL(theta)[56] * REAL(theta)[7] + REAL(theta)[57] * REAL(theta)[8] + REAL(theta)[58] * exp(REAL(theta)[9]);
      REAL(ret)[44] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[45] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[46] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[47] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[48] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + exp(2 * REAL(theta)[14]);
      REAL(ret)[49] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[50] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[51] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[52] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[53] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + REAL(theta)[49] * exp(REAL(theta)[14]);
      REAL(ret)[54] = REAL(theta)[10] * REAL(theta)[55] + REAL(theta)[11] * REAL(theta)[56] + REAL(theta)[12] * REAL(theta)[57] + REAL(theta)[13] * REAL(theta)[58] + REAL(theta)[59] * exp(REAL(theta)[14]);
      REAL(ret)[55] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[56] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[57] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[58] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[59] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[60] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + exp(2 * REAL(theta)[20]);
      REAL(ret)[61] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[62] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[63] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[64] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + REAL(theta)[50] * exp(REAL(theta)[20]);
      REAL(ret)[65] = REAL(theta)[15] * REAL(theta)[55] + REAL(theta)[16] * REAL(theta)[56] + REAL(theta)[17] * REAL(theta)[57] + REAL(theta)[18] * REAL(theta)[58] + REAL(theta)[19] * REAL(theta)[59] + REAL(theta)[60] * exp(REAL(theta)[20]);
      REAL(ret)[66] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[67] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[68] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[69] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[70] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[71] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[72] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + exp(2 * REAL(theta)[27]);
      REAL(ret)[73] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[74] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[75] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + REAL(theta)[51] * exp(REAL(theta)[27]);
      REAL(ret)[76] = REAL(theta)[21] * REAL(theta)[55] + REAL(theta)[22] * REAL(theta)[56] + REAL(theta)[23] * REAL(theta)[57] + REAL(theta)[24] * REAL(theta)[58] + REAL(theta)[25] * REAL(theta)[59] + REAL(theta)[26] * REAL(theta)[60] + REAL(theta)[61] * exp(REAL(theta)[27]);
      REAL(ret)[77] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[78] = REAL(theta)[1] * REAL(theta)[28] + REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[79] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[80] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[81] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[82] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[83] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[84] = R_pow_di(REAL(theta)[28], 2) + R_pow_di(REAL(theta)[29], 2) + R_pow_di(REAL(theta)[30], 2) + R_pow_di(REAL(theta)[31], 2) + R_pow_di(REAL(theta)[32], 2) + R_pow_di(REAL(theta)[33], 2) + R_pow_di(REAL(theta)[34], 2) + exp(2 * REAL(theta)[35]);
      REAL(ret)[85] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[86] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + REAL(theta)[52] * exp(REAL(theta)[35]);
      REAL(ret)[87] = REAL(theta)[28] * REAL(theta)[55] + REAL(theta)[29] * REAL(theta)[56] + REAL(theta)[30] * REAL(theta)[57] + REAL(theta)[31] * REAL(theta)[58] + REAL(theta)[32] * REAL(theta)[59] + REAL(theta)[33] * REAL(theta)[60] + REAL(theta)[34] * REAL(theta)[61] + REAL(theta)[62] * exp(REAL(theta)[35]);
      REAL(ret)[88] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[89] = REAL(theta)[1] * REAL(theta)[36] + REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[90] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[91] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[92] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[93] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[94] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[95] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[96] = R_pow_di(REAL(theta)[36], 2) + R_pow_di(REAL(theta)[37], 2) + R_pow_di(REAL(theta)[38], 2) + R_pow_di(REAL(theta)[39], 2) + R_pow_di(REAL(theta)[40], 2) + R_pow_di(REAL(theta)[41], 2) + R_pow_di(REAL(theta)[42], 2) + R_pow_di(REAL(theta)[43], 2) + exp(2 * REAL(theta)[44]);
      REAL(ret)[97] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + REAL(theta)[53] * exp(REAL(theta)[44]);
      REAL(ret)[98] = REAL(theta)[36] * REAL(theta)[55] + REAL(theta)[37] * REAL(theta)[56] + REAL(theta)[38] * REAL(theta)[57] + REAL(theta)[39] * REAL(theta)[58] + REAL(theta)[40] * REAL(theta)[59] + REAL(theta)[41] * REAL(theta)[60] + REAL(theta)[42] * REAL(theta)[61] + REAL(theta)[43] * REAL(theta)[62] + REAL(theta)[63] * exp(REAL(theta)[44]);
      REAL(ret)[99] = REAL(theta)[45] * exp(REAL(theta)[0]);
      REAL(ret)[100] = REAL(theta)[1] * REAL(theta)[45] + REAL(theta)[46] * exp(REAL(theta)[2]);
      REAL(ret)[101] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * exp(REAL(theta)[5]);
      REAL(ret)[102] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * exp(REAL(theta)[9]);
      REAL(ret)[103] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + REAL(theta)[49] * exp(REAL(theta)[14]);
      REAL(ret)[104] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + REAL(theta)[50] * exp(REAL(theta)[20]);
      REAL(ret)[105] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + REAL(theta)[51] * exp(REAL(theta)[27]);
      REAL(ret)[106] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + REAL(theta)[52] * exp(REAL(theta)[35]);
      REAL(ret)[107] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + REAL(theta)[53] * exp(REAL(theta)[44]);
      REAL(ret)[108] = R_pow_di(REAL(theta)[45], 2) + R_pow_di(REAL(theta)[46], 2) + R_pow_di(REAL(theta)[47], 2) + R_pow_di(REAL(theta)[48], 2) + R_pow_di(REAL(theta)[49], 2) + R_pow_di(REAL(theta)[50], 2) + R_pow_di(REAL(theta)[51], 2) + R_pow_di(REAL(theta)[52], 2) + R_pow_di(REAL(theta)[53], 2) + exp(2 * REAL(theta)[54]);
      REAL(ret)[109] = REAL(theta)[45] * REAL(theta)[55] + REAL(theta)[46] * REAL(theta)[56] + REAL(theta)[47] * REAL(theta)[57] + REAL(theta)[48] * REAL(theta)[58] + REAL(theta)[49] * REAL(theta)[59] + REAL(theta)[50] * REAL(theta)[60] + REAL(theta)[51] * REAL(theta)[61] + REAL(theta)[52] * REAL(theta)[62] + REAL(theta)[53] * REAL(theta)[63] + REAL(theta)[64] * exp(REAL(theta)[54]);
      REAL(ret)[110] = REAL(theta)[55] * exp(REAL(theta)[0]);
      REAL(ret)[111] = REAL(theta)[1] * REAL(theta)[55] + REAL(theta)[56] * exp(REAL(theta)[2]);
      REAL(ret)[112] = REAL(theta)[3] * REAL(theta)[55] + REAL(theta)[4] * REAL(theta)[56] + REAL(theta)[57] * exp(REAL(theta)[5]);
      REAL(ret)[113] = REAL(theta)[55] * REAL(theta)[6] + REAL(theta)[56] * REAL(theta)[7] + REAL(theta)[57] * REAL(theta)[8] + REAL(theta)[58] * exp(REAL(theta)[9]);
      REAL(ret)[114] = REAL(theta)[10] * REAL(theta)[55] + REAL(theta)[11] * REAL(theta)[56] + REAL(theta)[12] * REAL(theta)[57] + REAL(theta)[13] * REAL(theta)[58] + REAL(theta)[59] * exp(REAL(theta)[14]);
      REAL(ret)[115] = REAL(theta)[15] * REAL(theta)[55] + REAL(theta)[16] * REAL(theta)[56] + REAL(theta)[17] * REAL(theta)[57] + REAL(theta)[18] * REAL(theta)[58] + REAL(theta)[19] * REAL(theta)[59] + REAL(theta)[60] * exp(REAL(theta)[20]);
      REAL(ret)[116] = REAL(theta)[21] * REAL(theta)[55] + REAL(theta)[22] * REAL(theta)[56] + REAL(theta)[23] * REAL(theta)[57] + REAL(theta)[24] * REAL(theta)[58] + REAL(theta)[25] * REAL(theta)[59] + REAL(theta)[26] * REAL(theta)[60] + REAL(theta)[61] * exp(REAL(theta)[27]);
      REAL(ret)[117] = REAL(theta)[28] * REAL(theta)[55] + REAL(theta)[29] * REAL(theta)[56] + REAL(theta)[30] * REAL(theta)[57] + REAL(theta)[31] * REAL(theta)[58] + REAL(theta)[32] * REAL(theta)[59] + REAL(theta)[33] * REAL(theta)[60] + REAL(theta)[34] * REAL(theta)[61] + REAL(theta)[62] * exp(REAL(theta)[35]);
      REAL(ret)[118] = REAL(theta)[36] * REAL(theta)[55] + REAL(theta)[37] * REAL(theta)[56] + REAL(theta)[38] * REAL(theta)[57] + REAL(theta)[39] * REAL(theta)[58] + REAL(theta)[40] * REAL(theta)[59] + REAL(theta)[41] * REAL(theta)[60] + REAL(theta)[42] * REAL(theta)[61] + REAL(theta)[43] * REAL(theta)[62] + REAL(theta)[63] * exp(REAL(theta)[44]);
      REAL(ret)[119] = REAL(theta)[45] * REAL(theta)[55] + REAL(theta)[46] * REAL(theta)[56] + REAL(theta)[47] * REAL(theta)[57] + REAL(theta)[48] * REAL(theta)[58] + REAL(theta)[49] * REAL(theta)[59] + REAL(theta)[50] * REAL(theta)[60] + REAL(theta)[51] * REAL(theta)[61] + REAL(theta)[52] * REAL(theta)[62] + REAL(theta)[53] * REAL(theta)[63] + REAL(theta)[64] * exp(REAL(theta)[54]);
      REAL(ret)[120] = R_pow_di(REAL(theta)[55], 2) + R_pow_di(REAL(theta)[56], 2) + R_pow_di(REAL(theta)[57], 2) + R_pow_di(REAL(theta)[58], 2) + R_pow_di(REAL(theta)[59], 2) + R_pow_di(REAL(theta)[60], 2) + R_pow_di(REAL(theta)[61], 2) + R_pow_di(REAL(theta)[62], 2) + R_pow_di(REAL(theta)[63], 2) + R_pow_di(REAL(theta)[64], 2) + exp(2 * REAL(theta)[65]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[45] * exp(REAL(theta)[0]);
      REAL(ret)[10] = REAL(theta)[55] * exp(REAL(theta)[0]);
      REAL(ret)[11] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[22] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[33] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[44] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[55] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[66] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[77] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[88] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[99] = REAL(theta)[45] * exp(REAL(theta)[0]);
      REAL(ret)[110] = REAL(theta)[55] * exp(REAL(theta)[0]);
    }
    else if (theta_n == 2){
      REAL(ret)[1] = exp(REAL(theta)[0]);
      REAL(ret)[11] = exp(REAL(theta)[0]);
      REAL(ret)[12] = 2 * REAL(theta)[1];
      REAL(ret)[13] = REAL(theta)[3];
      REAL(ret)[14] = REAL(theta)[6];
      REAL(ret)[15] = REAL(theta)[10];
      REAL(ret)[16] = REAL(theta)[15];
      REAL(ret)[17] = REAL(theta)[21];
      REAL(ret)[18] = REAL(theta)[28];
      REAL(ret)[19] = REAL(theta)[36];
      REAL(ret)[20] = REAL(theta)[45];
      REAL(ret)[21] = REAL(theta)[55];
      REAL(ret)[23] = REAL(theta)[3];
      REAL(ret)[34] = REAL(theta)[6];
      REAL(ret)[45] = REAL(theta)[10];
      REAL(ret)[56] = REAL(theta)[15];
      REAL(ret)[67] = REAL(theta)[21];
      REAL(ret)[78] = REAL(theta)[28];
      REAL(ret)[89] = REAL(theta)[36];
      REAL(ret)[100] = REAL(theta)[45];
      REAL(ret)[111] = REAL(theta)[55];
    }
    else if (theta_n == 3){
      REAL(ret)[12] = 2 * exp(2 * REAL(theta)[2]);
      REAL(ret)[13] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[18] = REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[19] = REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[20] = REAL(theta)[46] * exp(REAL(theta)[2]);
      REAL(ret)[21] = REAL(theta)[56] * exp(REAL(theta)[2]);
      REAL(ret)[23] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[34] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[45] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[56] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[67] = REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[78] = REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[89] = REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[100] = REAL(theta)[46] * exp(REAL(theta)[2]);
      REAL(ret)[111] = REAL(theta)[56] * exp(REAL(theta)[2]);
    }
    else if (theta_n == 4){
      REAL(ret)[2] = exp(REAL(theta)[0]);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[22] = exp(REAL(theta)[0]);
      REAL(ret)[23] = REAL(theta)[1];
      REAL(ret)[24] = 2 * REAL(theta)[3];
      REAL(ret)[25] = REAL(theta)[6];
      REAL(ret)[26] = REAL(theta)[10];
      REAL(ret)[27] = REAL(theta)[15];
      REAL(ret)[28] = REAL(theta)[21];
      REAL(ret)[29] = REAL(theta)[28];
      REAL(ret)[30] = REAL(theta)[36];
      REAL(ret)[31] = REAL(theta)[45];
      REAL(ret)[32] = REAL(theta)[55];
      REAL(ret)[35] = REAL(theta)[6];
      REAL(ret)[46] = REAL(theta)[10];
      REAL(ret)[57] = REAL(theta)[15];
      REAL(ret)[68] = REAL(theta)[21];
      REAL(ret)[79] = REAL(theta)[28];
      REAL(ret)[90] = REAL(theta)[36];
      REAL(ret)[101] = REAL(theta)[45];
      REAL(ret)[112] = REAL(theta)[55];
    }
    else if (theta_n == 5){
      REAL(ret)[13] = exp(REAL(theta)[2]);
      REAL(ret)[23] = exp(REAL(theta)[2]);
      REAL(ret)[24] = 2 * REAL(theta)[4];
      REAL(ret)[25] = REAL(theta)[7];
      REAL(ret)[26] = REAL(theta)[11];
      REAL(ret)[27] = REAL(theta)[16];
      REAL(ret)[28] = REAL(theta)[22];
      REAL(ret)[29] = REAL(theta)[29];
      REAL(ret)[30] = REAL(theta)[37];
      REAL(ret)[31] = REAL(theta)[46];
      REAL(ret)[32] = REAL(theta)[56];
      REAL(ret)[35] = REAL(theta)[7];
      REAL(ret)[46] = REAL(theta)[11];
      REAL(ret)[57] = REAL(theta)[16];
      REAL(ret)[68] = REAL(theta)[22];
      REAL(ret)[79] = REAL(theta)[29];
      REAL(ret)[90] = REAL(theta)[37];
      REAL(ret)[101] = REAL(theta)[46];
      REAL(ret)[112] = REAL(theta)[56];
    }
    else if (theta_n == 6){
      REAL(ret)[24] = 2 * exp(2 * REAL(theta)[5]);
      REAL(ret)[25] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[26] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[27] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[28] = REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[29] = REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[30] = REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[31] = REAL(theta)[47] * exp(REAL(theta)[5]);
      REAL(ret)[32] = REAL(theta)[57] * exp(REAL(theta)[5]);
      REAL(ret)[35] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[46] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[57] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[68] = REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[79] = REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[90] = REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[101] = REAL(theta)[47] * exp(REAL(theta)[5]);
      REAL(ret)[112] = REAL(theta)[57] * exp(REAL(theta)[5]);
    }
    else if (theta_n == 7){
      REAL(ret)[3] = exp(REAL(theta)[0]);
      REAL(ret)[14] = REAL(theta)[1];
      REAL(ret)[25] = REAL(theta)[3];
      REAL(ret)[33] = exp(REAL(theta)[0]);
      REAL(ret)[34] = REAL(theta)[1];
      REAL(ret)[35] = REAL(theta)[3];
      REAL(ret)[36] = 2 * REAL(theta)[6];
      REAL(ret)[37] = REAL(theta)[10];
      REAL(ret)[38] = REAL(theta)[15];
      REAL(ret)[39] = REAL(theta)[21];
      REAL(ret)[40] = REAL(theta)[28];
      REAL(ret)[41] = REAL(theta)[36];
      REAL(ret)[42] = REAL(theta)[45];
      REAL(ret)[43] = REAL(theta)[55];
      REAL(ret)[47] = REAL(theta)[10];
      REAL(ret)[58] = REAL(theta)[15];
      REAL(ret)[69] = REAL(theta)[21];
      REAL(ret)[80] = REAL(theta)[28];
      REAL(ret)[91] = REAL(theta)[36];
      REAL(ret)[102] = REAL(theta)[45];
      REAL(ret)[113] = REAL(theta)[55];
    }
    else if (theta_n == 8){
      REAL(ret)[14] = exp(REAL(theta)[2]);
      REAL(ret)[25] = REAL(theta)[4];
      REAL(ret)[34] = exp(REAL(theta)[2]);
      REAL(ret)[35] = REAL(theta)[4];
      REAL(ret)[36] = 2 * REAL(theta)[7];
      REAL(ret)[37] = REAL(theta)[11];
      REAL(ret)[38] = REAL(theta)[16];
      REAL(ret)[39] = REAL(theta)[22];
      REAL(ret)[40] = REAL(theta)[29];
      REAL(ret)[41] = REAL(theta)[37];
      REAL(ret)[42] = REAL(theta)[46];
      REAL(ret)[43] = REAL(theta)[56];
      REAL(ret)[47] = REAL(theta)[11];
      REAL(ret)[58] = REAL(theta)[16];
      REAL(ret)[69] = REAL(theta)[22];
      REAL(ret)[80] = REAL(theta)[29];
      REAL(ret)[91] = REAL(theta)[37];
      REAL(ret)[102] = REAL(theta)[46];
      REAL(ret)[113] = REAL(theta)[56];
    }
    else if (theta_n == 9){
      REAL(ret)[25] = exp(REAL(theta)[5]);
      REAL(ret)[35] = exp(REAL(theta)[5]);
      REAL(ret)[36] = 2 * REAL(theta)[8];
      REAL(ret)[37] = REAL(theta)[12];
      REAL(ret)[38] = REAL(theta)[17];
      REAL(ret)[39] = REAL(theta)[23];
      REAL(ret)[40] = REAL(theta)[30];
      REAL(ret)[41] = REAL(theta)[38];
      REAL(ret)[42] = REAL(theta)[47];
      REAL(ret)[43] = REAL(theta)[57];
      REAL(ret)[47] = REAL(theta)[12];
      REAL(ret)[58] = REAL(theta)[17];
      REAL(ret)[69] = REAL(theta)[23];
      REAL(ret)[80] = REAL(theta)[30];
      REAL(ret)[91] = REAL(theta)[38];
      REAL(ret)[102] = REAL(theta)[47];
      REAL(ret)[113] = REAL(theta)[57];
    }
    else if (theta_n == 10){
      REAL(ret)[36] = 2 * exp(2 * REAL(theta)[9]);
      REAL(ret)[37] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[38] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[39] = REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[40] = REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[41] = REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[42] = REAL(theta)[48] * exp(REAL(theta)[9]);
      REAL(ret)[43] = REAL(theta)[58] * exp(REAL(theta)[9]);
      REAL(ret)[47] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[58] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[69] = REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[80] = REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[91] = REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[102] = REAL(theta)[48] * exp(REAL(theta)[9]);
      REAL(ret)[113] = REAL(theta)[58] * exp(REAL(theta)[9]);
    }
    else if (theta_n == 11){
      REAL(ret)[4] = exp(REAL(theta)[0]);
      REAL(ret)[15] = REAL(theta)[1];
      REAL(ret)[26] = REAL(theta)[3];
      REAL(ret)[37] = REAL(theta)[6];
      REAL(ret)[44] = exp(REAL(theta)[0]);
      REAL(ret)[45] = REAL(theta)[1];
      REAL(ret)[46] = REAL(theta)[3];
      REAL(ret)[47] = REAL(theta)[6];
      REAL(ret)[48] = 2 * REAL(theta)[10];
      REAL(ret)[49] = REAL(theta)[15];
      REAL(ret)[50] = REAL(theta)[21];
      REAL(ret)[51] = REAL(theta)[28];
      REAL(ret)[52] = REAL(theta)[36];
      REAL(ret)[53] = REAL(theta)[45];
      REAL(ret)[54] = REAL(theta)[55];
      REAL(ret)[59] = REAL(theta)[15];
      REAL(ret)[70] = REAL(theta)[21];
      REAL(ret)[81] = REAL(theta)[28];
      REAL(ret)[92] = REAL(theta)[36];
      REAL(ret)[103] = REAL(theta)[45];
      REAL(ret)[114] = REAL(theta)[55];
    }
    else if (theta_n == 12){
      REAL(ret)[15] = exp(REAL(theta)[2]);
      REAL(ret)[26] = REAL(theta)[4];
      REAL(ret)[37] = REAL(theta)[7];
      REAL(ret)[45] = exp(REAL(theta)[2]);
      REAL(ret)[46] = REAL(theta)[4];
      REAL(ret)[47] = REAL(theta)[7];
      REAL(ret)[48] = 2 * REAL(theta)[11];
      REAL(ret)[49] = REAL(theta)[16];
      REAL(ret)[50] = REAL(theta)[22];
      REAL(ret)[51] = REAL(theta)[29];
      REAL(ret)[52] = REAL(theta)[37];
      REAL(ret)[53] = REAL(theta)[46];
      REAL(ret)[54] = REAL(theta)[56];
      REAL(ret)[59] = REAL(theta)[16];
      REAL(ret)[70] = REAL(theta)[22];
      REAL(ret)[81] = REAL(theta)[29];
      REAL(ret)[92] = REAL(theta)[37];
      REAL(ret)[103] = REAL(theta)[46];
      REAL(ret)[114] = REAL(theta)[56];
    }
    else if (theta_n == 13){
      REAL(ret)[26] = exp(REAL(theta)[5]);
      REAL(ret)[37] = REAL(theta)[8];
      REAL(ret)[46] = exp(REAL(theta)[5]);
      REAL(ret)[47] = REAL(theta)[8];
      REAL(ret)[48] = 2 * REAL(theta)[12];
      REAL(ret)[49] = REAL(theta)[17];
      REAL(ret)[50] = REAL(theta)[23];
      REAL(ret)[51] = REAL(theta)[30];
      REAL(ret)[52] = REAL(theta)[38];
      REAL(ret)[53] = REAL(theta)[47];
      REAL(ret)[54] = REAL(theta)[57];
      REAL(ret)[59] = REAL(theta)[17];
      REAL(ret)[70] = REAL(theta)[23];
      REAL(ret)[81] = REAL(theta)[30];
      REAL(ret)[92] = REAL(theta)[38];
      REAL(ret)[103] = REAL(theta)[47];
      REAL(ret)[114] = REAL(theta)[57];
    }
    else if (theta_n == 14){
      REAL(ret)[37] = exp(REAL(theta)[9]);
      REAL(ret)[47] = exp(REAL(theta)[9]);
      REAL(ret)[48] = 2 * REAL(theta)[13];
      REAL(ret)[49] = REAL(theta)[18];
      REAL(ret)[50] = REAL(theta)[24];
      REAL(ret)[51] = REAL(theta)[31];
      REAL(ret)[52] = REAL(theta)[39];
      REAL(ret)[53] = REAL(theta)[48];
      REAL(ret)[54] = REAL(theta)[58];
      REAL(ret)[59] = REAL(theta)[18];
      REAL(ret)[70] = REAL(theta)[24];
      REAL(ret)[81] = REAL(theta)[31];
      REAL(ret)[92] = REAL(theta)[39];
      REAL(ret)[103] = REAL(theta)[48];
      REAL(ret)[114] = REAL(theta)[58];
    }
    else if (theta_n == 15){
      REAL(ret)[48] = 2 * exp(2 * REAL(theta)[14]);
      REAL(ret)[49] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[50] = REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[51] = REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[52] = REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[53] = REAL(theta)[49] * exp(REAL(theta)[14]);
      REAL(ret)[54] = REAL(theta)[59] * exp(REAL(theta)[14]);
      REAL(ret)[59] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[70] = REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[81] = REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[92] = REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[103] = REAL(theta)[49] * exp(REAL(theta)[14]);
      REAL(ret)[114] = REAL(theta)[59] * exp(REAL(theta)[14]);
    }
    else if (theta_n == 16){
      REAL(ret)[5] = exp(REAL(theta)[0]);
      REAL(ret)[16] = REAL(theta)[1];
      REAL(ret)[27] = REAL(theta)[3];
      REAL(ret)[38] = REAL(theta)[6];
      REAL(ret)[49] = REAL(theta)[10];
      REAL(ret)[55] = exp(REAL(theta)[0]);
      REAL(ret)[56] = REAL(theta)[1];
      REAL(ret)[57] = REAL(theta)[3];
      REAL(ret)[58] = REAL(theta)[6];
      REAL(ret)[59] = REAL(theta)[10];
      REAL(ret)[60] = 2 * REAL(theta)[15];
      REAL(ret)[61] = REAL(theta)[21];
      REAL(ret)[62] = REAL(theta)[28];
      REAL(ret)[63] = REAL(theta)[36];
      REAL(ret)[64] = REAL(theta)[45];
      REAL(ret)[65] = REAL(theta)[55];
      REAL(ret)[71] = REAL(theta)[21];
      REAL(ret)[82] = REAL(theta)[28];
      REAL(ret)[93] = REAL(theta)[36];
      REAL(ret)[104] = REAL(theta)[45];
      REAL(ret)[115] = REAL(theta)[55];
    }
    else if (theta_n == 17){
      REAL(ret)[16] = exp(REAL(theta)[2]);
      REAL(ret)[27] = REAL(theta)[4];
      REAL(ret)[38] = REAL(theta)[7];
      REAL(ret)[49] = REAL(theta)[11];
      REAL(ret)[56] = exp(REAL(theta)[2]);
      REAL(ret)[57] = REAL(theta)[4];
      REAL(ret)[58] = REAL(theta)[7];
      REAL(ret)[59] = REAL(theta)[11];
      REAL(ret)[60] = 2 * REAL(theta)[16];
      REAL(ret)[61] = REAL(theta)[22];
      REAL(ret)[62] = REAL(theta)[29];
      REAL(ret)[63] = REAL(theta)[37];
      REAL(ret)[64] = REAL(theta)[46];
      REAL(ret)[65] = REAL(theta)[56];
      REAL(ret)[71] = REAL(theta)[22];
      REAL(ret)[82] = REAL(theta)[29];
      REAL(ret)[93] = REAL(theta)[37];
      REAL(ret)[104] = REAL(theta)[46];
      REAL(ret)[115] = REAL(theta)[56];
    }
    else if (theta_n == 18){
      REAL(ret)[27] = exp(REAL(theta)[5]);
      REAL(ret)[38] = REAL(theta)[8];
      REAL(ret)[49] = REAL(theta)[12];
      REAL(ret)[57] = exp(REAL(theta)[5]);
      REAL(ret)[58] = REAL(theta)[8];
      REAL(ret)[59] = REAL(theta)[12];
      REAL(ret)[60] = 2 * REAL(theta)[17];
      REAL(ret)[61] = REAL(theta)[23];
      REAL(ret)[62] = REAL(theta)[30];
      REAL(ret)[63] = REAL(theta)[38];
      REAL(ret)[64] = REAL(theta)[47];
      REAL(ret)[65] = REAL(theta)[57];
      REAL(ret)[71] = REAL(theta)[23];
      REAL(ret)[82] = REAL(theta)[30];
      REAL(ret)[93] = REAL(theta)[38];
      REAL(ret)[104] = REAL(theta)[47];
      REAL(ret)[115] = REAL(theta)[57];
    }
    else if (theta_n == 19){
      REAL(ret)[38] = exp(REAL(theta)[9]);
      REAL(ret)[49] = REAL(theta)[13];
      REAL(ret)[58] = exp(REAL(theta)[9]);
      REAL(ret)[59] = REAL(theta)[13];
      REAL(ret)[60] = 2 * REAL(theta)[18];
      REAL(ret)[61] = REAL(theta)[24];
      REAL(ret)[62] = REAL(theta)[31];
      REAL(ret)[63] = REAL(theta)[39];
      REAL(ret)[64] = REAL(theta)[48];
      REAL(ret)[65] = REAL(theta)[58];
      REAL(ret)[71] = REAL(theta)[24];
      REAL(ret)[82] = REAL(theta)[31];
      REAL(ret)[93] = REAL(theta)[39];
      REAL(ret)[104] = REAL(theta)[48];
      REAL(ret)[115] = REAL(theta)[58];
    }
    else if (theta_n == 20){
      REAL(ret)[49] = exp(REAL(theta)[14]);
      REAL(ret)[59] = exp(REAL(theta)[14]);
      REAL(ret)[60] = 2 * REAL(theta)[19];
      REAL(ret)[61] = REAL(theta)[25];
      REAL(ret)[62] = REAL(theta)[32];
      REAL(ret)[63] = REAL(theta)[40];
      REAL(ret)[64] = REAL(theta)[49];
      REAL(ret)[65] = REAL(theta)[59];
      REAL(ret)[71] = REAL(theta)[25];
      REAL(ret)[82] = REAL(theta)[32];
      REAL(ret)[93] = REAL(theta)[40];
      REAL(ret)[104] = REAL(theta)[49];
      REAL(ret)[115] = REAL(theta)[59];
    }
    else if (theta_n == 21){
      REAL(ret)[60] = 2 * exp(2 * REAL(theta)[20]);
      REAL(ret)[61] = REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[62] = REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[63] = REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[64] = REAL(theta)[50] * exp(REAL(theta)[20]);
      REAL(ret)[65] = REAL(theta)[60] * exp(REAL(theta)[20]);
      REAL(ret)[71] = REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[82] = REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[93] = REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[104] = REAL(theta)[50] * exp(REAL(theta)[20]);
      REAL(ret)[115] = REAL(theta)[60] * exp(REAL(theta)[20]);
    }
    else if (theta_n == 22){
      REAL(ret)[6] = exp(REAL(theta)[0]);
      REAL(ret)[17] = REAL(theta)[1];
      REAL(ret)[28] = REAL(theta)[3];
      REAL(ret)[39] = REAL(theta)[6];
      REAL(ret)[50] = REAL(theta)[10];
      REAL(ret)[61] = REAL(theta)[15];
      REAL(ret)[66] = exp(REAL(theta)[0]);
      REAL(ret)[67] = REAL(theta)[1];
      REAL(ret)[68] = REAL(theta)[3];
      REAL(ret)[69] = REAL(theta)[6];
      REAL(ret)[70] = REAL(theta)[10];
      REAL(ret)[71] = REAL(theta)[15];
      REAL(ret)[72] = 2 * REAL(theta)[21];
      REAL(ret)[73] = REAL(theta)[28];
      REAL(ret)[74] = REAL(theta)[36];
      REAL(ret)[75] = REAL(theta)[45];
      REAL(ret)[76] = REAL(theta)[55];
      REAL(ret)[83] = REAL(theta)[28];
      REAL(ret)[94] = REAL(theta)[36];
      REAL(ret)[105] = REAL(theta)[45];
      REAL(ret)[116] = REAL(theta)[55];
    }
    else if (theta_n == 23){
      REAL(ret)[17] = exp(REAL(theta)[2]);
      REAL(ret)[28] = REAL(theta)[4];
      REAL(ret)[39] = REAL(theta)[7];
      REAL(ret)[50] = REAL(theta)[11];
      REAL(ret)[61] = REAL(theta)[16];
      REAL(ret)[67] = exp(REAL(theta)[2]);
      REAL(ret)[68] = REAL(theta)[4];
      REAL(ret)[69] = REAL(theta)[7];
      REAL(ret)[70] = REAL(theta)[11];
      REAL(ret)[71] = REAL(theta)[16];
      REAL(ret)[72] = 2 * REAL(theta)[22];
      REAL(ret)[73] = REAL(theta)[29];
      REAL(ret)[74] = REAL(theta)[37];
      REAL(ret)[75] = REAL(theta)[46];
      REAL(ret)[76] = REAL(theta)[56];
      REAL(ret)[83] = REAL(theta)[29];
      REAL(ret)[94] = REAL(theta)[37];
      REAL(ret)[105] = REAL(theta)[46];
      REAL(ret)[116] = REAL(theta)[56];
    }
    else if (theta_n == 24){
      REAL(ret)[28] = exp(REAL(theta)[5]);
      REAL(ret)[39] = REAL(theta)[8];
      REAL(ret)[50] = REAL(theta)[12];
      REAL(ret)[61] = REAL(theta)[17];
      REAL(ret)[68] = exp(REAL(theta)[5]);
      REAL(ret)[69] = REAL(theta)[8];
      REAL(ret)[70] = REAL(theta)[12];
      REAL(ret)[71] = REAL(theta)[17];
      REAL(ret)[72] = 2 * REAL(theta)[23];
      REAL(ret)[73] = REAL(theta)[30];
      REAL(ret)[74] = REAL(theta)[38];
      REAL(ret)[75] = REAL(theta)[47];
      REAL(ret)[76] = REAL(theta)[57];
      REAL(ret)[83] = REAL(theta)[30];
      REAL(ret)[94] = REAL(theta)[38];
      REAL(ret)[105] = REAL(theta)[47];
      REAL(ret)[116] = REAL(theta)[57];
    }
    else if (theta_n == 25){
      REAL(ret)[39] = exp(REAL(theta)[9]);
      REAL(ret)[50] = REAL(theta)[13];
      REAL(ret)[61] = REAL(theta)[18];
      REAL(ret)[69] = exp(REAL(theta)[9]);
      REAL(ret)[70] = REAL(theta)[13];
      REAL(ret)[71] = REAL(theta)[18];
      REAL(ret)[72] = 2 * REAL(theta)[24];
      REAL(ret)[73] = REAL(theta)[31];
      REAL(ret)[74] = REAL(theta)[39];
      REAL(ret)[75] = REAL(theta)[48];
      REAL(ret)[76] = REAL(theta)[58];
      REAL(ret)[83] = REAL(theta)[31];
      REAL(ret)[94] = REAL(theta)[39];
      REAL(ret)[105] = REAL(theta)[48];
      REAL(ret)[116] = REAL(theta)[58];
    }
    else if (theta_n == 26){
      REAL(ret)[50] = exp(REAL(theta)[14]);
      REAL(ret)[61] = REAL(theta)[19];
      REAL(ret)[70] = exp(REAL(theta)[14]);
      REAL(ret)[71] = REAL(theta)[19];
      REAL(ret)[72] = 2 * REAL(theta)[25];
      REAL(ret)[73] = REAL(theta)[32];
      REAL(ret)[74] = REAL(theta)[40];
      REAL(ret)[75] = REAL(theta)[49];
      REAL(ret)[76] = REAL(theta)[59];
      REAL(ret)[83] = REAL(theta)[32];
      REAL(ret)[94] = REAL(theta)[40];
      REAL(ret)[105] = REAL(theta)[49];
      REAL(ret)[116] = REAL(theta)[59];
    }
    else if (theta_n == 27){
      REAL(ret)[61] = exp(REAL(theta)[20]);
      REAL(ret)[71] = exp(REAL(theta)[20]);
      REAL(ret)[72] = 2 * REAL(theta)[26];
      REAL(ret)[73] = REAL(theta)[33];
      REAL(ret)[74] = REAL(theta)[41];
      REAL(ret)[75] = REAL(theta)[50];
      REAL(ret)[76] = REAL(theta)[60];
      REAL(ret)[83] = REAL(theta)[33];
      REAL(ret)[94] = REAL(theta)[41];
      REAL(ret)[105] = REAL(theta)[50];
      REAL(ret)[116] = REAL(theta)[60];
    }
    else if (theta_n == 28){
      REAL(ret)[72] = 2 * exp(2 * REAL(theta)[27]);
      REAL(ret)[73] = REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[74] = REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[75] = REAL(theta)[51] * exp(REAL(theta)[27]);
      REAL(ret)[76] = REAL(theta)[61] * exp(REAL(theta)[27]);
      REAL(ret)[83] = REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[94] = REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[105] = REAL(theta)[51] * exp(REAL(theta)[27]);
      REAL(ret)[116] = REAL(theta)[61] * exp(REAL(theta)[27]);
    }
    else if (theta_n == 29){
      REAL(ret)[7] = exp(REAL(theta)[0]);
      REAL(ret)[18] = REAL(theta)[1];
      REAL(ret)[29] = REAL(theta)[3];
      REAL(ret)[40] = REAL(theta)[6];
      REAL(ret)[51] = REAL(theta)[10];
      REAL(ret)[62] = REAL(theta)[15];
      REAL(ret)[73] = REAL(theta)[21];
      REAL(ret)[77] = exp(REAL(theta)[0]);
      REAL(ret)[78] = REAL(theta)[1];
      REAL(ret)[79] = REAL(theta)[3];
      REAL(ret)[80] = REAL(theta)[6];
      REAL(ret)[81] = REAL(theta)[10];
      REAL(ret)[82] = REAL(theta)[15];
      REAL(ret)[83] = REAL(theta)[21];
      REAL(ret)[84] = 2 * REAL(theta)[28];
      REAL(ret)[85] = REAL(theta)[36];
      REAL(ret)[86] = REAL(theta)[45];
      REAL(ret)[87] = REAL(theta)[55];
      REAL(ret)[95] = REAL(theta)[36];
      REAL(ret)[106] = REAL(theta)[45];
      REAL(ret)[117] = REAL(theta)[55];
    }
    else if (theta_n == 30){
      REAL(ret)[18] = exp(REAL(theta)[2]);
      REAL(ret)[29] = REAL(theta)[4];
      REAL(ret)[40] = REAL(theta)[7];
      REAL(ret)[51] = REAL(theta)[11];
      REAL(ret)[62] = REAL(theta)[16];
      REAL(ret)[73] = REAL(theta)[22];
      REAL(ret)[78] = exp(REAL(theta)[2]);
      REAL(ret)[79] = REAL(theta)[4];
      REAL(ret)[80] = REAL(theta)[7];
      REAL(ret)[81] = REAL(theta)[11];
      REAL(ret)[82] = REAL(theta)[16];
      REAL(ret)[83] = REAL(theta)[22];
      REAL(ret)[84] = 2 * REAL(theta)[29];
      REAL(ret)[85] = REAL(theta)[37];
      REAL(ret)[86] = REAL(theta)[46];
      REAL(ret)[87] = REAL(theta)[56];
      REAL(ret)[95] = REAL(theta)[37];
      REAL(ret)[106] = REAL(theta)[46];
      REAL(ret)[117] = REAL(theta)[56];
    }
    else if (theta_n == 31){
      REAL(ret)[29] = exp(REAL(theta)[5]);
      REAL(ret)[40] = REAL(theta)[8];
      REAL(ret)[51] = REAL(theta)[12];
      REAL(ret)[62] = REAL(theta)[17];
      REAL(ret)[73] = REAL(theta)[23];
      REAL(ret)[79] = exp(REAL(theta)[5]);
      REAL(ret)[80] = REAL(theta)[8];
      REAL(ret)[81] = REAL(theta)[12];
      REAL(ret)[82] = REAL(theta)[17];
      REAL(ret)[83] = REAL(theta)[23];
      REAL(ret)[84] = 2 * REAL(theta)[30];
      REAL(ret)[85] = REAL(theta)[38];
      REAL(ret)[86] = REAL(theta)[47];
      REAL(ret)[87] = REAL(theta)[57];
      REAL(ret)[95] = REAL(theta)[38];
      REAL(ret)[106] = REAL(theta)[47];
      REAL(ret)[117] = REAL(theta)[57];
    }
    else if (theta_n == 32){
      REAL(ret)[40] = exp(REAL(theta)[9]);
      REAL(ret)[51] = REAL(theta)[13];
      REAL(ret)[62] = REAL(theta)[18];
      REAL(ret)[73] = REAL(theta)[24];
      REAL(ret)[80] = exp(REAL(theta)[9]);
      REAL(ret)[81] = REAL(theta)[13];
      REAL(ret)[82] = REAL(theta)[18];
      REAL(ret)[83] = REAL(theta)[24];
      REAL(ret)[84] = 2 * REAL(theta)[31];
      REAL(ret)[85] = REAL(theta)[39];
      REAL(ret)[86] = REAL(theta)[48];
      REAL(ret)[87] = REAL(theta)[58];
      REAL(ret)[95] = REAL(theta)[39];
      REAL(ret)[106] = REAL(theta)[48];
      REAL(ret)[117] = REAL(theta)[58];
    }
    else if (theta_n == 33){
      REAL(ret)[51] = exp(REAL(theta)[14]);
      REAL(ret)[62] = REAL(theta)[19];
      REAL(ret)[73] = REAL(theta)[25];
      REAL(ret)[81] = exp(REAL(theta)[14]);
      REAL(ret)[82] = REAL(theta)[19];
      REAL(ret)[83] = REAL(theta)[25];
      REAL(ret)[84] = 2 * REAL(theta)[32];
      REAL(ret)[85] = REAL(theta)[40];
      REAL(ret)[86] = REAL(theta)[49];
      REAL(ret)[87] = REAL(theta)[59];
      REAL(ret)[95] = REAL(theta)[40];
      REAL(ret)[106] = REAL(theta)[49];
      REAL(ret)[117] = REAL(theta)[59];
    }
    else if (theta_n == 34){
      REAL(ret)[62] = exp(REAL(theta)[20]);
      REAL(ret)[73] = REAL(theta)[26];
      REAL(ret)[82] = exp(REAL(theta)[20]);
      REAL(ret)[83] = REAL(theta)[26];
      REAL(ret)[84] = 2 * REAL(theta)[33];
      REAL(ret)[85] = REAL(theta)[41];
      REAL(ret)[86] = REAL(theta)[50];
      REAL(ret)[87] = REAL(theta)[60];
      REAL(ret)[95] = REAL(theta)[41];
      REAL(ret)[106] = REAL(theta)[50];
      REAL(ret)[117] = REAL(theta)[60];
    }
    else if (theta_n == 35){
      REAL(ret)[73] = exp(REAL(theta)[27]);
      REAL(ret)[83] = exp(REAL(theta)[27]);
      REAL(ret)[84] = 2 * REAL(theta)[34];
      REAL(ret)[85] = REAL(theta)[42];
      REAL(ret)[86] = REAL(theta)[51];
      REAL(ret)[87] = REAL(theta)[61];
      REAL(ret)[95] = REAL(theta)[42];
      REAL(ret)[106] = REAL(theta)[51];
      REAL(ret)[117] = REAL(theta)[61];
    }
    else if (theta_n == 36){
      REAL(ret)[84] = 2 * exp(2 * REAL(theta)[35]);
      REAL(ret)[85] = REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[86] = REAL(theta)[52] * exp(REAL(theta)[35]);
      REAL(ret)[87] = REAL(theta)[62] * exp(REAL(theta)[35]);
      REAL(ret)[95] = REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[106] = REAL(theta)[52] * exp(REAL(theta)[35]);
      REAL(ret)[117] = REAL(theta)[62] * exp(REAL(theta)[35]);
    }
    else if (theta_n == 37){
      REAL(ret)[8] = exp(REAL(theta)[0]);
      REAL(ret)[19] = REAL(theta)[1];
      REAL(ret)[30] = REAL(theta)[3];
      REAL(ret)[41] = REAL(theta)[6];
      REAL(ret)[52] = REAL(theta)[10];
      REAL(ret)[63] = REAL(theta)[15];
      REAL(ret)[74] = REAL(theta)[21];
      REAL(ret)[85] = REAL(theta)[28];
      REAL(ret)[88] = exp(REAL(theta)[0]);
      REAL(ret)[89] = REAL(theta)[1];
      REAL(ret)[90] = REAL(theta)[3];
      REAL(ret)[91] = REAL(theta)[6];
      REAL(ret)[92] = REAL(theta)[10];
      REAL(ret)[93] = REAL(theta)[15];
      REAL(ret)[94] = REAL(theta)[21];
      REAL(ret)[95] = REAL(theta)[28];
      REAL(ret)[96] = 2 * REAL(theta)[36];
      REAL(ret)[97] = REAL(theta)[45];
      REAL(ret)[98] = REAL(theta)[55];
      REAL(ret)[107] = REAL(theta)[45];
      REAL(ret)[118] = REAL(theta)[55];
    }
    else if (theta_n == 38){
      REAL(ret)[19] = exp(REAL(theta)[2]);
      REAL(ret)[30] = REAL(theta)[4];
      REAL(ret)[41] = REAL(theta)[7];
      REAL(ret)[52] = REAL(theta)[11];
      REAL(ret)[63] = REAL(theta)[16];
      REAL(ret)[74] = REAL(theta)[22];
      REAL(ret)[85] = REAL(theta)[29];
      REAL(ret)[89] = exp(REAL(theta)[2]);
      REAL(ret)[90] = REAL(theta)[4];
      REAL(ret)[91] = REAL(theta)[7];
      REAL(ret)[92] = REAL(theta)[11];
      REAL(ret)[93] = REAL(theta)[16];
      REAL(ret)[94] = REAL(theta)[22];
      REAL(ret)[95] = REAL(theta)[29];
      REAL(ret)[96] = 2 * REAL(theta)[37];
      REAL(ret)[97] = REAL(theta)[46];
      REAL(ret)[98] = REAL(theta)[56];
      REAL(ret)[107] = REAL(theta)[46];
      REAL(ret)[118] = REAL(theta)[56];
    }
    else if (theta_n == 39){
      REAL(ret)[30] = exp(REAL(theta)[5]);
      REAL(ret)[41] = REAL(theta)[8];
      REAL(ret)[52] = REAL(theta)[12];
      REAL(ret)[63] = REAL(theta)[17];
      REAL(ret)[74] = REAL(theta)[23];
      REAL(ret)[85] = REAL(theta)[30];
      REAL(ret)[90] = exp(REAL(theta)[5]);
      REAL(ret)[91] = REAL(theta)[8];
      REAL(ret)[92] = REAL(theta)[12];
      REAL(ret)[93] = REAL(theta)[17];
      REAL(ret)[94] = REAL(theta)[23];
      REAL(ret)[95] = REAL(theta)[30];
      REAL(ret)[96] = 2 * REAL(theta)[38];
      REAL(ret)[97] = REAL(theta)[47];
      REAL(ret)[98] = REAL(theta)[57];
      REAL(ret)[107] = REAL(theta)[47];
      REAL(ret)[118] = REAL(theta)[57];
    }
    else if (theta_n == 40){
      REAL(ret)[41] = exp(REAL(theta)[9]);
      REAL(ret)[52] = REAL(theta)[13];
      REAL(ret)[63] = REAL(theta)[18];
      REAL(ret)[74] = REAL(theta)[24];
      REAL(ret)[85] = REAL(theta)[31];
      REAL(ret)[91] = exp(REAL(theta)[9]);
      REAL(ret)[92] = REAL(theta)[13];
      REAL(ret)[93] = REAL(theta)[18];
      REAL(ret)[94] = REAL(theta)[24];
      REAL(ret)[95] = REAL(theta)[31];
      REAL(ret)[96] = 2 * REAL(theta)[39];
      REAL(ret)[97] = REAL(theta)[48];
      REAL(ret)[98] = REAL(theta)[58];
      REAL(ret)[107] = REAL(theta)[48];
      REAL(ret)[118] = REAL(theta)[58];
    }
    else if (theta_n == 41){
      REAL(ret)[52] = exp(REAL(theta)[14]);
      REAL(ret)[63] = REAL(theta)[19];
      REAL(ret)[74] = REAL(theta)[25];
      REAL(ret)[85] = REAL(theta)[32];
      REAL(ret)[92] = exp(REAL(theta)[14]);
      REAL(ret)[93] = REAL(theta)[19];
      REAL(ret)[94] = REAL(theta)[25];
      REAL(ret)[95] = REAL(theta)[32];
      REAL(ret)[96] = 2 * REAL(theta)[40];
      REAL(ret)[97] = REAL(theta)[49];
      REAL(ret)[98] = REAL(theta)[59];
      REAL(ret)[107] = REAL(theta)[49];
      REAL(ret)[118] = REAL(theta)[59];
    }
    else if (theta_n == 42){
      REAL(ret)[63] = exp(REAL(theta)[20]);
      REAL(ret)[74] = REAL(theta)[26];
      REAL(ret)[85] = REAL(theta)[33];
      REAL(ret)[93] = exp(REAL(theta)[20]);
      REAL(ret)[94] = REAL(theta)[26];
      REAL(ret)[95] = REAL(theta)[33];
      REAL(ret)[96] = 2 * REAL(theta)[41];
      REAL(ret)[97] = REAL(theta)[50];
      REAL(ret)[98] = REAL(theta)[60];
      REAL(ret)[107] = REAL(theta)[50];
      REAL(ret)[118] = REAL(theta)[60];
    }
    else if (theta_n == 43){
      REAL(ret)[74] = exp(REAL(theta)[27]);
      REAL(ret)[85] = REAL(theta)[34];
      REAL(ret)[94] = exp(REAL(theta)[27]);
      REAL(ret)[95] = REAL(theta)[34];
      REAL(ret)[96] = 2 * REAL(theta)[42];
      REAL(ret)[97] = REAL(theta)[51];
      REAL(ret)[98] = REAL(theta)[61];
      REAL(ret)[107] = REAL(theta)[51];
      REAL(ret)[118] = REAL(theta)[61];
    }
    else if (theta_n == 44){
      REAL(ret)[85] = exp(REAL(theta)[35]);
      REAL(ret)[95] = exp(REAL(theta)[35]);
      REAL(ret)[96] = 2 * REAL(theta)[43];
      REAL(ret)[97] = REAL(theta)[52];
      REAL(ret)[98] = REAL(theta)[62];
      REAL(ret)[107] = REAL(theta)[52];
      REAL(ret)[118] = REAL(theta)[62];
    }
    else if (theta_n == 45){
      REAL(ret)[96] = 2 * exp(2 * REAL(theta)[44]);
      REAL(ret)[97] = REAL(theta)[53] * exp(REAL(theta)[44]);
      REAL(ret)[98] = REAL(theta)[63] * exp(REAL(theta)[44]);
      REAL(ret)[107] = REAL(theta)[53] * exp(REAL(theta)[44]);
      REAL(ret)[118] = REAL(theta)[63] * exp(REAL(theta)[44]);
    }
    else if (theta_n == 46){
      REAL(ret)[9] = exp(REAL(theta)[0]);
      REAL(ret)[20] = REAL(theta)[1];
      REAL(ret)[31] = REAL(theta)[3];
      REAL(ret)[42] = REAL(theta)[6];
      REAL(ret)[53] = REAL(theta)[10];
      REAL(ret)[64] = REAL(theta)[15];
      REAL(ret)[75] = REAL(theta)[21];
      REAL(ret)[86] = REAL(theta)[28];
      REAL(ret)[97] = REAL(theta)[36];
      REAL(ret)[99] = exp(REAL(theta)[0]);
      REAL(ret)[100] = REAL(theta)[1];
      REAL(ret)[101] = REAL(theta)[3];
      REAL(ret)[102] = REAL(theta)[6];
      REAL(ret)[103] = REAL(theta)[10];
      REAL(ret)[104] = REAL(theta)[15];
      REAL(ret)[105] = REAL(theta)[21];
      REAL(ret)[106] = REAL(theta)[28];
      REAL(ret)[107] = REAL(theta)[36];
      REAL(ret)[108] = 2 * REAL(theta)[45];
      REAL(ret)[109] = REAL(theta)[55];
      REAL(ret)[119] = REAL(theta)[55];
    }
    else if (theta_n == 47){
      REAL(ret)[20] = exp(REAL(theta)[2]);
      REAL(ret)[31] = REAL(theta)[4];
      REAL(ret)[42] = REAL(theta)[7];
      REAL(ret)[53] = REAL(theta)[11];
      REAL(ret)[64] = REAL(theta)[16];
      REAL(ret)[75] = REAL(theta)[22];
      REAL(ret)[86] = REAL(theta)[29];
      REAL(ret)[97] = REAL(theta)[37];
      REAL(ret)[100] = exp(REAL(theta)[2]);
      REAL(ret)[101] = REAL(theta)[4];
      REAL(ret)[102] = REAL(theta)[7];
      REAL(ret)[103] = REAL(theta)[11];
      REAL(ret)[104] = REAL(theta)[16];
      REAL(ret)[105] = REAL(theta)[22];
      REAL(ret)[106] = REAL(theta)[29];
      REAL(ret)[107] = REAL(theta)[37];
      REAL(ret)[108] = 2 * REAL(theta)[46];
      REAL(ret)[109] = REAL(theta)[56];
      REAL(ret)[119] = REAL(theta)[56];
    }
    else if (theta_n == 48){
      REAL(ret)[31] = exp(REAL(theta)[5]);
      REAL(ret)[42] = REAL(theta)[8];
      REAL(ret)[53] = REAL(theta)[12];
      REAL(ret)[64] = REAL(theta)[17];
      REAL(ret)[75] = REAL(theta)[23];
      REAL(ret)[86] = REAL(theta)[30];
      REAL(ret)[97] = REAL(theta)[38];
      REAL(ret)[101] = exp(REAL(theta)[5]);
      REAL(ret)[102] = REAL(theta)[8];
      REAL(ret)[103] = REAL(theta)[12];
      REAL(ret)[104] = REAL(theta)[17];
      REAL(ret)[105] = REAL(theta)[23];
      REAL(ret)[106] = REAL(theta)[30];
      REAL(ret)[107] = REAL(theta)[38];
      REAL(ret)[108] = 2 * REAL(theta)[47];
      REAL(ret)[109] = REAL(theta)[57];
      REAL(ret)[119] = REAL(theta)[57];
    }
    else if (theta_n == 49){
      REAL(ret)[42] = exp(REAL(theta)[9]);
      REAL(ret)[53] = REAL(theta)[13];
      REAL(ret)[64] = REAL(theta)[18];
      REAL(ret)[75] = REAL(theta)[24];
      REAL(ret)[86] = REAL(theta)[31];
      REAL(ret)[97] = REAL(theta)[39];
      REAL(ret)[102] = exp(REAL(theta)[9]);
      REAL(ret)[103] = REAL(theta)[13];
      REAL(ret)[104] = REAL(theta)[18];
      REAL(ret)[105] = REAL(theta)[24];
      REAL(ret)[106] = REAL(theta)[31];
      REAL(ret)[107] = REAL(theta)[39];
      REAL(ret)[108] = 2 * REAL(theta)[48];
      REAL(ret)[109] = REAL(theta)[58];
      REAL(ret)[119] = REAL(theta)[58];
    }
    else if (theta_n == 50){
      REAL(ret)[53] = exp(REAL(theta)[14]);
      REAL(ret)[64] = REAL(theta)[19];
      REAL(ret)[75] = REAL(theta)[25];
      REAL(ret)[86] = REAL(theta)[32];
      REAL(ret)[97] = REAL(theta)[40];
      REAL(ret)[103] = exp(REAL(theta)[14]);
      REAL(ret)[104] = REAL(theta)[19];
      REAL(ret)[105] = REAL(theta)[25];
      REAL(ret)[106] = REAL(theta)[32];
      REAL(ret)[107] = REAL(theta)[40];
      REAL(ret)[108] = 2 * REAL(theta)[49];
      REAL(ret)[109] = REAL(theta)[59];
      REAL(ret)[119] = REAL(theta)[59];
    }
    else if (theta_n == 51){
      REAL(ret)[64] = exp(REAL(theta)[20]);
      REAL(ret)[75] = REAL(theta)[26];
      REAL(ret)[86] = REAL(theta)[33];
      REAL(ret)[97] = REAL(theta)[41];
      REAL(ret)[104] = exp(REAL(theta)[20]);
      REAL(ret)[105] = REAL(theta)[26];
      REAL(ret)[106] = REAL(theta)[33];
      REAL(ret)[107] = REAL(theta)[41];
      REAL(ret)[108] = 2 * REAL(theta)[50];
      REAL(ret)[109] = REAL(theta)[60];
      REAL(ret)[119] = REAL(theta)[60];
    }
    else if (theta_n == 52){
      REAL(ret)[75] = exp(REAL(theta)[27]);
      REAL(ret)[86] = REAL(theta)[34];
      REAL(ret)[97] = REAL(theta)[42];
      REAL(ret)[105] = exp(REAL(theta)[27]);
      REAL(ret)[106] = REAL(theta)[34];
      REAL(ret)[107] = REAL(theta)[42];
      REAL(ret)[108] = 2 * REAL(theta)[51];
      REAL(ret)[109] = REAL(theta)[61];
      REAL(ret)[119] = REAL(theta)[61];
    }
    else if (theta_n == 53){
      REAL(ret)[86] = exp(REAL(theta)[35]);
      REAL(ret)[97] = REAL(theta)[43];
      REAL(ret)[106] = exp(REAL(theta)[35]);
      REAL(ret)[107] = REAL(theta)[43];
      REAL(ret)[108] = 2 * REAL(theta)[52];
      REAL(ret)[109] = REAL(theta)[62];
      REAL(ret)[119] = REAL(theta)[62];
    }
    else if (theta_n == 54){
      REAL(ret)[97] = exp(REAL(theta)[44]);
      REAL(ret)[107] = exp(REAL(theta)[44]);
      REAL(ret)[108] = 2 * REAL(theta)[53];
      REAL(ret)[109] = REAL(theta)[63];
      REAL(ret)[119] = REAL(theta)[63];
    }
    else if (theta_n == 55){
      REAL(ret)[108] = 2 * exp(2 * REAL(theta)[54]);
      REAL(ret)[109] = REAL(theta)[64] * exp(REAL(theta)[54]);
      REAL(ret)[119] = REAL(theta)[64] * exp(REAL(theta)[54]);
    }
    else if (theta_n == 56){
      REAL(ret)[10] = exp(REAL(theta)[0]);
      REAL(ret)[21] = REAL(theta)[1];
      REAL(ret)[32] = REAL(theta)[3];
      REAL(ret)[43] = REAL(theta)[6];
      REAL(ret)[54] = REAL(theta)[10];
      REAL(ret)[65] = REAL(theta)[15];
      REAL(ret)[76] = REAL(theta)[21];
      REAL(ret)[87] = REAL(theta)[28];
      REAL(ret)[98] = REAL(theta)[36];
      REAL(ret)[109] = REAL(theta)[45];
      REAL(ret)[110] = exp(REAL(theta)[0]);
      REAL(ret)[111] = REAL(theta)[1];
      REAL(ret)[112] = REAL(theta)[3];
      REAL(ret)[113] = REAL(theta)[6];
      REAL(ret)[114] = REAL(theta)[10];
      REAL(ret)[115] = REAL(theta)[15];
      REAL(ret)[116] = REAL(theta)[21];
      REAL(ret)[117] = REAL(theta)[28];
      REAL(ret)[118] = REAL(theta)[36];
      REAL(ret)[119] = REAL(theta)[45];
      REAL(ret)[120] = 2 * REAL(theta)[55];
    }
    else if (theta_n == 57){
      REAL(ret)[21] = exp(REAL(theta)[2]);
      REAL(ret)[32] = REAL(theta)[4];
      REAL(ret)[43] = REAL(theta)[7];
      REAL(ret)[54] = REAL(theta)[11];
      REAL(ret)[65] = REAL(theta)[16];
      REAL(ret)[76] = REAL(theta)[22];
      REAL(ret)[87] = REAL(theta)[29];
      REAL(ret)[98] = REAL(theta)[37];
      REAL(ret)[109] = REAL(theta)[46];
      REAL(ret)[111] = exp(REAL(theta)[2]);
      REAL(ret)[112] = REAL(theta)[4];
      REAL(ret)[113] = REAL(theta)[7];
      REAL(ret)[114] = REAL(theta)[11];
      REAL(ret)[115] = REAL(theta)[16];
      REAL(ret)[116] = REAL(theta)[22];
      REAL(ret)[117] = REAL(theta)[29];
      REAL(ret)[118] = REAL(theta)[37];
      REAL(ret)[119] = REAL(theta)[46];
      REAL(ret)[120] = 2 * REAL(theta)[56];
    }
    else if (theta_n == 58){
      REAL(ret)[32] = exp(REAL(theta)[5]);
      REAL(ret)[43] = REAL(theta)[8];
      REAL(ret)[54] = REAL(theta)[12];
      REAL(ret)[65] = REAL(theta)[17];
      REAL(ret)[76] = REAL(theta)[23];
      REAL(ret)[87] = REAL(theta)[30];
      REAL(ret)[98] = REAL(theta)[38];
      REAL(ret)[109] = REAL(theta)[47];
      REAL(ret)[112] = exp(REAL(theta)[5]);
      REAL(ret)[113] = REAL(theta)[8];
      REAL(ret)[114] = REAL(theta)[12];
      REAL(ret)[115] = REAL(theta)[17];
      REAL(ret)[116] = REAL(theta)[23];
      REAL(ret)[117] = REAL(theta)[30];
      REAL(ret)[118] = REAL(theta)[38];
      REAL(ret)[119] = REAL(theta)[47];
      REAL(ret)[120] = 2 * REAL(theta)[57];
    }
    else if (theta_n == 59){
      REAL(ret)[43] = exp(REAL(theta)[9]);
      REAL(ret)[54] = REAL(theta)[13];
      REAL(ret)[65] = REAL(theta)[18];
      REAL(ret)[76] = REAL(theta)[24];
      REAL(ret)[87] = REAL(theta)[31];
      REAL(ret)[98] = REAL(theta)[39];
      REAL(ret)[109] = REAL(theta)[48];
      REAL(ret)[113] = exp(REAL(theta)[9]);
      REAL(ret)[114] = REAL(theta)[13];
      REAL(ret)[115] = REAL(theta)[18];
      REAL(ret)[116] = REAL(theta)[24];
      REAL(ret)[117] = REAL(theta)[31];
      REAL(ret)[118] = REAL(theta)[39];
      REAL(ret)[119] = REAL(theta)[48];
      REAL(ret)[120] = 2 * REAL(theta)[58];
    }
    else if (theta_n == 60){
      REAL(ret)[54] = exp(REAL(theta)[14]);
      REAL(ret)[65] = REAL(theta)[19];
      REAL(ret)[76] = REAL(theta)[25];
      REAL(ret)[87] = REAL(theta)[32];
      REAL(ret)[98] = REAL(theta)[40];
      REAL(ret)[109] = REAL(theta)[49];
      REAL(ret)[114] = exp(REAL(theta)[14]);
      REAL(ret)[115] = REAL(theta)[19];
      REAL(ret)[116] = REAL(theta)[25];
      REAL(ret)[117] = REAL(theta)[32];
      REAL(ret)[118] = REAL(theta)[40];
      REAL(ret)[119] = REAL(theta)[49];
      REAL(ret)[120] = 2 * REAL(theta)[59];
    }
    else if (theta_n == 61){
      REAL(ret)[65] = exp(REAL(theta)[20]);
      REAL(ret)[76] = REAL(theta)[26];
      REAL(ret)[87] = REAL(theta)[33];
      REAL(ret)[98] = REAL(theta)[41];
      REAL(ret)[109] = REAL(theta)[50];
      REAL(ret)[115] = exp(REAL(theta)[20]);
      REAL(ret)[116] = REAL(theta)[26];
      REAL(ret)[117] = REAL(theta)[33];
      REAL(ret)[118] = REAL(theta)[41];
      REAL(ret)[119] = REAL(theta)[50];
      REAL(ret)[120] = 2 * REAL(theta)[60];
    }
    else if (theta_n == 62){
      REAL(ret)[76] = exp(REAL(theta)[27]);
      REAL(ret)[87] = REAL(theta)[34];
      REAL(ret)[98] = REAL(theta)[42];
      REAL(ret)[109] = REAL(theta)[51];
      REAL(ret)[116] = exp(REAL(theta)[27]);
      REAL(ret)[117] = REAL(theta)[34];
      REAL(ret)[118] = REAL(theta)[42];
      REAL(ret)[119] = REAL(theta)[51];
      REAL(ret)[120] = 2 * REAL(theta)[61];
    }
    else if (theta_n == 63){
      REAL(ret)[87] = exp(REAL(theta)[35]);
      REAL(ret)[98] = REAL(theta)[43];
      REAL(ret)[109] = REAL(theta)[52];
      REAL(ret)[117] = exp(REAL(theta)[35]);
      REAL(ret)[118] = REAL(theta)[43];
      REAL(ret)[119] = REAL(theta)[52];
      REAL(ret)[120] = 2 * REAL(theta)[62];
    }
    else if (theta_n == 64){
      REAL(ret)[98] = exp(REAL(theta)[44]);
      REAL(ret)[109] = REAL(theta)[53];
      REAL(ret)[118] = exp(REAL(theta)[44]);
      REAL(ret)[119] = REAL(theta)[53];
      REAL(ret)[120] = 2 * REAL(theta)[63];
    }
    else if (theta_n == 65){
      REAL(ret)[109] = exp(REAL(theta)[54]);
      REAL(ret)[119] = exp(REAL(theta)[54]);
      REAL(ret)[120] = 2 * REAL(theta)[64];
    }
    else if (theta_n == 66){
      REAL(ret)[120] = 2 * exp(2 * REAL(theta)[65]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 11));for(int i = 0; i < 11; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 2 * exp(2 * REAL(theta)[2]);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 2 * exp(2 * REAL(theta)[5]);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 2 * exp(2 * REAL(theta)[9]);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 2 * exp(2 * REAL(theta)[14]);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 2 * exp(2 * REAL(theta)[20]);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 2 * exp(2 * REAL(theta)[27]);
    }
    else if (theta_n == -38){
      REAL(ret)[7] = 2 * exp(2 * REAL(theta)[35]);
    }
    else if (theta_n == -47){
      REAL(ret)[8] = 2 * exp(2 * REAL(theta)[44]);
    }
    else if (theta_n == -57){
      REAL(ret)[9] = 2 * exp(2 * REAL(theta)[54]);
    }
    else if (theta_n == -68){
      REAL(ret)[10] = 2 * exp(2 * REAL(theta)[65]);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 12){
  int theta_n = INTEGER(tn)[0];
  if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 12;
    UNPROTECT(1);
    return ret;
  }
  else if (theta_n < -80 || theta_n > 78){
    error("d(Omega^-1) Derivative outside bounds.");
  }
  else if (length(theta) != 78){
    error("Requires vector with 78 arguments.");
  }
  if (theta_n >= -1){
    SEXP ret = PROTECT(allocMatrix(REALSXP, 12, 12));for (int i = 0; i < 144; i++){REAL(ret)[i]=0;}
    if (theta_n == 0){
      REAL(ret)[0] = exp(REAL(theta)[0]);
      REAL(ret)[12] = REAL(theta)[1];
      REAL(ret)[13] = exp(REAL(theta)[2]);
      REAL(ret)[24] = REAL(theta)[3];
      REAL(ret)[25] = REAL(theta)[4];
      REAL(ret)[26] = exp(REAL(theta)[5]);
      REAL(ret)[36] = REAL(theta)[6];
      REAL(ret)[37] = REAL(theta)[7];
      REAL(ret)[38] = REAL(theta)[8];
      REAL(ret)[39] = exp(REAL(theta)[9]);
      REAL(ret)[48] = REAL(theta)[10];
      REAL(ret)[49] = REAL(theta)[11];
      REAL(ret)[50] = REAL(theta)[12];
      REAL(ret)[51] = REAL(theta)[13];
      REAL(ret)[52] = exp(REAL(theta)[14]);
      REAL(ret)[60] = REAL(theta)[15];
      REAL(ret)[61] = REAL(theta)[16];
      REAL(ret)[62] = REAL(theta)[17];
      REAL(ret)[63] = REAL(theta)[18];
      REAL(ret)[64] = REAL(theta)[19];
      REAL(ret)[65] = exp(REAL(theta)[20]);
      REAL(ret)[72] = REAL(theta)[21];
      REAL(ret)[73] = REAL(theta)[22];
      REAL(ret)[74] = REAL(theta)[23];
      REAL(ret)[75] = REAL(theta)[24];
      REAL(ret)[76] = REAL(theta)[25];
      REAL(ret)[77] = REAL(theta)[26];
      REAL(ret)[78] = exp(REAL(theta)[27]);
      REAL(ret)[84] = REAL(theta)[28];
      REAL(ret)[85] = REAL(theta)[29];
      REAL(ret)[86] = REAL(theta)[30];
      REAL(ret)[87] = REAL(theta)[31];
      REAL(ret)[88] = REAL(theta)[32];
      REAL(ret)[89] = REAL(theta)[33];
      REAL(ret)[90] = REAL(theta)[34];
      REAL(ret)[91] = exp(REAL(theta)[35]);
      REAL(ret)[96] = REAL(theta)[36];
      REAL(ret)[97] = REAL(theta)[37];
      REAL(ret)[98] = REAL(theta)[38];
      REAL(ret)[99] = REAL(theta)[39];
      REAL(ret)[100] = REAL(theta)[40];
      REAL(ret)[101] = REAL(theta)[41];
      REAL(ret)[102] = REAL(theta)[42];
      REAL(ret)[103] = REAL(theta)[43];
      REAL(ret)[104] = exp(REAL(theta)[44]);
      REAL(ret)[108] = REAL(theta)[45];
      REAL(ret)[109] = REAL(theta)[46];
      REAL(ret)[110] = REAL(theta)[47];
      REAL(ret)[111] = REAL(theta)[48];
      REAL(ret)[112] = REAL(theta)[49];
      REAL(ret)[113] = REAL(theta)[50];
      REAL(ret)[114] = REAL(theta)[51];
      REAL(ret)[115] = REAL(theta)[52];
      REAL(ret)[116] = REAL(theta)[53];
      REAL(ret)[117] = exp(REAL(theta)[54]);
      REAL(ret)[120] = REAL(theta)[55];
      REAL(ret)[121] = REAL(theta)[56];
      REAL(ret)[122] = REAL(theta)[57];
      REAL(ret)[123] = REAL(theta)[58];
      REAL(ret)[124] = REAL(theta)[59];
      REAL(ret)[125] = REAL(theta)[60];
      REAL(ret)[126] = REAL(theta)[61];
      REAL(ret)[127] = REAL(theta)[62];
      REAL(ret)[128] = REAL(theta)[63];
      REAL(ret)[129] = REAL(theta)[64];
      REAL(ret)[130] = exp(REAL(theta)[65]);
      REAL(ret)[132] = REAL(theta)[66];
      REAL(ret)[133] = REAL(theta)[67];
      REAL(ret)[134] = REAL(theta)[68];
      REAL(ret)[135] = REAL(theta)[69];
      REAL(ret)[136] = REAL(theta)[70];
      REAL(ret)[137] = REAL(theta)[71];
      REAL(ret)[138] = REAL(theta)[72];
      REAL(ret)[139] = REAL(theta)[73];
      REAL(ret)[140] = REAL(theta)[74];
      REAL(ret)[141] = REAL(theta)[75];
      REAL(ret)[142] = REAL(theta)[76];
      REAL(ret)[143] = exp(REAL(theta)[77]);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[45] * exp(REAL(theta)[0]);
      REAL(ret)[10] = REAL(theta)[55] * exp(REAL(theta)[0]);
      REAL(ret)[11] = REAL(theta)[66] * exp(REAL(theta)[0]);
      REAL(ret)[12] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[13] = R_pow_di(REAL(theta)[1], 2) + exp(2 * REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[18] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[19] = REAL(theta)[1] * REAL(theta)[28] + REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[20] = REAL(theta)[1] * REAL(theta)[36] + REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[21] = REAL(theta)[1] * REAL(theta)[45] + REAL(theta)[46] * exp(REAL(theta)[2]);
      REAL(ret)[22] = REAL(theta)[1] * REAL(theta)[55] + REAL(theta)[56] * exp(REAL(theta)[2]);
      REAL(ret)[23] = REAL(theta)[1] * REAL(theta)[66] + REAL(theta)[67] * exp(REAL(theta)[2]);
      REAL(ret)[24] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[25] = REAL(theta)[1] * REAL(theta)[3] + REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[26] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + exp(2 * REAL(theta)[5]);
      REAL(ret)[27] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[28] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[29] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[30] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[31] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[32] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[33] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * exp(REAL(theta)[5]);
      REAL(ret)[34] = REAL(theta)[3] * REAL(theta)[55] + REAL(theta)[4] * REAL(theta)[56] + REAL(theta)[57] * exp(REAL(theta)[5]);
      REAL(ret)[35] = REAL(theta)[3] * REAL(theta)[66] + REAL(theta)[4] * REAL(theta)[67] + REAL(theta)[68] * exp(REAL(theta)[5]);
      REAL(ret)[36] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[37] = REAL(theta)[1] * REAL(theta)[6] + REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[38] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[39] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + exp(2 * REAL(theta)[9]);
      REAL(ret)[40] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[41] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[42] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[43] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[44] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[45] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * exp(REAL(theta)[9]);
      REAL(ret)[46] = REAL(theta)[55] * REAL(theta)[6] + REAL(theta)[56] * REAL(theta)[7] + REAL(theta)[57] * REAL(theta)[8] + REAL(theta)[58] * exp(REAL(theta)[9]);
      REAL(ret)[47] = REAL(theta)[6] * REAL(theta)[66] + REAL(theta)[67] * REAL(theta)[7] + REAL(theta)[68] * REAL(theta)[8] + REAL(theta)[69] * exp(REAL(theta)[9]);
      REAL(ret)[48] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[49] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[50] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[51] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[52] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + exp(2 * REAL(theta)[14]);
      REAL(ret)[53] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[54] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[55] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[56] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[57] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + REAL(theta)[49] * exp(REAL(theta)[14]);
      REAL(ret)[58] = REAL(theta)[10] * REAL(theta)[55] + REAL(theta)[11] * REAL(theta)[56] + REAL(theta)[12] * REAL(theta)[57] + REAL(theta)[13] * REAL(theta)[58] + REAL(theta)[59] * exp(REAL(theta)[14]);
      REAL(ret)[59] = REAL(theta)[10] * REAL(theta)[66] + REAL(theta)[11] * REAL(theta)[67] + REAL(theta)[12] * REAL(theta)[68] + REAL(theta)[13] * REAL(theta)[69] + REAL(theta)[70] * exp(REAL(theta)[14]);
      REAL(ret)[60] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[61] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[62] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[63] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[64] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[65] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + exp(2 * REAL(theta)[20]);
      REAL(ret)[66] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[67] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[68] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[69] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + REAL(theta)[50] * exp(REAL(theta)[20]);
      REAL(ret)[70] = REAL(theta)[15] * REAL(theta)[55] + REAL(theta)[16] * REAL(theta)[56] + REAL(theta)[17] * REAL(theta)[57] + REAL(theta)[18] * REAL(theta)[58] + REAL(theta)[19] * REAL(theta)[59] + REAL(theta)[60] * exp(REAL(theta)[20]);
      REAL(ret)[71] = REAL(theta)[15] * REAL(theta)[66] + REAL(theta)[16] * REAL(theta)[67] + REAL(theta)[17] * REAL(theta)[68] + REAL(theta)[18] * REAL(theta)[69] + REAL(theta)[19] * REAL(theta)[70] + REAL(theta)[71] * exp(REAL(theta)[20]);
      REAL(ret)[72] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[73] = REAL(theta)[1] * REAL(theta)[21] + REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[74] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[75] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[76] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[77] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[78] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + exp(2 * REAL(theta)[27]);
      REAL(ret)[79] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[80] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[81] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + REAL(theta)[51] * exp(REAL(theta)[27]);
      REAL(ret)[82] = REAL(theta)[21] * REAL(theta)[55] + REAL(theta)[22] * REAL(theta)[56] + REAL(theta)[23] * REAL(theta)[57] + REAL(theta)[24] * REAL(theta)[58] + REAL(theta)[25] * REAL(theta)[59] + REAL(theta)[26] * REAL(theta)[60] + REAL(theta)[61] * exp(REAL(theta)[27]);
      REAL(ret)[83] = REAL(theta)[21] * REAL(theta)[66] + REAL(theta)[22] * REAL(theta)[67] + REAL(theta)[23] * REAL(theta)[68] + REAL(theta)[24] * REAL(theta)[69] + REAL(theta)[25] * REAL(theta)[70] + REAL(theta)[26] * REAL(theta)[71] + REAL(theta)[72] * exp(REAL(theta)[27]);
      REAL(ret)[84] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[85] = REAL(theta)[1] * REAL(theta)[28] + REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[86] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[87] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[88] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[89] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[90] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[91] = R_pow_di(REAL(theta)[28], 2) + R_pow_di(REAL(theta)[29], 2) + R_pow_di(REAL(theta)[30], 2) + R_pow_di(REAL(theta)[31], 2) + R_pow_di(REAL(theta)[32], 2) + R_pow_di(REAL(theta)[33], 2) + R_pow_di(REAL(theta)[34], 2) + exp(2 * REAL(theta)[35]);
      REAL(ret)[92] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[93] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + REAL(theta)[52] * exp(REAL(theta)[35]);
      REAL(ret)[94] = REAL(theta)[28] * REAL(theta)[55] + REAL(theta)[29] * REAL(theta)[56] + REAL(theta)[30] * REAL(theta)[57] + REAL(theta)[31] * REAL(theta)[58] + REAL(theta)[32] * REAL(theta)[59] + REAL(theta)[33] * REAL(theta)[60] + REAL(theta)[34] * REAL(theta)[61] + REAL(theta)[62] * exp(REAL(theta)[35]);
      REAL(ret)[95] = REAL(theta)[28] * REAL(theta)[66] + REAL(theta)[29] * REAL(theta)[67] + REAL(theta)[30] * REAL(theta)[68] + REAL(theta)[31] * REAL(theta)[69] + REAL(theta)[32] * REAL(theta)[70] + REAL(theta)[33] * REAL(theta)[71] + REAL(theta)[34] * REAL(theta)[72] + REAL(theta)[73] * exp(REAL(theta)[35]);
      REAL(ret)[96] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[97] = REAL(theta)[1] * REAL(theta)[36] + REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[98] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[99] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[100] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[101] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[102] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[103] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[104] = R_pow_di(REAL(theta)[36], 2) + R_pow_di(REAL(theta)[37], 2) + R_pow_di(REAL(theta)[38], 2) + R_pow_di(REAL(theta)[39], 2) + R_pow_di(REAL(theta)[40], 2) + R_pow_di(REAL(theta)[41], 2) + R_pow_di(REAL(theta)[42], 2) + R_pow_di(REAL(theta)[43], 2) + exp(2 * REAL(theta)[44]);
      REAL(ret)[105] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + REAL(theta)[53] * exp(REAL(theta)[44]);
      REAL(ret)[106] = REAL(theta)[36] * REAL(theta)[55] + REAL(theta)[37] * REAL(theta)[56] + REAL(theta)[38] * REAL(theta)[57] + REAL(theta)[39] * REAL(theta)[58] + REAL(theta)[40] * REAL(theta)[59] + REAL(theta)[41] * REAL(theta)[60] + REAL(theta)[42] * REAL(theta)[61] + REAL(theta)[43] * REAL(theta)[62] + REAL(theta)[63] * exp(REAL(theta)[44]);
      REAL(ret)[107] = REAL(theta)[36] * REAL(theta)[66] + REAL(theta)[37] * REAL(theta)[67] + REAL(theta)[38] * REAL(theta)[68] + REAL(theta)[39] * REAL(theta)[69] + REAL(theta)[40] * REAL(theta)[70] + REAL(theta)[41] * REAL(theta)[71] + REAL(theta)[42] * REAL(theta)[72] + REAL(theta)[43] * REAL(theta)[73] + REAL(theta)[74] * exp(REAL(theta)[44]);
      REAL(ret)[108] = REAL(theta)[45] * exp(REAL(theta)[0]);
      REAL(ret)[109] = REAL(theta)[1] * REAL(theta)[45] + REAL(theta)[46] * exp(REAL(theta)[2]);
      REAL(ret)[110] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * exp(REAL(theta)[5]);
      REAL(ret)[111] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * exp(REAL(theta)[9]);
      REAL(ret)[112] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + REAL(theta)[49] * exp(REAL(theta)[14]);
      REAL(ret)[113] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + REAL(theta)[50] * exp(REAL(theta)[20]);
      REAL(ret)[114] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + REAL(theta)[51] * exp(REAL(theta)[27]);
      REAL(ret)[115] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + REAL(theta)[52] * exp(REAL(theta)[35]);
      REAL(ret)[116] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + REAL(theta)[53] * exp(REAL(theta)[44]);
      REAL(ret)[117] = R_pow_di(REAL(theta)[45], 2) + R_pow_di(REAL(theta)[46], 2) + R_pow_di(REAL(theta)[47], 2) + R_pow_di(REAL(theta)[48], 2) + R_pow_di(REAL(theta)[49], 2) + R_pow_di(REAL(theta)[50], 2) + R_pow_di(REAL(theta)[51], 2) + R_pow_di(REAL(theta)[52], 2) + R_pow_di(REAL(theta)[53], 2) + exp(2 * REAL(theta)[54]);
      REAL(ret)[118] = REAL(theta)[45] * REAL(theta)[55] + REAL(theta)[46] * REAL(theta)[56] + REAL(theta)[47] * REAL(theta)[57] + REAL(theta)[48] * REAL(theta)[58] + REAL(theta)[49] * REAL(theta)[59] + REAL(theta)[50] * REAL(theta)[60] + REAL(theta)[51] * REAL(theta)[61] + REAL(theta)[52] * REAL(theta)[62] + REAL(theta)[53] * REAL(theta)[63] + REAL(theta)[64] * exp(REAL(theta)[54]);
      REAL(ret)[119] = REAL(theta)[45] * REAL(theta)[66] + REAL(theta)[46] * REAL(theta)[67] + REAL(theta)[47] * REAL(theta)[68] + REAL(theta)[48] * REAL(theta)[69] + REAL(theta)[49] * REAL(theta)[70] + REAL(theta)[50] * REAL(theta)[71] + REAL(theta)[51] * REAL(theta)[72] + REAL(theta)[52] * REAL(theta)[73] + REAL(theta)[53] * REAL(theta)[74] + REAL(theta)[75] * exp(REAL(theta)[54]);
      REAL(ret)[120] = REAL(theta)[55] * exp(REAL(theta)[0]);
      REAL(ret)[121] = REAL(theta)[1] * REAL(theta)[55] + REAL(theta)[56] * exp(REAL(theta)[2]);
      REAL(ret)[122] = REAL(theta)[3] * REAL(theta)[55] + REAL(theta)[4] * REAL(theta)[56] + REAL(theta)[57] * exp(REAL(theta)[5]);
      REAL(ret)[123] = REAL(theta)[55] * REAL(theta)[6] + REAL(theta)[56] * REAL(theta)[7] + REAL(theta)[57] * REAL(theta)[8] + REAL(theta)[58] * exp(REAL(theta)[9]);
      REAL(ret)[124] = REAL(theta)[10] * REAL(theta)[55] + REAL(theta)[11] * REAL(theta)[56] + REAL(theta)[12] * REAL(theta)[57] + REAL(theta)[13] * REAL(theta)[58] + REAL(theta)[59] * exp(REAL(theta)[14]);
      REAL(ret)[125] = REAL(theta)[15] * REAL(theta)[55] + REAL(theta)[16] * REAL(theta)[56] + REAL(theta)[17] * REAL(theta)[57] + REAL(theta)[18] * REAL(theta)[58] + REAL(theta)[19] * REAL(theta)[59] + REAL(theta)[60] * exp(REAL(theta)[20]);
      REAL(ret)[126] = REAL(theta)[21] * REAL(theta)[55] + REAL(theta)[22] * REAL(theta)[56] + REAL(theta)[23] * REAL(theta)[57] + REAL(theta)[24] * REAL(theta)[58] + REAL(theta)[25] * REAL(theta)[59] + REAL(theta)[26] * REAL(theta)[60] + REAL(theta)[61] * exp(REAL(theta)[27]);
      REAL(ret)[127] = REAL(theta)[28] * REAL(theta)[55] + REAL(theta)[29] * REAL(theta)[56] + REAL(theta)[30] * REAL(theta)[57] + REAL(theta)[31] * REAL(theta)[58] + REAL(theta)[32] * REAL(theta)[59] + REAL(theta)[33] * REAL(theta)[60] + REAL(theta)[34] * REAL(theta)[61] + REAL(theta)[62] * exp(REAL(theta)[35]);
      REAL(ret)[128] = REAL(theta)[36] * REAL(theta)[55] + REAL(theta)[37] * REAL(theta)[56] + REAL(theta)[38] * REAL(theta)[57] + REAL(theta)[39] * REAL(theta)[58] + REAL(theta)[40] * REAL(theta)[59] + REAL(theta)[41] * REAL(theta)[60] + REAL(theta)[42] * REAL(theta)[61] + REAL(theta)[43] * REAL(theta)[62] + REAL(theta)[63] * exp(REAL(theta)[44]);
      REAL(ret)[129] = REAL(theta)[45] * REAL(theta)[55] + REAL(theta)[46] * REAL(theta)[56] + REAL(theta)[47] * REAL(theta)[57] + REAL(theta)[48] * REAL(theta)[58] + REAL(theta)[49] * REAL(theta)[59] + REAL(theta)[50] * REAL(theta)[60] + REAL(theta)[51] * REAL(theta)[61] + REAL(theta)[52] * REAL(theta)[62] + REAL(theta)[53] * REAL(theta)[63] + REAL(theta)[64] * exp(REAL(theta)[54]);
      REAL(ret)[130] = R_pow_di(REAL(theta)[55], 2) + R_pow_di(REAL(theta)[56], 2) + R_pow_di(REAL(theta)[57], 2) + R_pow_di(REAL(theta)[58], 2) + R_pow_di(REAL(theta)[59], 2) + R_pow_di(REAL(theta)[60], 2) + R_pow_di(REAL(theta)[61], 2) + R_pow_di(REAL(theta)[62], 2) + R_pow_di(REAL(theta)[63], 2) + R_pow_di(REAL(theta)[64], 2) + exp(2 * REAL(theta)[65]);
      REAL(ret)[131] = REAL(theta)[55] * REAL(theta)[66] + REAL(theta)[56] * REAL(theta)[67] + REAL(theta)[57] * REAL(theta)[68] + REAL(theta)[58] * REAL(theta)[69] + REAL(theta)[59] * REAL(theta)[70] + REAL(theta)[60] * REAL(theta)[71] + REAL(theta)[61] * REAL(theta)[72] + REAL(theta)[62] * REAL(theta)[73] + REAL(theta)[63] * REAL(theta)[74] + REAL(theta)[64] * REAL(theta)[75] + REAL(theta)[76] * exp(REAL(theta)[65]);
      REAL(ret)[132] = REAL(theta)[66] * exp(REAL(theta)[0]);
      REAL(ret)[133] = REAL(theta)[1] * REAL(theta)[66] + REAL(theta)[67] * exp(REAL(theta)[2]);
      REAL(ret)[134] = REAL(theta)[3] * REAL(theta)[66] + REAL(theta)[4] * REAL(theta)[67] + REAL(theta)[68] * exp(REAL(theta)[5]);
      REAL(ret)[135] = REAL(theta)[6] * REAL(theta)[66] + REAL(theta)[67] * REAL(theta)[7] + REAL(theta)[68] * REAL(theta)[8] + REAL(theta)[69] * exp(REAL(theta)[9]);
      REAL(ret)[136] = REAL(theta)[10] * REAL(theta)[66] + REAL(theta)[11] * REAL(theta)[67] + REAL(theta)[12] * REAL(theta)[68] + REAL(theta)[13] * REAL(theta)[69] + REAL(theta)[70] * exp(REAL(theta)[14]);
      REAL(ret)[137] = REAL(theta)[15] * REAL(theta)[66] + REAL(theta)[16] * REAL(theta)[67] + REAL(theta)[17] * REAL(theta)[68] + REAL(theta)[18] * REAL(theta)[69] + REAL(theta)[19] * REAL(theta)[70] + REAL(theta)[71] * exp(REAL(theta)[20]);
      REAL(ret)[138] = REAL(theta)[21] * REAL(theta)[66] + REAL(theta)[22] * REAL(theta)[67] + REAL(theta)[23] * REAL(theta)[68] + REAL(theta)[24] * REAL(theta)[69] + REAL(theta)[25] * REAL(theta)[70] + REAL(theta)[26] * REAL(theta)[71] + REAL(theta)[72] * exp(REAL(theta)[27]);
      REAL(ret)[139] = REAL(theta)[28] * REAL(theta)[66] + REAL(theta)[29] * REAL(theta)[67] + REAL(theta)[30] * REAL(theta)[68] + REAL(theta)[31] * REAL(theta)[69] + REAL(theta)[32] * REAL(theta)[70] + REAL(theta)[33] * REAL(theta)[71] + REAL(theta)[34] * REAL(theta)[72] + REAL(theta)[73] * exp(REAL(theta)[35]);
      REAL(ret)[140] = REAL(theta)[36] * REAL(theta)[66] + REAL(theta)[37] * REAL(theta)[67] + REAL(theta)[38] * REAL(theta)[68] + REAL(theta)[39] * REAL(theta)[69] + REAL(theta)[40] * REAL(theta)[70] + REAL(theta)[41] * REAL(theta)[71] + REAL(theta)[42] * REAL(theta)[72] + REAL(theta)[43] * REAL(theta)[73] + REAL(theta)[74] * exp(REAL(theta)[44]);
      REAL(ret)[141] = REAL(theta)[45] * REAL(theta)[66] + REAL(theta)[46] * REAL(theta)[67] + REAL(theta)[47] * REAL(theta)[68] + REAL(theta)[48] * REAL(theta)[69] + REAL(theta)[49] * REAL(theta)[70] + REAL(theta)[50] * REAL(theta)[71] + REAL(theta)[51] * REAL(theta)[72] + REAL(theta)[52] * REAL(theta)[73] + REAL(theta)[53] * REAL(theta)[74] + REAL(theta)[75] * exp(REAL(theta)[54]);
      REAL(ret)[142] = REAL(theta)[55] * REAL(theta)[66] + REAL(theta)[56] * REAL(theta)[67] + REAL(theta)[57] * REAL(theta)[68] + REAL(theta)[58] * REAL(theta)[69] + REAL(theta)[59] * REAL(theta)[70] + REAL(theta)[60] * REAL(theta)[71] + REAL(theta)[61] * REAL(theta)[72] + REAL(theta)[62] * REAL(theta)[73] + REAL(theta)[63] * REAL(theta)[74] + REAL(theta)[64] * REAL(theta)[75] + REAL(theta)[76] * exp(REAL(theta)[65]);
      REAL(ret)[143] = R_pow_di(REAL(theta)[66], 2) + R_pow_di(REAL(theta)[67], 2) + R_pow_di(REAL(theta)[68], 2) + R_pow_di(REAL(theta)[69], 2) + R_pow_di(REAL(theta)[70], 2) + R_pow_di(REAL(theta)[71], 2) + R_pow_di(REAL(theta)[72], 2) + R_pow_di(REAL(theta)[73], 2) + R_pow_di(REAL(theta)[74], 2) + R_pow_di(REAL(theta)[75], 2) + R_pow_di(REAL(theta)[76], 2) + exp(2 * REAL(theta)[77]);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
      REAL(ret)[1] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[2] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[3] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[4] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[5] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[6] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[7] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[8] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[9] = REAL(theta)[45] * exp(REAL(theta)[0]);
      REAL(ret)[10] = REAL(theta)[55] * exp(REAL(theta)[0]);
      REAL(ret)[11] = REAL(theta)[66] * exp(REAL(theta)[0]);
      REAL(ret)[12] = REAL(theta)[1] * exp(REAL(theta)[0]);
      REAL(ret)[24] = REAL(theta)[3] * exp(REAL(theta)[0]);
      REAL(ret)[36] = REAL(theta)[6] * exp(REAL(theta)[0]);
      REAL(ret)[48] = REAL(theta)[10] * exp(REAL(theta)[0]);
      REAL(ret)[60] = REAL(theta)[15] * exp(REAL(theta)[0]);
      REAL(ret)[72] = REAL(theta)[21] * exp(REAL(theta)[0]);
      REAL(ret)[84] = REAL(theta)[28] * exp(REAL(theta)[0]);
      REAL(ret)[96] = REAL(theta)[36] * exp(REAL(theta)[0]);
      REAL(ret)[108] = REAL(theta)[45] * exp(REAL(theta)[0]);
      REAL(ret)[120] = REAL(theta)[55] * exp(REAL(theta)[0]);
      REAL(ret)[132] = REAL(theta)[66] * exp(REAL(theta)[0]);
    }
    else if (theta_n == 2){
      REAL(ret)[1] = exp(REAL(theta)[0]);
      REAL(ret)[12] = exp(REAL(theta)[0]);
      REAL(ret)[13] = 2 * REAL(theta)[1];
      REAL(ret)[14] = REAL(theta)[3];
      REAL(ret)[15] = REAL(theta)[6];
      REAL(ret)[16] = REAL(theta)[10];
      REAL(ret)[17] = REAL(theta)[15];
      REAL(ret)[18] = REAL(theta)[21];
      REAL(ret)[19] = REAL(theta)[28];
      REAL(ret)[20] = REAL(theta)[36];
      REAL(ret)[21] = REAL(theta)[45];
      REAL(ret)[22] = REAL(theta)[55];
      REAL(ret)[23] = REAL(theta)[66];
      REAL(ret)[25] = REAL(theta)[3];
      REAL(ret)[37] = REAL(theta)[6];
      REAL(ret)[49] = REAL(theta)[10];
      REAL(ret)[61] = REAL(theta)[15];
      REAL(ret)[73] = REAL(theta)[21];
      REAL(ret)[85] = REAL(theta)[28];
      REAL(ret)[97] = REAL(theta)[36];
      REAL(ret)[109] = REAL(theta)[45];
      REAL(ret)[121] = REAL(theta)[55];
      REAL(ret)[133] = REAL(theta)[66];
    }
    else if (theta_n == 3){
      REAL(ret)[13] = 2 * exp(2 * REAL(theta)[2]);
      REAL(ret)[14] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[15] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[16] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[17] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[18] = REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[19] = REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[20] = REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[21] = REAL(theta)[46] * exp(REAL(theta)[2]);
      REAL(ret)[22] = REAL(theta)[56] * exp(REAL(theta)[2]);
      REAL(ret)[23] = REAL(theta)[67] * exp(REAL(theta)[2]);
      REAL(ret)[25] = REAL(theta)[4] * exp(REAL(theta)[2]);
      REAL(ret)[37] = REAL(theta)[7] * exp(REAL(theta)[2]);
      REAL(ret)[49] = REAL(theta)[11] * exp(REAL(theta)[2]);
      REAL(ret)[61] = REAL(theta)[16] * exp(REAL(theta)[2]);
      REAL(ret)[73] = REAL(theta)[22] * exp(REAL(theta)[2]);
      REAL(ret)[85] = REAL(theta)[29] * exp(REAL(theta)[2]);
      REAL(ret)[97] = REAL(theta)[37] * exp(REAL(theta)[2]);
      REAL(ret)[109] = REAL(theta)[46] * exp(REAL(theta)[2]);
      REAL(ret)[121] = REAL(theta)[56] * exp(REAL(theta)[2]);
      REAL(ret)[133] = REAL(theta)[67] * exp(REAL(theta)[2]);
    }
    else if (theta_n == 4){
      REAL(ret)[2] = exp(REAL(theta)[0]);
      REAL(ret)[14] = REAL(theta)[1];
      REAL(ret)[24] = exp(REAL(theta)[0]);
      REAL(ret)[25] = REAL(theta)[1];
      REAL(ret)[26] = 2 * REAL(theta)[3];
      REAL(ret)[27] = REAL(theta)[6];
      REAL(ret)[28] = REAL(theta)[10];
      REAL(ret)[29] = REAL(theta)[15];
      REAL(ret)[30] = REAL(theta)[21];
      REAL(ret)[31] = REAL(theta)[28];
      REAL(ret)[32] = REAL(theta)[36];
      REAL(ret)[33] = REAL(theta)[45];
      REAL(ret)[34] = REAL(theta)[55];
      REAL(ret)[35] = REAL(theta)[66];
      REAL(ret)[38] = REAL(theta)[6];
      REAL(ret)[50] = REAL(theta)[10];
      REAL(ret)[62] = REAL(theta)[15];
      REAL(ret)[74] = REAL(theta)[21];
      REAL(ret)[86] = REAL(theta)[28];
      REAL(ret)[98] = REAL(theta)[36];
      REAL(ret)[110] = REAL(theta)[45];
      REAL(ret)[122] = REAL(theta)[55];
      REAL(ret)[134] = REAL(theta)[66];
    }
    else if (theta_n == 5){
      REAL(ret)[14] = exp(REAL(theta)[2]);
      REAL(ret)[25] = exp(REAL(theta)[2]);
      REAL(ret)[26] = 2 * REAL(theta)[4];
      REAL(ret)[27] = REAL(theta)[7];
      REAL(ret)[28] = REAL(theta)[11];
      REAL(ret)[29] = REAL(theta)[16];
      REAL(ret)[30] = REAL(theta)[22];
      REAL(ret)[31] = REAL(theta)[29];
      REAL(ret)[32] = REAL(theta)[37];
      REAL(ret)[33] = REAL(theta)[46];
      REAL(ret)[34] = REAL(theta)[56];
      REAL(ret)[35] = REAL(theta)[67];
      REAL(ret)[38] = REAL(theta)[7];
      REAL(ret)[50] = REAL(theta)[11];
      REAL(ret)[62] = REAL(theta)[16];
      REAL(ret)[74] = REAL(theta)[22];
      REAL(ret)[86] = REAL(theta)[29];
      REAL(ret)[98] = REAL(theta)[37];
      REAL(ret)[110] = REAL(theta)[46];
      REAL(ret)[122] = REAL(theta)[56];
      REAL(ret)[134] = REAL(theta)[67];
    }
    else if (theta_n == 6){
      REAL(ret)[26] = 2 * exp(2 * REAL(theta)[5]);
      REAL(ret)[27] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[28] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[29] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[30] = REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[31] = REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[32] = REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[33] = REAL(theta)[47] * exp(REAL(theta)[5]);
      REAL(ret)[34] = REAL(theta)[57] * exp(REAL(theta)[5]);
      REAL(ret)[35] = REAL(theta)[68] * exp(REAL(theta)[5]);
      REAL(ret)[38] = REAL(theta)[8] * exp(REAL(theta)[5]);
      REAL(ret)[50] = REAL(theta)[12] * exp(REAL(theta)[5]);
      REAL(ret)[62] = REAL(theta)[17] * exp(REAL(theta)[5]);
      REAL(ret)[74] = REAL(theta)[23] * exp(REAL(theta)[5]);
      REAL(ret)[86] = REAL(theta)[30] * exp(REAL(theta)[5]);
      REAL(ret)[98] = REAL(theta)[38] * exp(REAL(theta)[5]);
      REAL(ret)[110] = REAL(theta)[47] * exp(REAL(theta)[5]);
      REAL(ret)[122] = REAL(theta)[57] * exp(REAL(theta)[5]);
      REAL(ret)[134] = REAL(theta)[68] * exp(REAL(theta)[5]);
    }
    else if (theta_n == 7){
      REAL(ret)[3] = exp(REAL(theta)[0]);
      REAL(ret)[15] = REAL(theta)[1];
      REAL(ret)[27] = REAL(theta)[3];
      REAL(ret)[36] = exp(REAL(theta)[0]);
      REAL(ret)[37] = REAL(theta)[1];
      REAL(ret)[38] = REAL(theta)[3];
      REAL(ret)[39] = 2 * REAL(theta)[6];
      REAL(ret)[40] = REAL(theta)[10];
      REAL(ret)[41] = REAL(theta)[15];
      REAL(ret)[42] = REAL(theta)[21];
      REAL(ret)[43] = REAL(theta)[28];
      REAL(ret)[44] = REAL(theta)[36];
      REAL(ret)[45] = REAL(theta)[45];
      REAL(ret)[46] = REAL(theta)[55];
      REAL(ret)[47] = REAL(theta)[66];
      REAL(ret)[51] = REAL(theta)[10];
      REAL(ret)[63] = REAL(theta)[15];
      REAL(ret)[75] = REAL(theta)[21];
      REAL(ret)[87] = REAL(theta)[28];
      REAL(ret)[99] = REAL(theta)[36];
      REAL(ret)[111] = REAL(theta)[45];
      REAL(ret)[123] = REAL(theta)[55];
      REAL(ret)[135] = REAL(theta)[66];
    }
    else if (theta_n == 8){
      REAL(ret)[15] = exp(REAL(theta)[2]);
      REAL(ret)[27] = REAL(theta)[4];
      REAL(ret)[37] = exp(REAL(theta)[2]);
      REAL(ret)[38] = REAL(theta)[4];
      REAL(ret)[39] = 2 * REAL(theta)[7];
      REAL(ret)[40] = REAL(theta)[11];
      REAL(ret)[41] = REAL(theta)[16];
      REAL(ret)[42] = REAL(theta)[22];
      REAL(ret)[43] = REAL(theta)[29];
      REAL(ret)[44] = REAL(theta)[37];
      REAL(ret)[45] = REAL(theta)[46];
      REAL(ret)[46] = REAL(theta)[56];
      REAL(ret)[47] = REAL(theta)[67];
      REAL(ret)[51] = REAL(theta)[11];
      REAL(ret)[63] = REAL(theta)[16];
      REAL(ret)[75] = REAL(theta)[22];
      REAL(ret)[87] = REAL(theta)[29];
      REAL(ret)[99] = REAL(theta)[37];
      REAL(ret)[111] = REAL(theta)[46];
      REAL(ret)[123] = REAL(theta)[56];
      REAL(ret)[135] = REAL(theta)[67];
    }
    else if (theta_n == 9){
      REAL(ret)[27] = exp(REAL(theta)[5]);
      REAL(ret)[38] = exp(REAL(theta)[5]);
      REAL(ret)[39] = 2 * REAL(theta)[8];
      REAL(ret)[40] = REAL(theta)[12];
      REAL(ret)[41] = REAL(theta)[17];
      REAL(ret)[42] = REAL(theta)[23];
      REAL(ret)[43] = REAL(theta)[30];
      REAL(ret)[44] = REAL(theta)[38];
      REAL(ret)[45] = REAL(theta)[47];
      REAL(ret)[46] = REAL(theta)[57];
      REAL(ret)[47] = REAL(theta)[68];
      REAL(ret)[51] = REAL(theta)[12];
      REAL(ret)[63] = REAL(theta)[17];
      REAL(ret)[75] = REAL(theta)[23];
      REAL(ret)[87] = REAL(theta)[30];
      REAL(ret)[99] = REAL(theta)[38];
      REAL(ret)[111] = REAL(theta)[47];
      REAL(ret)[123] = REAL(theta)[57];
      REAL(ret)[135] = REAL(theta)[68];
    }
    else if (theta_n == 10){
      REAL(ret)[39] = 2 * exp(2 * REAL(theta)[9]);
      REAL(ret)[40] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[41] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[42] = REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[43] = REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[44] = REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[45] = REAL(theta)[48] * exp(REAL(theta)[9]);
      REAL(ret)[46] = REAL(theta)[58] * exp(REAL(theta)[9]);
      REAL(ret)[47] = REAL(theta)[69] * exp(REAL(theta)[9]);
      REAL(ret)[51] = REAL(theta)[13] * exp(REAL(theta)[9]);
      REAL(ret)[63] = REAL(theta)[18] * exp(REAL(theta)[9]);
      REAL(ret)[75] = REAL(theta)[24] * exp(REAL(theta)[9]);
      REAL(ret)[87] = REAL(theta)[31] * exp(REAL(theta)[9]);
      REAL(ret)[99] = REAL(theta)[39] * exp(REAL(theta)[9]);
      REAL(ret)[111] = REAL(theta)[48] * exp(REAL(theta)[9]);
      REAL(ret)[123] = REAL(theta)[58] * exp(REAL(theta)[9]);
      REAL(ret)[135] = REAL(theta)[69] * exp(REAL(theta)[9]);
    }
    else if (theta_n == 11){
      REAL(ret)[4] = exp(REAL(theta)[0]);
      REAL(ret)[16] = REAL(theta)[1];
      REAL(ret)[28] = REAL(theta)[3];
      REAL(ret)[40] = REAL(theta)[6];
      REAL(ret)[48] = exp(REAL(theta)[0]);
      REAL(ret)[49] = REAL(theta)[1];
      REAL(ret)[50] = REAL(theta)[3];
      REAL(ret)[51] = REAL(theta)[6];
      REAL(ret)[52] = 2 * REAL(theta)[10];
      REAL(ret)[53] = REAL(theta)[15];
      REAL(ret)[54] = REAL(theta)[21];
      REAL(ret)[55] = REAL(theta)[28];
      REAL(ret)[56] = REAL(theta)[36];
      REAL(ret)[57] = REAL(theta)[45];
      REAL(ret)[58] = REAL(theta)[55];
      REAL(ret)[59] = REAL(theta)[66];
      REAL(ret)[64] = REAL(theta)[15];
      REAL(ret)[76] = REAL(theta)[21];
      REAL(ret)[88] = REAL(theta)[28];
      REAL(ret)[100] = REAL(theta)[36];
      REAL(ret)[112] = REAL(theta)[45];
      REAL(ret)[124] = REAL(theta)[55];
      REAL(ret)[136] = REAL(theta)[66];
    }
    else if (theta_n == 12){
      REAL(ret)[16] = exp(REAL(theta)[2]);
      REAL(ret)[28] = REAL(theta)[4];
      REAL(ret)[40] = REAL(theta)[7];
      REAL(ret)[49] = exp(REAL(theta)[2]);
      REAL(ret)[50] = REAL(theta)[4];
      REAL(ret)[51] = REAL(theta)[7];
      REAL(ret)[52] = 2 * REAL(theta)[11];
      REAL(ret)[53] = REAL(theta)[16];
      REAL(ret)[54] = REAL(theta)[22];
      REAL(ret)[55] = REAL(theta)[29];
      REAL(ret)[56] = REAL(theta)[37];
      REAL(ret)[57] = REAL(theta)[46];
      REAL(ret)[58] = REAL(theta)[56];
      REAL(ret)[59] = REAL(theta)[67];
      REAL(ret)[64] = REAL(theta)[16];
      REAL(ret)[76] = REAL(theta)[22];
      REAL(ret)[88] = REAL(theta)[29];
      REAL(ret)[100] = REAL(theta)[37];
      REAL(ret)[112] = REAL(theta)[46];
      REAL(ret)[124] = REAL(theta)[56];
      REAL(ret)[136] = REAL(theta)[67];
    }
    else if (theta_n == 13){
      REAL(ret)[28] = exp(REAL(theta)[5]);
      REAL(ret)[40] = REAL(theta)[8];
      REAL(ret)[50] = exp(REAL(theta)[5]);
      REAL(ret)[51] = REAL(theta)[8];
      REAL(ret)[52] = 2 * REAL(theta)[12];
      REAL(ret)[53] = REAL(theta)[17];
      REAL(ret)[54] = REAL(theta)[23];
      REAL(ret)[55] = REAL(theta)[30];
      REAL(ret)[56] = REAL(theta)[38];
      REAL(ret)[57] = REAL(theta)[47];
      REAL(ret)[58] = REAL(theta)[57];
      REAL(ret)[59] = REAL(theta)[68];
      REAL(ret)[64] = REAL(theta)[17];
      REAL(ret)[76] = REAL(theta)[23];
      REAL(ret)[88] = REAL(theta)[30];
      REAL(ret)[100] = REAL(theta)[38];
      REAL(ret)[112] = REAL(theta)[47];
      REAL(ret)[124] = REAL(theta)[57];
      REAL(ret)[136] = REAL(theta)[68];
    }
    else if (theta_n == 14){
      REAL(ret)[40] = exp(REAL(theta)[9]);
      REAL(ret)[51] = exp(REAL(theta)[9]);
      REAL(ret)[52] = 2 * REAL(theta)[13];
      REAL(ret)[53] = REAL(theta)[18];
      REAL(ret)[54] = REAL(theta)[24];
      REAL(ret)[55] = REAL(theta)[31];
      REAL(ret)[56] = REAL(theta)[39];
      REAL(ret)[57] = REAL(theta)[48];
      REAL(ret)[58] = REAL(theta)[58];
      REAL(ret)[59] = REAL(theta)[69];
      REAL(ret)[64] = REAL(theta)[18];
      REAL(ret)[76] = REAL(theta)[24];
      REAL(ret)[88] = REAL(theta)[31];
      REAL(ret)[100] = REAL(theta)[39];
      REAL(ret)[112] = REAL(theta)[48];
      REAL(ret)[124] = REAL(theta)[58];
      REAL(ret)[136] = REAL(theta)[69];
    }
    else if (theta_n == 15){
      REAL(ret)[52] = 2 * exp(2 * REAL(theta)[14]);
      REAL(ret)[53] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[54] = REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[55] = REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[56] = REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[57] = REAL(theta)[49] * exp(REAL(theta)[14]);
      REAL(ret)[58] = REAL(theta)[59] * exp(REAL(theta)[14]);
      REAL(ret)[59] = REAL(theta)[70] * exp(REAL(theta)[14]);
      REAL(ret)[64] = REAL(theta)[19] * exp(REAL(theta)[14]);
      REAL(ret)[76] = REAL(theta)[25] * exp(REAL(theta)[14]);
      REAL(ret)[88] = REAL(theta)[32] * exp(REAL(theta)[14]);
      REAL(ret)[100] = REAL(theta)[40] * exp(REAL(theta)[14]);
      REAL(ret)[112] = REAL(theta)[49] * exp(REAL(theta)[14]);
      REAL(ret)[124] = REAL(theta)[59] * exp(REAL(theta)[14]);
      REAL(ret)[136] = REAL(theta)[70] * exp(REAL(theta)[14]);
    }
    else if (theta_n == 16){
      REAL(ret)[5] = exp(REAL(theta)[0]);
      REAL(ret)[17] = REAL(theta)[1];
      REAL(ret)[29] = REAL(theta)[3];
      REAL(ret)[41] = REAL(theta)[6];
      REAL(ret)[53] = REAL(theta)[10];
      REAL(ret)[60] = exp(REAL(theta)[0]);
      REAL(ret)[61] = REAL(theta)[1];
      REAL(ret)[62] = REAL(theta)[3];
      REAL(ret)[63] = REAL(theta)[6];
      REAL(ret)[64] = REAL(theta)[10];
      REAL(ret)[65] = 2 * REAL(theta)[15];
      REAL(ret)[66] = REAL(theta)[21];
      REAL(ret)[67] = REAL(theta)[28];
      REAL(ret)[68] = REAL(theta)[36];
      REAL(ret)[69] = REAL(theta)[45];
      REAL(ret)[70] = REAL(theta)[55];
      REAL(ret)[71] = REAL(theta)[66];
      REAL(ret)[77] = REAL(theta)[21];
      REAL(ret)[89] = REAL(theta)[28];
      REAL(ret)[101] = REAL(theta)[36];
      REAL(ret)[113] = REAL(theta)[45];
      REAL(ret)[125] = REAL(theta)[55];
      REAL(ret)[137] = REAL(theta)[66];
    }
    else if (theta_n == 17){
      REAL(ret)[17] = exp(REAL(theta)[2]);
      REAL(ret)[29] = REAL(theta)[4];
      REAL(ret)[41] = REAL(theta)[7];
      REAL(ret)[53] = REAL(theta)[11];
      REAL(ret)[61] = exp(REAL(theta)[2]);
      REAL(ret)[62] = REAL(theta)[4];
      REAL(ret)[63] = REAL(theta)[7];
      REAL(ret)[64] = REAL(theta)[11];
      REAL(ret)[65] = 2 * REAL(theta)[16];
      REAL(ret)[66] = REAL(theta)[22];
      REAL(ret)[67] = REAL(theta)[29];
      REAL(ret)[68] = REAL(theta)[37];
      REAL(ret)[69] = REAL(theta)[46];
      REAL(ret)[70] = REAL(theta)[56];
      REAL(ret)[71] = REAL(theta)[67];
      REAL(ret)[77] = REAL(theta)[22];
      REAL(ret)[89] = REAL(theta)[29];
      REAL(ret)[101] = REAL(theta)[37];
      REAL(ret)[113] = REAL(theta)[46];
      REAL(ret)[125] = REAL(theta)[56];
      REAL(ret)[137] = REAL(theta)[67];
    }
    else if (theta_n == 18){
      REAL(ret)[29] = exp(REAL(theta)[5]);
      REAL(ret)[41] = REAL(theta)[8];
      REAL(ret)[53] = REAL(theta)[12];
      REAL(ret)[62] = exp(REAL(theta)[5]);
      REAL(ret)[63] = REAL(theta)[8];
      REAL(ret)[64] = REAL(theta)[12];
      REAL(ret)[65] = 2 * REAL(theta)[17];
      REAL(ret)[66] = REAL(theta)[23];
      REAL(ret)[67] = REAL(theta)[30];
      REAL(ret)[68] = REAL(theta)[38];
      REAL(ret)[69] = REAL(theta)[47];
      REAL(ret)[70] = REAL(theta)[57];
      REAL(ret)[71] = REAL(theta)[68];
      REAL(ret)[77] = REAL(theta)[23];
      REAL(ret)[89] = REAL(theta)[30];
      REAL(ret)[101] = REAL(theta)[38];
      REAL(ret)[113] = REAL(theta)[47];
      REAL(ret)[125] = REAL(theta)[57];
      REAL(ret)[137] = REAL(theta)[68];
    }
    else if (theta_n == 19){
      REAL(ret)[41] = exp(REAL(theta)[9]);
      REAL(ret)[53] = REAL(theta)[13];
      REAL(ret)[63] = exp(REAL(theta)[9]);
      REAL(ret)[64] = REAL(theta)[13];
      REAL(ret)[65] = 2 * REAL(theta)[18];
      REAL(ret)[66] = REAL(theta)[24];
      REAL(ret)[67] = REAL(theta)[31];
      REAL(ret)[68] = REAL(theta)[39];
      REAL(ret)[69] = REAL(theta)[48];
      REAL(ret)[70] = REAL(theta)[58];
      REAL(ret)[71] = REAL(theta)[69];
      REAL(ret)[77] = REAL(theta)[24];
      REAL(ret)[89] = REAL(theta)[31];
      REAL(ret)[101] = REAL(theta)[39];
      REAL(ret)[113] = REAL(theta)[48];
      REAL(ret)[125] = REAL(theta)[58];
      REAL(ret)[137] = REAL(theta)[69];
    }
    else if (theta_n == 20){
      REAL(ret)[53] = exp(REAL(theta)[14]);
      REAL(ret)[64] = exp(REAL(theta)[14]);
      REAL(ret)[65] = 2 * REAL(theta)[19];
      REAL(ret)[66] = REAL(theta)[25];
      REAL(ret)[67] = REAL(theta)[32];
      REAL(ret)[68] = REAL(theta)[40];
      REAL(ret)[69] = REAL(theta)[49];
      REAL(ret)[70] = REAL(theta)[59];
      REAL(ret)[71] = REAL(theta)[70];
      REAL(ret)[77] = REAL(theta)[25];
      REAL(ret)[89] = REAL(theta)[32];
      REAL(ret)[101] = REAL(theta)[40];
      REAL(ret)[113] = REAL(theta)[49];
      REAL(ret)[125] = REAL(theta)[59];
      REAL(ret)[137] = REAL(theta)[70];
    }
    else if (theta_n == 21){
      REAL(ret)[65] = 2 * exp(2 * REAL(theta)[20]);
      REAL(ret)[66] = REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[67] = REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[68] = REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[69] = REAL(theta)[50] * exp(REAL(theta)[20]);
      REAL(ret)[70] = REAL(theta)[60] * exp(REAL(theta)[20]);
      REAL(ret)[71] = REAL(theta)[71] * exp(REAL(theta)[20]);
      REAL(ret)[77] = REAL(theta)[26] * exp(REAL(theta)[20]);
      REAL(ret)[89] = REAL(theta)[33] * exp(REAL(theta)[20]);
      REAL(ret)[101] = REAL(theta)[41] * exp(REAL(theta)[20]);
      REAL(ret)[113] = REAL(theta)[50] * exp(REAL(theta)[20]);
      REAL(ret)[125] = REAL(theta)[60] * exp(REAL(theta)[20]);
      REAL(ret)[137] = REAL(theta)[71] * exp(REAL(theta)[20]);
    }
    else if (theta_n == 22){
      REAL(ret)[6] = exp(REAL(theta)[0]);
      REAL(ret)[18] = REAL(theta)[1];
      REAL(ret)[30] = REAL(theta)[3];
      REAL(ret)[42] = REAL(theta)[6];
      REAL(ret)[54] = REAL(theta)[10];
      REAL(ret)[66] = REAL(theta)[15];
      REAL(ret)[72] = exp(REAL(theta)[0]);
      REAL(ret)[73] = REAL(theta)[1];
      REAL(ret)[74] = REAL(theta)[3];
      REAL(ret)[75] = REAL(theta)[6];
      REAL(ret)[76] = REAL(theta)[10];
      REAL(ret)[77] = REAL(theta)[15];
      REAL(ret)[78] = 2 * REAL(theta)[21];
      REAL(ret)[79] = REAL(theta)[28];
      REAL(ret)[80] = REAL(theta)[36];
      REAL(ret)[81] = REAL(theta)[45];
      REAL(ret)[82] = REAL(theta)[55];
      REAL(ret)[83] = REAL(theta)[66];
      REAL(ret)[90] = REAL(theta)[28];
      REAL(ret)[102] = REAL(theta)[36];
      REAL(ret)[114] = REAL(theta)[45];
      REAL(ret)[126] = REAL(theta)[55];
      REAL(ret)[138] = REAL(theta)[66];
    }
    else if (theta_n == 23){
      REAL(ret)[18] = exp(REAL(theta)[2]);
      REAL(ret)[30] = REAL(theta)[4];
      REAL(ret)[42] = REAL(theta)[7];
      REAL(ret)[54] = REAL(theta)[11];
      REAL(ret)[66] = REAL(theta)[16];
      REAL(ret)[73] = exp(REAL(theta)[2]);
      REAL(ret)[74] = REAL(theta)[4];
      REAL(ret)[75] = REAL(theta)[7];
      REAL(ret)[76] = REAL(theta)[11];
      REAL(ret)[77] = REAL(theta)[16];
      REAL(ret)[78] = 2 * REAL(theta)[22];
      REAL(ret)[79] = REAL(theta)[29];
      REAL(ret)[80] = REAL(theta)[37];
      REAL(ret)[81] = REAL(theta)[46];
      REAL(ret)[82] = REAL(theta)[56];
      REAL(ret)[83] = REAL(theta)[67];
      REAL(ret)[90] = REAL(theta)[29];
      REAL(ret)[102] = REAL(theta)[37];
      REAL(ret)[114] = REAL(theta)[46];
      REAL(ret)[126] = REAL(theta)[56];
      REAL(ret)[138] = REAL(theta)[67];
    }
    else if (theta_n == 24){
      REAL(ret)[30] = exp(REAL(theta)[5]);
      REAL(ret)[42] = REAL(theta)[8];
      REAL(ret)[54] = REAL(theta)[12];
      REAL(ret)[66] = REAL(theta)[17];
      REAL(ret)[74] = exp(REAL(theta)[5]);
      REAL(ret)[75] = REAL(theta)[8];
      REAL(ret)[76] = REAL(theta)[12];
      REAL(ret)[77] = REAL(theta)[17];
      REAL(ret)[78] = 2 * REAL(theta)[23];
      REAL(ret)[79] = REAL(theta)[30];
      REAL(ret)[80] = REAL(theta)[38];
      REAL(ret)[81] = REAL(theta)[47];
      REAL(ret)[82] = REAL(theta)[57];
      REAL(ret)[83] = REAL(theta)[68];
      REAL(ret)[90] = REAL(theta)[30];
      REAL(ret)[102] = REAL(theta)[38];
      REAL(ret)[114] = REAL(theta)[47];
      REAL(ret)[126] = REAL(theta)[57];
      REAL(ret)[138] = REAL(theta)[68];
    }
    else if (theta_n == 25){
      REAL(ret)[42] = exp(REAL(theta)[9]);
      REAL(ret)[54] = REAL(theta)[13];
      REAL(ret)[66] = REAL(theta)[18];
      REAL(ret)[75] = exp(REAL(theta)[9]);
      REAL(ret)[76] = REAL(theta)[13];
      REAL(ret)[77] = REAL(theta)[18];
      REAL(ret)[78] = 2 * REAL(theta)[24];
      REAL(ret)[79] = REAL(theta)[31];
      REAL(ret)[80] = REAL(theta)[39];
      REAL(ret)[81] = REAL(theta)[48];
      REAL(ret)[82] = REAL(theta)[58];
      REAL(ret)[83] = REAL(theta)[69];
      REAL(ret)[90] = REAL(theta)[31];
      REAL(ret)[102] = REAL(theta)[39];
      REAL(ret)[114] = REAL(theta)[48];
      REAL(ret)[126] = REAL(theta)[58];
      REAL(ret)[138] = REAL(theta)[69];
    }
    else if (theta_n == 26){
      REAL(ret)[54] = exp(REAL(theta)[14]);
      REAL(ret)[66] = REAL(theta)[19];
      REAL(ret)[76] = exp(REAL(theta)[14]);
      REAL(ret)[77] = REAL(theta)[19];
      REAL(ret)[78] = 2 * REAL(theta)[25];
      REAL(ret)[79] = REAL(theta)[32];
      REAL(ret)[80] = REAL(theta)[40];
      REAL(ret)[81] = REAL(theta)[49];
      REAL(ret)[82] = REAL(theta)[59];
      REAL(ret)[83] = REAL(theta)[70];
      REAL(ret)[90] = REAL(theta)[32];
      REAL(ret)[102] = REAL(theta)[40];
      REAL(ret)[114] = REAL(theta)[49];
      REAL(ret)[126] = REAL(theta)[59];
      REAL(ret)[138] = REAL(theta)[70];
    }
    else if (theta_n == 27){
      REAL(ret)[66] = exp(REAL(theta)[20]);
      REAL(ret)[77] = exp(REAL(theta)[20]);
      REAL(ret)[78] = 2 * REAL(theta)[26];
      REAL(ret)[79] = REAL(theta)[33];
      REAL(ret)[80] = REAL(theta)[41];
      REAL(ret)[81] = REAL(theta)[50];
      REAL(ret)[82] = REAL(theta)[60];
      REAL(ret)[83] = REAL(theta)[71];
      REAL(ret)[90] = REAL(theta)[33];
      REAL(ret)[102] = REAL(theta)[41];
      REAL(ret)[114] = REAL(theta)[50];
      REAL(ret)[126] = REAL(theta)[60];
      REAL(ret)[138] = REAL(theta)[71];
    }
    else if (theta_n == 28){
      REAL(ret)[78] = 2 * exp(2 * REAL(theta)[27]);
      REAL(ret)[79] = REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[80] = REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[81] = REAL(theta)[51] * exp(REAL(theta)[27]);
      REAL(ret)[82] = REAL(theta)[61] * exp(REAL(theta)[27]);
      REAL(ret)[83] = REAL(theta)[72] * exp(REAL(theta)[27]);
      REAL(ret)[90] = REAL(theta)[34] * exp(REAL(theta)[27]);
      REAL(ret)[102] = REAL(theta)[42] * exp(REAL(theta)[27]);
      REAL(ret)[114] = REAL(theta)[51] * exp(REAL(theta)[27]);
      REAL(ret)[126] = REAL(theta)[61] * exp(REAL(theta)[27]);
      REAL(ret)[138] = REAL(theta)[72] * exp(REAL(theta)[27]);
    }
    else if (theta_n == 29){
      REAL(ret)[7] = exp(REAL(theta)[0]);
      REAL(ret)[19] = REAL(theta)[1];
      REAL(ret)[31] = REAL(theta)[3];
      REAL(ret)[43] = REAL(theta)[6];
      REAL(ret)[55] = REAL(theta)[10];
      REAL(ret)[67] = REAL(theta)[15];
      REAL(ret)[79] = REAL(theta)[21];
      REAL(ret)[84] = exp(REAL(theta)[0]);
      REAL(ret)[85] = REAL(theta)[1];
      REAL(ret)[86] = REAL(theta)[3];
      REAL(ret)[87] = REAL(theta)[6];
      REAL(ret)[88] = REAL(theta)[10];
      REAL(ret)[89] = REAL(theta)[15];
      REAL(ret)[90] = REAL(theta)[21];
      REAL(ret)[91] = 2 * REAL(theta)[28];
      REAL(ret)[92] = REAL(theta)[36];
      REAL(ret)[93] = REAL(theta)[45];
      REAL(ret)[94] = REAL(theta)[55];
      REAL(ret)[95] = REAL(theta)[66];
      REAL(ret)[103] = REAL(theta)[36];
      REAL(ret)[115] = REAL(theta)[45];
      REAL(ret)[127] = REAL(theta)[55];
      REAL(ret)[139] = REAL(theta)[66];
    }
    else if (theta_n == 30){
      REAL(ret)[19] = exp(REAL(theta)[2]);
      REAL(ret)[31] = REAL(theta)[4];
      REAL(ret)[43] = REAL(theta)[7];
      REAL(ret)[55] = REAL(theta)[11];
      REAL(ret)[67] = REAL(theta)[16];
      REAL(ret)[79] = REAL(theta)[22];
      REAL(ret)[85] = exp(REAL(theta)[2]);
      REAL(ret)[86] = REAL(theta)[4];
      REAL(ret)[87] = REAL(theta)[7];
      REAL(ret)[88] = REAL(theta)[11];
      REAL(ret)[89] = REAL(theta)[16];
      REAL(ret)[90] = REAL(theta)[22];
      REAL(ret)[91] = 2 * REAL(theta)[29];
      REAL(ret)[92] = REAL(theta)[37];
      REAL(ret)[93] = REAL(theta)[46];
      REAL(ret)[94] = REAL(theta)[56];
      REAL(ret)[95] = REAL(theta)[67];
      REAL(ret)[103] = REAL(theta)[37];
      REAL(ret)[115] = REAL(theta)[46];
      REAL(ret)[127] = REAL(theta)[56];
      REAL(ret)[139] = REAL(theta)[67];
    }
    else if (theta_n == 31){
      REAL(ret)[31] = exp(REAL(theta)[5]);
      REAL(ret)[43] = REAL(theta)[8];
      REAL(ret)[55] = REAL(theta)[12];
      REAL(ret)[67] = REAL(theta)[17];
      REAL(ret)[79] = REAL(theta)[23];
      REAL(ret)[86] = exp(REAL(theta)[5]);
      REAL(ret)[87] = REAL(theta)[8];
      REAL(ret)[88] = REAL(theta)[12];
      REAL(ret)[89] = REAL(theta)[17];
      REAL(ret)[90] = REAL(theta)[23];
      REAL(ret)[91] = 2 * REAL(theta)[30];
      REAL(ret)[92] = REAL(theta)[38];
      REAL(ret)[93] = REAL(theta)[47];
      REAL(ret)[94] = REAL(theta)[57];
      REAL(ret)[95] = REAL(theta)[68];
      REAL(ret)[103] = REAL(theta)[38];
      REAL(ret)[115] = REAL(theta)[47];
      REAL(ret)[127] = REAL(theta)[57];
      REAL(ret)[139] = REAL(theta)[68];
    }
    else if (theta_n == 32){
      REAL(ret)[43] = exp(REAL(theta)[9]);
      REAL(ret)[55] = REAL(theta)[13];
      REAL(ret)[67] = REAL(theta)[18];
      REAL(ret)[79] = REAL(theta)[24];
      REAL(ret)[87] = exp(REAL(theta)[9]);
      REAL(ret)[88] = REAL(theta)[13];
      REAL(ret)[89] = REAL(theta)[18];
      REAL(ret)[90] = REAL(theta)[24];
      REAL(ret)[91] = 2 * REAL(theta)[31];
      REAL(ret)[92] = REAL(theta)[39];
      REAL(ret)[93] = REAL(theta)[48];
      REAL(ret)[94] = REAL(theta)[58];
      REAL(ret)[95] = REAL(theta)[69];
      REAL(ret)[103] = REAL(theta)[39];
      REAL(ret)[115] = REAL(theta)[48];
      REAL(ret)[127] = REAL(theta)[58];
      REAL(ret)[139] = REAL(theta)[69];
    }
    else if (theta_n == 33){
      REAL(ret)[55] = exp(REAL(theta)[14]);
      REAL(ret)[67] = REAL(theta)[19];
      REAL(ret)[79] = REAL(theta)[25];
      REAL(ret)[88] = exp(REAL(theta)[14]);
      REAL(ret)[89] = REAL(theta)[19];
      REAL(ret)[90] = REAL(theta)[25];
      REAL(ret)[91] = 2 * REAL(theta)[32];
      REAL(ret)[92] = REAL(theta)[40];
      REAL(ret)[93] = REAL(theta)[49];
      REAL(ret)[94] = REAL(theta)[59];
      REAL(ret)[95] = REAL(theta)[70];
      REAL(ret)[103] = REAL(theta)[40];
      REAL(ret)[115] = REAL(theta)[49];
      REAL(ret)[127] = REAL(theta)[59];
      REAL(ret)[139] = REAL(theta)[70];
    }
    else if (theta_n == 34){
      REAL(ret)[67] = exp(REAL(theta)[20]);
      REAL(ret)[79] = REAL(theta)[26];
      REAL(ret)[89] = exp(REAL(theta)[20]);
      REAL(ret)[90] = REAL(theta)[26];
      REAL(ret)[91] = 2 * REAL(theta)[33];
      REAL(ret)[92] = REAL(theta)[41];
      REAL(ret)[93] = REAL(theta)[50];
      REAL(ret)[94] = REAL(theta)[60];
      REAL(ret)[95] = REAL(theta)[71];
      REAL(ret)[103] = REAL(theta)[41];
      REAL(ret)[115] = REAL(theta)[50];
      REAL(ret)[127] = REAL(theta)[60];
      REAL(ret)[139] = REAL(theta)[71];
    }
    else if (theta_n == 35){
      REAL(ret)[79] = exp(REAL(theta)[27]);
      REAL(ret)[90] = exp(REAL(theta)[27]);
      REAL(ret)[91] = 2 * REAL(theta)[34];
      REAL(ret)[92] = REAL(theta)[42];
      REAL(ret)[93] = REAL(theta)[51];
      REAL(ret)[94] = REAL(theta)[61];
      REAL(ret)[95] = REAL(theta)[72];
      REAL(ret)[103] = REAL(theta)[42];
      REAL(ret)[115] = REAL(theta)[51];
      REAL(ret)[127] = REAL(theta)[61];
      REAL(ret)[139] = REAL(theta)[72];
    }
    else if (theta_n == 36){
      REAL(ret)[91] = 2 * exp(2 * REAL(theta)[35]);
      REAL(ret)[92] = REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[93] = REAL(theta)[52] * exp(REAL(theta)[35]);
      REAL(ret)[94] = REAL(theta)[62] * exp(REAL(theta)[35]);
      REAL(ret)[95] = REAL(theta)[73] * exp(REAL(theta)[35]);
      REAL(ret)[103] = REAL(theta)[43] * exp(REAL(theta)[35]);
      REAL(ret)[115] = REAL(theta)[52] * exp(REAL(theta)[35]);
      REAL(ret)[127] = REAL(theta)[62] * exp(REAL(theta)[35]);
      REAL(ret)[139] = REAL(theta)[73] * exp(REAL(theta)[35]);
    }
    else if (theta_n == 37){
      REAL(ret)[8] = exp(REAL(theta)[0]);
      REAL(ret)[20] = REAL(theta)[1];
      REAL(ret)[32] = REAL(theta)[3];
      REAL(ret)[44] = REAL(theta)[6];
      REAL(ret)[56] = REAL(theta)[10];
      REAL(ret)[68] = REAL(theta)[15];
      REAL(ret)[80] = REAL(theta)[21];
      REAL(ret)[92] = REAL(theta)[28];
      REAL(ret)[96] = exp(REAL(theta)[0]);
      REAL(ret)[97] = REAL(theta)[1];
      REAL(ret)[98] = REAL(theta)[3];
      REAL(ret)[99] = REAL(theta)[6];
      REAL(ret)[100] = REAL(theta)[10];
      REAL(ret)[101] = REAL(theta)[15];
      REAL(ret)[102] = REAL(theta)[21];
      REAL(ret)[103] = REAL(theta)[28];
      REAL(ret)[104] = 2 * REAL(theta)[36];
      REAL(ret)[105] = REAL(theta)[45];
      REAL(ret)[106] = REAL(theta)[55];
      REAL(ret)[107] = REAL(theta)[66];
      REAL(ret)[116] = REAL(theta)[45];
      REAL(ret)[128] = REAL(theta)[55];
      REAL(ret)[140] = REAL(theta)[66];
    }
    else if (theta_n == 38){
      REAL(ret)[20] = exp(REAL(theta)[2]);
      REAL(ret)[32] = REAL(theta)[4];
      REAL(ret)[44] = REAL(theta)[7];
      REAL(ret)[56] = REAL(theta)[11];
      REAL(ret)[68] = REAL(theta)[16];
      REAL(ret)[80] = REAL(theta)[22];
      REAL(ret)[92] = REAL(theta)[29];
      REAL(ret)[97] = exp(REAL(theta)[2]);
      REAL(ret)[98] = REAL(theta)[4];
      REAL(ret)[99] = REAL(theta)[7];
      REAL(ret)[100] = REAL(theta)[11];
      REAL(ret)[101] = REAL(theta)[16];
      REAL(ret)[102] = REAL(theta)[22];
      REAL(ret)[103] = REAL(theta)[29];
      REAL(ret)[104] = 2 * REAL(theta)[37];
      REAL(ret)[105] = REAL(theta)[46];
      REAL(ret)[106] = REAL(theta)[56];
      REAL(ret)[107] = REAL(theta)[67];
      REAL(ret)[116] = REAL(theta)[46];
      REAL(ret)[128] = REAL(theta)[56];
      REAL(ret)[140] = REAL(theta)[67];
    }
    else if (theta_n == 39){
      REAL(ret)[32] = exp(REAL(theta)[5]);
      REAL(ret)[44] = REAL(theta)[8];
      REAL(ret)[56] = REAL(theta)[12];
      REAL(ret)[68] = REAL(theta)[17];
      REAL(ret)[80] = REAL(theta)[23];
      REAL(ret)[92] = REAL(theta)[30];
      REAL(ret)[98] = exp(REAL(theta)[5]);
      REAL(ret)[99] = REAL(theta)[8];
      REAL(ret)[100] = REAL(theta)[12];
      REAL(ret)[101] = REAL(theta)[17];
      REAL(ret)[102] = REAL(theta)[23];
      REAL(ret)[103] = REAL(theta)[30];
      REAL(ret)[104] = 2 * REAL(theta)[38];
      REAL(ret)[105] = REAL(theta)[47];
      REAL(ret)[106] = REAL(theta)[57];
      REAL(ret)[107] = REAL(theta)[68];
      REAL(ret)[116] = REAL(theta)[47];
      REAL(ret)[128] = REAL(theta)[57];
      REAL(ret)[140] = REAL(theta)[68];
    }
    else if (theta_n == 40){
      REAL(ret)[44] = exp(REAL(theta)[9]);
      REAL(ret)[56] = REAL(theta)[13];
      REAL(ret)[68] = REAL(theta)[18];
      REAL(ret)[80] = REAL(theta)[24];
      REAL(ret)[92] = REAL(theta)[31];
      REAL(ret)[99] = exp(REAL(theta)[9]);
      REAL(ret)[100] = REAL(theta)[13];
      REAL(ret)[101] = REAL(theta)[18];
      REAL(ret)[102] = REAL(theta)[24];
      REAL(ret)[103] = REAL(theta)[31];
      REAL(ret)[104] = 2 * REAL(theta)[39];
      REAL(ret)[105] = REAL(theta)[48];
      REAL(ret)[106] = REAL(theta)[58];
      REAL(ret)[107] = REAL(theta)[69];
      REAL(ret)[116] = REAL(theta)[48];
      REAL(ret)[128] = REAL(theta)[58];
      REAL(ret)[140] = REAL(theta)[69];
    }
    else if (theta_n == 41){
      REAL(ret)[56] = exp(REAL(theta)[14]);
      REAL(ret)[68] = REAL(theta)[19];
      REAL(ret)[80] = REAL(theta)[25];
      REAL(ret)[92] = REAL(theta)[32];
      REAL(ret)[100] = exp(REAL(theta)[14]);
      REAL(ret)[101] = REAL(theta)[19];
      REAL(ret)[102] = REAL(theta)[25];
      REAL(ret)[103] = REAL(theta)[32];
      REAL(ret)[104] = 2 * REAL(theta)[40];
      REAL(ret)[105] = REAL(theta)[49];
      REAL(ret)[106] = REAL(theta)[59];
      REAL(ret)[107] = REAL(theta)[70];
      REAL(ret)[116] = REAL(theta)[49];
      REAL(ret)[128] = REAL(theta)[59];
      REAL(ret)[140] = REAL(theta)[70];
    }
    else if (theta_n == 42){
      REAL(ret)[68] = exp(REAL(theta)[20]);
      REAL(ret)[80] = REAL(theta)[26];
      REAL(ret)[92] = REAL(theta)[33];
      REAL(ret)[101] = exp(REAL(theta)[20]);
      REAL(ret)[102] = REAL(theta)[26];
      REAL(ret)[103] = REAL(theta)[33];
      REAL(ret)[104] = 2 * REAL(theta)[41];
      REAL(ret)[105] = REAL(theta)[50];
      REAL(ret)[106] = REAL(theta)[60];
      REAL(ret)[107] = REAL(theta)[71];
      REAL(ret)[116] = REAL(theta)[50];
      REAL(ret)[128] = REAL(theta)[60];
      REAL(ret)[140] = REAL(theta)[71];
    }
    else if (theta_n == 43){
      REAL(ret)[80] = exp(REAL(theta)[27]);
      REAL(ret)[92] = REAL(theta)[34];
      REAL(ret)[102] = exp(REAL(theta)[27]);
      REAL(ret)[103] = REAL(theta)[34];
      REAL(ret)[104] = 2 * REAL(theta)[42];
      REAL(ret)[105] = REAL(theta)[51];
      REAL(ret)[106] = REAL(theta)[61];
      REAL(ret)[107] = REAL(theta)[72];
      REAL(ret)[116] = REAL(theta)[51];
      REAL(ret)[128] = REAL(theta)[61];
      REAL(ret)[140] = REAL(theta)[72];
    }
    else if (theta_n == 44){
      REAL(ret)[92] = exp(REAL(theta)[35]);
      REAL(ret)[103] = exp(REAL(theta)[35]);
      REAL(ret)[104] = 2 * REAL(theta)[43];
      REAL(ret)[105] = REAL(theta)[52];
      REAL(ret)[106] = REAL(theta)[62];
      REAL(ret)[107] = REAL(theta)[73];
      REAL(ret)[116] = REAL(theta)[52];
      REAL(ret)[128] = REAL(theta)[62];
      REAL(ret)[140] = REAL(theta)[73];
    }
    else if (theta_n == 45){
      REAL(ret)[104] = 2 * exp(2 * REAL(theta)[44]);
      REAL(ret)[105] = REAL(theta)[53] * exp(REAL(theta)[44]);
      REAL(ret)[106] = REAL(theta)[63] * exp(REAL(theta)[44]);
      REAL(ret)[107] = REAL(theta)[74] * exp(REAL(theta)[44]);
      REAL(ret)[116] = REAL(theta)[53] * exp(REAL(theta)[44]);
      REAL(ret)[128] = REAL(theta)[63] * exp(REAL(theta)[44]);
      REAL(ret)[140] = REAL(theta)[74] * exp(REAL(theta)[44]);
    }
    else if (theta_n == 46){
      REAL(ret)[9] = exp(REAL(theta)[0]);
      REAL(ret)[21] = REAL(theta)[1];
      REAL(ret)[33] = REAL(theta)[3];
      REAL(ret)[45] = REAL(theta)[6];
      REAL(ret)[57] = REAL(theta)[10];
      REAL(ret)[69] = REAL(theta)[15];
      REAL(ret)[81] = REAL(theta)[21];
      REAL(ret)[93] = REAL(theta)[28];
      REAL(ret)[105] = REAL(theta)[36];
      REAL(ret)[108] = exp(REAL(theta)[0]);
      REAL(ret)[109] = REAL(theta)[1];
      REAL(ret)[110] = REAL(theta)[3];
      REAL(ret)[111] = REAL(theta)[6];
      REAL(ret)[112] = REAL(theta)[10];
      REAL(ret)[113] = REAL(theta)[15];
      REAL(ret)[114] = REAL(theta)[21];
      REAL(ret)[115] = REAL(theta)[28];
      REAL(ret)[116] = REAL(theta)[36];
      REAL(ret)[117] = 2 * REAL(theta)[45];
      REAL(ret)[118] = REAL(theta)[55];
      REAL(ret)[119] = REAL(theta)[66];
      REAL(ret)[129] = REAL(theta)[55];
      REAL(ret)[141] = REAL(theta)[66];
    }
    else if (theta_n == 47){
      REAL(ret)[21] = exp(REAL(theta)[2]);
      REAL(ret)[33] = REAL(theta)[4];
      REAL(ret)[45] = REAL(theta)[7];
      REAL(ret)[57] = REAL(theta)[11];
      REAL(ret)[69] = REAL(theta)[16];
      REAL(ret)[81] = REAL(theta)[22];
      REAL(ret)[93] = REAL(theta)[29];
      REAL(ret)[105] = REAL(theta)[37];
      REAL(ret)[109] = exp(REAL(theta)[2]);
      REAL(ret)[110] = REAL(theta)[4];
      REAL(ret)[111] = REAL(theta)[7];
      REAL(ret)[112] = REAL(theta)[11];
      REAL(ret)[113] = REAL(theta)[16];
      REAL(ret)[114] = REAL(theta)[22];
      REAL(ret)[115] = REAL(theta)[29];
      REAL(ret)[116] = REAL(theta)[37];
      REAL(ret)[117] = 2 * REAL(theta)[46];
      REAL(ret)[118] = REAL(theta)[56];
      REAL(ret)[119] = REAL(theta)[67];
      REAL(ret)[129] = REAL(theta)[56];
      REAL(ret)[141] = REAL(theta)[67];
    }
    else if (theta_n == 48){
      REAL(ret)[33] = exp(REAL(theta)[5]);
      REAL(ret)[45] = REAL(theta)[8];
      REAL(ret)[57] = REAL(theta)[12];
      REAL(ret)[69] = REAL(theta)[17];
      REAL(ret)[81] = REAL(theta)[23];
      REAL(ret)[93] = REAL(theta)[30];
      REAL(ret)[105] = REAL(theta)[38];
      REAL(ret)[110] = exp(REAL(theta)[5]);
      REAL(ret)[111] = REAL(theta)[8];
      REAL(ret)[112] = REAL(theta)[12];
      REAL(ret)[113] = REAL(theta)[17];
      REAL(ret)[114] = REAL(theta)[23];
      REAL(ret)[115] = REAL(theta)[30];
      REAL(ret)[116] = REAL(theta)[38];
      REAL(ret)[117] = 2 * REAL(theta)[47];
      REAL(ret)[118] = REAL(theta)[57];
      REAL(ret)[119] = REAL(theta)[68];
      REAL(ret)[129] = REAL(theta)[57];
      REAL(ret)[141] = REAL(theta)[68];
    }
    else if (theta_n == 49){
      REAL(ret)[45] = exp(REAL(theta)[9]);
      REAL(ret)[57] = REAL(theta)[13];
      REAL(ret)[69] = REAL(theta)[18];
      REAL(ret)[81] = REAL(theta)[24];
      REAL(ret)[93] = REAL(theta)[31];
      REAL(ret)[105] = REAL(theta)[39];
      REAL(ret)[111] = exp(REAL(theta)[9]);
      REAL(ret)[112] = REAL(theta)[13];
      REAL(ret)[113] = REAL(theta)[18];
      REAL(ret)[114] = REAL(theta)[24];
      REAL(ret)[115] = REAL(theta)[31];
      REAL(ret)[116] = REAL(theta)[39];
      REAL(ret)[117] = 2 * REAL(theta)[48];
      REAL(ret)[118] = REAL(theta)[58];
      REAL(ret)[119] = REAL(theta)[69];
      REAL(ret)[129] = REAL(theta)[58];
      REAL(ret)[141] = REAL(theta)[69];
    }
    else if (theta_n == 50){
      REAL(ret)[57] = exp(REAL(theta)[14]);
      REAL(ret)[69] = REAL(theta)[19];
      REAL(ret)[81] = REAL(theta)[25];
      REAL(ret)[93] = REAL(theta)[32];
      REAL(ret)[105] = REAL(theta)[40];
      REAL(ret)[112] = exp(REAL(theta)[14]);
      REAL(ret)[113] = REAL(theta)[19];
      REAL(ret)[114] = REAL(theta)[25];
      REAL(ret)[115] = REAL(theta)[32];
      REAL(ret)[116] = REAL(theta)[40];
      REAL(ret)[117] = 2 * REAL(theta)[49];
      REAL(ret)[118] = REAL(theta)[59];
      REAL(ret)[119] = REAL(theta)[70];
      REAL(ret)[129] = REAL(theta)[59];
      REAL(ret)[141] = REAL(theta)[70];
    }
    else if (theta_n == 51){
      REAL(ret)[69] = exp(REAL(theta)[20]);
      REAL(ret)[81] = REAL(theta)[26];
      REAL(ret)[93] = REAL(theta)[33];
      REAL(ret)[105] = REAL(theta)[41];
      REAL(ret)[113] = exp(REAL(theta)[20]);
      REAL(ret)[114] = REAL(theta)[26];
      REAL(ret)[115] = REAL(theta)[33];
      REAL(ret)[116] = REAL(theta)[41];
      REAL(ret)[117] = 2 * REAL(theta)[50];
      REAL(ret)[118] = REAL(theta)[60];
      REAL(ret)[119] = REAL(theta)[71];
      REAL(ret)[129] = REAL(theta)[60];
      REAL(ret)[141] = REAL(theta)[71];
    }
    else if (theta_n == 52){
      REAL(ret)[81] = exp(REAL(theta)[27]);
      REAL(ret)[93] = REAL(theta)[34];
      REAL(ret)[105] = REAL(theta)[42];
      REAL(ret)[114] = exp(REAL(theta)[27]);
      REAL(ret)[115] = REAL(theta)[34];
      REAL(ret)[116] = REAL(theta)[42];
      REAL(ret)[117] = 2 * REAL(theta)[51];
      REAL(ret)[118] = REAL(theta)[61];
      REAL(ret)[119] = REAL(theta)[72];
      REAL(ret)[129] = REAL(theta)[61];
      REAL(ret)[141] = REAL(theta)[72];
    }
    else if (theta_n == 53){
      REAL(ret)[93] = exp(REAL(theta)[35]);
      REAL(ret)[105] = REAL(theta)[43];
      REAL(ret)[115] = exp(REAL(theta)[35]);
      REAL(ret)[116] = REAL(theta)[43];
      REAL(ret)[117] = 2 * REAL(theta)[52];
      REAL(ret)[118] = REAL(theta)[62];
      REAL(ret)[119] = REAL(theta)[73];
      REAL(ret)[129] = REAL(theta)[62];
      REAL(ret)[141] = REAL(theta)[73];
    }
    else if (theta_n == 54){
      REAL(ret)[105] = exp(REAL(theta)[44]);
      REAL(ret)[116] = exp(REAL(theta)[44]);
      REAL(ret)[117] = 2 * REAL(theta)[53];
      REAL(ret)[118] = REAL(theta)[63];
      REAL(ret)[119] = REAL(theta)[74];
      REAL(ret)[129] = REAL(theta)[63];
      REAL(ret)[141] = REAL(theta)[74];
    }
    else if (theta_n == 55){
      REAL(ret)[117] = 2 * exp(2 * REAL(theta)[54]);
      REAL(ret)[118] = REAL(theta)[64] * exp(REAL(theta)[54]);
      REAL(ret)[119] = REAL(theta)[75] * exp(REAL(theta)[54]);
      REAL(ret)[129] = REAL(theta)[64] * exp(REAL(theta)[54]);
      REAL(ret)[141] = REAL(theta)[75] * exp(REAL(theta)[54]);
    }
    else if (theta_n == 56){
      REAL(ret)[10] = exp(REAL(theta)[0]);
      REAL(ret)[22] = REAL(theta)[1];
      REAL(ret)[34] = REAL(theta)[3];
      REAL(ret)[46] = REAL(theta)[6];
      REAL(ret)[58] = REAL(theta)[10];
      REAL(ret)[70] = REAL(theta)[15];
      REAL(ret)[82] = REAL(theta)[21];
      REAL(ret)[94] = REAL(theta)[28];
      REAL(ret)[106] = REAL(theta)[36];
      REAL(ret)[118] = REAL(theta)[45];
      REAL(ret)[120] = exp(REAL(theta)[0]);
      REAL(ret)[121] = REAL(theta)[1];
      REAL(ret)[122] = REAL(theta)[3];
      REAL(ret)[123] = REAL(theta)[6];
      REAL(ret)[124] = REAL(theta)[10];
      REAL(ret)[125] = REAL(theta)[15];
      REAL(ret)[126] = REAL(theta)[21];
      REAL(ret)[127] = REAL(theta)[28];
      REAL(ret)[128] = REAL(theta)[36];
      REAL(ret)[129] = REAL(theta)[45];
      REAL(ret)[130] = 2 * REAL(theta)[55];
      REAL(ret)[131] = REAL(theta)[66];
      REAL(ret)[142] = REAL(theta)[66];
    }
    else if (theta_n == 57){
      REAL(ret)[22] = exp(REAL(theta)[2]);
      REAL(ret)[34] = REAL(theta)[4];
      REAL(ret)[46] = REAL(theta)[7];
      REAL(ret)[58] = REAL(theta)[11];
      REAL(ret)[70] = REAL(theta)[16];
      REAL(ret)[82] = REAL(theta)[22];
      REAL(ret)[94] = REAL(theta)[29];
      REAL(ret)[106] = REAL(theta)[37];
      REAL(ret)[118] = REAL(theta)[46];
      REAL(ret)[121] = exp(REAL(theta)[2]);
      REAL(ret)[122] = REAL(theta)[4];
      REAL(ret)[123] = REAL(theta)[7];
      REAL(ret)[124] = REAL(theta)[11];
      REAL(ret)[125] = REAL(theta)[16];
      REAL(ret)[126] = REAL(theta)[22];
      REAL(ret)[127] = REAL(theta)[29];
      REAL(ret)[128] = REAL(theta)[37];
      REAL(ret)[129] = REAL(theta)[46];
      REAL(ret)[130] = 2 * REAL(theta)[56];
      REAL(ret)[131] = REAL(theta)[67];
      REAL(ret)[142] = REAL(theta)[67];
    }
    else if (theta_n == 58){
      REAL(ret)[34] = exp(REAL(theta)[5]);
      REAL(ret)[46] = REAL(theta)[8];
      REAL(ret)[58] = REAL(theta)[12];
      REAL(ret)[70] = REAL(theta)[17];
      REAL(ret)[82] = REAL(theta)[23];
      REAL(ret)[94] = REAL(theta)[30];
      REAL(ret)[106] = REAL(theta)[38];
      REAL(ret)[118] = REAL(theta)[47];
      REAL(ret)[122] = exp(REAL(theta)[5]);
      REAL(ret)[123] = REAL(theta)[8];
      REAL(ret)[124] = REAL(theta)[12];
      REAL(ret)[125] = REAL(theta)[17];
      REAL(ret)[126] = REAL(theta)[23];
      REAL(ret)[127] = REAL(theta)[30];
      REAL(ret)[128] = REAL(theta)[38];
      REAL(ret)[129] = REAL(theta)[47];
      REAL(ret)[130] = 2 * REAL(theta)[57];
      REAL(ret)[131] = REAL(theta)[68];
      REAL(ret)[142] = REAL(theta)[68];
    }
    else if (theta_n == 59){
      REAL(ret)[46] = exp(REAL(theta)[9]);
      REAL(ret)[58] = REAL(theta)[13];
      REAL(ret)[70] = REAL(theta)[18];
      REAL(ret)[82] = REAL(theta)[24];
      REAL(ret)[94] = REAL(theta)[31];
      REAL(ret)[106] = REAL(theta)[39];
      REAL(ret)[118] = REAL(theta)[48];
      REAL(ret)[123] = exp(REAL(theta)[9]);
      REAL(ret)[124] = REAL(theta)[13];
      REAL(ret)[125] = REAL(theta)[18];
      REAL(ret)[126] = REAL(theta)[24];
      REAL(ret)[127] = REAL(theta)[31];
      REAL(ret)[128] = REAL(theta)[39];
      REAL(ret)[129] = REAL(theta)[48];
      REAL(ret)[130] = 2 * REAL(theta)[58];
      REAL(ret)[131] = REAL(theta)[69];
      REAL(ret)[142] = REAL(theta)[69];
    }
    else if (theta_n == 60){
      REAL(ret)[58] = exp(REAL(theta)[14]);
      REAL(ret)[70] = REAL(theta)[19];
      REAL(ret)[82] = REAL(theta)[25];
      REAL(ret)[94] = REAL(theta)[32];
      REAL(ret)[106] = REAL(theta)[40];
      REAL(ret)[118] = REAL(theta)[49];
      REAL(ret)[124] = exp(REAL(theta)[14]);
      REAL(ret)[125] = REAL(theta)[19];
      REAL(ret)[126] = REAL(theta)[25];
      REAL(ret)[127] = REAL(theta)[32];
      REAL(ret)[128] = REAL(theta)[40];
      REAL(ret)[129] = REAL(theta)[49];
      REAL(ret)[130] = 2 * REAL(theta)[59];
      REAL(ret)[131] = REAL(theta)[70];
      REAL(ret)[142] = REAL(theta)[70];
    }
    else if (theta_n == 61){
      REAL(ret)[70] = exp(REAL(theta)[20]);
      REAL(ret)[82] = REAL(theta)[26];
      REAL(ret)[94] = REAL(theta)[33];
      REAL(ret)[106] = REAL(theta)[41];
      REAL(ret)[118] = REAL(theta)[50];
      REAL(ret)[125] = exp(REAL(theta)[20]);
      REAL(ret)[126] = REAL(theta)[26];
      REAL(ret)[127] = REAL(theta)[33];
      REAL(ret)[128] = REAL(theta)[41];
      REAL(ret)[129] = REAL(theta)[50];
      REAL(ret)[130] = 2 * REAL(theta)[60];
      REAL(ret)[131] = REAL(theta)[71];
      REAL(ret)[142] = REAL(theta)[71];
    }
    else if (theta_n == 62){
      REAL(ret)[82] = exp(REAL(theta)[27]);
      REAL(ret)[94] = REAL(theta)[34];
      REAL(ret)[106] = REAL(theta)[42];
      REAL(ret)[118] = REAL(theta)[51];
      REAL(ret)[126] = exp(REAL(theta)[27]);
      REAL(ret)[127] = REAL(theta)[34];
      REAL(ret)[128] = REAL(theta)[42];
      REAL(ret)[129] = REAL(theta)[51];
      REAL(ret)[130] = 2 * REAL(theta)[61];
      REAL(ret)[131] = REAL(theta)[72];
      REAL(ret)[142] = REAL(theta)[72];
    }
    else if (theta_n == 63){
      REAL(ret)[94] = exp(REAL(theta)[35]);
      REAL(ret)[106] = REAL(theta)[43];
      REAL(ret)[118] = REAL(theta)[52];
      REAL(ret)[127] = exp(REAL(theta)[35]);
      REAL(ret)[128] = REAL(theta)[43];
      REAL(ret)[129] = REAL(theta)[52];
      REAL(ret)[130] = 2 * REAL(theta)[62];
      REAL(ret)[131] = REAL(theta)[73];
      REAL(ret)[142] = REAL(theta)[73];
    }
    else if (theta_n == 64){
      REAL(ret)[106] = exp(REAL(theta)[44]);
      REAL(ret)[118] = REAL(theta)[53];
      REAL(ret)[128] = exp(REAL(theta)[44]);
      REAL(ret)[129] = REAL(theta)[53];
      REAL(ret)[130] = 2 * REAL(theta)[63];
      REAL(ret)[131] = REAL(theta)[74];
      REAL(ret)[142] = REAL(theta)[74];
    }
    else if (theta_n == 65){
      REAL(ret)[118] = exp(REAL(theta)[54]);
      REAL(ret)[129] = exp(REAL(theta)[54]);
      REAL(ret)[130] = 2 * REAL(theta)[64];
      REAL(ret)[131] = REAL(theta)[75];
      REAL(ret)[142] = REAL(theta)[75];
    }
    else if (theta_n == 66){
      REAL(ret)[130] = 2 * exp(2 * REAL(theta)[65]);
      REAL(ret)[131] = REAL(theta)[76] * exp(REAL(theta)[65]);
      REAL(ret)[142] = REAL(theta)[76] * exp(REAL(theta)[65]);
    }
    else if (theta_n == 67){
      REAL(ret)[11] = exp(REAL(theta)[0]);
      REAL(ret)[23] = REAL(theta)[1];
      REAL(ret)[35] = REAL(theta)[3];
      REAL(ret)[47] = REAL(theta)[6];
      REAL(ret)[59] = REAL(theta)[10];
      REAL(ret)[71] = REAL(theta)[15];
      REAL(ret)[83] = REAL(theta)[21];
      REAL(ret)[95] = REAL(theta)[28];
      REAL(ret)[107] = REAL(theta)[36];
      REAL(ret)[119] = REAL(theta)[45];
      REAL(ret)[131] = REAL(theta)[55];
      REAL(ret)[132] = exp(REAL(theta)[0]);
      REAL(ret)[133] = REAL(theta)[1];
      REAL(ret)[134] = REAL(theta)[3];
      REAL(ret)[135] = REAL(theta)[6];
      REAL(ret)[136] = REAL(theta)[10];
      REAL(ret)[137] = REAL(theta)[15];
      REAL(ret)[138] = REAL(theta)[21];
      REAL(ret)[139] = REAL(theta)[28];
      REAL(ret)[140] = REAL(theta)[36];
      REAL(ret)[141] = REAL(theta)[45];
      REAL(ret)[142] = REAL(theta)[55];
      REAL(ret)[143] = 2 * REAL(theta)[66];
    }
    else if (theta_n == 68){
      REAL(ret)[23] = exp(REAL(theta)[2]);
      REAL(ret)[35] = REAL(theta)[4];
      REAL(ret)[47] = REAL(theta)[7];
      REAL(ret)[59] = REAL(theta)[11];
      REAL(ret)[71] = REAL(theta)[16];
      REAL(ret)[83] = REAL(theta)[22];
      REAL(ret)[95] = REAL(theta)[29];
      REAL(ret)[107] = REAL(theta)[37];
      REAL(ret)[119] = REAL(theta)[46];
      REAL(ret)[131] = REAL(theta)[56];
      REAL(ret)[133] = exp(REAL(theta)[2]);
      REAL(ret)[134] = REAL(theta)[4];
      REAL(ret)[135] = REAL(theta)[7];
      REAL(ret)[136] = REAL(theta)[11];
      REAL(ret)[137] = REAL(theta)[16];
      REAL(ret)[138] = REAL(theta)[22];
      REAL(ret)[139] = REAL(theta)[29];
      REAL(ret)[140] = REAL(theta)[37];
      REAL(ret)[141] = REAL(theta)[46];
      REAL(ret)[142] = REAL(theta)[56];
      REAL(ret)[143] = 2 * REAL(theta)[67];
    }
    else if (theta_n == 69){
      REAL(ret)[35] = exp(REAL(theta)[5]);
      REAL(ret)[47] = REAL(theta)[8];
      REAL(ret)[59] = REAL(theta)[12];
      REAL(ret)[71] = REAL(theta)[17];
      REAL(ret)[83] = REAL(theta)[23];
      REAL(ret)[95] = REAL(theta)[30];
      REAL(ret)[107] = REAL(theta)[38];
      REAL(ret)[119] = REAL(theta)[47];
      REAL(ret)[131] = REAL(theta)[57];
      REAL(ret)[134] = exp(REAL(theta)[5]);
      REAL(ret)[135] = REAL(theta)[8];
      REAL(ret)[136] = REAL(theta)[12];
      REAL(ret)[137] = REAL(theta)[17];
      REAL(ret)[138] = REAL(theta)[23];
      REAL(ret)[139] = REAL(theta)[30];
      REAL(ret)[140] = REAL(theta)[38];
      REAL(ret)[141] = REAL(theta)[47];
      REAL(ret)[142] = REAL(theta)[57];
      REAL(ret)[143] = 2 * REAL(theta)[68];
    }
    else if (theta_n == 70){
      REAL(ret)[47] = exp(REAL(theta)[9]);
      REAL(ret)[59] = REAL(theta)[13];
      REAL(ret)[71] = REAL(theta)[18];
      REAL(ret)[83] = REAL(theta)[24];
      REAL(ret)[95] = REAL(theta)[31];
      REAL(ret)[107] = REAL(theta)[39];
      REAL(ret)[119] = REAL(theta)[48];
      REAL(ret)[131] = REAL(theta)[58];
      REAL(ret)[135] = exp(REAL(theta)[9]);
      REAL(ret)[136] = REAL(theta)[13];
      REAL(ret)[137] = REAL(theta)[18];
      REAL(ret)[138] = REAL(theta)[24];
      REAL(ret)[139] = REAL(theta)[31];
      REAL(ret)[140] = REAL(theta)[39];
      REAL(ret)[141] = REAL(theta)[48];
      REAL(ret)[142] = REAL(theta)[58];
      REAL(ret)[143] = 2 * REAL(theta)[69];
    }
    else if (theta_n == 71){
      REAL(ret)[59] = exp(REAL(theta)[14]);
      REAL(ret)[71] = REAL(theta)[19];
      REAL(ret)[83] = REAL(theta)[25];
      REAL(ret)[95] = REAL(theta)[32];
      REAL(ret)[107] = REAL(theta)[40];
      REAL(ret)[119] = REAL(theta)[49];
      REAL(ret)[131] = REAL(theta)[59];
      REAL(ret)[136] = exp(REAL(theta)[14]);
      REAL(ret)[137] = REAL(theta)[19];
      REAL(ret)[138] = REAL(theta)[25];
      REAL(ret)[139] = REAL(theta)[32];
      REAL(ret)[140] = REAL(theta)[40];
      REAL(ret)[141] = REAL(theta)[49];
      REAL(ret)[142] = REAL(theta)[59];
      REAL(ret)[143] = 2 * REAL(theta)[70];
    }
    else if (theta_n == 72){
      REAL(ret)[71] = exp(REAL(theta)[20]);
      REAL(ret)[83] = REAL(theta)[26];
      REAL(ret)[95] = REAL(theta)[33];
      REAL(ret)[107] = REAL(theta)[41];
      REAL(ret)[119] = REAL(theta)[50];
      REAL(ret)[131] = REAL(theta)[60];
      REAL(ret)[137] = exp(REAL(theta)[20]);
      REAL(ret)[138] = REAL(theta)[26];
      REAL(ret)[139] = REAL(theta)[33];
      REAL(ret)[140] = REAL(theta)[41];
      REAL(ret)[141] = REAL(theta)[50];
      REAL(ret)[142] = REAL(theta)[60];
      REAL(ret)[143] = 2 * REAL(theta)[71];
    }
    else if (theta_n == 73){
      REAL(ret)[83] = exp(REAL(theta)[27]);
      REAL(ret)[95] = REAL(theta)[34];
      REAL(ret)[107] = REAL(theta)[42];
      REAL(ret)[119] = REAL(theta)[51];
      REAL(ret)[131] = REAL(theta)[61];
      REAL(ret)[138] = exp(REAL(theta)[27]);
      REAL(ret)[139] = REAL(theta)[34];
      REAL(ret)[140] = REAL(theta)[42];
      REAL(ret)[141] = REAL(theta)[51];
      REAL(ret)[142] = REAL(theta)[61];
      REAL(ret)[143] = 2 * REAL(theta)[72];
    }
    else if (theta_n == 74){
      REAL(ret)[95] = exp(REAL(theta)[35]);
      REAL(ret)[107] = REAL(theta)[43];
      REAL(ret)[119] = REAL(theta)[52];
      REAL(ret)[131] = REAL(theta)[62];
      REAL(ret)[139] = exp(REAL(theta)[35]);
      REAL(ret)[140] = REAL(theta)[43];
      REAL(ret)[141] = REAL(theta)[52];
      REAL(ret)[142] = REAL(theta)[62];
      REAL(ret)[143] = 2 * REAL(theta)[73];
    }
    else if (theta_n == 75){
      REAL(ret)[107] = exp(REAL(theta)[44]);
      REAL(ret)[119] = REAL(theta)[53];
      REAL(ret)[131] = REAL(theta)[63];
      REAL(ret)[140] = exp(REAL(theta)[44]);
      REAL(ret)[141] = REAL(theta)[53];
      REAL(ret)[142] = REAL(theta)[63];
      REAL(ret)[143] = 2 * REAL(theta)[74];
    }
    else if (theta_n == 76){
      REAL(ret)[119] = exp(REAL(theta)[54]);
      REAL(ret)[131] = REAL(theta)[64];
      REAL(ret)[141] = exp(REAL(theta)[54]);
      REAL(ret)[142] = REAL(theta)[64];
      REAL(ret)[143] = 2 * REAL(theta)[75];
    }
    else if (theta_n == 77){
      REAL(ret)[131] = exp(REAL(theta)[65]);
      REAL(ret)[142] = exp(REAL(theta)[65]);
      REAL(ret)[143] = 2 * REAL(theta)[76];
    }
    else if (theta_n == 78){
      REAL(ret)[143] = 2 * exp(2 * REAL(theta)[77]);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 12));for(int i = 0; i < 12; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 2 * exp(2 * REAL(theta)[0]);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 2 * exp(2 * REAL(theta)[2]);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 2 * exp(2 * REAL(theta)[5]);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 2 * exp(2 * REAL(theta)[9]);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 2 * exp(2 * REAL(theta)[14]);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 2 * exp(2 * REAL(theta)[20]);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 2 * exp(2 * REAL(theta)[27]);
    }
    else if (theta_n == -38){
      REAL(ret)[7] = 2 * exp(2 * REAL(theta)[35]);
    }
    else if (theta_n == -47){
      REAL(ret)[8] = 2 * exp(2 * REAL(theta)[44]);
    }
    else if (theta_n == -57){
      REAL(ret)[9] = 2 * exp(2 * REAL(theta)[54]);
    }
    else if (theta_n == -68){
      REAL(ret)[10] = 2 * exp(2 * REAL(theta)[65]);
    }
    else if (theta_n == -80){
      REAL(ret)[11] = 2 * exp(2 * REAL(theta)[77]);
    }
    UNPROTECT(1);
    return(ret);
  }
}

  return R_NilValue;
}
