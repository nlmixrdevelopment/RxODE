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
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,1));
    INTEGER(ret)[0]=2;
    UNPROTECT(1);
    return(ret);  
}

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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 1));for(int i = 0; i < 1; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 2){
  int theta_n = INTEGER(tn)[0];
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,3));
    INTEGER(ret)[0]=2;
    INTEGER(ret)[1]=5;
    INTEGER(ret)[2]=2;
    UNPROTECT(1);
    return(ret);  
}

if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 3;
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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[2] = REAL(theta)[1];
      REAL(ret)[3] = R_pow_di(REAL(theta)[2], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[3] = R_pow_di(REAL(theta)[1], 2) + R_pow_di(REAL(theta)[2], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
      REAL(ret)[1] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[2] = 2 * REAL(theta)[0] * REAL(theta)[1];
    }
    else if (theta_n == 2){
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[3] = 2 * REAL(theta)[1];
    }
    else if (theta_n == 3){
      REAL(ret)[3] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 2));for(int i = 0; i < 2; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 3){
  int theta_n = INTEGER(tn)[0];
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,6));
    INTEGER(ret)[0]=2;
    INTEGER(ret)[1]=5;
    INTEGER(ret)[2]=2;
    INTEGER(ret)[3]=5;
    INTEGER(ret)[4]=5;
    INTEGER(ret)[5]=2;
    UNPROTECT(1);
    return(ret);  
}

if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 6;
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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[3] = REAL(theta)[1];
      REAL(ret)[4] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[6] = REAL(theta)[3];
      REAL(ret)[7] = REAL(theta)[4];
      REAL(ret)[8] = R_pow_di(REAL(theta)[5], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[4] = R_pow_di(REAL(theta)[1], 2) + R_pow_di(REAL(theta)[2], 4);
      REAL(ret)[5] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[7] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[8] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + R_pow_di(REAL(theta)[5], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
      REAL(ret)[1] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[2] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[3] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[6] = 2 * REAL(theta)[0] * REAL(theta)[3];
    }
    else if (theta_n == 2){
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[4] = 2 * REAL(theta)[1];
      REAL(ret)[5] = REAL(theta)[3];
      REAL(ret)[7] = REAL(theta)[3];
    }
    else if (theta_n == 3){
      REAL(ret)[4] = 4 * R_pow_di(REAL(theta)[2], 3);
      REAL(ret)[5] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[7] = 2 * REAL(theta)[2] * REAL(theta)[4];
    }
    else if (theta_n == 4){
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[5] = REAL(theta)[1];
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[7] = REAL(theta)[1];
      REAL(ret)[8] = 2 * REAL(theta)[3];
    }
    else if (theta_n == 5){
      REAL(ret)[5] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[7] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[8] = 2 * REAL(theta)[4];
    }
    else if (theta_n == 6){
      REAL(ret)[8] = 4 * R_pow_di(REAL(theta)[5], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 3));for(int i = 0; i < 3; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 4 * R_pow_di(REAL(theta)[5], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 4){
  int theta_n = INTEGER(tn)[0];
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,10));
    INTEGER(ret)[0]=2;
    INTEGER(ret)[1]=5;
    INTEGER(ret)[2]=2;
    INTEGER(ret)[3]=5;
    INTEGER(ret)[4]=5;
    INTEGER(ret)[5]=2;
    INTEGER(ret)[6]=5;
    INTEGER(ret)[7]=5;
    INTEGER(ret)[8]=5;
    INTEGER(ret)[9]=2;
    UNPROTECT(1);
    return(ret);  
}

if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 10;
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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[4] = REAL(theta)[1];
      REAL(ret)[5] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[8] = REAL(theta)[3];
      REAL(ret)[9] = REAL(theta)[4];
      REAL(ret)[10] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[12] = REAL(theta)[6];
      REAL(ret)[13] = REAL(theta)[7];
      REAL(ret)[14] = REAL(theta)[8];
      REAL(ret)[15] = R_pow_di(REAL(theta)[9], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[5] = R_pow_di(REAL(theta)[1], 2) + R_pow_di(REAL(theta)[2], 4);
      REAL(ret)[6] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[7] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[9] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[10] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + R_pow_di(REAL(theta)[5], 4);
      REAL(ret)[11] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[12] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[14] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[15] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + R_pow_di(REAL(theta)[9], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
      REAL(ret)[1] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[2] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[3] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[4] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[8] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[12] = 2 * REAL(theta)[0] * REAL(theta)[6];
    }
    else if (theta_n == 2){
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[5] = 2 * REAL(theta)[1];
      REAL(ret)[6] = REAL(theta)[3];
      REAL(ret)[7] = REAL(theta)[6];
      REAL(ret)[9] = REAL(theta)[3];
      REAL(ret)[13] = REAL(theta)[6];
    }
    else if (theta_n == 3){
      REAL(ret)[5] = 4 * R_pow_di(REAL(theta)[2], 3);
      REAL(ret)[6] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[7] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[9] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[13] = 2 * REAL(theta)[2] * REAL(theta)[7];
    }
    else if (theta_n == 4){
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[6] = REAL(theta)[1];
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[9] = REAL(theta)[1];
      REAL(ret)[10] = 2 * REAL(theta)[3];
      REAL(ret)[11] = REAL(theta)[6];
      REAL(ret)[14] = REAL(theta)[6];
    }
    else if (theta_n == 5){
      REAL(ret)[6] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[9] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[10] = 2 * REAL(theta)[4];
      REAL(ret)[11] = REAL(theta)[7];
      REAL(ret)[14] = REAL(theta)[7];
    }
    else if (theta_n == 6){
      REAL(ret)[10] = 4 * R_pow_di(REAL(theta)[5], 3);
      REAL(ret)[11] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[14] = 2 * REAL(theta)[5] * REAL(theta)[8];
    }
    else if (theta_n == 7){
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[7] = REAL(theta)[1];
      REAL(ret)[11] = REAL(theta)[3];
      REAL(ret)[12] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[14] = REAL(theta)[3];
      REAL(ret)[15] = 2 * REAL(theta)[6];
    }
    else if (theta_n == 8){
      REAL(ret)[7] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[11] = REAL(theta)[4];
      REAL(ret)[13] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[14] = REAL(theta)[4];
      REAL(ret)[15] = 2 * REAL(theta)[7];
    }
    else if (theta_n == 9){
      REAL(ret)[11] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[14] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[15] = 2 * REAL(theta)[8];
    }
    else if (theta_n == 10){
      REAL(ret)[15] = 4 * R_pow_di(REAL(theta)[9], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 4));for(int i = 0; i < 4; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 4 * R_pow_di(REAL(theta)[5], 3);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 4 * R_pow_di(REAL(theta)[9], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 5){
  int theta_n = INTEGER(tn)[0];
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,15));
    INTEGER(ret)[0]=2;
    INTEGER(ret)[1]=5;
    INTEGER(ret)[2]=2;
    INTEGER(ret)[3]=5;
    INTEGER(ret)[4]=5;
    INTEGER(ret)[5]=2;
    INTEGER(ret)[6]=5;
    INTEGER(ret)[7]=5;
    INTEGER(ret)[8]=5;
    INTEGER(ret)[9]=2;
    INTEGER(ret)[10]=5;
    INTEGER(ret)[11]=5;
    INTEGER(ret)[12]=5;
    INTEGER(ret)[13]=5;
    INTEGER(ret)[14]=2;
    UNPROTECT(1);
    return(ret);  
}

if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 15;
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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[5] = REAL(theta)[1];
      REAL(ret)[6] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[10] = REAL(theta)[3];
      REAL(ret)[11] = REAL(theta)[4];
      REAL(ret)[12] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[15] = REAL(theta)[6];
      REAL(ret)[16] = REAL(theta)[7];
      REAL(ret)[17] = REAL(theta)[8];
      REAL(ret)[18] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[20] = REAL(theta)[10];
      REAL(ret)[21] = REAL(theta)[11];
      REAL(ret)[22] = REAL(theta)[12];
      REAL(ret)[23] = REAL(theta)[13];
      REAL(ret)[24] = R_pow_di(REAL(theta)[14], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[6] = R_pow_di(REAL(theta)[1], 2) + R_pow_di(REAL(theta)[2], 4);
      REAL(ret)[7] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[8] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[9] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[10] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[11] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[12] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + R_pow_di(REAL(theta)[5], 4);
      REAL(ret)[13] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[14] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[15] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[16] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[17] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[18] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + R_pow_di(REAL(theta)[9], 4);
      REAL(ret)[19] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[20] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[21] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[22] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[23] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[24] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + R_pow_di(REAL(theta)[14], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
      REAL(ret)[1] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[2] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[3] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[4] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[5] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[10] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[15] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[20] = 2 * REAL(theta)[0] * REAL(theta)[10];
    }
    else if (theta_n == 2){
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[6] = 2 * REAL(theta)[1];
      REAL(ret)[7] = REAL(theta)[3];
      REAL(ret)[8] = REAL(theta)[6];
      REAL(ret)[9] = REAL(theta)[10];
      REAL(ret)[11] = REAL(theta)[3];
      REAL(ret)[16] = REAL(theta)[6];
      REAL(ret)[21] = REAL(theta)[10];
    }
    else if (theta_n == 3){
      REAL(ret)[6] = 4 * R_pow_di(REAL(theta)[2], 3);
      REAL(ret)[7] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[8] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[9] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[11] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[16] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[21] = 2 * REAL(theta)[11] * REAL(theta)[2];
    }
    else if (theta_n == 4){
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[7] = REAL(theta)[1];
      REAL(ret)[10] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[12] = 2 * REAL(theta)[3];
      REAL(ret)[13] = REAL(theta)[6];
      REAL(ret)[14] = REAL(theta)[10];
      REAL(ret)[17] = REAL(theta)[6];
      REAL(ret)[22] = REAL(theta)[10];
    }
    else if (theta_n == 5){
      REAL(ret)[7] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[11] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[12] = 2 * REAL(theta)[4];
      REAL(ret)[13] = REAL(theta)[7];
      REAL(ret)[14] = REAL(theta)[11];
      REAL(ret)[17] = REAL(theta)[7];
      REAL(ret)[22] = REAL(theta)[11];
    }
    else if (theta_n == 6){
      REAL(ret)[12] = 4 * R_pow_di(REAL(theta)[5], 3);
      REAL(ret)[13] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[14] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[17] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[22] = 2 * REAL(theta)[12] * REAL(theta)[5];
    }
    else if (theta_n == 7){
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[8] = REAL(theta)[1];
      REAL(ret)[13] = REAL(theta)[3];
      REAL(ret)[15] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[16] = REAL(theta)[1];
      REAL(ret)[17] = REAL(theta)[3];
      REAL(ret)[18] = 2 * REAL(theta)[6];
      REAL(ret)[19] = REAL(theta)[10];
      REAL(ret)[23] = REAL(theta)[10];
    }
    else if (theta_n == 8){
      REAL(ret)[8] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[13] = REAL(theta)[4];
      REAL(ret)[16] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[17] = REAL(theta)[4];
      REAL(ret)[18] = 2 * REAL(theta)[7];
      REAL(ret)[19] = REAL(theta)[11];
      REAL(ret)[23] = REAL(theta)[11];
    }
    else if (theta_n == 9){
      REAL(ret)[13] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[17] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[18] = 2 * REAL(theta)[8];
      REAL(ret)[19] = REAL(theta)[12];
      REAL(ret)[23] = REAL(theta)[12];
    }
    else if (theta_n == 10){
      REAL(ret)[18] = 4 * R_pow_di(REAL(theta)[9], 3);
      REAL(ret)[19] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[23] = 2 * REAL(theta)[13] * REAL(theta)[9];
    }
    else if (theta_n == 11){
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[9] = REAL(theta)[1];
      REAL(ret)[14] = REAL(theta)[3];
      REAL(ret)[19] = REAL(theta)[6];
      REAL(ret)[20] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[21] = REAL(theta)[1];
      REAL(ret)[22] = REAL(theta)[3];
      REAL(ret)[23] = REAL(theta)[6];
      REAL(ret)[24] = 2 * REAL(theta)[10];
    }
    else if (theta_n == 12){
      REAL(ret)[9] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[14] = REAL(theta)[4];
      REAL(ret)[19] = REAL(theta)[7];
      REAL(ret)[21] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[22] = REAL(theta)[4];
      REAL(ret)[23] = REAL(theta)[7];
      REAL(ret)[24] = 2 * REAL(theta)[11];
    }
    else if (theta_n == 13){
      REAL(ret)[14] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[19] = REAL(theta)[8];
      REAL(ret)[22] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[23] = REAL(theta)[8];
      REAL(ret)[24] = 2 * REAL(theta)[12];
    }
    else if (theta_n == 14){
      REAL(ret)[19] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[23] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[24] = 2 * REAL(theta)[13];
    }
    else if (theta_n == 15){
      REAL(ret)[24] = 4 * R_pow_di(REAL(theta)[14], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 5));for(int i = 0; i < 5; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 4 * R_pow_di(REAL(theta)[5], 3);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 4 * R_pow_di(REAL(theta)[9], 3);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 4 * R_pow_di(REAL(theta)[14], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 6){
  int theta_n = INTEGER(tn)[0];
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,21));
    INTEGER(ret)[0]=2;
    INTEGER(ret)[1]=5;
    INTEGER(ret)[2]=2;
    INTEGER(ret)[3]=5;
    INTEGER(ret)[4]=5;
    INTEGER(ret)[5]=2;
    INTEGER(ret)[6]=5;
    INTEGER(ret)[7]=5;
    INTEGER(ret)[8]=5;
    INTEGER(ret)[9]=2;
    INTEGER(ret)[10]=5;
    INTEGER(ret)[11]=5;
    INTEGER(ret)[12]=5;
    INTEGER(ret)[13]=5;
    INTEGER(ret)[14]=2;
    INTEGER(ret)[15]=5;
    INTEGER(ret)[16]=5;
    INTEGER(ret)[17]=5;
    INTEGER(ret)[18]=5;
    INTEGER(ret)[19]=5;
    INTEGER(ret)[20]=2;
    UNPROTECT(1);
    return(ret);  
}

if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 21;
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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[6] = REAL(theta)[1];
      REAL(ret)[7] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[12] = REAL(theta)[3];
      REAL(ret)[13] = REAL(theta)[4];
      REAL(ret)[14] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[18] = REAL(theta)[6];
      REAL(ret)[19] = REAL(theta)[7];
      REAL(ret)[20] = REAL(theta)[8];
      REAL(ret)[21] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[24] = REAL(theta)[10];
      REAL(ret)[25] = REAL(theta)[11];
      REAL(ret)[26] = REAL(theta)[12];
      REAL(ret)[27] = REAL(theta)[13];
      REAL(ret)[28] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[30] = REAL(theta)[15];
      REAL(ret)[31] = REAL(theta)[16];
      REAL(ret)[32] = REAL(theta)[17];
      REAL(ret)[33] = REAL(theta)[18];
      REAL(ret)[34] = REAL(theta)[19];
      REAL(ret)[35] = R_pow_di(REAL(theta)[20], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[7] = R_pow_di(REAL(theta)[1], 2) + R_pow_di(REAL(theta)[2], 4);
      REAL(ret)[8] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[9] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[10] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[11] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[12] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[14] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + R_pow_di(REAL(theta)[5], 4);
      REAL(ret)[15] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[16] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[17] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[18] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[19] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[20] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[21] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + R_pow_di(REAL(theta)[9], 4);
      REAL(ret)[22] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[23] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[24] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[25] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[26] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[27] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[28] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + R_pow_di(REAL(theta)[14], 4);
      REAL(ret)[29] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[30] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[31] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[32] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[33] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[34] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[35] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + R_pow_di(REAL(theta)[20], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
      REAL(ret)[1] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[2] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[3] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[4] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[5] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[6] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[12] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[18] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[24] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[30] = 2 * REAL(theta)[0] * REAL(theta)[15];
    }
    else if (theta_n == 2){
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[7] = 4 * R_pow_di(REAL(theta)[2], 3);
      REAL(ret)[8] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[9] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[10] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[11] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[13] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[19] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[25] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[31] = 2 * REAL(theta)[16] * REAL(theta)[2];
    }
    else if (theta_n == 4){
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[8] = REAL(theta)[1];
      REAL(ret)[12] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[8] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[13] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[14] = 2 * REAL(theta)[4];
      REAL(ret)[15] = REAL(theta)[7];
      REAL(ret)[16] = REAL(theta)[11];
      REAL(ret)[17] = REAL(theta)[16];
      REAL(ret)[20] = REAL(theta)[7];
      REAL(ret)[26] = REAL(theta)[11];
      REAL(ret)[32] = REAL(theta)[16];
    }
    else if (theta_n == 6){
      REAL(ret)[14] = 4 * R_pow_di(REAL(theta)[5], 3);
      REAL(ret)[15] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[16] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[17] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[20] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[26] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[32] = 2 * REAL(theta)[17] * REAL(theta)[5];
    }
    else if (theta_n == 7){
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[9] = REAL(theta)[1];
      REAL(ret)[15] = REAL(theta)[3];
      REAL(ret)[18] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[19] = REAL(theta)[1];
      REAL(ret)[20] = REAL(theta)[3];
      REAL(ret)[21] = 2 * REAL(theta)[6];
      REAL(ret)[22] = REAL(theta)[10];
      REAL(ret)[23] = REAL(theta)[15];
      REAL(ret)[27] = REAL(theta)[10];
      REAL(ret)[33] = REAL(theta)[15];
    }
    else if (theta_n == 8){
      REAL(ret)[9] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[15] = REAL(theta)[4];
      REAL(ret)[19] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[20] = REAL(theta)[4];
      REAL(ret)[21] = 2 * REAL(theta)[7];
      REAL(ret)[22] = REAL(theta)[11];
      REAL(ret)[23] = REAL(theta)[16];
      REAL(ret)[27] = REAL(theta)[11];
      REAL(ret)[33] = REAL(theta)[16];
    }
    else if (theta_n == 9){
      REAL(ret)[15] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[20] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[21] = 2 * REAL(theta)[8];
      REAL(ret)[22] = REAL(theta)[12];
      REAL(ret)[23] = REAL(theta)[17];
      REAL(ret)[27] = REAL(theta)[12];
      REAL(ret)[33] = REAL(theta)[17];
    }
    else if (theta_n == 10){
      REAL(ret)[21] = 4 * R_pow_di(REAL(theta)[9], 3);
      REAL(ret)[22] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[23] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[27] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[33] = 2 * REAL(theta)[18] * REAL(theta)[9];
    }
    else if (theta_n == 11){
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[10] = REAL(theta)[1];
      REAL(ret)[16] = REAL(theta)[3];
      REAL(ret)[22] = REAL(theta)[6];
      REAL(ret)[24] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[25] = REAL(theta)[1];
      REAL(ret)[26] = REAL(theta)[3];
      REAL(ret)[27] = REAL(theta)[6];
      REAL(ret)[28] = 2 * REAL(theta)[10];
      REAL(ret)[29] = REAL(theta)[15];
      REAL(ret)[34] = REAL(theta)[15];
    }
    else if (theta_n == 12){
      REAL(ret)[10] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[16] = REAL(theta)[4];
      REAL(ret)[22] = REAL(theta)[7];
      REAL(ret)[25] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[26] = REAL(theta)[4];
      REAL(ret)[27] = REAL(theta)[7];
      REAL(ret)[28] = 2 * REAL(theta)[11];
      REAL(ret)[29] = REAL(theta)[16];
      REAL(ret)[34] = REAL(theta)[16];
    }
    else if (theta_n == 13){
      REAL(ret)[16] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[22] = REAL(theta)[8];
      REAL(ret)[26] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[27] = REAL(theta)[8];
      REAL(ret)[28] = 2 * REAL(theta)[12];
      REAL(ret)[29] = REAL(theta)[17];
      REAL(ret)[34] = REAL(theta)[17];
    }
    else if (theta_n == 14){
      REAL(ret)[22] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[27] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[28] = 2 * REAL(theta)[13];
      REAL(ret)[29] = REAL(theta)[18];
      REAL(ret)[34] = REAL(theta)[18];
    }
    else if (theta_n == 15){
      REAL(ret)[28] = 4 * R_pow_di(REAL(theta)[14], 3);
      REAL(ret)[29] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[34] = 2 * REAL(theta)[14] * REAL(theta)[19];
    }
    else if (theta_n == 16){
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[17] = REAL(theta)[3];
      REAL(ret)[23] = REAL(theta)[6];
      REAL(ret)[29] = REAL(theta)[10];
      REAL(ret)[30] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[31] = REAL(theta)[1];
      REAL(ret)[32] = REAL(theta)[3];
      REAL(ret)[33] = REAL(theta)[6];
      REAL(ret)[34] = REAL(theta)[10];
      REAL(ret)[35] = 2 * REAL(theta)[15];
    }
    else if (theta_n == 17){
      REAL(ret)[11] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[17] = REAL(theta)[4];
      REAL(ret)[23] = REAL(theta)[7];
      REAL(ret)[29] = REAL(theta)[11];
      REAL(ret)[31] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[32] = REAL(theta)[4];
      REAL(ret)[33] = REAL(theta)[7];
      REAL(ret)[34] = REAL(theta)[11];
      REAL(ret)[35] = 2 * REAL(theta)[16];
    }
    else if (theta_n == 18){
      REAL(ret)[17] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[23] = REAL(theta)[8];
      REAL(ret)[29] = REAL(theta)[12];
      REAL(ret)[32] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[33] = REAL(theta)[8];
      REAL(ret)[34] = REAL(theta)[12];
      REAL(ret)[35] = 2 * REAL(theta)[17];
    }
    else if (theta_n == 19){
      REAL(ret)[23] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[29] = REAL(theta)[13];
      REAL(ret)[33] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[34] = REAL(theta)[13];
      REAL(ret)[35] = 2 * REAL(theta)[18];
    }
    else if (theta_n == 20){
      REAL(ret)[29] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[34] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[35] = 2 * REAL(theta)[19];
    }
    else if (theta_n == 21){
      REAL(ret)[35] = 4 * R_pow_di(REAL(theta)[20], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 6));for(int i = 0; i < 6; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 4 * R_pow_di(REAL(theta)[5], 3);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 4 * R_pow_di(REAL(theta)[9], 3);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 4 * R_pow_di(REAL(theta)[14], 3);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 4 * R_pow_di(REAL(theta)[20], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 7){
  int theta_n = INTEGER(tn)[0];
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,28));
    INTEGER(ret)[0]=2;
    INTEGER(ret)[1]=5;
    INTEGER(ret)[2]=2;
    INTEGER(ret)[3]=5;
    INTEGER(ret)[4]=5;
    INTEGER(ret)[5]=2;
    INTEGER(ret)[6]=5;
    INTEGER(ret)[7]=5;
    INTEGER(ret)[8]=5;
    INTEGER(ret)[9]=2;
    INTEGER(ret)[10]=5;
    INTEGER(ret)[11]=5;
    INTEGER(ret)[12]=5;
    INTEGER(ret)[13]=5;
    INTEGER(ret)[14]=2;
    INTEGER(ret)[15]=5;
    INTEGER(ret)[16]=5;
    INTEGER(ret)[17]=5;
    INTEGER(ret)[18]=5;
    INTEGER(ret)[19]=5;
    INTEGER(ret)[20]=2;
    INTEGER(ret)[21]=5;
    INTEGER(ret)[22]=5;
    INTEGER(ret)[23]=5;
    INTEGER(ret)[24]=5;
    INTEGER(ret)[25]=5;
    INTEGER(ret)[26]=5;
    INTEGER(ret)[27]=2;
    UNPROTECT(1);
    return(ret);  
}

if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 28;
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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[7] = REAL(theta)[1];
      REAL(ret)[8] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[14] = REAL(theta)[3];
      REAL(ret)[15] = REAL(theta)[4];
      REAL(ret)[16] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[21] = REAL(theta)[6];
      REAL(ret)[22] = REAL(theta)[7];
      REAL(ret)[23] = REAL(theta)[8];
      REAL(ret)[24] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[28] = REAL(theta)[10];
      REAL(ret)[29] = REAL(theta)[11];
      REAL(ret)[30] = REAL(theta)[12];
      REAL(ret)[31] = REAL(theta)[13];
      REAL(ret)[32] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[35] = REAL(theta)[15];
      REAL(ret)[36] = REAL(theta)[16];
      REAL(ret)[37] = REAL(theta)[17];
      REAL(ret)[38] = REAL(theta)[18];
      REAL(ret)[39] = REAL(theta)[19];
      REAL(ret)[40] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[42] = REAL(theta)[21];
      REAL(ret)[43] = REAL(theta)[22];
      REAL(ret)[44] = REAL(theta)[23];
      REAL(ret)[45] = REAL(theta)[24];
      REAL(ret)[46] = REAL(theta)[25];
      REAL(ret)[47] = REAL(theta)[26];
      REAL(ret)[48] = R_pow_di(REAL(theta)[27], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[8] = R_pow_di(REAL(theta)[1], 2) + R_pow_di(REAL(theta)[2], 4);
      REAL(ret)[9] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[10] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[11] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[12] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[14] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[16] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + R_pow_di(REAL(theta)[5], 4);
      REAL(ret)[17] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[18] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[19] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[20] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[21] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[22] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[23] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[24] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + R_pow_di(REAL(theta)[9], 4);
      REAL(ret)[25] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[26] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[27] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[28] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[29] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[30] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[31] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[32] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + R_pow_di(REAL(theta)[14], 4);
      REAL(ret)[33] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[34] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[35] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[36] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[37] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[38] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[39] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[40] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + R_pow_di(REAL(theta)[20], 4);
      REAL(ret)[41] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[42] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[43] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[44] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[45] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[46] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[47] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[48] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + R_pow_di(REAL(theta)[27], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
      REAL(ret)[1] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[2] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[3] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[4] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[5] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[6] = 2 * REAL(theta)[0] * REAL(theta)[21];
      REAL(ret)[7] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[14] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[21] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[28] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[35] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[42] = 2 * REAL(theta)[0] * REAL(theta)[21];
    }
    else if (theta_n == 2){
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[8] = 4 * R_pow_di(REAL(theta)[2], 3);
      REAL(ret)[9] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[10] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[11] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[12] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[13] = 2 * REAL(theta)[2] * REAL(theta)[22];
      REAL(ret)[15] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[22] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[29] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[36] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[43] = 2 * REAL(theta)[2] * REAL(theta)[22];
    }
    else if (theta_n == 4){
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[9] = REAL(theta)[1];
      REAL(ret)[14] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[9] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[15] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[16] = 4 * R_pow_di(REAL(theta)[5], 3);
      REAL(ret)[17] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[18] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[19] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[20] = 2 * REAL(theta)[23] * REAL(theta)[5];
      REAL(ret)[23] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[30] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[37] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[44] = 2 * REAL(theta)[23] * REAL(theta)[5];
    }
    else if (theta_n == 7){
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[10] = REAL(theta)[1];
      REAL(ret)[17] = REAL(theta)[3];
      REAL(ret)[21] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[10] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[17] = REAL(theta)[4];
      REAL(ret)[22] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[17] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[23] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[24] = 2 * REAL(theta)[8];
      REAL(ret)[25] = REAL(theta)[12];
      REAL(ret)[26] = REAL(theta)[17];
      REAL(ret)[27] = REAL(theta)[23];
      REAL(ret)[31] = REAL(theta)[12];
      REAL(ret)[38] = REAL(theta)[17];
      REAL(ret)[45] = REAL(theta)[23];
    }
    else if (theta_n == 10){
      REAL(ret)[24] = 4 * R_pow_di(REAL(theta)[9], 3);
      REAL(ret)[25] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[26] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[27] = 2 * REAL(theta)[24] * REAL(theta)[9];
      REAL(ret)[31] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[38] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[45] = 2 * REAL(theta)[24] * REAL(theta)[9];
    }
    else if (theta_n == 11){
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[18] = REAL(theta)[3];
      REAL(ret)[25] = REAL(theta)[6];
      REAL(ret)[28] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[11] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[18] = REAL(theta)[4];
      REAL(ret)[25] = REAL(theta)[7];
      REAL(ret)[29] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[30] = REAL(theta)[4];
      REAL(ret)[31] = REAL(theta)[7];
      REAL(ret)[32] = 2 * REAL(theta)[11];
      REAL(ret)[33] = REAL(theta)[16];
      REAL(ret)[34] = REAL(theta)[22];
      REAL(ret)[39] = REAL(theta)[16];
      REAL(ret)[46] = REAL(theta)[22];
    }
    else if (theta_n == 13){
      REAL(ret)[18] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[25] = REAL(theta)[8];
      REAL(ret)[30] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[31] = REAL(theta)[8];
      REAL(ret)[32] = 2 * REAL(theta)[12];
      REAL(ret)[33] = REAL(theta)[17];
      REAL(ret)[34] = REAL(theta)[23];
      REAL(ret)[39] = REAL(theta)[17];
      REAL(ret)[46] = REAL(theta)[23];
    }
    else if (theta_n == 14){
      REAL(ret)[25] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[31] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[32] = 2 * REAL(theta)[13];
      REAL(ret)[33] = REAL(theta)[18];
      REAL(ret)[34] = REAL(theta)[24];
      REAL(ret)[39] = REAL(theta)[18];
      REAL(ret)[46] = REAL(theta)[24];
    }
    else if (theta_n == 15){
      REAL(ret)[32] = 4 * R_pow_di(REAL(theta)[14], 3);
      REAL(ret)[33] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[34] = 2 * REAL(theta)[14] * REAL(theta)[25];
      REAL(ret)[39] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[46] = 2 * REAL(theta)[14] * REAL(theta)[25];
    }
    else if (theta_n == 16){
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[12] = REAL(theta)[1];
      REAL(ret)[19] = REAL(theta)[3];
      REAL(ret)[26] = REAL(theta)[6];
      REAL(ret)[33] = REAL(theta)[10];
      REAL(ret)[35] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[36] = REAL(theta)[1];
      REAL(ret)[37] = REAL(theta)[3];
      REAL(ret)[38] = REAL(theta)[6];
      REAL(ret)[39] = REAL(theta)[10];
      REAL(ret)[40] = 2 * REAL(theta)[15];
      REAL(ret)[41] = REAL(theta)[21];
      REAL(ret)[47] = REAL(theta)[21];
    }
    else if (theta_n == 17){
      REAL(ret)[12] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[19] = REAL(theta)[4];
      REAL(ret)[26] = REAL(theta)[7];
      REAL(ret)[33] = REAL(theta)[11];
      REAL(ret)[36] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[37] = REAL(theta)[4];
      REAL(ret)[38] = REAL(theta)[7];
      REAL(ret)[39] = REAL(theta)[11];
      REAL(ret)[40] = 2 * REAL(theta)[16];
      REAL(ret)[41] = REAL(theta)[22];
      REAL(ret)[47] = REAL(theta)[22];
    }
    else if (theta_n == 18){
      REAL(ret)[19] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[26] = REAL(theta)[8];
      REAL(ret)[33] = REAL(theta)[12];
      REAL(ret)[37] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[38] = REAL(theta)[8];
      REAL(ret)[39] = REAL(theta)[12];
      REAL(ret)[40] = 2 * REAL(theta)[17];
      REAL(ret)[41] = REAL(theta)[23];
      REAL(ret)[47] = REAL(theta)[23];
    }
    else if (theta_n == 19){
      REAL(ret)[26] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[33] = REAL(theta)[13];
      REAL(ret)[38] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[39] = REAL(theta)[13];
      REAL(ret)[40] = 2 * REAL(theta)[18];
      REAL(ret)[41] = REAL(theta)[24];
      REAL(ret)[47] = REAL(theta)[24];
    }
    else if (theta_n == 20){
      REAL(ret)[33] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[39] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[40] = 2 * REAL(theta)[19];
      REAL(ret)[41] = REAL(theta)[25];
      REAL(ret)[47] = REAL(theta)[25];
    }
    else if (theta_n == 21){
      REAL(ret)[40] = 4 * R_pow_di(REAL(theta)[20], 3);
      REAL(ret)[41] = 2 * REAL(theta)[20] * REAL(theta)[26];
      REAL(ret)[47] = 2 * REAL(theta)[20] * REAL(theta)[26];
    }
    else if (theta_n == 22){
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[20] = REAL(theta)[3];
      REAL(ret)[27] = REAL(theta)[6];
      REAL(ret)[34] = REAL(theta)[10];
      REAL(ret)[41] = REAL(theta)[15];
      REAL(ret)[42] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[43] = REAL(theta)[1];
      REAL(ret)[44] = REAL(theta)[3];
      REAL(ret)[45] = REAL(theta)[6];
      REAL(ret)[46] = REAL(theta)[10];
      REAL(ret)[47] = REAL(theta)[15];
      REAL(ret)[48] = 2 * REAL(theta)[21];
    }
    else if (theta_n == 23){
      REAL(ret)[13] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[20] = REAL(theta)[4];
      REAL(ret)[27] = REAL(theta)[7];
      REAL(ret)[34] = REAL(theta)[11];
      REAL(ret)[41] = REAL(theta)[16];
      REAL(ret)[43] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[44] = REAL(theta)[4];
      REAL(ret)[45] = REAL(theta)[7];
      REAL(ret)[46] = REAL(theta)[11];
      REAL(ret)[47] = REAL(theta)[16];
      REAL(ret)[48] = 2 * REAL(theta)[22];
    }
    else if (theta_n == 24){
      REAL(ret)[20] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[27] = REAL(theta)[8];
      REAL(ret)[34] = REAL(theta)[12];
      REAL(ret)[41] = REAL(theta)[17];
      REAL(ret)[44] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[45] = REAL(theta)[8];
      REAL(ret)[46] = REAL(theta)[12];
      REAL(ret)[47] = REAL(theta)[17];
      REAL(ret)[48] = 2 * REAL(theta)[23];
    }
    else if (theta_n == 25){
      REAL(ret)[27] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[34] = REAL(theta)[13];
      REAL(ret)[41] = REAL(theta)[18];
      REAL(ret)[45] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[46] = REAL(theta)[13];
      REAL(ret)[47] = REAL(theta)[18];
      REAL(ret)[48] = 2 * REAL(theta)[24];
    }
    else if (theta_n == 26){
      REAL(ret)[34] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[41] = REAL(theta)[19];
      REAL(ret)[46] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[47] = REAL(theta)[19];
      REAL(ret)[48] = 2 * REAL(theta)[25];
    }
    else if (theta_n == 27){
      REAL(ret)[41] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[47] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[48] = 2 * REAL(theta)[26];
    }
    else if (theta_n == 28){
      REAL(ret)[48] = 4 * R_pow_di(REAL(theta)[27], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 7));for(int i = 0; i < 7; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 4 * R_pow_di(REAL(theta)[5], 3);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 4 * R_pow_di(REAL(theta)[9], 3);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 4 * R_pow_di(REAL(theta)[14], 3);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 4 * R_pow_di(REAL(theta)[20], 3);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 4 * R_pow_di(REAL(theta)[27], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 8){
  int theta_n = INTEGER(tn)[0];
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,36));
    INTEGER(ret)[0]=2;
    INTEGER(ret)[1]=5;
    INTEGER(ret)[2]=2;
    INTEGER(ret)[3]=5;
    INTEGER(ret)[4]=5;
    INTEGER(ret)[5]=2;
    INTEGER(ret)[6]=5;
    INTEGER(ret)[7]=5;
    INTEGER(ret)[8]=5;
    INTEGER(ret)[9]=2;
    INTEGER(ret)[10]=5;
    INTEGER(ret)[11]=5;
    INTEGER(ret)[12]=5;
    INTEGER(ret)[13]=5;
    INTEGER(ret)[14]=2;
    INTEGER(ret)[15]=5;
    INTEGER(ret)[16]=5;
    INTEGER(ret)[17]=5;
    INTEGER(ret)[18]=5;
    INTEGER(ret)[19]=5;
    INTEGER(ret)[20]=2;
    INTEGER(ret)[21]=5;
    INTEGER(ret)[22]=5;
    INTEGER(ret)[23]=5;
    INTEGER(ret)[24]=5;
    INTEGER(ret)[25]=5;
    INTEGER(ret)[26]=5;
    INTEGER(ret)[27]=2;
    INTEGER(ret)[28]=5;
    INTEGER(ret)[29]=5;
    INTEGER(ret)[30]=5;
    INTEGER(ret)[31]=5;
    INTEGER(ret)[32]=5;
    INTEGER(ret)[33]=5;
    INTEGER(ret)[34]=5;
    INTEGER(ret)[35]=2;
    UNPROTECT(1);
    return(ret);  
}

if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 36;
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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[8] = REAL(theta)[1];
      REAL(ret)[9] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[16] = REAL(theta)[3];
      REAL(ret)[17] = REAL(theta)[4];
      REAL(ret)[18] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[24] = REAL(theta)[6];
      REAL(ret)[25] = REAL(theta)[7];
      REAL(ret)[26] = REAL(theta)[8];
      REAL(ret)[27] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[32] = REAL(theta)[10];
      REAL(ret)[33] = REAL(theta)[11];
      REAL(ret)[34] = REAL(theta)[12];
      REAL(ret)[35] = REAL(theta)[13];
      REAL(ret)[36] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[40] = REAL(theta)[15];
      REAL(ret)[41] = REAL(theta)[16];
      REAL(ret)[42] = REAL(theta)[17];
      REAL(ret)[43] = REAL(theta)[18];
      REAL(ret)[44] = REAL(theta)[19];
      REAL(ret)[45] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[48] = REAL(theta)[21];
      REAL(ret)[49] = REAL(theta)[22];
      REAL(ret)[50] = REAL(theta)[23];
      REAL(ret)[51] = REAL(theta)[24];
      REAL(ret)[52] = REAL(theta)[25];
      REAL(ret)[53] = REAL(theta)[26];
      REAL(ret)[54] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[56] = REAL(theta)[28];
      REAL(ret)[57] = REAL(theta)[29];
      REAL(ret)[58] = REAL(theta)[30];
      REAL(ret)[59] = REAL(theta)[31];
      REAL(ret)[60] = REAL(theta)[32];
      REAL(ret)[61] = REAL(theta)[33];
      REAL(ret)[62] = REAL(theta)[34];
      REAL(ret)[63] = R_pow_di(REAL(theta)[35], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[28];
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[9] = R_pow_di(REAL(theta)[1], 2) + R_pow_di(REAL(theta)[2], 4);
      REAL(ret)[10] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[11] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[12] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[14] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[28] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[29];
      REAL(ret)[16] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[17] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[18] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + R_pow_di(REAL(theta)[5], 4);
      REAL(ret)[19] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[20] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[21] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[22] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[23] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[24] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[25] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[26] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[27] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + R_pow_di(REAL(theta)[9], 4);
      REAL(ret)[28] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[29] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[30] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[31] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[32] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[33] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[34] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[35] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[36] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + R_pow_di(REAL(theta)[14], 4);
      REAL(ret)[37] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[38] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[39] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[32];
      REAL(ret)[40] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[41] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[42] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[43] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[44] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[45] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + R_pow_di(REAL(theta)[20], 4);
      REAL(ret)[46] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[47] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[33];
      REAL(ret)[48] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[49] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[50] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[51] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[52] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[53] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[54] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + R_pow_di(REAL(theta)[27], 4);
      REAL(ret)[55] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[34];
      REAL(ret)[56] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[28];
      REAL(ret)[57] = REAL(theta)[1] * REAL(theta)[28] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[29];
      REAL(ret)[58] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[59] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[60] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[32];
      REAL(ret)[61] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[33];
      REAL(ret)[62] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[34];
      REAL(ret)[63] = R_pow_di(REAL(theta)[28], 2) + R_pow_di(REAL(theta)[29], 2) + R_pow_di(REAL(theta)[30], 2) + R_pow_di(REAL(theta)[31], 2) + R_pow_di(REAL(theta)[32], 2) + R_pow_di(REAL(theta)[33], 2) + R_pow_di(REAL(theta)[34], 2) + R_pow_di(REAL(theta)[35], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
      REAL(ret)[1] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[2] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[3] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[4] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[5] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[6] = 2 * REAL(theta)[0] * REAL(theta)[21];
      REAL(ret)[7] = 2 * REAL(theta)[0] * REAL(theta)[28];
      REAL(ret)[8] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[16] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[24] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[32] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[40] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[48] = 2 * REAL(theta)[0] * REAL(theta)[21];
      REAL(ret)[56] = 2 * REAL(theta)[0] * REAL(theta)[28];
    }
    else if (theta_n == 2){
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[9] = 4 * R_pow_di(REAL(theta)[2], 3);
      REAL(ret)[10] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[11] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[12] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[13] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[14] = 2 * REAL(theta)[2] * REAL(theta)[22];
      REAL(ret)[15] = 2 * REAL(theta)[2] * REAL(theta)[29];
      REAL(ret)[17] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[25] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[33] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[41] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[49] = 2 * REAL(theta)[2] * REAL(theta)[22];
      REAL(ret)[57] = 2 * REAL(theta)[2] * REAL(theta)[29];
    }
    else if (theta_n == 4){
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[10] = REAL(theta)[1];
      REAL(ret)[16] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[10] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[17] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[18] = 4 * R_pow_di(REAL(theta)[5], 3);
      REAL(ret)[19] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[20] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[21] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[22] = 2 * REAL(theta)[23] * REAL(theta)[5];
      REAL(ret)[23] = 2 * REAL(theta)[30] * REAL(theta)[5];
      REAL(ret)[26] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[34] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[42] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[50] = 2 * REAL(theta)[23] * REAL(theta)[5];
      REAL(ret)[58] = 2 * REAL(theta)[30] * REAL(theta)[5];
    }
    else if (theta_n == 7){
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[19] = REAL(theta)[3];
      REAL(ret)[24] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[11] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[19] = REAL(theta)[4];
      REAL(ret)[25] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[19] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[26] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[27] = 4 * R_pow_di(REAL(theta)[9], 3);
      REAL(ret)[28] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[29] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[30] = 2 * REAL(theta)[24] * REAL(theta)[9];
      REAL(ret)[31] = 2 * REAL(theta)[31] * REAL(theta)[9];
      REAL(ret)[35] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[43] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[51] = 2 * REAL(theta)[24] * REAL(theta)[9];
      REAL(ret)[59] = 2 * REAL(theta)[31] * REAL(theta)[9];
    }
    else if (theta_n == 11){
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[12] = REAL(theta)[1];
      REAL(ret)[20] = REAL(theta)[3];
      REAL(ret)[28] = REAL(theta)[6];
      REAL(ret)[32] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[12] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[20] = REAL(theta)[4];
      REAL(ret)[28] = REAL(theta)[7];
      REAL(ret)[33] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[20] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[28] = REAL(theta)[8];
      REAL(ret)[34] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[28] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[35] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[36] = 2 * REAL(theta)[13];
      REAL(ret)[37] = REAL(theta)[18];
      REAL(ret)[38] = REAL(theta)[24];
      REAL(ret)[39] = REAL(theta)[31];
      REAL(ret)[44] = REAL(theta)[18];
      REAL(ret)[52] = REAL(theta)[24];
      REAL(ret)[60] = REAL(theta)[31];
    }
    else if (theta_n == 15){
      REAL(ret)[36] = 4 * R_pow_di(REAL(theta)[14], 3);
      REAL(ret)[37] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[38] = 2 * REAL(theta)[14] * REAL(theta)[25];
      REAL(ret)[39] = 2 * REAL(theta)[14] * REAL(theta)[32];
      REAL(ret)[44] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[52] = 2 * REAL(theta)[14] * REAL(theta)[25];
      REAL(ret)[60] = 2 * REAL(theta)[14] * REAL(theta)[32];
    }
    else if (theta_n == 16){
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[21] = REAL(theta)[3];
      REAL(ret)[29] = REAL(theta)[6];
      REAL(ret)[37] = REAL(theta)[10];
      REAL(ret)[40] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[13] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[21] = REAL(theta)[4];
      REAL(ret)[29] = REAL(theta)[7];
      REAL(ret)[37] = REAL(theta)[11];
      REAL(ret)[41] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[21] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[29] = REAL(theta)[8];
      REAL(ret)[37] = REAL(theta)[12];
      REAL(ret)[42] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[43] = REAL(theta)[8];
      REAL(ret)[44] = REAL(theta)[12];
      REAL(ret)[45] = 2 * REAL(theta)[17];
      REAL(ret)[46] = REAL(theta)[23];
      REAL(ret)[47] = REAL(theta)[30];
      REAL(ret)[53] = REAL(theta)[23];
      REAL(ret)[61] = REAL(theta)[30];
    }
    else if (theta_n == 19){
      REAL(ret)[29] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[37] = REAL(theta)[13];
      REAL(ret)[43] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[44] = REAL(theta)[13];
      REAL(ret)[45] = 2 * REAL(theta)[18];
      REAL(ret)[46] = REAL(theta)[24];
      REAL(ret)[47] = REAL(theta)[31];
      REAL(ret)[53] = REAL(theta)[24];
      REAL(ret)[61] = REAL(theta)[31];
    }
    else if (theta_n == 20){
      REAL(ret)[37] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[44] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[45] = 2 * REAL(theta)[19];
      REAL(ret)[46] = REAL(theta)[25];
      REAL(ret)[47] = REAL(theta)[32];
      REAL(ret)[53] = REAL(theta)[25];
      REAL(ret)[61] = REAL(theta)[32];
    }
    else if (theta_n == 21){
      REAL(ret)[45] = 4 * R_pow_di(REAL(theta)[20], 3);
      REAL(ret)[46] = 2 * REAL(theta)[20] * REAL(theta)[26];
      REAL(ret)[47] = 2 * REAL(theta)[20] * REAL(theta)[33];
      REAL(ret)[53] = 2 * REAL(theta)[20] * REAL(theta)[26];
      REAL(ret)[61] = 2 * REAL(theta)[20] * REAL(theta)[33];
    }
    else if (theta_n == 22){
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[14] = REAL(theta)[1];
      REAL(ret)[22] = REAL(theta)[3];
      REAL(ret)[30] = REAL(theta)[6];
      REAL(ret)[38] = REAL(theta)[10];
      REAL(ret)[46] = REAL(theta)[15];
      REAL(ret)[48] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[14] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[22] = REAL(theta)[4];
      REAL(ret)[30] = REAL(theta)[7];
      REAL(ret)[38] = REAL(theta)[11];
      REAL(ret)[46] = REAL(theta)[16];
      REAL(ret)[49] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[50] = REAL(theta)[4];
      REAL(ret)[51] = REAL(theta)[7];
      REAL(ret)[52] = REAL(theta)[11];
      REAL(ret)[53] = REAL(theta)[16];
      REAL(ret)[54] = 2 * REAL(theta)[22];
      REAL(ret)[55] = REAL(theta)[29];
      REAL(ret)[62] = REAL(theta)[29];
    }
    else if (theta_n == 24){
      REAL(ret)[22] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[30] = REAL(theta)[8];
      REAL(ret)[38] = REAL(theta)[12];
      REAL(ret)[46] = REAL(theta)[17];
      REAL(ret)[50] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[51] = REAL(theta)[8];
      REAL(ret)[52] = REAL(theta)[12];
      REAL(ret)[53] = REAL(theta)[17];
      REAL(ret)[54] = 2 * REAL(theta)[23];
      REAL(ret)[55] = REAL(theta)[30];
      REAL(ret)[62] = REAL(theta)[30];
    }
    else if (theta_n == 25){
      REAL(ret)[30] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[38] = REAL(theta)[13];
      REAL(ret)[46] = REAL(theta)[18];
      REAL(ret)[51] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[52] = REAL(theta)[13];
      REAL(ret)[53] = REAL(theta)[18];
      REAL(ret)[54] = 2 * REAL(theta)[24];
      REAL(ret)[55] = REAL(theta)[31];
      REAL(ret)[62] = REAL(theta)[31];
    }
    else if (theta_n == 26){
      REAL(ret)[38] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[46] = REAL(theta)[19];
      REAL(ret)[52] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[53] = REAL(theta)[19];
      REAL(ret)[54] = 2 * REAL(theta)[25];
      REAL(ret)[55] = REAL(theta)[32];
      REAL(ret)[62] = REAL(theta)[32];
    }
    else if (theta_n == 27){
      REAL(ret)[46] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[53] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[54] = 2 * REAL(theta)[26];
      REAL(ret)[55] = REAL(theta)[33];
      REAL(ret)[62] = REAL(theta)[33];
    }
    else if (theta_n == 28){
      REAL(ret)[54] = 4 * R_pow_di(REAL(theta)[27], 3);
      REAL(ret)[55] = 2 * REAL(theta)[27] * REAL(theta)[34];
      REAL(ret)[62] = 2 * REAL(theta)[27] * REAL(theta)[34];
    }
    else if (theta_n == 29){
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[15] = REAL(theta)[1];
      REAL(ret)[23] = REAL(theta)[3];
      REAL(ret)[31] = REAL(theta)[6];
      REAL(ret)[39] = REAL(theta)[10];
      REAL(ret)[47] = REAL(theta)[15];
      REAL(ret)[55] = REAL(theta)[21];
      REAL(ret)[56] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[57] = REAL(theta)[1];
      REAL(ret)[58] = REAL(theta)[3];
      REAL(ret)[59] = REAL(theta)[6];
      REAL(ret)[60] = REAL(theta)[10];
      REAL(ret)[61] = REAL(theta)[15];
      REAL(ret)[62] = REAL(theta)[21];
      REAL(ret)[63] = 2 * REAL(theta)[28];
    }
    else if (theta_n == 30){
      REAL(ret)[15] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[23] = REAL(theta)[4];
      REAL(ret)[31] = REAL(theta)[7];
      REAL(ret)[39] = REAL(theta)[11];
      REAL(ret)[47] = REAL(theta)[16];
      REAL(ret)[55] = REAL(theta)[22];
      REAL(ret)[57] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[58] = REAL(theta)[4];
      REAL(ret)[59] = REAL(theta)[7];
      REAL(ret)[60] = REAL(theta)[11];
      REAL(ret)[61] = REAL(theta)[16];
      REAL(ret)[62] = REAL(theta)[22];
      REAL(ret)[63] = 2 * REAL(theta)[29];
    }
    else if (theta_n == 31){
      REAL(ret)[23] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[31] = REAL(theta)[8];
      REAL(ret)[39] = REAL(theta)[12];
      REAL(ret)[47] = REAL(theta)[17];
      REAL(ret)[55] = REAL(theta)[23];
      REAL(ret)[58] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[59] = REAL(theta)[8];
      REAL(ret)[60] = REAL(theta)[12];
      REAL(ret)[61] = REAL(theta)[17];
      REAL(ret)[62] = REAL(theta)[23];
      REAL(ret)[63] = 2 * REAL(theta)[30];
    }
    else if (theta_n == 32){
      REAL(ret)[31] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[39] = REAL(theta)[13];
      REAL(ret)[47] = REAL(theta)[18];
      REAL(ret)[55] = REAL(theta)[24];
      REAL(ret)[59] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[60] = REAL(theta)[13];
      REAL(ret)[61] = REAL(theta)[18];
      REAL(ret)[62] = REAL(theta)[24];
      REAL(ret)[63] = 2 * REAL(theta)[31];
    }
    else if (theta_n == 33){
      REAL(ret)[39] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[47] = REAL(theta)[19];
      REAL(ret)[55] = REAL(theta)[25];
      REAL(ret)[60] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[61] = REAL(theta)[19];
      REAL(ret)[62] = REAL(theta)[25];
      REAL(ret)[63] = 2 * REAL(theta)[32];
    }
    else if (theta_n == 34){
      REAL(ret)[47] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[55] = REAL(theta)[26];
      REAL(ret)[61] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[62] = REAL(theta)[26];
      REAL(ret)[63] = 2 * REAL(theta)[33];
    }
    else if (theta_n == 35){
      REAL(ret)[55] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[62] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[63] = 2 * REAL(theta)[34];
    }
    else if (theta_n == 36){
      REAL(ret)[63] = 4 * R_pow_di(REAL(theta)[35], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 8));for(int i = 0; i < 8; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 4 * R_pow_di(REAL(theta)[5], 3);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 4 * R_pow_di(REAL(theta)[9], 3);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 4 * R_pow_di(REAL(theta)[14], 3);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 4 * R_pow_di(REAL(theta)[20], 3);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 4 * R_pow_di(REAL(theta)[27], 3);
    }
    else if (theta_n == -38){
      REAL(ret)[7] = 4 * R_pow_di(REAL(theta)[35], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 9){
  int theta_n = INTEGER(tn)[0];
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,45));
    INTEGER(ret)[0]=2;
    INTEGER(ret)[1]=5;
    INTEGER(ret)[2]=2;
    INTEGER(ret)[3]=5;
    INTEGER(ret)[4]=5;
    INTEGER(ret)[5]=2;
    INTEGER(ret)[6]=5;
    INTEGER(ret)[7]=5;
    INTEGER(ret)[8]=5;
    INTEGER(ret)[9]=2;
    INTEGER(ret)[10]=5;
    INTEGER(ret)[11]=5;
    INTEGER(ret)[12]=5;
    INTEGER(ret)[13]=5;
    INTEGER(ret)[14]=2;
    INTEGER(ret)[15]=5;
    INTEGER(ret)[16]=5;
    INTEGER(ret)[17]=5;
    INTEGER(ret)[18]=5;
    INTEGER(ret)[19]=5;
    INTEGER(ret)[20]=2;
    INTEGER(ret)[21]=5;
    INTEGER(ret)[22]=5;
    INTEGER(ret)[23]=5;
    INTEGER(ret)[24]=5;
    INTEGER(ret)[25]=5;
    INTEGER(ret)[26]=5;
    INTEGER(ret)[27]=2;
    INTEGER(ret)[28]=5;
    INTEGER(ret)[29]=5;
    INTEGER(ret)[30]=5;
    INTEGER(ret)[31]=5;
    INTEGER(ret)[32]=5;
    INTEGER(ret)[33]=5;
    INTEGER(ret)[34]=5;
    INTEGER(ret)[35]=2;
    INTEGER(ret)[36]=5;
    INTEGER(ret)[37]=5;
    INTEGER(ret)[38]=5;
    INTEGER(ret)[39]=5;
    INTEGER(ret)[40]=5;
    INTEGER(ret)[41]=5;
    INTEGER(ret)[42]=5;
    INTEGER(ret)[43]=5;
    INTEGER(ret)[44]=2;
    UNPROTECT(1);
    return(ret);  
}

if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 45;
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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[9] = REAL(theta)[1];
      REAL(ret)[10] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[18] = REAL(theta)[3];
      REAL(ret)[19] = REAL(theta)[4];
      REAL(ret)[20] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[27] = REAL(theta)[6];
      REAL(ret)[28] = REAL(theta)[7];
      REAL(ret)[29] = REAL(theta)[8];
      REAL(ret)[30] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[36] = REAL(theta)[10];
      REAL(ret)[37] = REAL(theta)[11];
      REAL(ret)[38] = REAL(theta)[12];
      REAL(ret)[39] = REAL(theta)[13];
      REAL(ret)[40] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[45] = REAL(theta)[15];
      REAL(ret)[46] = REAL(theta)[16];
      REAL(ret)[47] = REAL(theta)[17];
      REAL(ret)[48] = REAL(theta)[18];
      REAL(ret)[49] = REAL(theta)[19];
      REAL(ret)[50] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[54] = REAL(theta)[21];
      REAL(ret)[55] = REAL(theta)[22];
      REAL(ret)[56] = REAL(theta)[23];
      REAL(ret)[57] = REAL(theta)[24];
      REAL(ret)[58] = REAL(theta)[25];
      REAL(ret)[59] = REAL(theta)[26];
      REAL(ret)[60] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[63] = REAL(theta)[28];
      REAL(ret)[64] = REAL(theta)[29];
      REAL(ret)[65] = REAL(theta)[30];
      REAL(ret)[66] = REAL(theta)[31];
      REAL(ret)[67] = REAL(theta)[32];
      REAL(ret)[68] = REAL(theta)[33];
      REAL(ret)[69] = REAL(theta)[34];
      REAL(ret)[70] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[72] = REAL(theta)[36];
      REAL(ret)[73] = REAL(theta)[37];
      REAL(ret)[74] = REAL(theta)[38];
      REAL(ret)[75] = REAL(theta)[39];
      REAL(ret)[76] = REAL(theta)[40];
      REAL(ret)[77] = REAL(theta)[41];
      REAL(ret)[78] = REAL(theta)[42];
      REAL(ret)[79] = REAL(theta)[43];
      REAL(ret)[80] = R_pow_di(REAL(theta)[44], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[28];
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[36];
      REAL(ret)[9] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[10] = R_pow_di(REAL(theta)[1], 2) + R_pow_di(REAL(theta)[2], 4);
      REAL(ret)[11] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[12] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[14] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[16] = REAL(theta)[1] * REAL(theta)[28] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[29];
      REAL(ret)[17] = REAL(theta)[1] * REAL(theta)[36] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[37];
      REAL(ret)[18] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[19] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[20] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + R_pow_di(REAL(theta)[5], 4);
      REAL(ret)[21] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[22] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[23] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[24] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[25] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[26] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[27] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[28] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[29] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[30] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + R_pow_di(REAL(theta)[9], 4);
      REAL(ret)[31] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[32] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[33] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[34] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[35] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[36] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[37] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[38] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[39] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[40] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + R_pow_di(REAL(theta)[14], 4);
      REAL(ret)[41] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[42] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[43] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[32];
      REAL(ret)[44] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[40];
      REAL(ret)[45] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[46] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[47] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[48] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[49] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[50] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + R_pow_di(REAL(theta)[20], 4);
      REAL(ret)[51] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[52] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[33];
      REAL(ret)[53] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[41];
      REAL(ret)[54] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[55] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[56] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[57] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[58] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[59] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[60] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + R_pow_di(REAL(theta)[27], 4);
      REAL(ret)[61] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[34];
      REAL(ret)[62] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[42];
      REAL(ret)[63] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[28];
      REAL(ret)[64] = REAL(theta)[1] * REAL(theta)[28] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[29];
      REAL(ret)[65] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[66] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[67] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[32];
      REAL(ret)[68] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[33];
      REAL(ret)[69] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[34];
      REAL(ret)[70] = R_pow_di(REAL(theta)[28], 2) + R_pow_di(REAL(theta)[29], 2) + R_pow_di(REAL(theta)[30], 2) + R_pow_di(REAL(theta)[31], 2) + R_pow_di(REAL(theta)[32], 2) + R_pow_di(REAL(theta)[33], 2) + R_pow_di(REAL(theta)[34], 2) + R_pow_di(REAL(theta)[35], 4);
      REAL(ret)[71] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[43];
      REAL(ret)[72] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[36];
      REAL(ret)[73] = REAL(theta)[1] * REAL(theta)[36] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[37];
      REAL(ret)[74] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[75] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[76] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[40];
      REAL(ret)[77] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[41];
      REAL(ret)[78] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[42];
      REAL(ret)[79] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[43];
      REAL(ret)[80] = R_pow_di(REAL(theta)[36], 2) + R_pow_di(REAL(theta)[37], 2) + R_pow_di(REAL(theta)[38], 2) + R_pow_di(REAL(theta)[39], 2) + R_pow_di(REAL(theta)[40], 2) + R_pow_di(REAL(theta)[41], 2) + R_pow_di(REAL(theta)[42], 2) + R_pow_di(REAL(theta)[43], 2) + R_pow_di(REAL(theta)[44], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
      REAL(ret)[1] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[2] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[3] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[4] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[5] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[6] = 2 * REAL(theta)[0] * REAL(theta)[21];
      REAL(ret)[7] = 2 * REAL(theta)[0] * REAL(theta)[28];
      REAL(ret)[8] = 2 * REAL(theta)[0] * REAL(theta)[36];
      REAL(ret)[9] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[18] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[27] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[36] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[45] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[54] = 2 * REAL(theta)[0] * REAL(theta)[21];
      REAL(ret)[63] = 2 * REAL(theta)[0] * REAL(theta)[28];
      REAL(ret)[72] = 2 * REAL(theta)[0] * REAL(theta)[36];
    }
    else if (theta_n == 2){
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[9] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[10] = 4 * R_pow_di(REAL(theta)[2], 3);
      REAL(ret)[11] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[12] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[13] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[14] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[15] = 2 * REAL(theta)[2] * REAL(theta)[22];
      REAL(ret)[16] = 2 * REAL(theta)[2] * REAL(theta)[29];
      REAL(ret)[17] = 2 * REAL(theta)[2] * REAL(theta)[37];
      REAL(ret)[19] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[28] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[37] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[46] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[55] = 2 * REAL(theta)[2] * REAL(theta)[22];
      REAL(ret)[64] = 2 * REAL(theta)[2] * REAL(theta)[29];
      REAL(ret)[73] = 2 * REAL(theta)[2] * REAL(theta)[37];
    }
    else if (theta_n == 4){
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[18] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[11] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[19] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[20] = 4 * R_pow_di(REAL(theta)[5], 3);
      REAL(ret)[21] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[22] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[23] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[24] = 2 * REAL(theta)[23] * REAL(theta)[5];
      REAL(ret)[25] = 2 * REAL(theta)[30] * REAL(theta)[5];
      REAL(ret)[26] = 2 * REAL(theta)[38] * REAL(theta)[5];
      REAL(ret)[29] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[38] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[47] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[56] = 2 * REAL(theta)[23] * REAL(theta)[5];
      REAL(ret)[65] = 2 * REAL(theta)[30] * REAL(theta)[5];
      REAL(ret)[74] = 2 * REAL(theta)[38] * REAL(theta)[5];
    }
    else if (theta_n == 7){
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[12] = REAL(theta)[1];
      REAL(ret)[21] = REAL(theta)[3];
      REAL(ret)[27] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[12] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[21] = REAL(theta)[4];
      REAL(ret)[28] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[21] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[29] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[30] = 4 * R_pow_di(REAL(theta)[9], 3);
      REAL(ret)[31] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[32] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[33] = 2 * REAL(theta)[24] * REAL(theta)[9];
      REAL(ret)[34] = 2 * REAL(theta)[31] * REAL(theta)[9];
      REAL(ret)[35] = 2 * REAL(theta)[39] * REAL(theta)[9];
      REAL(ret)[39] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[48] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[57] = 2 * REAL(theta)[24] * REAL(theta)[9];
      REAL(ret)[66] = 2 * REAL(theta)[31] * REAL(theta)[9];
      REAL(ret)[75] = 2 * REAL(theta)[39] * REAL(theta)[9];
    }
    else if (theta_n == 11){
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[22] = REAL(theta)[3];
      REAL(ret)[31] = REAL(theta)[6];
      REAL(ret)[36] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[13] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[22] = REAL(theta)[4];
      REAL(ret)[31] = REAL(theta)[7];
      REAL(ret)[37] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[22] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[31] = REAL(theta)[8];
      REAL(ret)[38] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[31] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[39] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[40] = 4 * R_pow_di(REAL(theta)[14], 3);
      REAL(ret)[41] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[42] = 2 * REAL(theta)[14] * REAL(theta)[25];
      REAL(ret)[43] = 2 * REAL(theta)[14] * REAL(theta)[32];
      REAL(ret)[44] = 2 * REAL(theta)[14] * REAL(theta)[40];
      REAL(ret)[49] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[58] = 2 * REAL(theta)[14] * REAL(theta)[25];
      REAL(ret)[67] = 2 * REAL(theta)[14] * REAL(theta)[32];
      REAL(ret)[76] = 2 * REAL(theta)[14] * REAL(theta)[40];
    }
    else if (theta_n == 16){
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[14] = REAL(theta)[1];
      REAL(ret)[23] = REAL(theta)[3];
      REAL(ret)[32] = REAL(theta)[6];
      REAL(ret)[41] = REAL(theta)[10];
      REAL(ret)[45] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[14] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[23] = REAL(theta)[4];
      REAL(ret)[32] = REAL(theta)[7];
      REAL(ret)[41] = REAL(theta)[11];
      REAL(ret)[46] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[23] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[32] = REAL(theta)[8];
      REAL(ret)[41] = REAL(theta)[12];
      REAL(ret)[47] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[32] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[41] = REAL(theta)[13];
      REAL(ret)[48] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[41] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[49] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[50] = 2 * REAL(theta)[19];
      REAL(ret)[51] = REAL(theta)[25];
      REAL(ret)[52] = REAL(theta)[32];
      REAL(ret)[53] = REAL(theta)[40];
      REAL(ret)[59] = REAL(theta)[25];
      REAL(ret)[68] = REAL(theta)[32];
      REAL(ret)[77] = REAL(theta)[40];
    }
    else if (theta_n == 21){
      REAL(ret)[50] = 4 * R_pow_di(REAL(theta)[20], 3);
      REAL(ret)[51] = 2 * REAL(theta)[20] * REAL(theta)[26];
      REAL(ret)[52] = 2 * REAL(theta)[20] * REAL(theta)[33];
      REAL(ret)[53] = 2 * REAL(theta)[20] * REAL(theta)[41];
      REAL(ret)[59] = 2 * REAL(theta)[20] * REAL(theta)[26];
      REAL(ret)[68] = 2 * REAL(theta)[20] * REAL(theta)[33];
      REAL(ret)[77] = 2 * REAL(theta)[20] * REAL(theta)[41];
    }
    else if (theta_n == 22){
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[15] = REAL(theta)[1];
      REAL(ret)[24] = REAL(theta)[3];
      REAL(ret)[33] = REAL(theta)[6];
      REAL(ret)[42] = REAL(theta)[10];
      REAL(ret)[51] = REAL(theta)[15];
      REAL(ret)[54] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[15] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[24] = REAL(theta)[4];
      REAL(ret)[33] = REAL(theta)[7];
      REAL(ret)[42] = REAL(theta)[11];
      REAL(ret)[51] = REAL(theta)[16];
      REAL(ret)[55] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[24] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[33] = REAL(theta)[8];
      REAL(ret)[42] = REAL(theta)[12];
      REAL(ret)[51] = REAL(theta)[17];
      REAL(ret)[56] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[33] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[42] = REAL(theta)[13];
      REAL(ret)[51] = REAL(theta)[18];
      REAL(ret)[57] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[58] = REAL(theta)[13];
      REAL(ret)[59] = REAL(theta)[18];
      REAL(ret)[60] = 2 * REAL(theta)[24];
      REAL(ret)[61] = REAL(theta)[31];
      REAL(ret)[62] = REAL(theta)[39];
      REAL(ret)[69] = REAL(theta)[31];
      REAL(ret)[78] = REAL(theta)[39];
    }
    else if (theta_n == 26){
      REAL(ret)[42] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[51] = REAL(theta)[19];
      REAL(ret)[58] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[59] = REAL(theta)[19];
      REAL(ret)[60] = 2 * REAL(theta)[25];
      REAL(ret)[61] = REAL(theta)[32];
      REAL(ret)[62] = REAL(theta)[40];
      REAL(ret)[69] = REAL(theta)[32];
      REAL(ret)[78] = REAL(theta)[40];
    }
    else if (theta_n == 27){
      REAL(ret)[51] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[59] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[60] = 2 * REAL(theta)[26];
      REAL(ret)[61] = REAL(theta)[33];
      REAL(ret)[62] = REAL(theta)[41];
      REAL(ret)[69] = REAL(theta)[33];
      REAL(ret)[78] = REAL(theta)[41];
    }
    else if (theta_n == 28){
      REAL(ret)[60] = 4 * R_pow_di(REAL(theta)[27], 3);
      REAL(ret)[61] = 2 * REAL(theta)[27] * REAL(theta)[34];
      REAL(ret)[62] = 2 * REAL(theta)[27] * REAL(theta)[42];
      REAL(ret)[69] = 2 * REAL(theta)[27] * REAL(theta)[34];
      REAL(ret)[78] = 2 * REAL(theta)[27] * REAL(theta)[42];
    }
    else if (theta_n == 29){
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[16] = REAL(theta)[1];
      REAL(ret)[25] = REAL(theta)[3];
      REAL(ret)[34] = REAL(theta)[6];
      REAL(ret)[43] = REAL(theta)[10];
      REAL(ret)[52] = REAL(theta)[15];
      REAL(ret)[61] = REAL(theta)[21];
      REAL(ret)[63] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[16] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[25] = REAL(theta)[4];
      REAL(ret)[34] = REAL(theta)[7];
      REAL(ret)[43] = REAL(theta)[11];
      REAL(ret)[52] = REAL(theta)[16];
      REAL(ret)[61] = REAL(theta)[22];
      REAL(ret)[64] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[25] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[34] = REAL(theta)[8];
      REAL(ret)[43] = REAL(theta)[12];
      REAL(ret)[52] = REAL(theta)[17];
      REAL(ret)[61] = REAL(theta)[23];
      REAL(ret)[65] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[66] = REAL(theta)[8];
      REAL(ret)[67] = REAL(theta)[12];
      REAL(ret)[68] = REAL(theta)[17];
      REAL(ret)[69] = REAL(theta)[23];
      REAL(ret)[70] = 2 * REAL(theta)[30];
      REAL(ret)[71] = REAL(theta)[38];
      REAL(ret)[79] = REAL(theta)[38];
    }
    else if (theta_n == 32){
      REAL(ret)[34] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[43] = REAL(theta)[13];
      REAL(ret)[52] = REAL(theta)[18];
      REAL(ret)[61] = REAL(theta)[24];
      REAL(ret)[66] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[67] = REAL(theta)[13];
      REAL(ret)[68] = REAL(theta)[18];
      REAL(ret)[69] = REAL(theta)[24];
      REAL(ret)[70] = 2 * REAL(theta)[31];
      REAL(ret)[71] = REAL(theta)[39];
      REAL(ret)[79] = REAL(theta)[39];
    }
    else if (theta_n == 33){
      REAL(ret)[43] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[52] = REAL(theta)[19];
      REAL(ret)[61] = REAL(theta)[25];
      REAL(ret)[67] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[68] = REAL(theta)[19];
      REAL(ret)[69] = REAL(theta)[25];
      REAL(ret)[70] = 2 * REAL(theta)[32];
      REAL(ret)[71] = REAL(theta)[40];
      REAL(ret)[79] = REAL(theta)[40];
    }
    else if (theta_n == 34){
      REAL(ret)[52] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[61] = REAL(theta)[26];
      REAL(ret)[68] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[69] = REAL(theta)[26];
      REAL(ret)[70] = 2 * REAL(theta)[33];
      REAL(ret)[71] = REAL(theta)[41];
      REAL(ret)[79] = REAL(theta)[41];
    }
    else if (theta_n == 35){
      REAL(ret)[61] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[69] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[70] = 2 * REAL(theta)[34];
      REAL(ret)[71] = REAL(theta)[42];
      REAL(ret)[79] = REAL(theta)[42];
    }
    else if (theta_n == 36){
      REAL(ret)[70] = 4 * R_pow_di(REAL(theta)[35], 3);
      REAL(ret)[71] = 2 * REAL(theta)[35] * REAL(theta)[43];
      REAL(ret)[79] = 2 * REAL(theta)[35] * REAL(theta)[43];
    }
    else if (theta_n == 37){
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[17] = REAL(theta)[1];
      REAL(ret)[26] = REAL(theta)[3];
      REAL(ret)[35] = REAL(theta)[6];
      REAL(ret)[44] = REAL(theta)[10];
      REAL(ret)[53] = REAL(theta)[15];
      REAL(ret)[62] = REAL(theta)[21];
      REAL(ret)[71] = REAL(theta)[28];
      REAL(ret)[72] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[17] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[26] = REAL(theta)[4];
      REAL(ret)[35] = REAL(theta)[7];
      REAL(ret)[44] = REAL(theta)[11];
      REAL(ret)[53] = REAL(theta)[16];
      REAL(ret)[62] = REAL(theta)[22];
      REAL(ret)[71] = REAL(theta)[29];
      REAL(ret)[73] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[74] = REAL(theta)[4];
      REAL(ret)[75] = REAL(theta)[7];
      REAL(ret)[76] = REAL(theta)[11];
      REAL(ret)[77] = REAL(theta)[16];
      REAL(ret)[78] = REAL(theta)[22];
      REAL(ret)[79] = REAL(theta)[29];
      REAL(ret)[80] = 2 * REAL(theta)[37];
    }
    else if (theta_n == 39){
      REAL(ret)[26] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[35] = REAL(theta)[8];
      REAL(ret)[44] = REAL(theta)[12];
      REAL(ret)[53] = REAL(theta)[17];
      REAL(ret)[62] = REAL(theta)[23];
      REAL(ret)[71] = REAL(theta)[30];
      REAL(ret)[74] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[75] = REAL(theta)[8];
      REAL(ret)[76] = REAL(theta)[12];
      REAL(ret)[77] = REAL(theta)[17];
      REAL(ret)[78] = REAL(theta)[23];
      REAL(ret)[79] = REAL(theta)[30];
      REAL(ret)[80] = 2 * REAL(theta)[38];
    }
    else if (theta_n == 40){
      REAL(ret)[35] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[44] = REAL(theta)[13];
      REAL(ret)[53] = REAL(theta)[18];
      REAL(ret)[62] = REAL(theta)[24];
      REAL(ret)[71] = REAL(theta)[31];
      REAL(ret)[75] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[76] = REAL(theta)[13];
      REAL(ret)[77] = REAL(theta)[18];
      REAL(ret)[78] = REAL(theta)[24];
      REAL(ret)[79] = REAL(theta)[31];
      REAL(ret)[80] = 2 * REAL(theta)[39];
    }
    else if (theta_n == 41){
      REAL(ret)[44] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[53] = REAL(theta)[19];
      REAL(ret)[62] = REAL(theta)[25];
      REAL(ret)[71] = REAL(theta)[32];
      REAL(ret)[76] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[77] = REAL(theta)[19];
      REAL(ret)[78] = REAL(theta)[25];
      REAL(ret)[79] = REAL(theta)[32];
      REAL(ret)[80] = 2 * REAL(theta)[40];
    }
    else if (theta_n == 42){
      REAL(ret)[53] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[62] = REAL(theta)[26];
      REAL(ret)[71] = REAL(theta)[33];
      REAL(ret)[77] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[78] = REAL(theta)[26];
      REAL(ret)[79] = REAL(theta)[33];
      REAL(ret)[80] = 2 * REAL(theta)[41];
    }
    else if (theta_n == 43){
      REAL(ret)[62] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[71] = REAL(theta)[34];
      REAL(ret)[78] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[79] = REAL(theta)[34];
      REAL(ret)[80] = 2 * REAL(theta)[42];
    }
    else if (theta_n == 44){
      REAL(ret)[71] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[79] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[80] = 2 * REAL(theta)[43];
    }
    else if (theta_n == 45){
      REAL(ret)[80] = 4 * R_pow_di(REAL(theta)[44], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 9));for(int i = 0; i < 9; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 4 * R_pow_di(REAL(theta)[5], 3);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 4 * R_pow_di(REAL(theta)[9], 3);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 4 * R_pow_di(REAL(theta)[14], 3);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 4 * R_pow_di(REAL(theta)[20], 3);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 4 * R_pow_di(REAL(theta)[27], 3);
    }
    else if (theta_n == -38){
      REAL(ret)[7] = 4 * R_pow_di(REAL(theta)[35], 3);
    }
    else if (theta_n == -47){
      REAL(ret)[8] = 4 * R_pow_di(REAL(theta)[44], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 10){
  int theta_n = INTEGER(tn)[0];
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,55));
    INTEGER(ret)[0]=2;
    INTEGER(ret)[1]=5;
    INTEGER(ret)[2]=2;
    INTEGER(ret)[3]=5;
    INTEGER(ret)[4]=5;
    INTEGER(ret)[5]=2;
    INTEGER(ret)[6]=5;
    INTEGER(ret)[7]=5;
    INTEGER(ret)[8]=5;
    INTEGER(ret)[9]=2;
    INTEGER(ret)[10]=5;
    INTEGER(ret)[11]=5;
    INTEGER(ret)[12]=5;
    INTEGER(ret)[13]=5;
    INTEGER(ret)[14]=2;
    INTEGER(ret)[15]=5;
    INTEGER(ret)[16]=5;
    INTEGER(ret)[17]=5;
    INTEGER(ret)[18]=5;
    INTEGER(ret)[19]=5;
    INTEGER(ret)[20]=2;
    INTEGER(ret)[21]=5;
    INTEGER(ret)[22]=5;
    INTEGER(ret)[23]=5;
    INTEGER(ret)[24]=5;
    INTEGER(ret)[25]=5;
    INTEGER(ret)[26]=5;
    INTEGER(ret)[27]=2;
    INTEGER(ret)[28]=5;
    INTEGER(ret)[29]=5;
    INTEGER(ret)[30]=5;
    INTEGER(ret)[31]=5;
    INTEGER(ret)[32]=5;
    INTEGER(ret)[33]=5;
    INTEGER(ret)[34]=5;
    INTEGER(ret)[35]=2;
    INTEGER(ret)[36]=5;
    INTEGER(ret)[37]=5;
    INTEGER(ret)[38]=5;
    INTEGER(ret)[39]=5;
    INTEGER(ret)[40]=5;
    INTEGER(ret)[41]=5;
    INTEGER(ret)[42]=5;
    INTEGER(ret)[43]=5;
    INTEGER(ret)[44]=2;
    INTEGER(ret)[45]=5;
    INTEGER(ret)[46]=5;
    INTEGER(ret)[47]=5;
    INTEGER(ret)[48]=5;
    INTEGER(ret)[49]=5;
    INTEGER(ret)[50]=5;
    INTEGER(ret)[51]=5;
    INTEGER(ret)[52]=5;
    INTEGER(ret)[53]=5;
    INTEGER(ret)[54]=2;
    UNPROTECT(1);
    return(ret);  
}

if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 55;
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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[10] = REAL(theta)[1];
      REAL(ret)[11] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[20] = REAL(theta)[3];
      REAL(ret)[21] = REAL(theta)[4];
      REAL(ret)[22] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[30] = REAL(theta)[6];
      REAL(ret)[31] = REAL(theta)[7];
      REAL(ret)[32] = REAL(theta)[8];
      REAL(ret)[33] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[40] = REAL(theta)[10];
      REAL(ret)[41] = REAL(theta)[11];
      REAL(ret)[42] = REAL(theta)[12];
      REAL(ret)[43] = REAL(theta)[13];
      REAL(ret)[44] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[50] = REAL(theta)[15];
      REAL(ret)[51] = REAL(theta)[16];
      REAL(ret)[52] = REAL(theta)[17];
      REAL(ret)[53] = REAL(theta)[18];
      REAL(ret)[54] = REAL(theta)[19];
      REAL(ret)[55] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[60] = REAL(theta)[21];
      REAL(ret)[61] = REAL(theta)[22];
      REAL(ret)[62] = REAL(theta)[23];
      REAL(ret)[63] = REAL(theta)[24];
      REAL(ret)[64] = REAL(theta)[25];
      REAL(ret)[65] = REAL(theta)[26];
      REAL(ret)[66] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[70] = REAL(theta)[28];
      REAL(ret)[71] = REAL(theta)[29];
      REAL(ret)[72] = REAL(theta)[30];
      REAL(ret)[73] = REAL(theta)[31];
      REAL(ret)[74] = REAL(theta)[32];
      REAL(ret)[75] = REAL(theta)[33];
      REAL(ret)[76] = REAL(theta)[34];
      REAL(ret)[77] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[80] = REAL(theta)[36];
      REAL(ret)[81] = REAL(theta)[37];
      REAL(ret)[82] = REAL(theta)[38];
      REAL(ret)[83] = REAL(theta)[39];
      REAL(ret)[84] = REAL(theta)[40];
      REAL(ret)[85] = REAL(theta)[41];
      REAL(ret)[86] = REAL(theta)[42];
      REAL(ret)[87] = REAL(theta)[43];
      REAL(ret)[88] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[90] = REAL(theta)[45];
      REAL(ret)[91] = REAL(theta)[46];
      REAL(ret)[92] = REAL(theta)[47];
      REAL(ret)[93] = REAL(theta)[48];
      REAL(ret)[94] = REAL(theta)[49];
      REAL(ret)[95] = REAL(theta)[50];
      REAL(ret)[96] = REAL(theta)[51];
      REAL(ret)[97] = REAL(theta)[52];
      REAL(ret)[98] = REAL(theta)[53];
      REAL(ret)[99] = R_pow_di(REAL(theta)[54], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[28];
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[36];
      REAL(ret)[9] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[45];
      REAL(ret)[10] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[11] = R_pow_di(REAL(theta)[1], 2) + R_pow_di(REAL(theta)[2], 4);
      REAL(ret)[12] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[14] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[16] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[17] = REAL(theta)[1] * REAL(theta)[28] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[29];
      REAL(ret)[18] = REAL(theta)[1] * REAL(theta)[36] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[37];
      REAL(ret)[19] = REAL(theta)[1] * REAL(theta)[45] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[46];
      REAL(ret)[20] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[21] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[22] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + R_pow_di(REAL(theta)[5], 4);
      REAL(ret)[23] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[24] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[25] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[26] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[27] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[28] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[29] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[30] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[31] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[32] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[33] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + R_pow_di(REAL(theta)[9], 4);
      REAL(ret)[34] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[35] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[36] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[37] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[38] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[39] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[40] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[41] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[42] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[43] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[44] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + R_pow_di(REAL(theta)[14], 4);
      REAL(ret)[45] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[46] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[47] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[32];
      REAL(ret)[48] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[40];
      REAL(ret)[49] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[49];
      REAL(ret)[50] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[51] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[52] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[53] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[54] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[55] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + R_pow_di(REAL(theta)[20], 4);
      REAL(ret)[56] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[57] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[33];
      REAL(ret)[58] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[41];
      REAL(ret)[59] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[50];
      REAL(ret)[60] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[61] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[62] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[63] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[64] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[65] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[66] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + R_pow_di(REAL(theta)[27], 4);
      REAL(ret)[67] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[34];
      REAL(ret)[68] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[42];
      REAL(ret)[69] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[51];
      REAL(ret)[70] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[28];
      REAL(ret)[71] = REAL(theta)[1] * REAL(theta)[28] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[29];
      REAL(ret)[72] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[73] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[74] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[32];
      REAL(ret)[75] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[33];
      REAL(ret)[76] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[34];
      REAL(ret)[77] = R_pow_di(REAL(theta)[28], 2) + R_pow_di(REAL(theta)[29], 2) + R_pow_di(REAL(theta)[30], 2) + R_pow_di(REAL(theta)[31], 2) + R_pow_di(REAL(theta)[32], 2) + R_pow_di(REAL(theta)[33], 2) + R_pow_di(REAL(theta)[34], 2) + R_pow_di(REAL(theta)[35], 4);
      REAL(ret)[78] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[43];
      REAL(ret)[79] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[52];
      REAL(ret)[80] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[36];
      REAL(ret)[81] = REAL(theta)[1] * REAL(theta)[36] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[37];
      REAL(ret)[82] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[83] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[84] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[40];
      REAL(ret)[85] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[41];
      REAL(ret)[86] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[42];
      REAL(ret)[87] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[43];
      REAL(ret)[88] = R_pow_di(REAL(theta)[36], 2) + R_pow_di(REAL(theta)[37], 2) + R_pow_di(REAL(theta)[38], 2) + R_pow_di(REAL(theta)[39], 2) + R_pow_di(REAL(theta)[40], 2) + R_pow_di(REAL(theta)[41], 2) + R_pow_di(REAL(theta)[42], 2) + R_pow_di(REAL(theta)[43], 2) + R_pow_di(REAL(theta)[44], 4);
      REAL(ret)[89] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[53];
      REAL(ret)[90] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[45];
      REAL(ret)[91] = REAL(theta)[1] * REAL(theta)[45] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[46];
      REAL(ret)[92] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[93] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[94] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[49];
      REAL(ret)[95] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[50];
      REAL(ret)[96] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[51];
      REAL(ret)[97] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[52];
      REAL(ret)[98] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[53];
      REAL(ret)[99] = R_pow_di(REAL(theta)[45], 2) + R_pow_di(REAL(theta)[46], 2) + R_pow_di(REAL(theta)[47], 2) + R_pow_di(REAL(theta)[48], 2) + R_pow_di(REAL(theta)[49], 2) + R_pow_di(REAL(theta)[50], 2) + R_pow_di(REAL(theta)[51], 2) + R_pow_di(REAL(theta)[52], 2) + R_pow_di(REAL(theta)[53], 2) + R_pow_di(REAL(theta)[54], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
      REAL(ret)[1] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[2] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[3] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[4] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[5] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[6] = 2 * REAL(theta)[0] * REAL(theta)[21];
      REAL(ret)[7] = 2 * REAL(theta)[0] * REAL(theta)[28];
      REAL(ret)[8] = 2 * REAL(theta)[0] * REAL(theta)[36];
      REAL(ret)[9] = 2 * REAL(theta)[0] * REAL(theta)[45];
      REAL(ret)[10] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[20] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[30] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[40] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[50] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[60] = 2 * REAL(theta)[0] * REAL(theta)[21];
      REAL(ret)[70] = 2 * REAL(theta)[0] * REAL(theta)[28];
      REAL(ret)[80] = 2 * REAL(theta)[0] * REAL(theta)[36];
      REAL(ret)[90] = 2 * REAL(theta)[0] * REAL(theta)[45];
    }
    else if (theta_n == 2){
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[10] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[11] = 4 * R_pow_di(REAL(theta)[2], 3);
      REAL(ret)[12] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[13] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[14] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[15] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[16] = 2 * REAL(theta)[2] * REAL(theta)[22];
      REAL(ret)[17] = 2 * REAL(theta)[2] * REAL(theta)[29];
      REAL(ret)[18] = 2 * REAL(theta)[2] * REAL(theta)[37];
      REAL(ret)[19] = 2 * REAL(theta)[2] * REAL(theta)[46];
      REAL(ret)[21] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[31] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[41] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[51] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[61] = 2 * REAL(theta)[2] * REAL(theta)[22];
      REAL(ret)[71] = 2 * REAL(theta)[2] * REAL(theta)[29];
      REAL(ret)[81] = 2 * REAL(theta)[2] * REAL(theta)[37];
      REAL(ret)[91] = 2 * REAL(theta)[2] * REAL(theta)[46];
    }
    else if (theta_n == 4){
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[12] = REAL(theta)[1];
      REAL(ret)[20] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[12] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[21] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[22] = 4 * R_pow_di(REAL(theta)[5], 3);
      REAL(ret)[23] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[24] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[25] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[26] = 2 * REAL(theta)[23] * REAL(theta)[5];
      REAL(ret)[27] = 2 * REAL(theta)[30] * REAL(theta)[5];
      REAL(ret)[28] = 2 * REAL(theta)[38] * REAL(theta)[5];
      REAL(ret)[29] = 2 * REAL(theta)[47] * REAL(theta)[5];
      REAL(ret)[32] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[42] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[52] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[62] = 2 * REAL(theta)[23] * REAL(theta)[5];
      REAL(ret)[72] = 2 * REAL(theta)[30] * REAL(theta)[5];
      REAL(ret)[82] = 2 * REAL(theta)[38] * REAL(theta)[5];
      REAL(ret)[92] = 2 * REAL(theta)[47] * REAL(theta)[5];
    }
    else if (theta_n == 7){
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[23] = REAL(theta)[3];
      REAL(ret)[30] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[13] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[23] = REAL(theta)[4];
      REAL(ret)[31] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[23] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[32] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[33] = 4 * R_pow_di(REAL(theta)[9], 3);
      REAL(ret)[34] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[35] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[36] = 2 * REAL(theta)[24] * REAL(theta)[9];
      REAL(ret)[37] = 2 * REAL(theta)[31] * REAL(theta)[9];
      REAL(ret)[38] = 2 * REAL(theta)[39] * REAL(theta)[9];
      REAL(ret)[39] = 2 * REAL(theta)[48] * REAL(theta)[9];
      REAL(ret)[43] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[53] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[63] = 2 * REAL(theta)[24] * REAL(theta)[9];
      REAL(ret)[73] = 2 * REAL(theta)[31] * REAL(theta)[9];
      REAL(ret)[83] = 2 * REAL(theta)[39] * REAL(theta)[9];
      REAL(ret)[93] = 2 * REAL(theta)[48] * REAL(theta)[9];
    }
    else if (theta_n == 11){
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[14] = REAL(theta)[1];
      REAL(ret)[24] = REAL(theta)[3];
      REAL(ret)[34] = REAL(theta)[6];
      REAL(ret)[40] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[14] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[24] = REAL(theta)[4];
      REAL(ret)[34] = REAL(theta)[7];
      REAL(ret)[41] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[24] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[34] = REAL(theta)[8];
      REAL(ret)[42] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[34] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[43] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[44] = 4 * R_pow_di(REAL(theta)[14], 3);
      REAL(ret)[45] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[46] = 2 * REAL(theta)[14] * REAL(theta)[25];
      REAL(ret)[47] = 2 * REAL(theta)[14] * REAL(theta)[32];
      REAL(ret)[48] = 2 * REAL(theta)[14] * REAL(theta)[40];
      REAL(ret)[49] = 2 * REAL(theta)[14] * REAL(theta)[49];
      REAL(ret)[54] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[64] = 2 * REAL(theta)[14] * REAL(theta)[25];
      REAL(ret)[74] = 2 * REAL(theta)[14] * REAL(theta)[32];
      REAL(ret)[84] = 2 * REAL(theta)[14] * REAL(theta)[40];
      REAL(ret)[94] = 2 * REAL(theta)[14] * REAL(theta)[49];
    }
    else if (theta_n == 16){
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[15] = REAL(theta)[1];
      REAL(ret)[25] = REAL(theta)[3];
      REAL(ret)[35] = REAL(theta)[6];
      REAL(ret)[45] = REAL(theta)[10];
      REAL(ret)[50] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[15] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[25] = REAL(theta)[4];
      REAL(ret)[35] = REAL(theta)[7];
      REAL(ret)[45] = REAL(theta)[11];
      REAL(ret)[51] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[25] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[35] = REAL(theta)[8];
      REAL(ret)[45] = REAL(theta)[12];
      REAL(ret)[52] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[35] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[45] = REAL(theta)[13];
      REAL(ret)[53] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[45] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[54] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[55] = 4 * R_pow_di(REAL(theta)[20], 3);
      REAL(ret)[56] = 2 * REAL(theta)[20] * REAL(theta)[26];
      REAL(ret)[57] = 2 * REAL(theta)[20] * REAL(theta)[33];
      REAL(ret)[58] = 2 * REAL(theta)[20] * REAL(theta)[41];
      REAL(ret)[59] = 2 * REAL(theta)[20] * REAL(theta)[50];
      REAL(ret)[65] = 2 * REAL(theta)[20] * REAL(theta)[26];
      REAL(ret)[75] = 2 * REAL(theta)[20] * REAL(theta)[33];
      REAL(ret)[85] = 2 * REAL(theta)[20] * REAL(theta)[41];
      REAL(ret)[95] = 2 * REAL(theta)[20] * REAL(theta)[50];
    }
    else if (theta_n == 22){
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[16] = REAL(theta)[1];
      REAL(ret)[26] = REAL(theta)[3];
      REAL(ret)[36] = REAL(theta)[6];
      REAL(ret)[46] = REAL(theta)[10];
      REAL(ret)[56] = REAL(theta)[15];
      REAL(ret)[60] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[16] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[26] = REAL(theta)[4];
      REAL(ret)[36] = REAL(theta)[7];
      REAL(ret)[46] = REAL(theta)[11];
      REAL(ret)[56] = REAL(theta)[16];
      REAL(ret)[61] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[26] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[36] = REAL(theta)[8];
      REAL(ret)[46] = REAL(theta)[12];
      REAL(ret)[56] = REAL(theta)[17];
      REAL(ret)[62] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[36] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[46] = REAL(theta)[13];
      REAL(ret)[56] = REAL(theta)[18];
      REAL(ret)[63] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[46] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[56] = REAL(theta)[19];
      REAL(ret)[64] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[56] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[65] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[66] = 2 * REAL(theta)[26];
      REAL(ret)[67] = REAL(theta)[33];
      REAL(ret)[68] = REAL(theta)[41];
      REAL(ret)[69] = REAL(theta)[50];
      REAL(ret)[76] = REAL(theta)[33];
      REAL(ret)[86] = REAL(theta)[41];
      REAL(ret)[96] = REAL(theta)[50];
    }
    else if (theta_n == 28){
      REAL(ret)[66] = 4 * R_pow_di(REAL(theta)[27], 3);
      REAL(ret)[67] = 2 * REAL(theta)[27] * REAL(theta)[34];
      REAL(ret)[68] = 2 * REAL(theta)[27] * REAL(theta)[42];
      REAL(ret)[69] = 2 * REAL(theta)[27] * REAL(theta)[51];
      REAL(ret)[76] = 2 * REAL(theta)[27] * REAL(theta)[34];
      REAL(ret)[86] = 2 * REAL(theta)[27] * REAL(theta)[42];
      REAL(ret)[96] = 2 * REAL(theta)[27] * REAL(theta)[51];
    }
    else if (theta_n == 29){
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[17] = REAL(theta)[1];
      REAL(ret)[27] = REAL(theta)[3];
      REAL(ret)[37] = REAL(theta)[6];
      REAL(ret)[47] = REAL(theta)[10];
      REAL(ret)[57] = REAL(theta)[15];
      REAL(ret)[67] = REAL(theta)[21];
      REAL(ret)[70] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[17] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[27] = REAL(theta)[4];
      REAL(ret)[37] = REAL(theta)[7];
      REAL(ret)[47] = REAL(theta)[11];
      REAL(ret)[57] = REAL(theta)[16];
      REAL(ret)[67] = REAL(theta)[22];
      REAL(ret)[71] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[27] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[37] = REAL(theta)[8];
      REAL(ret)[47] = REAL(theta)[12];
      REAL(ret)[57] = REAL(theta)[17];
      REAL(ret)[67] = REAL(theta)[23];
      REAL(ret)[72] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[37] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[47] = REAL(theta)[13];
      REAL(ret)[57] = REAL(theta)[18];
      REAL(ret)[67] = REAL(theta)[24];
      REAL(ret)[73] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[47] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[57] = REAL(theta)[19];
      REAL(ret)[67] = REAL(theta)[25];
      REAL(ret)[74] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[75] = REAL(theta)[19];
      REAL(ret)[76] = REAL(theta)[25];
      REAL(ret)[77] = 2 * REAL(theta)[32];
      REAL(ret)[78] = REAL(theta)[40];
      REAL(ret)[79] = REAL(theta)[49];
      REAL(ret)[87] = REAL(theta)[40];
      REAL(ret)[97] = REAL(theta)[49];
    }
    else if (theta_n == 34){
      REAL(ret)[57] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[67] = REAL(theta)[26];
      REAL(ret)[75] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[76] = REAL(theta)[26];
      REAL(ret)[77] = 2 * REAL(theta)[33];
      REAL(ret)[78] = REAL(theta)[41];
      REAL(ret)[79] = REAL(theta)[50];
      REAL(ret)[87] = REAL(theta)[41];
      REAL(ret)[97] = REAL(theta)[50];
    }
    else if (theta_n == 35){
      REAL(ret)[67] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[76] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[77] = 2 * REAL(theta)[34];
      REAL(ret)[78] = REAL(theta)[42];
      REAL(ret)[79] = REAL(theta)[51];
      REAL(ret)[87] = REAL(theta)[42];
      REAL(ret)[97] = REAL(theta)[51];
    }
    else if (theta_n == 36){
      REAL(ret)[77] = 4 * R_pow_di(REAL(theta)[35], 3);
      REAL(ret)[78] = 2 * REAL(theta)[35] * REAL(theta)[43];
      REAL(ret)[79] = 2 * REAL(theta)[35] * REAL(theta)[52];
      REAL(ret)[87] = 2 * REAL(theta)[35] * REAL(theta)[43];
      REAL(ret)[97] = 2 * REAL(theta)[35] * REAL(theta)[52];
    }
    else if (theta_n == 37){
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[18] = REAL(theta)[1];
      REAL(ret)[28] = REAL(theta)[3];
      REAL(ret)[38] = REAL(theta)[6];
      REAL(ret)[48] = REAL(theta)[10];
      REAL(ret)[58] = REAL(theta)[15];
      REAL(ret)[68] = REAL(theta)[21];
      REAL(ret)[78] = REAL(theta)[28];
      REAL(ret)[80] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[18] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[28] = REAL(theta)[4];
      REAL(ret)[38] = REAL(theta)[7];
      REAL(ret)[48] = REAL(theta)[11];
      REAL(ret)[58] = REAL(theta)[16];
      REAL(ret)[68] = REAL(theta)[22];
      REAL(ret)[78] = REAL(theta)[29];
      REAL(ret)[81] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[28] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[38] = REAL(theta)[8];
      REAL(ret)[48] = REAL(theta)[12];
      REAL(ret)[58] = REAL(theta)[17];
      REAL(ret)[68] = REAL(theta)[23];
      REAL(ret)[78] = REAL(theta)[30];
      REAL(ret)[82] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[38] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[48] = REAL(theta)[13];
      REAL(ret)[58] = REAL(theta)[18];
      REAL(ret)[68] = REAL(theta)[24];
      REAL(ret)[78] = REAL(theta)[31];
      REAL(ret)[83] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[84] = REAL(theta)[13];
      REAL(ret)[85] = REAL(theta)[18];
      REAL(ret)[86] = REAL(theta)[24];
      REAL(ret)[87] = REAL(theta)[31];
      REAL(ret)[88] = 2 * REAL(theta)[39];
      REAL(ret)[89] = REAL(theta)[48];
      REAL(ret)[98] = REAL(theta)[48];
    }
    else if (theta_n == 41){
      REAL(ret)[48] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[58] = REAL(theta)[19];
      REAL(ret)[68] = REAL(theta)[25];
      REAL(ret)[78] = REAL(theta)[32];
      REAL(ret)[84] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[85] = REAL(theta)[19];
      REAL(ret)[86] = REAL(theta)[25];
      REAL(ret)[87] = REAL(theta)[32];
      REAL(ret)[88] = 2 * REAL(theta)[40];
      REAL(ret)[89] = REAL(theta)[49];
      REAL(ret)[98] = REAL(theta)[49];
    }
    else if (theta_n == 42){
      REAL(ret)[58] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[68] = REAL(theta)[26];
      REAL(ret)[78] = REAL(theta)[33];
      REAL(ret)[85] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[86] = REAL(theta)[26];
      REAL(ret)[87] = REAL(theta)[33];
      REAL(ret)[88] = 2 * REAL(theta)[41];
      REAL(ret)[89] = REAL(theta)[50];
      REAL(ret)[98] = REAL(theta)[50];
    }
    else if (theta_n == 43){
      REAL(ret)[68] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[78] = REAL(theta)[34];
      REAL(ret)[86] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[87] = REAL(theta)[34];
      REAL(ret)[88] = 2 * REAL(theta)[42];
      REAL(ret)[89] = REAL(theta)[51];
      REAL(ret)[98] = REAL(theta)[51];
    }
    else if (theta_n == 44){
      REAL(ret)[78] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[87] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[88] = 2 * REAL(theta)[43];
      REAL(ret)[89] = REAL(theta)[52];
      REAL(ret)[98] = REAL(theta)[52];
    }
    else if (theta_n == 45){
      REAL(ret)[88] = 4 * R_pow_di(REAL(theta)[44], 3);
      REAL(ret)[89] = 2 * REAL(theta)[44] * REAL(theta)[53];
      REAL(ret)[98] = 2 * REAL(theta)[44] * REAL(theta)[53];
    }
    else if (theta_n == 46){
      REAL(ret)[9] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[19] = REAL(theta)[1];
      REAL(ret)[29] = REAL(theta)[3];
      REAL(ret)[39] = REAL(theta)[6];
      REAL(ret)[49] = REAL(theta)[10];
      REAL(ret)[59] = REAL(theta)[15];
      REAL(ret)[69] = REAL(theta)[21];
      REAL(ret)[79] = REAL(theta)[28];
      REAL(ret)[89] = REAL(theta)[36];
      REAL(ret)[90] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[19] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[29] = REAL(theta)[4];
      REAL(ret)[39] = REAL(theta)[7];
      REAL(ret)[49] = REAL(theta)[11];
      REAL(ret)[59] = REAL(theta)[16];
      REAL(ret)[69] = REAL(theta)[22];
      REAL(ret)[79] = REAL(theta)[29];
      REAL(ret)[89] = REAL(theta)[37];
      REAL(ret)[91] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[29] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[39] = REAL(theta)[8];
      REAL(ret)[49] = REAL(theta)[12];
      REAL(ret)[59] = REAL(theta)[17];
      REAL(ret)[69] = REAL(theta)[23];
      REAL(ret)[79] = REAL(theta)[30];
      REAL(ret)[89] = REAL(theta)[38];
      REAL(ret)[92] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[93] = REAL(theta)[8];
      REAL(ret)[94] = REAL(theta)[12];
      REAL(ret)[95] = REAL(theta)[17];
      REAL(ret)[96] = REAL(theta)[23];
      REAL(ret)[97] = REAL(theta)[30];
      REAL(ret)[98] = REAL(theta)[38];
      REAL(ret)[99] = 2 * REAL(theta)[47];
    }
    else if (theta_n == 49){
      REAL(ret)[39] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[49] = REAL(theta)[13];
      REAL(ret)[59] = REAL(theta)[18];
      REAL(ret)[69] = REAL(theta)[24];
      REAL(ret)[79] = REAL(theta)[31];
      REAL(ret)[89] = REAL(theta)[39];
      REAL(ret)[93] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[94] = REAL(theta)[13];
      REAL(ret)[95] = REAL(theta)[18];
      REAL(ret)[96] = REAL(theta)[24];
      REAL(ret)[97] = REAL(theta)[31];
      REAL(ret)[98] = REAL(theta)[39];
      REAL(ret)[99] = 2 * REAL(theta)[48];
    }
    else if (theta_n == 50){
      REAL(ret)[49] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[59] = REAL(theta)[19];
      REAL(ret)[69] = REAL(theta)[25];
      REAL(ret)[79] = REAL(theta)[32];
      REAL(ret)[89] = REAL(theta)[40];
      REAL(ret)[94] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[95] = REAL(theta)[19];
      REAL(ret)[96] = REAL(theta)[25];
      REAL(ret)[97] = REAL(theta)[32];
      REAL(ret)[98] = REAL(theta)[40];
      REAL(ret)[99] = 2 * REAL(theta)[49];
    }
    else if (theta_n == 51){
      REAL(ret)[59] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[69] = REAL(theta)[26];
      REAL(ret)[79] = REAL(theta)[33];
      REAL(ret)[89] = REAL(theta)[41];
      REAL(ret)[95] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[96] = REAL(theta)[26];
      REAL(ret)[97] = REAL(theta)[33];
      REAL(ret)[98] = REAL(theta)[41];
      REAL(ret)[99] = 2 * REAL(theta)[50];
    }
    else if (theta_n == 52){
      REAL(ret)[69] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[79] = REAL(theta)[34];
      REAL(ret)[89] = REAL(theta)[42];
      REAL(ret)[96] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[97] = REAL(theta)[34];
      REAL(ret)[98] = REAL(theta)[42];
      REAL(ret)[99] = 2 * REAL(theta)[51];
    }
    else if (theta_n == 53){
      REAL(ret)[79] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[89] = REAL(theta)[43];
      REAL(ret)[97] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[98] = REAL(theta)[43];
      REAL(ret)[99] = 2 * REAL(theta)[52];
    }
    else if (theta_n == 54){
      REAL(ret)[89] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[98] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[99] = 2 * REAL(theta)[53];
    }
    else if (theta_n == 55){
      REAL(ret)[99] = 4 * R_pow_di(REAL(theta)[54], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 10));for(int i = 0; i < 10; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 4 * R_pow_di(REAL(theta)[5], 3);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 4 * R_pow_di(REAL(theta)[9], 3);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 4 * R_pow_di(REAL(theta)[14], 3);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 4 * R_pow_di(REAL(theta)[20], 3);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 4 * R_pow_di(REAL(theta)[27], 3);
    }
    else if (theta_n == -38){
      REAL(ret)[7] = 4 * R_pow_di(REAL(theta)[35], 3);
    }
    else if (theta_n == -47){
      REAL(ret)[8] = 4 * R_pow_di(REAL(theta)[44], 3);
    }
    else if (theta_n == -57){
      REAL(ret)[9] = 4 * R_pow_di(REAL(theta)[54], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 11){
  int theta_n = INTEGER(tn)[0];
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,66));
    INTEGER(ret)[0]=2;
    INTEGER(ret)[1]=5;
    INTEGER(ret)[2]=2;
    INTEGER(ret)[3]=5;
    INTEGER(ret)[4]=5;
    INTEGER(ret)[5]=2;
    INTEGER(ret)[6]=5;
    INTEGER(ret)[7]=5;
    INTEGER(ret)[8]=5;
    INTEGER(ret)[9]=2;
    INTEGER(ret)[10]=5;
    INTEGER(ret)[11]=5;
    INTEGER(ret)[12]=5;
    INTEGER(ret)[13]=5;
    INTEGER(ret)[14]=2;
    INTEGER(ret)[15]=5;
    INTEGER(ret)[16]=5;
    INTEGER(ret)[17]=5;
    INTEGER(ret)[18]=5;
    INTEGER(ret)[19]=5;
    INTEGER(ret)[20]=2;
    INTEGER(ret)[21]=5;
    INTEGER(ret)[22]=5;
    INTEGER(ret)[23]=5;
    INTEGER(ret)[24]=5;
    INTEGER(ret)[25]=5;
    INTEGER(ret)[26]=5;
    INTEGER(ret)[27]=2;
    INTEGER(ret)[28]=5;
    INTEGER(ret)[29]=5;
    INTEGER(ret)[30]=5;
    INTEGER(ret)[31]=5;
    INTEGER(ret)[32]=5;
    INTEGER(ret)[33]=5;
    INTEGER(ret)[34]=5;
    INTEGER(ret)[35]=2;
    INTEGER(ret)[36]=5;
    INTEGER(ret)[37]=5;
    INTEGER(ret)[38]=5;
    INTEGER(ret)[39]=5;
    INTEGER(ret)[40]=5;
    INTEGER(ret)[41]=5;
    INTEGER(ret)[42]=5;
    INTEGER(ret)[43]=5;
    INTEGER(ret)[44]=2;
    INTEGER(ret)[45]=5;
    INTEGER(ret)[46]=5;
    INTEGER(ret)[47]=5;
    INTEGER(ret)[48]=5;
    INTEGER(ret)[49]=5;
    INTEGER(ret)[50]=5;
    INTEGER(ret)[51]=5;
    INTEGER(ret)[52]=5;
    INTEGER(ret)[53]=5;
    INTEGER(ret)[54]=2;
    INTEGER(ret)[55]=5;
    INTEGER(ret)[56]=5;
    INTEGER(ret)[57]=5;
    INTEGER(ret)[58]=5;
    INTEGER(ret)[59]=5;
    INTEGER(ret)[60]=5;
    INTEGER(ret)[61]=5;
    INTEGER(ret)[62]=5;
    INTEGER(ret)[63]=5;
    INTEGER(ret)[64]=5;
    INTEGER(ret)[65]=2;
    UNPROTECT(1);
    return(ret);  
}

if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 66;
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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[11] = REAL(theta)[1];
      REAL(ret)[12] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[22] = REAL(theta)[3];
      REAL(ret)[23] = REAL(theta)[4];
      REAL(ret)[24] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[33] = REAL(theta)[6];
      REAL(ret)[34] = REAL(theta)[7];
      REAL(ret)[35] = REAL(theta)[8];
      REAL(ret)[36] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[44] = REAL(theta)[10];
      REAL(ret)[45] = REAL(theta)[11];
      REAL(ret)[46] = REAL(theta)[12];
      REAL(ret)[47] = REAL(theta)[13];
      REAL(ret)[48] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[55] = REAL(theta)[15];
      REAL(ret)[56] = REAL(theta)[16];
      REAL(ret)[57] = REAL(theta)[17];
      REAL(ret)[58] = REAL(theta)[18];
      REAL(ret)[59] = REAL(theta)[19];
      REAL(ret)[60] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[66] = REAL(theta)[21];
      REAL(ret)[67] = REAL(theta)[22];
      REAL(ret)[68] = REAL(theta)[23];
      REAL(ret)[69] = REAL(theta)[24];
      REAL(ret)[70] = REAL(theta)[25];
      REAL(ret)[71] = REAL(theta)[26];
      REAL(ret)[72] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[77] = REAL(theta)[28];
      REAL(ret)[78] = REAL(theta)[29];
      REAL(ret)[79] = REAL(theta)[30];
      REAL(ret)[80] = REAL(theta)[31];
      REAL(ret)[81] = REAL(theta)[32];
      REAL(ret)[82] = REAL(theta)[33];
      REAL(ret)[83] = REAL(theta)[34];
      REAL(ret)[84] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[88] = REAL(theta)[36];
      REAL(ret)[89] = REAL(theta)[37];
      REAL(ret)[90] = REAL(theta)[38];
      REAL(ret)[91] = REAL(theta)[39];
      REAL(ret)[92] = REAL(theta)[40];
      REAL(ret)[93] = REAL(theta)[41];
      REAL(ret)[94] = REAL(theta)[42];
      REAL(ret)[95] = REAL(theta)[43];
      REAL(ret)[96] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[99] = REAL(theta)[45];
      REAL(ret)[100] = REAL(theta)[46];
      REAL(ret)[101] = REAL(theta)[47];
      REAL(ret)[102] = REAL(theta)[48];
      REAL(ret)[103] = REAL(theta)[49];
      REAL(ret)[104] = REAL(theta)[50];
      REAL(ret)[105] = REAL(theta)[51];
      REAL(ret)[106] = REAL(theta)[52];
      REAL(ret)[107] = REAL(theta)[53];
      REAL(ret)[108] = R_pow_di(REAL(theta)[54], 2);
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
      REAL(ret)[120] = R_pow_di(REAL(theta)[65], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[28];
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[36];
      REAL(ret)[9] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[45];
      REAL(ret)[10] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[55];
      REAL(ret)[11] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[12] = R_pow_di(REAL(theta)[1], 2) + R_pow_di(REAL(theta)[2], 4);
      REAL(ret)[13] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[14] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[16] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[17] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[18] = REAL(theta)[1] * REAL(theta)[28] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[29];
      REAL(ret)[19] = REAL(theta)[1] * REAL(theta)[36] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[37];
      REAL(ret)[20] = REAL(theta)[1] * REAL(theta)[45] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[46];
      REAL(ret)[21] = REAL(theta)[1] * REAL(theta)[55] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[56];
      REAL(ret)[22] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[23] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[24] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + R_pow_di(REAL(theta)[5], 4);
      REAL(ret)[25] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[26] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[27] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[28] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[29] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[30] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[31] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[32] = REAL(theta)[3] * REAL(theta)[55] + REAL(theta)[4] * REAL(theta)[56] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[57];
      REAL(ret)[33] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[34] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[35] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[36] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + R_pow_di(REAL(theta)[9], 4);
      REAL(ret)[37] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[38] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[39] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[40] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[41] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[42] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[43] = REAL(theta)[55] * REAL(theta)[6] + REAL(theta)[56] * REAL(theta)[7] + REAL(theta)[57] * REAL(theta)[8] + REAL(theta)[58] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[44] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[45] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[46] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[47] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[48] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + R_pow_di(REAL(theta)[14], 4);
      REAL(ret)[49] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[50] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[51] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[32];
      REAL(ret)[52] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[40];
      REAL(ret)[53] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[49];
      REAL(ret)[54] = REAL(theta)[10] * REAL(theta)[55] + REAL(theta)[11] * REAL(theta)[56] + REAL(theta)[12] * REAL(theta)[57] + REAL(theta)[13] * REAL(theta)[58] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[59];
      REAL(ret)[55] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[56] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[57] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[58] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[59] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[60] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + R_pow_di(REAL(theta)[20], 4);
      REAL(ret)[61] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[62] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[33];
      REAL(ret)[63] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[41];
      REAL(ret)[64] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[50];
      REAL(ret)[65] = REAL(theta)[15] * REAL(theta)[55] + REAL(theta)[16] * REAL(theta)[56] + REAL(theta)[17] * REAL(theta)[57] + REAL(theta)[18] * REAL(theta)[58] + REAL(theta)[19] * REAL(theta)[59] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[60];
      REAL(ret)[66] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[67] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[68] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[69] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[70] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[71] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[72] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + R_pow_di(REAL(theta)[27], 4);
      REAL(ret)[73] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[34];
      REAL(ret)[74] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[42];
      REAL(ret)[75] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[51];
      REAL(ret)[76] = REAL(theta)[21] * REAL(theta)[55] + REAL(theta)[22] * REAL(theta)[56] + REAL(theta)[23] * REAL(theta)[57] + REAL(theta)[24] * REAL(theta)[58] + REAL(theta)[25] * REAL(theta)[59] + REAL(theta)[26] * REAL(theta)[60] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[61];
      REAL(ret)[77] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[28];
      REAL(ret)[78] = REAL(theta)[1] * REAL(theta)[28] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[29];
      REAL(ret)[79] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[80] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[81] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[32];
      REAL(ret)[82] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[33];
      REAL(ret)[83] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[34];
      REAL(ret)[84] = R_pow_di(REAL(theta)[28], 2) + R_pow_di(REAL(theta)[29], 2) + R_pow_di(REAL(theta)[30], 2) + R_pow_di(REAL(theta)[31], 2) + R_pow_di(REAL(theta)[32], 2) + R_pow_di(REAL(theta)[33], 2) + R_pow_di(REAL(theta)[34], 2) + R_pow_di(REAL(theta)[35], 4);
      REAL(ret)[85] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[43];
      REAL(ret)[86] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[52];
      REAL(ret)[87] = REAL(theta)[28] * REAL(theta)[55] + REAL(theta)[29] * REAL(theta)[56] + REAL(theta)[30] * REAL(theta)[57] + REAL(theta)[31] * REAL(theta)[58] + REAL(theta)[32] * REAL(theta)[59] + REAL(theta)[33] * REAL(theta)[60] + REAL(theta)[34] * REAL(theta)[61] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[62];
      REAL(ret)[88] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[36];
      REAL(ret)[89] = REAL(theta)[1] * REAL(theta)[36] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[37];
      REAL(ret)[90] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[91] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[92] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[40];
      REAL(ret)[93] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[41];
      REAL(ret)[94] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[42];
      REAL(ret)[95] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[43];
      REAL(ret)[96] = R_pow_di(REAL(theta)[36], 2) + R_pow_di(REAL(theta)[37], 2) + R_pow_di(REAL(theta)[38], 2) + R_pow_di(REAL(theta)[39], 2) + R_pow_di(REAL(theta)[40], 2) + R_pow_di(REAL(theta)[41], 2) + R_pow_di(REAL(theta)[42], 2) + R_pow_di(REAL(theta)[43], 2) + R_pow_di(REAL(theta)[44], 4);
      REAL(ret)[97] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[53];
      REAL(ret)[98] = REAL(theta)[36] * REAL(theta)[55] + REAL(theta)[37] * REAL(theta)[56] + REAL(theta)[38] * REAL(theta)[57] + REAL(theta)[39] * REAL(theta)[58] + REAL(theta)[40] * REAL(theta)[59] + REAL(theta)[41] * REAL(theta)[60] + REAL(theta)[42] * REAL(theta)[61] + REAL(theta)[43] * REAL(theta)[62] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[63];
      REAL(ret)[99] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[45];
      REAL(ret)[100] = REAL(theta)[1] * REAL(theta)[45] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[46];
      REAL(ret)[101] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[102] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[103] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[49];
      REAL(ret)[104] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[50];
      REAL(ret)[105] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[51];
      REAL(ret)[106] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[52];
      REAL(ret)[107] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[53];
      REAL(ret)[108] = R_pow_di(REAL(theta)[45], 2) + R_pow_di(REAL(theta)[46], 2) + R_pow_di(REAL(theta)[47], 2) + R_pow_di(REAL(theta)[48], 2) + R_pow_di(REAL(theta)[49], 2) + R_pow_di(REAL(theta)[50], 2) + R_pow_di(REAL(theta)[51], 2) + R_pow_di(REAL(theta)[52], 2) + R_pow_di(REAL(theta)[53], 2) + R_pow_di(REAL(theta)[54], 4);
      REAL(ret)[109] = REAL(theta)[45] * REAL(theta)[55] + REAL(theta)[46] * REAL(theta)[56] + REAL(theta)[47] * REAL(theta)[57] + REAL(theta)[48] * REAL(theta)[58] + REAL(theta)[49] * REAL(theta)[59] + REAL(theta)[50] * REAL(theta)[60] + REAL(theta)[51] * REAL(theta)[61] + REAL(theta)[52] * REAL(theta)[62] + REAL(theta)[53] * REAL(theta)[63] + R_pow_di(REAL(theta)[54], 2) * REAL(theta)[64];
      REAL(ret)[110] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[55];
      REAL(ret)[111] = REAL(theta)[1] * REAL(theta)[55] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[56];
      REAL(ret)[112] = REAL(theta)[3] * REAL(theta)[55] + REAL(theta)[4] * REAL(theta)[56] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[57];
      REAL(ret)[113] = REAL(theta)[55] * REAL(theta)[6] + REAL(theta)[56] * REAL(theta)[7] + REAL(theta)[57] * REAL(theta)[8] + REAL(theta)[58] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[114] = REAL(theta)[10] * REAL(theta)[55] + REAL(theta)[11] * REAL(theta)[56] + REAL(theta)[12] * REAL(theta)[57] + REAL(theta)[13] * REAL(theta)[58] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[59];
      REAL(ret)[115] = REAL(theta)[15] * REAL(theta)[55] + REAL(theta)[16] * REAL(theta)[56] + REAL(theta)[17] * REAL(theta)[57] + REAL(theta)[18] * REAL(theta)[58] + REAL(theta)[19] * REAL(theta)[59] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[60];
      REAL(ret)[116] = REAL(theta)[21] * REAL(theta)[55] + REAL(theta)[22] * REAL(theta)[56] + REAL(theta)[23] * REAL(theta)[57] + REAL(theta)[24] * REAL(theta)[58] + REAL(theta)[25] * REAL(theta)[59] + REAL(theta)[26] * REAL(theta)[60] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[61];
      REAL(ret)[117] = REAL(theta)[28] * REAL(theta)[55] + REAL(theta)[29] * REAL(theta)[56] + REAL(theta)[30] * REAL(theta)[57] + REAL(theta)[31] * REAL(theta)[58] + REAL(theta)[32] * REAL(theta)[59] + REAL(theta)[33] * REAL(theta)[60] + REAL(theta)[34] * REAL(theta)[61] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[62];
      REAL(ret)[118] = REAL(theta)[36] * REAL(theta)[55] + REAL(theta)[37] * REAL(theta)[56] + REAL(theta)[38] * REAL(theta)[57] + REAL(theta)[39] * REAL(theta)[58] + REAL(theta)[40] * REAL(theta)[59] + REAL(theta)[41] * REAL(theta)[60] + REAL(theta)[42] * REAL(theta)[61] + REAL(theta)[43] * REAL(theta)[62] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[63];
      REAL(ret)[119] = REAL(theta)[45] * REAL(theta)[55] + REAL(theta)[46] * REAL(theta)[56] + REAL(theta)[47] * REAL(theta)[57] + REAL(theta)[48] * REAL(theta)[58] + REAL(theta)[49] * REAL(theta)[59] + REAL(theta)[50] * REAL(theta)[60] + REAL(theta)[51] * REAL(theta)[61] + REAL(theta)[52] * REAL(theta)[62] + REAL(theta)[53] * REAL(theta)[63] + R_pow_di(REAL(theta)[54], 2) * REAL(theta)[64];
      REAL(ret)[120] = R_pow_di(REAL(theta)[55], 2) + R_pow_di(REAL(theta)[56], 2) + R_pow_di(REAL(theta)[57], 2) + R_pow_di(REAL(theta)[58], 2) + R_pow_di(REAL(theta)[59], 2) + R_pow_di(REAL(theta)[60], 2) + R_pow_di(REAL(theta)[61], 2) + R_pow_di(REAL(theta)[62], 2) + R_pow_di(REAL(theta)[63], 2) + R_pow_di(REAL(theta)[64], 2) + R_pow_di(REAL(theta)[65], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
      REAL(ret)[1] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[2] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[3] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[4] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[5] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[6] = 2 * REAL(theta)[0] * REAL(theta)[21];
      REAL(ret)[7] = 2 * REAL(theta)[0] * REAL(theta)[28];
      REAL(ret)[8] = 2 * REAL(theta)[0] * REAL(theta)[36];
      REAL(ret)[9] = 2 * REAL(theta)[0] * REAL(theta)[45];
      REAL(ret)[10] = 2 * REAL(theta)[0] * REAL(theta)[55];
      REAL(ret)[11] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[22] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[33] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[44] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[55] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[66] = 2 * REAL(theta)[0] * REAL(theta)[21];
      REAL(ret)[77] = 2 * REAL(theta)[0] * REAL(theta)[28];
      REAL(ret)[88] = 2 * REAL(theta)[0] * REAL(theta)[36];
      REAL(ret)[99] = 2 * REAL(theta)[0] * REAL(theta)[45];
      REAL(ret)[110] = 2 * REAL(theta)[0] * REAL(theta)[55];
    }
    else if (theta_n == 2){
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[11] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[12] = 4 * R_pow_di(REAL(theta)[2], 3);
      REAL(ret)[13] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[14] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[15] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[16] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[17] = 2 * REAL(theta)[2] * REAL(theta)[22];
      REAL(ret)[18] = 2 * REAL(theta)[2] * REAL(theta)[29];
      REAL(ret)[19] = 2 * REAL(theta)[2] * REAL(theta)[37];
      REAL(ret)[20] = 2 * REAL(theta)[2] * REAL(theta)[46];
      REAL(ret)[21] = 2 * REAL(theta)[2] * REAL(theta)[56];
      REAL(ret)[23] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[34] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[45] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[56] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[67] = 2 * REAL(theta)[2] * REAL(theta)[22];
      REAL(ret)[78] = 2 * REAL(theta)[2] * REAL(theta)[29];
      REAL(ret)[89] = 2 * REAL(theta)[2] * REAL(theta)[37];
      REAL(ret)[100] = 2 * REAL(theta)[2] * REAL(theta)[46];
      REAL(ret)[111] = 2 * REAL(theta)[2] * REAL(theta)[56];
    }
    else if (theta_n == 4){
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[13] = REAL(theta)[1];
      REAL(ret)[22] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[13] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[23] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[24] = 4 * R_pow_di(REAL(theta)[5], 3);
      REAL(ret)[25] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[26] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[27] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[28] = 2 * REAL(theta)[23] * REAL(theta)[5];
      REAL(ret)[29] = 2 * REAL(theta)[30] * REAL(theta)[5];
      REAL(ret)[30] = 2 * REAL(theta)[38] * REAL(theta)[5];
      REAL(ret)[31] = 2 * REAL(theta)[47] * REAL(theta)[5];
      REAL(ret)[32] = 2 * REAL(theta)[5] * REAL(theta)[57];
      REAL(ret)[35] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[46] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[57] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[68] = 2 * REAL(theta)[23] * REAL(theta)[5];
      REAL(ret)[79] = 2 * REAL(theta)[30] * REAL(theta)[5];
      REAL(ret)[90] = 2 * REAL(theta)[38] * REAL(theta)[5];
      REAL(ret)[101] = 2 * REAL(theta)[47] * REAL(theta)[5];
      REAL(ret)[112] = 2 * REAL(theta)[5] * REAL(theta)[57];
    }
    else if (theta_n == 7){
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[14] = REAL(theta)[1];
      REAL(ret)[25] = REAL(theta)[3];
      REAL(ret)[33] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[14] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[25] = REAL(theta)[4];
      REAL(ret)[34] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[25] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[35] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[36] = 4 * R_pow_di(REAL(theta)[9], 3);
      REAL(ret)[37] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[38] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[39] = 2 * REAL(theta)[24] * REAL(theta)[9];
      REAL(ret)[40] = 2 * REAL(theta)[31] * REAL(theta)[9];
      REAL(ret)[41] = 2 * REAL(theta)[39] * REAL(theta)[9];
      REAL(ret)[42] = 2 * REAL(theta)[48] * REAL(theta)[9];
      REAL(ret)[43] = 2 * REAL(theta)[58] * REAL(theta)[9];
      REAL(ret)[47] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[58] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[69] = 2 * REAL(theta)[24] * REAL(theta)[9];
      REAL(ret)[80] = 2 * REAL(theta)[31] * REAL(theta)[9];
      REAL(ret)[91] = 2 * REAL(theta)[39] * REAL(theta)[9];
      REAL(ret)[102] = 2 * REAL(theta)[48] * REAL(theta)[9];
      REAL(ret)[113] = 2 * REAL(theta)[58] * REAL(theta)[9];
    }
    else if (theta_n == 11){
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[15] = REAL(theta)[1];
      REAL(ret)[26] = REAL(theta)[3];
      REAL(ret)[37] = REAL(theta)[6];
      REAL(ret)[44] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[15] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[26] = REAL(theta)[4];
      REAL(ret)[37] = REAL(theta)[7];
      REAL(ret)[45] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[26] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[37] = REAL(theta)[8];
      REAL(ret)[46] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[37] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[47] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[48] = 4 * R_pow_di(REAL(theta)[14], 3);
      REAL(ret)[49] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[50] = 2 * REAL(theta)[14] * REAL(theta)[25];
      REAL(ret)[51] = 2 * REAL(theta)[14] * REAL(theta)[32];
      REAL(ret)[52] = 2 * REAL(theta)[14] * REAL(theta)[40];
      REAL(ret)[53] = 2 * REAL(theta)[14] * REAL(theta)[49];
      REAL(ret)[54] = 2 * REAL(theta)[14] * REAL(theta)[59];
      REAL(ret)[59] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[70] = 2 * REAL(theta)[14] * REAL(theta)[25];
      REAL(ret)[81] = 2 * REAL(theta)[14] * REAL(theta)[32];
      REAL(ret)[92] = 2 * REAL(theta)[14] * REAL(theta)[40];
      REAL(ret)[103] = 2 * REAL(theta)[14] * REAL(theta)[49];
      REAL(ret)[114] = 2 * REAL(theta)[14] * REAL(theta)[59];
    }
    else if (theta_n == 16){
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[16] = REAL(theta)[1];
      REAL(ret)[27] = REAL(theta)[3];
      REAL(ret)[38] = REAL(theta)[6];
      REAL(ret)[49] = REAL(theta)[10];
      REAL(ret)[55] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[16] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[27] = REAL(theta)[4];
      REAL(ret)[38] = REAL(theta)[7];
      REAL(ret)[49] = REAL(theta)[11];
      REAL(ret)[56] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[27] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[38] = REAL(theta)[8];
      REAL(ret)[49] = REAL(theta)[12];
      REAL(ret)[57] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[38] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[49] = REAL(theta)[13];
      REAL(ret)[58] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[49] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[59] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[60] = 4 * R_pow_di(REAL(theta)[20], 3);
      REAL(ret)[61] = 2 * REAL(theta)[20] * REAL(theta)[26];
      REAL(ret)[62] = 2 * REAL(theta)[20] * REAL(theta)[33];
      REAL(ret)[63] = 2 * REAL(theta)[20] * REAL(theta)[41];
      REAL(ret)[64] = 2 * REAL(theta)[20] * REAL(theta)[50];
      REAL(ret)[65] = 2 * REAL(theta)[20] * REAL(theta)[60];
      REAL(ret)[71] = 2 * REAL(theta)[20] * REAL(theta)[26];
      REAL(ret)[82] = 2 * REAL(theta)[20] * REAL(theta)[33];
      REAL(ret)[93] = 2 * REAL(theta)[20] * REAL(theta)[41];
      REAL(ret)[104] = 2 * REAL(theta)[20] * REAL(theta)[50];
      REAL(ret)[115] = 2 * REAL(theta)[20] * REAL(theta)[60];
    }
    else if (theta_n == 22){
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[17] = REAL(theta)[1];
      REAL(ret)[28] = REAL(theta)[3];
      REAL(ret)[39] = REAL(theta)[6];
      REAL(ret)[50] = REAL(theta)[10];
      REAL(ret)[61] = REAL(theta)[15];
      REAL(ret)[66] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[17] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[28] = REAL(theta)[4];
      REAL(ret)[39] = REAL(theta)[7];
      REAL(ret)[50] = REAL(theta)[11];
      REAL(ret)[61] = REAL(theta)[16];
      REAL(ret)[67] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[28] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[39] = REAL(theta)[8];
      REAL(ret)[50] = REAL(theta)[12];
      REAL(ret)[61] = REAL(theta)[17];
      REAL(ret)[68] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[39] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[50] = REAL(theta)[13];
      REAL(ret)[61] = REAL(theta)[18];
      REAL(ret)[69] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[50] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[61] = REAL(theta)[19];
      REAL(ret)[70] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[61] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[71] = R_pow_di(REAL(theta)[20], 2);
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
      REAL(ret)[72] = 4 * R_pow_di(REAL(theta)[27], 3);
      REAL(ret)[73] = 2 * REAL(theta)[27] * REAL(theta)[34];
      REAL(ret)[74] = 2 * REAL(theta)[27] * REAL(theta)[42];
      REAL(ret)[75] = 2 * REAL(theta)[27] * REAL(theta)[51];
      REAL(ret)[76] = 2 * REAL(theta)[27] * REAL(theta)[61];
      REAL(ret)[83] = 2 * REAL(theta)[27] * REAL(theta)[34];
      REAL(ret)[94] = 2 * REAL(theta)[27] * REAL(theta)[42];
      REAL(ret)[105] = 2 * REAL(theta)[27] * REAL(theta)[51];
      REAL(ret)[116] = 2 * REAL(theta)[27] * REAL(theta)[61];
    }
    else if (theta_n == 29){
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[18] = REAL(theta)[1];
      REAL(ret)[29] = REAL(theta)[3];
      REAL(ret)[40] = REAL(theta)[6];
      REAL(ret)[51] = REAL(theta)[10];
      REAL(ret)[62] = REAL(theta)[15];
      REAL(ret)[73] = REAL(theta)[21];
      REAL(ret)[77] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[18] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[29] = REAL(theta)[4];
      REAL(ret)[40] = REAL(theta)[7];
      REAL(ret)[51] = REAL(theta)[11];
      REAL(ret)[62] = REAL(theta)[16];
      REAL(ret)[73] = REAL(theta)[22];
      REAL(ret)[78] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[29] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[40] = REAL(theta)[8];
      REAL(ret)[51] = REAL(theta)[12];
      REAL(ret)[62] = REAL(theta)[17];
      REAL(ret)[73] = REAL(theta)[23];
      REAL(ret)[79] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[40] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[51] = REAL(theta)[13];
      REAL(ret)[62] = REAL(theta)[18];
      REAL(ret)[73] = REAL(theta)[24];
      REAL(ret)[80] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[51] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[62] = REAL(theta)[19];
      REAL(ret)[73] = REAL(theta)[25];
      REAL(ret)[81] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[62] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[73] = REAL(theta)[26];
      REAL(ret)[82] = R_pow_di(REAL(theta)[20], 2);
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
      REAL(ret)[73] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[83] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[84] = 2 * REAL(theta)[34];
      REAL(ret)[85] = REAL(theta)[42];
      REAL(ret)[86] = REAL(theta)[51];
      REAL(ret)[87] = REAL(theta)[61];
      REAL(ret)[95] = REAL(theta)[42];
      REAL(ret)[106] = REAL(theta)[51];
      REAL(ret)[117] = REAL(theta)[61];
    }
    else if (theta_n == 36){
      REAL(ret)[84] = 4 * R_pow_di(REAL(theta)[35], 3);
      REAL(ret)[85] = 2 * REAL(theta)[35] * REAL(theta)[43];
      REAL(ret)[86] = 2 * REAL(theta)[35] * REAL(theta)[52];
      REAL(ret)[87] = 2 * REAL(theta)[35] * REAL(theta)[62];
      REAL(ret)[95] = 2 * REAL(theta)[35] * REAL(theta)[43];
      REAL(ret)[106] = 2 * REAL(theta)[35] * REAL(theta)[52];
      REAL(ret)[117] = 2 * REAL(theta)[35] * REAL(theta)[62];
    }
    else if (theta_n == 37){
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[19] = REAL(theta)[1];
      REAL(ret)[30] = REAL(theta)[3];
      REAL(ret)[41] = REAL(theta)[6];
      REAL(ret)[52] = REAL(theta)[10];
      REAL(ret)[63] = REAL(theta)[15];
      REAL(ret)[74] = REAL(theta)[21];
      REAL(ret)[85] = REAL(theta)[28];
      REAL(ret)[88] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[19] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[30] = REAL(theta)[4];
      REAL(ret)[41] = REAL(theta)[7];
      REAL(ret)[52] = REAL(theta)[11];
      REAL(ret)[63] = REAL(theta)[16];
      REAL(ret)[74] = REAL(theta)[22];
      REAL(ret)[85] = REAL(theta)[29];
      REAL(ret)[89] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[30] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[41] = REAL(theta)[8];
      REAL(ret)[52] = REAL(theta)[12];
      REAL(ret)[63] = REAL(theta)[17];
      REAL(ret)[74] = REAL(theta)[23];
      REAL(ret)[85] = REAL(theta)[30];
      REAL(ret)[90] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[41] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[52] = REAL(theta)[13];
      REAL(ret)[63] = REAL(theta)[18];
      REAL(ret)[74] = REAL(theta)[24];
      REAL(ret)[85] = REAL(theta)[31];
      REAL(ret)[91] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[52] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[63] = REAL(theta)[19];
      REAL(ret)[74] = REAL(theta)[25];
      REAL(ret)[85] = REAL(theta)[32];
      REAL(ret)[92] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[63] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[74] = REAL(theta)[26];
      REAL(ret)[85] = REAL(theta)[33];
      REAL(ret)[93] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[94] = REAL(theta)[26];
      REAL(ret)[95] = REAL(theta)[33];
      REAL(ret)[96] = 2 * REAL(theta)[41];
      REAL(ret)[97] = REAL(theta)[50];
      REAL(ret)[98] = REAL(theta)[60];
      REAL(ret)[107] = REAL(theta)[50];
      REAL(ret)[118] = REAL(theta)[60];
    }
    else if (theta_n == 43){
      REAL(ret)[74] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[85] = REAL(theta)[34];
      REAL(ret)[94] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[95] = REAL(theta)[34];
      REAL(ret)[96] = 2 * REAL(theta)[42];
      REAL(ret)[97] = REAL(theta)[51];
      REAL(ret)[98] = REAL(theta)[61];
      REAL(ret)[107] = REAL(theta)[51];
      REAL(ret)[118] = REAL(theta)[61];
    }
    else if (theta_n == 44){
      REAL(ret)[85] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[95] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[96] = 2 * REAL(theta)[43];
      REAL(ret)[97] = REAL(theta)[52];
      REAL(ret)[98] = REAL(theta)[62];
      REAL(ret)[107] = REAL(theta)[52];
      REAL(ret)[118] = REAL(theta)[62];
    }
    else if (theta_n == 45){
      REAL(ret)[96] = 4 * R_pow_di(REAL(theta)[44], 3);
      REAL(ret)[97] = 2 * REAL(theta)[44] * REAL(theta)[53];
      REAL(ret)[98] = 2 * REAL(theta)[44] * REAL(theta)[63];
      REAL(ret)[107] = 2 * REAL(theta)[44] * REAL(theta)[53];
      REAL(ret)[118] = 2 * REAL(theta)[44] * REAL(theta)[63];
    }
    else if (theta_n == 46){
      REAL(ret)[9] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[20] = REAL(theta)[1];
      REAL(ret)[31] = REAL(theta)[3];
      REAL(ret)[42] = REAL(theta)[6];
      REAL(ret)[53] = REAL(theta)[10];
      REAL(ret)[64] = REAL(theta)[15];
      REAL(ret)[75] = REAL(theta)[21];
      REAL(ret)[86] = REAL(theta)[28];
      REAL(ret)[97] = REAL(theta)[36];
      REAL(ret)[99] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[20] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[31] = REAL(theta)[4];
      REAL(ret)[42] = REAL(theta)[7];
      REAL(ret)[53] = REAL(theta)[11];
      REAL(ret)[64] = REAL(theta)[16];
      REAL(ret)[75] = REAL(theta)[22];
      REAL(ret)[86] = REAL(theta)[29];
      REAL(ret)[97] = REAL(theta)[37];
      REAL(ret)[100] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[31] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[42] = REAL(theta)[8];
      REAL(ret)[53] = REAL(theta)[12];
      REAL(ret)[64] = REAL(theta)[17];
      REAL(ret)[75] = REAL(theta)[23];
      REAL(ret)[86] = REAL(theta)[30];
      REAL(ret)[97] = REAL(theta)[38];
      REAL(ret)[101] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[42] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[53] = REAL(theta)[13];
      REAL(ret)[64] = REAL(theta)[18];
      REAL(ret)[75] = REAL(theta)[24];
      REAL(ret)[86] = REAL(theta)[31];
      REAL(ret)[97] = REAL(theta)[39];
      REAL(ret)[102] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[53] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[64] = REAL(theta)[19];
      REAL(ret)[75] = REAL(theta)[25];
      REAL(ret)[86] = REAL(theta)[32];
      REAL(ret)[97] = REAL(theta)[40];
      REAL(ret)[103] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[104] = REAL(theta)[19];
      REAL(ret)[105] = REAL(theta)[25];
      REAL(ret)[106] = REAL(theta)[32];
      REAL(ret)[107] = REAL(theta)[40];
      REAL(ret)[108] = 2 * REAL(theta)[49];
      REAL(ret)[109] = REAL(theta)[59];
      REAL(ret)[119] = REAL(theta)[59];
    }
    else if (theta_n == 51){
      REAL(ret)[64] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[75] = REAL(theta)[26];
      REAL(ret)[86] = REAL(theta)[33];
      REAL(ret)[97] = REAL(theta)[41];
      REAL(ret)[104] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[105] = REAL(theta)[26];
      REAL(ret)[106] = REAL(theta)[33];
      REAL(ret)[107] = REAL(theta)[41];
      REAL(ret)[108] = 2 * REAL(theta)[50];
      REAL(ret)[109] = REAL(theta)[60];
      REAL(ret)[119] = REAL(theta)[60];
    }
    else if (theta_n == 52){
      REAL(ret)[75] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[86] = REAL(theta)[34];
      REAL(ret)[97] = REAL(theta)[42];
      REAL(ret)[105] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[106] = REAL(theta)[34];
      REAL(ret)[107] = REAL(theta)[42];
      REAL(ret)[108] = 2 * REAL(theta)[51];
      REAL(ret)[109] = REAL(theta)[61];
      REAL(ret)[119] = REAL(theta)[61];
    }
    else if (theta_n == 53){
      REAL(ret)[86] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[97] = REAL(theta)[43];
      REAL(ret)[106] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[107] = REAL(theta)[43];
      REAL(ret)[108] = 2 * REAL(theta)[52];
      REAL(ret)[109] = REAL(theta)[62];
      REAL(ret)[119] = REAL(theta)[62];
    }
    else if (theta_n == 54){
      REAL(ret)[97] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[107] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[108] = 2 * REAL(theta)[53];
      REAL(ret)[109] = REAL(theta)[63];
      REAL(ret)[119] = REAL(theta)[63];
    }
    else if (theta_n == 55){
      REAL(ret)[108] = 4 * R_pow_di(REAL(theta)[54], 3);
      REAL(ret)[109] = 2 * REAL(theta)[54] * REAL(theta)[64];
      REAL(ret)[119] = 2 * REAL(theta)[54] * REAL(theta)[64];
    }
    else if (theta_n == 56){
      REAL(ret)[10] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[21] = REAL(theta)[1];
      REAL(ret)[32] = REAL(theta)[3];
      REAL(ret)[43] = REAL(theta)[6];
      REAL(ret)[54] = REAL(theta)[10];
      REAL(ret)[65] = REAL(theta)[15];
      REAL(ret)[76] = REAL(theta)[21];
      REAL(ret)[87] = REAL(theta)[28];
      REAL(ret)[98] = REAL(theta)[36];
      REAL(ret)[109] = REAL(theta)[45];
      REAL(ret)[110] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[21] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[32] = REAL(theta)[4];
      REAL(ret)[43] = REAL(theta)[7];
      REAL(ret)[54] = REAL(theta)[11];
      REAL(ret)[65] = REAL(theta)[16];
      REAL(ret)[76] = REAL(theta)[22];
      REAL(ret)[87] = REAL(theta)[29];
      REAL(ret)[98] = REAL(theta)[37];
      REAL(ret)[109] = REAL(theta)[46];
      REAL(ret)[111] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[32] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[43] = REAL(theta)[8];
      REAL(ret)[54] = REAL(theta)[12];
      REAL(ret)[65] = REAL(theta)[17];
      REAL(ret)[76] = REAL(theta)[23];
      REAL(ret)[87] = REAL(theta)[30];
      REAL(ret)[98] = REAL(theta)[38];
      REAL(ret)[109] = REAL(theta)[47];
      REAL(ret)[112] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[43] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[54] = REAL(theta)[13];
      REAL(ret)[65] = REAL(theta)[18];
      REAL(ret)[76] = REAL(theta)[24];
      REAL(ret)[87] = REAL(theta)[31];
      REAL(ret)[98] = REAL(theta)[39];
      REAL(ret)[109] = REAL(theta)[48];
      REAL(ret)[113] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[114] = REAL(theta)[13];
      REAL(ret)[115] = REAL(theta)[18];
      REAL(ret)[116] = REAL(theta)[24];
      REAL(ret)[117] = REAL(theta)[31];
      REAL(ret)[118] = REAL(theta)[39];
      REAL(ret)[119] = REAL(theta)[48];
      REAL(ret)[120] = 2 * REAL(theta)[58];
    }
    else if (theta_n == 60){
      REAL(ret)[54] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[65] = REAL(theta)[19];
      REAL(ret)[76] = REAL(theta)[25];
      REAL(ret)[87] = REAL(theta)[32];
      REAL(ret)[98] = REAL(theta)[40];
      REAL(ret)[109] = REAL(theta)[49];
      REAL(ret)[114] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[115] = REAL(theta)[19];
      REAL(ret)[116] = REAL(theta)[25];
      REAL(ret)[117] = REAL(theta)[32];
      REAL(ret)[118] = REAL(theta)[40];
      REAL(ret)[119] = REAL(theta)[49];
      REAL(ret)[120] = 2 * REAL(theta)[59];
    }
    else if (theta_n == 61){
      REAL(ret)[65] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[76] = REAL(theta)[26];
      REAL(ret)[87] = REAL(theta)[33];
      REAL(ret)[98] = REAL(theta)[41];
      REAL(ret)[109] = REAL(theta)[50];
      REAL(ret)[115] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[116] = REAL(theta)[26];
      REAL(ret)[117] = REAL(theta)[33];
      REAL(ret)[118] = REAL(theta)[41];
      REAL(ret)[119] = REAL(theta)[50];
      REAL(ret)[120] = 2 * REAL(theta)[60];
    }
    else if (theta_n == 62){
      REAL(ret)[76] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[87] = REAL(theta)[34];
      REAL(ret)[98] = REAL(theta)[42];
      REAL(ret)[109] = REAL(theta)[51];
      REAL(ret)[116] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[117] = REAL(theta)[34];
      REAL(ret)[118] = REAL(theta)[42];
      REAL(ret)[119] = REAL(theta)[51];
      REAL(ret)[120] = 2 * REAL(theta)[61];
    }
    else if (theta_n == 63){
      REAL(ret)[87] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[98] = REAL(theta)[43];
      REAL(ret)[109] = REAL(theta)[52];
      REAL(ret)[117] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[118] = REAL(theta)[43];
      REAL(ret)[119] = REAL(theta)[52];
      REAL(ret)[120] = 2 * REAL(theta)[62];
    }
    else if (theta_n == 64){
      REAL(ret)[98] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[109] = REAL(theta)[53];
      REAL(ret)[118] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[119] = REAL(theta)[53];
      REAL(ret)[120] = 2 * REAL(theta)[63];
    }
    else if (theta_n == 65){
      REAL(ret)[109] = R_pow_di(REAL(theta)[54], 2);
      REAL(ret)[119] = R_pow_di(REAL(theta)[54], 2);
      REAL(ret)[120] = 2 * REAL(theta)[64];
    }
    else if (theta_n == 66){
      REAL(ret)[120] = 4 * R_pow_di(REAL(theta)[65], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 11));for(int i = 0; i < 11; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 4 * R_pow_di(REAL(theta)[5], 3);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 4 * R_pow_di(REAL(theta)[9], 3);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 4 * R_pow_di(REAL(theta)[14], 3);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 4 * R_pow_di(REAL(theta)[20], 3);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 4 * R_pow_di(REAL(theta)[27], 3);
    }
    else if (theta_n == -38){
      REAL(ret)[7] = 4 * R_pow_di(REAL(theta)[35], 3);
    }
    else if (theta_n == -47){
      REAL(ret)[8] = 4 * R_pow_di(REAL(theta)[44], 3);
    }
    else if (theta_n == -57){
      REAL(ret)[9] = 4 * R_pow_di(REAL(theta)[54], 3);
    }
    else if (theta_n == -68){
      REAL(ret)[10] = 4 * R_pow_di(REAL(theta)[65], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}
else if (dm == 12){
  int theta_n = INTEGER(tn)[0];
  if (theta_n== NA_INTEGER){
    SEXP ret=  PROTECT(allocVector(INTSXP,78));
    INTEGER(ret)[0]=2;
    INTEGER(ret)[1]=5;
    INTEGER(ret)[2]=2;
    INTEGER(ret)[3]=5;
    INTEGER(ret)[4]=5;
    INTEGER(ret)[5]=2;
    INTEGER(ret)[6]=5;
    INTEGER(ret)[7]=5;
    INTEGER(ret)[8]=5;
    INTEGER(ret)[9]=2;
    INTEGER(ret)[10]=5;
    INTEGER(ret)[11]=5;
    INTEGER(ret)[12]=5;
    INTEGER(ret)[13]=5;
    INTEGER(ret)[14]=2;
    INTEGER(ret)[15]=5;
    INTEGER(ret)[16]=5;
    INTEGER(ret)[17]=5;
    INTEGER(ret)[18]=5;
    INTEGER(ret)[19]=5;
    INTEGER(ret)[20]=2;
    INTEGER(ret)[21]=5;
    INTEGER(ret)[22]=5;
    INTEGER(ret)[23]=5;
    INTEGER(ret)[24]=5;
    INTEGER(ret)[25]=5;
    INTEGER(ret)[26]=5;
    INTEGER(ret)[27]=2;
    INTEGER(ret)[28]=5;
    INTEGER(ret)[29]=5;
    INTEGER(ret)[30]=5;
    INTEGER(ret)[31]=5;
    INTEGER(ret)[32]=5;
    INTEGER(ret)[33]=5;
    INTEGER(ret)[34]=5;
    INTEGER(ret)[35]=2;
    INTEGER(ret)[36]=5;
    INTEGER(ret)[37]=5;
    INTEGER(ret)[38]=5;
    INTEGER(ret)[39]=5;
    INTEGER(ret)[40]=5;
    INTEGER(ret)[41]=5;
    INTEGER(ret)[42]=5;
    INTEGER(ret)[43]=5;
    INTEGER(ret)[44]=2;
    INTEGER(ret)[45]=5;
    INTEGER(ret)[46]=5;
    INTEGER(ret)[47]=5;
    INTEGER(ret)[48]=5;
    INTEGER(ret)[49]=5;
    INTEGER(ret)[50]=5;
    INTEGER(ret)[51]=5;
    INTEGER(ret)[52]=5;
    INTEGER(ret)[53]=5;
    INTEGER(ret)[54]=2;
    INTEGER(ret)[55]=5;
    INTEGER(ret)[56]=5;
    INTEGER(ret)[57]=5;
    INTEGER(ret)[58]=5;
    INTEGER(ret)[59]=5;
    INTEGER(ret)[60]=5;
    INTEGER(ret)[61]=5;
    INTEGER(ret)[62]=5;
    INTEGER(ret)[63]=5;
    INTEGER(ret)[64]=5;
    INTEGER(ret)[65]=2;
    INTEGER(ret)[66]=5;
    INTEGER(ret)[67]=5;
    INTEGER(ret)[68]=5;
    INTEGER(ret)[69]=5;
    INTEGER(ret)[70]=5;
    INTEGER(ret)[71]=5;
    INTEGER(ret)[72]=5;
    INTEGER(ret)[73]=5;
    INTEGER(ret)[74]=5;
    INTEGER(ret)[75]=5;
    INTEGER(ret)[76]=5;
    INTEGER(ret)[77]=2;
    UNPROTECT(1);
    return(ret);  
}

if (theta_n == -2){
    SEXP ret = PROTECT(allocVector(INTSXP, 1));
    INTEGER(ret)[0] = 78;
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
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[12] = REAL(theta)[1];
      REAL(ret)[13] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[24] = REAL(theta)[3];
      REAL(ret)[25] = REAL(theta)[4];
      REAL(ret)[26] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[36] = REAL(theta)[6];
      REAL(ret)[37] = REAL(theta)[7];
      REAL(ret)[38] = REAL(theta)[8];
      REAL(ret)[39] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[48] = REAL(theta)[10];
      REAL(ret)[49] = REAL(theta)[11];
      REAL(ret)[50] = REAL(theta)[12];
      REAL(ret)[51] = REAL(theta)[13];
      REAL(ret)[52] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[60] = REAL(theta)[15];
      REAL(ret)[61] = REAL(theta)[16];
      REAL(ret)[62] = REAL(theta)[17];
      REAL(ret)[63] = REAL(theta)[18];
      REAL(ret)[64] = REAL(theta)[19];
      REAL(ret)[65] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[72] = REAL(theta)[21];
      REAL(ret)[73] = REAL(theta)[22];
      REAL(ret)[74] = REAL(theta)[23];
      REAL(ret)[75] = REAL(theta)[24];
      REAL(ret)[76] = REAL(theta)[25];
      REAL(ret)[77] = REAL(theta)[26];
      REAL(ret)[78] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[84] = REAL(theta)[28];
      REAL(ret)[85] = REAL(theta)[29];
      REAL(ret)[86] = REAL(theta)[30];
      REAL(ret)[87] = REAL(theta)[31];
      REAL(ret)[88] = REAL(theta)[32];
      REAL(ret)[89] = REAL(theta)[33];
      REAL(ret)[90] = REAL(theta)[34];
      REAL(ret)[91] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[96] = REAL(theta)[36];
      REAL(ret)[97] = REAL(theta)[37];
      REAL(ret)[98] = REAL(theta)[38];
      REAL(ret)[99] = REAL(theta)[39];
      REAL(ret)[100] = REAL(theta)[40];
      REAL(ret)[101] = REAL(theta)[41];
      REAL(ret)[102] = REAL(theta)[42];
      REAL(ret)[103] = REAL(theta)[43];
      REAL(ret)[104] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[108] = REAL(theta)[45];
      REAL(ret)[109] = REAL(theta)[46];
      REAL(ret)[110] = REAL(theta)[47];
      REAL(ret)[111] = REAL(theta)[48];
      REAL(ret)[112] = REAL(theta)[49];
      REAL(ret)[113] = REAL(theta)[50];
      REAL(ret)[114] = REAL(theta)[51];
      REAL(ret)[115] = REAL(theta)[52];
      REAL(ret)[116] = REAL(theta)[53];
      REAL(ret)[117] = R_pow_di(REAL(theta)[54], 2);
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
      REAL(ret)[130] = R_pow_di(REAL(theta)[65], 2);
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
      REAL(ret)[143] = R_pow_di(REAL(theta)[77], 2);
    }
    else if (theta_n == -1){
      REAL(ret)[0] = R_pow_di(REAL(theta)[0], 4);
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[28];
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[36];
      REAL(ret)[9] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[45];
      REAL(ret)[10] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[55];
      REAL(ret)[11] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[66];
      REAL(ret)[12] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[1];
      REAL(ret)[13] = R_pow_di(REAL(theta)[1], 2) + R_pow_di(REAL(theta)[2], 4);
      REAL(ret)[14] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[15] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[16] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[17] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[18] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[19] = REAL(theta)[1] * REAL(theta)[28] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[29];
      REAL(ret)[20] = REAL(theta)[1] * REAL(theta)[36] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[37];
      REAL(ret)[21] = REAL(theta)[1] * REAL(theta)[45] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[46];
      REAL(ret)[22] = REAL(theta)[1] * REAL(theta)[55] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[56];
      REAL(ret)[23] = REAL(theta)[1] * REAL(theta)[66] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[67];
      REAL(ret)[24] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[3];
      REAL(ret)[25] = REAL(theta)[1] * REAL(theta)[3] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[4];
      REAL(ret)[26] = R_pow_di(REAL(theta)[3], 2) + R_pow_di(REAL(theta)[4], 2) + R_pow_di(REAL(theta)[5], 4);
      REAL(ret)[27] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[28] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[29] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[30] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[31] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[32] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[33] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[34] = REAL(theta)[3] * REAL(theta)[55] + REAL(theta)[4] * REAL(theta)[56] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[57];
      REAL(ret)[35] = REAL(theta)[3] * REAL(theta)[66] + REAL(theta)[4] * REAL(theta)[67] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[68];
      REAL(ret)[36] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[6];
      REAL(ret)[37] = REAL(theta)[1] * REAL(theta)[6] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[7];
      REAL(ret)[38] = REAL(theta)[3] * REAL(theta)[6] + REAL(theta)[4] * REAL(theta)[7] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[8];
      REAL(ret)[39] = R_pow_di(REAL(theta)[6], 2) + R_pow_di(REAL(theta)[7], 2) + R_pow_di(REAL(theta)[8], 2) + R_pow_di(REAL(theta)[9], 4);
      REAL(ret)[40] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[41] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[42] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[43] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[44] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[45] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[46] = REAL(theta)[55] * REAL(theta)[6] + REAL(theta)[56] * REAL(theta)[7] + REAL(theta)[57] * REAL(theta)[8] + REAL(theta)[58] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[47] = REAL(theta)[6] * REAL(theta)[66] + REAL(theta)[67] * REAL(theta)[7] + REAL(theta)[68] * REAL(theta)[8] + REAL(theta)[69] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[48] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[10];
      REAL(ret)[49] = REAL(theta)[1] * REAL(theta)[10] + REAL(theta)[11] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[50] = REAL(theta)[10] * REAL(theta)[3] + REAL(theta)[11] * REAL(theta)[4] + REAL(theta)[12] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[51] = REAL(theta)[10] * REAL(theta)[6] + REAL(theta)[11] * REAL(theta)[7] + REAL(theta)[12] * REAL(theta)[8] + REAL(theta)[13] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[52] = R_pow_di(REAL(theta)[10], 2) + R_pow_di(REAL(theta)[11], 2) + R_pow_di(REAL(theta)[12], 2) + R_pow_di(REAL(theta)[13], 2) + R_pow_di(REAL(theta)[14], 4);
      REAL(ret)[53] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[54] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[55] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[32];
      REAL(ret)[56] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[40];
      REAL(ret)[57] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[49];
      REAL(ret)[58] = REAL(theta)[10] * REAL(theta)[55] + REAL(theta)[11] * REAL(theta)[56] + REAL(theta)[12] * REAL(theta)[57] + REAL(theta)[13] * REAL(theta)[58] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[59];
      REAL(ret)[59] = REAL(theta)[10] * REAL(theta)[66] + REAL(theta)[11] * REAL(theta)[67] + REAL(theta)[12] * REAL(theta)[68] + REAL(theta)[13] * REAL(theta)[69] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[70];
      REAL(ret)[60] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[15];
      REAL(ret)[61] = REAL(theta)[1] * REAL(theta)[15] + REAL(theta)[16] * R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[62] = REAL(theta)[15] * REAL(theta)[3] + REAL(theta)[16] * REAL(theta)[4] + REAL(theta)[17] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[63] = REAL(theta)[15] * REAL(theta)[6] + REAL(theta)[16] * REAL(theta)[7] + REAL(theta)[17] * REAL(theta)[8] + REAL(theta)[18] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[64] = REAL(theta)[10] * REAL(theta)[15] + REAL(theta)[11] * REAL(theta)[16] + REAL(theta)[12] * REAL(theta)[17] + REAL(theta)[13] * REAL(theta)[18] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[19];
      REAL(ret)[65] = R_pow_di(REAL(theta)[15], 2) + R_pow_di(REAL(theta)[16], 2) + R_pow_di(REAL(theta)[17], 2) + R_pow_di(REAL(theta)[18], 2) + R_pow_di(REAL(theta)[19], 2) + R_pow_di(REAL(theta)[20], 4);
      REAL(ret)[66] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[67] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[33];
      REAL(ret)[68] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[41];
      REAL(ret)[69] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[50];
      REAL(ret)[70] = REAL(theta)[15] * REAL(theta)[55] + REAL(theta)[16] * REAL(theta)[56] + REAL(theta)[17] * REAL(theta)[57] + REAL(theta)[18] * REAL(theta)[58] + REAL(theta)[19] * REAL(theta)[59] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[60];
      REAL(ret)[71] = REAL(theta)[15] * REAL(theta)[66] + REAL(theta)[16] * REAL(theta)[67] + REAL(theta)[17] * REAL(theta)[68] + REAL(theta)[18] * REAL(theta)[69] + REAL(theta)[19] * REAL(theta)[70] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[71];
      REAL(ret)[72] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[21];
      REAL(ret)[73] = REAL(theta)[1] * REAL(theta)[21] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[22];
      REAL(ret)[74] = REAL(theta)[21] * REAL(theta)[3] + REAL(theta)[22] * REAL(theta)[4] + REAL(theta)[23] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[75] = REAL(theta)[21] * REAL(theta)[6] + REAL(theta)[22] * REAL(theta)[7] + REAL(theta)[23] * REAL(theta)[8] + REAL(theta)[24] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[76] = REAL(theta)[10] * REAL(theta)[21] + REAL(theta)[11] * REAL(theta)[22] + REAL(theta)[12] * REAL(theta)[23] + REAL(theta)[13] * REAL(theta)[24] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[25];
      REAL(ret)[77] = REAL(theta)[15] * REAL(theta)[21] + REAL(theta)[16] * REAL(theta)[22] + REAL(theta)[17] * REAL(theta)[23] + REAL(theta)[18] * REAL(theta)[24] + REAL(theta)[19] * REAL(theta)[25] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[26];
      REAL(ret)[78] = R_pow_di(REAL(theta)[21], 2) + R_pow_di(REAL(theta)[22], 2) + R_pow_di(REAL(theta)[23], 2) + R_pow_di(REAL(theta)[24], 2) + R_pow_di(REAL(theta)[25], 2) + R_pow_di(REAL(theta)[26], 2) + R_pow_di(REAL(theta)[27], 4);
      REAL(ret)[79] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[34];
      REAL(ret)[80] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[42];
      REAL(ret)[81] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[51];
      REAL(ret)[82] = REAL(theta)[21] * REAL(theta)[55] + REAL(theta)[22] * REAL(theta)[56] + REAL(theta)[23] * REAL(theta)[57] + REAL(theta)[24] * REAL(theta)[58] + REAL(theta)[25] * REAL(theta)[59] + REAL(theta)[26] * REAL(theta)[60] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[61];
      REAL(ret)[83] = REAL(theta)[21] * REAL(theta)[66] + REAL(theta)[22] * REAL(theta)[67] + REAL(theta)[23] * REAL(theta)[68] + REAL(theta)[24] * REAL(theta)[69] + REAL(theta)[25] * REAL(theta)[70] + REAL(theta)[26] * REAL(theta)[71] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[72];
      REAL(ret)[84] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[28];
      REAL(ret)[85] = REAL(theta)[1] * REAL(theta)[28] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[29];
      REAL(ret)[86] = REAL(theta)[28] * REAL(theta)[3] + REAL(theta)[29] * REAL(theta)[4] + REAL(theta)[30] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[87] = REAL(theta)[28] * REAL(theta)[6] + REAL(theta)[29] * REAL(theta)[7] + REAL(theta)[30] * REAL(theta)[8] + REAL(theta)[31] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[88] = REAL(theta)[10] * REAL(theta)[28] + REAL(theta)[11] * REAL(theta)[29] + REAL(theta)[12] * REAL(theta)[30] + REAL(theta)[13] * REAL(theta)[31] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[32];
      REAL(ret)[89] = REAL(theta)[15] * REAL(theta)[28] + REAL(theta)[16] * REAL(theta)[29] + REAL(theta)[17] * REAL(theta)[30] + REAL(theta)[18] * REAL(theta)[31] + REAL(theta)[19] * REAL(theta)[32] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[33];
      REAL(ret)[90] = REAL(theta)[21] * REAL(theta)[28] + REAL(theta)[22] * REAL(theta)[29] + REAL(theta)[23] * REAL(theta)[30] + REAL(theta)[24] * REAL(theta)[31] + REAL(theta)[25] * REAL(theta)[32] + REAL(theta)[26] * REAL(theta)[33] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[34];
      REAL(ret)[91] = R_pow_di(REAL(theta)[28], 2) + R_pow_di(REAL(theta)[29], 2) + R_pow_di(REAL(theta)[30], 2) + R_pow_di(REAL(theta)[31], 2) + R_pow_di(REAL(theta)[32], 2) + R_pow_di(REAL(theta)[33], 2) + R_pow_di(REAL(theta)[34], 2) + R_pow_di(REAL(theta)[35], 4);
      REAL(ret)[92] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[43];
      REAL(ret)[93] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[52];
      REAL(ret)[94] = REAL(theta)[28] * REAL(theta)[55] + REAL(theta)[29] * REAL(theta)[56] + REAL(theta)[30] * REAL(theta)[57] + REAL(theta)[31] * REAL(theta)[58] + REAL(theta)[32] * REAL(theta)[59] + REAL(theta)[33] * REAL(theta)[60] + REAL(theta)[34] * REAL(theta)[61] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[62];
      REAL(ret)[95] = REAL(theta)[28] * REAL(theta)[66] + REAL(theta)[29] * REAL(theta)[67] + REAL(theta)[30] * REAL(theta)[68] + REAL(theta)[31] * REAL(theta)[69] + REAL(theta)[32] * REAL(theta)[70] + REAL(theta)[33] * REAL(theta)[71] + REAL(theta)[34] * REAL(theta)[72] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[73];
      REAL(ret)[96] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[36];
      REAL(ret)[97] = REAL(theta)[1] * REAL(theta)[36] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[37];
      REAL(ret)[98] = REAL(theta)[3] * REAL(theta)[36] + REAL(theta)[37] * REAL(theta)[4] + REAL(theta)[38] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[99] = REAL(theta)[36] * REAL(theta)[6] + REAL(theta)[37] * REAL(theta)[7] + REAL(theta)[38] * REAL(theta)[8] + REAL(theta)[39] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[100] = REAL(theta)[10] * REAL(theta)[36] + REAL(theta)[11] * REAL(theta)[37] + REAL(theta)[12] * REAL(theta)[38] + REAL(theta)[13] * REAL(theta)[39] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[40];
      REAL(ret)[101] = REAL(theta)[15] * REAL(theta)[36] + REAL(theta)[16] * REAL(theta)[37] + REAL(theta)[17] * REAL(theta)[38] + REAL(theta)[18] * REAL(theta)[39] + REAL(theta)[19] * REAL(theta)[40] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[41];
      REAL(ret)[102] = REAL(theta)[21] * REAL(theta)[36] + REAL(theta)[22] * REAL(theta)[37] + REAL(theta)[23] * REAL(theta)[38] + REAL(theta)[24] * REAL(theta)[39] + REAL(theta)[25] * REAL(theta)[40] + REAL(theta)[26] * REAL(theta)[41] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[42];
      REAL(ret)[103] = REAL(theta)[28] * REAL(theta)[36] + REAL(theta)[29] * REAL(theta)[37] + REAL(theta)[30] * REAL(theta)[38] + REAL(theta)[31] * REAL(theta)[39] + REAL(theta)[32] * REAL(theta)[40] + REAL(theta)[33] * REAL(theta)[41] + REAL(theta)[34] * REAL(theta)[42] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[43];
      REAL(ret)[104] = R_pow_di(REAL(theta)[36], 2) + R_pow_di(REAL(theta)[37], 2) + R_pow_di(REAL(theta)[38], 2) + R_pow_di(REAL(theta)[39], 2) + R_pow_di(REAL(theta)[40], 2) + R_pow_di(REAL(theta)[41], 2) + R_pow_di(REAL(theta)[42], 2) + R_pow_di(REAL(theta)[43], 2) + R_pow_di(REAL(theta)[44], 4);
      REAL(ret)[105] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[53];
      REAL(ret)[106] = REAL(theta)[36] * REAL(theta)[55] + REAL(theta)[37] * REAL(theta)[56] + REAL(theta)[38] * REAL(theta)[57] + REAL(theta)[39] * REAL(theta)[58] + REAL(theta)[40] * REAL(theta)[59] + REAL(theta)[41] * REAL(theta)[60] + REAL(theta)[42] * REAL(theta)[61] + REAL(theta)[43] * REAL(theta)[62] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[63];
      REAL(ret)[107] = REAL(theta)[36] * REAL(theta)[66] + REAL(theta)[37] * REAL(theta)[67] + REAL(theta)[38] * REAL(theta)[68] + REAL(theta)[39] * REAL(theta)[69] + REAL(theta)[40] * REAL(theta)[70] + REAL(theta)[41] * REAL(theta)[71] + REAL(theta)[42] * REAL(theta)[72] + REAL(theta)[43] * REAL(theta)[73] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[74];
      REAL(ret)[108] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[45];
      REAL(ret)[109] = REAL(theta)[1] * REAL(theta)[45] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[46];
      REAL(ret)[110] = REAL(theta)[3] * REAL(theta)[45] + REAL(theta)[4] * REAL(theta)[46] + REAL(theta)[47] * R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[111] = REAL(theta)[45] * REAL(theta)[6] + REAL(theta)[46] * REAL(theta)[7] + REAL(theta)[47] * REAL(theta)[8] + REAL(theta)[48] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[112] = REAL(theta)[10] * REAL(theta)[45] + REAL(theta)[11] * REAL(theta)[46] + REAL(theta)[12] * REAL(theta)[47] + REAL(theta)[13] * REAL(theta)[48] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[49];
      REAL(ret)[113] = REAL(theta)[15] * REAL(theta)[45] + REAL(theta)[16] * REAL(theta)[46] + REAL(theta)[17] * REAL(theta)[47] + REAL(theta)[18] * REAL(theta)[48] + REAL(theta)[19] * REAL(theta)[49] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[50];
      REAL(ret)[114] = REAL(theta)[21] * REAL(theta)[45] + REAL(theta)[22] * REAL(theta)[46] + REAL(theta)[23] * REAL(theta)[47] + REAL(theta)[24] * REAL(theta)[48] + REAL(theta)[25] * REAL(theta)[49] + REAL(theta)[26] * REAL(theta)[50] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[51];
      REAL(ret)[115] = REAL(theta)[28] * REAL(theta)[45] + REAL(theta)[29] * REAL(theta)[46] + REAL(theta)[30] * REAL(theta)[47] + REAL(theta)[31] * REAL(theta)[48] + REAL(theta)[32] * REAL(theta)[49] + REAL(theta)[33] * REAL(theta)[50] + REAL(theta)[34] * REAL(theta)[51] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[52];
      REAL(ret)[116] = REAL(theta)[36] * REAL(theta)[45] + REAL(theta)[37] * REAL(theta)[46] + REAL(theta)[38] * REAL(theta)[47] + REAL(theta)[39] * REAL(theta)[48] + REAL(theta)[40] * REAL(theta)[49] + REAL(theta)[41] * REAL(theta)[50] + REAL(theta)[42] * REAL(theta)[51] + REAL(theta)[43] * REAL(theta)[52] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[53];
      REAL(ret)[117] = R_pow_di(REAL(theta)[45], 2) + R_pow_di(REAL(theta)[46], 2) + R_pow_di(REAL(theta)[47], 2) + R_pow_di(REAL(theta)[48], 2) + R_pow_di(REAL(theta)[49], 2) + R_pow_di(REAL(theta)[50], 2) + R_pow_di(REAL(theta)[51], 2) + R_pow_di(REAL(theta)[52], 2) + R_pow_di(REAL(theta)[53], 2) + R_pow_di(REAL(theta)[54], 4);
      REAL(ret)[118] = REAL(theta)[45] * REAL(theta)[55] + REAL(theta)[46] * REAL(theta)[56] + REAL(theta)[47] * REAL(theta)[57] + REAL(theta)[48] * REAL(theta)[58] + REAL(theta)[49] * REAL(theta)[59] + REAL(theta)[50] * REAL(theta)[60] + REAL(theta)[51] * REAL(theta)[61] + REAL(theta)[52] * REAL(theta)[62] + REAL(theta)[53] * REAL(theta)[63] + R_pow_di(REAL(theta)[54], 2) * REAL(theta)[64];
      REAL(ret)[119] = REAL(theta)[45] * REAL(theta)[66] + REAL(theta)[46] * REAL(theta)[67] + REAL(theta)[47] * REAL(theta)[68] + REAL(theta)[48] * REAL(theta)[69] + REAL(theta)[49] * REAL(theta)[70] + REAL(theta)[50] * REAL(theta)[71] + REAL(theta)[51] * REAL(theta)[72] + REAL(theta)[52] * REAL(theta)[73] + REAL(theta)[53] * REAL(theta)[74] + R_pow_di(REAL(theta)[54], 2) * REAL(theta)[75];
      REAL(ret)[120] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[55];
      REAL(ret)[121] = REAL(theta)[1] * REAL(theta)[55] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[56];
      REAL(ret)[122] = REAL(theta)[3] * REAL(theta)[55] + REAL(theta)[4] * REAL(theta)[56] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[57];
      REAL(ret)[123] = REAL(theta)[55] * REAL(theta)[6] + REAL(theta)[56] * REAL(theta)[7] + REAL(theta)[57] * REAL(theta)[8] + REAL(theta)[58] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[124] = REAL(theta)[10] * REAL(theta)[55] + REAL(theta)[11] * REAL(theta)[56] + REAL(theta)[12] * REAL(theta)[57] + REAL(theta)[13] * REAL(theta)[58] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[59];
      REAL(ret)[125] = REAL(theta)[15] * REAL(theta)[55] + REAL(theta)[16] * REAL(theta)[56] + REAL(theta)[17] * REAL(theta)[57] + REAL(theta)[18] * REAL(theta)[58] + REAL(theta)[19] * REAL(theta)[59] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[60];
      REAL(ret)[126] = REAL(theta)[21] * REAL(theta)[55] + REAL(theta)[22] * REAL(theta)[56] + REAL(theta)[23] * REAL(theta)[57] + REAL(theta)[24] * REAL(theta)[58] + REAL(theta)[25] * REAL(theta)[59] + REAL(theta)[26] * REAL(theta)[60] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[61];
      REAL(ret)[127] = REAL(theta)[28] * REAL(theta)[55] + REAL(theta)[29] * REAL(theta)[56] + REAL(theta)[30] * REAL(theta)[57] + REAL(theta)[31] * REAL(theta)[58] + REAL(theta)[32] * REAL(theta)[59] + REAL(theta)[33] * REAL(theta)[60] + REAL(theta)[34] * REAL(theta)[61] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[62];
      REAL(ret)[128] = REAL(theta)[36] * REAL(theta)[55] + REAL(theta)[37] * REAL(theta)[56] + REAL(theta)[38] * REAL(theta)[57] + REAL(theta)[39] * REAL(theta)[58] + REAL(theta)[40] * REAL(theta)[59] + REAL(theta)[41] * REAL(theta)[60] + REAL(theta)[42] * REAL(theta)[61] + REAL(theta)[43] * REAL(theta)[62] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[63];
      REAL(ret)[129] = REAL(theta)[45] * REAL(theta)[55] + REAL(theta)[46] * REAL(theta)[56] + REAL(theta)[47] * REAL(theta)[57] + REAL(theta)[48] * REAL(theta)[58] + REAL(theta)[49] * REAL(theta)[59] + REAL(theta)[50] * REAL(theta)[60] + REAL(theta)[51] * REAL(theta)[61] + REAL(theta)[52] * REAL(theta)[62] + REAL(theta)[53] * REAL(theta)[63] + R_pow_di(REAL(theta)[54], 2) * REAL(theta)[64];
      REAL(ret)[130] = R_pow_di(REAL(theta)[55], 2) + R_pow_di(REAL(theta)[56], 2) + R_pow_di(REAL(theta)[57], 2) + R_pow_di(REAL(theta)[58], 2) + R_pow_di(REAL(theta)[59], 2) + R_pow_di(REAL(theta)[60], 2) + R_pow_di(REAL(theta)[61], 2) + R_pow_di(REAL(theta)[62], 2) + R_pow_di(REAL(theta)[63], 2) + R_pow_di(REAL(theta)[64], 2) + R_pow_di(REAL(theta)[65], 4);
      REAL(ret)[131] = REAL(theta)[55] * REAL(theta)[66] + REAL(theta)[56] * REAL(theta)[67] + REAL(theta)[57] * REAL(theta)[68] + REAL(theta)[58] * REAL(theta)[69] + REAL(theta)[59] * REAL(theta)[70] + REAL(theta)[60] * REAL(theta)[71] + REAL(theta)[61] * REAL(theta)[72] + REAL(theta)[62] * REAL(theta)[73] + REAL(theta)[63] * REAL(theta)[74] + REAL(theta)[64] * REAL(theta)[75] + R_pow_di(REAL(theta)[65], 2) * REAL(theta)[76];
      REAL(ret)[132] = R_pow_di(REAL(theta)[0], 2) * REAL(theta)[66];
      REAL(ret)[133] = REAL(theta)[1] * REAL(theta)[66] + R_pow_di(REAL(theta)[2], 2) * REAL(theta)[67];
      REAL(ret)[134] = REAL(theta)[3] * REAL(theta)[66] + REAL(theta)[4] * REAL(theta)[67] + R_pow_di(REAL(theta)[5], 2) * REAL(theta)[68];
      REAL(ret)[135] = REAL(theta)[6] * REAL(theta)[66] + REAL(theta)[67] * REAL(theta)[7] + REAL(theta)[68] * REAL(theta)[8] + REAL(theta)[69] * R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[136] = REAL(theta)[10] * REAL(theta)[66] + REAL(theta)[11] * REAL(theta)[67] + REAL(theta)[12] * REAL(theta)[68] + REAL(theta)[13] * REAL(theta)[69] + R_pow_di(REAL(theta)[14], 2) * REAL(theta)[70];
      REAL(ret)[137] = REAL(theta)[15] * REAL(theta)[66] + REAL(theta)[16] * REAL(theta)[67] + REAL(theta)[17] * REAL(theta)[68] + REAL(theta)[18] * REAL(theta)[69] + REAL(theta)[19] * REAL(theta)[70] + R_pow_di(REAL(theta)[20], 2) * REAL(theta)[71];
      REAL(ret)[138] = REAL(theta)[21] * REAL(theta)[66] + REAL(theta)[22] * REAL(theta)[67] + REAL(theta)[23] * REAL(theta)[68] + REAL(theta)[24] * REAL(theta)[69] + REAL(theta)[25] * REAL(theta)[70] + REAL(theta)[26] * REAL(theta)[71] + R_pow_di(REAL(theta)[27], 2) * REAL(theta)[72];
      REAL(ret)[139] = REAL(theta)[28] * REAL(theta)[66] + REAL(theta)[29] * REAL(theta)[67] + REAL(theta)[30] * REAL(theta)[68] + REAL(theta)[31] * REAL(theta)[69] + REAL(theta)[32] * REAL(theta)[70] + REAL(theta)[33] * REAL(theta)[71] + REAL(theta)[34] * REAL(theta)[72] + R_pow_di(REAL(theta)[35], 2) * REAL(theta)[73];
      REAL(ret)[140] = REAL(theta)[36] * REAL(theta)[66] + REAL(theta)[37] * REAL(theta)[67] + REAL(theta)[38] * REAL(theta)[68] + REAL(theta)[39] * REAL(theta)[69] + REAL(theta)[40] * REAL(theta)[70] + REAL(theta)[41] * REAL(theta)[71] + REAL(theta)[42] * REAL(theta)[72] + REAL(theta)[43] * REAL(theta)[73] + R_pow_di(REAL(theta)[44], 2) * REAL(theta)[74];
      REAL(ret)[141] = REAL(theta)[45] * REAL(theta)[66] + REAL(theta)[46] * REAL(theta)[67] + REAL(theta)[47] * REAL(theta)[68] + REAL(theta)[48] * REAL(theta)[69] + REAL(theta)[49] * REAL(theta)[70] + REAL(theta)[50] * REAL(theta)[71] + REAL(theta)[51] * REAL(theta)[72] + REAL(theta)[52] * REAL(theta)[73] + REAL(theta)[53] * REAL(theta)[74] + R_pow_di(REAL(theta)[54], 2) * REAL(theta)[75];
      REAL(ret)[142] = REAL(theta)[55] * REAL(theta)[66] + REAL(theta)[56] * REAL(theta)[67] + REAL(theta)[57] * REAL(theta)[68] + REAL(theta)[58] * REAL(theta)[69] + REAL(theta)[59] * REAL(theta)[70] + REAL(theta)[60] * REAL(theta)[71] + REAL(theta)[61] * REAL(theta)[72] + REAL(theta)[62] * REAL(theta)[73] + REAL(theta)[63] * REAL(theta)[74] + REAL(theta)[64] * REAL(theta)[75] + R_pow_di(REAL(theta)[65], 2) * REAL(theta)[76];
      REAL(ret)[143] = R_pow_di(REAL(theta)[66], 2) + R_pow_di(REAL(theta)[67], 2) + R_pow_di(REAL(theta)[68], 2) + R_pow_di(REAL(theta)[69], 2) + R_pow_di(REAL(theta)[70], 2) + R_pow_di(REAL(theta)[71], 2) + R_pow_di(REAL(theta)[72], 2) + R_pow_di(REAL(theta)[73], 2) + R_pow_di(REAL(theta)[74], 2) + R_pow_di(REAL(theta)[75], 2) + R_pow_di(REAL(theta)[76], 2) + R_pow_di(REAL(theta)[77], 4);
    }
    else if (theta_n == 1){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
      REAL(ret)[1] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[2] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[3] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[4] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[5] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[6] = 2 * REAL(theta)[0] * REAL(theta)[21];
      REAL(ret)[7] = 2 * REAL(theta)[0] * REAL(theta)[28];
      REAL(ret)[8] = 2 * REAL(theta)[0] * REAL(theta)[36];
      REAL(ret)[9] = 2 * REAL(theta)[0] * REAL(theta)[45];
      REAL(ret)[10] = 2 * REAL(theta)[0] * REAL(theta)[55];
      REAL(ret)[11] = 2 * REAL(theta)[0] * REAL(theta)[66];
      REAL(ret)[12] = 2 * REAL(theta)[0] * REAL(theta)[1];
      REAL(ret)[24] = 2 * REAL(theta)[0] * REAL(theta)[3];
      REAL(ret)[36] = 2 * REAL(theta)[0] * REAL(theta)[6];
      REAL(ret)[48] = 2 * REAL(theta)[0] * REAL(theta)[10];
      REAL(ret)[60] = 2 * REAL(theta)[0] * REAL(theta)[15];
      REAL(ret)[72] = 2 * REAL(theta)[0] * REAL(theta)[21];
      REAL(ret)[84] = 2 * REAL(theta)[0] * REAL(theta)[28];
      REAL(ret)[96] = 2 * REAL(theta)[0] * REAL(theta)[36];
      REAL(ret)[108] = 2 * REAL(theta)[0] * REAL(theta)[45];
      REAL(ret)[120] = 2 * REAL(theta)[0] * REAL(theta)[55];
      REAL(ret)[132] = 2 * REAL(theta)[0] * REAL(theta)[66];
    }
    else if (theta_n == 2){
      REAL(ret)[1] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[12] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[13] = 4 * R_pow_di(REAL(theta)[2], 3);
      REAL(ret)[14] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[15] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[16] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[17] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[18] = 2 * REAL(theta)[2] * REAL(theta)[22];
      REAL(ret)[19] = 2 * REAL(theta)[2] * REAL(theta)[29];
      REAL(ret)[20] = 2 * REAL(theta)[2] * REAL(theta)[37];
      REAL(ret)[21] = 2 * REAL(theta)[2] * REAL(theta)[46];
      REAL(ret)[22] = 2 * REAL(theta)[2] * REAL(theta)[56];
      REAL(ret)[23] = 2 * REAL(theta)[2] * REAL(theta)[67];
      REAL(ret)[25] = 2 * REAL(theta)[2] * REAL(theta)[4];
      REAL(ret)[37] = 2 * REAL(theta)[2] * REAL(theta)[7];
      REAL(ret)[49] = 2 * REAL(theta)[11] * REAL(theta)[2];
      REAL(ret)[61] = 2 * REAL(theta)[16] * REAL(theta)[2];
      REAL(ret)[73] = 2 * REAL(theta)[2] * REAL(theta)[22];
      REAL(ret)[85] = 2 * REAL(theta)[2] * REAL(theta)[29];
      REAL(ret)[97] = 2 * REAL(theta)[2] * REAL(theta)[37];
      REAL(ret)[109] = 2 * REAL(theta)[2] * REAL(theta)[46];
      REAL(ret)[121] = 2 * REAL(theta)[2] * REAL(theta)[56];
      REAL(ret)[133] = 2 * REAL(theta)[2] * REAL(theta)[67];
    }
    else if (theta_n == 4){
      REAL(ret)[2] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[14] = REAL(theta)[1];
      REAL(ret)[24] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[14] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[25] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[26] = 4 * R_pow_di(REAL(theta)[5], 3);
      REAL(ret)[27] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[28] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[29] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[30] = 2 * REAL(theta)[23] * REAL(theta)[5];
      REAL(ret)[31] = 2 * REAL(theta)[30] * REAL(theta)[5];
      REAL(ret)[32] = 2 * REAL(theta)[38] * REAL(theta)[5];
      REAL(ret)[33] = 2 * REAL(theta)[47] * REAL(theta)[5];
      REAL(ret)[34] = 2 * REAL(theta)[5] * REAL(theta)[57];
      REAL(ret)[35] = 2 * REAL(theta)[5] * REAL(theta)[68];
      REAL(ret)[38] = 2 * REAL(theta)[5] * REAL(theta)[8];
      REAL(ret)[50] = 2 * REAL(theta)[12] * REAL(theta)[5];
      REAL(ret)[62] = 2 * REAL(theta)[17] * REAL(theta)[5];
      REAL(ret)[74] = 2 * REAL(theta)[23] * REAL(theta)[5];
      REAL(ret)[86] = 2 * REAL(theta)[30] * REAL(theta)[5];
      REAL(ret)[98] = 2 * REAL(theta)[38] * REAL(theta)[5];
      REAL(ret)[110] = 2 * REAL(theta)[47] * REAL(theta)[5];
      REAL(ret)[122] = 2 * REAL(theta)[5] * REAL(theta)[57];
      REAL(ret)[134] = 2 * REAL(theta)[5] * REAL(theta)[68];
    }
    else if (theta_n == 7){
      REAL(ret)[3] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[15] = REAL(theta)[1];
      REAL(ret)[27] = REAL(theta)[3];
      REAL(ret)[36] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[15] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[27] = REAL(theta)[4];
      REAL(ret)[37] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[27] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[38] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[39] = 4 * R_pow_di(REAL(theta)[9], 3);
      REAL(ret)[40] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[41] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[42] = 2 * REAL(theta)[24] * REAL(theta)[9];
      REAL(ret)[43] = 2 * REAL(theta)[31] * REAL(theta)[9];
      REAL(ret)[44] = 2 * REAL(theta)[39] * REAL(theta)[9];
      REAL(ret)[45] = 2 * REAL(theta)[48] * REAL(theta)[9];
      REAL(ret)[46] = 2 * REAL(theta)[58] * REAL(theta)[9];
      REAL(ret)[47] = 2 * REAL(theta)[69] * REAL(theta)[9];
      REAL(ret)[51] = 2 * REAL(theta)[13] * REAL(theta)[9];
      REAL(ret)[63] = 2 * REAL(theta)[18] * REAL(theta)[9];
      REAL(ret)[75] = 2 * REAL(theta)[24] * REAL(theta)[9];
      REAL(ret)[87] = 2 * REAL(theta)[31] * REAL(theta)[9];
      REAL(ret)[99] = 2 * REAL(theta)[39] * REAL(theta)[9];
      REAL(ret)[111] = 2 * REAL(theta)[48] * REAL(theta)[9];
      REAL(ret)[123] = 2 * REAL(theta)[58] * REAL(theta)[9];
      REAL(ret)[135] = 2 * REAL(theta)[69] * REAL(theta)[9];
    }
    else if (theta_n == 11){
      REAL(ret)[4] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[16] = REAL(theta)[1];
      REAL(ret)[28] = REAL(theta)[3];
      REAL(ret)[40] = REAL(theta)[6];
      REAL(ret)[48] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[16] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[28] = REAL(theta)[4];
      REAL(ret)[40] = REAL(theta)[7];
      REAL(ret)[49] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[28] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[40] = REAL(theta)[8];
      REAL(ret)[50] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[40] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[51] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[52] = 4 * R_pow_di(REAL(theta)[14], 3);
      REAL(ret)[53] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[54] = 2 * REAL(theta)[14] * REAL(theta)[25];
      REAL(ret)[55] = 2 * REAL(theta)[14] * REAL(theta)[32];
      REAL(ret)[56] = 2 * REAL(theta)[14] * REAL(theta)[40];
      REAL(ret)[57] = 2 * REAL(theta)[14] * REAL(theta)[49];
      REAL(ret)[58] = 2 * REAL(theta)[14] * REAL(theta)[59];
      REAL(ret)[59] = 2 * REAL(theta)[14] * REAL(theta)[70];
      REAL(ret)[64] = 2 * REAL(theta)[14] * REAL(theta)[19];
      REAL(ret)[76] = 2 * REAL(theta)[14] * REAL(theta)[25];
      REAL(ret)[88] = 2 * REAL(theta)[14] * REAL(theta)[32];
      REAL(ret)[100] = 2 * REAL(theta)[14] * REAL(theta)[40];
      REAL(ret)[112] = 2 * REAL(theta)[14] * REAL(theta)[49];
      REAL(ret)[124] = 2 * REAL(theta)[14] * REAL(theta)[59];
      REAL(ret)[136] = 2 * REAL(theta)[14] * REAL(theta)[70];
    }
    else if (theta_n == 16){
      REAL(ret)[5] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[17] = REAL(theta)[1];
      REAL(ret)[29] = REAL(theta)[3];
      REAL(ret)[41] = REAL(theta)[6];
      REAL(ret)[53] = REAL(theta)[10];
      REAL(ret)[60] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[17] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[29] = REAL(theta)[4];
      REAL(ret)[41] = REAL(theta)[7];
      REAL(ret)[53] = REAL(theta)[11];
      REAL(ret)[61] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[29] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[41] = REAL(theta)[8];
      REAL(ret)[53] = REAL(theta)[12];
      REAL(ret)[62] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[41] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[53] = REAL(theta)[13];
      REAL(ret)[63] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[53] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[64] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[65] = 4 * R_pow_di(REAL(theta)[20], 3);
      REAL(ret)[66] = 2 * REAL(theta)[20] * REAL(theta)[26];
      REAL(ret)[67] = 2 * REAL(theta)[20] * REAL(theta)[33];
      REAL(ret)[68] = 2 * REAL(theta)[20] * REAL(theta)[41];
      REAL(ret)[69] = 2 * REAL(theta)[20] * REAL(theta)[50];
      REAL(ret)[70] = 2 * REAL(theta)[20] * REAL(theta)[60];
      REAL(ret)[71] = 2 * REAL(theta)[20] * REAL(theta)[71];
      REAL(ret)[77] = 2 * REAL(theta)[20] * REAL(theta)[26];
      REAL(ret)[89] = 2 * REAL(theta)[20] * REAL(theta)[33];
      REAL(ret)[101] = 2 * REAL(theta)[20] * REAL(theta)[41];
      REAL(ret)[113] = 2 * REAL(theta)[20] * REAL(theta)[50];
      REAL(ret)[125] = 2 * REAL(theta)[20] * REAL(theta)[60];
      REAL(ret)[137] = 2 * REAL(theta)[20] * REAL(theta)[71];
    }
    else if (theta_n == 22){
      REAL(ret)[6] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[18] = REAL(theta)[1];
      REAL(ret)[30] = REAL(theta)[3];
      REAL(ret)[42] = REAL(theta)[6];
      REAL(ret)[54] = REAL(theta)[10];
      REAL(ret)[66] = REAL(theta)[15];
      REAL(ret)[72] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[18] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[30] = REAL(theta)[4];
      REAL(ret)[42] = REAL(theta)[7];
      REAL(ret)[54] = REAL(theta)[11];
      REAL(ret)[66] = REAL(theta)[16];
      REAL(ret)[73] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[30] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[42] = REAL(theta)[8];
      REAL(ret)[54] = REAL(theta)[12];
      REAL(ret)[66] = REAL(theta)[17];
      REAL(ret)[74] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[42] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[54] = REAL(theta)[13];
      REAL(ret)[66] = REAL(theta)[18];
      REAL(ret)[75] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[54] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[66] = REAL(theta)[19];
      REAL(ret)[76] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[66] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[77] = R_pow_di(REAL(theta)[20], 2);
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
      REAL(ret)[78] = 4 * R_pow_di(REAL(theta)[27], 3);
      REAL(ret)[79] = 2 * REAL(theta)[27] * REAL(theta)[34];
      REAL(ret)[80] = 2 * REAL(theta)[27] * REAL(theta)[42];
      REAL(ret)[81] = 2 * REAL(theta)[27] * REAL(theta)[51];
      REAL(ret)[82] = 2 * REAL(theta)[27] * REAL(theta)[61];
      REAL(ret)[83] = 2 * REAL(theta)[27] * REAL(theta)[72];
      REAL(ret)[90] = 2 * REAL(theta)[27] * REAL(theta)[34];
      REAL(ret)[102] = 2 * REAL(theta)[27] * REAL(theta)[42];
      REAL(ret)[114] = 2 * REAL(theta)[27] * REAL(theta)[51];
      REAL(ret)[126] = 2 * REAL(theta)[27] * REAL(theta)[61];
      REAL(ret)[138] = 2 * REAL(theta)[27] * REAL(theta)[72];
    }
    else if (theta_n == 29){
      REAL(ret)[7] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[19] = REAL(theta)[1];
      REAL(ret)[31] = REAL(theta)[3];
      REAL(ret)[43] = REAL(theta)[6];
      REAL(ret)[55] = REAL(theta)[10];
      REAL(ret)[67] = REAL(theta)[15];
      REAL(ret)[79] = REAL(theta)[21];
      REAL(ret)[84] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[19] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[31] = REAL(theta)[4];
      REAL(ret)[43] = REAL(theta)[7];
      REAL(ret)[55] = REAL(theta)[11];
      REAL(ret)[67] = REAL(theta)[16];
      REAL(ret)[79] = REAL(theta)[22];
      REAL(ret)[85] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[31] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[43] = REAL(theta)[8];
      REAL(ret)[55] = REAL(theta)[12];
      REAL(ret)[67] = REAL(theta)[17];
      REAL(ret)[79] = REAL(theta)[23];
      REAL(ret)[86] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[43] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[55] = REAL(theta)[13];
      REAL(ret)[67] = REAL(theta)[18];
      REAL(ret)[79] = REAL(theta)[24];
      REAL(ret)[87] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[55] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[67] = REAL(theta)[19];
      REAL(ret)[79] = REAL(theta)[25];
      REAL(ret)[88] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[67] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[79] = REAL(theta)[26];
      REAL(ret)[89] = R_pow_di(REAL(theta)[20], 2);
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
      REAL(ret)[79] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[90] = R_pow_di(REAL(theta)[27], 2);
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
      REAL(ret)[91] = 4 * R_pow_di(REAL(theta)[35], 3);
      REAL(ret)[92] = 2 * REAL(theta)[35] * REAL(theta)[43];
      REAL(ret)[93] = 2 * REAL(theta)[35] * REAL(theta)[52];
      REAL(ret)[94] = 2 * REAL(theta)[35] * REAL(theta)[62];
      REAL(ret)[95] = 2 * REAL(theta)[35] * REAL(theta)[73];
      REAL(ret)[103] = 2 * REAL(theta)[35] * REAL(theta)[43];
      REAL(ret)[115] = 2 * REAL(theta)[35] * REAL(theta)[52];
      REAL(ret)[127] = 2 * REAL(theta)[35] * REAL(theta)[62];
      REAL(ret)[139] = 2 * REAL(theta)[35] * REAL(theta)[73];
    }
    else if (theta_n == 37){
      REAL(ret)[8] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[20] = REAL(theta)[1];
      REAL(ret)[32] = REAL(theta)[3];
      REAL(ret)[44] = REAL(theta)[6];
      REAL(ret)[56] = REAL(theta)[10];
      REAL(ret)[68] = REAL(theta)[15];
      REAL(ret)[80] = REAL(theta)[21];
      REAL(ret)[92] = REAL(theta)[28];
      REAL(ret)[96] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[20] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[32] = REAL(theta)[4];
      REAL(ret)[44] = REAL(theta)[7];
      REAL(ret)[56] = REAL(theta)[11];
      REAL(ret)[68] = REAL(theta)[16];
      REAL(ret)[80] = REAL(theta)[22];
      REAL(ret)[92] = REAL(theta)[29];
      REAL(ret)[97] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[32] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[44] = REAL(theta)[8];
      REAL(ret)[56] = REAL(theta)[12];
      REAL(ret)[68] = REAL(theta)[17];
      REAL(ret)[80] = REAL(theta)[23];
      REAL(ret)[92] = REAL(theta)[30];
      REAL(ret)[98] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[44] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[56] = REAL(theta)[13];
      REAL(ret)[68] = REAL(theta)[18];
      REAL(ret)[80] = REAL(theta)[24];
      REAL(ret)[92] = REAL(theta)[31];
      REAL(ret)[99] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[56] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[68] = REAL(theta)[19];
      REAL(ret)[80] = REAL(theta)[25];
      REAL(ret)[92] = REAL(theta)[32];
      REAL(ret)[100] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[68] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[80] = REAL(theta)[26];
      REAL(ret)[92] = REAL(theta)[33];
      REAL(ret)[101] = R_pow_di(REAL(theta)[20], 2);
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
      REAL(ret)[80] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[92] = REAL(theta)[34];
      REAL(ret)[102] = R_pow_di(REAL(theta)[27], 2);
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
      REAL(ret)[92] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[103] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[104] = 2 * REAL(theta)[43];
      REAL(ret)[105] = REAL(theta)[52];
      REAL(ret)[106] = REAL(theta)[62];
      REAL(ret)[107] = REAL(theta)[73];
      REAL(ret)[116] = REAL(theta)[52];
      REAL(ret)[128] = REAL(theta)[62];
      REAL(ret)[140] = REAL(theta)[73];
    }
    else if (theta_n == 45){
      REAL(ret)[104] = 4 * R_pow_di(REAL(theta)[44], 3);
      REAL(ret)[105] = 2 * REAL(theta)[44] * REAL(theta)[53];
      REAL(ret)[106] = 2 * REAL(theta)[44] * REAL(theta)[63];
      REAL(ret)[107] = 2 * REAL(theta)[44] * REAL(theta)[74];
      REAL(ret)[116] = 2 * REAL(theta)[44] * REAL(theta)[53];
      REAL(ret)[128] = 2 * REAL(theta)[44] * REAL(theta)[63];
      REAL(ret)[140] = 2 * REAL(theta)[44] * REAL(theta)[74];
    }
    else if (theta_n == 46){
      REAL(ret)[9] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[21] = REAL(theta)[1];
      REAL(ret)[33] = REAL(theta)[3];
      REAL(ret)[45] = REAL(theta)[6];
      REAL(ret)[57] = REAL(theta)[10];
      REAL(ret)[69] = REAL(theta)[15];
      REAL(ret)[81] = REAL(theta)[21];
      REAL(ret)[93] = REAL(theta)[28];
      REAL(ret)[105] = REAL(theta)[36];
      REAL(ret)[108] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[21] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[33] = REAL(theta)[4];
      REAL(ret)[45] = REAL(theta)[7];
      REAL(ret)[57] = REAL(theta)[11];
      REAL(ret)[69] = REAL(theta)[16];
      REAL(ret)[81] = REAL(theta)[22];
      REAL(ret)[93] = REAL(theta)[29];
      REAL(ret)[105] = REAL(theta)[37];
      REAL(ret)[109] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[33] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[45] = REAL(theta)[8];
      REAL(ret)[57] = REAL(theta)[12];
      REAL(ret)[69] = REAL(theta)[17];
      REAL(ret)[81] = REAL(theta)[23];
      REAL(ret)[93] = REAL(theta)[30];
      REAL(ret)[105] = REAL(theta)[38];
      REAL(ret)[110] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[45] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[57] = REAL(theta)[13];
      REAL(ret)[69] = REAL(theta)[18];
      REAL(ret)[81] = REAL(theta)[24];
      REAL(ret)[93] = REAL(theta)[31];
      REAL(ret)[105] = REAL(theta)[39];
      REAL(ret)[111] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[57] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[69] = REAL(theta)[19];
      REAL(ret)[81] = REAL(theta)[25];
      REAL(ret)[93] = REAL(theta)[32];
      REAL(ret)[105] = REAL(theta)[40];
      REAL(ret)[112] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[69] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[81] = REAL(theta)[26];
      REAL(ret)[93] = REAL(theta)[33];
      REAL(ret)[105] = REAL(theta)[41];
      REAL(ret)[113] = R_pow_di(REAL(theta)[20], 2);
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
      REAL(ret)[81] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[93] = REAL(theta)[34];
      REAL(ret)[105] = REAL(theta)[42];
      REAL(ret)[114] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[115] = REAL(theta)[34];
      REAL(ret)[116] = REAL(theta)[42];
      REAL(ret)[117] = 2 * REAL(theta)[51];
      REAL(ret)[118] = REAL(theta)[61];
      REAL(ret)[119] = REAL(theta)[72];
      REAL(ret)[129] = REAL(theta)[61];
      REAL(ret)[141] = REAL(theta)[72];
    }
    else if (theta_n == 53){
      REAL(ret)[93] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[105] = REAL(theta)[43];
      REAL(ret)[115] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[116] = REAL(theta)[43];
      REAL(ret)[117] = 2 * REAL(theta)[52];
      REAL(ret)[118] = REAL(theta)[62];
      REAL(ret)[119] = REAL(theta)[73];
      REAL(ret)[129] = REAL(theta)[62];
      REAL(ret)[141] = REAL(theta)[73];
    }
    else if (theta_n == 54){
      REAL(ret)[105] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[116] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[117] = 2 * REAL(theta)[53];
      REAL(ret)[118] = REAL(theta)[63];
      REAL(ret)[119] = REAL(theta)[74];
      REAL(ret)[129] = REAL(theta)[63];
      REAL(ret)[141] = REAL(theta)[74];
    }
    else if (theta_n == 55){
      REAL(ret)[117] = 4 * R_pow_di(REAL(theta)[54], 3);
      REAL(ret)[118] = 2 * REAL(theta)[54] * REAL(theta)[64];
      REAL(ret)[119] = 2 * REAL(theta)[54] * REAL(theta)[75];
      REAL(ret)[129] = 2 * REAL(theta)[54] * REAL(theta)[64];
      REAL(ret)[141] = 2 * REAL(theta)[54] * REAL(theta)[75];
    }
    else if (theta_n == 56){
      REAL(ret)[10] = R_pow_di(REAL(theta)[0], 2);
      REAL(ret)[22] = REAL(theta)[1];
      REAL(ret)[34] = REAL(theta)[3];
      REAL(ret)[46] = REAL(theta)[6];
      REAL(ret)[58] = REAL(theta)[10];
      REAL(ret)[70] = REAL(theta)[15];
      REAL(ret)[82] = REAL(theta)[21];
      REAL(ret)[94] = REAL(theta)[28];
      REAL(ret)[106] = REAL(theta)[36];
      REAL(ret)[118] = REAL(theta)[45];
      REAL(ret)[120] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[22] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[34] = REAL(theta)[4];
      REAL(ret)[46] = REAL(theta)[7];
      REAL(ret)[58] = REAL(theta)[11];
      REAL(ret)[70] = REAL(theta)[16];
      REAL(ret)[82] = REAL(theta)[22];
      REAL(ret)[94] = REAL(theta)[29];
      REAL(ret)[106] = REAL(theta)[37];
      REAL(ret)[118] = REAL(theta)[46];
      REAL(ret)[121] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[34] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[46] = REAL(theta)[8];
      REAL(ret)[58] = REAL(theta)[12];
      REAL(ret)[70] = REAL(theta)[17];
      REAL(ret)[82] = REAL(theta)[23];
      REAL(ret)[94] = REAL(theta)[30];
      REAL(ret)[106] = REAL(theta)[38];
      REAL(ret)[118] = REAL(theta)[47];
      REAL(ret)[122] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[46] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[58] = REAL(theta)[13];
      REAL(ret)[70] = REAL(theta)[18];
      REAL(ret)[82] = REAL(theta)[24];
      REAL(ret)[94] = REAL(theta)[31];
      REAL(ret)[106] = REAL(theta)[39];
      REAL(ret)[118] = REAL(theta)[48];
      REAL(ret)[123] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[58] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[70] = REAL(theta)[19];
      REAL(ret)[82] = REAL(theta)[25];
      REAL(ret)[94] = REAL(theta)[32];
      REAL(ret)[106] = REAL(theta)[40];
      REAL(ret)[118] = REAL(theta)[49];
      REAL(ret)[124] = R_pow_di(REAL(theta)[14], 2);
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
      REAL(ret)[70] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[82] = REAL(theta)[26];
      REAL(ret)[94] = REAL(theta)[33];
      REAL(ret)[106] = REAL(theta)[41];
      REAL(ret)[118] = REAL(theta)[50];
      REAL(ret)[125] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[126] = REAL(theta)[26];
      REAL(ret)[127] = REAL(theta)[33];
      REAL(ret)[128] = REAL(theta)[41];
      REAL(ret)[129] = REAL(theta)[50];
      REAL(ret)[130] = 2 * REAL(theta)[60];
      REAL(ret)[131] = REAL(theta)[71];
      REAL(ret)[142] = REAL(theta)[71];
    }
    else if (theta_n == 62){
      REAL(ret)[82] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[94] = REAL(theta)[34];
      REAL(ret)[106] = REAL(theta)[42];
      REAL(ret)[118] = REAL(theta)[51];
      REAL(ret)[126] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[127] = REAL(theta)[34];
      REAL(ret)[128] = REAL(theta)[42];
      REAL(ret)[129] = REAL(theta)[51];
      REAL(ret)[130] = 2 * REAL(theta)[61];
      REAL(ret)[131] = REAL(theta)[72];
      REAL(ret)[142] = REAL(theta)[72];
    }
    else if (theta_n == 63){
      REAL(ret)[94] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[106] = REAL(theta)[43];
      REAL(ret)[118] = REAL(theta)[52];
      REAL(ret)[127] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[128] = REAL(theta)[43];
      REAL(ret)[129] = REAL(theta)[52];
      REAL(ret)[130] = 2 * REAL(theta)[62];
      REAL(ret)[131] = REAL(theta)[73];
      REAL(ret)[142] = REAL(theta)[73];
    }
    else if (theta_n == 64){
      REAL(ret)[106] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[118] = REAL(theta)[53];
      REAL(ret)[128] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[129] = REAL(theta)[53];
      REAL(ret)[130] = 2 * REAL(theta)[63];
      REAL(ret)[131] = REAL(theta)[74];
      REAL(ret)[142] = REAL(theta)[74];
    }
    else if (theta_n == 65){
      REAL(ret)[118] = R_pow_di(REAL(theta)[54], 2);
      REAL(ret)[129] = R_pow_di(REAL(theta)[54], 2);
      REAL(ret)[130] = 2 * REAL(theta)[64];
      REAL(ret)[131] = REAL(theta)[75];
      REAL(ret)[142] = REAL(theta)[75];
    }
    else if (theta_n == 66){
      REAL(ret)[130] = 4 * R_pow_di(REAL(theta)[65], 3);
      REAL(ret)[131] = 2 * REAL(theta)[65] * REAL(theta)[76];
      REAL(ret)[142] = 2 * REAL(theta)[65] * REAL(theta)[76];
    }
    else if (theta_n == 67){
      REAL(ret)[11] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[132] = R_pow_di(REAL(theta)[0], 2);
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
      REAL(ret)[23] = R_pow_di(REAL(theta)[2], 2);
      REAL(ret)[35] = REAL(theta)[4];
      REAL(ret)[47] = REAL(theta)[7];
      REAL(ret)[59] = REAL(theta)[11];
      REAL(ret)[71] = REAL(theta)[16];
      REAL(ret)[83] = REAL(theta)[22];
      REAL(ret)[95] = REAL(theta)[29];
      REAL(ret)[107] = REAL(theta)[37];
      REAL(ret)[119] = REAL(theta)[46];
      REAL(ret)[131] = REAL(theta)[56];
      REAL(ret)[133] = R_pow_di(REAL(theta)[2], 2);
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
      REAL(ret)[35] = R_pow_di(REAL(theta)[5], 2);
      REAL(ret)[47] = REAL(theta)[8];
      REAL(ret)[59] = REAL(theta)[12];
      REAL(ret)[71] = REAL(theta)[17];
      REAL(ret)[83] = REAL(theta)[23];
      REAL(ret)[95] = REAL(theta)[30];
      REAL(ret)[107] = REAL(theta)[38];
      REAL(ret)[119] = REAL(theta)[47];
      REAL(ret)[131] = REAL(theta)[57];
      REAL(ret)[134] = R_pow_di(REAL(theta)[5], 2);
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
      REAL(ret)[47] = R_pow_di(REAL(theta)[9], 2);
      REAL(ret)[59] = REAL(theta)[13];
      REAL(ret)[71] = REAL(theta)[18];
      REAL(ret)[83] = REAL(theta)[24];
      REAL(ret)[95] = REAL(theta)[31];
      REAL(ret)[107] = REAL(theta)[39];
      REAL(ret)[119] = REAL(theta)[48];
      REAL(ret)[131] = REAL(theta)[58];
      REAL(ret)[135] = R_pow_di(REAL(theta)[9], 2);
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
      REAL(ret)[59] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[71] = REAL(theta)[19];
      REAL(ret)[83] = REAL(theta)[25];
      REAL(ret)[95] = REAL(theta)[32];
      REAL(ret)[107] = REAL(theta)[40];
      REAL(ret)[119] = REAL(theta)[49];
      REAL(ret)[131] = REAL(theta)[59];
      REAL(ret)[136] = R_pow_di(REAL(theta)[14], 2);
      REAL(ret)[137] = REAL(theta)[19];
      REAL(ret)[138] = REAL(theta)[25];
      REAL(ret)[139] = REAL(theta)[32];
      REAL(ret)[140] = REAL(theta)[40];
      REAL(ret)[141] = REAL(theta)[49];
      REAL(ret)[142] = REAL(theta)[59];
      REAL(ret)[143] = 2 * REAL(theta)[70];
    }
    else if (theta_n == 72){
      REAL(ret)[71] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[83] = REAL(theta)[26];
      REAL(ret)[95] = REAL(theta)[33];
      REAL(ret)[107] = REAL(theta)[41];
      REAL(ret)[119] = REAL(theta)[50];
      REAL(ret)[131] = REAL(theta)[60];
      REAL(ret)[137] = R_pow_di(REAL(theta)[20], 2);
      REAL(ret)[138] = REAL(theta)[26];
      REAL(ret)[139] = REAL(theta)[33];
      REAL(ret)[140] = REAL(theta)[41];
      REAL(ret)[141] = REAL(theta)[50];
      REAL(ret)[142] = REAL(theta)[60];
      REAL(ret)[143] = 2 * REAL(theta)[71];
    }
    else if (theta_n == 73){
      REAL(ret)[83] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[95] = REAL(theta)[34];
      REAL(ret)[107] = REAL(theta)[42];
      REAL(ret)[119] = REAL(theta)[51];
      REAL(ret)[131] = REAL(theta)[61];
      REAL(ret)[138] = R_pow_di(REAL(theta)[27], 2);
      REAL(ret)[139] = REAL(theta)[34];
      REAL(ret)[140] = REAL(theta)[42];
      REAL(ret)[141] = REAL(theta)[51];
      REAL(ret)[142] = REAL(theta)[61];
      REAL(ret)[143] = 2 * REAL(theta)[72];
    }
    else if (theta_n == 74){
      REAL(ret)[95] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[107] = REAL(theta)[43];
      REAL(ret)[119] = REAL(theta)[52];
      REAL(ret)[131] = REAL(theta)[62];
      REAL(ret)[139] = R_pow_di(REAL(theta)[35], 2);
      REAL(ret)[140] = REAL(theta)[43];
      REAL(ret)[141] = REAL(theta)[52];
      REAL(ret)[142] = REAL(theta)[62];
      REAL(ret)[143] = 2 * REAL(theta)[73];
    }
    else if (theta_n == 75){
      REAL(ret)[107] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[119] = REAL(theta)[53];
      REAL(ret)[131] = REAL(theta)[63];
      REAL(ret)[140] = R_pow_di(REAL(theta)[44], 2);
      REAL(ret)[141] = REAL(theta)[53];
      REAL(ret)[142] = REAL(theta)[63];
      REAL(ret)[143] = 2 * REAL(theta)[74];
    }
    else if (theta_n == 76){
      REAL(ret)[119] = R_pow_di(REAL(theta)[54], 2);
      REAL(ret)[131] = REAL(theta)[64];
      REAL(ret)[141] = R_pow_di(REAL(theta)[54], 2);
      REAL(ret)[142] = REAL(theta)[64];
      REAL(ret)[143] = 2 * REAL(theta)[75];
    }
    else if (theta_n == 77){
      REAL(ret)[131] = R_pow_di(REAL(theta)[65], 2);
      REAL(ret)[142] = R_pow_di(REAL(theta)[65], 2);
      REAL(ret)[143] = 2 * REAL(theta)[76];
    }
    else if (theta_n == 78){
      REAL(ret)[143] = 4 * R_pow_di(REAL(theta)[77], 3);
    }
    UNPROTECT(1);
    return(ret);
  } else {
    SEXP ret = PROTECT(allocVector(REALSXP, 12));for(int i = 0; i < 12; i++){REAL(ret)[i]=0;}
    if (theta_n == -3){
      REAL(ret)[0] = 4 * R_pow_di(REAL(theta)[0], 3);
    }
    else if (theta_n == -5){
      REAL(ret)[1] = 4 * R_pow_di(REAL(theta)[2], 3);
    }
    else if (theta_n == -8){
      REAL(ret)[2] = 4 * R_pow_di(REAL(theta)[5], 3);
    }
    else if (theta_n == -12){
      REAL(ret)[3] = 4 * R_pow_di(REAL(theta)[9], 3);
    }
    else if (theta_n == -17){
      REAL(ret)[4] = 4 * R_pow_di(REAL(theta)[14], 3);
    }
    else if (theta_n == -23){
      REAL(ret)[5] = 4 * R_pow_di(REAL(theta)[20], 3);
    }
    else if (theta_n == -30){
      REAL(ret)[6] = 4 * R_pow_di(REAL(theta)[27], 3);
    }
    else if (theta_n == -38){
      REAL(ret)[7] = 4 * R_pow_di(REAL(theta)[35], 3);
    }
    else if (theta_n == -47){
      REAL(ret)[8] = 4 * R_pow_di(REAL(theta)[44], 3);
    }
    else if (theta_n == -57){
      REAL(ret)[9] = 4 * R_pow_di(REAL(theta)[54], 3);
    }
    else if (theta_n == -68){
      REAL(ret)[10] = 4 * R_pow_di(REAL(theta)[65], 3);
    }
    else if (theta_n == -80){
      REAL(ret)[11] = 4 * R_pow_di(REAL(theta)[77], 3);
    }
    UNPROTECT(1);
    return(ret);
  }
}

  return R_NilValue;
}
