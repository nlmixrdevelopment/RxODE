#include <sys/stat.h> 
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>   /* dj: import intptr_t */
#include <errno.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif


SEXP _vecDF(SEXP cv, SEXP n_) {
  int n=0;
  int typ = TYPEOF(n_);
  if (typ == REALSXP) {
    n = (int)(REAL(n_)[0]);
  } else if (typ == INTSXP) {
    n = INTEGER(n_)[0];
  }
  if (n <= 0) error(_("'n' must be greater than 0"));
  int pro = 0;
  int len = length(cv);
  SEXP ret = PROTECT(allocVector(VECSXP, len)); pro++;
  SEXP retN = PROTECT(allocVector(STRSXP, len)); pro++;
  SEXP cvN = getAttrib(cv, R_NamesSymbol);
  for (int i = len; i--;) {
    SEXP tmp = PROTECT(allocVector(REALSXP, n)); pro++;
    for (int j = n; j--;) {
      REAL(tmp)[j] = REAL(cv)[i];
    }
    SET_VECTOR_ELT(ret, i, tmp);
    SET_STRING_ELT(retN, i, STRING_ELT(cvN, i));
  }
  SEXP sexp_rownames = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(sexp_rownames)[0] = NA_INTEGER;
  INTEGER(sexp_rownames)[1] = -n;
  setAttrib(ret, R_RowNamesSymbol, sexp_rownames);
  SEXP sexp_class = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(sexp_class,0,mkChar("data.frame"));
  setAttrib(ret, R_ClassSymbol, sexp_class);
  setAttrib(ret, R_NamesSymbol, retN);
  UNPROTECT(pro);
  return ret;
}

SEXP _cbindOme(SEXP et_, SEXP mat_, SEXP n_) {
  int n = INTEGER(n_)[0];
  if (n <= 0) error(_("'n' must be greater than 0"));
  
  int len1 = length(et_);
  int len1a = 0;
  if (len1 > 0) {
    len1a = length(VECTOR_ELT(et_,0));
  }
  SEXP etN = getAttrib(et_, R_NamesSymbol);

  SEXP matD;
  SEXP matDN;
  int len2;
  int lenOut;
  int lenItem;
  if (!Rf_isNull(mat_) && !Rf_isNull(et_)) {
    matD = getAttrib(mat_, install("dim"));
    matDN = VECTOR_ELT(getAttrib(mat_, R_DimNamesSymbol), 1);
    len2 = INTEGER(matD)[1];
    lenOut = INTEGER(matD)[0];
    lenItem = lenOut/len1a;
  } else if (!Rf_isNull(et_)) {
    len2 = 0;
    lenOut = n*len1a;
    lenItem = n;
  } else {
    matD = getAttrib(mat_, install("dim"));;
    matDN = VECTOR_ELT(getAttrib(mat_, R_DimNamesSymbol), 1);
    len2 = INTEGER(matD)[1];
    lenOut = INTEGER(matD)[0];
    lenItem = n;
  }
  int pro = 0;
  SEXP ret = PROTECT(allocVector(VECSXP, len1+len2)); pro++;
  SEXP retN = PROTECT(allocVector(STRSXP, len1+len2)); pro++;
  for (int i = len1; i--; ) {
    SEXP tmp = PROTECT(allocVector(REALSXP, lenOut)); pro++;
    SEXP in = VECTOR_ELT(et_, i);
    int l = lenOut;
    for (int j = len1a; j--;) {
      for (int k = lenItem; k--; ) {
	REAL(tmp)[--l] = REAL(in)[j];
      }
    }
    SET_VECTOR_ELT(ret, i, tmp);
    SET_STRING_ELT(retN, i, STRING_ELT(etN, i));
  }
  for (int i = len2; i--; ) {
    SEXP tmp = PROTECT(allocVector(REALSXP, lenOut)); pro++;
    memcpy(&(REAL(tmp)[0]), &(REAL(mat_)[lenOut*i]), lenOut*sizeof(double));
    SET_VECTOR_ELT(ret, i+len1, tmp);
    SET_STRING_ELT(retN, i+len1, STRING_ELT(matDN, i));
  }
  SEXP sexp_rownames = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(sexp_rownames)[0] = NA_INTEGER;
  INTEGER(sexp_rownames)[1] = -lenOut;
  setAttrib(ret, R_RowNamesSymbol, sexp_rownames);
  SEXP sexp_class = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(sexp_class,0,mkChar("data.frame"));
  setAttrib(ret, R_ClassSymbol, sexp_class);
  setAttrib(ret, R_NamesSymbol, retN);
  UNPROTECT(pro);
  return ret;
}
