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
  int n = INTEGER(n_)[0];
  if (n <= 0) error(_("'n' must be greater than 0"));
  int pro = 0;
  int len = length(cv);
  SEXP ret = PROTECT(allocVector(VECSXP, len)); pro++;
  SEXP retN = PROTECT(allocVector(STRSXP, len)); pro++;
  SEXP cvN = getAttrib(cv, install("names"));
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
  setAttrib(ret, install("class"), sexp_class);
  setAttrib(ret, install("names"), retN);
  UNPROTECT(pro);
  return ret;
}
