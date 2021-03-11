#ifndef __checkmate_H__
#define __checkmate_H__

#if defined(__cplusplus)
SEXP qassertS(SEXP in, const char *test, const char *what);
bool qtest(SEXP in, const char *test);
SEXP qstrictS(SEXP nn, const char *what);
SEXP qstrictSn(SEXP x_, const char *what);
SEXP qstrictSdn(SEXP x_, const char *what);
extern "C" {
#endif
  
  R_xlen_t check_strict_names(SEXP x);
  R_xlen_t find_missing_string(SEXP x);
  void qstrict(SEXP x, const char *what);

#if defined(__cplusplus)
}
#endif

#endif
