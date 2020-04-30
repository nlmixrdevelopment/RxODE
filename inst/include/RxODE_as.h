#ifndef RxODE_as
#define RxODE_as
#include <Rcpp.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

static inline bool rxIsNull(RObject obj) {
  return obj.sexp_type() == 0;
}
static inline bool rxIsNum(RObject obj) {
  if (obj.sexp_type() == REALSXP) {
    return (!obj.hasAttribute("dim"));
  }
  return false;
}
static inline bool rxIsNum1(RObject obj) {
  if (obj.sexp_type() == REALSXP) {
    if (!obj.hasAttribute("dim")){
      return (Rf_length(obj) == 1);
    }
  }
  return false;
}
static inline bool rxIsInt(RObject obj) {
  if (obj.sexp_type() == 13) {
    return (!obj.hasAttribute("dim"));
  }
  return false;
}

static inline bool rxIsNumInt(RObject obj) {
  int type = obj.sexp_type();
  if (type == REALSXP || type == 13) {
    return (!obj.hasAttribute("dim"));
  }
  return false;
}
static inline bool rxIsChar(RObject obj) {
  int type = obj.sexp_type();
  if (type == STRSXP) {
    return (!obj.hasAttribute("dim"));
  }
  return false;
}

static inline int asInt(SEXP in, const char* what) {
  if (Rf_length(in) != 1 || !qtest(in,"x")) {
    Rcpp::stop(_("'%s' needs to be an integer"), what);
  }
  return as<int>(in);
}

static inline unsigned int asUnsignedInt(SEXP in, const char* what) {
  if (Rf_length(in) != 1 || !qtest(in, "x[0,)")) {
    Rcpp::stop(_("'%s' needs to be an integer greater than 0"), what);
  }
  return as<unsigned int>(in);
}

static inline double asDouble(SEXP in, const char* what) {
  int type = TYPEOF(in);
  if (Rf_length(in) != 1 || type != REALSXP) {
    Rcpp::stop(_("'%s' needs to be an double"), what);
  } 
  return as<double>(in);
}

static inline bool asBool(SEXP in, const char *what) {
  int type = TYPEOF(in);
  if (Rf_length(in) != 1 || type != LGLSXP){
    Rcpp::stop(_("'%s' needs to be a boolean"), what);
  }
  return as<bool>(in);
}

static inline std::string asStr(SEXP in, const char *what) {
  int type = TYPEOF(in);
  if (Rf_length(in) != 1 || type != CHARSXP){
    Rcpp::stop(_("'%s' needs to be a string"), what);
  }
  return as<std::string>(in);
}

static inline CharacterVector asCv(SEXP in, const char *what) {
  int type = TYPEOF(in);
  if (type != STRSXP){
    Rcpp::stop(_("'%s' needs to be a vector of strings"), what);
  } 
  return as<CharacterVector>(in);
}

static inline Nullable<LogicalVector> asNLv(SEXP in, const char* what) {
  int type = TYPEOF(in);
  if (!(type == 0 || type == LGLSXP)){
    Rcpp::stop(_("'%s' needs to be an logical vector or NULL"), what);
  }
  return as<Nullable<LogicalVector>>(in);
}

static inline NumericVector asNv(SEXP in, const char* what) {
  if (TYPEOF(in) != REALSXP) {
    Rcpp::stop(_("'%s' needs to be a numeric vector"), what);
  }
  return as<NumericVector>(in);
}

static inline Nullable<NumericVector> asNNv(SEXP in, const char* what) {
  int type = TYPEOF(in);
  if (!(type == 0 || type == REALSXP)){
    Rcpp::stop(_("'%s' needs to be an numeric vector or NULL"), what);
  }
  return as<Nullable<NumericVector>>(in);
}

static inline IntegerVector asIv(SEXP in, const char* what) {
  if (!qtest(in, "x")) {
    Rcpp::stop(_("'%s' needs to be a integer vector"), what);
  }
  return as<IntegerVector>(in);
}

#endif
