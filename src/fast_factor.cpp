// Author: Kevin Ushey, David Cooley

#include <Rcpp.h>
using namespace Rcpp;
#include <checkmate.h>
#include "../inst/include/RxODE_as.h"
//
// The unique ordered by occurrence comes from David Cooley:
// It is found https://stackoverflow.com/questions/44697544/rcpp-unique-order-output
template < typename T, int RTYPE >
inline SEXP sexp_unique( Rcpp::Vector< RTYPE > x ) {
  std::set< T > seen;
  auto newEnd = std::remove_if( x.begin(), x.end(), [&seen]( const T value ) {
						      if ( seen.find( value ) != std::end( seen ) ) {
							return true;
						      }
						      seen.insert( value );
						      return false;
						    });
  x.erase( newEnd, x.end() );
  return x;
}

// returns unique values in their original input order
inline SEXP get_sexp_unique( SEXP s ) {

  SEXP s2 = Rcpp::clone( s );

  switch( TYPEOF( s2 ) ) {
  case LGLSXP: {
    return sexp_unique< bool, LGLSXP >( s2 );
  }
  case REALSXP: {
    return sexp_unique< double, REALSXP >( s2 );
  }
  case INTSXP: {
    return sexp_unique< int, INTSXP >( s2 );
  }
  case STRSXP: {
    return sexp_unique< char* , STRSXP >( s2 );
  }
  default: Rcpp::stop("unknown vector type");
  }
  return 0;
}

// adapted from https://gallery.rcpp.org/articles/fast-factor-generation/
// This was modified to use symbols for levels and class.
template <int RTYPE>
SEXP fast_factor_unsorted( const Vector<RTYPE>& x, SEXP oldLvl) {
  Vector<RTYPE> levs = get_sexp_unique(x); 
  IntegerVector out = match(x, levs);
  SEXP outS = wrap(out);
  SEXP lvl = R_NilValue;
  SEXP fac = wrap(CharacterVector("factor"));
  if (Rf_isNull(oldLvl)) {
    lvl = wrap(as<CharacterVector>(levs));
  } else {
    // RTYPE should be INTSXP
    SEXP levsSEXP = wrap(levs);
    IntegerVector lvlI = as<IntegerVector>(levsSEXP);
    int hasNa = 0;
    for (int i = lvlI.size(); i--;) {
      if (lvlI[i] == NA_INTEGER) {
	hasNa = 1;
	break;
      }
    }
    CharacterVector lvlC(lvlI.size()-hasNa);
    lvl = wrap(lvlC);
    int cur=0;
    int j=0;
    for (int i = 0; i < lvlI.size(); ++i) {
      cur = lvlI[i];
      if (cur != NA_INTEGER) {
	SET_STRING_ELT(lvl, j++, STRING_ELT(oldLvl, lvlI[i]-1));
      }
    }
  }
  Rf_setAttrib(outS, R_LevelsSymbol, lvl);
  Rf_setAttrib(outS, R_ClassSymbol, fac);
  return outS;
}
//[[Rcpp::export]]
SEXP convertId_(SEXP x) {
  SEXP oldLvl = R_NilValue;
  switch( TYPEOF(x) ) {
  case INTSXP: {
    oldLvl = Rf_getAttrib(x, R_LevelsSymbol);
    return fast_factor_unsorted<INTSXP>(x, oldLvl);
  }
  case REALSXP: return fast_factor_unsorted<REALSXP>(x, oldLvl);
  case STRSXP: return fast_factor_unsorted<STRSXP>(x, oldLvl);
  }
  return R_NilValue;
}

