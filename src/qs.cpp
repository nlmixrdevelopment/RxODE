#include <Rcpp.h>
#include <qs.h>


//[[Rcpp::export]]
Rcpp::CharacterVector rxQs(SEXP const x) {
  Rcpp::CharacterVector ret(1);
  ret[0] = qs::base91_encode(qs::c_qserialize(x, "high", "zstd", 4, 15, false));
  return ret;
}

//[[Rcpp::export]]
SEXP rxQr(const std::string& encoded_string) {
  return qs::c_qdeserialize(qs::base91_decode(encoded_string), true, false);
}
