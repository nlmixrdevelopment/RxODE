#include <Rcpp.h>

Rcpp::Function loadNamespaceQs("loadNamespace", R_BaseNamespace);
Rcpp::Environment qsNs;
bool loadQsC = false;

static void loadQs() {
  if (!loadQsC) {
    qsNs = loadNamespaceQs("qs");
    loadQsC = true;
  }
}

//[[Rcpp::export]]
Rcpp::CharacterVector rxQs(SEXP const x) {
  Rcpp::CharacterVector ret(1);
  loadQs();
  Rcpp::Function base91_encode = Rcpp::as<Rcpp::Function>(qsNs["base91_encode"]);
  Rcpp::Function qserialize = Rcpp::as<Rcpp::Function>(qsNs["qserialize"]);
  ret[0] = base91_encode(qserialize(x, "high", "zstd", 22, 15, false));
  return ret;
}

//[[Rcpp::export]]
SEXP rxQr(const std::string& encoded_string) {
  loadQs();
  Rcpp::Function base91_decode = Rcpp::as<Rcpp::Function>(qsNs["base91_decode"]);
  Rcpp::Function qdeserialize = Rcpp::as<Rcpp::Function>(qsNs["qdeserialize"]);
  return qdeserialize(base91_decode(encoded_string), true, false);
}
