// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>
using namespace Rcpp;
bool rxIs(const RObject &obj, std::string cls);

//[[Rcpp::export]]
List rxExpandGrid_(RObject &c1, RObject &c2){
  if (rxIs(c1, "character") && rxIs(c2, "character")){
    CharacterVector in1 = as<CharacterVector>(c1);
    CharacterVector in2 = as<CharacterVector>(c2);
    int len1 = in1.size();
    int len2 = in2.size();
    int lenF = len1*len2;
    CharacterVector out1(lenF);
    CharacterVector out2(lenF);
    int i1, i2;
    for (int i = lenF; i--;){
      i1 = i % len1;
      i2 = floor(i / len1);
      out1[i] = in1[i1];
      out2[i] = in2[i2];
    }
    List out(2);
    out[0] = out1;
    out[1] = out2;
    out.attr("class") = "data.frame";
    out.attr("row.names") = IntegerVector::create(NA_INTEGER, -lenF);
    out.attr("names") = CharacterVector::create("Var1","Var2");
    return out;
  } else {
    stop("Unanticipated input for rxExpandGrid_");
  }
}
  
