// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>
using namespace Rcpp;
bool rxIs(const RObject &obj, std::string cls);

//[[Rcpp::export]]
List rxExpandGrid_(RObject &c1, RObject &c2, RObject &type){
  if (rxIs(c1, "character") && rxIs(c2, "character")){
    CharacterVector in1 = as<CharacterVector>(c1);
    CharacterVector in2 = as<CharacterVector>(c2);
    int len1 = in1.size();
    int len2 = in2.size();
    int lenF = len1*len2;
    CharacterVector out1(lenF);
    CharacterVector out2(lenF);
    int i1, i2;
    int iType = as<int>(type);
    if (iType == 0){
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
    } else if (iType == 1) {
      CharacterVector out3(lenF);
      CharacterVector out4(lenF);
      CharacterVector out5(lenF);
      for (int i = lenF; i--;){
	i1 = i % len1;
	i2 = floor(i / len1);
	std::string s1 = as<std::string>(in1[i1]);
	std::string s2 = as<std::string>(in2[i2]);
	out1[i] = s1;
	out2[i] = s2;
	out3[i] = "df(" + s1 + ")/dy(" + s2 + ")";
	std::string sDf = "__d_df_" + s1 + "_dy_" + s2 + "__";
	out4[i] = sDf;
	out5[i] = sDf + " = diff(rx__d_dt_" + s1 + "__, " +
	  s2 +")";
      }
      List out(5);
      out[0] = out1;
      out[1] = out2;
      out[2] = out3;
      out[3] = out4;
      out[4] = out5;
      out.attr("class") = "data.frame";
      out.attr("row.names") = IntegerVector::create(NA_INTEGER, -lenF);
      out.attr("names") = CharacterVector::create("s1", "s2", "rx", "sym", "line");
      return out;
    }
  } else {
    stop("Unanticipated input for rxExpandGrid_");
  }
  List tmp = List(0);
  return tmp;
}
  
