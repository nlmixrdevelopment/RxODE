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
	std::string sDf = "rx__df_" + s1 + "_dy_" + s2 + "__";
	out4[i] = sDf;
	out5[i] = "assign(\"" + sDf + "\",with(model,D(rx__d_dt_" +
	  s1 + "__, " + s2 +")), envir=model)";
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

//[[Rcpp::export]]
List rxExpandSens_(CharacterVector state, CharacterVector calcSens){
  // rx__d_dt_rx__sens_V1_BY_V2____
  int len1 = state.size();
  int len2 = calcSens.size();
  int lenF = len1*len2;
  int i1, i2;
  CharacterVector ddt(lenF);
  CharacterVector ddtS(lenF);
  CharacterVector ddtS2(lenF);
  CharacterVector line(lenF);
  CharacterVector state0(lenF);
  CharacterVector state0r(lenF);
  CharacterVector stateD(lenF);
  CharacterVector rateS(lenF);
  CharacterVector rateR(lenF);
  CharacterVector durS(lenF);
  CharacterVector durR(lenF);
  CharacterVector lagS(lenF);
  CharacterVector lagR(lenF);
  CharacterVector FS(lenF);
  CharacterVector FR(lenF);
  for (int i = lenF; i--;){
    //"d/dt(rx__sens_%s_BY_%s__)"
    i1 = i % len1;
    i2 = floor(i / len1);
    std::string curState = as<std::string>(state[i1]);
    std::string curSens = as<std::string>(calcSens[i2]);
    std::string sensSp = "rx__sens_" + curState + "_BY_" + curSens + "__";
    ddt[i] = "d/dt("+ sensSp + ")";
    std::string curLine= "rx__d_dt_rx__sens_" + curState + "_BY_" +
      curSens + "____";
    ddtS[i] = curLine;
    ddtS2[i] = sensSp;
    curLine = "assign(\"" + curLine + "\", with(model,";
    for (int j = len1; j--;){
      std::string curState2=as<std::string>(state[j]);
      // sprintf("df(%s)/dy(%s)*rx__sens_%s_BY_%s__", s1, rxToSymPy(s2), s2, rxToSymPy(sns))
      curLine += "rx__df_"+curState+"_dy_"+curState2+
	"__*rx__sens_"+curState2+"_BY_"+curSens+"__+";
    }
    curLine += "rx__df_"+curState+ "_dy_"+curSens+"__),envir=model)";
    line[i] = curLine;
    std::string tmp = "rx_" + curState  + "_ini_0__";
    state0[i] = tmp;
    state0r[i] = sensSp + "(0)";
    stateD[i] = "assign(\"rx_"+sensSp+"_ini_0__\",with(model,D(" +
      tmp + "," + curSens + ")),envir=model)";
    rateS[i] = "rx_rate_" + curState + "_";
    rateR[i] = "rate(" + curState + ")";
    durS[i] = "rx_dur_" + curState + "_";
    durR[i] = "dur(" + curState + ")";
    lagS[i] = "rx_lag_" + curState + "_";
    lagR[i] = "lag(" + curState + ")";
    FS[i] = "rx_f_" + curState + "_";
    FR[i] = "f(" + curState + ")";
  }
  List out(15);
  out[0] = ddt;
  out[1] = ddtS;
  out[2] = line;
  out[3] = state0;
  out[4] = stateD;
  out[5] = state0r;
  out[6] = ddtS2;
  out[7] = rateS;
  out[8] = rateR;
  out[9] = durS;
  out[10] = durR;
  out[11] = lagS;
  out[12] = lagR;
  out[13] = FS;
  out[14] = FR;
  out.attr("names") = CharacterVector::create("ddt","ddtS","line","s0", "s0D","s0r", "ddS2",
					      "rateS","rateR", "durS","durR","lagS", "lagR",
					      "fS","fR");
  out.attr("class") = "data.frame";
  out.attr("row.names") = IntegerVector::create(NA_INTEGER, -lenF);
  return out;
}


//[[Rcpp::export]]
List rxExpandSens2_(CharacterVector state, CharacterVector s1, CharacterVector s2){
  // rx__d_dt_rx__sens_V1_BY_V2____
  int len1 = state.size();
  int len2 = s1.size();
  int len3 = s2.size();
  int lenF = len1*len2*len3;
  CharacterVector ddt(lenF);
  CharacterVector ddtS(lenF);
  CharacterVector ddtS2(lenF);
  CharacterVector line(lenF);
  CharacterVector state0r(lenF);
  CharacterVector stateD(lenF);
  int i1, i2, i3;
  for (int i = lenF; i--;){
    i1 = i % len1;
    i2 = (int)(floor(i / len1)) % len2;
    i3 = floor(i / (len1*len2));
    std::string cS = as<std::string>(state[i1]);
    std::string cS1 = as<std::string>(s1[i2]);
    std::string cS2 = as<std::string>(s2[i3]);
    std::string sensSp = "rx__sens_" + cS + "_BY_" + cS1 + "_BY_" + cS2 + "__";
    ddt[i] = "d/dt("+ sensSp + ")";
    std::string curLine= "rx__d_dt_rx__sens_" + cS + "_BY_" + cS1 + "_BY_" + cS2 + "____";
    ddtS[i] = curLine;
    ddtS2[i] = sensSp;
    curLine = "assign(\"" + curLine + "\", with(model,";
    std::string v1 = "rx__d_dt_rx__sens_"+cS+"_BY_"+cS2+"____";
    curLine += "D("+v1+","+cS1+")";
    for (int j = len1; j--;){
      std::string s2 = as<std::string>(state[j]);
      curLine += "+D("+v1+","+s2+")*rx__sens_"+s2+"_BY_"+cS1+
	"__+rx__sens_"+s2+"_BY_"+cS1+"_BY_"+cS2+"__*rx__df_"+
	cS + "_dy_"+s2+"__";
    }
    curLine += "),envir=model)";
    line[i] = curLine;
    state0r[i] = sensSp + "(0)";
    stateD[i] = "assign(\"rx_"+sensSp+"_ini_0__\",with(model,D(D(rx_" +
      cS + "_ini_0__," + cS1 + "),"+cS2+")),envir=model)";
  }
  List out(6);
  out[0] = ddt;
  out[1] = ddtS;
  out[2] = ddtS2;
  out[3] = line;
  out[4] = state0r;
  out[5] = stateD;
  out.attr("names") = CharacterVector::create("ddt","ddtS","ddS2","line","s0r","s0D");
  out.attr("class") = "data.frame";
  out.attr("row.names") = IntegerVector::create(NA_INTEGER, -lenF);
  return out;  
}
