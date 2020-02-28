// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>
using namespace Rcpp;

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif
bool rxIs(const RObject &obj, std::string cls);
Function getRxFn(std::string name);

std::string symengineRes(std::string val);

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
	std::string sDf = "rx__df_" + symengineRes(s1) + "_dy_" + symengineRes(s2) + "__";
	out4[i] = sDf;
	out5[i] = "assign(\"" + sDf + "\",with(model,D(rx__d_dt_" +
	  symengineRes(s1) + "__, \"" + symengineRes(s2) + "\")), envir=model)";
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
    stop(_("unanticipated input for rxExpandGrid_"));
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
    std::string curLine= "rx__d_dt_rx__sens_" + symengineRes(curState) + "_BY_" +
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
      tmp + ",\"" + curSens + "\")),envir=model)";
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
  CharacterVector s0(lenF);
  int i1, i2, i3;
  for (int i = lenF; i--;){
    i1 = i % len1;
    i2 = (int)(floor(i / len1)) % len2;
    i3 = floor(i / (len1*len2));
    std::string cS = as<std::string>(state[i1]);
    std::string cS1 = as<std::string>(s1[i2]);
    std::string cS2 = as<std::string>(s2[i3]);
    std::string sensSp = "rx__sens_" + symengineRes(cS) + "_BY_" +
      symengineRes(cS1) + "_BY_" + symengineRes(cS2) + "__";
    ddt[i] = "d/dt("+ sensSp + ")";
    std::string curLine= "rx__d_dt_rx__sens_" + symengineRes(cS) +
      "_BY_" + symengineRes(cS1) + "_BY_" + symengineRes(cS2) + "____";
    ddtS[i] = curLine;
    ddtS2[i] = sensSp;
    curLine = "assign(\"" + curLine + "\", with(model,";
    std::string v1 = "rx__d_dt_rx__sens_"+symengineRes(cS)+"_BY_"+
      symengineRes(cS2)+"____";
    curLine += "D("+v1+",\""+cS1+"\")";
    for (int j = len1; j--;){
      std::string s2 = as<std::string>(state[j]);
      curLine += "+D("+v1+",\""+s2+"\")*rx__sens_"+symengineRes(s2)+"_BY_"+symengineRes(cS1)+
	"__+rx__sens_"+symengineRes(s2)+"_BY_"+symengineRes(cS1)+"_BY_"+symengineRes(cS2)+
	"__*rx__df_"+
	symengineRes(cS) + "_dy_"+symengineRes(s2)+"__";
    }
    curLine += "),envir=model)";
    line[i] = curLine;
    state0r[i] = sensSp + "(0)";
    s0[i] = "rx_" + symengineRes(cS) + "_ini_0__";
    stateD[i] = "assign(\"rx_"+sensSp+"_ini_0__\",with(model,D(D(rx_" +
      symengineRes(cS) + "_ini_0__,\"" + symengineRes(cS1) + "\"),\""+symengineRes(cS2)+
      "\")),envir=model)";
  }
  List out(7);
  out[0] = ddt;
  out[1] = ddtS;
  out[2] = ddtS2;
  out[3] = line;
  out[4] = state0r;
  out[5] = stateD;
  out[6] = s0;
  out.attr("names") = CharacterVector::create("ddt","ddtS","ddS2","line","s0r","s0D","s0");
  out.attr("class") = "data.frame";
  out.attr("row.names") = IntegerVector::create(NA_INTEGER, -lenF);
  return out;  
}

//[[Rcpp::export]]
List rxExpandFEta_(CharacterVector state, int neta, int pred){
  CharacterVector fe(neta);
  CharacterVector calcS(neta);
  int nstate = state.size();
  for (int i = 0; i < neta; i++){
    std::string etaN = std::to_string(i+1);
    std::string feta;
    std::string calc;

    switch(pred){
    case 2:
      feta = "rx__sens_rx_pred__BY_ETA_" +
	etaN + "___";
      calc = "assign(\"" + feta + "\",with(.s,-D(rx_pred_, ETA_" +
	etaN + "_)";
      for (int j = nstate; j--;){
	calc += "-rx__sens_" + symengineRes(as<std::string>(state[j])) + "_BY_ETA_" + etaN +
	  "___*D(rx_pred_,\""+ symengineRes(as<std::string>(state[j])) + "\")";
      }
      calc += "), envir=.s)";
      break;
    case 1:
      feta = "rx__sens_rx_pred__BY_ETA_" +
	etaN + "___";
      calc = "assign(\"" + feta + "\",with(.s,D(rx_pred_, ETA_" +
	etaN + "_)";
      for (int j = nstate; j--;){
	calc += "+rx__sens_" + symengineRes(as<std::string>(state[j])) + "_BY_ETA_" + etaN +
	  "___*D(rx_pred_,"+ symengineRes(as<std::string>(state[j])) + ")";
      }
      calc += "), envir=.s)";
      break;
    case 0:
      feta = "rx__sens_rx_r__BY_ETA_" +
	etaN + "___";
      calc = "assign(\"" + feta + "\",with(.s,D(rx_r_,ETA_" +
	etaN + "_)";
      for (int j = nstate; j--;){
	calc += "+rx__sens_" + symengineRes(as<std::string>(state[j])) + "_BY_ETA_" + etaN +
	  "___*D(rx_r_,"+ symengineRes(as<std::string>(state[j])) + ")";
      }
      calc += "), envir=.s)";
      break;
    }
    fe[i] = feta;
    calcS[i] = calc;
  }
  List out(2);
  out[0] = fe;
  out[1] = calcS;
  out.attr("names") = CharacterVector::create("dfe","calc");
  out.attr("class") = "data.frame";
  out.attr("row.names") = IntegerVector::create(NA_INTEGER, -neta);
  return out;
}

//[[Rcpp::export]]
std::string rxRepR0_(int neta){
  std::string cur = "rx_r_";
  for (int j = neta; j--;){
    cur = "symengine::subs(" + cur+ ",ETA_" + std::to_string(j+1) + "_,..s0)";
  }
  cur = "assign(\"rx_r_\",with(.s,"+cur+"),envir=.s)";
  return cur;
}

List rxModelVars_(const RObject &obj);

// Expands nesting for theta/eta;
//
// @param thetaNest
// 
// @noRd
void rxExpandNestingRep(CharacterVector &thetaNest,
			CharacterVector &thetaNestTran,
			CharacterVector &thetaNestFull,
			int &thCnt, int &curtheta,
			List &aboveVars, NumericVector& above,
			std::string &retS, std::string &retF,
			std::string &namesAbove,
			List &data, std::string thetaVar = "THETA[") {
  std::string theta;
  int lastTheta;
  int firstTheta = curtheta;
  int nameTheta = curtheta;
  for (int j = 0; j < thetaNest.size(); ++j) {
    std::string curNest = as<std::string>(thetaNest[j]);
    CharacterVector nestVars = as<CharacterVector>(aboveVars[curNest]);
    RObject curNestV = as<RObject>(data[curNest]);
    CharacterVector curNestLvl = curNestV.attr("levels");
    int nnest = as<int>(above[curNest]);
    if (nnest > 1){
      retF += "rep(list(cvPost(omega$nu$" +curNest + ",omega$" +curNest +
	",type=\"invWishart\")),"+std::to_string(nnest)+")";
    } else {
      retF += "list(cvPost(omega$nu$" +curNest + ",omega$" +curNest +",type=\"invWishart\"))";
    }
    if (j != thetaNest.size()-1) retF += ",";
    // This is the base theta count
    lastTheta = thCnt;
    for (int i = 0; i < nestVars.size(); ++i) {
      std::string curPar = as<std::string>(nestVars[i]);
      retS += curPar + "=";
      for (int k = 0; k < nnest; ++k) {
	theta = thetaVar + std::to_string(lastTheta+i+k*nnest+firstTheta) + "]";
	namesAbove += "\"" +thetaVar + std::to_string(nameTheta++) + "]\"";
	namesAbove += ",";
	retS += "(" + curNest + "==" + std::to_string(k+1)+")*" + theta;
	thetaNestTran[thCnt] = curNest + "("+as<std::string>(curNestLvl[k])+")";
	thetaNestFull[thCnt] = theta;
	curtheta++; thCnt++;
	if (k != nnest-1) retS += "+";
	else retS += ";\n";
      }
    }
  }
  thetaNestTran.attr("names") = thetaNestFull;
}

RObject evalFun(std::string str, bool pretty = true){
  Function deparse2("deparse", R_BaseNamespace);
  Function parse2("parse", R_BaseNamespace);
  Function eval2("eval", R_BaseNamespace);
  RObject ret = eval2(parse2(_["text"]=str));
  if (pretty){
    CharacterVector retC = deparse2(ret);
    std::string str2 = "";
    for (int i = 0; i < retC.size(); ++i){
      str2+=as<std::string>(retC[i]) + "\n";
    }
    ret = eval2(parse2(_["text"]=str2));
  }
  return ret;
}

//[[Rcpp::export]]
List rxExpandNesting(const RObject& obj, List& nestingInfo,
		     bool compile=false){
  std::string retS="";
  std::string aboveF="";
  std::string namesBelow="";
  std::string namesAbove="";
  std::string belowF="";
  List mv = rxModelVars_(obj);
  IntegerVector flags = as<IntegerVector>(mv["flags"]);
  int cureta = as<int>(flags["maxeta"])+1;
  int curtheta = as<int>(flags["maxtheta"])+1;
  int extraTheta = as<int>(nestingInfo["extraTheta"]);
  int extraEta = as<int>(nestingInfo["extraEta"]);
  CharacterVector thetaNames(extraTheta);
  CharacterVector etaNames(extraEta);
  List aboveVars = as<List>(nestingInfo["aboveVars"]);
  NumericVector above = as<NumericVector>(nestingInfo["above"]);
  List data = as<List>(nestingInfo["data"]);
  CharacterVector thetaNest = aboveVars.attr("names");
  CharacterVector thetaNestTran(extraTheta);
  CharacterVector thetaNestFull(extraTheta);

  List belowVars = as<List>(nestingInfo["belowVars"]);
  NumericVector below = as<NumericVector>(nestingInfo["below"]);
  CharacterVector etaNest = belowVars.attr("names");
  CharacterVector etaNestTran(extraEta);
  CharacterVector etaNestFull(extraEta);
  int thCnt=0;
  std::string idName=nestingInfo["idName"];

  Function deparse2("deparse", R_BaseNamespace);
  CharacterVector dpO = deparse2(nestingInfo["omega"]);
  std::string dp0="";
  for (int i = 0; i < dpO.size(); ++i) {
    dp0 += as<std::string>(dpO[i]);
  }
  
  aboveF += ".theta <- as.matrix(Matrix::bdiag(list(";
  rxExpandNestingRep(thetaNest, thetaNestTran, thetaNestFull,
		     thCnt, curtheta,
		     aboveVars, above, retS, aboveF, namesAbove, data,
		     "THETA[");
  aboveF +=")));\n";
  
  namesBelow = ".dim <- c(dimnames(omega$" + idName + ")[[1]],";
  belowF += "function(omega=" + dp0 + "){";
  belowF += ".omega <- as.matrix(Matrix::bdiag(c(list(cvPost(omega$nu$" + idName + ",omega$" +
    idName + ",type=\"invWishart\"))";
  if (etaNest.size() > 0) belowF += ",";
  int etCnt = 0;
  rxExpandNestingRep(etaNest, etaNestTran, etaNestFull,
		     etCnt, cureta,
		     belowVars, below, retS, belowF, namesBelow, data,
		     "ETA[");
  belowF +=")));" + namesBelow + 
    "NULL);dimnames(.omega) <- list(.dim,.dim);return(.omega);}";
  
  CharacterVector mod = mv["model"];
  List ret(5);
  retS += as<std::string>(mod[0]);
  if (compile){
    Function rxode = getRxFn("RxODE");
    RObject retO = rxode(retS);
    ret[0] = retO;
  } else {
    ret[0] = retS;
  }
  ret[1] = thetaNestTran;
  ret[2] = etaNestTran;
  ret[3] = aboveF;
  ret[4] = evalFun(belowF, true);
  ret.attr("names") = CharacterVector::create("mod","theta","eta",
					      "aboveF", "belowF");
  return(ret);
}
