#include <RcppArmadillo.h>
#include <algorithm>
#include "../inst/include/RxODE.h"

#define rxModelVars(a) rxModelVars_(a)
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
using namespace Rcpp;

List rxModelVars_(const RObject &obj);
bool rxIs(const RObject &obj, std::string cls);
Environment RxODEenv();

IntegerVector curDvid;

IntegerVector toCmt(RObject inCmt, CharacterVector state, bool isDvid=false,
		    int stateSize = 0, int sensSize=0){
  RObject cmtInfo = R_NilValue;
  List extraCmt;
  if (rxIs(inCmt, "numeric") || rxIs(inCmt, "integer")){
    if (rxIs(inCmt, "factor")){
      CharacterVector lvl = inCmt.attr("levels");
      IntegerVector lvlI(lvl.size());
      int i, j, k=0;
      std::string curLvl, curState,negSub;
      bool foundState = false;
      bool isNeg = false;
      for (i = 0; i < lvlI.size(); i++){
	curLvl = as<std::string>(lvl[i]);
	negSub = curLvl.substr(0,1);
	if (negSub == "-"){
	  isNeg = true;
	  curLvl = curLvl.substr(1, std::string::npos);
	} else {
	  isNeg = false;
	}
	foundState=false;
	for (j = state.size(); j--;){
	  curState = as<std::string>(state[j]);
	  if (curState == curLvl){
	    if (isNeg){
	      lvlI[i] = -j-1;
	    } else {
	      lvlI[i] = j+1;
	    }
	    foundState = true;
	    break;
	  }
	}
	if (!foundState){
	  if (curLvl == "(default)" || curLvl == "(obs)"){
	    lvlI[i] = 1;
	  } else {
	    k++;
	    if (isNeg){
	      stop("Negative compartments on non-ode cmt (%s) do not make sense.", curLvl.c_str());
	    } else {
	      List tmpList(extraCmt.size()+1);
	      for (int i = extraCmt.size(); i--;) tmpList[i] = extraCmt[i];
	      tmpList[extraCmt.size()]= curLvl;
	      extraCmt = tmpList;
	      lvlI[i] = state.size() + k;
	    }
	  }
	}
      }
      IntegerVector cmtIn = IntegerVector(inCmt);
      IntegerVector ret(cmtIn.size());
      for (j=ret.size(); j--;){
	if (IntegerVector::is_na(cmtIn[j])){
	  ret[j] = NA_INTEGER;
	} else {
	  ret[j] = lvlI[cmtIn[j]-1];
	}
      }
      CharacterVector newCmt(state.size()+extraCmt.size());
      for (int j = state.size(); j--;) newCmt[j]=state[j];
      for (int j = extraCmt.size(); j--;)
	newCmt[j+state.size()] = as<std::string>(extraCmt[j]);
      ret.attr("cmtNames") = newCmt;
      return ret;
    } else {
      if (isDvid){
	// This converts DVID to cmt; Things that don't match become -9999
	Environment rx = RxODEenv();
	Function convertDvid = rx[".convertDvid"];
	IntegerVector in = convertDvid(inCmt, curDvid.length());
	IntegerVector out(in.size());
	IntegerVector conv = curDvid;
	std::vector<int> warnDvid;
	std::vector<int> warnConvertDvid;
	std::string warnC = "'dvid' were not numbered 1, 2, etc\nThey were converted from->to:";
	IntegerVector cmtIn = IntegerVector(inCmt);
	int curConv;
	for (int i = in.size(); i--;){
	  if (in[i] != cmtIn[i]){
	    if (std::find(warnConvertDvid.begin(), warnConvertDvid.end(), in[i]) == warnConvertDvid.end()){
	      warnConvertDvid.push_back(in[i]);
	      warnC = warnC + " " + std::to_string(cmtIn[i]) + "->" + std::to_string(in[i]);
	    }
	  }
	  if (in[i] > 0 && in[i] <= conv.size()){
	    curConv= conv[in[i]-1];
	    if (curConv > 0){
	      out[i] = curConv;
	    } else {
	      out[i] = curConv +state.size()+1;
	    }
	  } else {
	    if (std::find(warnDvid.begin(), warnDvid.end(), in[i]) == warnDvid.end()){
	      warnDvid.push_back(in[i]);
	    }
	    out[i] = -9999;
	  }
	}
	if (warnDvid.size() > 1){
	  std::string warn = "Undefined 'dvid' integer values in data: ";
	  std::sort(warnDvid.begin(), warnDvid.end());
	  for (int i = 0; i < (int)(warnDvid.size()-1); i++){
	    warn = warn + std::to_string(warnDvid[i]) + ", ";
	  }
	  warn = warn + std::to_string(warnDvid[warnDvid.size()-1]);
	  warning(warn);
	}
	if (warnConvertDvid.size() > 0){
	  warning(warnC);
	}
	return out;
      } else {
	IntegerVector in = as<IntegerVector>(inCmt);
	IntegerVector out(in.size());
	int baseSize = stateSize - sensSize;
	// Sensitivity equations are ignored in CMT data items.
	for (int i = in.size(); i--;){
	  if (in[i] > 0 && in[i] <= baseSize){
	    out[i] = in[i];
	  } else if (in[i] > 0) {
	    out[i] = in[i]+sensSize;
	  } else if (in[i] < 0 && in[i] >= -baseSize){
	    out[i] = in[i];
	  } else if (in[i] < 0){
	    out[i] = in[i] - sensSize;
	  } else {
	    out[i] = 0;
	  }
	}
	return out;
      }
    }
  } else if (rxIs(inCmt, "character")) {
    CharacterVector iCmt = as<CharacterVector>(inCmt);
    std::vector<int> newCmt;
    newCmt.reserve(iCmt.size());
    std::string strCmt, negSub;
    bool foundState=false;
    bool isNeg = false;
    int i, j, k = 0;
    for (i = 0; i < iCmt.size(); i++){
      strCmt = as<std::string>(iCmt[i]);
      negSub = strCmt.substr(0,1);
      if (negSub == "-"){
	isNeg = true;
	strCmt = strCmt.substr(1, std::string::npos);
      } else {
	isNeg = false;
      }
      foundState=false;
      if (strCmt == "(default)" || strCmt == "(obs)" || CharacterVector::is_na(iCmt[i])){
	foundState=true;
	newCmt.push_back(1);
      } else {
	for (j = state.size(); j--;){
	  if (as<std::string>(state[j]) == strCmt){
	    foundState = true;
	    if (isNeg){
	      newCmt.push_back(-j-1);
	    } else {
	      newCmt.push_back(j+1);
	    }
	    break;
	  }
	}
	if (!foundState){
	  for (j = k; j--;){
	    CharacterVector cur = extraCmt[j];
	    if (as<std::string>(cur) == strCmt){
	      foundState = true;
	      if (isNeg){
		stop("Negative compartments on non-ode cmt (%s) do not make sense.", strCmt.c_str());
	      } else {
		newCmt.push_back(state.size()+j+1);
	      }
	      break;
	    }
	  }
	  if (!foundState){
	    if (isNeg){
	      stop("Negative compartments on non-ode cmt (%s) do not make sense.", strCmt.c_str());
	    } else {
	      List tmpList(extraCmt.size()+1);
	      for (int i = extraCmt.size(); i--;) tmpList[i] = extraCmt[i];
	      extraCmt = tmpList;	    
	      newCmt.push_back(state.size()+k+1);
	      extraCmt[k++] = strCmt;
	    }
	  }
	}
      }
      CharacterVector newCmt(state.size()+extraCmt.size());
      for (int j = state.size(); j--;) newCmt[j]=state[j];
      for (int j = extraCmt.size(); j--;) newCmt[j+state.size()] = as<std::string>(extraCmt[j]);
      cmtInfo = as<RObject>(newCmt);
    }
    IntegerVector ret = wrap(newCmt);
    ret.attr("cmtNames") = cmtInfo;
    return ret;
  }
  stop("Should not reach here.");
  return IntegerVector::create(0);
}

bool _ini0=true;
//' Set Initial conditions to time zero instead of the first observed/dosed time
//'
//' @param ini0 When TRUE (default), set initial conditions to time
//'   zero. Otherwise the initial conditions are the first observed
//'   time.
//'
//' @export
//[[Rcpp::export]]
bool rxSetIni0(bool ini0 = true){
  _ini0=ini0;
  return _ini0;
}

extern void setFkeep(List keep);

//' Event translation for RxODE
//'
//' @param inData Data frame to translate
//' @param obj Model to translate data 
//' @param addCmt Add compartment to data frame (default code{FALSE}).
//' @param dropUnits Boolean to drop the units (default \code{FALSE}).
//' @param allTimeVar Treat all covariates as if they were time-varying
//' @param keepDosingOnly keep the individuals who only have dosing records and any
//'   trailing dosing records after the last observation.
//' @param combineDvid is a boolean indicating if RxODE will use DVID on observation
//'     records to change the cmt value; Useful for multiple-endpoint nlmixr models.  By default
//'     this is determined by code{option("RxODE.combine.dvid")} and if the option has not been set,
//'     this is \code{TRUE}. This typically does not affect RxODE simulations.
//' @param keep This is a named vector of items you want to keep in the final RxODE dataset.
//'     For added RxODE event records (if seen), last observation carried forward will be used.
//' @return Object for solving in RxODE
//' @keywords internal
//' @export
//[[Rcpp::export]]
List etTrans(List inData, const RObject &obj, bool addCmt=false,
	     bool dropUnits=false, bool allTimeVar=false,
	     bool keepDosingOnly=false, Nullable<LogicalVector> combineDvid=R_NilValue,
	     CharacterVector keep = CharacterVector(0)){
  Environment rx = RxODEenv();
  bool combineDvidB = false;
  if (!combineDvid.isNull()){
    combineDvidB = (as<LogicalVector>(combineDvid))[1];
  } else {
    Environment b=Rcpp::Environment::base_namespace();
    Function getOption = b["getOption"];
    combineDvidB = as<bool>(getOption("RxODE.combine.dvid", true));
  }
  List mv = rxModelVars_(obj);
  curDvid = clone(as<IntegerVector>(mv["dvid"]));
  CharacterVector trans = mv["trans"];
  if (rxIs(inData,"rxEtTran")){
    CharacterVector cls = inData.attr("class");
    List e0 = cls.attr(".RxODE.lst");
    if (as<std::string>(trans["lib.name"]) == as<std::string>(e0["lib.name"])){
      if (as<bool>(e0["allTimeVar"]) && !allTimeVar){
	LogicalVector sub0 = as<LogicalVector>(e0["sub0"]);
	int baseSize = as<int>(e0["baseSize"]);
	int nTv = as<int>(e0["nTv"]);
	List lst = as<List>(e0["lst"]);
	CharacterVector nme = as<CharacterVector>(e0["nme"]);
	List e = clone(e0);
	e["baseSize"] = R_NilValue;
	e["nTv"] = R_NilValue;
	e["lst"] = R_NilValue;
	e["nme"] = R_NilValue;	
	e["sub0"] = R_NilValue;
	e["allTimeVar"] = false;
	cls.attr(".RxODE.lst") = e;
	List lstF = List(baseSize+nTv);
	CharacterVector nmeF = CharacterVector(baseSize+nTv);
	int j=0;
	for (int i = 0; i < lst.size();i++){
	  if (sub0[i]){
	    lstF[j]=lst[i];
	    nmeF[j]=nme[i];
	    j++;
	  }
	}
	lstF.attr("names") = nmeF;
	lstF.attr("class") = cls;
	IntegerVector tmp = lstF[0];
	lstF.attr("row.names") = IntegerVector::create(NA_INTEGER,-tmp.size());
	return(lstF);
      } else {
	return inData;
      }
    }
    // stop("This dataset was prepared for another model.");
  }
  // Translates events + model into translated events
  CharacterVector dName = as<CharacterVector>(inData.attr("names"));
  CharacterVector lName = clone(dName);

  int i, idCol = -1, evidCol=-1, timeCol=-1, amtCol=-1, cmtCol=-1,
    dvCol=-1, ssCol=-1, rateCol=-1, addlCol=-1, iiCol=-1, durCol=-1, j,
    mdvCol=-1, dvidCol=-1;
  std::string tmpS;
  
  CharacterVector pars = as<CharacterVector>(mv["params"]);
  std::vector<int> covCol;
  std::vector<int> covParPos;
  std::vector<int> keepCol;
  std::string tmpS0;
  bool allBolus = true;
  bool allInf = true;
  int mxCmt = 0;
  for (i = lName.size(); i--;){
    tmpS0= as<std::string>(lName[i]);
    tmpS = as<std::string>(lName[i]);
    std::transform(tmpS.begin(), tmpS.end(), tmpS.begin(), ::tolower);
    lName[i] = tmpS;
    if (tmpS == "id") idCol=i;
    else if (tmpS == "evid") evidCol=i;
    else if (tmpS == "time") timeCol=i;
    else if (tmpS == "amt") amtCol=i;
    else if (tmpS == "cmt") cmtCol=i;
    else if (tmpS == "dv") dvCol=i;
    else if (tmpS == "ss")   ssCol=i;
    else if (tmpS == "rate") rateCol=i;
    else if (tmpS == "dur") durCol=i;
    else if (tmpS == "addl") addlCol=i;
    else if (tmpS == "ii")   iiCol=i;
    else if (tmpS == "mdv") mdvCol=i;
    else if (tmpS == "dvid") dvidCol=i;
    for (j = keep.size(); j--;){
      if (as<std::string>(dName[i]) == as<std::string>(keep[j])){
	if (tmpS == "evid") stop("Cannot keep 'evid'; try 'addDosing'");
	keepCol.push_back(i);
	break;
      }
    }
    for (j = pars.size(); j--;){
      // Check lower case
      if (tmpS == as<std::string>(pars[j])){
	// Covariate found.
	covCol.push_back(i);
	covParPos.push_back(j);
	break;
      }
      if (tmpS0 == as<std::string>(pars[j])){
	// Covariate found.
	covCol.push_back(i);
	covParPos.push_back(j);
	break;
      }
      // Check upper case.
      std::transform(tmpS.begin(), tmpS.end(), tmpS.begin(), ::toupper);
      if (tmpS == as<std::string>(pars[j])){
	// Covariate found.
	covCol.push_back(i);
	covParPos.push_back(j);
	break;
      }
    }
  }
  List covUnits(covCol.size());
  CharacterVector covUnitsN(covCol.size());
  NumericVector nvTmp, nvTmp2;
  bool hasCmt = false;
  for (i = covCol.size(); i--;){
    covUnitsN[i] = lName[covCol[i]];
    nvTmp2 = NumericVector::create(1.0);
    if (as<std::string>(lName[covCol[i]]) != "cmt"){
      nvTmp = as<NumericVector>(inData[covCol[i]]);
      if (!dropUnits && rxIs(nvTmp, "units")){
	nvTmp2.attr("class") = "units";
	nvTmp2.attr("units") = nvTmp.attr("units");
      }
    } else {
      hasCmt=true;
    }
    covUnits[i] = nvTmp2;
  }
  covUnits.attr("names") = covUnitsN;
  // EVID = 0; Observations
  // EVID = 1; is illegal, but converted from NONMEM
  // EVID = 2; Non-observation, possibly covariate
  // EVID = 3; Reset ODE states to zero; Non-observation event
  // EVID = 4; Reset and then dose event;  Illegal
  // EVID = 9; Non-observation event to ini system at time zero; This is to set the INIs at the correct place.
  // EVID = 10-99; mtime events (from ODE system)
  // When EVID > 100
  // EVID: ## # ## ##
  //       c2 I c1 xx
  // c2 = Compartment numbers over 100
  //  I = Infusion Flag
  //      0 = no Infusion
  //      1 = Infusion, AMT=rate (mg/hr for instance)
  //      2 = Infusion, duration is fixed
  //      9 = Rate is modeled, AMT=dose; Duration = AMT/(Modeled Rate) NONMEM RATE=-1
  //      8 = Duration is modeled, AMT=dose; Rate = AMT/(Modeled Duration) NONMEM RATE=-2
  //      7 = Turn off modeled rate compartment
  //      6 = Turn off modeled duration
  // c1 = Compartment numbers below 99
  // xx = 1, regular event
  // xx = 10, steady state event SS=1
  // xx = 20, steady state event + last observed info.
  // xx = 30, Turn off compartment
  // Steady state events need a II data item > 0
  
  CharacterVector state0 = as<CharacterVector>(mv["state"]);
  CharacterVector stateE = as<CharacterVector>(mv["stateExtra"]);
  CharacterVector stateS = as<CharacterVector>(mv["sens"]);
  int extraCmt  = as<int>(mv["extraCmt"]);
  // Enlarge compartments
  if (extraCmt == 2){
    CharacterVector newState(state0.size()+2);
    for (int j = state0.size();j--;) newState[j] = state0[j];
    newState[state0.size()] = "depot";
    newState[state0.size()+1] = "central";
    state0 = newState;
  } else if (extraCmt==1){
    CharacterVector newState(state0.size()+1);
    for (int j = state0.size();j--;) newState[j] = state0[j];
    newState[state0.size()] = "central";
    state0 = newState;
  }
  CharacterVector state(state0.size() + stateE.size());
  for (int i = 0; i < state0.size(); i++){
    state[i] = state0[i];
  }
  for (int i = 0; i < stateE.size(); i++){
    state[i+state0.size()] = stateE[i];
  }
  // Now adjust the compartment numbers if needed.
  int baseSize = state0.size() - stateS.size();
  for (int i = curDvid.size(); i--;){
    if (curDvid[i] > baseSize){
      curDvid[i] = stateS.size()+curDvid[i];
    }
  }
  std::vector<int> id;
  std::vector<int> allId;
  std::vector<int> obsId;
  std::vector<int> zeroId;
  std::vector<int> doseId;
  std::vector<int> evid;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> ii;
  std::vector<int> idx;
  std::vector<int> cmtF;
  std::vector<int> dvidF;
  std::vector<double> dv;
  std::vector<int> idxO;
  if (timeCol== -1){
    stop("time is required in dataset.");
  }
  NumericVector inTime;
  if (rxIs(inData[timeCol], "numeric") || rxIs(inData[timeCol], "integer")){
    inTime = as<NumericVector>(inData[timeCol]);
  } else {
    List newInData = clone(inData);
    Function convDate = rx[".convertExtra"];
    newInData =  convDate(newInData);
    return etTrans(newInData, obj, addCmt, dropUnits, allTimeVar, keepDosingOnly, combineDvid);
  }
  // save units information
  bool addTimeUnits = false;
  RObject timeUnits;
  if (rxIs(inTime, "units")){
    addTimeUnits=true;
    timeUnits=inTime.attr("units");
  }
  IntegerVector inCmt;
  RObject cmtInfo = R_NilValue;
  if (cmtCol != -1){
    inCmt = as<IntegerVector>(toCmt(inData[cmtCol], state, false,
				    state0.size(), stateS.size()));//as<IntegerVector>();
    cmtInfo = inCmt.attr("cmtNames");
    inCmt.attr("cmtNames") = R_NilValue;
  }
  IntegerVector inDvid;
  if (dvidCol != -1){
    inDvid = as<IntegerVector>(toCmt(inData[dvidCol], state, true,
				     state0.size(), stateS.size()));//as<IntegerVector>();
    inDvid.attr("cmtNames") = R_NilValue;
  }
  int tmpCmt = 1;
  IntegerVector inId;
  CharacterVector idLvl;
  if (idCol != -1){
    Function convId = rx[".convertId"];
    inId = convId(inData[idCol]);//as<IntegerVector>();
    idLvl = inId.attr("levels");
  } else {
    // warning("ID=1 added to dataset");
    idLvl = CharacterVector::create("1");
  }
  IntegerVector inSs;
  if (ssCol != -1){
    if (rxIs(inData[ssCol], "integer") || rxIs(inData[ssCol], "numeric") ||
	rxIs(inData[ssCol], "logical")){
      // NA by default is NA_logical
      inSs = as<IntegerVector>(inData[ssCol]);
    } else {
      stop("Steady state column (ss) needs to be an integer");
    }
  }
  IntegerVector inEvid;
  if (evidCol != -1){
    if (rxIs(inData[evidCol], "integer") || rxIs(inData[evidCol], "numeric") ||
	rxIs(inData[evidCol], "logical")){
      inEvid = as<IntegerVector>(inData[evidCol]);
    } else {
      stop("Event id (evid) needs to be an integer");
    }
  }
  IntegerVector inMdv;
  if (mdvCol != -1){
    if (rxIs(inData[mdvCol], "integer") || rxIs(inData[mdvCol], "numeric") ||
	rxIs(inData[mdvCol], "logical")){
      inMdv = as<IntegerVector>(inData[mdvCol]);
    } else {
      stop("Missing dependent variable (mdv) needs to be an integer");
    }
  }
  NumericVector inRate;
  if (rateCol != -1){
    if (rxIs(inData[rateCol], "integer") || rxIs(inData[rateCol], "numeric") ||
	rxIs(inData[rateCol], "logical")){
      inRate = as<NumericVector>(inData[rateCol]);
    } else {
      stop("'rate' needs to be a number");
    }
  }

  NumericVector inDur;
  if (durCol != -1){
    if (rxIs(inData[durCol], "integer") || rxIs(inData[durCol], "numeric") ||
	rxIs(inData[durCol], "logical")){
      inDur = as<NumericVector>(inData[durCol]);
    } else {
      stop("'dur' needs to be a number");
    }
  }
  
  bool addAmtUnits = false;
  RObject amtUnits;
  NumericVector inAmt;
  if (amtCol != -1){
    if (rxIs(inData[amtCol], "integer") || rxIs(inData[amtCol], "numeric") ||
	rxIs(inData[amtCol], "logical")){
      inAmt = as<NumericVector>(inData[amtCol]);
      if (rxIs(inAmt, "units")){
	addAmtUnits=true;
	amtUnits=inAmt.attr("units");
      }
    } else {
      stop("Amount (amt) needs to be a number");
    }
  }
  NumericVector inIi;
  if (iiCol != -1){
    if (rxIs(inData[iiCol], "integer") || rxIs(inData[iiCol], "numeric") ||
	rxIs(inData[iiCol], "logical")){
      inIi = as<NumericVector>(inData[iiCol]);
    } else {
      stop("Inter-dose interval (ii) needs to be a number.");
    }
  }
  IntegerVector inAddl;
  if (addlCol != -1){
    if (rxIs(inData[addlCol], "integer") || rxIs(inData[addlCol], "numeric")||
	rxIs(inData[addlCol], "level")){
      inAddl = as<IntegerVector>(inData[addlCol]);
    } else {
      stop("Number of additional doses (addl) needs to be an integer");
    }
  }
  NumericVector inDv;
  if (dvCol != -1){
    if (rxIs(inData[dvCol], "integer") || rxIs(inData[dvCol], "numeric") ||
	rxIs(inData[dvCol], "logical")){
      inDv = as<NumericVector>(inData[dvCol]);
    } else {
      stop("The dependent variable (dv) needs to be a number");
    }
  }
  int flg = 0;
  int cid = 0;
  int nMtime = as<int>(mv["nMtime"]);
  double rate = 0.0;
  int nid=0;
  int cmt = 0;
  int rateI = 0;
  int cmt100; //= amt[i]/100;
  int cmt99;  //= amt[i]-amt100*100;
  int cevid;
  int nevid;
  int caddl;
  double ctime;
  double cii;
  double dur =0.0;
  double camt;
  int curIdx=0;
  double cdv;
  int nobs=0, ndose=0;
  for (int i = 0; i < inTime.size(); i++){
    if (idCol == -1) cid = 1;
    else cid = inId[i];
    if (dvCol == -1) cdv = NA_REAL;
    else cdv = inDv[i];
    ctime=inTime[i];
    if (std::isinf(ctime)){
      stop("Infinite times are not allowed");
    }
    if (ctime < 0){
      stop("Negative times are not allowed");
    }
    if (iiCol == -1) cii = 0;
    else cii = inIi[i];
    
    if (std::find(allId.begin(), allId.end(), cid) == allId.end()){
      allId.push_back(cid);
      // New ID
      // Add mtime records
      for (j = nMtime; j--;){
	id.push_back(cid);
	evid.push_back(j+10);
	cmtF.push_back(0);
	time.push_back(0.0);
	amt.push_back(NA_REAL);
	ii.push_back(0.0);
	dv.push_back(NA_REAL);
	idx.push_back(-1);
	idxO.push_back(curIdx);curIdx++;
      }
      nid++;
    }
    
    // SS flag
    flg=1;
    if (ssCol == -1) flg=1;
    else if (inSs[i] == 0) flg=1;
    else if (IntegerVector::is_na(inSs[i])) flg=1;
    else if (inSs[i] == 1 && cii > 0) flg=10;
    else if (inSs[i] == 2 && cii > 0) flg=20;
    if (cmtCol != -1){
      tmpCmt = inCmt[i];
      if (inCmt[i] == 0){
	if (evidCol == -1){
	  tmpCmt=1;
	} else if (inEvid[i] == 0){
	  tmpCmt=1;
	} else {
	  tmpCmt=1;
	}
      }      
      if (IntegerVector::is_na(inCmt[i])){
	tmpCmt = 1;
      } else if (inCmt[i] < 0){
	if (flg != 1) stop("Steady state records cannot be on negative compartments.");
	flg = 30;
	tmpCmt = -tmpCmt;
      }
      cmt = tmpCmt;
    } else if (cmtCol == -1) {
      cmt = 1;
    }
    // CMT flag
    else cmt = tmpCmt;
    if (cmt <= 99){
      cmt100=0;
      cmt99=cmt;
    } else {
      cmt100=cmt/100; 
      cmt99=cmt-cmt100*100; 
    }
    mxCmt = max2(cmt,mxCmt);

    // Amt
    if (amtCol == -1) camt = 0.0;
    else camt = inAmt[i];
    
    rateI = 0;
    // Rate
    if (durCol == -1 || inDur[i] == 0 || ISNA(inDur[i])){
      if (rateCol == -1 || inRate[i] == 0 || ISNA(inRate[i])) rate = 0.0;
      else rate = inRate[i];
      if (rate == -1.0){
	// rate is modeled
	rateI = 9;
      } else if (rate == -2.0){
	// duration is modeled
	rateI = 8;
      } else if (rate > 0){
	// Rate is fixed
	if (evidCol == -1 || inEvid[i] == 1 || inEvid[i] == 4){
	  rateI = 1;
	} else {
	  rateI = 0;
	  rate = 0.0;
	}
      }
    } else if (rateCol == -1 || inRate[i] == 0 || ISNA(inRate[i])) {
      if (durCol == -1) rate = 0.0;
      if (inDur[i] == 0) rate = 0;
      // if (inDur[i] > 0)
      if (inDur[i] == -1.0){
	// rate is modeled
	rateI = 9;
      } else if (inDur[i] == -2.0){
	// duration is modeled
	rateI = 8;
      } else if (inDur[i] > 0){
	// Duration is fixed
	if (evidCol == -1 || inEvid[i] == 1 || inEvid[i] == 4){
	  rateI = 2;
	  rate = camt/inDur[i];
	} else if (inEvid[i] > 4){
	  rateI=0;
	  rate = 0.0;
	}  
      }
    } else {
      stop("'rate' and/or 'dur' are not specified correctly.");
    }
    if (addlCol == -1) caddl=0;
    else caddl = inAddl[i];
    // EVID flag
    if (evidCol == -1){
      // Missing EVID
      if (rateI == 0 && (ISNA(camt) || camt == 0.0)){
	cevid = 0;
	if (mdvCol != -1 && inMdv[i] == 1){
	  cevid=2;
	}
	if (std::find(obsId.begin(), obsId.end(), cid) == obsId.end()){
	  obsId.push_back(cid);
	}
      } else {
	if (mdvCol != -1 && inMdv[i] == 0){
	  stop("amt or dur/rate are non-zero therefore MDV cannot = 0.");
	}
	// For Rates and non-zero amts, assume dosing event
	cevid = cmt100*100000+rateI*10000+cmt99*100+flg;
	allBolus=false;
      }
    } else {
      cevid = inEvid[i];
    }
    switch(cevid){
    case 0:
      // Observation
      nobs++;
      cevid = 0;
      if (mdvCol != -1 && inMdv[i] == 1){
	cevid = 2;
      }
      if (dvCol != -1 && ISNA(inDv[i])){
	if (amtCol==-1){
	  cevid=2;
	} else if (ISNA(inAmt[i]) || inAmt[i] == 0){
	  cevid=2;
	}
      }
      if (std::find(obsId.begin(), obsId.end(), cid) == obsId.end()){
	obsId.push_back(cid);
      }
      if (caddl > 0){
	warning("addl is ignored with observations.");
      }
      if (flg != 1){
	// warning("ss is ignored with observations.");
	flg=1;
      }
      if (ISNA(ctime)){
	id.push_back(cid);
	evid.push_back(2);
	cmtF.push_back(cmt);
	time.push_back(ctime);
	amt.push_back(NA_REAL);
	ii.push_back(0.0);
	idx.push_back(i);
	dv.push_back(cdv);
	idxO.push_back(curIdx);curIdx++;
	cevid = -1;
      } else {
	bool goodCmt = false;
	int cmpCmt;
	if ((curDvid.size()) > 1){
	  for (j=curDvid.size();j--;){
	    if (curDvid[j] > 0){
	      if (curDvid[j] == cmt || curDvid[j] == -cmt){
		goodCmt=true;
		break;
	      }
	    } else {
	      cmpCmt = state.size()+curDvid[j]+1;
	      if (cmpCmt == cmt || cmpCmt == -cmt){
		goodCmt=true;
		break;
	      }
	    }
	  }
	  if (combineDvidB && dvidCol != -1 && !IntegerVector::is_na(inDvid[i]) &&
	      inDvid[i]>0){
	    if (goodCmt && cmt != inDvid[i] && cmt != 1 && cmt != 0){
	      stop("'cmt' and 'dvid' specify different compartments; Please correct.");
	    }
	    cmt = inDvid[i];
	    goodCmt=true;
	  }
	} else {
	  goodCmt = true;
	  if (combineDvidB && dvidCol != -1 && !IntegerVector::is_na(inDvid[i]) &&
	      inDvid[i]>0){
	    cmt = inDvid[i];
	  }
	}
	if (!goodCmt){
	  IntegerVector dvidDF(curDvid.size());
	  for (i = dvidDF.size(); i--;){
	    dvidDF[i] = i+1;
	  }
	  List dvidTrans = List::create(_["dvid"]=dvidDF, _["modeled cmt"]=curDvid);
	  dvidTrans.attr("class") = "data.frame";
	  dvidTrans.attr("row.names") = IntegerVector::create(NA_INTEGER, -dvidDF.size());
	  Rprintf("DVID/CMT translation:\n");
	  print(dvidTrans);
	  if (dvidCol != -1){
	    Rprintf("DVID: %d\t", inDvid[i]);
	  }
	  Rprintf("CMT: %d\n", cmt);
	  stop("'dvid'->'cmt' or 'cmt' on observation record on a undefined compartment (use `cmt()` `dvid()`).");
	}
	id.push_back(cid);
	evid.push_back(cevid);
	cmtF.push_back(cmt);
	time.push_back(ctime);
	if (ctime == 0){
	  if (std::find(zeroId.begin(), zeroId.end(), cid) == zeroId.end()){
	    zeroId.push_back(cid);
	  }
	}
	amt.push_back(NA_REAL);
	ii.push_back(0.0);
	idx.push_back(i);
	dv.push_back(cdv);
	idxO.push_back(curIdx);curIdx++;
	cevid = -1;
      }
      break;
    case 1:
      if (mdvCol != -1 && inMdv[i] == 0){
	stop("MDV cannot be 0 when EVID=1");
      }
      cevid = cmt100*100000+rateI*10000+cmt99*100+flg;
      if (rateI == 0) allInf=false;
      else allBolus=false;
      break;
    case 2:
      cevid = 2;
      if (flg == 30){
	cevid = cmt100*100000+rateI*10000+cmt99*100+flg;
	if (rateI == 0) allInf=false;
	else allBolus=false;
      } else {
	cevid = cmt100*100000+rateI*10000+cmt99*100+1;
	if (rateI == 0) allInf=false;
	else allBolus=false;
	if (caddl > 0){
	  warning("addl is ignored with EVID=2.");
	}
	if (flg != 1){
	  warning("ss is ignored with EVID=2.");
	}	
	id.push_back(cid);
	evid.push_back(2);
	cmtF.push_back(cmt);
	time.push_back(ctime);
	if (ctime == 0){
	  if (std::find(zeroId.begin(), zeroId.end(), cid) == zeroId.end()){
	    zeroId.push_back(cid);
	  }
	}
	amt.push_back(NA_REAL);
	ii.push_back(0.0);
	idx.push_back(i);
	dv.push_back(NA_REAL);
	idxO.push_back(curIdx);curIdx++;
	ndose++;
	// + cmt needs to turn on cmts.
	// This gives a zero dose to cmt
	if (cmtCol != -1 && cmt > 0 && cmt <= baseSize){
	  // Turn on state with dose
	  id.push_back(cid);
	  evid.push_back(cevid);
	  cmtF.push_back(cmt);
	  time.push_back(ctime);
	  amt.push_back(0.0);
	  ii.push_back(0.0);
	  idx.push_back(i);
	  dv.push_back(NA_REAL);
	  idxO.push_back(curIdx);curIdx++;
	  ndose++;
	}
	
	cevid = -1;
      }
      break;
    case 3:
      cevid = 3;
      if (caddl > 0){
	warning("addl is ignored with EVID=3.");
      }
      if (flg != 1){
	warning("ss is ignored with EVID=3.");
      }
      id.push_back(cid);
      evid.push_back(3);
      cmtF.push_back(cmt);
      time.push_back(ctime);
      if (ctime == 0){
	if (std::find(zeroId.begin(), zeroId.end(), cid) == zeroId.end()){
	  zeroId.push_back(cid);
	}
      }
      amt.push_back(NA_REAL);
      ii.push_back(0.0);
      idx.push_back(i);
      dv.push_back(NA_REAL);
      idxO.push_back(curIdx);curIdx++;
      ndose++;
      cevid = -1;
      break;
    case 4:
      if (mdvCol != -1 && inMdv[i] == 0){
	stop("MDV cannot be 0 when EVID=4");
      }
      id.push_back(cid);
      evid.push_back(3);
      cmtF.push_back(cmt);
      time.push_back(ctime);
      if (ctime == 0){
	if (std::find(zeroId.begin(), zeroId.end(), cid) == zeroId.end()){
	  zeroId.push_back(cid);
	}
      }
      amt.push_back(NA_REAL);
      ii.push_back(0.0);
      idx.push_back(-1);
      dv.push_back(NA_REAL);
      idxO.push_back(curIdx);curIdx++;
      ndose++;
      // Now use the transformed compartment
      cevid = cmt100*100000+rateI*10000+cmt99*100+flg;
      if (rateI == 0) allInf=false;
      else allBolus=false;
      break;
    default:
      if (cevid > 10 && cevid < 99){
	continue;
      }
      if (rateI != 0){
	warning("'rate' or 'dur' is ignored with classic RxODE EVIDs");
	rateI = 0;
      }
      if (flg!=1){ // ss=1 is the same as ss=0 for NONMEM
	warning("'ss' is ignored with classic RxODE EVIDs.");
	flg=1;
      }
    }
    if (cevid != -1){
      if (rateI == 9){
	nevid = cmt100*100000+70001+cmt99*100;
      } else if (rateI == 8) {
	nevid = cmt100*100000+60001+cmt99*100;
      }
      id.push_back(cid);
      evid.push_back(cevid);
      cmtF.push_back(cmt);
      time.push_back(ctime);
      if (ctime == 0){
	if (std::find(zeroId.begin(), zeroId.end(), cid) == zeroId.end()){
	  zeroId.push_back(cid);
	}
      }
      if (flg >= 10){
	ii.push_back(cii);
	if (caddl > 0){
	  stop("ss with addl not supported yet.");
	}
      } else {
	ii.push_back(0.0);
      }
      idx.push_back(i);
      dv.push_back(NA_REAL);
      idxO.push_back(curIdx);curIdx++;
      ndose++;
      if (rateI > 2){
	amt.push_back(camt);
	// turn off
	id.push_back(cid);
	evid.push_back(nevid);
	cmtF.push_back(cmt);
	time.push_back(ctime);
	amt.push_back(camt);
	ii.push_back(0.0);
	idx.push_back(-1);
	dv.push_back(NA_REAL);
	idxO.push_back(curIdx);curIdx++;
	ndose++;
      } else if (rateI == 1 || rateI == 2){
	// In this case amt needs to be changed.
	dur = camt/rate;
	amt.push_back(rate); // turn on
	// turn off
	id.push_back(cid);
	evid.push_back(cevid);
	cmtF.push_back(cmt);
	time.push_back(ctime+dur);
	amt.push_back(-rate);
	ii.push_back(0.0);
	idx.push_back(-1);
	dv.push_back(NA_REAL);
	idxO.push_back(curIdx);curIdx++;
	ndose++;
      } else {
	amt.push_back(camt);
      }
      if (cii > 0 && caddl > 0 && flg < 10){
	for (j=caddl;j--;){
	  ctime+=cii;
	  id.push_back(cid);
	  evid.push_back(cevid);
	  cmtF.push_back(cmt);
	  time.push_back(ctime);
	  ii.push_back(0.0);
	  idx.push_back(-1);
	  dv.push_back(NA_REAL);
	  idxO.push_back(curIdx);curIdx++;
	  ndose++;
	
	  if (rateI > 2){
	    amt.push_back(camt);
	    // turn off
	    id.push_back(cid);
	    evid.push_back(nevid);
	    cmtF.push_back(cmt);
	    time.push_back(ctime);
	    amt.push_back(camt);
	    ii.push_back(0.0);
	    idx.push_back(-1);
	    dv.push_back(NA_REAL);
	    idxO.push_back(curIdx);curIdx++;
	    ndose++;
	  } else if (rateI == 1 || rateI == 2){
	    amt.push_back(rate);
	    // turn off
	    id.push_back(cid);
	    evid.push_back(cevid);
	    cmtF.push_back(cmt);
	    time.push_back(ctime+dur);
	    amt.push_back(-rate);
	    ii.push_back(0.0);
	    idx.push_back(-1);
	    dv.push_back(NA_REAL);
	    idxO.push_back(curIdx);curIdx++;
	    ndose++;
	  } else {
	    amt.push_back(camt);
	  }
	}
      }
    }
  }
  bool redoId=false;
  if (!keepDosingOnly){
    if (obsId.size() != allId.size()){
      std::string idWarn = "IDs without observations dropped:";
      for (j = allId.size(); j--;){
	if (std::find(obsId.begin(), obsId.end(), allId[j]) == obsId.end()){
	  idWarn = idWarn + " " +as<std::string>(idLvl[allId[j]-1]);
	  doseId.push_back(allId[j]);
	}
      }
      warning(idWarn.c_str());
      redoId=true;
    }
  }
  if (zeroId.size() != allId.size()){
    std::string idWarn = "IDs without zero-time start at the first observed time:";
    for (j = allId.size(); j--;){
      if (std::find(zeroId.begin(), zeroId.end(), allId[j]) == zeroId.end()){
	bool skipIt=false;
	if (!keepDosingOnly){
	  // Excluded from list
	  if (std::find(obsId.begin(), obsId.end(), allId[j]) == obsId.end()){
	    skipIt=true;
	  }
	}
	if (!skipIt){
	  if (!_ini0){
	    idWarn = idWarn + " " +as<std::string>(idLvl[allId[j]-1]);
	  } else {
	    id.push_back(allId[j]);
	    evid.push_back(9);
	    cmtF.push_back(0);
	    time.push_back(0.0);
	    amt.push_back(NA_REAL);
	    ii.push_back(0.0);
	    dv.push_back(NA_REAL);
	    idx.push_back(-1);
	    idxO.push_back(curIdx);curIdx++;	  
	  }
	}
      }
    }
    if (!_ini0) warning(idWarn.c_str());
  }
  std::sort(idxO.begin(),idxO.end(),
	    [id,time,evid,amt,doseId,keepDosingOnly](int a, int b){
	      // Bad IDs are pushed to the end to be popped off.
	      if (!keepDosingOnly){
		if (doseId.size() > 0 && !(std::find(doseId.begin(), doseId.end(), id[a]) == doseId.end())){
		  return false;
		}
		if (doseId.size() > 0 && !(std::find(doseId.begin(), doseId.end(), id[b]) == doseId.end())){
		  return true;
		}
	      }
	      if (id[a] == id[b]){
		if (time[a] == time[b]){
		  if (evid[a] == evid[b]){
		    return a < b;
		  }
		  // Reset should be before all other events.
		  if (evid[a] == 3){
		    return true;
		  }
		  if (evid[b] == 3){
		    return false;
		  }
		  // Zero amts turn on and off compartments and should be first.
		  if (amt[a] == 0){
		    return true;
		  }
		  if (amt[b] == 0){
		    return false;
		  }
		  return evid[a] > evid[b];
		}
		return time[a] < time[b];
	      }
	      return id[a] < id[b];
	    });
  if (!keepDosingOnly && doseId.size() > 0){
    while (idxO.size() > 0 && std::find(doseId.begin(), doseId.end(), id[idxO.back()]) != doseId.end()){
      idxO.pop_back();
    }
  }
  if (idxO.size()==0) stop("Empty data.");
  int lastId = id[idxO.back()]+42;
  int rmAmt = 0;
  // Remove trailing doses
  if (!keepDosingOnly){
    for (j =idxO.size(); j--; ){
      if (lastId != id[idxO[j]]){
	// New id
	lastId = id[idxO[j]];
	while (isDose(evid[idxO[j]]) && j--){
	  idxO[j+1] = -1;
	  rmAmt++;
	}
	if (j <= 0) break;
      }
    }
    // Remove trailing dose from idO
    while(idxO.size() > 0 && idxO.back() == -1){
      idxO.pop_back();
      rmAmt--;
    }
  }
  if (idxO.size()-rmAmt <= 0) stop("Empty data.");
  nid = obsId.size();
  NumericVector fPars = NumericVector(pars.size()*nid, NA_REAL);
  // sorted create the vectors/list
  if (addCmt && !hasCmt){
    baseSize = 7;
  } else {
    baseSize = 6;
  }
  List lst = List(baseSize+covCol.size());
  std::vector<bool> sub0(baseSize+covCol.size(), true);
  CharacterVector nme(baseSize+covCol.size());
  
  lst[0] = IntegerVector(idxO.size()-rmAmt);
  nme[0] = "ID";
  
  lst[1] = NumericVector(idxO.size()-rmAmt);
  nme[1] = "TIME";
  
  lst[2] = IntegerVector(idxO.size()-rmAmt);
  nme[2] = "EVID";
  
  lst[3] = NumericVector(idxO.size()-rmAmt);
  nme[3] = "AMT";
  
  lst[4] = NumericVector(idxO.size()-rmAmt);
  nme[4] = "II";
  
  lst[5] = NumericVector(idxO.size()-rmAmt);
  nme[5] = "DV";

  if (baseSize == 7){
    lst[6] = IntegerVector(idxO.size()-rmAmt);
    nme[6] = "CMT";
  }
  

  List lst1(1+covCol.size());
  CharacterVector nme1(1+covCol.size());
  std::vector<bool> sub1(1+covCol.size(), true);

  lst1[0] = IntegerVector(nid);
  nme1[0] = "ID";
  
  for (j = 0; j < (int)(covCol.size()); j++){
    if (as<std::string>(lName[covCol[j]]) == "cmt"){
      lst[baseSize+j] = IntegerVector(idxO.size()-rmAmt);
      nme[baseSize+j] = pars[covParPos[j]];
      sub0[baseSize+j] = false;
      lst1[1+j] = IntegerVector(nid);
      nme1[1+j] = nme[baseSize+j];
      sub1[1+j] = true;
    } else {
      lst[baseSize+j] = NumericVector(idxO.size()-rmAmt);
      nme[baseSize+j] = pars[covParPos[j]];
      sub0[baseSize+j] = false;
      lst1[1+j] = NumericVector(nid);
      nme1[1+j] = nme[baseSize+j];
      sub1[1+j] = true;
    }
  }
  List keepL = List(keepCol.size());
  CharacterVector keepN(keepCol.size());
  for (j = 0; j < (int)(keepCol.size()); j++){
    keepL[j] = NumericVector(idxO.size()-rmAmt);
    keepN[j] = dName[keepCol[j]];
  }

  IntegerVector ivTmp;
  // Since we removed the -1 in idx, you can get the last id here.
  lastId = id[idxO.back()]+1;
  bool addId = false, added=false;
  int idx1=nid, nTv=0;
  std::vector<int> covParPosTV;
  bool cmtFadd = false;
  int jj = idxO.size()-rmAmt;
  for (i =idxO.size(); i--;){
    if (idxO[i] != -1){
      jj--;
      ivTmp = as<IntegerVector>(lst[0]);
      ivTmp[jj] = id[idxO[i]];
      if (lastId != id[idxO[i]]){
	addId=true;
	idx1--;
	// Add ID
	ivTmp = as<IntegerVector>(lst1[0]);
	ivTmp[idx1] = id[idxO[i]];
	lastId=id[idxO[i]];
      }
      // retId[i]=id[idxO[i]];
      nvTmp = as<NumericVector>(lst[1]);
      // retTime[i]=time[idxO[i]];
      nvTmp[jj] = time[idxO[i]];
      ivTmp = as<IntegerVector>(lst[2]);
      ivTmp[jj] = evid[idxO[i]];
      // retEvid[i]=evid[idxO[i]];
      nvTmp = as<NumericVector>(lst[3]);
      // retAmt[i]=amt[idxO[i]];
      nvTmp[jj] = amt[idxO[i]];
      nvTmp = as<NumericVector>(lst[4]);
      nvTmp[jj]=ii[idxO[i]];
      nvTmp = as<NumericVector>(lst[5]);
      nvTmp[jj]=dv[idxO[i]];
      if (baseSize == 7){
	ivTmp = as<IntegerVector>(lst[6]);
	ivTmp[jj] = cmtF[idxO[i]];
      }
      // Now add the other items.
      added=false;
      for (j = 0; j < (int)(keepCol.size()); j++){
	nvTmp = as<NumericVector>(keepL[j]);
	if (idx[idxO[i]] == -1){
	  // Implement LOCF
	  if (addId){
	    nvTmp[jj] = NA_REAL;
	  } else {
	    // LOCF
	    nvTmp[jj] = nvTmp[jj-1];
	  }
	} else {
	  // These keepers are added.
	  nvTmp2   = as<NumericVector>(inData[keepCol[j]]);
	  nvTmp[jj] = nvTmp2[idx[idxO[i]]];
	}
      }
      for (j = 0; j < (int)(covCol.size()); j++){
	if (as<std::string>(lName[covCol[j]]) == "cmt"){
	  ivTmp = as<IntegerVector>(lst[baseSize+j]);
	  ivTmp[jj] = cmtF[idxO[i]];
	  if (!cmtFadd){
	    sub0[baseSize+j] = true;
	    sub1[1+j] = false;
	    covParPosTV.push_back(covParPos[j]);
	    cmtFadd=true;
	    nTv++;
	  }
	} else {
	  nvTmp = as<NumericVector>(lst[baseSize+j]);
	  if (idx[idxO[i]] == -1){
	    // These should be ignored for interpolation.
	    nvTmp[jj] = NA_REAL;
	  } else {
	    // These covariates are added.
	    nvTmp2   = as<NumericVector>(inData[covCol[j]]);
	    nvTmp[jj] = nvTmp2[idx[idxO[i]]];
	    if (addId){
	      nvTmp = as<NumericVector>(lst1[1+j]);
	      nvTmp[idx1] = nvTmp2[idx[idxO[i]]];
	      fPars[idx1*pars.size()+covParPos[j]] = nvTmp[idx1];
	      added = true;
	    } else if (sub1[1+j]) {
	      nvTmp = as<NumericVector>(lst1[1+j]);
	      if (nvTmp[idx1] != nvTmp2[idx[idxO[i]]]){
		sub0[baseSize+j] = true;
		sub1[1+j] = false;
		fPars[idx1*pars.size()+covParPos[j]] = NA_REAL;
		if (std::find(covParPosTV.begin(), covParPosTV.end(), covParPos[j]) == covParPosTV.end()){
		  covParPosTV.push_back(covParPos[j]);
		}
		nTv++;
	      }
	    }
	  }
	}
      }
      if (added && addId){
	addId=false;
	added=false;
      }
    }
  }
  if (!dropUnits && addTimeUnits){
    NumericVector tmpN = as<NumericVector>(lst[1]);
    tmpN.attr("class") = "units";
    tmpN.attr("units") = timeUnits;
  }
  if (!dropUnits && addAmtUnits){
    NumericVector tmpN = as<NumericVector>(lst[3]);
    tmpN.attr("class") = "units";
    tmpN.attr("units") = amtUnits;
  }
  // Now subset based on time-varying covariates
  List lstF;
  CharacterVector nmeF;
  if (allTimeVar){
    lstF = List(baseSize+covCol.size());
    nmeF = CharacterVector(baseSize+covCol.size());
  } else {
    lstF = List(baseSize+nTv);
    nmeF = CharacterVector(baseSize+nTv);
  }
  j=0;
  for (i = 0; i < lst.size();i++){
    if (allTimeVar || sub0[i]){
      lstF[j]=lst[i];
      nmeF[j]=nme[i];
      j++;
    }
  }
  j=0;
  List lst1F(1+covCol.size()-nTv);
  CharacterVector nme1F(1+covCol.size()-nTv);
  for (i = 0; i < lst1.size();i++){
    if (sub1[i]){
      lst1F[j]=lst1[i];
      nme1F[j]=nme1[i];
      j++;
    }
  }

  CharacterVector cls = CharacterVector::create("rxEtTran","data.frame");
  
  IntegerVector tmp = lst1F[0];
  if (redoId){
    Function convId = rx[".convertId"];
    tmp = convId(tmp);//as<IntegerVector>();
    idLvl = tmp.attr("levels");
    tmp.attr("class")  = R_NilValue;
    tmp.attr("levels") = R_NilValue;
    lst1F[0] = tmp;
  }
  
  lst1F.attr("names") = nme1F;
  lst1F.attr("class") = CharacterVector::create("data.frame");
  lst1F.attr("row.names") = IntegerVector::create(NA_INTEGER, -nid);
  List e;
  e["ndose"] = ndose;
  e["nobs"]  = nobs;
  e["nid"]   = nid;
  e["cov1"] = lst1F;
  e["covParPos"]  = wrap(covParPos);
  e["covParPosTV"] = wrap(covParPosTV); // Time-varying pos
  if (allTimeVar){
    e["sub0"] = wrap(sub0);
    e["baseSize"] = baseSize;
    e["nTv"] = nTv;
    e["lst"] = lst;
    e["nme"] = nme;
  }
  std::vector<int> covParPos0;
  for (j = covParPos.size();j--;){
    if (std::find(covParPosTV.begin(), covParPosTV.end(), covParPos[j]) == covParPosTV.end()){
      covParPos0.push_back(covParPos[j]);
    }
  }
  e["covParPos0"] = wrap(covParPos0);
  e["covUnits"] = covUnits;
  fPars.attr("dim")= IntegerVector::create(pars.size(), nid);
  fPars.attr("dimnames") = List::create(pars, R_NilValue);
  e["pars"] = fPars;
  e["allBolus"] = allBolus;
  e["allInf"] = allInf;
  e["mxCmt"] = mxCmt;
  e["lib.name"] = trans["lib.name"];
  e["addCmt"] = addCmt;
  e["cmtInfo"] = cmtInfo;
  e["idLvl"] = idLvl;
  e["allTimeVar"] = allTimeVar;
  e["keepDosingOnly"] = true;
  keepL.attr("names") = keepN;
  keepL.attr("class") = CharacterVector::create("data.frame");
  keepL.attr("row.names") = IntegerVector::create(NA_INTEGER,-idxO.size()+rmAmt);
  setFkeep(keepL);
  e.attr("class") = "rxHidden";
  cls.attr(".RxODE.lst") = e;
  tmp = lstF[0];
  if (redoId){
    Function convId = rx[".convertId"];
    tmp = convId(tmp);//as<IntegerVector>();
    idLvl = tmp.attr("levels");
    tmp.attr("class")  = R_NilValue;
    tmp.attr("levels") = R_NilValue;
    lstF[0]=tmp;
  }
  if (!dropUnits){
    tmp.attr("class") = "factor";
    tmp.attr("levels") = idLvl;
  }
  lstF.attr("names") = nmeF;
  lstF.attr("class") = cls;
  lstF.attr("row.names") = IntegerVector::create(NA_INTEGER,-idxO.size()+rmAmt);
  return lstF;
}
