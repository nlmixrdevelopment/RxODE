#include <Rcpp.h>
#include <algorithm>
#define rxModelVars(a) rxModelVars_(a)
using namespace Rcpp;

List rxModelVars_(const RObject &obj);
bool rxIs(const RObject &obj, std::string cls);
Environment RxODEenv();

IntegerVector toCmt(RObject inCmt, CharacterVector state){
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
	  k++;
	  lvlI[i] = state.size() + k;
	}
      }
      IntegerVector cmtIn = IntegerVector(inCmt);
      IntegerVector ret(cmtIn.size());
      for (j=ret.size(); j--;){
	ret[j] = lvlI[cmtIn[j]-1];
      }
      return ret;
    } else {
      return as<IntegerVector>(inCmt);
    }
  } else if (rxIs(inCmt, "character")) {
    CharacterVector iCmt = as<CharacterVector>(inCmt);
    std::vector<int> newCmt;
    newCmt.reserve(iCmt.size());
    std::string strCmt, negSub;
    bool foundState=false;
    bool isNeg = false;
    List extraCmt;
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
      if (strCmt == "(default)" || strCmt == "(obs)"){
	newCmt.push_back(1);
      } else {
	foundState=false;
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
		warning("Negative compartments on non-ode cmt (%s) do not make sense.", strCmt.c_str());
		newCmt.push_back(-state.size()-j-1);
	      } else {
		newCmt.push_back(state.size()+j+1);
	      }
	      break;
	    }
	  }
	  if (!foundState){
	    extraCmt[k++] = CharacterVector::create(strCmt);
	    newCmt.push_back(state.size()+k);
	  }
	}
      }
    }
    return wrap(newCmt);
  }
  stop("Should not reach here.");
  return IntegerVector::create(0);
}

//[[Rcpp::export]]
List etTrans(List inData, const RObject &obj, bool addCmt=false){
  // Translates events + model into translated events
  CharacterVector lName = clone(as<CharacterVector>(inData.attr("names")));
  int i, idCol = -1, evidCol=-1, timeCol=-1, amtCol=-1, cmtCol=-1,
    dvCol=-1, ssCol=-1, rateCol=-1, addlCol=-1, iiCol=-1, durCol=-1, j;
  std::string tmpS;
  List mv = rxModelVars_(obj);
  CharacterVector pars = as<CharacterVector>(mv["params"]);
  std::vector<int> covCol;
  std::vector<int> covParPos;
  std::string tmpS0;
  for (i = lName.size(); i--;){
    tmpS0= as<std::string>(lName[i]);
    tmpS = as<std::string>(lName[i]);
    std::transform(tmpS.begin(), tmpS.end(), tmpS.begin(), ::tolower);
    lName[i] = tmpS;
    if (tmpS == "id") idCol=i;
    else if (tmpS == "evid") evidCol=i;
    else if (tmpS == "time") timeCol=i;
    else if (tmpS == "amt") amtCol=i;
    else if (tmpS == "cmt" || tmpS == "ytype") cmtCol=i;
    else if (tmpS == "dv" || tmpS == "y") dvCol=i;
    else if (tmpS == "ss")   ssCol=i;
    else if (tmpS == "rate") rateCol=i;
    else if (tmpS == "dur") durCol=i;
    else if (tmpS == "addl") addlCol=i;
    else if (tmpS == "ii")   iiCol=i;
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
      if (rxIs(nvTmp, "units")){
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
  // EVID = 10-99; mtime events (from ODE system)
  // When EVID > 100
  // EVID: ## # ## ##
  //       c2 I c1 xx
  // c2 = Compartment numbers over 100
  //  I = Infusion Flag
  //      0 = no Infusion
  //      1 = Infusion, AMT=rate (mg/hr for instance)
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
  
  CharacterVector state = as<CharacterVector>(mv["state"]);
  std::vector<int> id;
  std::vector<int> allId;
  std::vector<int> evid;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> ii;
  std::vector<int> idx;
  std::vector<int> cmtF;
  std::vector<double> dv;
  std::vector<int> idxO;
  if (timeCol== -1){
    stop("time is required in dataset.");
  }
  NumericVector inTime = as<NumericVector>(inData[timeCol]);
  // save units information
  bool addTimeUnits = false;
  RObject timeUnits;
  if (rxIs(inTime, "units")){
    addTimeUnits=true;
    timeUnits=inTime.attr("units");
  }
  IntegerVector inCmt;
  if (cmtCol != -1){
    inCmt = toCmt(inData[cmtCol], mv["state"]);//as<IntegerVector>();
  }
  int tmpCmt;
  IntegerVector inId;
  if (idCol != -1){
    inId = as<IntegerVector>(inData[idCol]);
  }
  IntegerVector inSs;
  if (ssCol != -1){
    inSs = as<IntegerVector>(inData[ssCol]);
  }
  IntegerVector inEvid;
  if (evidCol != -1){
    inEvid = as<IntegerVector>(inData[evidCol]);
  }
  NumericVector inRate;
  if (rateCol != -1){
    inRate = as<NumericVector>(inData[rateCol]);
  }

  NumericVector inDur;
  if (durCol != -1){
    inDur = as<NumericVector>(inData[durCol]);
  }
  
  bool addAmtUnits = false;
  RObject amtUnits;
  NumericVector inAmt;
  if (amtCol != -1){
    inAmt = as<NumericVector>(inData[amtCol]);
    if (rxIs(inAmt, "units")){
      addAmtUnits=true;
      amtUnits=inAmt.attr("units");
    }
  }
  NumericVector inIi;
  if (iiCol != -1){
    inIi = as<NumericVector>(inData[iiCol]);
  }
  IntegerVector inAddl;
  if (addlCol != -1){
    inAddl = as<IntegerVector>(inData[addlCol]);
  }
  NumericVector inDv;
  if (dvCol != -1){
    inDv = as<NumericVector>(inData[dvCol]);
  }
  int ss = 0;
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
    if (iiCol == -1) cii = 0;
    else cii = inIi[i];
    if (std::find(allId.begin(), allId.end(), cid) == allId.end()){
      allId.push_back(cid);
      // New ID
      // Add mtime records
      for (j = nMtime; j--;){
	id.push_back(cid);
	evid.push_back(j+10);
	cmtF.push_back(cmt);
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
    ss=1;
    if (ssCol == -1) ss=1;
    else if (inSs[i] == 0) ss=1;
    else if (IntegerVector::is_na(inSs[i])) ss=1;
    else if (inSs[i] == 1 && cii > 0) ss=10;
    else if (inSs[i] == 2 && cii > 0) ss=20;
    tmpCmt = inCmt[i];
    if (inCmt[i] < 0){
      if (ss != 1) stop("Steady state records cannot be on negative compartments.");
      ss = 30;
      tmpCmt = -tmpCmt;
    }
    
    // CMT flag
    if (cmtCol == -1) cmt = 1;
    else cmt = tmpCmt;
    if (IntegerVector::is_na(tmpCmt)) cmt=1;
    if (cmt <= 99){
      cmt100=0;
      cmt99=cmt;
    } else {
      cmt100=cmt/100; 
      cmt99=cmt-cmt100*100; 
    }

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
	rateI = 1;
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
	rateI = 2;
	rate = camt/inDur[i];
      }
    } else {
      stop("'rate' and/or 'dur' are not specified correctly.");
    }
    
    if (addlCol == -1) caddl=0;
    else caddl = inAddl[i];

    // EVID flag
    if (evidCol == -1){
      // Missing EVID
      if (rateI != 0 || !ISNA(camt)){
	// For Rates and non-zero amts, assume dosing event
	cevid = cmt100*100000+rateI*10000+cmt99*100+ss;
      } else {
	cevid = 0;
      }
    } else {
      switch(inEvid[i]){
      case 0:
	nobs++;
	cevid = 0;
	if (caddl > 0){
	  warning("addl is ignored with observations.");
	}
	if (ss != 1){
	  // warning("ss is ignored with observations.");
	  ss=1;
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
	  id.push_back(cid);
	  evid.push_back(0);
	  cmtF.push_back(cmt);
	  time.push_back(ctime);
	  amt.push_back(NA_REAL);
	  ii.push_back(0.0);
	  idx.push_back(i);
	  dv.push_back(cdv);
	  idxO.push_back(curIdx);curIdx++;
	  cevid = -1;
	}
	break;
      case 1:
	cevid = cmt100*100000+rateI*10000+cmt99*100+ss;
	break;
      case 2:
	cevid = 2;
	if (ss == 30){
	  cevid = cmt100*100000+rateI*10000+cmt99*100+ss;
	} else {
	  if (caddl > 0){
	    warning("addl is ignored with EVID=2.");
	  }
	  if (ss != 1){
	    warning("ss is ignored with EVID=2.");
	  }	
	  id.push_back(cid);
	  evid.push_back(2);
	  cmtF.push_back(cmt);
	  time.push_back(ctime);
	  amt.push_back(NA_REAL);
	  ii.push_back(0.0);
	  idx.push_back(i);
	  dv.push_back(NA_REAL);
	  idxO.push_back(curIdx);curIdx++;
	  ndose++;
	  cevid = -1;
	}
	break;
      case 3:
	cevid = 3;
	if (caddl > 0){
	  warning("addl is ignored with EVID=3.");
	}
	if (ss != 1){
	  warning("ss is ignored with EVID=3.");
	}
	id.push_back(cid);
	evid.push_back(3);
	cmtF.push_back(cmt);
	time.push_back(ctime);
	amt.push_back(NA_REAL);
	ii.push_back(0.0);
	idx.push_back(i);
	dv.push_back(NA_REAL);
	idxO.push_back(curIdx);curIdx++;
	ndose++;
	cevid = -1;
	break;
      case 4:
	id.push_back(cid);
	evid.push_back(3);
	cmtF.push_back(cmt);
	time.push_back(ctime);
	amt.push_back(NA_REAL);
	ii.push_back(0.0);
	idx.push_back(-1);
	dv.push_back(NA_REAL);
	idxO.push_back(curIdx);curIdx++;
	ndose++;
	// Now use the transformed compartment
	cevid = cmt100*100000+rateI*10000+cmt99*100+ss;
	break;
      default:
	cevid = inEvid[i];
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
      if (ss >= 10){
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
      if (cii > 0 && caddl > 0 && ss < 10){
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
  NumericVector fPars = NumericVector(pars.size()*nid, NA_REAL);
  std::sort(idxO.begin(),idxO.end(),
	    [id,time,evid](int a, int b){
	      if (id[a] == id[b]){
		if (time[a] == time[b]){
		  if (evid[a] == evid[b]){
		    return a < b;
	  }
		  return evid[a] > evid[b];
		}
		return time[a] < time[b];
	      }
	      return id[a] < id[b];
	    });
  // sorted create the vectors/list
  int baseSize;
  if (addCmt && !hasCmt){
    baseSize = 7;
  } else {
    baseSize = 6;
  }
  List lst = List(baseSize+covCol.size());
  std::vector<bool> sub0(baseSize+covCol.size(), true);
  CharacterVector nme(baseSize+covCol.size());
  
  lst[0] = IntegerVector(idxO.size());
  nme[0] = "id";
  
  lst[1] = NumericVector(idxO.size());
  nme[1] = "time";
  
  lst[2] = IntegerVector(idxO.size());
  nme[2] = "evid";
  
  lst[3] = NumericVector(idxO.size());
  nme[3] = "amt";
  
  lst[4] = NumericVector(idxO.size());
  nme[4] = "ii";
  
  lst[5] = NumericVector(idxO.size());
  nme[5] = "dv";

  if (baseSize == 7){
    lst[6] = IntegerVector(idxO.size());
    nme[6] = "cmt";
  }
  

  List lst1(1+covCol.size());
  CharacterVector nme1(1+covCol.size());
  std::vector<bool> sub1(1+covCol.size(), true);

  lst1[0] = IntegerVector(nid);
  nme1[0] = "id";
  
  for (j = 0; j < (int)(covCol.size()); j++){
    if (as<std::string>(lName[covCol[j]]) == "cmt"){
      lst[baseSize+j] = IntegerVector(idxO.size());
      nme[baseSize+j] = pars[covParPos[j]];
      sub0[baseSize+j] = false;
      lst1[1+j] = IntegerVector(nid);
      nme1[1+j] = nme[baseSize+j];
      sub1[1+j] = true;
    } else {
      lst[baseSize+j] = NumericVector(idxO.size());
      nme[baseSize+j] = pars[covParPos[j]];
      sub0[baseSize+j] = false;
      lst1[1+j] = NumericVector(nid);
      nme1[1+j] = nme[baseSize+j];
      sub1[1+j] = true;
    }
  }

  IntegerVector ivTmp;
  int lastId = id[idxO[idxO.size()-1]]+1;
  bool addId = false, added=false;
  int idx1=nid, nTv=0;
  std::vector<int> covParPosTV;
  bool cmtFadd = false;
  for (i =idxO.size(); i--;){
    ivTmp = as<IntegerVector>(lst[0]);
    ivTmp[i] = id[idxO[i]];
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
    nvTmp[i] = time[idxO[i]];
    ivTmp = as<IntegerVector>(lst[2]);
    ivTmp[i] = evid[idxO[i]];
    // retEvid[i]=evid[idxO[i]];
    nvTmp = as<NumericVector>(lst[3]);
    // retAmt[i]=amt[idxO[i]];
    nvTmp[i] = amt[idxO[i]];
    nvTmp = as<NumericVector>(lst[4]);
    nvTmp[i]=ii[idxO[i]];
    nvTmp = as<NumericVector>(lst[5]);
    nvTmp[i]=dv[idxO[i]];
    if (baseSize == 7){
      ivTmp = as<IntegerVector>(lst[6]);
      ivTmp[i] = cmtF[idxO[i]];
    }
    // Now add the other items.
    added=false;
    for (j = 0; j < (int)(covCol.size()); j++){
      if (as<std::string>(lName[covCol[j]]) == "cmt"){
	ivTmp = as<IntegerVector>(lst[baseSize+j]);
	ivTmp[i] = cmtF[idxO[i]];
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
	  nvTmp[i] = NA_REAL;
	} else {
	  // These covariates are added.
	  nvTmp2   = as<NumericVector>(inData[covCol[j]]);
	  nvTmp[i] = nvTmp2[idx[idxO[i]]];
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
  if (addTimeUnits){
    NumericVector tmpN = as<NumericVector>(lst[1]);
    tmpN.attr("class") = "units";
    tmpN.attr("units") = timeUnits;
  }
  if (addAmtUnits){
    NumericVector tmpN = as<NumericVector>(lst[3]);
    tmpN.attr("class") = "units";
    tmpN.attr("units") = amtUnits;
  }
  // Now subset based on time-varying covariates
  List lstF(baseSize+nTv);
  CharacterVector nmeF(baseSize+nTv);
  j=0;
  for (i = 0; i < lst.size();i++){
    if (sub0[i]){
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
  e.attr("class") = "rxHidden";
  cls.attr(".RxODE.lst") = e;
  lstF.attr("names") = nmeF;
  lstF.attr("class") = cls;
  lstF.attr("row.names") = IntegerVector::create(NA_INTEGER,-idxO.size());
  return lstF;
}
