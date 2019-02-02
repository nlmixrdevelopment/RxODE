#include <Rcpp.h>
#include <algorithm>
#define rxModelVars(a) rxModelVars_(a)
using namespace Rcpp;

List rxModelVars_(const RObject &obj);
Environment RxODEenv();

//[[Rcpp::export]]
List evTrans(List inData, const RObject &obj){
  // Translates events + model into translated events
  CharacterVector lName = as<CharacterVector>(inData.attr("names"));
  int i, idCol = -1, evidCol=-1, timeCol=-1, amtCol=-1, cmtCol=-1,
    dvCol=-1, ssCol=-1, rateCol=-1, addlCol=-1, iiCol=-1;
  std::string tmpS;
  for (i = lName.size(); i--;){
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
    else if (tmpS == "addl") addlCol=i;
    else if (tmpS == "ii")   iiCol=i;
  }
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
  // Steady state events need a II data item > 0
  List mv = rxModelVars_(obj);
  CharacterVector state = as<CharacterVector>(mv["state"]);
  std::vector<int> id;
  std::vector<int> allId;
  std::vector<int> evid;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> ii;
  std::vector<int> idx;
  std::vector<double> dv;
  std::vector<int> idxO;
  if (timeCol== -1){
    stop("time is required in dataset.");
  }
  NumericVector inTime = as<NumericVector>(inData[timeCol]);
  IntegerVector inCmt;
  if (timeCol != -1){
    inCmt = as<IntegerVector>(inData[cmtCol]);
  }
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
  NumericVector inAmt;
  if (amtCol != -1){
    inAmt = as<NumericVector>(inData[amtCol]);
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
  int j;
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
  double dur;
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
    if (iiCol == -1) cii = 0;
    else cii = inIi[i];
    if (std::find(allId.begin(), allId.end(), cid) == allId.end()){
      allId.push_back(cid);
      // New ID
      // Add mtime records
      for (j = nMtime; j--;){
	id.push_back(cid);
	evid.push_back(j+10);
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
    if (ssCol == -1) ss=1;
    else if (inSs[i] == 0) ss=1;
    else if (inSs[i] == 1 && cii > 0) ss=10;
    else if (inSs[i] == 2 && cii > 0) ss=20;
    
    // CMT flag
    if (cmtCol == -1) cmt = 1;
    else cmt = inCmt[i];
    if (cmt > 99){
      cmt100=0;
      cmt99=cmt;
    } else {
      cmt100=cmt/100; 
      cmt99=cmt-cmt100*100; 
    }

    // Rate
    rateI = 0;
    if (rateCol == -1) rate = 0.0;
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

    // Amt
    if (amtCol == -1) camt = 0.0;
    else camt = inAmt[i];

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
	  warning("ss is ignored with observations.");
	  ss=1;
	}
	id.push_back(cid);
	evid.push_back(0);
	time.push_back(ctime);
	amt.push_back(NA_REAL);
	ii.push_back(0.0);
	idx.push_back(i);
	dv.push_back(cdv);
	idxO.push_back(curIdx);curIdx++;
	cevid = -1;
	break;
      case 1:
	cevid = cmt100*100000+rateI*10000+cmt99*100+ss;
	break;
      case 2:
	cevid = 2;
	if (caddl > 0){
	  warning("addl is ignored with EVID=2.");
	}
	if (ss != 1){
	  warning("ss is ignored with EVID=2.");
	}
	id.push_back(cid);
	evid.push_back(2);
	time.push_back(ctime);
	amt.push_back(NA_REAL);
	ii.push_back(0.0);
	idx.push_back(i);
	dv.push_back(NA_REAL);
	idxO.push_back(curIdx);curIdx++;
	ndose++;
	cevid = -1;
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
      if (rateI > 1){
	amt.push_back(camt);
      
	// turn off
	id.push_back(cid);
	evid.push_back(nevid);
	time.push_back(ctime);
	amt.push_back(camt);
	ii.push_back(0.0);
	idx.push_back(-1);
	dv.push_back(NA_REAL);
	idxO.push_back(curIdx);curIdx++;
	ndose++;
      } else if (rateI == 1){
	// In this case amt needs to be changed.
	dur = camt/rate;
	amt.push_back(rate); // turn on

	// turn off
	id.push_back(cid);
	evid.push_back(cevid);
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
	  time.push_back(ctime);
	  ii.push_back(0.0);
	  idx.push_back(-1);
	  dv.push_back(NA_REAL);
	  idxO.push_back(curIdx);curIdx++;
	  ndose++;
	
	  if (rateI > 1){
	    amt.push_back(camt);
	    // turn off
	    id.push_back(cid);
	    evid.push_back(nevid);
	    time.push_back(ctime);
	    amt.push_back(camt);
	    ii.push_back(0.0);
	    idx.push_back(-1);
	    dv.push_back(NA_REAL);
	    idxO.push_back(curIdx);curIdx++;
	    ndose++;
	  } else if (rateI == 1){
	    amt.push_back(rate);
	    // turn off
	    id.push_back(cid);
	    evid.push_back(cevid);
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
  IntegerVector retId(idxO.size());
  NumericVector retTime(idxO.size());
  IntegerVector retEvid(idxO.size());
  NumericVector retAmt(idxO.size());
  NumericVector retIi(idxO.size());
  NumericVector retDv(idxO.size());
  for (i =idxO.size(); i--;){
    retId[i]=id[idxO[i]];
    retTime[i]=time[idxO[i]];
    retEvid[i]=evid[idxO[i]];
    retAmt[i]=amt[idxO[i]];
    retIi[i]=ii[idxO[i]];
    retDv[i]=dv[idxO[i]];
  }
  List lst(6);
  CharacterVector nme(6);
  lst[0] = retId;
  nme[0] = "id";
  
  lst[1] = retTime;
  nme[1] = "time";
  
  lst[2] = retEvid;
  nme[2] = "evid";
  
  lst[3] = retAmt;
  nme[3] = "amt";
  
  lst[4] = retIi;
  nme[4] = "ii";
  
  lst[5] = retDv;
  nme[5] = "dv";

  CharacterVector cls = CharacterVector::create("rxEvTran","data.frame");
  Function newEnv("new.env", R_BaseNamespace);
  Environment e = newEnv(_["size"] = 29, _["parent"] = RxODEenv());
  e["ndose"] = ndose;
  e["nobs"] = nobs;
  e["nid"] = nid;
  cls.attr(".RxODE.env") = e;

  lst.attr("names") = nme;
  lst.attr("class") = cls;
  lst.attr("row.names") = IntegerVector::create(NA_INTEGER,-retId.size());
  return lst;
}
