#include <Rcpp.h>
using namespace Rcpp;

bool rxIs(const RObject &obj, std::string cls);
Environment RxODEenv();

RObject evCur;

Function loadNamespace2("loadNamespace", R_BaseNamespace);
Environment unitsPkg = loadNamespace2("units");
//(x, value, ..., mode = units_options("set_units_mode"))
NumericVector setUnits(NumericVector obj, std::string unit){
  Function f = as<Function>(unitsPkg["set_units"]);
  return as<NumericVector>(f(_["x"] = obj, _["value"] = unit, _["mode"] = "standard"));
}

//[[Rcpp::export]]
RObject etUpdate(RObject obj,
		 RObject arg = R_NilValue,
		 RObject value = R_NilValue){
  if (rxIs(obj,"rxEt")){
    evCur = obj;
    if (rxIs(value, "NULL")){
      CharacterVector cls = obj.attr("class");
      List e = clone(as<List>(cls.attr(".RxODE.lst")));
      if (rxIs(arg, "character")){
	CharacterVector carg = as<CharacterVector>(arg);
	std::string sarg = as<std::string>(carg[0]);
	if (sarg == "env"){
	  return as<RObject>(e);
	} else if (e.containsElementNamed(sarg.c_str())){
	  return as<RObject>(e[sarg]);
	}
      }
    } else {
      // Assign
    }
  }
  return R_NilValue;
}

List etEmpty(CharacterVector units){
  CharacterVector cls = CharacterVector::create("rxEt","data.frame");
  List e;
  e[".units"] = units;
  Function parse2("parse", R_BaseNamespace);
  Function eval2("eval", R_BaseNamespace);
  // eventTable style methods
  e["get.units"] = eval2(_["expr"]   = parse2(_["text"]="function() .units"),
			 _["envir"]  = e);
  
  e["getUnits"] = eval2(_["expr"]   = parse2(_["text"]="function() .units"),
			_["envir"]  = e);
  e["get_units"] = eval2(_["expr"]   = parse2(_["text"]="function() .units"),
			 _["envir"]  = e);

  std::string addDosing= "function(dose, nbr.doses = 1, dosing.interval = 24, dosing.to=1, rate=NULL, amount.units = NA, start.time, do.sampling=FALSE, time.units = NA, ...) {invisible(.Call(`_RxODE_et_`, as.list(match.call())[-1] ,list('last')))}";

  e["add.dosing"] = eval2(_["expr"] = parse2(_["text"] = addDosing),
			  _["envir"]  = e);
  e["add_dosing"] = eval2(_["expr"] = parse2(_["text"] = addDosing),
			  _["envir"]  = e);
  e["addDosing"] = eval2(_["expr"] = parse2(_["text"] = addDosing),
			 _["envir"]  = e);

  std::string addSampling="function(time, time.units = NA) invisible(.Call(`_RxODE_et_`, as.list(match.call())[-1], list('last')))";

  e["add.sampling"] = eval2(_["expr"] = parse2(_["text"] = addSampling),
			    _["envir"]  = e);
  e["add_sampling"] = eval2(_["expr"] = parse2(_["text"] = addSampling),
			    _["envir"]  = e);
  e["addSampling"] = eval2(_["expr"] = parse2(_["text"] = addSampling),
			   _["envir"]  = e);
  
  e["clear.sampling"] = eval2(_["expr"] = parse2(_["text"] = "function() invisible(.Call(`_RxODE_et_`, list(clearSampling=TRUE),list('last')))"),
			      _["envir"]  = e);

  e["clear_sampling"] = eval2(_["expr"] = parse2(_["text"] = "function() invisible(.Call(`_RxODE_et_`, list(clearSampling=TRUE),list('last')))"),
			      _["envir"]  = e);

  e["clearSampling"] = eval2(_["expr"] = parse2(_["text"] = "function() invisible(.Call(`_RxODE_et_`, list(clearSampling=TRUE),list('last')))"),
			     _["envir"]  = e);


  e["clear.dosing"] = eval2(_["expr"] = parse2(_["text"] = "function() invisible(.Call(`_RxODE_et_`, list(clearDosing=TRUE),list('last')))"),
			      _["envir"]  = e);

  e["clear_dosing"] = eval2(_["expr"] = parse2(_["text"] = "function() invisible(.Call(`_RxODE_et_`, list(clearDosing=TRUE),list('last')))"),
			      _["envir"]  = e);

  e["clearDosing"] = eval2(_["expr"] = parse2(_["text"] = "function() invisible(.Call(`_RxODE_et_`, list(clearDosing=TRUE),list('last')))"),
			     _["envir"]  = e);

  e["import.EventTable"] = eval2(_["expr"] = parse2(_["text"] = "function() invisible(.Call(`_RxODE_et_`, list(...),list('last')))"),
				 _["envir"]  = e);

  e["importEventTable"] = eval2(_["expr"] = parse2(_["text"] = "function() .Call(`_RxODE_et_`, list(...),list('last'))"),
				_["envir"]  = e);

  e["copy"] = eval2(_["expr"] = parse2(_["text"] = "function() .Call(`_RxODE_et_`, list(copy=TRUE),list('last'))"),
		    _["envir"]  = e);
  
  e["get.EventTable"] = eval2(_["expr"] = parse2(_["text"] = "function() .Call(`_RxODE_et_`, list(get.EventTable=TRUE),list('last'))"),
			      _["envir"]  = e);
  e["getEventTable"] = eval2(_["expr"] = parse2(_["text"] = "function() .Call(`_RxODE_et_`, list(get.EventTable=TRUE),list('last'))"),
			     _["envir"]  = e);
  e["get.obs.rec"] = eval2(_["expr"] = parse2(_["text"] = "function() .Call(`_RxODE_et_`, list(get.obs.rec=TRUE),list('last'))"),
			   _["envir"]  = e);

  e["get.nobs"] = eval2(_["expr"] = parse2(_["text"] = "function() .Call(`_RxODE_et_`, list(get.nobs=TRUE),list('last'))"),
			   _["envir"]  = e);


  e["get.dosing"] = eval2(_["expr"] = parse2(_["text"] = "function() .Call(`_RxODE_et_`, list(get.dosing=TRUE),list('last'))"),
			  _["envir"]  = e);
  e["getDosing"] = eval2(_["expr"] = parse2(_["text"] = "function() .Call(`_RxODE_et_`, list(get.dosing=TRUE),list('last'))"),
			 _["envir"]  = e);
  
  e["get.sampling"] = eval2(_["expr"] = parse2(_["text"] = "function() .Call(`_RxODE_et_`, list(get.sampling=TRUE),list('last'))"),
			    _["envir"]  = e);

  e["getSampling"] = eval2(_["expr"] = parse2(_["text"] = "function() .Call(`_RxODE_et_`, list(get.sampling=TRUE),list('last'))"),
			   _["envir"]  = e);
  e["nobs"] = 0;
  e["ndose"] = 0;
  e[".show"] = LogicalVector::create(_["id"] = false, _["low"] = false,_["time"] = true, _["high"] = false,
				     _["cmt"] =false, _["amt"]=false, _["rate"] = false,
				     _["ii"] = false, _["addl"] = false,
				     _["evid"] = true, _["ss"] = false);
  e[".maxId"] = 1;

  // Return an empty data frame.
  List lst(11);
  CharacterVector nme(11);
  nme[0] = "id";
  lst[0] = IntegerVector(0);
      
  nme[2] = "time";
  lst[2] = NumericVector(0);
      
  nme[1] = "low";
  lst[1] = NumericVector(0);
      
  nme[3] = "high";
  lst[3] = NumericVector(0);
      
  nme[4] = "cmt";
  lst[4] = CharacterVector(0);
      
  nme[5] = "amt";
  lst[5] = NumericVector(0);

  nme[6] = "rate";
  lst[6] = NumericVector(0);
      
  nme[7] = "ii";
  lst[7] = NumericVector(0);
      
  nme[8] = "addl";
  lst[8] = IntegerVector(0);
      
  nme[9] = "evid";
  lst[9] = IntegerVector(0);
      
  nme[10] = "ss";
  lst[10] = IntegerVector(0);

  e.attr("class") = "rxHidden";

  cls.attr(".RxODE.lst") = e;

  lst.attr("names") = nme;
  lst.attr("class") = cls;
  lst.attr("row.names") = IntegerVector::create(NA_INTEGER, 0);
  if (!CharacterVector::is_na(units[1])){
    lst["low"]  = setUnits(lst["low"],  as<std::string>(units[1]));
    lst["time"] = setUnits(lst["time"], as<std::string>(units[1]));
    lst["high"] = setUnits(lst["high"], as<std::string>(units[1]));
  }
  if (!CharacterVector::is_na(units[0])){
    lst["amt"] = setUnits(lst["amt"], as<std::string>(units[0]));
  }
  return lst;
}

List etSort(List curEt){
  std::vector<double> time = as<std::vector<double>>(curEt["time"]);
  std::vector<int> id = as<std::vector<int>>(curEt["id"]);
  std::vector<int> evid = as<std::vector<int>>(curEt["evid"]);
  std::vector<int> idx(id.size());
  std::iota(idx.begin(),idx.end(),0);
  std::sort(idx.begin(),idx.end(),
	    [id,time,evid](int a, int b){
	      if (id[a] == id[b]){
		if (time[a] == time[b]){
		  if (evid[a] == evid[b]){
		    return a < b;
		  }
		  return evid[a] < evid[b];
		}
		return time[a] < time[b];
	      }
	      return id[a] < id[b];
	    });
  List newEt(curEt.size());
  int i, j, newSize = time.size();
  IntegerVector tmpI, tmpI2;
  CharacterVector tmpC, tmpC2;
  NumericVector tmpN, tmpN2;
  for (j = newEt.size(); j--;){
    for (i = newSize;i--;){
      if (rxIs(curEt[j], "numeric")) {
	if (i == newSize-1) newEt[j] = NumericVector(newSize);
	tmpN=newEt[j];
	tmpN2 = curEt[j];
	tmpN[i] = tmpN2[idx[i]];
      } else if (rxIs(curEt[j], "integer")) {
	if (i == newSize-1) newEt[j] = IntegerVector(newSize);
	tmpI=newEt[j];
	tmpI2 = curEt[j];
	tmpI[i] = tmpI2[idx[i]];
      } else if (rxIs(curEt[j], "character")){
	if (i == newSize-1) newEt[j] = CharacterVector(newSize);
	tmpC=newEt[j];
	tmpC2 = curEt[j];
	tmpC[i] = tmpC2[idx[i]];
      } 
    }
  }
  newEt.attr("class") = curEt.attr("class");
  newEt.attr("names") = curEt.attr("names");
  return newEt;
}


List etAddWindow(List windowLst, int idMax, RObject cmt, bool turnOnShowCmt, List curEt){
  std::vector<double> time = as<std::vector<double>>(curEt["time"]);
  std::vector<double> low = as<std::vector<double>>(curEt["low"]);
  std::vector<double> high = as<std::vector<double>>(curEt["high"]);
  std::vector<int> id = as<std::vector<int>>(curEt["id"]);
  int oldSize =id.size();
  std::vector<int> idx(oldSize+windowLst.size()*idMax);
  std::vector<int> evid = as<std::vector<int>>(curEt["evid"]);
  std::iota(idx.begin(),idx.end(),0);
  double c = 0;
  for (int j = idMax; j--;){
    for (int i = windowLst.size(); i--;){
      NumericVector cur = as<NumericVector>(windowLst[i]);
      if (cur.size() != 2)
	stop("Windows need to be a list of observation windows, each of 2 elements e.g. list(c(0,2), c(2,7)).");
      if (cur[0]> cur[1])
	stop("Windows need to be ordered list(c(2,0)) is invalid.");
      id.push_back(j+1);
      low.push_back(cur[0]);
      high.push_back(cur[1]);
      c = Rf_runif(cur[0], cur[1]);
      time.push_back(c);
      evid.push_back(0);
    }
  }
  std::sort(idx.begin(),idx.end(),
	    [id,time,evid](int a, int b){
	      if (id[a] == id[b]){
		if (time[a] == time[b]){
		  if (evid[a] == evid[b]){
		    return a < b;
		  }
		  return evid[a] < evid[b];
		}
		return time[a] < time[b];
	      }
	      return id[a] < id[b];
	    });
  List lst(curEt.size());
  IntegerVector tmpI = as<IntegerVector>(curEt["id"]), tmpI2;
  NumericVector tmpN, tmpN2;
  CharacterVector tmpC, tmpC2;
  lst.attr("names") = curEt.attr("names");
  int i, j;
  // nme[0] = "id";
  lst[0] = IntegerVector(id.size());
      
  // nme[1] = "time";
  lst[2] = NumericVector(id.size());
      
  // nme[2] = "low";
  lst[1] = NumericVector(id.size());
      
  // nme[3] = "high";
  lst[3] = NumericVector(id.size());
      
  // nme[4] = "cmt";
  bool isCmtInt=false;
  if (!rxIs(cmt, "integer")){
    lst[4] = CharacterVector(id.size());
  } else { 
    lst[4] = IntegerVector(id.size());
    isCmtInt=true;
  }
  // nme[5] = "amt";
  lst[5] = NumericVector(id.size());

  // nme[6] = "rate";
  lst[6] = NumericVector(id.size());
      
  // nme[7] = "ii";
  lst[7] = NumericVector(id.size());
      
  // nme[8] = "addl";
  lst[8] = IntegerVector(id.size());
  
  // nme[9] = "evid";
  lst[9] = IntegerVector(id.size());
      
  // nme[10] = "ss";
  lst[10] = IntegerVector(id.size());
  for (i = idx.size(); i--;){
    tmpI = as<IntegerVector>(lst[0]); // id
    tmpI[i] = id[idx[i]];
    
    tmpI = as<IntegerVector>(lst[9]); // evid
    tmpI[i] = evid[idx[i]];

    tmpN = as<NumericVector>(lst[2]); // time
    tmpN[i] = time[idx[i]];

    // low
    tmpN = as<NumericVector>(lst[1]);
    tmpN[i] = low[idx[i]];

    // hi
    tmpN = as<NumericVector>(lst[3]);
    tmpN[i] = high[idx[i]];
    if (idx[i] >= oldSize){
      if (isCmtInt){
	tmpI = as<IntegerVector>(lst[4]);
	tmpI[i] = as<int>(cmt);
      } else {
	tmpC = as<CharacterVector>(lst[4]);
	tmpC2 = as<CharacterVector>(cmt);
	tmpC[i] = tmpC2[0];
      }
      
      // nme[5] = "amt";
      tmpN = as<NumericVector>(lst[5]);
      tmpN[i] = NA_REAL;

      // nme[6] = "rate";
      tmpN = as<NumericVector>(lst[6]);
      tmpN[i] = NA_REAL;
      
      // nme[7] = "ii";
      tmpN = as<NumericVector>(lst[7]);
      tmpN[i] = NA_REAL;
      
      // nme[8] = "addl";
      tmpI = as<IntegerVector>(lst[8]); // id
      tmpI[i] = NA_REAL;
  
      // nme[10] = "ss";
      tmpI = as<IntegerVector>(lst[10]); // id
      tmpI[i] = NA_REAL;
    } else {
      for (j = 11; j--;){
	if (rxIs(curEt[j], "numeric")){
	  tmpN = as<NumericVector>(lst[j]);
	  tmpN2 = as<NumericVector>(curEt[j]);
	  tmpN[i] = tmpN2[idx[i]];
	} else if (rxIs(curEt[j], "integer")) {
	  tmpI = as<IntegerVector>(lst[j]);
	  tmpI2 = as<IntegerVector>(curEt[j]);
	  tmpI[i] = tmpI2[idx[i]];
	} else if (rxIs(curEt[j], "character")){
	  // Char
	  tmpC = as<CharacterVector>(lst[j]);
	  tmpC2 = as<CharacterVector>(curEt[j]);
	  tmpC[i] = tmpC2[idx[i]];	
	}
      }
    }
  }
  CharacterVector cls = clone(as<CharacterVector>(curEt.attr("class")));
  List eOld = cls.attr(".RxODE.lst");
  List e = clone(eOld);
  e["nobs"] = as<int>(e["nobs"]) + windowLst.size()*idMax;
  LogicalVector show = e[".show"];
  if (turnOnShowCmt){
    show["cmt"] = true;
  }
  if (idMax > 1){
    show["id"] = true;
    e[".maxId"] = idMax;
  }
  show["low"] = true;
  show["high"] = true;
  e[".show"] = show;
  e.attr("names") = eOld.attr("names");
  e.attr("class") = eOld.attr("class");
  cls.attr(".RxODE.lst") = e;
  lst.attr("class") = cls;
  int len = as<int>(e["nobs"]) +as<int>(e["ndose"]);
  lst.attr("row.names") = IntegerVector::create(NA_INTEGER, -len);
  return lst;
}

List etAddTimes(NumericVector newTimes, int idMax, RObject cmt, bool turnOnShowCmt, List curEt){
  std::vector<double> time = as<std::vector<double>>(curEt["time"]);
  std::vector<int> id = as<std::vector<int>>(curEt["id"]);
  int oldSize =id.size();
  std::vector<int> idx(oldSize+newTimes.size()*idMax);
  std::vector<int> evid = as<std::vector<int>>(curEt["evid"]);
  std::iota(idx.begin(),idx.end(),0);
  for (int j = idMax; j--;){
    for (int i = newTimes.size(); i--;){
      id.push_back(j+1);
      time.push_back(newTimes[i]);
      evid.push_back(0);
    }
  }
  std::sort(idx.begin(),idx.end(),
	    [id,time,evid](int a, int b){
	      if (id[a] == id[b]){
		if (time[a] == time[b]){
		  if (evid[a] == evid[b]){
		    return a < b;
		  }
		  return evid[a] < evid[b];
		}
		return time[a] < time[b];
	      }
	      return id[a] < id[b];
	    });
  List lst(curEt.size());
  IntegerVector tmpI = as<IntegerVector>(curEt["id"]), tmpI2;
  NumericVector tmpN, tmpN2;
  CharacterVector tmpC, tmpC2;
  lst.attr("names") = curEt.attr("names");
  int i, j;
  // nme[0] = "id";
  lst[0] = IntegerVector(id.size());
      
  // nme[1] = "time";
  lst[2] = NumericVector(id.size());
      
  // nme[2] = "low";
  lst[1] = NumericVector(id.size());
      
  // nme[3] = "high";
  lst[3] = NumericVector(id.size());
      
  // nme[4] = "cmt";
  bool isCmtInt = false;
  if (rxIs(cmt, "integer")){
    lst[4] = IntegerVector(id.size());
    isCmtInt=true;
  } else {
    lst[4] = CharacterVector(id.size());
  }
      
  // nme[5] = "amt";
  lst[5] = NumericVector(id.size());

  // nme[6] = "rate";
  lst[6] = NumericVector(id.size());
      
  // nme[7] = "ii";
  lst[7] = NumericVector(id.size());
      
  // nme[8] = "addl";
  lst[8] = IntegerVector(id.size());
  
  // nme[9] = "evid";
  lst[9] = IntegerVector(id.size());
      
  // nme[10] = "ss";
  lst[10] = IntegerVector(id.size());
  for (i = idx.size(); i--;){
    tmpI = as<IntegerVector>(lst[0]); // id
    tmpI[i] = id[idx[i]];
    
    tmpI = as<IntegerVector>(lst[9]); // evid
    tmpI[i] = evid[idx[i]];

    tmpN = as<NumericVector>(lst[2]); // time
    tmpN[i] = time[idx[i]];
    if (idx[i] >= oldSize){
      // low
      tmpN = as<NumericVector>(lst[1]);
      tmpN[i] = NA_REAL;
      
      // hi
      tmpN = as<NumericVector>(lst[3]);
      tmpN[i] = NA_REAL;

      if (isCmtInt){
	tmpI = as<IntegerVector>(lst[4]);
	tmpI[i] = as<int>(cmt);
      } else {
	tmpC = as<CharacterVector>(lst[4]);
	tmpC2 = as<CharacterVector>(cmt);
	tmpC[i] = tmpC2[0];
      }

      // nme[5] = "amt";
      tmpN = as<NumericVector>(lst[5]);
      tmpN[i] = NA_REAL;

      // nme[6] = "rate";
      tmpN = as<NumericVector>(lst[6]);
      tmpN[i] = NA_REAL;
      
      // nme[7] = "ii";
      tmpN = as<NumericVector>(lst[7]);
      tmpN[i] = NA_REAL;
      
      // nme[8] = "addl";
      tmpI = as<IntegerVector>(lst[8]); // id
      tmpI[i] = NA_REAL;
  
      // nme[10] = "ss";
      tmpI = as<IntegerVector>(lst[10]); // id
      tmpI[i] = NA_REAL;
    } else {
      // low
      for (j = 11; j--;){
	if (rxIs(curEt[j], "numeric")){
	  tmpN = as<NumericVector>(lst[j]);
	  tmpN2 = as<NumericVector>(curEt[j]);
	  tmpN[i] = tmpN2[idx[i]];
	} else if (rxIs(curEt[j], "integer")) {
	  tmpI = as<IntegerVector>(lst[j]);
	  tmpI2 = as<IntegerVector>(curEt[j]);
	  tmpI[i] = tmpI2[idx[i]];
	} else if (rxIs(curEt[j], "character")){
	  // Char
	  tmpC = as<CharacterVector>(lst[j]);
	  tmpC2 = as<CharacterVector>(curEt[j]);
	  tmpC[i] = tmpC2[idx[i]];	
	}
      }
    }
  }
  CharacterVector cls = clone(as<CharacterVector>(curEt.attr("class")));
  List eOld = cls.attr(".RxODE.lst");
  List e = clone(eOld);
  e["nobs"] = as<int>(e["nobs"]) + newTimes.size()*idMax;
  LogicalVector show = e[".show"];
  if (turnOnShowCmt){
    show["cmt"] = true;
  }
  if (idMax > 1){
    show["id"] = true;
    e[".maxId"] = idMax;
  }
  e[".show"] = show;
  e.attr("names") = eOld.attr("names");
  e.attr("class") = eOld.attr("class");
  cls.attr(".RxODE.lst") = e;
  lst.attr("class") = cls;
  int len = as<int>(e["nobs"]) +as<int>(e["ndose"]);
  lst.attr("row.names") = IntegerVector::create(NA_INTEGER, -len);
  return lst;
}


List etResizeId(int maxId, List curEt){
  // Calculate size
  CharacterVector cls = clone(as<CharacterVector>(curEt.attr("class")));
  List eOld = cls.attr(".RxODE.lst");
  List e = clone(eOld);
  int oldMaxId = as<int>(e[".maxId"]);
  if (maxId == oldMaxId) return curEt;
  double c = (double)(maxId)/(double)(oldMaxId);
  int oldSize = as<int>(e["nobs"]) + as<int>(e["ndose"]);
  int newSize = (int)(oldSize*c);
  List newEt(curEt.size());
  IntegerVector tmpI, tmpI2;
  CharacterVector tmpC, tmpC2;
  NumericVector tmpN, tmpN2;
  int i, j;
  bool recalcTime=false;
  if (maxId < oldMaxId){
    // Reducing the number of IDs
    for (j = newEt.size(); j--;){
      if (rxIs(curEt[j], "integer")) {
	tmpI = IntegerVector(newSize);
	tmpI2 = as<IntegerVector>(curEt[j]);
	std::copy(tmpI2.begin(), tmpI2.begin()+newSize, tmpI.begin());
	newEt[j] = tmpI;
      } else if (rxIs(curEt[j], "character")){
	// Char
	tmpC = CharacterVector(newSize);
	tmpC2 = as<CharacterVector>(curEt[j]);
	std::copy(tmpC2.begin(), tmpC2.begin()+newSize, tmpC.begin());
	newEt[j] = tmpC;
      } else {
	tmpN = NumericVector(newSize);
	tmpN2 = as<NumericVector>(curEt[j]);
	std::copy(tmpN2.begin(), tmpN2.begin()+newSize, tmpN.begin());
	newEt[j] = tmpN;
      }
    }
  } else {
    // Enlarge data-set
    int idSize = (int)((double)(oldSize)/(double)(oldMaxId));
    for (j = newEt.size(); j--;){
      if (rxIs(curEt[j], "integer")) {
	tmpI = IntegerVector(newSize);
	tmpI2 = as<IntegerVector>(curEt[j]);
	std::copy(tmpI2.begin(), tmpI2.end(), tmpI.begin());
	if (j == 0){
	  for (i = oldMaxId+1; i <= maxId; i++){
	    std::fill_n(tmpI.begin() + oldSize + (i-oldMaxId-1)*idSize, idSize, i);
	  }
	} else {
	  for (i = newSize - oldSize; i--;){
	    tmpI[oldSize+i] = tmpI2[i % oldSize];
	  }
	}
	newEt[j] = tmpI;
      } else if (rxIs(curEt[j], "character")){
	// Char
	tmpC = CharacterVector(newSize);
	tmpC2 = as<CharacterVector>(curEt[j]);
	std::copy(tmpC2.begin(), tmpC2.end(), tmpC.begin());
	for (i = newSize - oldSize; i--;){
	  tmpC[oldSize+i] = tmpC2[i % oldSize];
	}
	newEt[j] = tmpC;
      } else {
	tmpN = NumericVector(newSize);
	tmpN2 = as<NumericVector>(curEt[j]);
	std::copy(tmpN2.begin(), tmpN2.end(), tmpN.begin());
	for (i = newSize - oldSize; i--;){
	  tmpN[oldSize+i] = tmpN2[i % oldSize];
	}
	newEt[j] = tmpN;
      }
    }
    recalcTime=true;
  }
  newEt.attr("names")     = curEt.attr("names");
  if (recalcTime){
    tmpN = as<NumericVector>(newEt["time"]);
    NumericVector tmpN1 = as<NumericVector>(newEt["low"]);
    tmpN2 = as<NumericVector>(newEt["high"]);
    // Update new observations with recalculated windows
    recalcTime=false;
    for (i = newSize - oldSize; i--;){
      if (!ISNA(tmpN1[oldSize+i]) && !ISNA(tmpN2[oldSize+i])){
	tmpN[oldSize+i] = Rf_runif(tmpN1[oldSize+i], tmpN2[oldSize+i]);
	recalcTime=true;
      }
    }
    curEt = etSort(curEt);
  }
  // Update new windows
  e["nobs"]   = (int)(as<double>(e["nobs"])*c);
  e["ndose"]  = (int)(as<double>(e["ndose"])*c);
  e[".maxId"] = maxId;
  LogicalVector show = e[".show"];
  if (maxId > 1){
    show["id"] = true;
  } else {
    show["id"] = false;
  }
  e[".maxId"]             = maxId;    
  e.attr("class")         = eOld.attr("class");
  cls.attr(".RxODE.lst")  = e;
  newEt.attr("class")     = cls;
  int len = as<int>(e["nobs"]) +as<int>(e["ndose"]);
  newEt.attr("row.names") = IntegerVector::create(NA_INTEGER, -len);
  return newEt;
}

List etAddDose(NumericVector curTime, RObject cmt,  double amt, double rate, double ii,
	       int addl, int curEvid, int ss,
	       int maxId, bool turnOnShowCmt, List curEt){
  std::vector<double> time = as<std::vector<double>>(curEt["time"]);
  std::vector<int> id = as<std::vector<int>>(curEt["id"]);
  std::vector<int> evid = as<std::vector<int>>(curEt["evid"]);
  std::vector<double> low = as<std::vector<double>>(curEt["low"]);
  std::vector<double> high = as<std::vector<double>>(curEt["high"]);
  int oldSize = id.size();
  int i, j;
  double a, b, c;
  for (j = maxId; j--;){
    if (curTime.size() == 1){
      id.push_back(j+1);
      evid.push_back(curEvid);
      time.push_back(curTime[0]);
      low.push_back(NA_REAL);
      high.push_back(NA_REAL);
    } else if (curTime.size() == 2) {
	if (curTime[0] < curTime[1]){
	  id.push_back(j+1);
	  evid.push_back(curEvid);
	  low.push_back(curTime[0]);
	  high.push_back(curTime[1]);
	  c = Rf_runif(curTime[0], curTime[1]);
	  time.push_back(c);
	  for (i = addl; i--;){
	    id.push_back(j+1);
	    evid.push_back(curEvid);
	    a = curTime[0]+ (i+1)*ii;
	    b = curTime[1]+ (i+1)*ii;
	    low.push_back(a);
	    high.push_back(b);
	    c = Rf_runif(a, b);
	    time.push_back(c);
	  }
	} else {
	  stop("For dosing window you need to specify window in order, e.g. et(time=c(0,2),amt=3).");
	}
    } else {
      stop("Time windows must only be 2 elements for dosing.");
    }
    
  }
  std::vector<int> idx(time.size());
  std::iota(idx.begin(),idx.end(),0);
  std::sort(idx.begin(),idx.end(),
	    [id,time,evid](int a, int b){
	      if (id[a] == id[b]){
		if (time[a] == time[b]){
		  if (evid[a] == evid[b]){
		    return a < b;
		  }
		  return evid[a] < evid[b];
		}
		return time[a] < time[b];
	      }
	      return id[a] < id[b];
	    });

  List lst(curEt.size());
  IntegerVector tmpI = as<IntegerVector>(curEt["id"]), tmpI2;
  NumericVector tmpN, tmpN2;
  CharacterVector tmpC, tmpC2;
  lst.attr("names") = curEt.attr("names");
  // nme[0] = "id";
  lst[0] = IntegerVector(id.size());
      
  // nme[1] = "time";
  lst[2] = NumericVector(id.size());
      
  // nme[2] = "low";
  lst[1] = NumericVector(id.size());
      
  // nme[3] = "high";
  lst[3] = NumericVector(id.size());
      
  // nme[4] = "cmt";
  bool isCmtInt=false;
  if (rxIs(cmt,"integer")){
      lst[4] = IntegerVector(id.size());
      isCmtInt=true;
  } else {
      lst[4] = CharacterVector(id.size());
  }
      
  // nme[5] = "amt";
  lst[5] = NumericVector(id.size());

  // nme[6] = "rate";
  lst[6] = NumericVector(id.size());
      
  // nme[7] = "ii";
  lst[7] = NumericVector(id.size());
      
  // nme[8] = "addl";
  lst[8] = IntegerVector(id.size());
  
  // nme[9] = "evid";
  lst[9] = IntegerVector(id.size());
      
  // nme[10] = "ss";
  lst[10] = IntegerVector(id.size());
  for (i = idx.size(); i--;){
    tmpI = as<IntegerVector>(lst[0]); // id
    tmpI[i] = id[idx[i]];
    
    tmpI = as<IntegerVector>(lst[9]); // evid
    tmpI[i] = evid[idx[i]];

    tmpN = as<NumericVector>(lst[1]); // low
    tmpN[i] = low[idx[i]];

    tmpN = as<NumericVector>(lst[2]); // time
    tmpN[i] = time[idx[i]];
    
    tmpN = as<NumericVector>(lst[3]); // high
    tmpN[i] = high[idx[i]];

    if (idx[i] >= oldSize){
      if (isCmtInt){
	tmpI = as<IntegerVector>(lst[4]);
	tmpI[i] = as<int>(cmt);
      } else {
	tmpC = as<CharacterVector>(lst[4]);
	tmpC2 = as<CharacterVector>(cmt);
	tmpC[i] = tmpC2[0];
      }

      // nme[5] = "amt";
      tmpN = as<NumericVector>(lst[5]);
      tmpN[i] = amt;

      // nme[6] = "rate";
      tmpN = as<NumericVector>(lst[6]);
      tmpN[i] = rate;
      
      // nme[7] = "ii";
      tmpN = as<NumericVector>(lst[7]);
      tmpN[i] = ii;
      
      // nme[8] = "addl";
      tmpI = as<IntegerVector>(lst[8]); // id
      tmpI[i] = addl;
  
      // nme[10] = "ss";
      tmpI = as<IntegerVector>(lst[10]); // id
      tmpI[i] = ss;
    } else {
      // low
      for (j = 11; j--;){
	if (rxIs(curEt[j], "numeric")){
	  tmpN = as<NumericVector>(lst[j]);
	  tmpN2 = as<NumericVector>(curEt[j]);
	  tmpN[i] = tmpN2[idx[i]];
	} else if (rxIs(curEt[j], "integer")) {
	  tmpI = as<IntegerVector>(lst[j]);
	  tmpI2 = as<IntegerVector>(curEt[j]);
	  tmpI[i] = tmpI2[idx[i]];
	} else if (rxIs(curEt[j], "character")){
	  // Char
	  tmpC = as<CharacterVector>(lst[j]);
	  tmpC2 = as<CharacterVector>(curEt[j]);
	  tmpC[i] = tmpC2[idx[i]];	
	}
      }
    }
  }
  CharacterVector cls = clone(as<CharacterVector>(curEt.attr("class")));
  List eOld = cls.attr(".RxODE.lst");
  List e = clone(eOld);
  LogicalVector show = e[".show"];
  if (turnOnShowCmt){
    show["cmt"] = true;
  }
  if (maxId > 1){
    show["id"] = true;
    e[".maxId"] = maxId;
  }
  show["amt"] = true;
  if (rate != 0){
    show["rate"] = true;
  }
  if (ss != 0){
    show["ss"] = true;
  }
  if (curTime.size() == 2){
    show["low"] = true;
    show["high"] = true;
    if (addl != 0){
      e["ndose"] = as<int>(e["ndose"]) + maxId*(addl+1);
    }
  } else {
    e["ndose"] = as<int>(e["ndose"]) + maxId;
    if (ii != 0){
      show["ii"] = true;
    }
    if (addl != 0){
      show["addl"] = true;
    }  
  }
  e[".show"] = show;
  e.attr("names") = eOld.attr("names");
  cls.attr(".RxODE.lst") = e;
  lst.attr("class") = cls;
  int len = as<int>(e["nobs"]) +as<int>(e["ndose"]);
  lst.attr("row.names") = IntegerVector::create(NA_INTEGER, -len);
  return lst;
}

RObject etUpdateObj(List curEt, bool update){
  CharacterVector cls=curEt.attr("class");
  List e = clone(as<List>(cls.attr(".RxODE.lst")));
  CharacterVector units = e[".units"];
  if (!CharacterVector::is_na(units[1])){
    curEt["low"]  = setUnits(curEt["low"],  as<std::string>(units[1]));
    curEt["time"] = setUnits(curEt["time"], as<std::string>(units[1]));
    curEt["high"] = setUnits(curEt["high"], as<std::string>(units[1]));
  }
  if (!CharacterVector::is_na(units[0])){
    curEt["amt"] = setUnits(curEt["amt"], as<std::string>(units[0]));
  }
  if (update){
    List cmp = as<List>(evCur);
    for (int j = curEt.size(); j--;){
      cmp[j] = curEt[j];
    }
    cmp.attr("class") = curEt.attr("class");
    cmp.attr("names") = curEt.attr("names");
    cmp.attr("row.names") = curEt.attr("row.names");
  }
  return as<RObject>(curEt);
}

RObject etCmtInt(RObject et){
  List cur = as<List>(et);
  List newEt;
  if (rxIs(cur[4], "character")){
    newEt = clone(cur);
    CharacterVector oldCmt = CharacterVector(cur[4]);
    IntegerVector newCmt = IntegerVector(oldCmt.size());
    for (int j = newCmt.size();j--;){
      if (oldCmt[j] == "(default)") newCmt[j] = 1;
      else if (oldCmt[j] == "(obs)") newCmt[j] = NA_INTEGER;
      else stop("Cannot mix named compartments and integer compartments.");
    }
    warning("Using numbered compartments is discouraged with RxODE simulations.");
    newEt[4] = newCmt;
  } else {
    newEt = cur;
  }
  return (as<RObject>(newEt));
}


//[[Rcpp::export]]
RObject et_(List input, List et__){
  // Create or modify new event table
  double doRet = false;
  bool turnOnShowCmt = false;
  bool doUpdateObj = false;
  RObject curEt;
  if (et__.size() > 0){
    if (rxIs(et__[0],"character")){
      if (as<std::string>(et__[0]) == "last"){
	doUpdateObj=true;
	curEt = evCur;
      }
    } else if (rxIs(et__, "rxEt")) {
      curEt = as<RObject>(et__);
    }
  }
  CharacterVector inN;
  if (input.hasAttribute("names")){
    inN = input.attr("names");
  }
  int i, amtIx = -1, iiIx = -1, addlIx = -1,
    untilIx = -1, evidIx=-1, idIx=-1, cmtIx=-1, 
    amtUnitIx=-1, timeUnitIx=-1, doSamplingIdx=-1, timeIx=-1,
    rateIx = -1, nbrIx=-1, ssIx=-1;
  // Wait should be in sequences and rep
  for (i = (int)inN.size(); i--;){
    if (inN[i] == "amt" || inN[i] == "dose") amtIx=i;
    else if (inN[i] == "ii" || inN[i] == "dosing.interval" || inN[i] == "dosingInterval" || inN[i] == "dosing_interval") iiIx=i;
    else if (inN[i] == "addl") addlIx = i;
    else if (inN[i] == "until") untilIx = i;
    else if (inN[i] == "evid") evidIx = i;
    else if (inN[i] == "ID" || inN[i] == "id") idIx=i;
    else if (inN[i] == "cmt" || inN[i] == "dosing.to" || inN[i] == "dosingTo" || inN[i] =="dosing_to" ||
	     inN[i] == "dose.to" || inN[i] == "doseTo" || inN[i] == "dose_to") cmtIx=i;
    else if (inN[i] == "amount.units" || inN[i] == "amountUnits" || inN[i] == "amount_units") amtUnitIx=i;
    else if (inN[i] == "time.units" || inN[i] == "timeUnits" || inN[i] == "time_units") timeUnitIx=i;
    else if (inN[i] == "do.sampling" || inN[i] == "doSampling" || inN[i] == "do_sampling") doSamplingIdx=i;
    else if (inN[i] == "start.time" || inN[i] == "startTime" || inN[i] == "start_time" ||
	     inN[i] == "start" || inN[i] == "time") timeIx = i;
    else if (inN[i] == "nbr.doses" || inN[i] == "nbrDoses" || inN[i] == "nbr") nbrIx=i;
    else if (inN[i] == "ss") ssIx = i;
    else if (inN[i] == "rate") rateIx = i;
  }
  // missing argument name handling.
  for (i = 0; i <(int)inN.size(); i++){
    if (inN[i] == ""){
      if (rxIs(input[i], "character") && cmtIx == -1) cmtIx = i;
      if (rxIs(input[i], "numeric") && timeIx == -1) timeIx = i;
      if (rxIs(input[i], "integer") && timeIx == -1) timeIx = i;
      if (rxIs(input[i], "list") && timeIx == -1) timeIx=i;
      if (rxIs(input[i], "rxEt")) curEt = input[i];
    }
  }
  if (inN.size() == 0 && input.size() != 0){
    for (i = 0; i <(int)input.size(); i++){
      if (rxIs(input[i], "character")) cmtIx = i;
      if (rxIs(input[i], "numeric") && timeIx == -1) timeIx = i;
      if (rxIs(input[i], "integer") && timeIx == -1) timeIx = i;
      if (rxIs(input[i], "list") && timeIx == -1) timeIx=i;
      if (rxIs(input[i], "rxEt")) curEt = input[i];
    }
  }
  if (rxIs(curEt, "rxEt")){
    // This is a modification to an existing event table.
    if (rxIs(input[0], "logical")){
      LogicalVector in0 = as<LogicalVector>(input[0]);
      if (in0[0]){
	CharacterVector cls = curEt.attr("class");
	List e = cls.attr(".RxODE.lst");
	CharacterVector nm = input.attr("names");
	if (nm[0] == "copy"){
	  // There is no need to copy any more.
	  return curEt;
	} else if (nm[0] == "get.EventTable"){
	  e.attr("class") = R_NilValue;
	  if (as<int>(e["nobs"]) == 0 && as<int>(e["ndose"]) == 0){
	    return R_NilValue;
	  } else {
	    List ret = clone(as<List>(curEt));
	    ret.attr("class") = "data.frame";
	    return as<RObject>(ret);
	  }
	} else if (nm[0] == "get.obs.rec"){
	} else if (nm[0] == "get.sampling" || nm[0] == "get.dosing" ||
		   nm[0] == "clearSampling" || nm[0] == "clearDosing"){
	  // Need to update
	  bool doDose = (nm[0] == "get.dosing" || nm[0] == "clearSampling");
	  bool updateObj = (nm[0] == "clearSampling" || nm[0] == "clearDosing");
	  int n = doDose ? as<int>(e["ndose"]) : as<int>(e["nobs"]);
	  List cmp;
	  if (n == 0){
	    if (updateObj){
	      cmp = as<List>(curEt);
	      List em =etEmpty(as<CharacterVector>(e[".units"]));
	      cls = em.attr("class");
	      for (int j = cmp.size(); j--;){
		if (rxIs(cmp[j], "numeric")){
		  cmp[j] = NumericVector(0);
		} else if (rxIs(cmp[j], "integer")){
		  cmp[j] = IntegerVector(0);
		} else if (rxIs(cmp[j], "character")) {
		  cmp[j] = CharacterVector(0);
		}
	      }
	      cmp.attr("row.names") = IntegerVector::create(NA_INTEGER, -n);
	      cmp.attr("class") = cls;
	      return R_NilValue;
	    } else {
	      return R_NilValue;
	    }
	  } else {
	    cmp = as<List>(curEt);
	    List ret(cmp.size());
	    IntegerVector evid = as<IntegerVector>(cmp["evid"]);
	    int i, j, k;
	    for (j = ret.size(); j--;){
	      if (rxIs(cmp[j], "numeric")){
		ret[j] = NumericVector(n);
	      } else if (rxIs(cmp[j], "integer")){
		ret[j] = IntegerVector(n);
	      } else if (rxIs(cmp[j], "character")){
		ret[j] = CharacterVector(n);
	      }
	    }
	    k = 0;
	    NumericVector tmpN, tmpN2;
	    IntegerVector tmpI, tmpI2;
	    CharacterVector tmpC, tmpC2;
	    for (i = 0; i < evid.size(); i++){
	      if ((doDose && evid[i] != 0) || (!doDose && evid[i] == 0)){
		for (j = ret.size(); j--;){
		  if (rxIs(cmp[j], "numeric")){
		    tmpN = ret[j];
		    tmpN2 = cmp[j];
		    tmpN[k] = tmpN2[i];
		  } else if (rxIs(cmp[j], "integer")){
		    tmpI = ret[j];
		    tmpI2 = cmp[j];
		    tmpI[k] = tmpI2[i];
		  } else if (rxIs(cmp[j], "character")){
		    tmpC = ret[j];
		    tmpC2 = cmp[j];
		    tmpC[k] = tmpC2[i];
		  }
		}
		k++;
	      }
	    }
	    ret.attr("names") = cmp.attr("names");
	    ret.attr("row.names") = IntegerVector::create(NA_INTEGER, -n);
	    if (updateObj){
	      // This updates the object in place.
	      if (doDose){
		e["nobs"] = 0;
	      } else {
		e["ndose"] = 0;
	      }
	      cls.attr(".RxODE.lst") = e;
	      for (j = cmp.size(); j--;){
		cmp[j] = ret[j];
	      }
	      cmp.attr("row.names") = IntegerVector::create(NA_INTEGER, -n);
	      cmp.attr("class") = cls;
	      return as<RObject>(ret);
	    } else {
	      ret.attr("class") = "data.frame";
	    }
	    return as<RObject>(ret);
	  }
	} else {
	  warning("Nothing done");
	  return  as<List>(curEt);
	}
      } else {
	warning("Nothing done");
	return  as<List>(curEt);
      }
    } else {
      // We are updating the event table
      CharacterVector cls = curEt.attr("class");
      List e = cls.attr(".RxODE.lst");
      int id=1;
      if (idIx != -1){
	id    = as<int>(input[idIx]);
	if ((int)(e["nobs"]) != 0 || (int)(e["ndose"]) != 0){
	  if ((int)(e[".maxId"]) != id){
	    curEt = as<RObject>(etResizeId(id, as<List>(curEt)));
	    doRet=true;
	  }
	}
      } else {
	id = (int)(e[".maxId"]);
      }
      CharacterVector cmtS;
      IntegerVector cmtI;
      RObject cmt;
      // Dose
      if (cmtIx == -1){
	List tmp = as<List>(curEt);
	if (rxIs(tmp[4], "integer")){
	  cmtI = IntegerVector(1);
	  if (amtIx == -1){
	    cmtI[0] = NA_INTEGER;
	  } else {
	    cmtI[0] = 1;
	  }
	  cmt = as<RObject>(cmtI);
	} else {
	  cmtS = CharacterVector(1);
	  if (amtIx == -1){
	    cmtS[0] = "(obs)";
	  } else {
	    cmtS[0] = "(default)";
	  }
	  cmt = as<RObject>(cmtS);
	}
      } else {
	if (rxIs(input[cmtIx], "character")){
	  cmtS = as<CharacterVector>(input[cmtIx]);
	  if (cmtS.size() == 1){
	    List tmp = as<List>(curEt);
	    if (rxIs(tmp[4], "integer")){
	      stop("Cannot mix named and integer compartments");
	    }
	    cmt = as<RObject>(cmtS);
	  } else {
	    stop("The compartment cannot be a vector.");
	  }
	} else if (rxIs(input[cmtIx], "integer") || rxIs(input[cmtIx], "numeric")){
	  cmtI = as<IntegerVector>(input[cmtIx]);
	  if (cmtI.size() == 1){
	    curEt=etCmtInt(curEt);
	    cls = curEt.attr("class");
	    e = cls.attr(".RxODE.lst");
	    cmt = as<RObject>(cmtI);
	  } else {
	    stop("The compartment cannot be an integer.");
	  }
	} else {
	  stop("The compartment must be an integer or a character.");
	}
	turnOnShowCmt=true;
      }
      NumericVector amt;
      bool isObs=false;
      if (amtIx == -1){
	isObs=true;
      } else {
	amt = as<NumericVector>(input[amtIx]);
	if (amt.size() != 1){
	  stop("Dose amount cannot be a vector.");
	}
	isObs = false;
      }
      if (isObs){
	IntegerVector addl;// = 0;
	if (addlIx != -1){
	  addl = as<IntegerVector>(input[addlIx]);
	  if (addl.size() != 1) stop("addl cannot be a vector.");
	  if (addl[0] != 0){
	    stop("addl needs a dose/amt.");
	  }
	} else {
	  addl = IntegerVector(1);
	  addl[0] = 0;
	}
	if (untilIx != -1){
	  stop("until needs a dose/amt.");
	}
	if (nbrIx != -1){
	  stop("nbr.doses needs a dose/amt.");
	}
	NumericVector rate;
	if (rateIx != -1){
	  rate = as<NumericVector>(input[rateIx]);
	  if (rate.size() != 1) stop("rate cannot be a vector");
	  if (rate[0] != 0.0){
	    stop("rate needs a dose/amt.");
	  }
	} else {
	  rate = NumericVector(1);
	  rate[0]=0.0;
	}
	NumericVector ii;
	if (iiIx != -1){
	  ii = as<NumericVector>(input[iiIx]);
	  if (ii.size() != 1) stop("ii cannot be a vector.");
	  if (ii[0] != 0.0){
	    stop("ii needs a dose/amt.");
	  }
	} else {
	  ii = NumericVector(1);
	  ii[0] = 0.0;
	}
	IntegerVector evid;
	if (evidIx != -1){
	  evid = as<IntegerVector>(input[evidIx]);
	  if (evid.size() != 1){
	    stop("evid cannot be a vector.");
	  }
	  if (evid[0] != 0){
	    stop("non-zero evid needs a dose/amt.");
	  }
	} else {
	  evid = IntegerVector(1);
	  evid[0]=0;
	}
	IntegerVector ss;// = 0;
	if (ssIx != -1){
	  ss = as<IntegerVector>(input[ssIx]);
	  if (ss.size() != 1) stop("ss cannot be a vector.");
	  if (ss[0] != 0){
	    stop("non-zero ss needs a dose/amt.");
	  }
	}
	if (timeIx != -1) {
	  if (rxIs(input[timeIx], "numeric") || rxIs(input[timeIx], "integer")){
	    NumericVector time = as<NumericVector>(input[timeIx]);
	    return etUpdateObj(etAddTimes(as<NumericVector>(input[timeIx]), id, cmt, turnOnShowCmt, as<List>(curEt)),
			       doUpdateObj);
	  } else if (rxIs(input[timeIx], "list")){
	    return etUpdateObj(etAddWindow(as<List>(input[timeIx]), id, cmt, turnOnShowCmt, as<List>(curEt)),
			       doUpdateObj);
	  }
	}
      } else {
	////////////////////////////////////////////////////////////////////////////////
	// Dose
	////////////////////////////////////////////////////////////////////////////////
	if (addlIx != -1 && untilIx != -1){
	  stop("Can only specify until or addl, not both.");
	}
	if (addlIx != -1 && nbrIx != -1){
	  stop("Can only specify addl or nbr.doses, not both.");
	}
	if (nbrIx != -1 && untilIx != -1){
	  stop("Can only specify nbr.doses or until, not both.");
	}
	NumericVector rate;
	if (rateIx != -1){
	  rate = as<NumericVector>(input[rateIx]);
	  if (rate.size() != 1) stop("rate cannot be a vector.");
	} else {
	  rate = NumericVector(1);
	  rate[0] = 0.0;
	}
	NumericVector ii;// =0.0;
	if (iiIx != -1){
	  ii = as<NumericVector>(input[iiIx]);
	  if (ii.size() != 1) stop("ii cannot be a vector.");
	} else {
	  ii = NumericVector(1);
	  ii[0] = 0.0;
	}
	IntegerVector evid;// =1;
	if (evidIx != -1){
	  evid = as<IntegerVector>(input[evidIx]);
	  if (evid.size()!= 1){
	    stop("evid cannot be a vector");
	  }
	  if (evid[0] == 0){
	    stop("zero evid cannot be used with dose/amt.");
	  }
	} else {
	  evid = IntegerVector(1);
	  evid[0] = 1;
	}
	IntegerVector ss;// = 0;
	if (ssIx != -1){
	  ss = as<IntegerVector>(input[ssIx]);
	  if (ss.size() != 1) stop("ss cannot be a vector.");
	} else {
	  ss = IntegerVector(1);
	  ss[0] = 0.0;
	}
	NumericVector time;
	if (timeIx != -1){
	  time = as<NumericVector>(input[timeIx]);
	} else {
	  time = NumericVector(1);
	  time[0] = 0;
	}
	IntegerVector addl;//=0;
	if (addlIx != -1){
	  addl = as<IntegerVector>(input[addlIx]);
	  if (addl.size() != 1) stop("addl cannot be a vector.");
	} else if (nbrIx != -1){
	  addl = as<IntegerVector>(input[nbrIx]) - 1;
	  if (addl.size() != 1) stop("Number of doses cannot be a vector.");
	  if (addl[0] < 0){
	    stop("Number of Doses must be at least one.");
	  }
	} else if (untilIx != -1){
	  // Need time for this
	  NumericVector until = as<NumericVector>(input[untilIx]);
	  if (until.size() != 1) stop("Until cannot be a vector.");
	  if (ii[0] < 0){
	    stop("'until' can only be used with positive inter-dose intervals (ii).");
	  }
	  if (time.size() == 1){
	    double tmp = until[0] - time[0] - ii[0];
	    if (tmp > 0){
	      addl[0] = ceil(tmp/ii[0]);
	    } else {
	      addl[0] = 0;
	    }
	  } else if (time.size() == 2){
	    double tmp = until[0] - time[1] - ii[0];
	    if (tmp > 0){
	      addl[0] = ceil(tmp/ii[0]);
	    } else {
	      addl[0] = 0;
	    }
	  }
	}	
	if (ii[0] > 0 && ss[0] == 0 && addl[0] == 0){
	  stop("ii requires non zero additional doses (addl/until/nbr.doses) or steady state dosing.");
	}
	if (ss[0] < 0 || ss[0] > 2){
	  stop("ss must be 0, 1 or 2.");
	}
	if (ss[0] > 1 && time.size() > 1){
	  stop("Steady state (ss) is not supported with dosing windows.");
	}
	if (addl[0] < 0){
	  stop("Additional doses must be positive.");
	}
	return etUpdateObj(etAddDose(time, cmt, amt[0], rate[0], ii[0], addl[0], evid[0], ss[0],
				     id, turnOnShowCmt, as<List>(curEt)),doUpdateObj);
      }
    }
  } else {
    // This is a new event table.
    // Get the amtUnit
    CharacterVector units(2);
    int foundArgs = 0;
    if (amtUnitIx == -1){
      units[0] = NA_STRING;
    } else  {
      units[0] = as<std::string>(input[amtUnitIx]);
      foundArgs++;
    }
    if (timeUnitIx == -1){
      units[1] = "hours";
    } else  {
      units[1] = as<std::string>(input[timeUnitIx]);
      foundArgs++;
    }
    units.attr("names") = CharacterVector::create("dosing", "time");
    List et = etEmpty(units);
    if (input.size() == foundArgs){
      return et;
    } else {
      return et_(input, et);
    }
  }
  if (doRet) return curEt;
  stop("Cannot figure out what type of EventTable you are trying to create.");
  // Should never get here...
  List ret(0);
  return ret;
}
