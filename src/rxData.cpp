// [[Rcpp::depends(RcppArmadillo)]]
#define NCMT 100
// NONMEM 7.1 has a max of 50 obesrvations/individual
#define MAXIDS 500
#define NALL 500
#define NDOSES 50
// NONMEM nTHETA=20
// NONMEM nETA=30
// NONMEM nSIGMA=10
#define NPARS 60
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <thread>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <climits>
extern "C" {
#include "RxODE.h"
}
#include "ode.h"
#define rxModelVars(a) rxModelVars_(a)
using namespace Rcpp;
using namespace arma;

int rxcEvid = -1;
int rxcTime = -1;
int rxcAmt  = -1;
int rxcId   = -1;
int rxcDv   = -1;
int rxcLen  = -1;
bool resetCache = true;
bool rxHasEventNames(CharacterVector &nm){
  int len = nm.size();
  bool reset  = resetCache;
  if (reset || len != rxcLen){
    reset   = resetCache;
    rxcEvid = -1;
    rxcTime = -1;
    rxcAmt  = -1;
    rxcId   = -1;
    rxcDv   = -1;
    rxcLen  = len;
    for (unsigned int i = len; i--;){
      if (as<std::string>(nm[i]) == "evid" || as<std::string>(nm[i]) == "EVID" || as<std::string>(nm[i]) == "Evid"){
        rxcEvid = i;
      } else if (as<std::string>(nm[i]) == "time" || as<std::string>(nm[i]) == "TIME" || as<std::string>(nm[i]) == "Time"){
        rxcTime = i;
      } else if (as<std::string>(nm[i]) == "amt" || as<std::string>(nm[i]) == "AMT" || as<std::string>(nm[i]) == "Amt"){
        rxcAmt = i;
      } else if (as<std::string>(nm[i]) == "id" || as<std::string>(nm[i]) == "ID" || as<std::string>(nm[i]) == "Id"){
        rxcId = i;
      } else if (as<std::string>(nm[i]) == "dv" || as<std::string>(nm[i]) == "DV" || as<std::string>(nm[i]) == "Dv"){
        rxcDv = i;
      }
    }
  }
  resetCache = true;
  if (rxcEvid >= 0 && rxcTime >= 0 && rxcAmt >= 0){
    return true;
  } else {
    return false;
  }
}

//' Check the type of an object using Rcpp
//'
//' @param obj Object to check
//' @param cls Type of class.  Only s3 classes for lists/environments and primitive classes are checked.
//'    For matrix types they are distinguished as \code{numeric.matrix}, \code{integer.matrix},
//'    \code{logical.matrix}, and \code{character.matrix} as well as the traditional \code{matrix}
//'    class. Additionally checks for \code{event.data.frame} which is an \code{data.frame} object
//'    with \code{time},  \code{evid} and \code{amt}. (UPPER, lower or Title cases accepted)
//'
//' @return A boolean indicating if the object is a member of the class.
//' @keywords internal
//' @author Matthew L. Fidler
//' @export
// [[Rcpp::export]]
bool rxIs(const RObject &obj, std::string cls){
  int type = obj.sexp_type();
  bool hasDim = false;
  bool hasCls = false;
  switch (type){
  case 0: return (cls == "NULL");
  case REALSXP: 
    hasDim = obj.hasAttribute("dim");
    if (hasDim){
      if (cls == "event.matrix" || cls ==  "rx.event"){
	if (obj.hasAttribute("dimnames")){
	  List dn = as<List>(obj.attr("dimnames"));
          if (dn.size() == 2){
            CharacterVector cv = as<CharacterVector>(dn[1]);
            return rxHasEventNames(cv);
          } else {
	    return false;
	  }
	} else {
	  return false;
	}
      } else {
	return (cls == "matrix" || cls == "numeric.matrix");
      }
    } else {
      return (cls == "numeric");
    }
  case 13: // integer vectors
    // An integer vector cannot be an event matrix.
    hasDim = obj.hasAttribute("dim");
    if (hasDim){
      return (cls == "matrix" || cls == "integer.matrix");
    } else {
      return (cls == "integer");
    }
  case LGLSXP:
    hasDim = obj.hasAttribute("dim");
    if (hasDim){
      return (cls == "matrix" ||  cls == "logical.matrix");
    } else {
      return (cls == "logical");
    }
  case STRSXP:
    hasDim = obj.hasAttribute("dim");
    if (hasDim){
      return (cls == "matrix" || cls == "character.matrix");
    } else {
      return (cls == "character");
    }
  case VECSXP:
    hasCls = obj.hasAttribute("class");
    if (hasCls){
      CharacterVector classattr = obj.attr("class");
      bool hasDf = false;
      bool hasEt = false;
      std::string cur;
      for (unsigned int i = classattr.size(); i--; ){
	cur = as<std::string>(classattr[i]);
	if (cur == cls){
	  if (cls == "rxSolve"){
	    Environment e = as<Environment>(classattr.attr(".RxODE.env"));
	    List lobj = List(obj);
	    CharacterVector cls2= CharacterVector::create("data.frame");
	    if (as<int>(e["check.ncol"]) != lobj.size()){
	      lobj.attr("class") = cls2;
	      return false;
	    }
	    int nrow = (as<NumericVector>(lobj[0])).size();
	    if (as<int>(e["check.nrow"]) != nrow){
	      lobj.attr("class") = cls2;
	      return false;
	    }
	    CharacterVector cn = CharacterVector(e["check.names"]);
	    if (cn.size() != lobj.size()){
	      lobj.attr("class") = cls2;
	      return false;
	    }
	    CharacterVector cn2 = CharacterVector(lobj.names());
	    for (int j = 0; j < cn.size();j++){
	      if (cn[j] != cn2[j]){
		lobj.attr("class") = cls2;
		return false;
	      }
	    }
	    return true;
	  } else {
	    return true;
	  }
	} else if (cur == "data.frame"){
	  hasDf=true;
        } else if (cur == "EventTable"){
	  hasEt=true;
	}
      }
      if (hasDf && (cls == "rx.event" || cls == "event.data.frame")){
	// Check for event.data.frame
	CharacterVector cv =as<CharacterVector>((as<DataFrame>(obj)).names());
	return rxHasEventNames(cv);
      } else if (hasEt) {
	return (cls == "rx.event");
      } else {
	return false;
      }
    } else {
      return (cls == "list");
    }
  case 4: // environment
    hasCls = obj.hasAttribute("class");
    if (hasCls){
      CharacterVector classattr = obj.attr("class");
      std::string cur;
      for (unsigned int i = classattr.size(); i--; ){
        cur = as<std::string>(classattr[i]);
	if (cur == cls) return true;
      }
    } else if (cls == "environment"){
       return true;
    }
    return false;
  case 22: // external pointer
    return (cls == "externalptr" || cls == "refObject");
  }
  return false;
}

Function loadNamespace("loadNamespace", R_BaseNamespace);
Environment mvnfast = loadNamespace("mvnfast");

RObject rxSimSigma(const RObject &sigma,
		   const RObject &df,
		   int ncores,
		   const bool &isChol,
		   int nObs,
		   const bool checkNames = true){
  if (nObs < 1) stop("Refusing to simulate %d items",nObs); 
  if (rxIs(sigma, "numeric.matrix")){
    // FIXME more distributions
    NumericMatrix sigmaM(sigma);
    if (sigmaM.nrow() != sigmaM.ncol()){
      stop("The matrix must be a square matrix.");
    }
    List dimnames;
    StringVector simNames;
    bool addNames = false;
    if (checkNames){
      if (!sigmaM.hasAttribute("dimnames")){
        stop("The matrix must have named dimensions.");
      }
      dimnames = sigmaM.attr("dimnames");
      simNames = as<StringVector>(dimnames[1]);
      addNames = true;
    } else if (sigmaM.hasAttribute("dimnames")){
      dimnames = sigmaM.attr("dimnames");
      simNames = as<StringVector>(dimnames[1]);
      addNames = true;
    }
    NumericMatrix simMat(nObs,sigmaM.ncol());
    NumericVector m(sigmaM.ncol());
    // Ncores = 1?  Should it be parallelized when it can be...?
    // Note that if so, the number of cores also affects the output.
    if (df.isNULL()){
      Function rmvn = as<Function>(mvnfast["rmvn"]);
      rmvn(_["n"]=nObs, _["mu"]=m, _["sigma"]=sigmaM, _["ncores"]=ncores,
	   _["isChol"]=isChol, _["A"] = simMat); // simMat is updated with the random deviates
    } else {
      double df2 = as<double>(df);
      if (R_FINITE(df2)){
        Function rmvt = as<Function>(mvnfast["rmvt"]);
        rmvt(_["n"]=nObs, _["mu"]=m, _["sigma"]=sigmaM, _["df"] = df,
	     _["ncores"]=ncores, _["isChol"]=isChol, _["A"] = simMat);
      } else {
        Function rmvn = as<Function>(mvnfast["rmvn"]);
        rmvn(_["n"]=nObs, _["mu"]=m, _["sigma"]=sigmaM, _["ncores"]=ncores,
	     _["isChol"]=isChol, _["A"] = simMat);
      }
    }
    if (addNames){
      simMat.attr("dimnames") = List::create(R_NilValue, simNames);
    }
    return wrap(simMat);
  } else {
    return R_NilValue;
  }
}

bool foundEnv = false;
Environment _rxModels;
bool _RxODE_found = false;
Environment _RxODE;

Environment RxODE(){
  if (_RxODE_found){
    return _RxODE;
  } else {
    Function loadNamespace("loadNamespace", R_BaseNamespace);
    _RxODE = loadNamespace("RxODE");
    _RxODE_found = true;
    return _RxODE;
  }
}
Function getRxFn(std::string name){
  Environment rx = RxODE();
  return as<Function>(rx[name]);
}

void getRxModels(){
  if (!foundEnv){ // minimize R call
    Function f = getRxFn("rxModels_");
    _rxModels = f();
    foundEnv = true;
  }
}

// 
inline bool fileExists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

// Use this function to keep dynLoad options consistent.
//[[Rcpp::export]]
SEXP dynLoad(std::string dll){
  Function dl("dyn.load", R_BaseNamespace);
  SEXP ret = dl(dll, _["local"]=false, _["now"]=true);
  return ret;
}


// [[Rcpp::export]]
List rxModelVars_(const RObject &obj){
  getRxModels();
  if (rxIs(obj, "rxModelVars")){
    List ret(obj);
    return ret;
  } else if (rxIs(obj,"RxODE")) {
    Environment e = as<Environment>(obj);
    List rxDll = e["rxDll"];
    List ret = rxDll["modVars"];
    return ret;
  } else if (rxIs(obj,"rxSolve")){
    CharacterVector cls = obj.attr("class");
    Environment e = as<Environment>(cls.attr(".RxODE.env"));
    return  rxModelVars_(as<RObject>(e["args.object"]));
  } else if (rxIs(obj,"rxDll")){
    List lobj = (as<List>(obj))["modVars"];
    return lobj;
  } else if (rxIs(obj, "environment")){
    Environment e = as<Environment>(obj);
    if (e.exists("args.object")){
      return rxModelVars_(e["args.object"]);
    } else {
      CharacterVector cls = obj.attr("class");
      int i = 0;
      Rprintf("Class:\t");
      for (i = 0; i < cls.size(); i++){
        Rprintf("%s\t", (as<std::string>(cls[i])).c_str());
      }
      Rprintf("\n");
      stop("Need an RxODE-type object to extract model variables from.");
    }
  } else if (rxIs(obj, "character")){
    CharacterVector modList = as<CharacterVector>(obj);
    if (modList.size() == 1){
      std::string sobj =as<std::string>(obj);
      if ((sobj.find("=") == std::string::npos) &&
	  (sobj.find("<-") == std::string::npos) &&
          (sobj.find("~") == std::string::npos)){
        if (_rxModels.exists(sobj)){
          RObject obj1 = _rxModels.get(sobj);
          if (rxIs(obj1, "rxModelVars")){
            return as<List>(obj1);
          } else if (rxIs(obj1, "RxODE")){
            return rxModelVars_(obj1);
          }
        }
        std::string sobj1 = sobj + "_model_vars";
        if (_rxModels.exists(sobj1)){
          RObject obj1 = _rxModels.get(sobj1);
          if (rxIs(obj1, "rxModelVars")){
            return as<List>(obj1);
          }
        }
        Function get("get",R_BaseNamespace);
        List platform = get(_["x"]=".Platform", _["envir"] = R_BaseEnv);
        sobj1 = sobj + "_" + as<std::string>(platform["r_arch"]) + "_model_vars";
        if (_rxModels.exists(sobj1)){
          RObject obj1 = _rxModels.get(sobj1);
          if (rxIs(obj1, "rxModelVars")){
            return as<List>(obj1);
          }
        }
        Function filePath("file.path", R_BaseNamespace);
        Function getwd("getwd", R_BaseNamespace);
        sobj1 = as<std::string>(getwd());
        std::string sobj2 = sobj + ".d";
        std::string sobj3 = sobj + "_" + as<std::string>(platform["r_arch"]) +
          as<std::string>(platform["dynlib.ext"]);
        sobj1 = as<std::string>(filePath(sobj1,sobj2, sobj3));
        if (fileExists(sobj1)){
          Rcout << "Path: " << sobj1 << "\n";
          dynLoad(sobj1);
          sobj1 = sobj + "_" + as<std::string>(platform["r_arch"]) +
            "_model_vars";
          Function call(".Call", R_BaseNamespace);
          List ret = as<List>(call(sobj1));
          return ret;
        }
      }
    } else if (modList.hasAttribute("names")){
      bool containsPrefix = false;
      CharacterVector modListNames = modList.names();
      for (int i = 0; i < modListNames.size(); i++){
	if (modListNames[i] == "prefix"){
	  containsPrefix=true;
	  break;
	}
      }
      if (containsPrefix){
	std::string mvstr = as<std::string>(modList["prefix"]) + "model_vars";
        if(_rxModels.exists(mvstr)){
          RObject obj1 = _rxModels.get(mvstr);
          if (rxIs(obj1, "rxModelVars")){
            return as<List>(obj1);
          }
        }
      }
    }
    // fileExists(const std::string& name)
    Function f = getRxFn(".rxModelVarsCharacter");
    return f(obj);
  } else if (rxIs(obj,"list")){
    bool params=false, lhs=false, state=false, trans=false, ini=false, model=false, md5=false, podo=false, dfdy=false;
    List lobj  = as<List>(obj);
    CharacterVector nobj = lobj.names();
    for (int i = 0; i < nobj.size(); i++){
      if (nobj[i] == "modVars"){
	return(rxModelVars_(lobj["modVars"]));
      } else if (!params && nobj[i]== "params"){
	params=true;
      } else if (!lhs && nobj[i] == "lhs"){
	lhs=true;
      } else if (!state && nobj[i] == "state"){
	state=true;
      } else if (!trans && nobj[i] == "trans"){
	trans=true;
      } else if (!ini && nobj[i] == "ini"){
	ini = true;
      } else if (!model && nobj[i] == "model"){
	model = true;
      } else if (!md5 && nobj[i] == "md5"){
	md5 = true;
      } else if (!podo && nobj[i] == "podo"){
	podo=true;
      } else if (!dfdy && nobj[i] == "dfdy"){
	dfdy = true;
      } else {
        return lobj;
      }
    }
    stop("Cannot figure out the model variables.");
  } else if (rxIs(obj,"NULL")) {
      stop("A NULL object does not have any RxODE model variables");
  } else {
    CharacterVector cls = obj.attr("class");
    int i = 0;
    Rprintf("Class:\t");
    for (i = 0; i < cls.size(); i++){
      Rprintf("%s\t", (as<std::string>(cls[i])).c_str());
    }
    Rprintf("\n");
    stop("Need an RxODE-type object to extract model variables from.");
  }
}


//' State variables
//'
//' This returns the model's compartments or states.
//'
//' @inheritParams rxModelVars
//'
//' @param state is a string indicating the state or compartment that
//'     you would like to lookup.
//'
//' @return If state is missing, return a character vector of all the states.
//'
//' If state is a string, return the compartment number of the named state.
//'
//' @seealso \code{\link{RxODE}}
//'
//' @author Matthew L.Fidler
//' @export
// [[Rcpp::export]]
RObject rxState(const RObject &obj = R_NilValue, RObject state = R_NilValue){
  List modVar = rxModelVars(obj);
  CharacterVector states = modVar["state"];
  if (state.isNULL()){
    return states;
  }
  else if (rxIs(state,"character")){
    CharacterVector lookup = as<CharacterVector>(state);
    if (lookup.size() > 1){
      // Fixme?
      stop("Can only lookup one state at a time.");
    }
    if (states.size() == 1){
      warning("Only one state variable should be input.");
    }
    IntegerVector ret(1);
    for (int i = 0; i < states.size(); i++){
      if (states[i] == lookup[0]){
	ret[0] = i+1;
	return ret;
      }
    }
    stop("Cannot locate compartment \"%s\".",as<std::string>(lookup[0]).c_str());
  }
  return R_NilValue;
}

//[[Rcpp::export]]
CharacterVector rxParams_(const RObject &obj){
  List modVar = rxModelVars(obj);
  CharacterVector ret = modVar["params"];
  return ret;
}


//' Jacobian and parameter derivatives
//'
//' Return Jacobain and parameter derivatives
//'
//' @inheritParams rxModelVars
//'
//' @return A list of the jacobian parameters defined in this RxODE
//'     object.
//' @author Matthew L. Fidler
//' @export
//[[Rcpp::export]]
CharacterVector rxDfdy(const RObject &obj){
  List modVar = rxModelVars(obj);
  CharacterVector ret = modVar["dfdy"];
  return ret;
}

//' Left handed Variables
//'
//' This returns the model calculated variables
//'
//' @inheritParams rxModelVars
//'
//' @return a character vector listing the calculated parameters
//' @seealso \code{\link{RxODE}}
//'
//' @author Matthew L.Fidler
//' @export
//[[Rcpp::export]]
CharacterVector rxLhs(const RObject &obj){
  List modVar = rxModelVars(obj);
  CharacterVector ret = modVar["lhs"];
  return ret;
}
NumericVector rxInits0(const RObject &obj,
		       Nullable<NumericVector> vec = R_NilValue,
		       Nullable<CharacterVector> req = R_NilValue,
		       double defaultValue = 0,
		       bool noerror = false,
		       bool noini=false){
  NumericVector oini;
  CharacterVector cini;
  List modVar = rxModelVars(obj);
  if (!noini){
    oini = (modVar["ini"]);
    cini = oini.names();
  }
  int i, j, k;
  CharacterVector nreq;
  NumericVector miss;
  if (!req.isNull()){
    nreq = CharacterVector(req);
    if ((ISNA(defaultValue) && noerror) || !ISNA(defaultValue)){
      miss = NumericVector(nreq.size());
      for (i = 0; i < nreq.size(); i++) {
	miss[i] = defaultValue;
      }
      miss.attr("names") = CharacterVector(nreq);
    }
  }
  NumericVector nvec;
  CharacterVector nvecNames;
  if (!vec.isNull()){
    nvec = NumericVector(vec);
    if (nvec.size() > 0){
      if (!nvec.hasAttribute("names")){
	if (!req.isNull() && nreq.size() == nvec.size()){
	  nvec.attr("names") = req;
	  nvecNames = req;
	  std::string wstr = "Assumed order of inputs: ";
	  for (i = 0; i < nreq.size(); i++){
	    wstr += (i == 0 ? "" : ", ") + nreq[i];
	  }
	  warning(wstr);
	} else {
	  std::string sstr = "Length mismatch\nreq: c(";
	  for (i = 0; i < nreq.size(); i++){
	    sstr += (i == 0 ? "" : ", ") + nreq[i];
	  }
	  sstr += ")\nvec: c(";
	  for (i = 0; i < nvec.size(); i++){
            sstr += (i == 0 ? "" : ", ") + std::to_string((double)(nvec[i]));
          }
	  sstr += ")";
	  stop(sstr);
	}
      } else {
	nvecNames = nvec.names();
      }
    }
  }
  // Prefer c(vec, ini, miss)
  NumericVector ret;
  CharacterVector nret;
  if (!req.isNull()){
    ret =  NumericVector(nreq.size());
    bool found = false;
    for (i = 0; i < nreq.size(); i++){
      found = false;
      for (j = 0; !found && j < nvec.size(); j++){
	if (nreq[i] == nvecNames[j]){
	  found =  true;
	  ret[i] = nvec[j];
	  break;
	}
      }
      for (j = 0; !found && j < cini.size(); j++){
	if (nreq[i] == cini[j]){
          found =  true;
          ret[i] = oini[j];
          break;
        }
      }
      if (!found)
	ret[i] = miss[i];
    }
    ret.attr("names")= nreq;
  } else {
    // In this case
    // vec <- c(vec, ini);
    // vec <- vec[!duplicated(names(vec))]
    CharacterVector dupnames(nvec.size()+oini.size()+miss.size());
    j = 0;
    for (i = 0; i < nvec.size(); i++){
      dupnames[j] = nvecNames[i];
      j++;
    }
    for (i = 0; i < oini.size(); i++){
      dupnames[j] = cini[i];
      j++;
    }
    LogicalVector dups = duplicated(dupnames);
    j = 0;
    for (i = 0; i < dups.size(); i++){
      if (!dups[i]) j++;
    }
    ret = NumericVector(j);
    CharacterVector retn(j);
    k = 0, j = 0;
    for (i = 0; i < nvec.size(); i++){
      if (!dups[j]){
	ret[k] = nvec[i];
	retn[k] = nvecNames[i];
	k++;
      }
      j++;
    }
    for (i = 0; i < oini.size(); i++){
      if (!dups[j]){
	ret[k] = oini[i];
        retn[k] = cini[i];
        k++;
      }
      j++;
    }
    ret.attr("names") = retn;
  }
  return ret;
}
//' Initial Values and State values for a RxODE object
//'
//' Returns the initial values of the rxDll object
//'
//' @param obj rxDll, RxODE, or named vector representing default
//'     initial arguments
//'
//' @param vec If supplied, named vector for the model.
//'
//' @param req Required names, and the required order for the ODE solver
//'
//' @param defaultValue a number or NA representing the default value for
//'     parameters missing in \code{vec}, but required in \code{req}.
//'
//' @param noerror is a boolean specifying if an error should be thrown
//'     for missing parameter values when \code{default} = \code{NA}
//'
//' @keywords internal
//' @author Matthew L.Fidler
//' @export
//[[Rcpp::export]]
SEXP rxInits(const RObject &obj,
	     RObject vec = R_NilValue,
	     Nullable<CharacterVector> req = R_NilValue,
	     double defaultValue = 0,
	     bool noerror = false,
	     bool noini=false,
	     bool rxLines=false){
  if (rxLines){
    if (rxIs(obj, "NULL")){
      CharacterVector ret = "";
      return ret;
    }
    NumericVector inits = rxInits(obj, vec, req, defaultValue, noerror,noini,false);
    CharacterVector nms = inits.names();
    List mv = rxModelVars(obj);
    CharacterVector state = mv["state"];
    std::string ret="";
    bool isState;
    for (unsigned int j=inits.size(); j--;){
      isState=false;
      for (unsigned int k=state.size(); k--;){
	if (nms[j] == state[k]){
	  isState=true;
	  break;
	}
      }
      ret += as<std::string>(nms[j]);
      if (isState) ret += "(0)";
      ret += "=" + std::to_string(inits[j]) + ";\n";
    }
    return wrap(ret);
  } else if (vec.isNULL()){
    return wrap(rxInits0(obj, R_NilValue, req, defaultValue, noerror,noini));
  } else if (rxIs(vec, "list")){
    List vecL = as<List>(vec);
    Function unlist("unlist", R_BaseNamespace);
    NumericVector vec2 = as<NumericVector>(unlist(vec));
    // if (!vec2.hasAttribute("names")){
    //   stop("When using a list for inits or scales, the list must be named. list(depot=1)");
    // }
    if (vec2.size() != vecL.size()){
      stop("Only one estimate per named list item; i.e. list(x=1) instead of list(x=1:2).");
    }
    return wrap(rxInits0(obj, vec2, req, defaultValue, noerror,noini));
  } else if (rxIs(vec, "integer") || rxIs(vec, "numeric")){
    return wrap(rxInits0(obj, as<NumericVector>(vec), req, defaultValue, noerror,noini));
  } else {
    stop("Incompatible initial estimate type.");
  }
}

//' Setup the initial conditions.
//'
//' @param obj RxODE object
//' @param inits A numeric vector of initial conditions.
//' @author Matthew L. Fidler
//' @keywords internal
//' @export
//[[Rcpp::export]]
NumericVector rxSetupIni(const RObject &obj,
			   RObject inits = R_NilValue){
  List modVars = rxModelVars(obj);
  CharacterVector state = modVars["state"];
  return rxInits(obj, inits, state, 0.0);
}

//' Setup the initial conditions.
//'
//' @param obj RxODE object
//' @param inits A numeric vector of initial conditions.
//' @param extraArgs A list of extra args to parse for initial conditions.
//' @author Matthew L. Fidler
//' @keywords internal
//' @export
//[[Rcpp::export]]
NumericVector rxSetupScale(const RObject &obj,
			   RObject scale = R_NilValue,
			   Nullable<List> extraArgs = R_NilValue){
  List modVars = rxModelVars(obj);
  CharacterVector state = modVars["state"];
  NumericVector ret = rxInits(obj, scale, state, 1.0, false, true);
  unsigned int usedS = 0, foundS = 0;
  if (!extraArgs.isNull()){
    List xtra = as<List>(extraArgs);
    int i, n=state.size();
    std::string cur;
    for (i = 0; i < n; i++){
      cur = "S" + std::to_string(i+1);
      if (xtra.containsElementNamed(cur.c_str())){
	if (ret[i] == 1.0){
	  ret[i] = as<double>(xtra[cur]);
	  usedS++;
	} else {
	  stop("Trying to scale the same compartment by scale=c(%s=%f,...) and S%d=%f;  Cannot do both.",
	       (as<std::string>(state[i])).c_str(), ret[i], i+1,as<double>(xtra[i]));
	}
      } else {
	cur = "s" + std::to_string(i+1);
        if (xtra.containsElementNamed(cur.c_str())){
          if (ret[i] == 1.0){
            ret[i] = as<double>(xtra[cur]);
	    usedS++;
          } else {
            stop("Trying to scale the same compartment by scale=c(%s=%f,...) and s%d=%f;  Cannot do both.",
                 (as<std::string>(state[i])).c_str(), ret[i], i+1,as<double>(xtra[i]));
          }
        }
      }
    }
    Nullable<StringVector> xtraN2 = xtra.names();
    if (!xtraN2.isNull()){
      StringVector xtraN = StringVector(xtraN2);
      n = xtraN.size();
      std::string tmpstr;
      for (i = 0; i < n; i++){
	tmpstr = (as<std::string>(xtraN[i])).substr(0,1);
	if (tmpstr == "S" || tmpstr == "s"){
	  tmpstr = (as<std::string>(xtraN[i])).substr(1,1);
	  if (tmpstr == "0" || tmpstr == "1" || tmpstr == "2" || tmpstr == "3" || tmpstr == "4" ||
	      tmpstr == "5" || tmpstr == "6" || tmpstr == "7" || tmpstr == "8" || tmpstr == "9"){
	    foundS++;
          }
        }
      }
    }
    if (foundS > usedS){
      warning("Scaled a compartment that is not defined by the RxODE model.");
    }
  }
  return ret;
}

NumericVector rxSetupParamsThetaEtaThetaN(const RObject &theta = R_NilValue, std::string thetaTxt = "theta"){
  NumericVector thetaN;
  if (rxIs(theta,"numeric") || rxIs(theta,"integer")){
    thetaN = as<NumericVector>(theta);
  } else if (rxIs(theta, "matrix")){
    NumericMatrix thetaM = as<NumericMatrix>(theta);
    if (thetaM.nrow() == 1){
      thetaN = NumericVector(thetaM.ncol());
      for (unsigned int i = thetaM.ncol() ; i--;){
        thetaN[i] = thetaM(1,i);
      }
    } else if (thetaM.ncol() == 1){
      thetaN = NumericVector(thetaM.nrow());
      for (unsigned int i = thetaM.ncol() ; i-- ;){
        thetaN[i] = thetaM(i, i);
      }
    } else {
      stop("'%s' is not compatible with params, check dimensions to make sure they are compatible.", thetaTxt.c_str());
    }
  } else if (!theta.isNULL()){
    stop("'%s' is not compatible with params, check dimensions to make sure they are compatible.", thetaTxt.c_str());
  }
  return thetaN;
}

NumericMatrix rxSetupParamsThetaEtaNullParams(const RObject &theta = R_NilValue, const RObject &eta = R_NilValue){
  // Create the matrix
  NumericVector thetaN = rxSetupParamsThetaEtaThetaN(theta,"theta");
  // Now eta
  NumericVector etaN = rxSetupParamsThetaEtaThetaN(eta, "eta");
  NumericMatrix tmp1(1, thetaN.size()+etaN.size());
  CharacterVector tmpN = CharacterVector(tmp1.size());
  unsigned int i;
  for (i = thetaN.size(); i--;){
    tmp1(0, i) = thetaN[i];
    tmpN[i] = "THETA[" + std::to_string(i + 1) + "]";
  }
  i = thetaN.size();
  for (; i < thetaN.size()+ etaN.size(); i++){
    tmp1(0, i) = etaN[i - thetaN.size()];
    tmpN[i] = "ETA[" + std::to_string(i - thetaN.size() + 1) + "]";
  }
  tmp1.attr("dimnames") = List::create(R_NilValue, tmpN);
  return tmp1;
}

RObject rxSetupParamsThetaEta(const RObject &params = R_NilValue,
                                    const RObject &theta = R_NilValue,
				    const RObject &eta = R_NilValue){
  // Now get the parameters as a data.frame
  if (rxIs(params,"list")) {
    Function asDf("as.data.frame", R_BaseNamespace);
    return rxSetupParamsThetaEta(asDf(params), theta, eta);
  }
  NumericMatrix parMat;
  int i;
  if (params.isNULL()){
    if (!theta.isNULL() || !eta.isNULL()){
      parMat = rxSetupParamsThetaEtaNullParams(theta, eta);
    }
  } else if (rxIs(params, "data.frame") || rxIs(params, "matrix")){
    if (rxIs(params,"data.frame")){
      DataFrame tmp = as<DataFrame>(params);
      int nr = tmp.nrows();
      NumericMatrix tmpM(nr,tmp.size());
      for (i = tmp.size(); i-- ;){
        tmpM(_,i) = NumericVector(tmp[i]);
      }
      tmpM.attr("dimnames") = List::create(R_NilValue,tmp.names());
      parMat=tmpM;
    } else {
      parMat = as<NumericMatrix>(params);
    }
  } else if (rxIs(params, "numeric") || rxIs(params, "integer")){
    // Create the matrix
    NumericVector thetaN;
    if (rxIs(theta,"numeric") || rxIs(theta,"integer")){
      thetaN = as<NumericVector>(theta);
    } else if (rxIs(theta, "matrix")){
      NumericMatrix thetaM = as<NumericMatrix>(theta);
      if (thetaM.nrow() == 1){
        thetaN = NumericVector(thetaM.ncol());
        for (i = thetaM.ncol() ; i-- ;){
          thetaN[i] = thetaM(1,i);
        }
      } else if (thetaM.ncol() == 1){
        thetaN = NumericVector(thetaM.nrow());
        for (i = thetaM.ncol(); i--;){
          thetaN[i] = thetaM(i, i);
        }
      } else {
        stop("'theta' is not compatible with params, check dimensions to make sure they are compatible.");
      }
    } else if (!theta.isNULL()){
      stop("'theta' is not compatible with params, check dimensions to make sure they are compatible.");
    }
    // Now eta
    NumericVector etaN;
    if (rxIs(eta,"numeric") || rxIs(eta,"integer")){
      etaN = as<NumericVector>(eta);
    } else if (rxIs(eta, "matrix")){
      NumericMatrix etaM = as<NumericMatrix>(eta);
      if (etaM.nrow() == 1){
        etaN = NumericVector(etaM.ncol());
        for (i = etaM.ncol() ; i-- ;){
          etaN[i] = etaM(0, i);
        }
      } else if (etaM.ncol() == 1){
        etaN = NumericVector(etaM.nrow());
        for (i = etaM.ncol() ; i-- ;){
          etaN[i] = etaM(i, 0);
        }
      } else {
        stop("'eta' is not compatible with params, check dimensions to make sure they are compatible.");
      }
    } else if (!eta.isNULL()){
      stop("'eta' is not compatible with params, check dimensions to make sure they are compatible.");
    }
    NumericVector tmp0 = as<NumericVector>(params);
    NumericMatrix tmp1(1, tmp0.size()+thetaN.size()+etaN.size());
    CharacterVector tmp0N ;
    NumericVector pars = as<NumericVector>(params);
    if (tmp0.hasAttribute("names")){
      tmp0N = tmp0.names();
    } else if (tmp0.size() == pars.size()){
      tmp0N = pars;
    } else if (tmp0.size() > 0){
      // In this case there are no names
      stop("The parameter names must be specified.");
    }
    CharacterVector tmpN = CharacterVector(tmp1.size());
    for (i = 0; i < tmp0.size(); i++){
      tmp1(0, i) = tmp0[i];
      tmpN[i]   = tmp0N[i];
    }
    for (; i < tmp0.size()+thetaN.size(); i++){
      tmp1(0, i) = thetaN[i - tmp0.size()];
      tmpN[i] = "THETA[" + std::to_string(i - tmp0.size() + 1) + "]";
    }
    for (; i < tmp0.size()+thetaN.size()+ etaN.size(); i++){
      tmp1(0, i) = etaN[i - tmp0.size() - thetaN.size()];
      tmpN[i] = "ETA[" + std::to_string(i - tmp0.size() - thetaN.size() + 1) + "]";
    }
    tmp1.attr("dimnames") = List::create(R_NilValue, tmpN);
    parMat = tmp1;
  }
  return as<RObject>(parMat);
}


typedef struct {
  double *gsolve;
  int gsolven;
  double *gInfusionRate;
  int gInfusionRaten;
  double *gall_times;
  int gall_timesn;
  double *gdv;
  int gdvn;
  double *gamt;
  int gamtn;
  double *glhs;
  int glhsn;
  double *gcov;
  int gcovn;
  double *ginits;
  int ginitsn;
  double *gscale;
  int gscalen;
  double *gatol2;
  int gatol2n;
  double *grtol2;
  int grtol2n;
  double *gpars;
  int gparsn;
  //ints
  int *gevid;
  int gevidn;
  int *gBadDose;
  int gBadDosen;
  int *grc;
  int grcn;
  int *gidose;
  int gidosen;
  int *gpar_cov;
  int gpar_covn;
  int *gParPos;
  int gParPosn;
  int *gsvar;
  int gsvarn;
  int *gsiV;
  int gsiVn;
  //
  int *slvr_counter;
  int slvr_countern;
  int *dadt_counter;
  int dadt_countern;
  int *jac_counter;
  int jac_countern;
} rx_globals;

rx_globals _globals;


extern "C" void rxOptionsIniData(){
  _globals.gsolve = NULL;//Calloc(NCMT*NALL,double);
  _globals.gsolven=0;//NCMT*NALL;
  _globals.gInfusionRate = NULL;//Calloc(NCMT,double);
  _globals.gInfusionRaten=0;//NCMT;
  _globals.gall_times = NULL;//Calloc(NALL,double);
  _globals.gall_timesn=0;//NALL;
  _globals.gdv = NULL;//Calloc(NALL,double);
  _globals.gdvn=0;//NALL;
  _globals.gamt = NULL;//Calloc(NDOSES,double);
  _globals.gamtn=0;//NDOSES;
  _globals.glhs = NULL;//Calloc(NPARS,double);
  _globals.glhsn=0;//NPARS;
  _globals.gcov = NULL;//Calloc(NALL*10,double);
  _globals.gcovn=0;//NALL*10;
  _globals.ginits = NULL;//Calloc(NCMT,double);
  _globals.ginitsn=0;//NCMT;
  _globals.gscale = NULL;//Calloc(NCMT,double);
  _globals.gscalen=0;//NCMT;
  _globals.gatol2 = NULL;//Calloc(NCMT,double);
  _globals.gatol2n=0;//NCMT;
  _globals.grtol2 = NULL;//Calloc(NCMT,double);
  _globals.grtol2n=0;//NCMT;
  _globals.gpars = NULL;//Calloc(NPARS,double);
  _globals.gparsn=0;//NPARS;
  //ints
  _globals.gevid = NULL;//Calloc(NALL, int);
  _globals.gevidn = 0;//NALL;
  _globals.gBadDose = NULL;//Calloc(NCMT, int);
  _globals.gBadDosen = 0;//NCMT;
  _globals.grc = NULL;//Calloc(MAXIDS, int);
  _globals.grcn = 0;//MAXIDS;
  _globals.gidose = NULL;//Calloc(NALL, int);
  _globals.gidosen = 0;//NALL;
  _globals.gpar_cov = NULL;//Calloc(NCMT, int);
  _globals.gpar_covn = 0;//NCMT;
  _globals.gParPos = NULL;//Calloc(NCMT, int);
  _globals.gParPosn = 0;//NCMT;
  _globals.gsvar = NULL;//Calloc(NPARS, int);
  _globals.gsvarn = 0;//NPARS;
  _globals.gsiV = NULL;//Calloc(NCMT, int);
  _globals.gsiVn = 0;//NCMT;
  _globals.slvr_counter = NULL;//Calloc(MAXIDS, int);
  _globals.slvr_countern = 0;//MAXIDS;
  _globals.dadt_counter = NULL;//Calloc(MAXIDS, int);
  _globals.dadt_countern = 0;//MAXIDS;
  _globals.jac_counter = NULL;//Calloc(MAXIDS, int);
  _globals.jac_countern = 0;//MAXIDS;
}

void gsolveSetup(int n){
  if (_globals.gsolven < 0){
    _globals.gsolven=0;
    _globals.gsolve=NULL;
  }
  if (_globals.gsolven < n){
    int cur = n;
    Free( _globals.gsolve);
    _globals.gsolve = Calloc(cur, double);
    _globals.gsolven=cur;
  }
}

void gInfusionRateSetup(int n){
  if (_globals.gInfusionRaten < n){
    int cur = n;
    Free(_globals.gInfusionRate);
    _globals.gInfusionRate = Calloc(cur, double);
    _globals.gInfusionRaten=cur;
  }
}

void gall_timesSetup(int n){
  if (_globals.gall_timesn < 0){
    _globals.gall_timesn=0;
    _globals.gall_times=NULL;
  }
  if (_globals.gall_timesn < n){
    int cur = n;
    Free(_globals.gall_times);
    _globals.gall_times = Calloc(cur, double);
    _globals.gall_timesn=cur;
  }
}

void gdvSetup(int n){
  if (_globals.gdvn < n){
    int cur = n;
    Free(_globals.gdv);
    _globals.gdv = Calloc(cur, double);
    _globals.gdvn=cur;
  }
}

void gamtSetup(int n){
  if (_globals.gamtn < 0){
    _globals.gamtn=0;
    _globals.gamt=NULL;
  }
  if (_globals.gamtn < n){
    int cur = n;
    Free(_globals.gamt);
    _globals.gamt = Calloc(cur, double);
    _globals.gamtn = cur;
  }
}

void glhsSetup(int n){
  if (_globals.glhsn < n){
    _globals.glhsn=0;
    _globals.glhs=NULL;
  }
  if (_globals.glhsn < n){
    int cur = n;
    Free(_globals.glhs);
    _globals.glhs = Calloc(cur, double);
    _globals.glhsn =cur;
  }
}

void gcovSetup(int n){
  if (_globals.gcovn < n){
    int cur = n;
    Free(_globals.gcov);
    _globals.gcov = Calloc(cur, double);
    _globals.gcovn = cur;
  }
}

void ginitsSetup(int n){
  if (_globals.ginitsn < 0){
    _globals.ginits = NULL;
    _globals.ginitsn = 0;
  }
  if (_globals.ginitsn < n){
    int cur = n;
    Free(_globals.ginits);
    _globals.ginits = Calloc(cur, double);
    _globals.ginitsn = cur;
  }
}

void gscaleSetup(int n){
  if (_globals.gscalen < n){
    int cur = n;
    Free(_globals.gscale);
    _globals.gscale = Calloc(cur, double);
    _globals.gscalen = cur;
  }
}

void gatol2Setup(int n){
  if (_globals.gatol2n < n){
    int cur = n;
    Free(_globals.gatol2);
    _globals.gatol2 = Calloc(cur, double);
    _globals.gatol2n = cur;
  }
}

void grtol2Setup(int n){
  if (_globals.grtol2n < n){
    int cur = n;
    Free(_globals.grtol2);
    _globals.grtol2 = Calloc(cur, double);
    _globals.grtol2n = cur;
  }
}


void gparsSetup(int n){
  if (_globals.gparsn < 0){
    _globals.gparsn=0;
    _globals.gpars=NULL;
  }
  if (_globals.gparsn < n){
    int cur = n;
    Free(_globals.gpars);
    _globals.gpars = Calloc(cur, double);
    cur = _globals.gparsn;
  }
}

void gevidSetup(int n){
  if (_globals.gevidn < 0){
    _globals.gevidn = 0;
    _globals.gevid = NULL;
  }
  if (_globals.gevidn < n){
    int cur = n;
    Free(_globals.gevid);
    _globals.gevid = Calloc(cur, int);
    _globals.gevidn  = cur;
  }
}

void gBadDoseSetup(int n){
  if (_globals.gBadDosen < n){
    int cur = n;
    Free(_globals.gBadDose);
    _globals.gBadDose = Calloc(cur, int);
    _globals.gBadDosen  = cur;
  }
}

void grcSetup(int n){
  if (_globals.grcn < 0){
    _globals.grcn=0;
    _globals.grc=NULL;
  }
  if (_globals.grcn < n){
    int cur = n;
    Free(_globals.grc);
    _globals.grc = Calloc(cur, int);
    _globals.grcn  = cur;
  }
}

extern "C" int *gslvr_counterSetup(int n){
  if (_globals.slvr_countern < n){
    int cur = n;
    Free(_globals.slvr_counter);
    _globals.slvr_counter = Calloc(cur , int);
    _globals.slvr_countern = cur;
  }
  return _globals.slvr_counter;
}

extern "C" int *gdadt_counterSetup(int n){
  if (_globals.dadt_countern < n){
    int cur = n;
    Free(_globals.dadt_counter);
    _globals.dadt_counter = Calloc(cur, int);
    _globals.dadt_countern = cur;
  }
  return _globals.dadt_counter;
}

extern "C" int *gjac_counterSetup(int n){
  if (_globals.jac_countern < n){
    int cur = n;
    Free(_globals.jac_counter);
    _globals.jac_counter = Calloc(cur, int);
    _globals.jac_countern = cur;
  }
  return _globals.jac_counter;
}

extern "C" int *gidoseSetup(int n){
  if (_globals.gidosen < n){
    int cur = n;
    Free(_globals.gidose);
    _globals.gidose = Calloc(cur, int);
    _globals.gidosen = cur;
  }
  return _globals.gidose;
}

void gpar_covSetup(int n){
  if (_globals.gpar_covn < n){
    int cur = n;
    Free(_globals.gpar_cov);
    _globals.gpar_cov = Calloc(cur, int);
    _globals.gpar_covn = cur;
  }
}

void gParPosSetup(int n){
  if (_globals.gParPosn < n){
    int cur = n;
    Free(_globals.gParPos);
    _globals.gParPos = Calloc(cur, int);
    _globals.gParPosn = cur;
  }
}

void gsvarSetup(int n){
  if (_globals.gsvarn < n){
    int cur = n;
    Free(_globals.gsvar);
    _globals.gsvar = Calloc(cur, int);
    _globals.gsvarn = cur;
  }
}

extern "C" void protectOld(){
  _globals.gparsn=-1;
  _globals.gamtn=-1;
  _globals.gsolven=-1;
  _globals.glhsn=-1;
  _globals.gevidn=-1;
  _globals.grcn=-1;
  _globals.gall_timesn=-1;
  _globals.ginitsn=-1;
}


extern "C" int *gsiVSetup(int n){
  if (_globals.gsiVn < n){
    int cur = n;
    Free(_globals.gsiV);
    _globals.gsiV = Calloc(cur, int);
    _globals.gsiVn = cur;
  }
  return _globals.gsiV;
}

extern "C" void gFree(){
  if (_globals.gsiV != NULL) Free(_globals.gsiV);
  _globals.gsiVn=0;
  if (_globals.gsvar != NULL) Free(_globals.gsvar);
  _globals.gsvarn=0;
  if (_globals.gpar_cov != NULL) Free(_globals.gpar_cov);
  _globals.gpar_covn=0;
  if (_globals.gidose != NULL) Free(_globals.gidose);
  _globals.gidosen=0;
  if (_globals.grc != NULL && _globals.grcn > 0) Free(_globals.grc);
  _globals.grc=NULL;
  _globals.grcn=0;
  if (_globals.gBadDose != NULL) Free(_globals.gBadDose);
  _globals.gBadDosen=0;
  if (_globals.gevid != NULL && _globals.gevidn > 0) Free(_globals.gevid);
  _globals.gevid=NULL;
  _globals.gevidn=0;
  if (_globals.gpars != NULL && _globals.gparsn>0) Free(_globals.gpars);
  _globals.gpars=NULL;
  _globals.gparsn=0;
  if (_globals.grtol2 != NULL) Free(_globals.grtol2);
  _globals.grtol2n=0;
  if (_globals.gatol2 != NULL) Free(_globals.gatol2);
  _globals.gatol2n=0;
  if (_globals.gscale != NULL) Free(_globals.gscale);
  _globals.gscalen=0;
  if (_globals.ginits != NULL && _globals.ginitsn > 0) Free(_globals.ginits);
  _globals.ginits=NULL;
  _globals.ginitsn=0;
  if (_globals.gcov != NULL) Free(_globals.gcov);
  _globals.gcovn=0;
  if (_globals.glhs != NULL && _globals.glhsn > 0) Free(_globals.glhs);
  _globals.glhs=NULL;
  _globals.glhsn=0;
  if (_globals.gamt != NULL && _globals.gamtn > 0) Free(_globals.gamt);
  _globals.gamt=NULL;
  _globals.gamtn=0;
  if (_globals.gall_times != NULL && _globals.gall_timesn>0) Free(_globals.gall_times);
  _globals.gall_times=NULL;
  _globals.gall_timesn=0;
  if (_globals.gdv != NULL) Free(_globals.gdv);
  _globals.gdvn=0;
  if (_globals.gInfusionRate != NULL) Free(_globals.gInfusionRate);
  _globals.gInfusionRaten=0;
  if (_globals.gsolve != NULL&& _globals.gsolven>0) Free(_globals.gsolve);
  _globals.gsolve=NULL;
  _globals.gsolven=0;
  if (_globals.gParPos != NULL) Free(_globals.gParPos);
  _globals.gParPosn = 0;
}

arma::mat rwish5(double nu, int p){
  GetRNGstate();
  arma::mat Z(p,p, fill::zeros);
  double curp = nu;
  double tmp =sqrt(Rf_rchisq(curp--));
  Z(0,0) = (tmp < 1e-100) ? 1e-100 : tmp;
  int i, j;
  if (p > 1){
    for (i = 1; i < (int)p; i++){
      tmp = sqrt(Rf_rchisq(curp--));
      Z(i,i) = (tmp < 1e-100) ? 1e-100 : tmp;
      for (j = 0; j < i; j++){
        // row,col
        Z(j,i) = norm_rand();
      }
    }
  }
  PutRNGstate();
  return Z;
}

NumericMatrix cvPost0(double nu, NumericMatrix omega, bool omegaIsChol = false,
                      bool returnChol = false){
  arma::mat S =as<arma::mat>(omega);
  int p = S.n_rows;
  if (p == 1){
    GetRNGstate();
    NumericMatrix ret(1,1);
    if (omegaIsChol){
      ret[0] = nu*omega[0]*omega[0]/(Rf_rgamma(nu/2.0,2.0));
    } else {
      ret[0] = nu*omega[0]/(Rf_rgamma(nu/2.0,2.0));
    }
    if (returnChol) ret[0] = sqrt(ret[0]);
    PutRNGstate();
    return ret;
  } else {
    arma::mat Z = rwish5(nu, p);
    // Backsolve isn't available in armadillo
    arma::mat Z2 = arma::trans(arma::inv(trimatu(Z)));
    arma::mat cv5;
    if (omegaIsChol){
      cv5 = S;
    } else {
      cv5 = arma::chol(S);
    }
    arma::mat mat1 = Z2 * cv5;
    mat1 = mat1.t() * mat1;
    mat1 = mat1 * nu;
    if (returnChol) mat1 = arma::chol(mat1);
    return wrap(mat1);
  }
}

//' Sample a covariance Matrix from the Posteior Inverse Wishart distribution.
//'
//' Note this Inverse wishart rescaled to match the original scale of the covariance matrix.
//'
//' If your covariance matrix is a 1x1 matrix, this uses an scaled inverse chi-squared which 
//' is equivalent to the Inverse Wishart distribution in the uni-directional case.
//'
//' @param nu Degrees of Freedom (Number of Observations) for 
//'        covariance matrix simulation.
//' @param omega Estimate of Covariance matrix.
//' @param n Number of Matricies to sample.  By default this is 1.
//' @param omegaIsChol is an indicator of if the omega matrix is in the cholesky decomposition. 
//' @param returnChol Return the cholesky decomposition of the covariance matrix sample.
//'
//' @return a matrix (n=1) or a list of matricies (n > 1)
//'
//' @author Matthew L.Fidler & Wenping Wang
//' 
//' @export
//[[Rcpp::export]]
RObject cvPost(double nu, RObject omega, int n = 1, bool omegaIsChol = false, bool returnChol = false){
  if (n == 1){
    if (rxIs(omega,"numeric.matrix") || rxIs(omega,"integer.matrix")){
      return as<RObject>(cvPost0(nu, as<NumericMatrix>(omega), omegaIsChol));
    } else if (rxIs(omega, "numeric") || rxIs(omega, "integer")){
      NumericVector om1 = as<NumericVector>(omega);
      if (om1.size() % 2 == 0){
        int n1 = om1.size()/2;
        NumericMatrix om2(n1,n1);
        for (int i = 0; i < om1.size();i++){
          om2[i] = om1[i];
        }
        return as<RObject>(cvPost0(nu, om2, omegaIsChol, returnChol));
      }
    }
  } else {
    List ret(n);
    for (int i = 0; i < n; i++){
      ret[i] = cvPost(nu, omega, 1, omegaIsChol, returnChol);
    }
    return(as<RObject>(ret));
  }
  stop("omega needs to be a matrix or a numberic vector that can be converted to a matrix.");
  return R_NilValue;
}

//' Scaled Inverse Chi Squared distribution
//'
//' @param n Number of random samples
//' @param nu degrees of freedom of inverse chi square
//' @param scale  Scale of inverse chi squared distribution 
//'         (default is 1).
//' @return a vector of inverse chi squared deviates .
//' @export
//[[Rcpp::export]]
NumericVector rinvchisq(const int n = 1, const double &nu = 1.0, const double &scale = 1){
  NumericVector ret(n);
  GetRNGstate();
  for (int i = 0; i < n; i++){
    ret[i] = nu*scale/(Rf_rgamma(nu/2.0,2.0));
  }
  PutRNGstate();
  return ret;
}

extern "C" double *rxGetErrs(){
  getRxModels();
  if (_rxModels.exists(".sigma")){
    NumericMatrix sigma = _rxModels[".sigma"];
    return &sigma[0];
  }
  return NULL;
}

extern "C" int rxGetErrsNcol(){
  getRxModels();
  if (_rxModels.exists(".sigma")){
    NumericMatrix sigma = _rxModels[".sigma"];
    int ret = sigma.ncol();
    return ret;
  } 
  return 0;
}
  
SEXP rxGetFromChar(char *ptr, std::string var){
  std::string str(ptr);
  // Rcout << str << "\n";
  CharacterVector cv(1);
  cv[0] = str;
  List mv = rxModelVars(as<RObject>(cv));
  if (var == ""){
    return wrap(mv);
  } else {
    return wrap(mv[var]);
  }
}

//' Simulate Parameters from a Theta/Omega specification
//'
//' @param params Named Vector of RxODE model parameters
//'
//' @param thetaMat Named theta matrix.
//'
//' @param thetaDf The degrees of freedom of a t-distribution for
//'     simulation.  By default this is \code{NULL} which is
//'     equivalent to \code{Inf} degrees, or to simulate from a normal
//'     distribution instead of a t-distribution.
//'
//' @param thetaIsChol Indicates if the \code{theta} supplied is a
//'     Cholesky decomposed matrix instead of the traditional
//'     symmetric matrix.
//'
//' @param nSub Number between subject variabilities (ETAs) simulated for every 
//'        realization of the parameters.
//'
//' @param omega Named omega matrix.
//'
//' @param omegaDf The degrees of freedom of a t-distribution for
//'     simulation.  By default this is \code{NULL} which is
//'     equivalent to \code{Inf} degrees, or to simulate from a normal
//'     distribution instead of a t-distribution.
//'
//' @param omegaIsChol Indicates if the \code{omega} supplied is a
//'     Cholesky decomposed matrix instead of the traditional
//'     symmetric matrix.
//'
//' @param nStud Number virtual studies to characterize uncertainty in estimated 
//'        parameters.
//'
//' @param nObs Number of observations to simulate (with \code{sigma} matrix)
//'
//' @param sigma Matrix for residual variation.  Adds a "NA" value for each of the 
//'     indivdual parameters, residuals are updated after solve is completed. 
//'
//' @inheritParams rxSolve
//'
//' @param dfSub Degrees of freedom to sample the between subject variaiblity matrix from the 
//'        inverse Wishart distribution (scaled) or scaled inverse chi squared distribution.
//'
//' @param dfObs Degrees of freedom to sample the unexplained variaiblity matrix from the 
//'        inverse Wishart distribution (scaled) or scaled inverse chi squared distribution. 
//'
//' @param simSubjects boolean indicated RxODE should simulate subjects in studies (\code{TRUE}, 
//'         default) or studies (\code{FALSE})
//'
//' @author Matthew L.Fidler
//'
//' @export
//[[Rcpp::export]]
List rxSimThetaOmega(const Nullable<NumericVector> &params    = R_NilValue,
                     const Nullable<NumericMatrix> &omega= R_NilValue,
                     const Nullable<NumericVector> &omegaDf= R_NilValue,
                     const bool &omegaIsChol = false,
                     int nSub = 1,
                     const Nullable<NumericMatrix> &thetaMat = R_NilValue,
                     const Nullable<NumericVector> &thetaDf  = R_NilValue,
                     const bool &thetaIsChol = false,
                     int nStud = 1,
                     const Nullable<NumericMatrix> sigma = R_NilValue,
                     const Nullable<NumericVector> &sigmaDf= R_NilValue,
                     const bool &sigmaIsChol = false,
                     int nCoresRV = 1,
                     int nObs = 1,
                     double dfSub = 0,
                     double dfObs = 0,
		     bool simSubjects=true){
  NumericVector par;
  if (params.isNull()){
    stop("This function requires overall parameters.");
  } else {
    par = NumericVector(params);
    if (!par.hasAttribute("names")){
      stop("'params' must be a named vector.");
    }
  }
  bool simSigma = false;
  NumericMatrix sigmaM;
  CharacterVector sigmaN;
  NumericMatrix sigmaMC;
  if (!sigma.isNull() && nObs > 1){
    simSigma = true;
    sigmaM = as<NumericMatrix>(sigma);
    if (!sigmaM.hasAttribute("dimnames")){
      stop("'sigma' must be a named Matrix.");
    }
    if (sigmaIsChol){
      sigmaMC = sigmaM;
    } else {
      sigmaMC = wrap(arma::chol(as<arma::mat>(sigmaM)));
    }
    sigmaN = as<CharacterVector>((as<List>(sigmaM.attr("dimnames")))[1]);
  }
  int scol = 0;
  if (simSigma){
    scol = sigmaMC.ncol();
    if (simSubjects){
      if (nObs*nStud*nSub*scol < 0){
        // nStud = INT_MAX/(nObs*nSub*scol)*0.25;
        stop("Simulation Overflow; Reduce the number of observations, number of subjects or number of studies.");
      }
    } else {
      if (nObs*nStud*scol < 0){
        // nStud = INT_MAX/(nObs*nSub*scol)*0.25;
        stop("Simulation Overflow; Reduce the number of observations or number of studies.");
      }
    }
  }
  NumericMatrix thetaM;
  CharacterVector thetaN;
  bool simTheta = false;
  CharacterVector parN = CharacterVector(par.attr("names"));
  IntegerVector thetaPar(parN.size());
  int i, j, k;
  if (!thetaMat.isNull() && nStud > 1){
    thetaM = as<NumericMatrix>(thetaMat);
    if (!thetaM.hasAttribute("dimnames")){
      stop("'thetaMat' must be a named Matrix.");
    }
    thetaM = as<NumericMatrix>(rxSimSigma(as<RObject>(thetaMat), as<RObject>(thetaDf), nCoresRV, thetaIsChol, nStud));
    thetaN = as<CharacterVector>((as<List>(thetaM.attr("dimnames")))[1]);
    for (i = 0; i < parN.size(); i++){
      thetaPar[i] = -1;
      for (j = 0; j < thetaN.size(); j++){
        if (parN[i] == thetaN[j]){
          thetaPar[i] = j;
          break;
        }
      }
    }
    simTheta = true;
  } else if (!thetaMat.isNull() && nStud <= 1){
    warning("'thetaMat' is ignored since nStud <= 1.");
  }
  bool simOmega = false;
  NumericMatrix omegaM;
  CharacterVector omegaN;
  NumericMatrix omegaMC;
  if (!omega.isNull() && nSub > 1){
    simOmega = true;
    omegaM = as<NumericMatrix>(omega);
    if (!omegaM.hasAttribute("dimnames")){
      stop("'omega' must be a named Matrix.");
    }
    if (omegaIsChol){
      omegaMC = omegaM;
    } else {
      omegaMC = wrap(arma::chol(as<arma::mat>(omegaM)));
    }
    omegaN = as<CharacterVector>((as<List>(omegaM.attr("dimnames")))[1]);
  } else if (nSub > 1){
    stop("'omega' is required for multi-subject simulations.");
  }
  // Now create data frame of parameter values
  List omegaList;
  List sigmaList;  
  if (nStud > 1){
    if (dfSub > 0 && simOmega) {
      omegaList = cvPost(dfSub, as<RObject>(omegaMC), nStud,  true, false);
    }
    if (dfObs > 0 && simSigma){
      sigmaList = cvPost(dfObs, as<RObject>(sigmaMC), nStud,  true, false);
    }
  }
  int pcol = par.size();
  int ocol = 0;
  int ncol = pcol;
  if (simOmega){
    ocol = omegaMC.ncol();
    ncol += ocol;
  }
  NumericMatrix ret1;
  if (simSigma){
    ncol += scol;
    if (simSubjects){
      ret1 = NumericMatrix(nObs*nStud*nSub, scol);
    } else {
      ret1 = NumericMatrix(nObs*nStud, scol);
    }
  }
  List ret0(ncol);
  NumericVector nm;
  NumericMatrix nm1;
  for (i = 0; i < ncol; i++){
    nm = NumericVector(nSub*nStud);
    ret0[i] = nm;
  }
  for (i = 0; i < nStud; i++){
    for (j = 0; j < pcol; j++){
      nm = ret0[j];
      for (k = 0; k < nSub; k++){
        nm[nSub*i + k] = par[j];
      }
      if (simTheta){
        if(thetaPar[j] != -1){
          for (k = 0; k < nSub; k++){
            nm[nSub*i + k] += thetaM(i, thetaPar[j]);
          }
        }
      }
      ret0[j] = nm;
    }
    // Now Omega Covariates
    if (ocol > 0){
      if (dfSub > 0 && nStud > 1){
        // nm = ret0[j]; // parameter column
        nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(omegaList[i]), as<RObject>(omegaDf), nCoresRV, false, nSub,false));
      } else {
        nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(omegaMC), as<RObject>(omegaDf), nCoresRV, true, nSub,false));
      }
      for (j=pcol; j < pcol+ocol; j++){
        nm = ret0[j];
        for (k = 0; k < nSub; k++){
          nm[nSub*i + k] = nm1(k, j-pcol);
        }
        ret0[j] = nm;
      }
    }
    if (scol > 0){
      if (simSubjects){
        if (dfObs > 0  && nStud > 1){
          nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(sigmaList[i]), as<RObject>(sigmaDf), nCoresRV, false, nObs*nSub, false));
        } else {
          nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(sigmaMC), as<RObject>(sigmaDf), nCoresRV, true, nObs*nSub, false));
        }
        for (j = 0; j < scol; j++){
          for (k = 0; k < nObs*nSub; k++){
            // ret1 = NumericMatrix(nObs*nStud, scol);
            ret1(nObs*nSub*i+k, j) = nm1(k, j);
          }
        }
      } else {
        if (dfObs > 0  && nStud > 1){
          nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(sigmaList[i]), as<RObject>(sigmaDf), nCoresRV, false, nObs, false));
        } else {
          nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(sigmaMC), as<RObject>(sigmaDf), nCoresRV, true, nObs, false));
        }
        for (j = 0; j < scol; j++){
          for (k = 0; k < nObs; k++){
            // ret1 = NumericMatrix(nObs*nStud, scol);
            ret1(nObs*i+k, j) = nm1(k, j);
          }
        }
      }
    }
  }
  CharacterVector dfName(ncol);
  for (i = 0; i < pcol; i++){
    dfName[i] = parN[i];
  }
  for (i = pcol; i < pcol+ocol; i++){
    dfName[i] = omegaN[i-pcol];
  }
  for (i = pcol+ocol; i < ncol; i++){
    dfName[i] = sigmaN[i-pcol-ocol];
  }
  ret0.attr("names") = dfName;
  ret0.attr("class") = "data.frame";
  ret0.attr("row.names") = IntegerVector::create(NA_INTEGER,-nSub*nStud);
  getRxModels();
  if (ret1.nrow() > 1){
    ret1.attr("dimnames") = List::create(R_NilValue, sigmaN);
    _rxModels[".sigma"] = ret1;
  }
  if (simTheta){
    _rxModels[".theta"] = thetaM;
  }
  if (dfSub > 0 && nStud > 1){
    _rxModels[".omegaL"] = omegaList;
  }
  if (dfObs > 0 && nStud > 1){
    _rxModels[".sigmaL"] =sigmaList;
  }
  return ret0;
}


extern "C" double *global_InfusionRate(unsigned int mx);

#define defrx_params R_NilValue
#define defrx_events R_NilValue
#define defrx_inits R_NilValue
#define defrx_covs R_NilValue
#define defrx_method 2
#define defrx_transit_abs R_NilValue
#define defrx_atol 1.0e-8
#define defrx_rtol 1.0e-8
#define defrx_maxsteps 5000
#define defrx_hmin 0
#define defrx_hmax R_NilValue
#define defrx_hini 0
#define defrx_maxordn 12
#define defrx_maxords 5
#define defrx_cores 1
#define defrx_covs_interpolation 0
#define defrx_addCov false
#define defrx_matrix 0
#define defrx_sigma  R_NilValue
#define defrx_sigmaDf R_NilValue
#define defrx_nCoresRV 1
#define defrx_sigmaIsChol false
#define defrx_nDisplayProgress 10000
#define defrx_amountUnits NA_STRING
#define defrx_timeUnits "hours"
#define defrx_addDosing false
#define defrx_stateTrim R_PosInf


RObject rxCurObj;

Nullable<Environment> rxRxODEenv(RObject obj);

std::string rxDll(RObject obj);

bool rxDynLoad(RObject obj);

void updateSolveEnvPost(Environment e){
  if (!e.exists("params.dat")){
    List mv = rxModelVars(as<RObject>(e));
    NumericVector mvIni = mv["ini"];
    CharacterVector pars = mv["params"];
    RObject parso = e["args.params"];
    IntegerVector ppos = e[".par.pos"];
    if (rxIs(parso, "numeric") || rxIs(parso, "integer") ||
	rxIs(parso, "NULL")){
      double *tmp=Calloc(ppos.size(),double);
      // NumericVector   prs(ppos.size()-nrm);
      // CharacterVector prsn(ppos.size()-nrm+1);
      NumericVector parNumeric;
      if (!rxIs(parso, "NULL")){
        parNumeric= as<NumericVector>(parso);
      }
      unsigned int i, j=0;
      for (i = 0; i < ppos.size(); i++){
	if (ppos[i] > 0){ // User specified parameter
          tmp[j] = parNumeric[ppos[i]-1];
          // prsn[j] = pars[i];
	  j++;
        } else if (ppos[i] < 0) { // ini specified parameter.
          tmp[j] = mvIni[-ppos[i]-1];
          // prsn[j] = pars[i];
	  j++;
        }
      }
      NumericVector prs(j);
      CharacterVector prsn(j);
      std::copy(&tmp[0],&tmp[0]+j,prs.begin());
      Free(tmp);
      j=0;
      for (i = 0; i < ppos.size(); i++){
        if (ppos[i] != 0){ // User specified parameter
          prsn[j] = pars[i];
          j++;
        }
      }
      prs.names() = prsn;
      e["params.single"] = prs;
      List pd(prs.size());
      for (unsigned int j = prs.size();j--;){
	pd[j] = NumericVector::create(prs[j]);
      }
      pd.names() = prsn;
      pd.attr("class") = "data.frame";
      pd.attr("row.names") = IntegerVector::create(NA_INTEGER,-1);
      e["params.dat"] = pd;
      e["counts"] = DataFrame::create(_["slvr"]=e[".slvr.counter"],
				      _["dadt"]=e[".dadt.counter"],
				      _["jac"]=e[".jac.counter"]);
    } else {
      DataFrame parsdf = as<DataFrame>(parso);
      int extran = 0;
      int nsub = e[".nsub"], nsim=e[".nsim"];
      int nrm=0;
      for (unsigned int i = ppos.size(); i--;){
	if (ppos[i] == 0){ // Covariate or simulated variable.
	  nrm++;
	}
      }
      if (nsub > 1) extran++;
      if (nsim > 1) extran++;
      CharacterVector prsn(ppos.size()-nrm+extran);
      List prsl(ppos.size()-nrm+extran);
      unsigned int i, j=0;
      if (nsim > 1) {
	IntegerVector tmp(parsdf.nrow());
	for (unsigned int k = parsdf.nrow(); k--;){
	  tmp[k] = (k / nsub)+1;
	}
        prsn[j]="sim.id";
	prsl[j] = tmp;
	j++;
      }
      if (nsub > 1) {
        IntegerVector tmp(parsdf.nrow());
	for (unsigned int k = parsdf.nrow(); k--;){
	  tmp[k] = (k % nsub)+1;
	}
	prsl[j] = tmp;
	prsn[j]="id";
	j++;
      }
      for (i = 0; i < ppos.size();i++){
        if (ppos[i] > 0){ // User specified parameter
          prsl[j] = parsdf[ppos[i]-1];
          prsn[j] = pars[i];
          j++;
        } else if (ppos[i] < 0) { // ini specified parameter.
	  NumericVector tmp(parsdf.nrow(), mvIni[-ppos[i]-1]);
	  prsl[j] = tmp;
          prsn[j] = pars[i];
          j++;
        }
      }
      prsl.names() = prsn;
      prsl.attr("class") = "data.frame";
      prsl.attr("row.names") = IntegerVector::create(NA_INTEGER,-parsdf.nrow());
      e["params.dat"] = prsl;
      if (parsdf.nrow() == 1){
	NumericVector prsnv(prsl.size());
	for (j = prsl.size(); j--;){
	  prsnv[j] = (as<NumericVector>(prsl[j]))[0];
	}
	prsnv.names() = prsn;
        e["params.single"] = prsnv;
      } else {
	e["params.single"] = R_NilValue;
      }
      List cnt(extran+3);
      CharacterVector cntn(extran+3);
      j = 0;
      if (nsim > 1) {
	cntn[j]="sim.id";
        cnt[j] = prsl["sim.id"];
	j++;
      }
      if (nsub > 1) {
        cntn[j]="id";
        cnt[j] = prsl["id"];
	j++;
      }
      cntn[j] = "slvr";
      cnt[j++] = e[".slvr.counter"];
      cntn[j] = "dadt";
      cnt[j++] = e[".dadt.counter"];
      cntn[j] = "jac";
      cnt[j++] = e[".jac.counter"];
      cnt.names() = cntn;
      cnt.attr("class") = "data.frame";
      cnt.attr("row.names") = IntegerVector::create(NA_INTEGER,-parsdf.nrow());
      e["counts"] = cnt;
    }
  }
  if (!e.exists("EventTable")){
    RObject eventso = e["args.events"];
    if (rxIs(eventso, "EventTable")){
      List et = as<List>(eventso);
      Function f = et["get.EventTable"];
      e["EventTable"] = f();
      f = et["get.obs.rec"];
      e["obs.rec"] = f();
      f = et["get.dosing"];
      e["dosing"] = f();
      f = et["get.sampling"];
      e["sampling"] = f();
      f = et["get.units"];
      e["units"] = f();
    }
    e["covs"] = as<RObject>(e["args.covs"]);
    Function parse2("parse", R_BaseNamespace);
    Function eval2("eval", R_BaseNamespace);
    // eventTable style methods
    e["get.EventTable"] = eval2(_["expr"]   = parse2(_["text"]="function() EventTable"),
				_["envir"]  = e);
    e["get.obs.rec"] = eval2(_["expr"]   = parse2(_["text"]="function() obs.rec"),
			     _["envir"]  = e);
    e["get.nobs"] = eval2(_["expr"]   = parse2(_["text"]="function() nobs"),
			  _["envir"]  = e);
    e["add.dosing"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {.et <- create.eventTable(); .et$add.dosing(...); invisible(rxSolve(args.object,events=.et,updateObject=TRUE))}"),
			    _["envir"]  = e);
    e["clear.dosing"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {.et <- create.eventTable(); .et$clear.dosing(...); invisible(rxSolve(args.object,events=.et,updateObject=TRUE))}"),
			      _["envir"]  = e);
    e["get.dosing"] = eval2(_["expr"]   = parse2(_["text"]="function() dosing"),
			    _["envir"]  = e);

    e["add.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {.et <- create.eventTable(); .et$add.sampling(...); invisible(rxSolve(args.object,events=.et,updateObject=TRUE))}"),
			      _["envir"]  = e);
      
    e["clear.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {.et <- create.eventTable(); .et$clear.sampling(...); invisible(rxSolve(args.object,events=.et,updateObject=TRUE))}"),
				_["envir"]  = e);

    e["replace.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {.et <- create.eventTable(); .et$clear.sampling(); .et$add.sampling(...); invisible(rxSolve(args.object,events=.et,updateObject=TRUE))}"),
				  _["envir"]  = e);

    e["get.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function() sampling"),
			      _["envir"]  = e);
      
    e["get.units"] = eval2(_["expr"]   = parse2(_["text"]="function() units"),
			   _["envir"]  = e);

    e["import.EventTable"] = eval2(_["expr"]   = parse2(_["text"]="function(imp) {.et <- create.eventTable(imp); invisible(rxSolve(args.object,events=.et,updateObject=TRUE))}"),
				   _["envir"]  = e);
      
    e["create.eventTable"] = eval2(_["expr"]   = parse2(_["text"]="function(new.event) {.et <- eventTable(amount.units=units[1],time.units=units[2]);if (missing(new.event)) {.nev <- EventTable; } else {.nev <- new.event;}; .et$import.EventTable(.nev); return(.et);}"),
				   _["envir"]  = e);
    // Note event.copy doesn't really make sense...?  The create.eventTable does basically the same thing.
  }
}

extern "C" void rxOptionsFree();
extern "C" void rxOptionsIni();
extern "C" rx_solving_options_ind *rxOptionsIniEnsure(int mx);
extern "C" void rxOptionsFreeFocei();

//' Free the C solving information.
//'
//' Take the ODE C system and free it.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
LogicalVector rxSolveFree(){
  gFree();
  rxOptionsFree();
  rxOptionsIni();
  rxOptionsIniData();
  rxOptionsFreeFocei();
  return LogicalVector::create(true);
}
extern "C" void RxODE_assign_fn_pointers(SEXP);
SEXP rxSolveC(const RObject &obj,
              const Nullable<CharacterVector> &specParams = R_NilValue,
              const Nullable<List> &extraArgs = R_NilValue,
              const RObject &params = R_NilValue,
              const RObject &events = R_NilValue,
              const RObject &inits = R_NilValue,
              const RObject &scale = R_NilValue,
              const RObject &covs  = R_NilValue,
              const int method = 2, // 0
	      const Nullable<LogicalVector> &transit_abs = R_NilValue, //1
	      const double atol = 1.0e-6, //2
	      const double rtol = 1.0e-4, //3
	      const int maxsteps = 5000, //4
	      const double hmin = 0, //5
	      const Nullable<NumericVector> &hmax = R_NilValue, //6
	      const double hini = 0, //7
	      const int maxordn = 12, //8
	      const int maxords = 5, //9
	      const unsigned int cores = 1, //10
	      const int covs_interpolation = 0, //11
	      bool addCov = false, //12
	      int matrix = 0, //13
	      const Nullable<NumericMatrix> &sigma= R_NilValue, //14
	      const Nullable<NumericVector> &sigmaDf= R_NilValue, //15
	      const int &nCoresRV= 1, //16
	      const bool &sigmaIsChol= false,
	      const int &nDisplayProgress = 10000,
	      const CharacterVector &amountUnits = NA_STRING,
	      const CharacterVector &timeUnits = "hours",
	      const bool addDosing = false,
	      const double stateTrim = R_PosInf,
	      const RObject &theta = R_NilValue,
	      const RObject &eta = R_NilValue,
	      const bool updateObject = false,
	      const bool doSolve = true,
              const Nullable<NumericMatrix> &omega = R_NilValue, 
	      const Nullable<NumericVector> &omegaDf = R_NilValue, 
	      const bool &omegaIsChol = false,
              const unsigned int nSub = 1, 
	      const Nullable<NumericMatrix> &thetaMat = R_NilValue, 
	      const Nullable<NumericVector> &thetaDf = R_NilValue, 
	      const bool &thetaIsChol = false,
              const unsigned int nStud = 1, 
	      const double dfSub=0.0,
	      const double dfObs=0.0,
	      const int setupOnly = 0){
  RObject object;
  bool isRxSolve = rxIs(obj, "rxSolve");
  bool isEnvironment = rxIs(obj, "environment");
  if (updateObject && !isRxSolve && !isEnvironment){
    if (rxIs(rxCurObj, "rxSolve")){
      object = rxCurObj;
      isRxSolve = true;
    } else {
      stop("Cannot update this object.");
    }
  } else {
    object =obj;
  }
  if (isRxSolve || isEnvironment){
    bool update_params = false,
      update_events = false,
      update_inits = false,
      update_covs = false,
      update_method = false,
      update_transit_abs = false,
      update_atol = false,
      update_rtol = false,
      update_maxsteps = false,
      update_hini = false,
      update_hmin = false,
      update_hmax = false,
      update_maxordn = false,
      update_maxords = false,
      update_cores = false,
      update_covs_interpolation = false,
      update_addCov = false,
      update_matrix = false,
      update_sigma  = false,
      update_sigmaDf = false,
      update_nCoresRV = false,
      update_sigmaIsChol = false,
      update_amountUnits = false,
      update_timeUnits = false,
      update_scale = false,
      update_dosing = false,
      update_trim=false;
    if (specParams.isNull()){
      warning("No additional parameters were specified; Returning fit.");
      return object;
    }
    CharacterVector specs = CharacterVector(specParams);
    int n = specs.size(), i;
    for (i = n; i--;){
      if (as<std::string>(specs[i]) == "params")
	update_params = true;
      else if (as<std::string>(specs[i]) == "events")
	update_events = true;
      else if (as<std::string>(specs[i]) == "inits")
	update_inits = true;
      else if (as<std::string>(specs[i]) == "covs")
	update_covs = true;
      else if (as<std::string>(specs[i]) == "method")
	update_method = true;
      else if (as<std::string>(specs[i]) == "transit_abs")
	update_transit_abs = true;
      else if (as<std::string>(specs[i]) == "atol")
	update_atol = true;
      else if (as<std::string>(specs[i]) == "rtol")
	update_rtol = true;
      else if (as<std::string>(specs[i]) == "maxsteps")
	update_maxsteps = true;
      else if (as<std::string>(specs[i]) == "hmin")
	update_hmin = true;
      else if (as<std::string>(specs[i]) == "hmax")
	update_hmax = true;
      else if (as<std::string>(specs[i]) == "maxordn")
	update_maxordn = true;
      else if (as<std::string>(specs[i]) == "maxords")
	update_maxords = true;
      else if (as<std::string>(specs[i]) == "cores")
	update_cores = true;
      else if (as<std::string>(specs[i]) == "covs_interpolation")
	update_covs_interpolation = true;
      else if (as<std::string>(specs[i]) == "addCov")
	update_addCov = true;
      else if (as<std::string>(specs[i]) == "matrix")
	update_matrix = true;
      else if (as<std::string>(specs[i]) == "sigma")
	update_sigma  = true;
      else if (as<std::string>(specs[i]) == "sigmaDf")
	update_sigmaDf = true;
      else if (as<std::string>(specs[i]) == "nCoresRV")
	update_nCoresRV = true;
      else if (as<std::string>(specs[i]) == "sigmaIsChol")
	update_sigmaIsChol = true;
      else if (as<std::string>(specs[i]) == "amountUnits")
	update_amountUnits = true;
      else if (as<std::string>(specs[i]) == "timeUnits")
	update_timeUnits = true;
      else if (as<std::string>(specs[i]) == "hini")
	update_hini = true;
      else if (as<std::string>(specs[i]) == "scale")
	update_scale = true;
      else if (as<std::string>(specs[i]) == "addDosing")
	update_dosing = true;
      else if (as<std::string>(specs[i]) == "stateTrim")
	update_trim = true;
    }
    // Now update
    Environment e;
    List obj;
    if (isRxSolve){
      obj = as<List>(obj);
      CharacterVector classattr = object.attr("class");
      e = as<Environment>(classattr.attr(".RxODE.env"));
    } else  { // if (rxIs(object, "environment")) 
      e = as<Environment>(object);
      obj = as<List>(e["obj"]);
    }
    getRxModels();
    // if (e.exists("params.dat")){
    //   e.remove("params.dat");
    // }
    // if (e.exists("EventTable")){
    //   e.remove("EventTable");
    // }
    if(e.exists(".sigma")){
      _rxModels[".sigma"]=as<NumericMatrix>(e[".sigma"]);
    }
    if(e.exists(".sigmaL")){
      _rxModels[".sigmaL"]=as<List>(e[".sigmaL"]);
    }
    if(e.exists(".omegaL")){
      _rxModels[".omegaL"] = as<List>(e[".omegaL"]);
    }
    if(e.exists(".theta")){
      _rxModels[".theta"] = as<NumericMatrix>(e[".theta"]);
    }
    RObject new_params = update_params ? params : e["args.params"];
    RObject new_events = update_events ? events : e["args.events"];
    RObject new_inits = update_inits ? inits : e["args.inits"];
    RObject new_covs  = update_covs  ? covs  : e["args.covs"];
    int new_method = update_method ? method : e["args.method"];
    Nullable<LogicalVector> new_transit_abs = update_transit_abs ? transit_abs : e["args.transit_abs"];
    double new_atol = update_atol ? atol : e["args.atol"];
    double new_rtol = update_rtol ? rtol : e["args.rtol"];
    int new_maxsteps = update_maxsteps ? maxsteps : e["args.maxsteps"];
    int new_hmin = update_hmin ? hmin : e["args.hmin"];
    Nullable<NumericVector> new_hmax = update_hmax ? hmax : e["args.hmax"];
    int new_hini = update_hini ? hini : e["args.hini"];
    int new_maxordn = update_maxordn ? maxordn : e["args.maxordn"];
    int new_maxords = update_maxords ? maxords : e["args.maxords"];
    int new_cores = update_cores ? cores : e["args.cores"];
    int new_covs_interpolation = update_covs_interpolation ? covs_interpolation : e["args.covs_interpolation"];
    bool new_addCov = update_addCov ? addCov : e["args.addCov"];
    int new_matrix = update_matrix ? matrix : e["args.matrix"];
    Nullable<NumericMatrix> new_sigma = update_sigma ? sigma : e["args.sigma"];
    Nullable<NumericVector> new_sigmaDf = update_sigmaDf ? sigmaDf : e["args.sigmaDf"];
    int new_nCoresRV = update_nCoresRV ? nCoresRV : e["args.nCoresRV"];
    bool new_sigmaIsChol = update_sigmaIsChol ? sigmaIsChol : e["args.sigmaIsChol"];
    int new_nDisplayProgress = e["args.nDisplayProgress"];
    CharacterVector new_amountUnits = update_amountUnits ? amountUnits : e["args.amountUnits"];
    CharacterVector new_timeUnits = update_timeUnits ? timeUnits : e["args.timeUnits"];
    RObject new_scale = update_scale ? scale : e["args.scale"];
    bool new_addDosing = update_dosing ? addDosing : e["args.addDosing"];
    double new_stateTrim = update_trim ? stateTrim : e["args.stateTrim"];
    RObject new_object = as<RObject>(e["args.object"]);
    List dat = as<List>(rxSolveC(new_object, R_NilValue, extraArgs, new_params, new_events, new_inits, new_scale, new_covs,
				 new_method, new_transit_abs, new_atol, new_rtol, new_maxsteps, new_hmin,
				 new_hmax, new_hini,new_maxordn, new_maxords, new_cores, new_covs_interpolation,
				 new_addCov, new_matrix, new_sigma, new_sigmaDf, new_nCoresRV, new_sigmaIsChol,
				 new_nDisplayProgress, new_amountUnits, new_timeUnits, new_addDosing, new_stateTrim));
    if (updateObject && as<bool>(e[".real.update"])){
      List old = as<List>(rxCurObj);
      //Should I zero out the List...?
      CharacterVector oldNms = old.names();
      CharacterVector nms = dat.names();
      if (oldNms.size() == nms.size()){
	int i;
	for (i = 0; i < nms.size(); i++){
	  old[as<std::string>(nms[i])] = as<SEXP>(dat[as<std::string>(nms[i])]);
	}
	old.attr("class") = dat.attr("class");
	old.attr("row.names") = dat.attr("row.names");
	return old;
      } else {
	warning("Cannot update object...");
	return dat;
      }
    }
    e[".real.update"] = true;
    return dat;
  } else {
    // Load model
    if (!rxDynLoad(object)){
      stop("Cannot load RxODE dlls for this model.");
    }
    // Get model 
    List mv = rxModelVars(object);
    // Assign Pointers
    RxODE_assign_fn_pointers(as<SEXP>(mv));
    // Get the C solve object
    rx_solve* rx = getRxSolve_();
    rx_solving_options* op = rx->op;
    rx_solving_options_ind* ind;

    rx->stateTrim = stateTrim;
    rx->matrix = matrix;
    rx->add_cov = (int)(addCov);
    op->stiff = method;
    if (method != 2){
      op->cores =1;
    } else {
      op->cores=cores;
    }
    // Now set up events and parameters
    RObject par0 = params;
    RObject ev0  = events;
    RObject ev1;
    RObject par1;
    bool swappedEvents = false;
    NumericVector initsC;
    if (rxIs(par0, "rx.event")){
      // Swapped events and parameters
      swappedEvents=true;
      ev1 = par0;
      par1 = ev0;
    } else if (rxIs(ev0, "rx.event")) {
      ev1 = ev0;
      par1 = par0;
    } else {
      Rcout << "Events:\n";
      Rcout << "Parameters:\n";
      stop("Need some event information (observation/dosing times) to solve.\nYou can use either 'eventTable' or an RxODE compatible data.frame/matrix.");
    }
    // Now get the parameters (and covariates)
    //
    // Unspecified parameters can be found in the modVars["ini"]
    NumericVector mvIni = mv["ini"];
    // The event table can contain covariate information, if it is acutally a data frame or matrix.
    Nullable<CharacterVector> covnames0, simnames0;
    CharacterVector covnames, simnames;
    CharacterVector pars = mv["params"];
    int npars = pars.size();
    IntegerVector eGparPos(npars);
    CharacterVector state = mv["state"];
    CharacterVector lhs = mv["lhs"];
    op->neq = state.size();
    op->badSolve = 0;
    op->abort = 0;
    op->ATOL = atol;          //absolute error
    op->RTOL = rtol;          //relative error
    gatol2Setup(op->neq);
    grtol2Setup(op->neq);
    std::fill_n(&_globals.gatol2[0],op->neq, atol);
    std::fill_n(&_globals.grtol2[0],op->neq, rtol);
    op->atol2 = &_globals.gatol2[0];
    op->rtol2 = &_globals.grtol2[0];
    op->H0 = hini;
    op->HMIN = hmin;
    op->mxstep = maxsteps;
    op->MXORDN = maxordn;
    op->MXORDS = maxords;
    // The initial conditions cannot be changed for each individual; If
    // they do they need to be a parameter.
    initsC = rxInits(object, inits, state, 0.0);
    ginitsSetup(initsC.size());
    std::copy(initsC.begin(), initsC.end(), &_globals.ginits[0]);
    op->inits = &_globals.ginits[0];
    NumericVector scaleC = rxSetupScale(object, scale, extraArgs);
    gscaleSetup(scaleC.size());
    std::copy(scaleC.begin(),scaleC.end(),&_globals.gscale[0]);
    op->scale = &_globals.gscale[0];
    //
    int transit = 0;
    if (transit_abs.isNull()){
      transit = mv["podo"];
      if (transit){
        warning("Assumed transit compartment model since 'podo' is in the model.");
      }
    }  else {
      LogicalVector tr = LogicalVector(transit_abs);
      if (tr[0]){
        transit=  1;
      }
    }
    op->do_transit_abs = transit;
    op->nlhs = lhs.size();
    CharacterVector trans = mv["trans"];
    // Make sure the model variables are assigned...
    // This fixes random issues on windows where the solves are done and the data set cannot be solved.
    std::string ptrS = (as<std::string>(trans["model_vars"]));
    getRxModels();
    _rxModels[ptrS] = mv;
    sprintf(op->modNamePtr, "%s", ptrS.c_str());
    // approx fun options
    op->is_locf = covs_interpolation;
    if (op->is_locf == 0){//linear
      op->f1 = 1.0;
      op->f2 = 0.0;
      op->kind=1;
    } else if (op->is_locf == 1){ // locf
      op->f2 = 0.0;
      op->f1 = 1.0;
      op->kind = 0;
      op->is_locf=1;
    } else if (op->is_locf == 2){ //nocb
      op->f2 = 1.0;
      op->f1 = 0.0;
      op->kind = 0;
    }  else if (op->is_locf == 3){ // midpoint
      op->f1 = op->f2 = 0.5;
      op->kind = 0;
    } else {
      stop("Unknown covariate interpolation specified.");
    }
    op->extraCmt=0;
    op->nDisplayProgress = nDisplayProgress;
    op->ncoresRV = nCoresRV;
    op->isChol = (int)(sigmaIsChol);
    unsigned int nsub = 0;
    unsigned int nobs = 0, ndoses = 0;
    unsigned int i, j, k = 0;
    int ncov =0, curcovi = 0;
    double tmp, hmax1 = 0.0, hmax2 = 0.0, tlast;
    // Covariate options
    // Simulation variabiles
    // int *svar;
    CharacterVector sigmaN;
    bool usePar1 = false;
    bool simSubjects = false;
    if (!thetaMat.isNull() || !omega.isNull() || !sigma.isNull()){
      // Simulated Variable3
      if (!rxIs(par1, "numeric")){
        stop("When specifying 'thetaMat', 'omega', or 'sigma' the parameters cannot be a data.frame/matrix.");
      }
      unsigned int nSub0 = 0;
      int curObs = 0;
      rx->nall = 0;
      rx->nobs = 0;
      if(rxIs(ev1, "EventTable")){
	nSub0 = 1;
	List et = List(ev1);
	Function f = et["get.EventTable"];
	DataFrame dataf = f();
	IntegerVector evid = as<IntegerVector>(dataf[1]);
	rx->nall = evid.size();
	for (unsigned int j = rx->nall;j--;){
	  if (evid[j] == 0) rx->nobs++;
	}
      } else if (rxIs(ev1,"event.data.frame")||
		 rxIs(ev1,"event.matrix")){
	if (rxcId > -1){
	  DataFrame dataf = as<DataFrame>(ev1);
          IntegerVector id = as<IntegerVector>(dataf[rxcId]);
	  IntegerVector evid  = as<IntegerVector>(dataf[rxcEvid]);
	  int lastid= id[id.size()-1]+42;
	  rx->nall = evid.size();
	  for (unsigned int j = rx->nall; j--;){
	    if (lastid != id[j]){
	      lastid=id[j];
	      nSub0++;
	    }
	    if (evid[j] == 0) rx->nobs++;
	  }
	} else {
	  nSub0 =1;
	  DataFrame dataf = as<DataFrame>(ev1);
          IntegerVector evid  = as<IntegerVector>(dataf[rxcEvid]);
          rx->nall = evid.size();
          for (unsigned int j =rx->nall; j--;){
            if (evid[j] == 0) rx->nobs++;
          }
	}
      }
      if (nSub > 1 && nSub0 > 1 && nSub != nSub0){
        stop("You provided multi-subject data and asked to simulate a different number of subjects;  I don't know what to do.");
      } else if (nSub > 1 && nSub0 == 1) {
	nSub0 = nSub;
        simSubjects = true;
      } 
      curObs = addDosing ? rx->nall : rx->nobs;
      par1 =  as<RObject>(rxSimThetaOmega(as<Nullable<NumericVector>>(par1), omega, omegaDf, omegaIsChol, nSub0, thetaMat, thetaDf, thetaIsChol, nStud,
                                          sigma, sigmaDf, sigmaIsChol, nCoresRV, curObs, dfSub, dfObs, simSubjects));
      usePar1=true;
      
      // The parameters are in the same format as they would be if they were specified as part of the original dataset.
    }
    // .sigma could be reassigned in an update, so check outside simulation function.
    if (_rxModels.exists(".sigma")){
      sigmaN = as<CharacterVector>((as<List>((as<NumericMatrix>(_rxModels[".sigma"])).attr("dimnames")))[1]);
    }
    int parType = 1;
    NumericVector parNumeric;
    DataFrame parDf;
    NumericMatrix parMat;
    CharacterVector nmP;
    int nPopPar = 1;
    if (!theta.isNULL() || !eta.isNULL()){
      par1 = rxSetupParamsThetaEta(par1, theta, eta);
    }
    if (rxIs(par1, "numeric") || rxIs(par1, "integer")){
      parNumeric = as<NumericVector>(par1);
      if (parNumeric.hasAttribute("names")){
        nmP = parNumeric.names();
      } else if (parNumeric.size() == pars.size()) {
        nmP = pars;
      } else {
        stop("If parameters are not named, they must match the order and size of the parameters in the model.");
      }
    } else if (rxIs(par1, "data.frame")){
      parDf = as<DataFrame>(par1);
      parType = 2;
      nmP = parDf.names();
      nPopPar = parDf.nrows();
    } else if (rxIs(par1, "matrix")){
      parMat = as<NumericMatrix>(par1);
      nPopPar = parMat.nrow();
      parType = 3;
      if (parMat.hasAttribute("dimnames")){
        Nullable<CharacterVector> colnames0 = as<Nullable<CharacterVector>>((as<List>(parMat.attr("dimnames")))[1]);
        if (colnames0.isNull()){
          if (parMat.ncol() == pars.size()){
            nmP = pars;
          } else {
            stop("If parameters are not named, they must match the order and size of the parameters in the model.");
          }
        } else {
          nmP = CharacterVector(colnames0);
        }
      } else if (parMat.ncol() == pars.size()) {
        nmP = pars;
      } else {
        stop("If parameters are not named, they must match the order and size of the parameters in the model.");
      }
    }
    rxOptionsIniEnsure(nPopPar); // 1 simulation per parameter specifcation
    if (rxIs(ev1, "EventTable")){
      List et = List(ev1);
      Function f = et["get.EventTable"];
      DataFrame dataf = f();
      NumericVector time = as<NumericVector>(dataf[0]);
      IntegerVector evid = as<IntegerVector>(dataf[1]);
      NumericVector amt  = as<NumericVector>(dataf[2]);
      ind = &(rx->subjects[0]);
      ind->id=1;
      ind->lambda=1.0;
      ind->yj = 0;
      rx->nsub= 1;
      // Time copy
      ind->n_all_times   = time.size();
      gall_timesSetup(ind->n_all_times);
      std::copy(time.begin(), time.end(), &_globals.gall_times[0]);
      ind->all_times     = &_globals.gall_times[0];
      // EVID copy
      gevidSetup(ind->n_all_times);
      std::copy(evid.begin(),evid.end(), &_globals.gevid[0]);
      ind->evid     = &_globals.gevid[0];
      gamtSetup(ind->n_all_times);
      ind->dose = &_globals.gamt[0];
      // Slower; These need to be built.
      tlast = time[0];
      // hmax1 = hmax2 = 0;
      gidoseSetup(ind->n_all_times);
      ind->idose = &_globals.gidose[0];
      j=0;
      for (i =0; i != (unsigned int)(ind->n_all_times); i++){
        if (ind->evid[i]){
          ndoses++;
	  _globals.gidose[j] = (int)i;
          _globals.gamt[j++] = amt[i];
	} else {
	  nobs++;
	  tmp = time[i]-tlast;
	  if (tmp < 0) stop("Dataset must be ordered by ID and TIME variables .");
	  if (tmp > hmax1){
	    hmax1 = tmp;
	  }
	}
      }
      ind->ndoses = ndoses;
      rx->nobs = nobs;
      rx->nall = nobs+ndoses;
      if (!hmax.isNull()){
	NumericVector hmax0(hmax);
	ind->HMAX = hmax0[0];
	op->hmax2 = hmax0[0];
      } else {
        ind->HMAX = hmax1;
	op->hmax2 = hmax1;
      }
      nsub=1;
      if (!covs.isNULL()){
        // op->do_par_cov = 1;
	op->do_par_cov = 1;
	ncov = 0;
        if (rxIs(covs,"data.frame")){
	  List df = as<List>(covs);
          CharacterVector dfNames = df.names();
	  int dfN = dfNames.size();
	  gcovSetup(dfN * ind->n_all_times);
	  gpar_covSetup(dfN);
	  k = 0;
	  for (i = dfN; i--;){
	    for (j = npars; j--;){
	      if (pars[j] == dfNames[i]){
		_globals.gpar_cov[k] = j+1;
		// Not clear if this is an integer/real.  Copy the values.
		NumericVector cur = as<NumericVector>(df[i]);
		std::copy(cur.begin(), cur.end(), _globals.gcov+curcovi);
		curcovi += ind->n_all_times;
                ncov++;
                k++;
                break;
              }
            }
          }
          op->ncov=ncov;
          ind->cov_ptr = &(_globals.gcov[0]);
	  op->par_cov=&(_globals.gpar_cov[0]);
        } else if (rxIs(covs, "matrix")){
	  // FIXME
	  stop("Covariates must be supplied as a data.frame.");
	} 
      } else {
	op->ncov = 0;
        // int *par_cov;
        op->do_par_cov = 0;
      }
    } else if (rxIs(ev1,"event.data.frame")||
               rxIs(ev1,"event.matrix")){
      // data.frame or matrix
      double hmax0 = 0.0;
      if (!hmax.isNull()){
        NumericVector hmax00(hmax);
        hmax0 = hmax00[0];
      }
      DataFrame dataf = as<DataFrame>(ev1);
      // Copy the information and make sure there is enough room.
      // - evid
      IntegerVector evid  = as<IntegerVector>(dataf[rxcEvid]);
      gevidSetup(evid.size());
      std::copy(evid.begin(),evid.end(), &_globals.gevid[0]);
      IntegerVector id(evid.size(), 1);
      if (rxcId > -1){
        id    = as<IntegerVector>(dataf[rxcId]);
      }
      NumericVector time0 = dataf[rxcTime];
      tlast = time0[0];
      // - all_times
      gall_timesSetup(time0.size());
      std::copy(time0.begin(), time0.end(), &_globals.gall_times[0]);
      NumericVector dv;
      if (rxcDv > -1){
	dv = as<NumericVector>(dataf[rxcDv]);
        gdvSetup(dv.size());
        std::copy(dv.begin(), dv.end(), &_globals.gdv[0]);
      }
      NumericVector amt   = dataf[rxcAmt];
      // Make sure that everything can take the correct sizes
      // - amt
      gamtSetup(amt.size());
      // - idose
      gidoseSetup(amt.size());
      // Get covariates
      CharacterVector dfNames = dataf.names();
      int dfN = dfNames.size();
      // - par cov needs to be at lest the size of the dataframe names.
      gpar_covSetup(dfN);
      ncov = 0;
      std::vector<int> covPos(dfN);
      curcovi=0;
      for (i = dfN; i--;){
        for (j = npars; j--;){
          if (pars[j] == dfNames[i]){
            _globals.gpar_cov[ncov] = j+1;
	    covPos[ncov] = i;
            ncov++;
	  }
	}
      }
      op->par_cov=&(_globals.gpar_cov[0]);
      op->ncov=ncov;
      op->do_par_cov = (ncov > 0);
      // Make sure the covariates are a #ncov * all times size
      gcovSetup(ncov * amt.size());
      unsigned int ids = id.size();
      // Get the number of subjects
      // Get the number of observations
      // Get the number of doses
      unsigned int nall = 0, nobst=0, lasti =0, ii=0;
      nsub = 0;
      ind = &(rx->subjects[0]);
      j=0;
      int lastId = id[0]-42;
      for (i = 0; i < ids; i++){
        if (lastId != id[i]){
	  if (nall != 0){
            // Finalize last solve.
            ind->n_all_times    = ndoses+nobs;
            ind->cov_ptr = &(_globals.gcov[curcovi]);
	    for (ii = 0; ii < (unsigned int)ncov; ii++){
	      NumericVector cur = as<NumericVector>(dataf[covPos[ii]]);
	      std::copy(cur.begin()+lasti, cur.begin()+lasti+ind->n_all_times,
                        _globals.gcov+curcovi);
              curcovi += ind->n_all_times;
	    }
            nsub++;
            ind = &(rx->subjects[nsub]);
          }
	  // Setup the pointers.
          ind->id             = nsub+1;
          ind->lambda         =1.0;
          ind->yj             = 0;
          ind->all_times      = &_globals.gall_times[i];
	  if (rxcDv > -1){
	    ind->dv = &_globals.gdv[i];
	  }
          ind->evid           = &_globals.gevid[i];
	  ind->idose          = &_globals.gidose[i];
          ind->dose           = &_globals.gamt[i];
	  lasti = i;
	  if (hmax0 == 0.0){
            ind->HMAX = hmax1;
	  } else {
	    ind->HMAX = hmax0;
	  }
	  hmax1 = 0.0;
          lastId=id[i];
	  j=i;
	  ndoses=0;
	  nobs=0;
	  tlast = NA_REAL;
        }
        if (_globals.gevid[i]){
          _globals.gidose[j] = i-lasti;
          _globals.gamt[j] = amt[i];
	  ind->ndoses++;
          ndoses++; nall++; j++;
        } else {
          nobs++; nobst++; nall++;
	  if (!ISNA(tlast)){
            tmp = time0[i]-tlast;
            if (tmp < 0) stop("Dataset must be ordered by ID and TIME variables", ind->id, tlast, time0[i]);
            if (tmp > hmax1){
              hmax1 = tmp;
              if (hmax1 > hmax2){
                hmax2=hmax1;
              }
            }
	  }
        }
	tlast = time0[i];
      }
      rx->nobs = nobst;
      rx->nall = nall;
      // Finalize the prior individual
      ind->n_all_times    = ndoses+nobs;
      ind->cov_ptr = &(_globals.gcov[curcovi]);
      for (ii = 0; ii < (unsigned int)ncov; ii++){
        NumericVector cur = as<NumericVector>(dataf[covPos[ii]]);
        std::copy(cur.begin()+lasti, cur.begin()+lasti+ind->n_all_times,
                  _globals.gcov+curcovi);
        curcovi += ind->n_all_times;
      }
      if (hmax0 == 0.0){
        op->hmax2 = hmax2;
      } else {
        op->hmax2 = hmax0;
      }
      nsub++;
      rx->nsub= nsub;
    }
    // Make sure the user input all the parmeters.
    gParPosSetup(npars);
    std::string errStr = "";
    bool allPars = true;
    bool curPar = false;
    CharacterVector mvIniN = mvIni.names();
    gsvarSetup(sigmaN.size());
    for (i = npars; i--;){
      curPar = false;
      // Check to see if this is a covariate.
      for (j = op->ncov; j--;){
	if (_globals.gpar_cov[j] == (int)(i + 1)){
	  _globals.gParPos[i] = 0; // These are set at run-time and "dont" matter.
	  curPar = true;
          eGparPos[i]=_globals.gParPos[i];
	  break;
	}
      }
      // Check for the sigma-style simulated parameters.
      if (!curPar){
	for (j = sigmaN.size(); j--;){
          if (sigmaN[j] == pars[i]){
	    _globals.gsvar[j] = i;
	    _globals.gParPos[i] = 0; // These are set at run-time and "dont" matter.
	    curPar = true;
            eGparPos[i]=_globals.gParPos[i];
            break;
	  }
	}
      }
      // Next, check to see if this is a user-specified parameter
      if (!curPar){
        for (j = nmP.size(); j--;){
          if (nmP[j] == pars[i]){
            curPar = true;
            _globals.gParPos[i] = j + 1;
            eGparPos[i]=_globals.gParPos[i];
            break;
          }
        }
      }
      // last, check for $ini values
      if (!curPar){
        for (j = mvIniN.size(); j--;){
          if (mvIniN[j] == pars[i]){
            curPar = true;
            _globals.gParPos[i] = -j - 1;
            eGparPos[i]=_globals.gParPos[i];
            break;
          }
        }
      }
      if (!curPar){
        if (errStr == ""){
          errStr = "The following parameter(s) are required for solving: " + pars[i];
        } else {
          errStr = errStr + ", " + pars[i];
        }
        allPars = false;
      }
    }
    if (!allPars){
      CharacterVector modSyntax = mv["model"];
      Rcout << "Model:\n\n" + modSyntax[0] + "\n";
      stop(errStr);
    }
    op->svar = &_globals.gsvar[0];
    // Now setup the rest of the rx_solve object
    if (nPopPar != 1 && nPopPar % rx->nsub != 0){
      stop("The number of parameters (%d) solved by RxODE for multi-subject data needs to be a multiple of the number of subjects (%d).",nPopPar, rx->nsub);
    }
    int nSize = nPopPar*rx->nsub;
    if (nPopPar == 1) nSize = rx->nsub;

    gInfusionRateSetup(op->neq*nSize);
    std::fill_n(&_globals.gInfusionRate[0], op->neq*nSize, 0.0);

    gBadDoseSetup(op->neq*nSize);
    std::fill_n(&_globals.gBadDose[0], op->neq*nSize, 0);
    
    grcSetup(nSize);
    std::fill_n(&_globals.grc[0], nSize, 0);

    gslvr_counterSetup(nSize);
    std::fill_n(&_globals.slvr_counter[0], nSize, 0);
    
    gdadt_counterSetup(nSize);
    std::fill_n(&_globals.dadt_counter[0], nSize, 0);

    gjac_counterSetup(nSize);
    std::fill_n(&_globals.jac_counter[0], nSize, 0);

    glhsSetup(lhs.size()*nSize);
    std::fill_n(&_globals.glhs[0],lhs.size()*nSize,0.0);

    rx->nsim = nPopPar / rx->nsub;
    if (rx->nsim < 1) rx->nsim=1;

    gsolveSetup(rx->nall*state.size()*rx->nsim);
    std::fill_n(&_globals.gsolve[0], rx->nall*state.size()*rx->nsim, 0.0);
    int curEvent = 0;
    
    switch(parType){
    case 1: // NumericVector
      if (nPopPar != 1) stop("Something is wrong... nPopPar != 1 but parameters are specified as a NumericVector.");
      gparsSetup(npars);
      for (i = npars; i--;){
	if (_globals.gParPos[i] == 0){ // Covariate or simulated variable.
	  _globals.gpars[i] = 0;//NA_REAL;
	} else if (_globals.gParPos[i] > 0){ // User specified parameter
	  _globals.gpars[i] = parNumeric[_globals.gParPos[i]-1];
	} else { // ini specified parameter.
	  _globals.gpars[i] = mvIni[-_globals.gParPos[i]-1];
	}
      }
      for (i = rx->nsub; i--;){
	ind = &(rx->subjects[i]);
	ind->par_ptr = &_globals.gpars[0];
	ind->InfusionRate = &_globals.gInfusionRate[op->neq*i];
        ind->BadDose = &_globals.gBadDose[op->neq*i];
        ind->nBadDose = 0;
	// Hmax defined above.
	ind->tlast=0.0;
	ind->podo = 0.0;
	ind->ixds =  0;
	ind->sim = i+1;
        ind->solve = &_globals.gsolve[curEvent];
        curEvent += op->neq*ind->n_all_times;
        ind->lhs = &_globals.glhs[i*lhs.size()];
	ind->rc = &_globals.grc[i];
        ind->slvr_counter = &_globals.slvr_counter[i];
        ind->dadt_counter = &_globals.dadt_counter[i];
        ind->jac_counter  = &_globals.jac_counter[i];
      }
      break;
    case 2: // DataFrame
      // Convert to NumericMatrix
      {
	unsigned int parDfi = parDf.size();
	parMat = NumericMatrix(nPopPar,parDfi);
	while (parDfi--){
	  parMat(_,parDfi)=NumericVector(parDf[parDfi]);
	}
	parMat.attr("dimnames") = List::create(R_NilValue, parDf.names());
      }
    case 3: // NumericMatrix
      {
	gparsSetup(npars*nPopPar);
        for (unsigned int j = 0; j < (unsigned int)nPopPar; j++){
	  for (unsigned int k = 0; k < (unsigned int)npars; k++){
	    i = k+npars*j;
	    if (_globals.gParPos[k] == 0){
	      _globals.gpars[i] = 0;
	    } else if (_globals.gParPos[k] > 0){
	      // posPar[i] = j + 1;
	      _globals.gpars[i] = parMat(j, _globals.gParPos[k]-1);
	    } else {
	      // posPar[i] = -j - 1;
	      _globals.gpars[i] = mvIni[-_globals.gParPos[k]-1];
	    }
	  }
	}
	curEvent=0;
	rx_solving_options_ind indS;
	for (unsigned int simNum = rx->nsim; simNum--;){
	  for (unsigned int id = rx->nsub; id--;){
	    unsigned int cid = id+simNum*rx->nsub;
	    ind = &(rx->subjects[cid]);
	    ind->par_ptr = &_globals.gpars[cid*npars];
	    ind->InfusionRate = &_globals.gInfusionRate[op->neq*cid];
	    ind->BadDose = &_globals.gBadDose[op->neq*cid];
	    ind->nBadDose = 0;
	    // Hmax defined above.
	    ind->tlast=0.0;
	    ind->podo = 0.0;
	    ind->ixds =  0;
	    ind->sim = simNum+1;
	    ind->lhs = &_globals.glhs[cid];
	    ind->rc = &_globals.grc[cid];
	    ind->slvr_counter = &_globals.slvr_counter[cid];
	    ind->dadt_counter = &_globals.dadt_counter[cid];
	    ind->jac_counter = &_globals.jac_counter[cid];
	    if (simNum){
	      // Assign the pointers to the shared data
	      indS = rx->subjects[id];
	      if (op->do_par_cov){
		ind->cov_ptr = indS.cov_ptr;
	      }
	      ind->n_all_times =indS.n_all_times;
	      ind->HMAX = indS.HMAX;
	      ind->idose = &(indS.idose[0]);
	      ind->ndoses = indS.ndoses;
	      ind->dose = &(indS.dose[0]);
	      ind->evid =&(indS.evid[0]);
	      ind->all_times = &(indS.all_times[0]);
	      ind->id=id+1;
              ind->lambda=1.0;
              ind->yj = 0;
	    }
	    int eLen = op->neq*ind->n_all_times;
	    ind->solve = &_globals.gsolve[curEvent];
	    curEvent += eLen;
	  }
	}
      }
      break;
    default: 
      stop("Something is wrong here.");
    }
    if (setupOnly){
      return as<SEXP>(LogicalVector::create(true));
    }
    par_solve(rx);
    if (op->abort){
      rxSolveFree();
      stop("Aborted solve.");
    }
    int doDose = 0;
    if (addDosing){
      doDose = 1;
    } else {
      doDose = 0;
    }
    IntegerVector si = mv["state.ignore"];
    rx->stateIgnore = &si[0];
    int doTBS = (rx->matrix == 3);
    if (doTBS) rx->matrix=2;
    List dat = RxODE_df(doDose, doTBS);
    dat.attr("class") = CharacterVector::create("data.frame");
    List xtra;
    int nr = rx->nr;
    int nc = dat.size();
    if (rx->matrix){
      rxSolveFree();
      if(_rxModels.exists(".sigma")){
      	_rxModels.remove(".sigma");
      }
      if(_rxModels.exists(".sigmaL")){
      	_rxModels.remove(".sigmaL");
      }
      if(_rxModels.exists(".omegaL")){
      	_rxModels.remove(".omegaL");
      }
      if(_rxModels.exists(".theta")){
      	_rxModels.remove(".theta");
      }
      if (rx->matrix == 2){
        dat.attr("class") = "data.frame";
	return dat;
      } else {
        NumericMatrix tmpM(nr,nc);
        for (unsigned int i = dat.size(); i--;){
          tmpM(_,i) = as<NumericVector>(dat[i]);
        }
        tmpM.attr("dimnames") = List::create(R_NilValue,dat.names());
        return tmpM;
      }
    } else {
      IntegerVector slvr_counterIv(nSize);
      IntegerVector dadt_counterIv(nSize);
      IntegerVector  jac_counterIv(nSize);
      std::copy(&_globals.slvr_counter[0], &_globals.slvr_counter[0] + nSize, slvr_counterIv.begin());
      std::copy(&_globals.dadt_counter[0], &_globals.dadt_counter[0] + nSize, dadt_counterIv.begin());
      std::copy(&_globals.jac_counter[0], &_globals.jac_counter[0] + nSize, jac_counterIv.begin());

      rxSolveFree();
      
      Function newEnv("new.env", R_BaseNamespace);
      Environment e = newEnv(_["size"] = 29, _["parent"] = RxODE());
      getRxModels();
      if(_rxModels.exists(".theta")){
  	e[".theta"] = as<NumericMatrix>(_rxModels[".theta"]);
  	_rxModels.remove(".theta");
      }
      if(_rxModels.exists(".sigma")){
  	e[".sigma"] = as<NumericMatrix>(_rxModels[".sigma"]);
  	_rxModels.remove(".sigma");
      }
      if(_rxModels.exists(".omegaL")){
  	e[".omegaL"] = as<List>(_rxModels[".omegaL"]);
  	_rxModels.remove(".omegaL");
      }
      if(_rxModels.exists(".sigmaL")){
  	e[".sigmaL"] = as<List>(_rxModels[".sigmaL"]);
  	_rxModels.remove(".sigmaL");
      }
      e["check.nrow"] = nr;
      e["check.ncol"] = nc;
      e["check.names"] = dat.names();
      
      e[".par.pos"] = eGparPos;
      e[".slvr.counter"] = slvr_counterIv;
      e[".dadt.counter"] = dadt_counterIv;
      e[".jac.counter"] = jac_counterIv;
      e[".nsub"] = rx->nsub;
      e[".nsim"] = rx->nsim;
      e["inits.dat"] = initsC;
      CharacterVector units = CharacterVector::create(amountUnits[0], timeUnits[0]);
      units.names() = CharacterVector::create("dosing","time");
      e["units"] = units;
      e["nobs"] = rx->nobs;
    
      e["args.object"] = object;
      e["dll"] = rxDll(object);
      if (!swappedEvents){
	if (usePar1){
          e["args.params"] = par1;
	} else {
	  e["args.params"] = params;
	}
  	e["args.events"] = events;
      } else {
        if (usePar1){
          e["args.params"] = par1;
        } else {
          e["args.params"] = events;
        }
  	e["args.events"] = params;
      }
      e["args.inits"] = inits;
      e["args.covs"] = covs;
      e["args.method"] = method;
      e["args.transit_abs"] = transit_abs;
      e["args.atol"] = atol;
      e["args.rtol"] = rtol;
      e["args.maxsteps"] = maxsteps;
      e["args.hmin"] = hmin;
      e["args.hmax"] = hmax;
      e["args.hini"] = hini;
      e["args.maxordn"] = maxordn;
      e["args.maxords"] = maxords;
      e["args.cores"] = cores;
      e["args.covs_interpolation"] = covs_interpolation;
      e["args.addCov"] = addCov;
      e["args.matrix"] = matrix;
      e["args.sigma"] = sigma;
      e["args.sigmaDf"] = sigmaDf;
      e["args.nCoresRV"] = nCoresRV;
      e["args.sigmaIsChol"] = sigmaIsChol;
      e["args.nDisplayProgress"] = nDisplayProgress;
      e["args.amountUnits"] = amountUnits;
      e["args.timeUnits"] = timeUnits;
      e["args.addDosing"] = addDosing;
      e["args.stateTrim"] = stateTrim;
      e[".real.update"] = true;
      CharacterVector cls= CharacterVector::create("rxSolve", "data.frame");
      cls.attr(".RxODE.env") = e;    
      dat.attr("class") = cls;
      return(dat);
    }
  }
  return R_NilValue;
}

//[[Rcpp::export]]
SEXP rxSolveCsmall(const RObject &object,
                   const Nullable<CharacterVector> &specParams = R_NilValue,
                   const Nullable<List> &extraArgs = R_NilValue,
                   const RObject &params = R_NilValue,
                   const RObject &events = R_NilValue,
                   const RObject &inits = R_NilValue,
		   const RObject &scale = R_NilValue,
                   const RObject &covs  = R_NilValue,
                   const Nullable<List> &optsL = R_NilValue){
  if (optsL.isNull()){
    stop("Not meant to be called directly.  Needs options setup.");
  }
  List opts = List(optsL);
  return rxSolveC(object, specParams, extraArgs, params, events, inits, scale, covs,
                  (int)opts[0], // const int method = 2,
                  as<Nullable<LogicalVector>>(opts[1]), // const Nullable<LogicalVector> &transit_abs = R_NilValue,
                  (double)opts[2], //const double atol = 1.0e-8,
                  (double)opts[3], // const double rtol = 1.0e-6,
                  (int)opts[4], //const int maxsteps = 5000,
                  (double)opts[5], //const double hmin = 0,
                  as<Nullable<NumericVector>>(opts[6]), //const Nullable<NumericVector> &hmax = R_NilValue,
                  (double)opts[7], //const double hini = 0,
                  (int)opts[8], //const int maxordn = 12,
                  opts[9], //const int maxords = 5,
                  (unsigned int)opts[10], //const int cores = 1,
                  (int)opts[11], //const int covs_interpolation = 0,
                  opts[12], //bool addCov = false,
                  opts[13], //bool matrix = false,
                  opts[14], //const RObject &sigma= R_NilValue,
                  opts[15], //const RObject &sigmaDf= R_NilValue,
                  (int)opts[16], //const int &nCoresRV= 1,
                  opts[17], //const bool &sigmaIsChol= false,
                  opts[18], // nDisplayProgress
                  opts[19], //const CharacterVector &amountUnits = NA_STRING,
                  opts[20], //const CharacterVector &timeUnits = "hours",
                  opts[21], //const bool updateObject = false
                  opts[22], //const RObject &eta = R_NilValue,
                  opts[23], //const bool addDosing = false
		  opts[24],//const double stateTrim = R_PosInf
		  opts[25], // 
		  opts[26], //
		  opts[27], // const RObject &omega = R_NilValue, 
		  opts[28], // const RObject &omegaDf = R_NilValue, 
                  opts[29], // const bool &omegaIsChol = false,
                  opts[30], // const int nSub = 1, 
                  opts[31], // const RObject &thetaMat = R_NilValue, 
                  opts[32], // const RObject &thetaDf = R_NilValue, 
                  opts[33], // const bool &thetaIsChol = false,
                  opts[34], // const int nStud = 1, 
                  opts[35], // const dfSub
                  opts[36], // dfObs
		  opts[37] // setupOnly
		  );
}

//[[Rcpp::export]]
RObject rxSolveGet(RObject obj, RObject arg, LogicalVector exact = true){
  std::string sarg;
  int i, j, n;
  if (rxIs(obj, "data.frame")){
    List lst = as<List>(obj);
    if (rxIs(arg, "character")){
      sarg = as<std::string>(arg);
      CharacterVector nm = lst.names();
      n = nm.size();
      unsigned int slen = strlen(sarg.c_str());
      int dexact = -1;
      if (exact[0] == TRUE){
	dexact = 1;
      } else if (exact[0] == FALSE){
	dexact = 0;
      }
      unsigned int slen2;
      for (i = 0; i < n; i++){
	slen2 = strlen((as<std::string>(nm[i])).c_str());
	if (slen <= slen2 &&
	    (strncmp((as<std::string>(nm[i])).c_str(), sarg.c_str(), slen)  == 0 ) &&
	    (dexact != 1 || (dexact == 1 && slen == slen2))){
	  if (dexact == -1){
	    warning("partial match of '%s' to '%s'",sarg.c_str(), (as<std::string>(nm[i])).c_str());
	  }
	  return lst[i];
	}
      }
      if (rxIs(obj, "rxSolve")){
	rxCurObj = obj;
	CharacterVector cls = lst.attr("class");
	Environment e = as<Environment>(cls.attr(".RxODE.env"));
	if (sarg == "env"){
	  return as<RObject>(e);
	}
	if (sarg == "model"){
	  List mv = rxModelVars(obj);
	  CharacterVector mods = mv["model"];
	  CharacterVector retS = as<std::string>(mods["model"]);
	  retS.attr("class") = "RxODE.modeltext";
	  return(retS);
	}
        updateSolveEnvPost(e);
        if (e.exists(sarg)){
          return e[sarg];
        }
        if (sarg == "params" || sarg == "par" || sarg == "pars" || sarg == "param"){
          return e["params.dat"];
	} else if (sarg == "inits" || sarg == "init"){
	  return e["inits.dat"];
	} else if (sarg == "t"){
	  return lst["time"];
	} else if ((sarg == "theta.mat" || sarg == "thetaMat") && e.exists(".theta")){
	  return e[".theta"];
	} else if ((sarg == "sigma.list" || sarg == "sigmaList") && e.exists(".sigmaL")){
	  return e[".sigmaL"];
	} else if ((sarg == "omega.list" || sarg == "omegaList") && e.exists(".omegaL")){
          return e[".omegaL"];
	}
	// Now parameters
	List pars = List(e["params.dat"]);
	CharacterVector nmp = pars.names();
	n = pars.size();
	for (i = n; i--;){
	  if (nmp[i] == sarg){
	    return pars[sarg];
	  }
	}
	// // Now inis.
	// Function sub("sub", R_BaseNamespace);
	NumericVector ini = NumericVector(e["inits.dat"]);
	CharacterVector nmi = ini.names();
	n = ini.size();
        std::string cur;
        NumericVector retN(1);
        for (i = n; i--; ){
	  cur = as<std::string>(nmi[i]) + "0";
	  if (cur == sarg){
	    retN = ini[i];
	    return as<RObject>(retN);
	  }
	  cur = as<std::string>(nmi[i]) + "_0";
          if (cur == sarg){
	    retN = ini[i];
            return as<RObject>(retN);
          }
          cur = as<std::string>(nmi[i]) + ".0";
          if (cur == sarg){
            retN = ini[i];
            return as<RObject>(retN);
          }
          cur = as<std::string>(nmi[i]) + "[0]";
          if (cur == sarg){
	    retN = ini[i];
            return as<RObject>(retN);
          }
          cur = as<std::string>(nmi[i]) + "(0)";
          if (cur == sarg){
	    retN = ini[i];
	    return as<RObject>(retN);
          }
          cur = as<std::string>(nmi[i]) + "{0}";
          if (cur == sarg){
	    retN = ini[i];
            return as<RObject>(retN);
          }
	}
	List mv = rxModelVars(obj);
	if (sarg == "rx" || sarg == "rxode" || sarg == "RxODE"){
          CharacterVector trans = mv["trans"];
	  getRxModels();
	  std::string pre = as<std::string>(trans["prefix"]);
          if (_rxModels.exists(pre)){
	    return as<RObject>(_rxModels[pre]);
	  }
	}
	CharacterVector normState = mv["normal.state"];;
	CharacterVector parsC = mv["params"];
        CharacterVector lhsC = mv["lhs"];
	for (i = normState.size(); i--;){
	  for (j = parsC.size(); j--; ){
	    std::string test = "_sens_" + as<std::string>(normState[i]) + "_" + as<std::string>(parsC[j]);
	    if (test == sarg){
	      test = "rx__sens_" + as<std::string>(normState[i]) + "_BY_" + as<std::string>(parsC[j]) + "__";
	      return lst[test];
	    }
            test = as<std::string>(normState[i]) + "_" + as<std::string>(parsC[j]);
            if (test == sarg){
              test = "rx__sens_" + as<std::string>(normState[i]) + "_BY_" + as<std::string>(parsC[j]) + "__";
	      return lst[test];
	    }
            test = as<std::string>(normState[i]) + "." + as<std::string>(parsC[j]);
            if (test == sarg){
              test = "rx__sens_" + as<std::string>(normState[i]) + "_BY_" + as<std::string>(parsC[j]) + "__";
              return lst[test];
            }
	  }
          for (j = lhsC.size(); j--;){
            std::string test = "_sens_" + as<std::string>(normState[i]) + "_" + as<std::string>(lhsC[j]);
            if (test == sarg){
              test = "rx__sens_" + as<std::string>(normState[i]) + "_BY_" + as<std::string>(lhsC[j]) + "__";
              return lst[test];
            }
            test = as<std::string>(normState[i]) + "_" + as<std::string>(lhsC[j]);
            if (test == sarg){
              test = "rx__sens_" + as<std::string>(normState[i]) + "_BY_" + as<std::string>(lhsC[j]) + "__";
              return lst[test];
            }
            test = as<std::string>(normState[i]) + "." + as<std::string>(lhsC[j]);
            if (test == sarg){
              test = "rx__sens_" + as<std::string>(normState[i]) + "_BY_" + as<std::string>(lhsC[j]) + "__";
              return lst[test];
            }
          }
	}
      }
    } else {
      if (rxIs(arg, "integer") || rxIs(arg, "numeric")){
	int iarg = as<int>(arg);
	if (iarg <= lst.size()){
	  return lst[iarg-1];
	}
      }
    }
  }
  return R_NilValue;
}

//[[Rcpp::export]]
RObject rxSolveUpdate(RObject obj,
		      RObject arg = R_NilValue,
		      RObject value = R_NilValue){
  if (rxIs(obj,"rxSolve")){
    rxCurObj = obj;
    if (rxIs(arg,"character")){
      CharacterVector what = CharacterVector(arg);
      // CharacterVector classattr = obj.attr("class");
      // Environment e = as<Environment>(classattr.attr(".RxODE.env"));
      // updateSolveEnvPost(e);
      if (what.size() == 1){
	std::string sarg = as<std::string>(what[0]);
	// Now check to see if this is something that can be updated...
	if (sarg == "params"){
	  return rxSolveC(obj,
                          CharacterVector::create("params"),
			  R_NilValue,
                          value, //defrx_params,
                          defrx_events,
                          defrx_inits,
			  R_NilValue, // scale (cannot be updated currently.)
                          defrx_covs,
                          defrx_method,
                          defrx_transit_abs,
                          defrx_atol,
                          defrx_rtol,
                          defrx_maxsteps,
                          defrx_hmin,
                          defrx_hmax,
                          defrx_hini,
                          defrx_maxordn,
                          defrx_maxords,
                          defrx_cores,
                          defrx_covs_interpolation,
                          defrx_addCov,
                          defrx_matrix,
                          defrx_sigma,
                          defrx_sigmaDf,
                          defrx_nCoresRV,
                          defrx_sigmaIsChol,
                          defrx_nDisplayProgress,
                          defrx_amountUnits,
                          defrx_timeUnits, defrx_addDosing, defrx_stateTrim);
	} else if (sarg == "events"){
	  return rxSolveC(obj,
			  CharacterVector::create("events"),
			  R_NilValue,
			  defrx_params,
			  value, // defrx_events,
			  defrx_inits,
			  R_NilValue, // scale
			  defrx_covs,
			  defrx_method,
			  defrx_transit_abs,
			  defrx_atol,
			  defrx_rtol,
			  defrx_maxsteps,
			  defrx_hmin,
			  defrx_hmax,
			  defrx_hini,
			  defrx_maxordn,
			  defrx_maxords,
			  defrx_cores,
			  defrx_covs_interpolation,
			  defrx_addCov,
			  defrx_matrix,
			  defrx_sigma,
			  defrx_sigmaDf,
			  defrx_nCoresRV,
			  defrx_sigmaIsChol,
                          defrx_nDisplayProgress,
			  defrx_amountUnits,
			  defrx_timeUnits, 
			  defrx_addDosing, defrx_stateTrim);
	} else if (sarg == "inits"){
	  return rxSolveC(obj,
                          CharacterVector::create("inits"),
			  R_NilValue,
                          defrx_params,
                          defrx_events,
                          as<RObject>(value), //defrx_inits,
			  R_NilValue, // scale
			  defrx_covs,
                          defrx_method,
                          defrx_transit_abs,
                          defrx_atol,
                          defrx_rtol,
                          defrx_maxsteps,
                          defrx_hmin,
                          defrx_hmax,
                          defrx_hini,
                          defrx_maxordn,
                          defrx_maxords,
                          defrx_cores,
                          defrx_covs_interpolation,
                          defrx_addCov,
                          defrx_matrix,
                          defrx_sigma,
                          defrx_sigmaDf,
                          defrx_nCoresRV,
                          defrx_sigmaIsChol,
                          defrx_nDisplayProgress,
                          defrx_amountUnits,
                          defrx_timeUnits, 
			  defrx_addDosing, defrx_stateTrim);
	} else if (sarg == "covs"){
	  return rxSolveC(obj,
                          CharacterVector::create("covs"),
			  R_NilValue,
                          defrx_params,
                          defrx_events,
                          defrx_inits,
			  R_NilValue,
                          value,// defrx_covs,
                          defrx_method,
                          defrx_transit_abs,
                          defrx_atol,
                          defrx_rtol,
                          defrx_maxsteps,
                          defrx_hmin,
                          defrx_hmax,
                          defrx_hini,
                          defrx_maxordn,
                          defrx_maxords,
                          defrx_cores,
                          defrx_covs_interpolation,
                          defrx_addCov,
                          defrx_matrix,
                          defrx_sigma,
                          defrx_sigmaDf,
                          defrx_nCoresRV,
                          defrx_sigmaIsChol,
                          defrx_nDisplayProgress,
                          defrx_amountUnits,
                          defrx_timeUnits, 
			  defrx_addDosing, defrx_stateTrim);
	} else if (sarg == "t" || sarg == "time"){
	  CharacterVector classattr = obj.attr("class");
          Environment e = as<Environment>(classattr.attr(".RxODE.env"));
          updateSolveEnvPost(e);
	  Function f = as<Function>(e["replace.sampling"]);
	  return f(value);
        } else {
	  CharacterVector classattr = obj.attr("class");
	  Environment e = as<Environment>(classattr.attr(".RxODE.env"));
          updateSolveEnvPost(e);
          if (rxIs(e["params.dat"], "NULL")) stop("Cannot update nonexistent parameters.");
          List pars = List(e["params.dat"]);
	  CharacterVector nmp = pars.names();
	  int i, n, np, nc, j;
	  np = (as<NumericVector>(pars[0])).size();
	  RObject covsR = e["covs"];
	  List covs;
	  if (!covsR.isNULL()){
	    covs = List(covsR);
	  }
	  CharacterVector nmc;
	  if (covs.hasAttribute("names")){
	    nmc = covs.names();
	    nc = (as<NumericVector>(covs[0])).size();
	  } else {
	    nc = as<int>(e["nobs"]);
	  }
	  //////////////////////////////////////////////////////////////////////////////
	  // Update Parameters by name
	  n = pars.size();
	  for (i = n; i--; ){
	    if (nmp[i] == sarg){
	      // Update solve.
	      NumericVector val = NumericVector(value);
	      if (val.size() == np){
		// Update Parameter
		pars[i] = val;
		return rxSolveC(obj,
				CharacterVector::create("params"),
				R_NilValue,
				pars, //defrx_params,
				defrx_events,
				defrx_inits,
				R_NilValue,
				defrx_covs,
				defrx_method,
				defrx_transit_abs,
				defrx_atol,
				defrx_rtol,
				defrx_maxsteps,
				defrx_hmin,
				defrx_hmax,
				defrx_hini,
				defrx_maxordn,
				defrx_maxords,
				defrx_cores,
				defrx_covs_interpolation,
				defrx_addCov,
				defrx_matrix,
				defrx_sigma,
				defrx_sigmaDf,
				defrx_nCoresRV,
				defrx_sigmaIsChol,
                                defrx_nDisplayProgress,
				defrx_amountUnits,
				defrx_timeUnits, 
				defrx_addDosing, defrx_stateTrim);
	      } else if (val.size() == nc){
		// Change parameter -> Covariate
		List newPars(pars.size()-1);
		CharacterVector newParNames(pars.size()-1);
		for (j = i; j--;){
		  newPars[j]     = pars[j];
		  newParNames[j] = nmp[j];
		}
		for (j=i+1; j < pars.size(); j++){
		  newPars[j-1]     = pars[j];
		  newParNames[j-1] = nmp[j];
		}
		newPars.attr("names") = newParNames;
		newPars.attr("class") = "data.frame";
		newPars.attr("row.names") = IntegerVector::create(NA_INTEGER,-np);
		List newCovs(covs.size()+1);
		CharacterVector newCovsNames(covs.size()+1);
		for (j = covs.size(); j--;){
		  newCovs[j]      = covs[j];
		  newCovsNames[j] = nmc[j];
		}
		j = covs.size();
		newCovs[j]      = val;
		newCovsNames[j] = nmp[i];
		newCovs.attr("names") = newCovsNames;
		newCovs.attr("class") = "data.frame";
		newCovs.attr("row.names") = IntegerVector::create(NA_INTEGER,-nc);
		return rxSolveC(obj,
				CharacterVector::create("params","covs"),
				R_NilValue,
				newPars, //defrx_params,
				defrx_events,
				defrx_inits,
				R_NilValue,
				newCovs, //defrx_covs
				defrx_method,
				defrx_transit_abs,
				defrx_atol,
				defrx_rtol,
				defrx_maxsteps,
				defrx_hmin,
				defrx_hmax,
				defrx_hini,
				defrx_maxordn,
				defrx_maxords,
				defrx_cores,
				defrx_covs_interpolation,
				defrx_addCov,
				defrx_matrix,
				defrx_sigma,
				defrx_sigmaDf,
				defrx_nCoresRV,
				defrx_sigmaIsChol,
                                defrx_nDisplayProgress,
				defrx_amountUnits,
				defrx_timeUnits, 
				defrx_addDosing, defrx_stateTrim);
	      }
	      return R_NilValue;
	    }
	  }
	  ///////////////////////////////////////////////////////////////////////////////
	  // Update Covariates by covariate name
	  n = covs.size();
	  for (i = n; i--;){
	    if (nmc[i] == sarg){
	      // Update solve.
	      NumericVector val = NumericVector(value);
	      if (val.size() == nc){
		// Update Covariate
		covs[i]=val;
		return rxSolveC(obj,
				CharacterVector::create("covs"),
				R_NilValue,
				defrx_params,
				defrx_events,
				defrx_inits,
				R_NilValue,
				covs, // defrx_covs,
				defrx_method,
				defrx_transit_abs,
				defrx_atol,
				defrx_rtol,
				defrx_maxsteps,
				defrx_hmin,
				defrx_hmax,
				defrx_hini,
				defrx_maxordn,
				defrx_maxords,
				defrx_cores,
				defrx_covs_interpolation,
				defrx_addCov,
				defrx_matrix,
				defrx_sigma,
				defrx_sigmaDf,
				defrx_nCoresRV,
				defrx_sigmaIsChol,
                                defrx_nDisplayProgress,
				defrx_amountUnits,
				defrx_timeUnits, 
				defrx_addDosing, defrx_stateTrim);
	      } else if (val.size() == np){
		// Change Covariate -> Parameter
		List newPars(pars.size()+1);
		CharacterVector newParNames(pars.size()+1);
		for (j = 0; j < pars.size(); j++){
		  newPars[j]     = pars[j];
		  newParNames[j] = nmp[j];
		}
		newPars[j]     = val;
		newParNames[j] = nmc[i];
		newPars.attr("names") = newParNames;
		newPars.attr("class") = "data.frame";
		newPars.attr("row.names") = IntegerVector::create(NA_INTEGER,-np);
		// if ()
		List newCovs(covs.size()-1);
		CharacterVector newCovsNames(covs.size()-1);
		for (j = 0; j < i; j++){
		  newCovs[j]      = covs[j];
		  newCovsNames[j] = nmc[j];
		}
		for (j=i+1; j < covs.size(); j++){
		  newCovs[j-1]      = covs[j];
		  newCovsNames[j-1] = nmc[j];
		}
		newCovs.attr("names") = newCovsNames;
		newCovs.attr("class") = "data.frame";
		newCovs.attr("row.names") = IntegerVector::create(NA_INTEGER,-nc);
		return rxSolveC(obj,
				CharacterVector::create("covs", "params"),
				R_NilValue,
				newPars,//defrx_params,
				defrx_events,
				defrx_inits,
				R_NilValue,
				newCovs, // defrx_covs,
				defrx_method,
				defrx_transit_abs,
				defrx_atol,
				defrx_rtol,
				defrx_maxsteps,
				defrx_hmin,
				defrx_hmax,
				defrx_hini,
				defrx_maxordn,
				defrx_maxords,
				defrx_cores,
				defrx_covs_interpolation,
				defrx_addCov,
				defrx_matrix,
				defrx_sigma,
				defrx_sigmaDf,
				defrx_nCoresRV,
				defrx_sigmaIsChol,
                                defrx_nDisplayProgress,
				defrx_amountUnits,
				defrx_timeUnits, 
				defrx_addDosing, defrx_stateTrim);
	      }
	    }
	  }
	  ////////////////////////////////////////////////////////////////////////////////
          // Update Initial Conditions
	  NumericVector ini = NumericVector(e["inits.dat"]);
          CharacterVector nmi = ini.names();
          n = ini.size();
          std::string cur;
	  bool doIt = false;
          for (i = 0; i < n; i++){
            cur = as<std::string>(nmi[i]) + "0";
	    if (cur == sarg){
	      doIt = true;
	    } else {
              cur = as<std::string>(nmi[i]) + ".0";
              if (cur == sarg){
                doIt = true;
	      } else {
		cur = as<std::string>(nmi[i]) + "_0";
                if (cur == sarg){
		  doIt = true;
		} else {
		  cur = as<std::string>(nmi[i]) + "(0)";
                  if (cur == sarg){
                    doIt = true;
                  } else {
		    cur = as<std::string>(nmi[i]) + "[0]";
                    if (cur == sarg){
		      doIt = true;
		    } 
		  }
		}
	      }
	    }
	    if (doIt){
	      cur=as<std::string>(nmi[i]);
              NumericVector ini = NumericVector(e["inits.dat"]);
	      double v = as<double>(value);
	      for (j = 0; j < n; j++){
		if (cur == as<std::string>(nmi[j])){
		  ini[j] = v;
		}
	      }
              return rxSolveC(obj,
			      CharacterVector::create("inits"),
			      R_NilValue,
			      defrx_params,
			      defrx_events,
			      ini,
			      R_NilValue,
			      defrx_covs,
			      defrx_method,
			      defrx_transit_abs,
			      defrx_atol,
			      defrx_rtol,
			      defrx_maxsteps,
			      defrx_hmin,
			      defrx_hmax,
			      defrx_hini,
			      defrx_maxordn,
			      defrx_maxords,
			      defrx_cores,
			      defrx_covs_interpolation,
			      defrx_addCov,
			      defrx_matrix,
			      defrx_sigma,
			      defrx_sigmaDf,
			      defrx_nCoresRV,
			      defrx_sigmaIsChol,
                              defrx_nDisplayProgress,
			      defrx_amountUnits,
			      defrx_timeUnits, 
			      defrx_addDosing, defrx_stateTrim);
	    }
	  }
	  return R_NilValue;
	}
      }
    }
  }
  return R_NilValue;
}

extern "C" SEXP rxGetModelLib(const char *s){
  std::string str(s);
  getRxModels();
  if (_rxModels.exists(str)){
    return wrap(_rxModels.get(str));
  } else {
    return R_NilValue;
  }
}

//[[Rcpp::export]]
void rxRmModelLib_(std::string str){
  getRxModels();
  if (_rxModels.exists(str)){
    List trans =(as<List>(as<List>(_rxModels[str]))["trans"]);
    std::string rxlib = as<std::string>(trans["prefix"]);
    _rxModels.remove(str);
    if (_rxModels.exists(rxlib)){
      _rxModels.remove(rxlib);
    }
  }  
}
extern "C" void rxClearFuns();
extern "C" void rxRmModelLib(const char* s){
  std::string str(s);
  rxClearFuns();
  rxRmModelLib_(str);
}

Nullable<Environment> rxRxODEenv(RObject obj){
  if (rxIs(obj, "RxODE")){
    return(as<Environment>(obj));
  } else if (rxIs(obj, "rxSolve")){
    CharacterVector cls = obj.attr("class");
    Environment e = as<Environment>(cls.attr(".RxODE.env"));
    return rxRxODEenv(as<RObject>(e["args.object"]));
  } else if (rxIs(obj, "rxModelVars")){
    List mv = as<List>(obj);
    CharacterVector trans = mv["trans"];
    getRxModels();
    std::string prefix = as<std::string>(trans["prefix"]);
    if (_rxModels.exists(prefix)){
      return as<Environment>(_rxModels[prefix]);
    } else {
      return R_NilValue;
    }
  } else {
    return rxRxODEenv(as<RObject>(rxModelVars(obj)));
  }
}

//' Get RxODE model from object
//' @param obj RxODE family of objects
//' @export
//[[Rcpp::export]]
RObject rxGetRxODE(RObject obj){
  Nullable<Environment> rxode1 = rxRxODEenv(obj);
  if (rxode1.isNull()){
    // FIXME compile if needed.
    stop("Can't figure out the RxODE object");
  } else {
    Environment e = as<Environment>(rxode1);
    e.attr("class") = "RxODE";
    return as<RObject>(e);
  }
}

//' Checks if the RxODE object was built with the current build
//'
//' @inheritParams rxModelVars
//'
//' @return boolean indicating if this was built with current RxODE
//'
//' @export
//[[Rcpp::export]]
bool rxIsCurrent(RObject obj){
  List mv = rxModelVars(obj);
  if (mv.containsElementNamed("version")){
    CharacterVector version = mv["version"];
    std::string str = __VER_md5__;
    std::string str2 = as<std::string>(version["md5"]);
    if (str2 == str){
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

extern "C" void RxODE_assign_fn_pointers_(const char *mv);

//' Assign pointer based on model variables
//' @param object RxODE family of objects
//' @export
//[[Rcpp::export]]
void rxAssignPtr(SEXP object = R_NilValue){
  List mv=rxModelVars(as<RObject>(object));
  CharacterVector trans = mv["trans"];
  RxODE_assign_fn_pointers_((as<std::string>(trans["model_vars"])).c_str());
  rxUpdateFuns(as<SEXP>(trans));
  getRxSolve_();
  // Update rxModels environment.
  getRxModels();
  
  std::string ptr = as<std::string>(trans["model_vars"]); 
  if (!_rxModels.exists(ptr)){
    _rxModels[ptr] = mv;
  } else if (!rxIsCurrent(as<RObject>(_rxModels[ptr]))) {
    _rxModels[ptr] = mv;
  }
  Nullable<Environment> e1 = rxRxODEenv(object);
  if (!e1.isNull()){
    std::string prefix = as<std::string>(trans["prefix"]);
    if (!_rxModels.exists(prefix)){
      Environment e = as<Environment>(e1);
      _rxModels[prefix] = e;
    }
  }
}

extern "C" void rxAssignPtrC(SEXP obj){
  rxAssignPtr(obj);
}

//' Get the number of cores in a system
//' @export
//[[Rcpp::export]]
IntegerVector rxCores(){
  unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
  return IntegerVector::create((int)(concurentThreadsSupported));
}

//' Return the DLL associated with the RxODE object
//'
//' This will return the dynamic load library or shared object used to
//' run the C code for RxODE.
//'
//' @param obj A RxODE family of objects or a character string of the
//'     model specification or location of a file with a model
//'     specification.
//'
//' @return a path of the library
//'
//' @keywords internal
//' @author Matthew L.Fidler
//' @export
//[[Rcpp::export]]
std::string rxDll(RObject obj){
  if (rxIs(obj,"RxODE")){
    Environment e = as<Environment>(obj);
    return as<std::string>((as<List>(e["rxDll"]))["dll"]);
  } else if (rxIs(obj,"rxSolve")) {
    CharacterVector cls = obj.attr("class");
    Environment e = as<Environment>(cls.attr(".RxODE.env"));
    return(as<std::string>(e["dll"]));
  } else if (rxIs(obj, "rxDll")){
    return as<std::string>(as<List>(obj)["dll"]);
  } else if (rxIs(obj, "character")){
    Function f = getRxFn("rxCompile.character");
    RObject newO = f(as<std::string>(obj));
    return(rxDll(newO));
  } else {
    List mv = rxModelVars(obj);
    Nullable<Environment> en = rxRxODEenv(mv);
    if (en.isNull()){
      stop("Can't figure out the DLL for this object.");
    } else {
      Environment e = as<Environment>(en);
      return as<std::string>((as<List>(e["rxDll"]))["dll"]);
    }
  }
}

//' Return the C file associated with the RxODE object
//'
//' This will return C code for generating the RxODE DLL.
//'
//' @param obj A RxODE family of objects or a character string of the
//'     model specification or location of a file with a model
//'     specification.
//'
//' @return a path of the library
//'
//' @keywords internal
//' @author Matthew L.Fidler
//' @export
//[[Rcpp::export]]
CharacterVector rxC(RObject obj){
  std::string rets;
  CharacterVector ret(1);
  if (rxIs(obj,"RxODE")){
    Environment e = as<Environment>(obj);
    rets = as<std::string>((as<List>(e["rxDll"]))["c"]);
  } else if (rxIs(obj,"rxSolve")) {
    CharacterVector cls = obj.attr("class");
    Environment e = as<Environment>(cls.attr(".RxODE.env"));
    rets = as<std::string>(e["c"]);
  } else if (rxIs(obj, "rxDll")){
    rets = as<std::string>(as<List>(obj)["c"]);
  } else if (rxIs(obj, "character")){
    Function f = getRxFn("rxCompile.character");
    RObject newO = f(as<std::string>(obj));
    rets = rxDll(newO);
  } else {
    List mv = rxModelVars(obj);
    Nullable<Environment> en = rxRxODEenv(mv);
    if (en.isNull()){
      stop("Can't figure out the DLL for this object");
    } else {
      Environment e = as<Environment>(en);
      rets = as<std::string>((as<List>(e["rxDll"]))["dll"]);
    }
  }
  ret[0] = rets;
  ret.attr("class") = "rxC";
  return ret;
}

//' Determine if the DLL associated with the RxODE object is loaded
//'
//' @param obj A RxODE family of objects 
//'
//' @return Boolean returning if the RxODE library is loaded.
//'
//' @keywords internal
//' @author Matthew L.Fidler
//' @export
//[[Rcpp::export]]
bool rxIsLoaded(RObject obj){
  if (obj.isNULL()) return false;
  Function isLoaded("is.loaded", R_BaseNamespace);
  List mv = rxModelVars(obj);
  CharacterVector trans = mv["trans"];
  std::string dydt = as<std::string>(trans["ode_solver"]);
  return as<bool>(isLoaded(dydt));
}

//' Load RxODE object
//'
//' @param obj A RxODE family of objects 
//'
//' @return Boolean returning if the RxODE library is loaded.
//'
//' @keywords internal
//' @author Matthew L.Fidler
//' @export
//[[Rcpp::export]]
bool rxDynLoad(RObject obj){
  if (!rxIsLoaded(obj)){
    std::string file = rxDll(obj);
    if (fileExists(file)){
      dynLoad(file);
    } else {
      Nullable<Environment> e1 = rxRxODEenv(obj);
      if (!e1.isNull()){
	Environment e = as<Environment>(e1);
	Function compile = as<Function>(e["compile"]);
	compile();
      }
    }
  }
  if (rxIsLoaded(obj)){
    rxAssignPtr(obj);
    return true;
  } else {
    return false;
  }
}

//' Unload RxODE object
//'
//' @param obj A RxODE family of objects 
//'
//' @return Boolean returning if the RxODE library is loaded.
//'
//' @keywords internal
//' @author Matthew L.Fidler
//' @export
//[[Rcpp::export]]
bool rxDynUnload(RObject obj){
  List mv = rxModelVars(obj);
  CharacterVector trans = mv["trans"];
  std::string ptr = as<std::string>(trans["model_vars"]);
  if (rxIsLoaded(obj)){
    Function dynUnload("dyn.unload", R_BaseNamespace);
    std::string file = rxDll(obj);
    dynUnload(file);
  } 
  rxRmModelLib_(ptr);
  return !(rxIsLoaded(obj));
}

//' Delete the DLL for the model
//'
//' This function deletes the DLL, but doesn't delete the model
//' information in the object.
//'
//' @param obj RxODE family of objects
//'
//' @return A boolean stating if the operation was successful.
//'
//' @author Matthew L.Fidler
//' @export
//[[Rcpp::export]]
bool rxDelete(RObject obj){
  std::string file = rxDll(obj);
  if (rxDynUnload(obj)){
    CharacterVector cfileV = rxC(obj);
    std::string cfile = as<std::string>(cfileV[0]);
    if (fileExists(cfile)) remove(cfile.c_str());
    if (!fileExists(file)) return true;
    if (remove(file.c_str()) == 0) return true;
  }
  return false;
}

extern "C" SEXP rxModelVarsC(char *ptr){
  // Rcout << "rxModelVars C: ";
  return rxGetFromChar(ptr, "");
}

extern "C" SEXP rxStateNames(char *ptr){
  // Rcout << "State: ";
  return rxGetFromChar(ptr, "state");
}

extern "C" SEXP rxLhsNames(char *ptr){
  // Rcout << "Lhs: ";
  return rxGetFromChar(ptr, "lhs");
}

extern "C" SEXP rxParamNames(char *ptr){
  // Rcout << "Param: ";
  return rxGetFromChar(ptr, "params");
}

extern "C" int rxIsCurrentC(SEXP obj){
  RObject robj = as<RObject>(obj);
  if (robj.isNULL()) return 0;
  bool ret = rxIsCurrent(robj);
  if (ret) return 1;
  return 0;
}

extern "C" int Rcat(char *msg){
  std::string str(msg);
  Rcpp::Function msgF("message"); 
  msgF(str, _["appendLF"]=false);
  return 1; 
}

int isRstudioI = 0;

//[[Rcpp::export]]
SEXP setRstudio(bool isRstudio=false){
  if (isRstudio) isRstudioI=1; 
  else isRstudioI=0;
  return wrap(isRstudioI);
}

extern "C" int isRstudio(){
  return isRstudioI;
}
