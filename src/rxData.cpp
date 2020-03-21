// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#define NCMT 100
// NONMEM 7.1 has a max of 50 obesrvations/individual
#define MAXIDS 500
#define NALL 500
#define NDOSES 50
//#define rxSolveT 1
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
#include "../inst/include/RxODE.h"
#include "ode.h"
#include "timsort.h"	
#ifdef rxSortStd	
#define SORT std::sort	
#else	
#define SORT gfx::timsort	
#endif
#define rxModelVars(a) rxModelVars_(a)
#define min2( a , b )  ( (a) < (b) ? (a) : (b) )
void resetSolveLinB();
using namespace Rcpp;
using namespace arma;

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

LogicalVector rxSolveFree();

List etTrans(List inData, const RObject &obj, bool addCmt=false,
	     bool dropUnits=false, bool allTimeVar=false,
	     bool keepDosingOnly=false, Nullable<LogicalVector> combineDvid=R_NilValue,
	     CharacterVector keep = CharacterVector(0));
List etImportEventTable(List inData);
RObject et_(List input, List et__);
void setEvCur(RObject cur);

RObject cvPost_(double nu, RObject omega, int n = 1, bool omegaIsChol = false, bool returnChol = false,
	       int type=1, int diagXformType=1);

int rxcEvid  = -1;
int rxcTime  = -1;
int rxcAmt   = -1;
int rxcId    = -1;
int rxcDv    = -1;
int rxcLimit = -1;
int rxcCens  = -1;
int rxcLen   = -1;
int rxcIi    = -1;
bool resetCache = true;
bool rxHasEventNames(CharacterVector &nm){
  int len = nm.size();
  bool reset  = resetCache;
  if (reset || len != rxcLen){
    reset    = resetCache;
    rxcEvid  = -1;
    rxcTime  = -1;
    rxcAmt   = -1;
    rxcId    = -1;
    rxcDv    = -1;
    rxcIi    = -1;
    rxcLimit = -1;
    rxcCens  = -1;
    rxcLen   = len;
    for (unsigned int i = len; i--;){
      if (as<std::string>(nm[i]) == "evid" ||
	  as<std::string>(nm[i]) == "EVID" ||
	  as<std::string>(nm[i]) == "Evid"){
        rxcEvid = i;
      } else if (as<std::string>(nm[i]) == "time" ||
		 as<std::string>(nm[i]) == "TIME" ||
		 as<std::string>(nm[i]) == "Time"){
        rxcTime = i;
      } else if (as<std::string>(nm[i]) == "amt" ||
		 as<std::string>(nm[i]) == "AMT" ||
		 as<std::string>(nm[i]) == "Amt"){
        rxcAmt = i;
      } else if (as<std::string>(nm[i]) == "id" ||
		 as<std::string>(nm[i]) == "ID" ||
		 as<std::string>(nm[i]) == "Id"){
        rxcId = i;
      } else if (as<std::string>(nm[i]) == "dv" ||
		 as<std::string>(nm[i]) == "DV" ||
		 as<std::string>(nm[i]) == "Dv"){
        rxcDv = i;
      } else if (as<std::string>(nm[i]) == "ii" ||
		 as<std::string>(nm[i]) == "II" ||
		 as<std::string>(nm[i]) == "Ii"){
	rxcIi = i;
      } else if (as<std::string>(nm[i]) == "cens" ||
		 as<std::string>(nm[i]) == "CENS" ||
		 as<std::string>(nm[i]) == "Cens"){
	rxcCens = i;
      } else if (as<std::string>(nm[i]) == "limit" ||
		 as<std::string>(nm[i]) == "LIMIT" ||
		 as<std::string>(nm[i]) == "Limit"){
	rxcLimit = i;
      }
    }
  }
  resetCache = true;
  if (rxcTime >= 0){
    return true;
  } else {
    return false;
  }
}

bool rxIs_list(const RObject &obj, std::string cls){
  bool hasCls = obj.hasAttribute("class");
  if (hasCls){
    CharacterVector classattr = obj.attr("class");
    bool hasDf = false;
    bool hasEt = false;
    std::string cur;
    for (unsigned int i = classattr.size(); i--; ){
      cur = as<std::string>(classattr[i]);
      if (cur == cls){
	if (cls == "rxEt"){
	  List ce = as<List>(classattr.attr(".RxODE.lst"));
	  List lobj = List(obj);
	  int nobs = ce["nobs"];
	  int ndose = ce["ndose"];
	  if (lobj.size() != 12){
	    lobj.attr("class") = CharacterVector::create("data.frame");
	    return false;
	  }
	  if ( (as<IntegerVector>(lobj[0])).size() != ndose + nobs){
	    lobj.attr("class") = CharacterVector::create("data.frame");
	    return false;
	  }
	  return true;
	} else if (cls == "rxSolve"){
	  Environment e = as<Environment>(classattr.attr(".RxODE.env"));
	  List lobj = List(obj);
	  CharacterVector cls2= CharacterVector::create("data.frame");
	  if (as<int>(e[".check.ncol"]) != lobj.size()){
	    lobj.attr("class") = cls2;
	    return false;
	  }
	  int nrow = (as<NumericVector>(lobj[0])).size();
	  if (as<int>(e[".check.nrow"]) != nrow){
	    lobj.attr("class") = cls2;
	    return false;
	  }
	  CharacterVector cn = CharacterVector(e[".check.names"]);
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
      }
    }
    if (hasDf && (cls == "rx.event" || cls == "event.data.frame")){
      if (classattr[0] == "rxEtTran"){
	rxcEvid = 2;
	rxcTime = 1;
	rxcAmt  = 3;
	rxcId   = 0;
	rxcDv   = 5;
	rxcIi   = 4;
	List e = as<List>(classattr.attr(".RxODE.lst"));
	int censAdd = as<int>(e["censAdd"]);
	int limitAdd = as<int>(e["limitAdd"]);
	if (censAdd == 1 && limitAdd == 1) {
	  rxcCens = 6;
	  rxcLimit = 7;
	} else if (censAdd == 1){
	  rxcCens = 6;
	  rxcLimit = -1;
	} else if (limitAdd == 1){
	  rxcCens  = -1;
	  rxcLimit = 6;
	} else {
	  rxcCens  = -1;
	  rxcLimit = -1;
	}
	return true;
      } else {
	// Check for event.data.frame
	CharacterVector cv =as<CharacterVector>((as<DataFrame>(obj)).names());
	return rxHasEventNames(cv);
      }
    } else if (hasEt) {
      return (cls == "rx.event");
    } else {
      return false;
    }
  } else {
    return (cls == "list");
  }
}

bool rxDropB = false;

List rxDrop(CharacterVector drop, List input, bool warnDrop) {
  rxDropB=false;
  CharacterVector inNames = input.attr("names");
  std::vector<int> keepI;
  int ndrop=0;
  for (int i = 0; i < inNames.size(); ++i) {
    std::string curName = as<std::string>(inNames[i]);
    bool dropCur = false;
    for (int j = drop.size();j--;){
      if (as<std::string>(drop[j]) == curName){
	dropCur = true;
	break;
      }
    }
    if (dropCur) ndrop++;
    else keepI.push_back(i);
    if (dropCur && !rxDropB && i < 10){
      if (curName == "time" ||
	  curName == "sim.id" ||
	  curName == "id") {
	rxDropB=true;
      }
    }
  }
  if (warnDrop && ndrop != drop.size()) {
    warning(_("column(s) in 'drop' were not in solved data"));
  }
  List ret(keepI.size());
  CharacterVector retN(keepI.size());
  for (int i = keepI.size(); i--;){
    ret[i] = input[keepI[i]];
    retN[i] = inNames[keepI[i]];
  }
  ret.attr("names") = retN;
  ret.attr("row.names")=input.attr("row.names");
  return ret;
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
  if (obj == NULL) return false;
  if (cls == "units"){
    if (obj.hasAttribute("class")){
      CharacterVector cls = obj.attr("class");
      return as<std::string>(cls[0]) == "units";
    }
  }
  if (obj == NULL) return false;
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
      if (cls == "integer") return true;
      hasCls = obj.hasAttribute("class");
      if (hasCls){
	CharacterVector classattr = obj.attr("class");
	if (as<std::string>(classattr[0]) == cls){
	  return true;
	}
	return false;
      }
      return false;
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
    return rxIs_list(obj, cls);
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
bool _mvnfast=false;

SEXP rxRmvn0(NumericMatrix& A_, arma::rowvec mu, arma::mat sigma,
	     arma::vec lower, arma::vec upper, int ncores=1, bool isChol=false,
	     double a=0.4, double tol = 2.05, double nlTol=1e-10, int nlMaxiter=100);

Function getRxFn(std::string name);
RObject rxSimSigma(const RObject &sigma,
		   const RObject &df,
		   int ncores,
		   const bool &isChol,
		   int nObs,
		   const bool checkNames = true,
		   NumericVector lowerIn =NumericVector::create(R_NegInf),
		   NumericVector upperIn = NumericVector::create(R_PosInf),
		   double a=0.4, double tol = 2.05, double nlTol=1e-10, int nlMaxiter=100){
  if (nObs < 1){
    rxSolveFree();
    stop(_("refusing to simulate %d items"),nObs);
  }
  if (rxIs(sigma, "numeric.matrix")){
    // FIXME more distributions
    NumericMatrix sigmaM(sigma);
    if (sigmaM.nrow() != sigmaM.ncol()){
      rxSolveFree();
      stop(_("matrix must be a square matrix"));
    }
    List dimnames;
    StringVector simNames;
    bool addNames = false;
    if (checkNames){
      if (!sigmaM.hasAttribute("dimnames")){
	rxSolveFree();
        stop(_("matrix must have named dimensions (try 'lotri')"));
      }
      dimnames = sigmaM.attr("dimnames");
      simNames = as<StringVector>(dimnames[1]);
      addNames=true;
    } else if (sigmaM.hasAttribute("dimnames")){
      dimnames = sigmaM.attr("dimnames");
      simNames = as<StringVector>(dimnames[1]);
      addNames = true;
    }
    NumericMatrix simMat(nObs,sigmaM.ncol());
    arma::rowvec m(sigmaM.ncol(),arma::fill::zeros);
    NumericVector lower(sigmaM.ncol());
    NumericVector upper(sigmaM.ncol());
    if (lowerIn.hasAttribute("names")){
      // Named
      CharacterVector lowN = lowerIn.attr("names");
      for (int j = simNames.size(); j--;){
	lower[j] = R_NegInf;
	for (int k = lowN.size(); k--;){
	  if (lowN[k] == simNames[j]){
	    lower[j] = lowerIn[k];
	    break;
	  }
	}
      }
    } else if (lowerIn.size() == 1){
      std::fill_n(lower.begin(), lower.size(), lowerIn[0]);
    } else if (lowerIn.size() == lower.size()){
      lower = lowerIn;
    } else {
      rxSolveFree();
      stop(_("lower bounds needs to be a named vector, a single value or exactly the same size"));
    }
    if (upperIn.hasAttribute("names")){
      // Named
      CharacterVector upN = upperIn.attr("names");
      for (int j = simNames.size(); j--;){
	upper[j] = R_NegInf;
	for (int k = upN.size(); k--;){
	  if (upN[k] == simNames[j]){
	    upper[j] = upperIn[k];
	    break;
	  }
	}
      }
    } else if (upperIn.size() == 1){
      std::fill_n(upper.begin(), upper.size(), upperIn[0]);
    } else if (!upperIn.hasAttribute("names") && upperIn.size() == upper.size()){
      upper = upperIn;
    } else {
      rxSolveFree();
      stop(_("upper bounds needs to be a named vector, a single value or exactly the same size"));
    }

    // Ncores = 1?  Should it be parallelized when it can be...?
    // Note that if so, the number of cores also affects the output.
    // while (totSim < nObs){
    if (_mvnfast){
      int totSim = 0;
      while (totSim < nObs){
	int curSimN = nObs - totSim;
	NumericMatrix simMat0(curSimN,sigmaM.ncol());
	if (df.isNULL()){
	  Function rmvn = getRxFn(".rmvn");
	  rmvn(_["n"]=curSimN, _["mu"]=m, _["sigma"]=sigmaM, _["ncores"]=ncores,
	       _["isChol"]=isChol, _["A"] = simMat0); // simMat is updated with the random deviates
	} else {
	  double df2 = as<double>(df);
	  if (R_FINITE(df2)){
	    Function rmvt = getRxFn(".rmvt");
	    rmvt(_["n"]=curSimN, _["mu"]=m, _["sigma"]=sigmaM, _["df"] = df,
		 _["ncores"]=ncores, _["isChol"]=isChol, _["A"] = simMat0);
	  } else {
	    Function rmvn = getRxFn(".rmvn");
	    rmvn(_["n"]=curSimN, _["mu"]=m, _["sigma"]=sigmaM, _["ncores"]=ncores,
		 _["isChol"]=isChol, _["A"] = simMat0);
	  }
	}
	// Reject any bad simulations if needed
	for (int i = 0; i < curSimN; i++){
	  bool goodSim = true;
	  for (int j = sigmaM.ncol(); j--;){
	    double cur = simMat0(i, j);
	    if (cur <= lower[j] || cur >= upper[j]){
	      goodSim=false;
	      break;
	    }
	  }
	  if (goodSim){
	    simMat(totSim, _) = simMat0(i, _);
	    totSim++;
	  }
	}
      }
    } else if (df.isNULL()){
      rxRmvn0(simMat, m, as<arma::mat>(sigmaM), as<arma::vec>(lower),
	      as<arma::vec>(upper), ncores, isChol, a, tol, nlTol, nlMaxiter);
	// Function rmvn = as<Function>(mvnfast["rmvn"]);
	// rmvn(_["n"]=curSimN, _["mu"]=m, _["sigma"]=sigmaM, _["ncores"]=ncores,
	//      _["isChol"]=isChol, _["A"] = simMat0); // simMat is updated with the random deviates
    } else {
      double df2 = as<double>(df);
      if (R_FINITE(df2)){
	stop(_("t distribution not yet supported"));
	  // Function rmvt = as<Function>(mvnfast["rmvt"]);
	// rmvt(_["n"]=curSimN, _["mu"]=m, _["sigma"]=sigmaM, _["df"] = df,
	//      _["ncores"]=ncores, _["isChol"]=isChol, _["A"] = simMat0);
      } else {
	// Function rmvn = as<Function>(mvnfast["rmvn"]);
	// rmvn(_["n"]=curSimN, _["mu"]=m, _["sigma"]=sigmaM, _["ncores"]=ncores,
	//      _["isChol"]=isChol, _["A"] = simMat0);
	rxRmvn0(simMat, m, as<arma::mat>(sigmaM), as<arma::vec>(lower),
	      as<arma::vec>(upper), ncores, isChol, a, tol, nlTol, nlMaxiter);
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

Environment RxODEenv(){
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  _RxODE = loadNamespace("RxODE");
  _RxODE_found = true;
  return _RxODE;
}
// Export for C.
//[[Rcpp::export]]
Function getRxFn(std::string name){
  Environment rx = RxODEenv();
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

List rxModelVars_(const RObject &obj); // model variables section

List rxModelVars_RxODE(const RObject &obj){
  Environment e = as<Environment>(obj);
  List rxDll = as<List>(e["rxDll"]);
  List ret = as<List>(rxDll["modVars"]);
  RObject pkgR = e["package"];
  if (rxIs(pkgR,"NULL")){
    return ret;
  } else {
    bool isV;
    Function isValid = e["isValid"];
    isV = as<bool>(isValid());
    if (isV){
      return ret;
    } else {
      std::string sobj = as<std::string>(e["modName"]);
      if (sobj.find("_new")!=std::string::npos){
	return ret;
      }
      Function pkgLoaded = getRxFn(".rxPkgLoaded");
      isV = as<bool>(pkgLoaded(pkgR));
      if (isV){
	Environment ns;
	if (as<std::string>(pkgR) == "RxODE"){
	  return ret;
	} else {
	  ns=loadNamespace(pkgR);
	}
	if (ns.exists(".rxUpdated")){
	  Environment rxU = ns[".rxUpdated"];
	  if (rxU.exists(sobj)){
	    e = rxU[sobj];
	    rxDll = e["rxDll"];
	    ret = rxDll["modVars"];
	    return ret;
	  } else {
	    Function rxode = getRxFn("RxODE");
	    e["modName"] = sobj + "_new"; // For RxODE parsing add "new"
	    Environment newRx = rxode(e);
	    e["modName"] = sobj; // Put back
	    rxU[sobj] = newRx;
	    return rxModelVars_(as<RObject>(newRx));
	  }
	} else {
	  return ret;
	}
      } else {
	return ret;
      }
    }
  }
}
List rxModelVars_blank(){
  List ret(21);
  CharacterVector retN(21);
  ret[0]  = CharacterVector::create(); // params
  retN[0] = "params";
  ret[1]  = CharacterVector::create(); // lhs
  retN[1] = "lhs";
  ret[2]  = CharacterVector::create(); // state
  retN[2] = "state";
  ret[3]  = CharacterVector::create(); // trans
  retN[3] = "trans";
  ret[4]  = CharacterVector::create(_["normModel"] = ""); // model
  retN[4] = "model";
  ret[5]  = NumericVector::create(); // ini
  retN[5] = "ini";
  ret[6]  = LogicalVector::create(false); // podo
  retN[6] = "podo";
  ret[7]  = LogicalVector::create(); // dfdy
  retN[7] = "dfdy";
  ret[8]  = LogicalVector::create(); // sens
  retN[8] = "sens";
  ret[9] = IntegerVector::create(); // state.ignore
  retN[9] = "state.ignore";
  ret[10] = CharacterVector::create(); // version
  retN[10] = "version";
  ret[11] = CharacterVector::create(); // normal.state
  retN[11] = "normal.state";
  ret[12] = IntegerVector::create(0); // need sort
  retN[12] = "needSort";
  ret[13] = IntegerVector::create(0); // nMtime
  retN[13] = "nMtime";
  ret[14] = IntegerVector::create(0); // extraCmt
  retN[14] = "extraCmt";
  ret[15] = CharacterVector::create(); // stateExtra
  retN[15] = "stateExtra";
  ret[16] = IntegerVector::create(); // dvid
  retN[16] = "dvid";
  ret[17] = List::create();
  retN[17] = "indLin";  
  ret[18] = IntegerVector::create(0); // timeId
  retN[18] = "timeId";
  ret[19] =CharacterVector::create(_["file_md5"] = "", _["parsed_md5"] = ""); // md5
  retN[19] = "md5";
  ret.attr("names") = retN;
  ret.attr("class") = "rxModelVars";
  return ret;
}

List rxModelVars_character(const RObject &obj){
  CharacterVector modList = as<CharacterVector>(obj);
  if (modList.size() == 1){
    std::string sobj =as<std::string>(obj);
    if (sobj == ""){
      // Blank RxODE model
      return rxModelVars_blank();
    } else if (fileExists(sobj)){
      // From file
      Function f = getRxFn(".rxModelVarsCharacter");
      return f(obj);
    } else if ((sobj.find("=") == std::string::npos) &&
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
  if (obj == NULL) {
    return rxModelVars_blank();
  }
  // fileExists(const std::string& name)
  Function f = getRxFn(".rxModelVarsCharacter");
  return f(obj);
}

List rxModelVars_list(const RObject &obj){
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
    rxSolveFree();
    stop(_("cannot figure out the model variables"));
}

// [[Rcpp::export]]
List rxModelVars_(const RObject &obj){
  if (obj == NULL) return rxModelVars_blank();
  getRxModels();
  if (rxIs(obj, "rxModelVars")){
    List ret(obj);
    return ret;
  } else if (rxIs(obj,"RxODE")) {
    return rxModelVars_RxODE(obj);
  } else if (rxIs(obj, "rxS")){
    Environment e = as<Environment>(obj);
    List ret = as<List>(e["..mv"]);
    return ret;
  } else if (rxIs(obj,"rxSolve")){
    CharacterVector cls = obj.attr("class");
    Environment e = as<Environment>(cls.attr(".RxODE.env"));
    return  rxModelVars_(as<RObject>(e[".args.object"]));
  } else if (rxIs(obj,"rxDll")){
    List lobj = (as<List>(obj))["modVars"];
    return lobj;
  } else if (rxIs(obj, "environment")){
    Environment e = as<Environment>(obj);
    if (e.exists(".args.object")){
      return rxModelVars_(e[".args.object"]);
    } else {
      CharacterVector cls = obj.attr("class");
      int i = 0;
      Rprintf(_("class:\t"));
      for (i = 0; i < cls.size(); i++){
        Rprintf("%s\t", (as<std::string>(cls[i])).c_str());
      }
      Rprintf("\n");
      rxSolveFree();
      stop(_("need an RxODE-type object to extract model variables"));
    }
  } else if (rxIs(obj, "character")){
    return rxModelVars_character(obj);
  } else if (rxIs(obj,"list")){
    return rxModelVars_list(obj);
  } else if (rxIs(obj,"NULL")) {
    rxSolveFree();
    stop(_("a NULL object does not have any RxODE model variables"));
  } else {
    CharacterVector cls = obj.attr("class");
    int i = 0;
    Rprintf(_("class:\t"));
    for (i = 0; i < cls.size(); i++){
      Rprintf("%s\t", (as<std::string>(cls[i])).c_str());
    }
    Rprintf("\n");
    rxSolveFree();
    stop(_("need an RxODE-type object to extract model variables"));
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
      rxSolveFree();
      stop(_("can only lookup one state at a time"));
    }
    if (states.size() == 1){
      warning(_("only one state variable should be input"));
    }
    IntegerVector ret(1);
    for (int i = 0; i < states.size(); i++){
      if (states[i] == lookup[0]){
	ret[0] = i+1;
	return ret;
      }
    }
    rxSolveFree();
    stop(_("cannot locate compartment \"%s\""),as<std::string>(lookup[0]).c_str());
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
	  rxSolveFree();
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
    for (int j=inits.size(); j--;){
      isState=false;
      for (int k=state.size(); k--;){
	if (nms[j] == state[k]){
	  isState=true;
	  break;
	}
      }
      if (as<std::string>(nms[j]) != ""){
	ret += as<std::string>(nms[j]);
	if (isState) ret += "(0)";
	ret += "=" + std::to_string(inits[j]) + ";\n";
      }
    }
    return wrap(ret);
  } else if (vec.isNULL()){
    return wrap(rxInits0(obj, R_NilValue, req, defaultValue, noerror,noini));
  } else if (rxIs(vec, "list")){
    List vecL = as<List>(vec);
    Function unlist("unlist", R_BaseNamespace);
    NumericVector vec2 = as<NumericVector>(unlist(vec));
    if (vec2.size() != vecL.size()){
      rxSolveFree();
      stop(_("only one estimate per named list item; use 'list(x=1)' instead of 'list(x=1:2)'"));
    }
    return wrap(rxInits0(obj, vec2, req, defaultValue, noerror,noini));
  } else if (rxIs(vec, "integer") || rxIs(vec, "numeric")){
    return wrap(rxInits0(obj, as<NumericVector>(vec), req, defaultValue, noerror,noini));
  } else {
    rxSolveFree();
    stop(_("incompatible initial estimate"));
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
//' 
//' @param inits A numeric vector of initial conditions.
//' 
//' @param extraArgs A list of extra args to parse for initial conditions.
//' 
//' @author Matthew L. Fidler
//' 
//' @keywords internal
//' 
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
	  rxSolveFree();
	  stop(_("trying to scale the same compartment by 'scale=c(%s=%f,...)' and 'S%d=%f' choose one"),
	       (as<std::string>(state[i])).c_str(), ret[i], i+1,as<double>(xtra[i]));
	}
      } else {
	cur = "s" + std::to_string(i+1);
        if (xtra.containsElementNamed(cur.c_str())){
          if (ret[i] == 1.0){
            ret[i] = as<double>(xtra[cur]);
	    usedS++;
          } else {
	    rxSolveFree();
            stop(_("trying to scale the same compartment by 'scale=c(%s=%f,...)' and 's%d=%f' choose one"),
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
      warning(_("scaled a compartment that is not defined by the RxODE model"));
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
      rxSolveFree();
      stop(_("'%s' is not compatible with 'params'"), thetaTxt.c_str());
    }
  } else if (!theta.isNULL()){
    rxSolveFree();
    stop(_("'%s' is not compatible with 'params'"), thetaTxt.c_str());
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
	rxSolveFree();
        stop(_("'theta' is not compatible with 'params'"));
      }
    } else if (!theta.isNULL()){
      rxSolveFree();
      stop(_("'theta' is not compatible with 'params'"));
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
	rxSolveFree();
        stop(_("'eta' is not compatible with 'params'"));
      }
    } else if (!eta.isNULL()){
      rxSolveFree();
      stop(_("'eta' is not compatible with 'params'"));
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
      rxSolveFree();
      stop(_("parameter names must be specified"));
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
  int *gon;
  double *gsolve;
  double *gadvan;
  double *gInfusionRate;
  double *gall_times;
  int *gix;
  double *gdv;
  double *glimit;
  int *gcens;
  double *gamt;
  double *gii;
  double *glhs;
  double *gcov;
  double *ginits;
  double *gscale;
  double *gatol2;
  double *grtol2;
  double *gssAtol;
  double *gssRtol;
  double *gpars;
  //ints
  int *gevid;
  int *gBadDose;
  int *grc;
  int *gidose;
  int *gpar_cov;
  
  int *gParPos;
  int *gParPos2;
  
  int *gsvar;
  int *gsiV;
  //
  int *slvr_counter;
  int *dadt_counter;
  int *jac_counter;
  double *gmtime;
} rx_globals;

rx_globals _globals;

double maxAtolRtolFactor = 0.1;

//[[Rcpp::export]]
void atolRtolFactor_(double factor){
  rx_solve* rx = getRxSolve_();
  rx_solving_options* op = rx->op;
  for (int i = op->neq;i--;){
    _globals.grtol2[i] = min2(_globals.grtol2[i]*factor, maxAtolRtolFactor);
    _globals.gatol2[i] = min2(_globals.gatol2[i]*factor, maxAtolRtolFactor);
  }
  op->ATOL = min2(op->ATOL*factor, maxAtolRtolFactor);
  op->RTOL = min2(op->RTOL*factor, maxAtolRtolFactor);
}

extern "C" double * getAol(int n, double atol){
  return _globals.gatol2;
}

extern "C" double * getRol(int n, double rtol){
  return _globals.grtol2;
}
// This sets up the constant covariates if ev1 is rxEtTran
static inline void gparsCovSetupConstant(RObject &ev1, int npars){
  if (rxIs(ev1, "rxEtTran")) {
    rx_solve* rx = getRxSolve_();
    CharacterVector tmpCls = ev1.attr("class");
    List envCls = tmpCls.attr(".RxODE.lst");
    NumericMatrix iniPars = envCls["pars"];
    // Copy the pre-filled covariates into the parameter values.
    for (int j = rx->nsim;j--;){
      std::copy(iniPars.begin(), iniPars.end(), &_globals.gpars[0]+rx->nsub*npars*j);
    }
    IntegerVector parPos = envCls["covParPos0"];
    std::copy(parPos.begin(),parPos.end(), &_globals.gParPos2[0]);
    rx->nCov0 = parPos.size();
    rx->cov0 = _globals.gParPos2;
  }
}
void gparsCovSetup(int npars, int nPopPar, RObject ev1,rx_solve* rx){
  if (_globals.gpars != NULL) free(_globals.gpars);
  _globals.gpars = (double*)calloc(npars*nPopPar, sizeof(double));
  if (_globals.gpars == NULL){
    rxSolveFree();
    stop(_("could not allocate memory for solving parameters"));
  }
  // Fill the parameters with NA.
  std::fill_n(&_globals.gpars[0], npars*nPopPar, NA_REAL);
  gparsCovSetupConstant(ev1, npars);
}

double *_rxGetErrs = NULL;
void rxFreeErrs(){
  free(_rxGetErrs);
  _rxGetErrs=NULL;
}

extern "C" void gFree(){
  if (_globals.gsolve != NULL) free(_globals.gsolve);
  _globals.gsolve=NULL;
  if (_globals.gon != NULL) free(_globals.gon);
  _globals.gon=NULL;
  if (_globals.gpars != NULL) free(_globals.gpars);
  _globals.gpars=NULL;
  if (_globals.gParPos != NULL) free(_globals.gParPos);
  _globals.gParPos=NULL;
  if (_globals.gevid != NULL) free(_globals.gevid);
  _globals.gevid=NULL;
  if (_globals.gall_times != NULL) free(_globals.gall_times);
  _globals.gall_times=NULL;
  if (_globals.gcov != NULL) free(_globals.gcov);
  _globals.gcov=NULL;
  if (_rxGetErrs != NULL) free(_rxGetErrs);
  _rxGetErrs=NULL;
}

extern "C" double *rxGetErrs(){
  getRxModels();
  if (_rxModels.exists(".sigma")){
    NumericMatrix sigma = _rxModels[".sigma"];
    if (_rxGetErrs == NULL){
      if (_rxGetErrs != NULL) free(_rxGetErrs);
      _rxGetErrs = (double*)calloc(sigma.ncol()*sigma.nrow(), sizeof(double));
      if (_rxGetErrs == NULL) {
	rxSolveFree();
	stop(_("memory for residual errors could not be allocated"));
      }
    }
    else {
      double *tmpErr = (double*)realloc(_rxGetErrs, sigma.ncol()*sigma.nrow()*sizeof(double));
      if (tmpErr == NULL){
	rxSolveFree();
	stop(_("cannot allocate memory to simulate the residuals"));
      }
      _rxGetErrs  = tmpErr;
    }
    std::copy(sigma.begin(),sigma.end(),&_rxGetErrs[0]);
    return _rxGetErrs;
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

void rxSimTheta(CharacterVector &thetaN,
		CharacterVector& parN,
		IntegerVector &thetaPar,
		NumericMatrix &thetaM,
		bool &simTheta,
		const Nullable<NumericMatrix> &thetaMat = R_NilValue,
		const NumericVector &thetaLower = NumericVector::create(R_NegInf),
		const NumericVector &thetaUpper = NumericVector::create(R_PosInf),
		const Nullable<NumericVector> &thetaDf  = R_NilValue,
		const bool &thetaIsChol = false,
		int nStud = 1,
		int nCoresRV = 1){
  int i, j;
  if (!thetaMat.isNull() && nStud > 1){
    thetaM = as<NumericMatrix>(thetaMat);
    if (!thetaM.hasAttribute("dimnames")){
      rxSolveFree();
      stop(_("'thetaMat' must be a named matrix"));
    }
    if (!thetaIsChol){
      arma::mat tmpM = as<arma::mat>(thetaMat);
      if (!tmpM.is_sympd()){
	rxSolveFree();
	stop(_("'thetaMat' must be symmetric"));
      }
    }
    thetaM = as<NumericMatrix>(rxSimSigma(as<RObject>(thetaMat), as<RObject>(thetaDf),
					  nCoresRV, thetaIsChol, nStud, true,
					  thetaLower, thetaUpper));
    thetaN = as<CharacterVector>((as<List>(thetaM.attr("dimnames")))[1]);
    // Now put in theta parameter position; This will be used for
    // alternate types of simulation strategies A position of -1 means
    // that the parameter is not defined in the theta matrix.
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
    warning(_("'thetaMat' is ignored since nStud <= 1"));
  }
}


void rxSimOmega(bool &simOmega,
		bool &omegaSep,
		NumericMatrix &omegaM,
		CharacterVector &omegaN,
		NumericMatrix &omegaMC,
		List &omegaList,
		CharacterVector &thetaN,
		NumericMatrix &thetaM,
		std::string omegatxt="omega",
		const RObject &omega= R_NilValue,
		const Nullable<NumericVector> &omegaDf= R_NilValue,
		const NumericVector &omegaLower = NumericVector::create(R_NegInf),
		const NumericVector &omegaUpper = NumericVector::create(R_PosInf),
		const bool &omegaIsChol = false,
		std::string omegaSeparation= "auto",//("lkj", "separation")
		const int omegaXform = 1,
		double dfSub = 0,
		int nStud = 1,
		int nSub = 1){
  int j;
  if (rxIs(omega, "NULL")){
  } else if (rxIs(omega,"character")){
    // Create a matrix in order of the names.
    omegaN = as<CharacterVector>(omega);
    omegaSep=true;
    omegaMC = NumericMatrix(nStud, omegaN.size());
    j=0;
    for (int i = 0; i < omegaN.size(); i++){
      for (j = 0; j < thetaN.size(); j++){
	if (omegaN[i] == thetaN[j]){
	  omegaMC(_,i) = thetaM(_,j);
	  break;
	}
      }
      if (j == thetaN.size()){
	stop(_("parameter '%s' was not simulated in 'thetaMat'"), (as<std::string>(omegaN[i])).c_str());
      }      
    }
    simOmega = true;
  } else if (rxIs(omega,"matrix") && nSub > 1){
    simOmega = true;
    omegaM = as<NumericMatrix>(omega);
    if (!omegaM.hasAttribute("dimnames")){
      rxSolveFree();
      stop(_("'%s' must be a named matrix"),omegatxt.c_str());
    }
    if (omegaIsChol){
      omegaMC = omegaM;
    } else {
      arma::mat tmpM = as<arma::mat>(omegaM);
      if (!tmpM.is_sympd()){
	rxSolveFree();
	stop(_("'%s' must be symmetric"),omegatxt.c_str());
      }
      omegaMC = wrap(arma::chol(as<arma::mat>(omegaM)));
    }
    omegaN = as<CharacterVector>((as<List>(omegaM.attr("dimnames")))[1]);
  }
  if (nStud > 1){
    if (dfSub > 0 && simOmega) {
      if (omegaSep){
	int defaultType = 2;
	if (omegaSeparation == "auto"){
	  if (omegaN.size() >= 10){
	    defaultType = 3;
	  }
	} else if (omegaSeparation == "separation") {
	  defaultType = 3;
	}
	omegaList = cvPost_(dfSub, as<RObject>(omegaMC), 1,  false, false, defaultType,
			    omegaXform);
      } else {
	omegaList = cvPost_(dfSub, as<RObject>(omegaMC), nStud,  true, false);
      }
    }
  }
}


//' Simulate Parameters from a Theta/Omega specification
//'
//' @param params Named Vector of RxODE model parameters
//'
//' @param thetaMat Named theta matrix.
//' 
//' @param thetaLower Lower bounds for simulated population parameter variability (by default -Inf)
//'
//' @param thetaUpper Upper bounds for simulated population unexplained variability (by default Inf)
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
//' @param omegaLower Lower bounds for simulated ETAs (by default -Inf)
//'
//' @param omegaUpper Upper bounds for simulated ETAs (by default Inf)
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
//' @param omegaSeparation @template separation
//'
//' @param nStud Number virtual studies to characterize uncertainty in estimated 
//'        parameters.
//'
//' @param nObs Number of observations to simulate (with \code{sigma} matrix)
//'
//' @param sigma Matrix for residual variation.  Adds a "NA" value for each of the 
//'     individual parameters, residuals are updated after solve is completed.
//'
//' @param sigmaSeparation @template separation
//'
//' @param sigmaLower Lower bounds for simulated unexplained variability (by default -Inf)
//'
//' @param sigmaUpper Upper bounds for simulated unexplained variability (by default Inf)
//'
//' @inheritParams rxSolve
//'
//' @param dfSub Degrees of freedom to sample the between subject variability matrix from the 
//'        inverse Wishart distribution (scaled) or scaled inverse chi squared distribution.
//'
//' @param dfObs Degrees of freedom to sample the unexplained variability matrix from the 
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
                     const RObject &omega= R_NilValue,
                     const Nullable<NumericVector> &omegaDf= R_NilValue,
		     const NumericVector &omegaLower = NumericVector::create(R_NegInf),
		     const NumericVector &omegaUpper = NumericVector::create(R_PosInf),
                     const bool &omegaIsChol = false,
		     std::string omegaSeparation= "auto",//("lkj", "separation")
		     const int omegaXform = 1,
                     int nSub = 1,
                     const Nullable<NumericMatrix> &thetaMat = R_NilValue,
		     const NumericVector &thetaLower = NumericVector::create(R_NegInf),
		     const NumericVector &thetaUpper = NumericVector::create(R_PosInf),
                     const Nullable<NumericVector> &thetaDf  = R_NilValue,
                     const bool &thetaIsChol = false,
                     int nStud = 1,
                     const RObject sigma = R_NilValue,
		     const NumericVector &sigmaLower = NumericVector::create(R_NegInf),
		     const NumericVector &sigmaUpper = NumericVector::create(R_PosInf),
                     const Nullable<NumericVector> &sigmaDf= R_NilValue,
                     const bool &sigmaIsChol = false,
		     std::string sigmaSeparation= "auto",//("lkj", "separation")
		     const int sigmaXform = 1,
                     int nCoresRV = 1,
                     int nObs = 1,
                     double dfSub = 0,
                     double dfObs = 0,
		     bool simSubjects=true){
  NumericVector par;
  if (params.isNull()){
    rxSolveFree();
    stop(_("requires 'params'"));
  } else {
    par = NumericVector(params);
    if (!par.hasAttribute("names")){
      rxSolveFree();
      stop(_("'params' must be a named vector"));
    }
  }
  NumericMatrix thetaM;
  CharacterVector thetaN;
  bool simTheta = false;
  CharacterVector parN = CharacterVector(par.attr("names"));
  IntegerVector thetaPar(parN.size());
  int i, j, k;
  rxSimTheta(thetaN, parN, thetaPar, thetaM, simTheta,
	     thetaMat, thetaLower, thetaUpper, thetaDf,
	     thetaIsChol, nStud, nCoresRV);
  bool simOmega = false;
  bool omegaSep=false;
  NumericMatrix omegaM;
  CharacterVector omegaN;
  NumericMatrix omegaMC;
  List omegaList;

  if (nSub > 1 && rxIs(omega, "NULL")){
    // rxSolveFree();
    warning(_("multi-subject simulation without without 'omega'"));
  }    
  rxSimOmega(simOmega, omegaSep, omegaM, omegaN, omegaMC,
	     omegaList, thetaN, thetaM, "omega", omega, omegaDf,
	     omegaLower, omegaUpper, omegaIsChol,
	     omegaSeparation, omegaXform, dfSub, nStud, nSub);
  
  bool simSigma = false;
  bool sigmaSep=false;
  NumericMatrix sigmaM;
  CharacterVector sigmaN;
  NumericMatrix sigmaMC;
  List sigmaList;
  
  rxSimOmega(simSigma, sigmaSep, sigmaM, sigmaN, sigmaMC,
	     sigmaList, thetaN, thetaM, "sigma", sigma, sigmaDf,
	     sigmaLower, sigmaUpper, sigmaIsChol,
	     sigmaSeparation, sigmaXform, dfObs, nStud, nObs);
  
  // Now create data frame of parameter values
  List iovList;
  
  int pcol = par.size();
  int ocol = 0;
  int ncol = pcol;
  if (simOmega){
    ocol = omegaMC.ncol();
    ncol += ocol;
  }
  int scol = 0;
  if (simSigma){
    scol = sigmaMC.ncol();
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
        // nm = ret0[j]; // parameter columnem
        nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(omegaList[i]), as<RObject>(omegaDf), nCoresRV, false, nSub,
					   false, omegaLower, omegaUpper));
      } else {
        nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(omegaMC), as<RObject>(omegaDf), nCoresRV, true, nSub,
					   false, omegaLower, omegaUpper));
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
          nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(sigmaList[i]), as<RObject>(sigmaDf), nCoresRV, false, nObs*nSub,
					     false, sigmaLower, sigmaUpper));
        } else {
          nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(sigmaMC), as<RObject>(sigmaDf), nCoresRV, true, nObs*nSub,
					     false, sigmaLower, sigmaUpper));
        }
        for (j = 0; j < scol; j++){
          for (k = 0; k < nObs*nSub; k++){
            // ret1 = NumericMatrix(nObs*nStud, scol);
            ret1(nObs*nSub*i+k, j) = nm1(k, j);
          }
        }
      } else {
        if (dfObs > 0  && nStud > 1){
          nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(sigmaList[i]), as<RObject>(sigmaDf), nCoresRV, false, nObs,
					     false, sigmaLower, sigmaUpper));
        } else {
          nm1 = as<NumericMatrix>(rxSimSigma(as<RObject>(sigmaMC), as<RObject>(sigmaDf), nCoresRV, true, nObs,
					     false, sigmaLower, sigmaUpper));
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

List getEtRxsolve(Environment e){
  if (!e.exists(".et")){
    RObject eventso = e[".args.events"];
    List emptyLst(0);
    RObject et = et_(emptyLst, emptyLst);
    setEvCur(et);
    et_(List::create(_["data"] = eventso), List::create("import"));
    e[".et"] = et;
    Function parse2("parse", R_BaseNamespace);
    Function eval2("eval", R_BaseNamespace);
    // eventTable style methods
    e["get.EventTable"] = eval2(_["expr"]   = parse2(_["text"]="function() .et$get.EventTable()"),
				_["envir"]  = e);
    e["get.obs.rec"] = eval2(_["expr"]   = parse2(_["text"]="function() .et$get.obs.rec()"),
			     _["envir"]  = e);
    e["get.nobs"] = eval2(_["expr"]   = parse2(_["text"]="function() .et$get.nobs()"),
			  _["envir"]  = e);
    e["add.dosing"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {.newEt <- .et; .newEt$add.dosing(...); invisible(rxSolve(.args.object,events=.newEt, updateObject=TRUE))}"),
			    _["envir"]  = e);
    e["clear.dosing"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {.newEt <- .et; .newEt$clear.dosing(...); invisible(rxSolve(.args.object,events=.newEt, updateObject=TRUE))}"),
			      _["envir"]  = e);
    e["get.dosing"] = eval2(_["expr"]   = parse2(_["text"]="function() .et$get.dosing()"),
			    _["envir"]  = e);

    e["add.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {.newEt <- .et; .newEt$add.sampling(...); invisible(rxSolve(.args.object,events=.et,updateObject=TRUE))}"),
			      _["envir"]  = e);
      
    e["clear.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {.newEt <- .et; .newEt$clear.sampling(...); invisible(rxSolve(.args.object,events=.newEt,updateObject=TRUE))}"),
				_["envir"]  = e);

    e[".replace.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {.newEt <- .et; .newEt$clear.sampling(); .newEt$add.sampling(...); invisible(rxSolve(.args.object,events=.newEt,updateObject=TRUE))}"),
				  _["envir"]  = e);

    e["get.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function() .et$get.sampling()"),
			      _["envir"]  = e);
      
    e["get.units"] = eval2(_["expr"]   = parse2(_["text"]="function() .et$get.units()"),
			   _["envir"]  = e);

    e["import.EventTable"] = eval2(_["expr"]   = parse2(_["text"]="function(imp) {.et <- as.et(imp); invisible(rxSolve(.args.object,events=.et,updateObject=TRUE))}"),
				   _["envir"]  = e);
      
    // Note event.copy doesn't really make sense...?  The create.eventTable does basically the same thing.
  }
  return e[".et"];
}

// This updates the evironment post solve after running.  Defers some
// computational cost until the rxSolve is looked at by the user.
void updateSolveEnvPost(Environment e){
  if (!e.exists(".params.dat")){
    List mv = rxModelVars(as<RObject>(e));
    NumericVector mvIni = mv["ini"];
    CharacterVector pars = mv["params"];
    RObject parso = e[".args.params"];
    IntegerVector ppos = e[".par.pos"];
    bool IsIni = e[".par.pos.ini"];
    CharacterVector idLevels = as<CharacterVector>(e[".idLevels"]);
    if (rxIs(parso, "numeric") || rxIs(parso, "integer") ||
	rxIs(parso, "NULL")){
      double *tmp=(double*)calloc(ppos.size(),sizeof(double));
      if (tmp == NULL){
	rxSolveFree();
	stop(_("ran out of memory during 'updateSolveEnvPost'"));
      }
      // NumericVector   prs(ppos.size()-nrm);
      // CharacterVector prsn(ppos.size()-nrm+1);
      NumericVector parNumeric;
      if (!rxIs(parso, "NULL")){
        parNumeric= as<NumericVector>(parso);
      }
      unsigned int i, j=0;
      for (i = 0; i < ppos.size(); i++){
	if (ppos[i] > 0){ // User specified parameter
	  if (IsIni){
	    tmp[j] = mvIni[ppos[i]-1];
	  } else {
	    tmp[j] = parNumeric[ppos[i]-1];
	  }
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
      free(tmp);
      j=0;
      for (i = 0; i < ppos.size(); i++){
        if (ppos[i] != 0){ // User specified parameter
          prsn[j] = pars[i];
          j++;
        }
      }
      prs.names() = prsn;
      e[".params.single"] = prs;
      List pd(prs.size());
      for (unsigned int j = prs.size();j--;){
	pd[j] = NumericVector::create(prs[j]);
      }
      pd.names() = prsn;
      pd.attr("class") = "data.frame";
      pd.attr("row.names") = IntegerVector::create(NA_INTEGER,-1);
      e[".params.dat"] = pd;
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
	if (idLevels.size() > 0){	
	  tmp.attr("class") = CharacterVector::create("factor");	
	  tmp.attr("levels") = idLevels;	
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
      e[".params.dat"] = prsl;
      if (parsdf.nrow() == 1){
	NumericVector prsnv(prsl.size());
	for (j = prsl.size(); j--;){
	  prsnv[j] = (as<NumericVector>(prsl[j]))[0];
	}
	prsnv.names() = prsn;
        e[".params.single"] = prsnv;
      } else {
	e[".params.single"] = R_NilValue;
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
  getEtRxsolve(e);
}

extern "C" void rxOptionsFree();
extern "C" void rxOptionsIni();
extern "C" void rxOptionsIniEnsure(int mx);
extern "C" void parseFree(int last);
extern "C" void rxClearFuns();
extern "C" void rxFreeLast();	
//' Free the C solving/parsing information.
//'
//' Take the ODE C system and free it.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
LogicalVector rxSolveFree(){
  rxOptionsFree(); // f77 losda free
  rxOptionsIni();// realloc f77 lsoda cache
  parseFree(0); //free parser
  rxClearFuns(); // Assign all the global ODE solving functions to NULL pointers
  gFree();// Frees all the global pointers
  return LogicalVector::create(true);
}

extern "C" void RxODE_assign_fn_pointers(SEXP);

List keepIcov;
List keepFcov;
extern void setFkeep(List keep){
  keepFcov=keep;
}

extern "C" double get_ikeep(int col, int id){
  NumericVector nv = as<NumericVector>(keepIcov[col]);
  return nv[id];
}

extern "C" double get_fkeep(int col, int id){
  NumericVector nv = as<NumericVector>(keepFcov[col]);
  return nv[id];
}
extern "C" SEXP get_ikeepn(){
  return as<SEXP>(keepIcov.attr("names"));
}

extern "C" SEXP get_fkeepn(){
  return as<SEXP>(keepFcov.attr("names"));
}
typedef struct{
  bool updateObject;
  bool isRxSolve;
  bool isEnvironment;
  bool hasCmt;
  List mv;
  Nullable<LogicalVector> addDosing;
  RObject timeUnitsU;
  bool addTimeUnits;
  List covUnits;
  RObject par1;
  bool usePar1 = false;
  bool swappedEvents = false;
  RObject par1ini;
  NumericVector initsC;
  int nSize;
  bool fromIni = false;
  IntegerVector eGparPos;
  CharacterVector sigmaN;
  NumericVector parNumeric;
  DataFrame parDf;
  NumericMatrix parMat;
  int parType = 1;
  int nPopPar = 1;
  CharacterVector nmP;
  int npars;
  NumericVector mvIni;
  int nsvar;
  bool idFactor;
  bool labelID;
  bool warnIdSort;
  CharacterVector idLevels;
} rxSolve_t;

SEXP rxSolve_(const RObject &obj, const List &rxControl, const Nullable<CharacterVector> &specParams,
	      const Nullable<List> &extraArgs, const RObject &params, const RObject &events,
	      const RObject &inits, const int setupOnly);


// This will update the rxSolve object with any new features
//
static inline SEXP rxSolve_update(const RObject &object,
				  const List &rxControl,
				  const Nullable<CharacterVector> &specParams,
				  const Nullable<List> &extraArgs,
				  const RObject &params,
				  const RObject &events,
				  const RObject &inits,
				  rxSolve_t* rxSolveDat){
  bool update_params = false,
    update_events = false,
    update_inits = false;
  if (specParams.isNull()){
    warning(_("No additional parameters were specified, returning original object"));
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
  }
  // Now update
  Environment e;
  List lobj;
  if (rxSolveDat->isRxSolve){
    lobj = as<List>(object);
    CharacterVector classattr = object.attr("class");
    e = as<Environment>(classattr.attr(".RxODE.env"));
  } else  { // if (rxIs(object, "environment")) 
    e = as<Environment>(object);
    lobj = as<List>(e["obj"]);
  }
  getRxModels();
  if (e.exists(".params.dat")){
    e.remove(".params.dat");
  }
  if (e.exists(".et")){
    e.remove(".et");
  }
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
  RObject new_params;
  RObject new_events;
  if (rxIs(params, "rx.event")){
    new_events = params;
    new_params = events;
  } else {
    new_params = update_params ? params : e[".args.params"];
    new_events = update_events ? events : e[".args.events"];
  }
    
  RObject new_inits = update_inits ? inits : e[".args.inits"];
  List newRxControl = clone(rxControl);
  RObject new_object = as<RObject>(e[".args.object"]);
  newRxControl["updateObject"] = false;
  List dat = as<List>(rxSolve_(new_object, newRxControl, R_NilValue, extraArgs,
			       new_params, new_events, new_inits, 0));
  if (rxSolveDat->updateObject && as<bool>(e[".real.update"])){
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
      warning(_("can not update object"));
      return dat;
    }
  }
  e[".real.update"] = true;
  return dat;
}

// This updates the event table with sequences and other similar items
// This is useful when there is an event table without any
// observations in it.  This is the way that deSolve/mrgsolve handles
// events.
static inline void rxSolve_ev1Update(const RObject &obj,
				     const List &rxControl,
				     const Nullable<CharacterVector> &specParams,
				     const Nullable<List> &extraArgs,
				     const RObject &params,
				     RObject &ev1,
				     const RObject &inits,
				     rxSolve_t* rxSolveDat){
  rx_solve* rx = getRxSolve_();
  if (rxIs(ev1, "rxEt")){
    CharacterVector cls = ev1.attr("class");
    List etE = cls.attr(".RxODE.lst");
    int nobs = etE["nobs"];
    if (nobs == 0){
      // KEEP/DROP?
      List ev1a = etTrans(as<List>(ev1), obj, rxSolveDat->hasCmt,
			  false, false, true, R_NilValue,
			  rxControl["keepF"]);
      rxSolveDat->labelID=true;
      CharacterVector tmpC = ev1a.attr("class");
      List tmpL = tmpC.attr(".RxODE.lst");
      rxSolveDat->idLevels = as<CharacterVector>(tmpL["idLvl"]);
      rx->nKeepF = keepFcov.size();
      int lenOut = 200;
      double by = NA_REAL;
      double to;
      double from = 0.0;
      NumericVector tmp;
      IntegerVector tmpI;
      if (rxIs(rxControl["from"], "integer") ||
	  rxIs(rxControl["from"], "numeric")){
	tmp = as<NumericVector>(rxControl["from"]);
	if (tmp.size() != 1){
	  rxSolveFree();
	  stop(_("'from' must be of length 1"));
	}
	from = tmp[0];
      }
      if (rxIs(rxControl["to"], "integer") || rxIs(rxControl["to"], "numeric")){
	tmp = as<NumericVector>(rxControl["to"]);
	if (tmp.size() != 1){
	  rxSolveFree();
	  stop(_("'to' must be of length 1"));
	}
	to = tmp[0];
      } else {
	to = (max(as<NumericVector>(ev1a["TIME"]))+24);
      }
      if (rxIs(rxControl["by"], "integer") || rxIs(rxControl["by"], "numeric")){
	tmp = as<NumericVector>(rxControl["by"]);
	if (tmp.size() != 1){
	  rxSolveFree();
	  stop(_("'by' must be of length 1"));
	}
	by = tmp[0];
      }
      if (rxIs(rxControl["length.out"], "integer") || rxIs(rxControl["length.out"], "numeric")){
	tmpI = as<IntegerVector>(rxControl["length.out"]);
	if (tmpI.size() != 1){
	  rxSolveFree();
	  stop(_("'length.out' must be of length 1"));
	}
	lenOut = tmpI[0];
	if (!ISNA(by)){
	  // Matches seq(0,1,by=0.1,length.out=3)
	  rxSolveFree();
	  stop(_("cannot use both 'by' and 'length.out' for RxODE simulations"));
	}
	by = (to-from)/(lenOut-1);
      } else if (ISNA(by)) {
	lenOut=200;
	by = (to-from)/(lenOut-1);
      } else {
	lenOut= (int)((to-from)/by+1.0);
      }
      NumericVector newObs(lenOut);
      // ((to - from)/(length.out - 1))
      List et = as<List>(ev1);
      for (int i = lenOut; i--;){
	newObs[i]=by*i+from;
      }
      rx->nobs2 = lenOut;
      ev1 = et_(List::create(newObs), as<List>(ev1));
    }
  }
  if (rxIs(ev1, "data.frame") && !rxIs(ev1, "rxEtTrans")){
    ev1 = as<List>(etTrans(as<List>(ev1), obj, rxSolveDat->hasCmt,
			   false, false, true, R_NilValue,
			   rxControl["keepF"]));
    rxSolveDat->labelID=true;
    CharacterVector tmpC = ev1.attr("class");
    List tmpL = tmpC.attr(".RxODE.lst");
    rxSolveDat->idLevels = as<CharacterVector>(tmpL["idLvl"]);
    rx->nKeepF = keepFcov.size();
    rxcEvid = 2;
    rxcTime = 1;
    rxcAmt  = 3;
    rxcId   = 0;
    rxcDv   = 5;
    rxcIi   = 4;
    int censAdd = as<int>(tmpL["censAdd"]);
    int limitAdd = as<int>(tmpL["limitAdd"]);
    if (censAdd == 1 && limitAdd == 1) {
      rxcCens = 6;
      rxcLimit = 7;
    } else if (censAdd == 1){
      rxcCens = 6;
      rxcLimit = -1;
    } else if (limitAdd == 1){
      rxcCens  = -1;
      rxcLimit = 6;
    } else {
      rxcCens  = -1;
      rxcLimit = -1;
    }
  }
  if (rxIs(ev1, "rxEtTrans")){
    CharacterVector cls = ev1.attr("class");
    List tmpL = cls.attr(".RxODE.lst");
    rx->nobs2 = as<int>(tmpL["nobs"]);
  }
  _rxModels[".lastEv1"] = ev1;
}


// This functin simulates individual parameter values and residual
// parameter values and then converts them to a data.frame.  This
// allows rxSolve_ to solve as if the user specified these parameters
// directly
static inline void rxSolve_simulate(const RObject &obj,
				    const List &rxControl,
				    const Nullable<CharacterVector> &specParams,
				    const Nullable<List> &extraArgs,
				    const RObject &params,
				    const RObject &ev1,
				    const RObject &inits,
				    rxSolve_t* rxSolveDat){
  rx_solve* rx = getRxSolve_();
  rx_solving_options* op = rx->op;
  
  RObject omega = rxControl["omega"];
  Nullable<NumericVector> omegaDf = as<Nullable<NumericVector>>(rxControl["omegaDf"]);
  bool omegaIsChol = as<bool>(rxControl["omegaIsChol"]);

  Nullable<NumericMatrix> thetaMat = as<Nullable<NumericMatrix>>(rxControl["thetaMat"]);
  Nullable<NumericVector> thetaDf = as<Nullable<NumericVector>>(rxControl["thetaDf"]);
  bool thetaIsChol = as<bool>(rxControl["thetaIsChol"]);
  
  RObject sigma= rxControl["sigma"];
  Nullable<NumericVector> sigmaDf= (rxControl["sigmaDf"]);
  bool sigmaIsChol= as<bool>(rxControl["sigmaIsChol"]);
  op->isChol = (int)(sigmaIsChol);

  unsigned int nSub = as<unsigned int>(rxControl["nSub"]);
  unsigned int nStud = as<unsigned int>(rxControl["nStud"]);
  double dfSub=as<double>(rxControl["dfSub"]);
  double dfObs=as<double>(rxControl["dfObs"]);
  
  int nCoresRV = as<int>(rxControl["nCoresRV"]);
  
  bool simSubjects = false;
  
  op->ncoresRV = nCoresRV;

  if (!thetaMat.isNull() || !rxIs(omega, "NULL") || !rxIs(sigma, "NULL")){
    // Simulated Variable3
    if (!rxIs(rxSolveDat->par1, "numeric")){
      rxSolveFree();
      stop(_("when specifying 'thetaMat', 'omega', or 'sigma' the parameters cannot be a 'data.frame'/'matrix'."));
    }
    unsigned int nSub0 = 0;
    int curObs = 0;
    rx->nevid9 = 0;
    rx->nall = 0;
    rx->nobs = 0;
    rx->nobs2 = 0;
    if (rxIs(ev1,"event.data.frame")||
	rxIs(ev1,"event.matrix")){
      if (rxcId > -1){
	DataFrame dataf = as<DataFrame>(ev1);
	IntegerVector id = as<IntegerVector>(dataf[rxcId]);
	IntegerVector evid  = as<IntegerVector>(dataf[rxcEvid]);
	int lastid= id[id.size()-1]+42;
	rx->nall = evid.size();
	int evid9=0;
	for (unsigned int j = rx->nall; j--;){
	  if (lastid != id[j]){
	    lastid=id[j];
	    nSub0++;
	  }
	  if (isObs(evid[j])) rx->nobs++;
	  if (evid[j] == 0) rx->nobs2++;
	  if (evid[j] == 9) evid9++;
	}
	rx->nevid9 = evid9;
      } else {
	nSub0 =1;
	DataFrame dataf = as<DataFrame>(ev1);
	IntegerVector evid  = as<IntegerVector>(dataf[rxcEvid]);
	rx->nall = evid.size();
	int evid9=0;
	for (unsigned int j =rx->nall; j--;){
	  if (isObs(evid[j])) rx->nobs++;
	  if (evid[j] == 0) rx->nobs2++;
	  if (evid[j] == 9) evid9++;
	}
	rx->nevid9= evid9;
      }
    }
    if (nStud > 1 && nSub0 == nSub*nStud){
    } else if (nSub > 1 && nSub0 > 1 && nSub != nSub0){
      rxSolveFree();
      stop(_("provided multi-subject data (n=%d) trying to simulate a different number of subjects (n=%d)"), nSub0, nSub);
    } else if (nSub > 1 && nSub0 == 1) {
      nSub0 = nSub;
      simSubjects = true;
    }
    if (rxSolveDat->addDosing.isNull()){
      // only evid=0
      curObs= rx->nobs2;
    } else {
      LogicalVector addDosing1 = as<LogicalVector>(rxSolveDat->addDosing);
      if (LogicalVector::is_na(addDosing1[0])){
	curObs = rx->nall - rx->nevid9;
      } if (addDosing1[0]){
	curObs = rx->nall - rx->nevid9;
      } else {
	curObs = rx->nobs - rx->nevid9;
      }
    }
    if (rxIs(as<RObject>(thetaMat), "matrix")){
      if (!thetaIsChol){
	arma::mat tmpM = as<arma::mat>(thetaMat);
	if (!tmpM.is_sympd()){
	  rxSolveFree();
	  stop(_("'thetaMat' must be symmetric"));
	}
      }
    }
    List lst = rxSimThetaOmega(as<Nullable<NumericVector>>(rxSolveDat->par1),
			       omega,
			       omegaDf,
			       as<NumericVector>(rxControl["omegaLower"]),
			       as<NumericVector>(rxControl["omegaUpper"]),
			       omegaIsChol,
			       as<std::string>(rxControl["omegaSeparation"]),
			       as<int>(rxControl["omegaXform"]),
			       nSub0, thetaMat,
			       as<NumericVector>(rxControl["thetaLower"]),
			       as<NumericVector>(rxControl["thetaUpper"]),
			       thetaDf, thetaIsChol, nStud,
			       sigma,
			       as<NumericVector>(rxControl["sigmaLower"]),
			       as<NumericVector>(rxControl["sigmaUpper"]),
			       sigmaDf, sigmaIsChol,
			       as<std::string>(rxControl["sigmaSeparation"]),
			       as<int>(rxControl["sigmaXform"]),
			       nCoresRV, curObs,
			       dfSub, dfObs, simSubjects);
    rxSolveDat->warnIdSort = false;
    rxSolveDat->par1 =  as<RObject>(lst);
    rxSolveDat->usePar1=true;
    // The parameters are in the same format as they would be if they were
    // specified as part of the original dataset.
  }
  // .sigma could be reassigned in an update, so check outside simulation function.
  if (_rxModels.exists(".sigma")){
    rxSolveDat->sigmaN= as<CharacterVector>((as<List>((as<NumericMatrix>(_rxModels[".sigma"])).attr("dimnames")))[1]);
  }
}

// This will setup the parNumeric, parDf, or parMat for solving. It
// will also set the parameter type.
static inline void rxSolve_parSetup(const RObject &obj,
				    const List &rxControl,
				    const Nullable<CharacterVector> &specParams,
				    const Nullable<List> &extraArgs,
				    const CharacterVector& pars,
				    const RObject &ev1,
				    const RObject &inits,
				    rxSolve_t* rxSolveDat){
  rx_solve* rx = getRxSolve_();
  RObject theta = rxControl["theta"];
  RObject eta = rxControl["eta"];
  if (!theta.isNULL() || !eta.isNULL()){
    rxSolveDat->usePar1=true;
    rxSolveDat->par1 = rxSetupParamsThetaEta(rxSolveDat->par1, theta, eta);
  }
  if (rxIs(rxSolveDat->par1, "numeric") || rxIs(rxSolveDat->par1, "integer")){
    rxSolveDat->parNumeric = as<NumericVector>(rxSolveDat->par1);
    if (rxSolveDat->parNumeric.hasAttribute("names")){
      rxSolveDat->nmP = rxSolveDat->parNumeric.names();
    } else if (rxSolveDat->parNumeric.size() == pars.size()) {
      rxSolveDat->nmP = pars;
    } else {
      rxSolveFree();
      stop(_("if parameters are not named, they must match the order and size of the parameters in the model"));
    }
    RObject iCov = rxControl["iCov"];
    if (!rxIs(iCov, "NULL")){
      // Create a data frame
      Function sortId = getRxFn(".sortId");
      iCov = clone(sortId(iCov, rxSolveDat->idLevels, "iCov", rxSolveDat->warnIdSort));
      CharacterVector keepC, keepCf;
      if (rxIs(rxControl["keepI"], "character")){
	keepC = as<CharacterVector>(rxControl["keepI"]);
      }
      IntegerVector keepIv(keepC.size());
      std::fill_n(keepIv.begin(), keepC.size(), -1);
      List lstT=as<List>(iCov);
      rxSolveDat->parDf = as<DataFrame>(iCov);
      int nr = rxSolveDat->parDf.nrows();
      List lstF(rxSolveDat->parNumeric.size()+lstT.size());
      CharacterVector nmF(lstF.size());
      CharacterVector nmL = rxSolveDat->nmP;
      CharacterVector nmT = lstT.attr("names");
      for (int ii = rxSolveDat->parNumeric.size(); ii--;){
	NumericVector tmp(nr);
	std::fill_n(tmp.begin(), nr, rxSolveDat->parNumeric[ii]);
	lstF[ii] = tmp;
	nmF[ii] = nmL[ii];
      }
      int nKeepi=0;
      for (int ii=lstT.size(); ii--;){
	lstF[ii+rxSolveDat->parNumeric.size()] = lstT[ii];
	nmF[ii+rxSolveDat->parNumeric.size()] = nmT[ii];
	for (int jj = keepC.size(); jj--;){
	  if (nmT[ii] == keepC[jj]){
	    keepIv[jj] = ii;
	    nKeepi++;
	    break;
	  }
	}
      }
      keepIcov=List(nKeepi);
      keepCf = CharacterVector(nKeepi);
      int iii=0;
      for (int ii=keepC.size(); ii--;){
	if (keepIv[ii] != -1){
	  keepIcov[iii] = lstT[keepIv[ii]];
	  keepCf[iii++] = nmT[keepIv[ii]];
	}
      }
      keepIcov.attr("names") = keepCf;
      keepIcov.attr("class") = "data.frame";
      keepIcov.attr("row.names") = lstT.attr("row.names");
      rx->nKeep0 = nKeepi;
      lstF.attr("names") = nmF;
      lstF.attr("class") = "data.frame";
      lstF.attr("row.names") = lstT.attr("row.names");
      rxSolveDat->par1 = as<RObject>(lstF);
      rxSolveDat->parDf = as<DataFrame>(rxSolveDat->par1);
      rxSolveDat->parType = 2;
      rxSolveDat->nmP = rxSolveDat->parDf.names();
      rxSolveDat->nPopPar = rxSolveDat->parDf.nrows();
      rxSolveDat->usePar1=true;
    }
  } else if (rxIs(rxSolveDat->par1, "data.frame")){
    Function sortId = getRxFn(".sortId");
    if (rxSolveDat->idLevels.size() > 0){
      rxSolveDat->par1 = clone(sortId(rxSolveDat->par1, rxSolveDat->idLevels, "parameters", rxSolveDat->warnIdSort));
      rxSolveDat->usePar1=true;
    }
    RObject iCov = rxControl["iCov"];
    if (!rxIs(iCov, "NULL")){
      Function sortId = getRxFn(".sortId");
      iCov = clone(sortId(iCov, rxSolveDat->idLevels, "iCov", rxSolveDat->warnIdSort));
      List lstT=as<List>(iCov);
      List lst = as<List>(rxSolveDat->par1);
      List lstF(lst.size()+lstT.size());
      CharacterVector nmF(lstF.size());
      CharacterVector nmL = lst.attr("names");
      CharacterVector nmT = lstT.attr("names");
      for (int ii = lst.size(); ii--;){
	lstF[ii] = lst[ii];
	nmF[ii] = nmL[ii];
      }
      int nKeepi=0;
      CharacterVector keepC, keepCf;
      if (rxIs(rxControl["keepI"], "character")){
	keepC = as<CharacterVector>(rxControl["keepI"]);
      }
      IntegerVector keepIv(keepC.size());
      for (int ii=lstT.size(); ii--;){
	lstF[ii+lst.size()] = lstT[ii];
	nmF[ii+lst.size()] = nmT[ii];
	for (int jj = keepC.size(); jj--;){
	  if (nmT[ii] == keepC[jj]){
	    keepIv[jj] = ii;
	    nKeepi++;
	    break;
	  }
	}
      }
      keepIcov=List(nKeepi);
      keepCf = CharacterVector(nKeepi);
      int iii=0;
      for (int ii=keepC.size(); ii--;){
	if (keepIv[ii] != -1){
	  keepIcov[iii] = lstT[keepIv[ii]];
	  keepCf[iii++] = nmT[keepIv[ii]];
	}
      }
      keepIcov.attr("names") = keepCf;
      keepIcov.attr("class") = "data.frame";
      keepIcov.attr("row.names") = lstT.attr("row.names");
      rx->nKeep0 = nKeepi;
      lstF.attr("names") = nmF;
      lstF.attr("class") = "data.frame";
      lstF.attr("row.names") = lst.attr("row.names");
      rxSolveDat->par1 = as<RObject>(lstF);
    }
    rxSolveDat->parDf = as<DataFrame>(rxSolveDat->par1);
    rxSolveDat->parType = 2;
    rxSolveDat->nmP = rxSolveDat->parDf.names();
    rxSolveDat->nPopPar = rxSolveDat->parDf.nrows();
  } else if (rxIs(rxSolveDat->par1, "matrix")){
    RObject iCov = rxControl["iCov"];
    if (!rxIs(iCov, "NULL")){
      rxSolveFree();
      stop(_("matrix parameters with 'iCov' 'data.frame' is not supported"));
    }
    rxSolveDat->parMat = as<NumericMatrix>(rxSolveDat->par1);
    rxSolveDat->nPopPar = rxSolveDat->parMat.nrow();
    rxSolveDat->parType = 3;
    if (rxSolveDat->parMat.hasAttribute("dimnames")){
      Nullable<CharacterVector> colnames0 = as<Nullable<CharacterVector>>((as<List>(rxSolveDat->parMat.attr("dimnames")))[1]);
      if (colnames0.isNull()){
	if (rxSolveDat->parMat.ncol() == pars.size()){
	  rxSolveDat->nmP = pars;
	} else {
	  rxSolveFree();
	  stop(_("if parameters are not named, they must match the order and size of the parameters in the model"));
	}
      } else {
	rxSolveDat->nmP = CharacterVector(colnames0);
      }
    } else if (rxSolveDat->parMat.ncol() == pars.size()) {
      rxSolveDat->nmP = pars;
    } else {
      rxSolveFree();
      stop(_("if parameters are not named, they must match the order and size of the parameters in the model"));
    }
  }
}

// This loops through the data to put each individual into the
// approiate data structure.
// At the same time calculate hmax per individual as well
// in the overall population.
static inline void rxSolve_datSetupHmax(const RObject &obj, const List &rxControl,
					const Nullable<CharacterVector> &specParams,
					const Nullable<List> &extraArgs,
					CharacterVector& pars,
					const RObject &ev1,
					const RObject &inits,
					rxSolve_t* rxSolveDat){
  Nullable<NumericVector> hmax = rxControl["hmax"];
  bool doMean=true;
  rx_solve* rx = getRxSolve_();
  rx_solving_options* op = rx->op;
  rx_solving_options_ind* ind;
  double tmp, hmax1 = 0.0, hmax1m = 0.0, hmax1mo, hmax1s=0.0,
      hmax1n=0.0, hmax2 = 0.0, hmax2m = 0.0, hmax2mo, hmax2s=0.0,
      hmax2n=0.0, tlast;
  int i, j;
  unsigned int nsub = 0;
  unsigned int nobs = 0, ndoses = 0;

  int ncov =0, curcovi = 0;

  if (rxIs(ev1,"event.data.frame")||
      rxIs(ev1,"event.matrix")){
    // data.frame or matrix
    double hmax0 = 0.0;
    if (!hmax.isNull()){
      NumericVector hmax00(hmax);
      hmax0 = hmax00[0];
      if (NumericVector::is_na(hmax00[0])){
	hmax0 = 0.0;
	doMean = true;
      } else {
	hmax0 = hmax00[0];
      }
    }
    DataFrame dataf = as<DataFrame>(ev1);
    CharacterVector dfNames = dataf.names();
    int dfN = dfNames.size();
    IntegerVector evid  = as<IntegerVector>(dataf[rxcEvid]);
    if (_globals.gevid != NULL) free(_globals.gevid);
    _globals.gevid = (int*)calloc(3*evid.size()+dfN, sizeof(int));
    if (_globals.gevid == NULL){
      rxSolveFree();
      stop(_("can not allocate enough memory to load 'evid'"));
    }
    std::copy(evid.begin(),evid.end(), &_globals.gevid[0]);
    _globals.gidose = _globals.gevid + evid.size();
    _globals.gcens = _globals.gidose + evid.size();
    _globals.gpar_cov = _globals.gidose + evid.size();//[dfN];

    int ntot = 1;

    IntegerVector id(evid.size(), 1);
    if (rxcId > -1){
      id    = as<IntegerVector>(dataf[rxcId]);
      int lastid = NA_INTEGER;
      int nid = 0;
      for (int ii = 0; ii < id.size(); ii++){
	if (id[ii] != lastid){
	  lastid = id[ii];
	  nid++;
	}
      }
      // if (nid == 0){
      // } else 
      if (nid == rxSolveDat->nPopPar || rxSolveDat->nPopPar == 1){
	ntot = nid;
      } else {
	ntot = nid*rxSolveDat->nPopPar;
      }
    } else {
      if (rxSolveDat->nPopPar != 0) ntot = rxSolveDat->nPopPar;
    }
    
    NumericVector time0 = dataf[rxcTime];
    if (rxIs(time0, "units")){
      rxSolveDat->addTimeUnits=true;
      rxSolveDat->timeUnitsU=time0.attr("units");
    }

    tlast = time0[0];
    // - all_times
    if (ntot <= 0) {
      rxSolveFree();
      stop(_("nothing to solve"));
    }
    rxOptionsIniEnsure(ntot);
    if (_globals.gall_times != NULL) free(_globals.gall_times);
    _globals.gall_times = (double*)calloc(5*time0.size(), sizeof(double));
    std::copy(time0.begin(), time0.end(), &_globals.gall_times[0]);
    _globals.gdv = _globals.gall_times + time0.size(); // Perhaps allocate zero size if missing?
    _globals.gamt = _globals.gdv + time0.size(); 
    _globals.gii = _globals.gamt + time0.size();
    _globals.glimit=_globals.gii + time0.size();
    NumericVector dv;
    if (rxcDv > -1){
      dv = as<NumericVector>(dataf[rxcDv]);
      std::copy(dv.begin(), dv.end(), &_globals.gdv[0]);
    }
    IntegerVector cens;
    rx->cens=0;
    if (rxcCens > -1){
      cens = as<IntegerVector>(dataf[rxcCens]);
      std::copy(cens.begin(), cens.end(), &_globals.gcens[0]);
      rx->cens=1;
    }
    NumericVector limit;
    rx->limit=0;
    if (rxcLimit > -1){
      rx->limit=1;
      limit = as<NumericVector>(dataf[rxcLimit]);
      std::copy(limit.begin(), limit.end(), &_globals.glimit[0]);
    }
    NumericVector amt   = dataf[rxcAmt];
    NumericVector datIi(amt.size());
    if (rxcIi > -1){
      datIi = as<NumericVector>(dataf[rxcIi]);
    } else {
      std::fill_n(datIi.begin(), datIi.size(), 0.0);
    }
    ncov = 0;
    std::vector<int> covPos(dfN);
    curcovi=0;
    for (i = dfN; i--;){
      for (j = rxSolveDat->npars; j--;){
	if (pars[j] == dfNames[i]){
	  _globals.gpar_cov[ncov] = j+1;
	  covPos[ncov] = i;
	  ncov++;
	}
      }
    }
    // Make sure the covariates are a #ncov * all times size
    if (_globals.gcov != NULL) free(_globals.gcov);
    _globals.gcov = (double*)calloc(ncov * amt.size(), sizeof(double));
    if (_globals.gcov == NULL){
      rxSolveFree();
      stop(_("can not allocate memory for the covariates"));
    }
    op->par_cov=&(_globals.gpar_cov[0]);
    op->ncov=ncov;
    op->do_par_cov = (ncov > 0);
    // Make sure the covariates are a #ncov * all times size
    int ids = id.size();
    // Get the number of subjects
    // Get the number of observations
    // Get the number of doses
    int nall = 0, nobst=0, lasti =0, ii=0, nobs2t=0, nevid9=0;
    nsub = 0;
    ind = &(rx->subjects[0]);
    ind->idx=0;
    j=0;
    int lastId = id[0]-42;
    for (i = 0; i < ids; i++){
      if (lastId != id[i]){
	if (nall != 0){
	  // Finalize last solve.
	  ind->n_all_times    = ndoses+nobs;
	  ind->cov_ptr = &(_globals.gcov[curcovi]);
	  for (ii = 0; ii < ncov; ii++){
	    NumericVector cur = as<NumericVector>(dataf[covPos[ii]]);
	    std::copy(cur.begin()+lasti, cur.begin()+lasti+ind->n_all_times,
		      _globals.gcov+curcovi);
	    curcovi += ind->n_all_times;
	  }
	  nsub++;
	  if (doMean){
	    hmax1s = hmax1s/(hmax1n-1);
	    hmax1  = hmax1m;
	    ind->HMAX = hmax1;
	  } else if (hmax0 == 0.0){
	    ind->HMAX = hmax1;
	  } else {
	    ind->HMAX = hmax0;
	  }
	  ind = &(rx->subjects[nsub]);
	  ind->idx=0;
	}
	// Setup the pointers.
	ind->id             = nsub+1;
	ind->bT = NA_REAL;
	ind->allCovWarn = 0;
	ind->wrongSSDur=0;
	ind->err = 0;
	ind->linCmtT = NA_REAL;
	ind->linCmtAdvanSetup=0;
	ind->cacheME=0;
	ind->timeReset=1;
	ind->lambda         =1.0;
	ind->yj             = 0;
	ind->doSS = 0;
	ind->all_times      = &_globals.gall_times[i];
	ind->dv = &_globals.gdv[i];
	ind->limit = &_globals.glimit[i];
	ind->cens = &_globals.gcens[i];
	ind->evid           = &_globals.gevid[i];
	ind->idose          = &_globals.gidose[i];
	ind->dose           = &_globals.gamt[i];
	ind->ii             = &_globals.gii[i];
	lasti = i;

	hmax1m=0.0;
	hmax1s=0.0;
	hmax1n=0.0;
	hmax1 = 0.0;
	lastId=id[i];
	j=i;
	ind->ndoses=0;
	ind->nevid2=0;
	ndoses=0;
	nobs=0;
	tlast = NA_REAL;
      }
      // Create index
      if (isDose(_globals.gevid[i])){
	_globals.gidose[j] = i-lasti;
	_globals.gii[j] = datIi[i];
	_globals.gamt[j] = amt[i];
	ind->ndoses++;
	ndoses++; nall++; j++;
      } else {
	nobs++; nobst++; nall++;
	if (_globals.gevid[i] == 2) ind->nevid2++;
	if (_globals.gevid[i] == 0) nobs2t++;
	if (_globals.gevid[i] == 9) nevid9++;
	if (!ISNA(tlast)){
	  tmp = time0[i]-tlast;
	  if (tmp < 0){
	    rxSolveFree();
	    stop(_("data must be ordered by 'ID' and 'TIME' variables"));
	  }
	  hmax1n++;
	  hmax1mo = hmax1m;
	  hmax1m += (tmp-hmax1m)/hmax1n;
	  hmax1s += (tmp-hmax1m)*(tmp-hmax1mo);
	  hmax2n++;
	  hmax2mo = hmax2m;
	  hmax2m += (tmp-hmax2m)/hmax2n;
	  hmax2s += (tmp-hmax2m)*(tmp-hmax2mo);
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
    if (doMean){
      hmax2  = hmax2m;
    }
    rx->nobs = nobst;
    rx->nobs2 = nobs2t;
    rx->nall = nall;
    rx->nevid9 = nevid9;
    // Finalize the prior individual
    ind->n_all_times    = ndoses+nobs;
    ind->cov_ptr = &(_globals.gcov[curcovi]);
    for (ii = 0; ii < ncov; ii++){
      NumericVector cur = as<NumericVector>(dataf[covPos[ii]]);
      std::copy(cur.begin()+lasti, cur.begin()+lasti+ind->n_all_times,
		_globals.gcov+curcovi);
      curcovi += ind->n_all_times;
    }
    if (doMean || hmax0 == 0.0){
      op->hmax2 = hmax2;
    } else {
      op->hmax2 = hmax0;
    }
    nsub++;
    rx->nsub= nsub;
    if (doMean){
      hmax1s = hmax1s/(hmax1n-1);
      hmax1  = hmax1m;
      ind->HMAX = hmax1 + as<double>(rxControl["hmaxSd"])*sqrt(hmax1s);
    } else if (hmax0 == 0.0){
      ind->HMAX = hmax1;
    } else {
      ind->HMAX = hmax0;
    }
  }
}

// Setup the parameter order and covariate from the dataset and input parameters
static inline void rxSolve_parOrder(const RObject &obj, const List &rxControl,
				    const Nullable<CharacterVector> &specParams,
				    const Nullable<List> &extraArgs,
				    CharacterVector& pars,
				    const RObject &ev1,
				    const RObject &inits,
				    rxSolve_t* rxSolveDat){
  rx_solve* rx = getRxSolve_();
  rx_solving_options* op = rx->op;
  if (_globals.gParPos != NULL) free(_globals.gParPos);
  _globals.gParPos = (int*)calloc(rxSolveDat->npars*2 + rxSolveDat->sigmaN.size(), sizeof(int));// [npars]
  if (_globals.gParPos == NULL){
    rxSolveFree();
    stop(_("cannot allocate enough memory to sort input parameters"));
  }
  _globals.gParPos2 =  _globals.gParPos + rxSolveDat->npars; // [npars]
  _globals.gsvar =  _globals.gParPos2 + rxSolveDat->npars;//[sigmaN.size()]
  std::string errStr = "";
  bool allPars = true;
  bool curPar = false;
  CharacterVector mvIniN = rxSolveDat->mvIni.names();
  CharacterVector mvCov1N;
  if (rxIs(ev1, "rxEtTran")){
    CharacterVector tmpCls = ev1.attr("class");
    List e = tmpCls.attr(".RxODE.lst");
    List tmpCov1 = e["cov1"];
    mvCov1N = tmpCov1.attr("names");
  }
  int i, j;
  for (i = rxSolveDat->npars; i--;){
    curPar = false;
    // Check to see if this is a covariate.
    for (j = op->ncov; j--;){
      if (_globals.gpar_cov[j] == (int)(i + 1)){
	_globals.gParPos[i] = 0; // These are set at run-time and "dont" matter.
	curPar = true;
	rxSolveDat->eGparPos[i]=_globals.gParPos[i];
	break;
      }
    }
    // Check for the sigma-style simulated parameters.
    if (!curPar){
      for (j = rxSolveDat->sigmaN.size(); j--;){
	if (rxSolveDat->sigmaN[j] == pars[i]){
	  _globals.gsvar[j] = i;
	  rxSolveDat->nsvar++;
	  _globals.gParPos[i] = 0; // These are set at run-time and "dont" matter.
	  curPar = true;
	  rxSolveDat->eGparPos[i]=_globals.gParPos[i];
	  break;
	}
      }
    }
    // Next, check to see if this is a user-specified parameter
    if (!curPar){
      for (j = rxSolveDat->nmP.size(); j--;){
	if (rxSolveDat->nmP[j] == pars[i]){
	  curPar = true;
	  _globals.gParPos[i] = j + 1;
	  rxSolveDat->eGparPos[i]=_globals.gParPos[i];
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
	  rxSolveDat->eGparPos[i]=_globals.gParPos[i];
	  break;
	}
      }
    }
    if (!curPar){
      for (j = 1; j < mvCov1N.size(); j++){
	if (mvCov1N[j] == pars[i]){
	  // These are setup once and don't need to be updated
	  _globals.gParPos[i] = 0; // These are set at run-time and "dont" matter.
	  curPar = true;
	  rxSolveDat->eGparPos[i]=_globals.gParPos[i];
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
    CharacterVector modSyntax = rxSolveDat->mv["model"];
    Rcout << "Model:\n\n" + modSyntax[0] + "\n";
    rxSolveFree();
    stop(errStr);
  }
}

static inline void rxSolve_assignGpars(rxSolve_t* rxSolveDat);

// Normalizes parameter input
// It will also load the parameters into the solving data structure.
static inline void rxSolve_normalizeParms(const RObject &obj, const List &rxControl,
					  const Nullable<CharacterVector> &specParams,
					  const Nullable<List> &extraArgs,
					  const CharacterVector& pars,
					  const RObject &ev1,
					  const RObject &inits,
					  rxSolve_t* rxSolveDat){
  int i;
  rx_solve* rx = getRxSolve_();
  rx_solving_options* op = rx->op;
  rx_solving_options_ind* ind;
  int curEvent = 0, curIdx = 0, curSolve=0, curLin=0;
  switch(rxSolveDat->parType){
  case 1: // NumericVector
    {
      if (rxSolveDat->nPopPar != 1) {
	rxSolveFree();
	stop(_("nPopPar != 1 but parameters are specified as a NumericVector"));
      }
      // Convert to DataFrame to simplify code.
      List parDfL(rxSolveDat->parNumeric.size());
      for (i = rxSolveDat->parNumeric.size(); i--;){
	parDfL[i] = NumericVector(rx->nsub,rxSolveDat->parNumeric[i]);
      }
      parDfL.attr("names") = rxSolveDat->nmP;
      parDfL.attr("row.names") = IntegerVector::create(NA_INTEGER,-rx->nsub);
      parDfL.attr("class") = CharacterVector::create("data.frame");
      rxSolveDat->parDf = as<DataFrame>(parDfL);
      rxSolveDat->parType = 2;
      rxSolveDat->nPopPar = rxSolveDat->parDf.nrows();
      rx->nsim=1;
    }
  case 2: // DataFrame
    // Convert to NumericMatrix
    {
      unsigned int parDfi = rxSolveDat->parDf.size();
      rxSolveDat->parMat = NumericMatrix(rxSolveDat->nPopPar,parDfi);
      while (parDfi--){
	rxSolveDat->parMat(_,parDfi)=NumericVector(rxSolveDat->parDf[parDfi]);
      }
      rxSolveDat->parMat.attr("dimnames") = List::create(R_NilValue, rxSolveDat->parDf.names());
    }
  case 3: // NumericMatrix
    {
      gparsCovSetup(rxSolveDat->npars, rxSolveDat->nPopPar, ev1, rx);
      rxSolve_assignGpars(rxSolveDat);
      curSolve=0;
      curLin=0;
      curEvent=0;
      curIdx=0;
      int curOn=0;
      rx_solving_options_ind indS;
      for (unsigned int simNum = rx->nsim; simNum--;){
	for (unsigned int id = rx->nsub; id--;){
	  unsigned int cid = id+simNum*rx->nsub;
	  ind = &(rx->subjects[cid]);
	  ind->idx=0;
	  ind->par_ptr = &_globals.gpars[cid*rxSolveDat->npars];
	  ind->mtime   = &_globals.gmtime[rx->nMtime*cid];
	  if (rx->nMtime > 0) ind->mtime[0]=-1;
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
	    ind->nevid2 = indS.nevid2;
	    ind->dose = &(indS.dose[0]);
	    ind->ii   = &(indS.ii[0]);
	    ind->evid =&(indS.evid[0]);
	    ind->all_times = &(indS.all_times[0]);
	    ind->id=id+1;
	    ind->lambda=1.0;
	    ind->yj = 0;
	    ind->doSS = 0;
	  }
	  int eLen = op->neq*ind->n_all_times;
	  int aLen = op->nlin;
	  if (aLen != 0) aLen = (aLen+2)*ind->n_all_times;
	  ind->linCmtAdvan = &_globals.gadvan[curLin];
	  curLin += aLen;
	  ind->solve = &_globals.gsolve[curSolve];
	  curSolve += eLen;
	  ind->solveLast = &_globals.gsolve[curSolve];
	  curSolve += op->neq;
	  ind->solveLast2 = &_globals.gsolve[curSolve];
	  curSolve += op->neq;
	  ind->solveSave = &_globals.gsolve[curSolve];
	  curSolve += op->neq;
	  ind->ix = &_globals.gix[curIdx];
	  std::iota(ind->ix,ind->ix+ind->n_all_times,0);
	  curEvent += eLen;
	  ind->on=&_globals.gon[curOn];
	  curOn +=op->neq;
	  curIdx += ind->n_all_times;
	  ind->_newind = -1;
	}
      }
    }
    break;
  default:
    rxSolveFree();
    stop(_("Something is wrong"));
  }
}

// This creates the final dataset from the currently solved object.
// Most of this is a direct C call, but some items are done in C++
static inline List rxSolve_df(const RObject &obj,
			      const List &rxControl,
			      const Nullable<CharacterVector> &specParams,
			      const Nullable<List> &extraArgs,
			      const RObject &params,
			      const RObject &events,
			      const RObject &inits,
			      rxSolve_t* rxSolveDat){
  rx_solve* rx = getRxSolve_();
  rx_solving_options* op = rx->op;
  if (op->abort){
    rxSolveFree();
    stop(_("aborted solve"));
  }
  int doDose = 0;
  if (rxSolveDat->addDosing.isNull()){
    // only evid=0
    doDose=-1;
  } else {
    LogicalVector addDosing1 = as<LogicalVector>(rxSolveDat->addDosing);
    if (LogicalVector::is_na(addDosing1[0])){
      doDose = 1;
    } else if (addDosing1[0]){
      doDose = 2;
      if (as<bool>(rxControl["subsetNonmem"])) doDose = 3;
    } else {
      doDose = 0;
    }
  }
  IntegerVector si = rxSolveDat->mv["state.ignore"];
  rx->stateIgnore = &si[0];
  int doTBS = (rx->matrix == 3);
  if (doTBS) rx->matrix=2;
  if (rx->matrix == 4 || rx->matrix == 5) rx->matrix=2;
  List dat = RxODE_df(doDose, doTBS);
  if (!rxIs(rxControl["drop"], "NULL")) {
    dat = rxDrop(as<CharacterVector>(rxControl["drop"]), dat, as<bool>(rxControl["warnDrop"]));
  }
  if (rxSolveDat->idFactor && rxSolveDat->labelID && rx->nsub > 1){
    IntegerVector did = as<IntegerVector>(dat["id"]);
    did.attr("class") = "factor";
    did.attr("levels") = rxSolveDat->idLevels;
  }
  if (rxSolveDat->addTimeUnits){
    NumericVector tmpN = as<NumericVector>(dat["time"]);
    tmpN.attr("class") = "units";
    tmpN.attr("units") = rxSolveDat->timeUnitsU;
  }
  dat.attr("class") = CharacterVector::create("data.frame");
  if (rx->add_cov && (rx->matrix == 2 || rx->matrix == 0) &&
      rxSolveDat->covUnits.hasAttribute("names")){
    CharacterVector nmC = rxSolveDat->covUnits.attr("names");
    NumericVector tmpN, tmpN2;
    for (int i = nmC.size(); i--;){
      tmpN = rxSolveDat->covUnits[i];
      if (rxIs(tmpN, "units")){
	tmpN2 = dat[as<std::string>(nmC[i])];
	tmpN2.attr("class") = "units";
	tmpN2.attr("units") = tmpN.attr("units");
      }
    }
  }
  return dat;
}

// Generate environment for saving solving information
static inline Environment rxSolve_genenv(const RObject &object,
					 const List &rxControl,
					 const Nullable<CharacterVector> &specParams,
					 const List &dat,
					 const RObject &params,
					 const RObject &events,
					 const RObject &inits,
					 rxSolve_t* rxSolveDat){
  rx_solve* rx = getRxSolve_();
  IntegerVector slvr_counterIv(rxSolveDat->nSize);
  IntegerVector dadt_counterIv(rxSolveDat->nSize);
  IntegerVector  jac_counterIv(rxSolveDat->nSize);
  CharacterVector timeUnits = as<CharacterVector>(rxControl["timeUnits"]);
  CharacterVector amountUnits = as<CharacterVector>(rxControl["amountUnits"]);

  std::copy(&_globals.slvr_counter[0], &_globals.slvr_counter[0] + rxSolveDat->nSize, slvr_counterIv.begin());
  std::copy(&_globals.dadt_counter[0], &_globals.dadt_counter[0] + rxSolveDat->nSize, dadt_counterIv.begin());
  std::copy(&_globals.jac_counter[0], &_globals.jac_counter[0] + rxSolveDat->nSize, jac_counterIv.begin());
  Function newEnv("new.env", R_BaseNamespace);
  Environment e = newEnv(_["size"] = 29, _["parent"] = RxODEenv());
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
  e[".check.nrow"] = rx->nr;
  e[".check.ncol"] = dat.size();
  e[".check.names"] = dat.names();
  e[".idLevels"] = as<CharacterVector>(rxSolveDat->idLevels);
  e[".par.pos"] = rxSolveDat->eGparPos;
  e[".par.pos.ini"] = rxSolveDat->fromIni;
  e[".slvr.counter"] = slvr_counterIv;
  e[".dadt.counter"] = dadt_counterIv;
  e[".jac.counter"] = jac_counterIv;
  e[".nsub"] = rx->nsub;
  e[".nsim"] = rx->nsim;
  e[".init.dat"] = rxSolveDat->initsC;
  CharacterVector units = CharacterVector::create(amountUnits[0], timeUnits[0]);
  units.names() = CharacterVector::create("dosing","time");
  e["units"] = units;
  e["nobs"] = rx->nobs - rx->nevid9;
  e[".args.object"] = object;
  e["dll"] = rxDll(object);
  e[".args.par0"] = rxSolveDat->par1ini;
  if (!rxSolveDat->swappedEvents){
    if (rxSolveDat->usePar1){
      e[".args.params"] = rxSolveDat->par1;
    } else {
      e[".args.params"] = params;
    }
    e[".args.events"] = events;
  } else {
    if (rxSolveDat->usePar1){
      e[".args.params"] = rxSolveDat->par1;
    } else {
      e[".args.params"] = events;
    }
    e[".args.events"] = params;
  }
  e[".args.inits"] = inits;
  e[".args"] = rxControl;
  e[".real.update"] = true;
  return e;
}

extern "C" void seedEng(int ncores);
void rxAssignPtr(SEXP object);

int getNRows(RObject obj){
  switch (obj.sexp_type()){
  case 13: // Integer Vectors
  case REALSXP:
    {
      bool hasDim = obj.hasAttribute("dim");
      if (hasDim){
	// matrix
	IntegerVector ret0  = obj.attr("dim");
	return ret0[0];
      } else {
	return 1;
      }
    }
  case VECSXP:
    {
      bool rn = obj.hasAttribute("row.names");
      if (rn){
	RObject ret0 = obj.attr("row.names");
	if (rxIs(ret0, "integer")){
	  IntegerVector ret1 = IntegerVector(ret0);
	  if (ret1.size() == 2 && IntegerVector::is_na(ret1[0])){
	    return(ret1[1]);
	  } else {
	    return ret1.size();
	  }
	} else if (rxIs(ret0, "character")){
	  CharacterVector ret1 = CharacterVector(ret0);
	  return ret1.size();
	} else {
	  return NA_INTEGER;
	}
      } else {
	return NA_INTEGER;
      }
    }
  default:
    return NA_INTEGER;
  }
}

rxSolve_t rxSolveDatLast;

static inline void rxSolveSaveRxSolve(rxSolve_t* rxSolveDat){
  rxSolveDatLast = rxSolveDat[0];
  // Assigned as a precaution against R's gc()
  _rxModels[".rxSolveDat.mv"] = rxSolveDat->mv;
  _rxModels[".rxSolveDat.addDosing"] = rxSolveDat->addDosing;
  _rxModels[".rxSolveDat.timeUnitsU"] = rxSolveDat->timeUnitsU;
  _rxModels[".rxSolveDat.covUnits"] = rxSolveDat->covUnits;
  _rxModels[".rxSolveDat.par1"] = rxSolveDat->par1;
  _rxModels[".rxSolveDat.par1ini"] = rxSolveDat->par1ini;
  _rxModels[".rxSolveDat.initsC"] = rxSolveDat->initsC;
  _rxModels[".rxSolveDat.eGparPos"] = rxSolveDat->eGparPos;
  _rxModels[".rxSolveDat.sigmaN"] = rxSolveDat->sigmaN;
  _rxModels[".rxSolveDat.parNumeric"] = rxSolveDat->parNumeric;
  _rxModels[".rxSolveDat.parDf"] = rxSolveDat->parDf;
  _rxModels[".rxSolveDat.parMat"] = rxSolveDat->parMat;
  _rxModels[".rxSolveDat.nmP"] = rxSolveDat->nmP;
  _rxModels[".rxSolveDat.mvIni"] = rxSolveDat->mvIni;
  _rxModels[".rxSolveDat.idLevels"] = rxSolveDat->idLevels;
}

static inline void rxSolveSaveRxRestore(rxSolve_t* rxSolveDat){
  // Assigned as a precaution against R's gc()
  rxSolveDat->mv = as<List>(_rxModels[".rxSolveDat.mv"]);
  rxSolveDat->addDosing = as<Nullable<LogicalVector>>(_rxModels[".rxSolveDat.addDosing"]);
  rxSolveDat->timeUnitsU = as<RObject>(_rxModels[".rxSolveDat.timeUnitsU"]);
  rxSolveDat->covUnits = as<List>(_rxModels[".rxSolveDat.covUnits"]);
  rxSolveDat->par1 = as<RObject>(_rxModels[".rxSolveDat.par1"]);
  rxSolveDat->par1ini = as<RObject>(_rxModels[".rxSolveDat.par1ini"]);
  rxSolveDat->initsC = as<NumericVector>(_rxModels[".rxSolveDat.initsC"]);
  rxSolveDat->eGparPos = as<IntegerVector>(_rxModels[".rxSolveDat.eGparPos"]);
  rxSolveDat->sigmaN = as<CharacterVector>(_rxModels[".rxSolveDat.sigmaN"]);
  rxSolveDat->parNumeric = as<NumericVector>(_rxModels[".rxSolveDat.parNumeric"]);
  rxSolveDat->parDf = as<DataFrame>(_rxModels[".rxSolveDat.parDf"]);
  rxSolveDat->parMat = as<NumericMatrix>(_rxModels[".rxSolveDat.parMat"]);
  rxSolveDat->nmP = as<CharacterVector>(_rxModels[".rxSolveDat.nmP"]);
  rxSolveDat->mvIni = as<NumericVector>(_rxModels[".rxSolveDat.mvIni"]);
  rxSolveDat->idLevels = as<CharacterVector>(_rxModels[".rxSolveDat.idLevels"]);
}

RObject _curPar;
static inline bool rxSolve_loaded(const RObject &trueEvents,
				  RObject &trueParams,
				  const List &rxControl,
				  const SEXP &mv,
				  const RObject &inits){
  RObject iCov = rxControl["iCov"];
  RObject thetaMat = rxControl["thetaMat"];
  RObject omega = rxControl["omega"];
  RObject sigma = rxControl["sigma"];
  // print(trueParams);
  bool ret = _globals.gsolve != NULL &&
    as<bool>(rxControl["cacheEvent"]) &&
    _rxModels.exists(".lastEvents") &&
    // FIXME:
    // - Inductive linearization requires a different setup
    //  - Check to make sure the models are the same.
    //  - Check oldModel = newModel
    //  - Options are the same
    rxIs(iCov,"NULL") && rxIs(omega, "NULL") &&
    rxIs(sigma, "NULL") && rxIs(thetaMat, "NULL") &&
    // the #rows of input should be equal.
    getNRows(trueParams) == as<int>(_rxModels[".lastNrow"]) &&
    R_compute_identical(as<SEXP>(_rxModels[".lastEvents"]),
			as<SEXP>(trueEvents), 16) &&
    R_compute_identical(as<SEXP>(_rxModels[".lastParNames"]),
			as<SEXP>(trueParams.attr("names")), 16) &&
    R_compute_identical(as<SEXP>(_rxModels[".lastMv"]), mv, 16) &&
    R_compute_identical(as<SEXP>(_rxModels[".lastInits"]),as<SEXP>(inits), 16) &&
    R_compute_identical(as<SEXP>(_rxModels[".lastControl"]), as<SEXP>(rxControl), 16);
  if (ret && rxIs(trueParams, "data.frame")) {
    rxSolve_t* rxSolveDat = &rxSolveDatLast;
    rxSolveSaveRxRestore(rxSolveDat);
    if (rxSolveDat->idLevels.size() > 0){
      // FIXME: check to see if IDs are dropped.
      Function sortId = getRxFn(".sortId");
      // print(rxSolveDat->idLevels);
      // print(trueParams);
      RObject cur = clone(sortId(trueParams,
  				 rxSolveDat->idLevels,
  				 "parameters",
  				 false, false));
      rx_solve* rx = getRxSolve_();
      // the #subs of curPar should equal rx->nsub
      if (getNRows(cur) == rx->nsub){
  	// Rprintf("1\n");
  	_rxModels[".saveDf"] = cur;
  	// print(_rxModels[".saveDf"]);
      } else {
  	ret = false;
      }
    } else {
      // Rprintf("2\n");
      _rxModels[".saveDf"] = trueParams;
      // print(_rxModels[".saveDf"]);
    }
  }
  // if (ret) Rprintf("cached\n");
  return ret;
}

static inline void rxSolve_assignGpars(rxSolve_t* rxSolveDat){
  unsigned int i;
  for (unsigned int j = 0; j < (unsigned int)rxSolveDat->nPopPar; j++){
    for (unsigned int k = 0; k < (unsigned int)rxSolveDat->npars; k++){
      i = k+rxSolveDat->npars*j;
      if (ISNA(_globals.gpars[i])){
	if (_globals.gParPos[k] == 0){
	  _globals.gpars[i] = 0;
	} else if (_globals.gParPos[k] > 0){
	  // posPar[i] = j + 1;
	  _globals.gpars[i] = rxSolveDat->parMat(j, _globals.gParPos[k]-1);
	} else {
	  // posPar[i] = -j - 1;
	  _globals.gpars[i] = rxSolveDat->mvIni[-_globals.gParPos[k]-1];
	}
      }
    }
  }
}

// This updates the parameters based on the already loaded event table.
static inline void rxSolve_updateParams(RObject &trueParams,
					const RObject &obj,
					const List &rxControl,
					const Nullable<CharacterVector> &specParams,
					const Nullable<List> &extraArgs,
					const RObject &params,
					const RObject &events,
					const RObject &inits,
					rxSolve_t* rxSolveDat){
  rxSolveDat->par1 = trueParams;
  rxSolveDat->usePar1=true;
  rx_solve* rx = getRxSolve_();
  rx_solving_options* op = rx->op;
  if (rxIs(rxSolveDat->par1, "data.frame")){
    if (rxSolveDat->idLevels.size() > 0){
      // FIXME: check to see if IDs are dropped.
      Function sortId = getRxFn(".sortId");
      rxSolveDat->parDf = clone(sortId(rxSolveDat->par1, rxSolveDat->idLevels, "parameters",
				      rxSolveDat->warnIdSort));
      rxSolveDat->par1 = rxSolveDat->parDf;
    } else {
      rxSolveDat->parDf = as<DataFrame>(rxSolveDat->par1);
    }  
    rxSolveDat->parType = 2;
    rxSolveDat->nPopPar = rxSolveDat->parDf.nrows();
    unsigned int parDfi = rxSolveDat->parDf.size();
    rxSolveDat->parMat = NumericMatrix(rxSolveDat->nPopPar, parDfi);
    while (parDfi--){
      rxSolveDat->parMat(_,parDfi)=NumericVector(rxSolveDat->parDf[parDfi]);
    }
    rxSolveDat->parMat.attr("dimnames") = List::create(R_NilValue, rxSolveDat->parDf.names());
  } else if (rxIs(rxSolveDat->par1, "matrix")){
    rxSolveDat->parMat = as<NumericMatrix>(rxSolveDat->par1);
    rxSolveDat->nPopPar = rxSolveDat->parMat.nrow();
    rxSolveDat->parType = 3;
  } else if (!rxIs(rxSolveDat->par1, "NULL")) {
    rxSolveDat->parNumeric = rxSolveDat->par1;
    // Convert NumericVector to a matrix
    
    rxSolveDat->parMat = NumericMatrix(rx->nsub, rxSolveDat->parNumeric.size());
    unsigned int parDfi = rxSolveDat->parNumeric.size();
    while (parDfi--){
      std::fill_n(rxSolveDat->parMat.begin()+parDfi*rx->nsub, rx->nsub, rxSolveDat->parNumeric[parDfi]);
    }
    rxSolveDat->parMat.attr("dimnames") = List::create(R_NilValue, rxSolveDat->parNumeric.attr("names"));
  }
  RObject ev1=_rxModels[".lastEv1"];
  std::fill_n(&_globals.gpars[0], rxSolveDat->npars*rxSolveDat->nPopPar, NA_REAL);
  gparsCovSetupConstant(ev1, rxSolveDat->npars);
  // Setup a possibly new scale.
  RObject scale = rxControl["scale"];
  NumericVector scaleC = rxSetupScale(obj, scale, extraArgs);
  std::copy(scaleC.begin(),scaleC.end(),&_globals.gscale[0]);
  rxSolve_assignGpars(rxSolveDat);
  // std::fill_n(&_globals.gsolve[0], rx->nall*op->neq*rx->nsim, 0.0);
  seedEng(op->cores); // Reseed
}

static inline SEXP rxSolve_finalize(const RObject &obj,
				    const List &rxControl,
				    const Nullable<CharacterVector> &specParams,
				    const Nullable<List> &extraArgs,
				    const RObject &params, const RObject &events,
				    const RObject &inits,
				    rxSolve_t* rxSolveDat){
  rxSolveSaveRxSolve(rxSolveDat);
  rx_solve* rx = getRxSolve_();
  par_solve(rx);
  List dat = rxSolve_df(obj, rxControl, specParams, extraArgs,
			params, events, inits, rxSolveDat);
  if (rx->matrix == 0 && rxDropB){
    rx->matrix=2;
    warning(_("dropped key column, returning data.frame instead of special solved data.frame"));
  }
  if (rx->matrix){
    // rxSolveFree();
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
      // Free(op->indLin);
      return dat;
    } else {
      NumericMatrix tmpM(rx->nr, dat.size());
      for (unsigned int i = dat.size(); i--;){
	tmpM(_,i) = as<NumericVector>(dat[i]);
      }
      tmpM.attr("dimnames") = List::create(R_NilValue,dat.names());
      // Free(op->indLin);
      return tmpM;
    }
  } else {
    CharacterVector cls= CharacterVector::create("rxSolve", "rxSolveParams","rxSolveCovs",
						 "rxSolveInits", "rxSolveSimType",
						 "data.frame");
    cls.attr(".RxODE.env") = rxSolve_genenv(obj, rxControl, specParams,
					    dat, params, events,
					    inits, rxSolveDat);
    // eGparPos
    dat.attr("class") = cls;
    // rxSolveFree();
    // Free(op->indLin);
    return(dat);
  }    
}

//[[Rcpp::export]]
SEXP rxSolve_(const RObject &obj, const List &rxControl,
	      const Nullable<CharacterVector> &specParams,
	      const Nullable<List> &extraArgs,
	      const RObject &params, const RObject &events, const RObject &inits,
	      const int setupOnly){
  rxDropB = false;
#ifdef rxSolveT
  clock_t _lastT0 = clock();
#endif
  if (rxIs(rxControl,"rxControl")){
    rxSolveFree();
    stop(_("control list not setup correctly"));
  }
  maxAtolRtolFactor = as<double>(rxControl["maxAtolRtolFactor"]);
  _mvnfast = as<bool>(rxControl["mvnfast"]);
  RObject scale = rxControl["scale"];
  int method = as<int>(rxControl["method"]);
  Nullable<LogicalVector> transit_abs = rxControl["transitAbs"];
  NumericVector atolNV = as<NumericVector>(rxControl["atol"]);
  NumericVector rtolNV = as<NumericVector>(rxControl["rtol"]);
  NumericVector atolNVss = as<NumericVector>(rxControl["ssAtol"]);
  NumericVector rtolNVss = as<NumericVector>(rxControl["ssRtol"]);
  int maxsteps = as<int>(rxControl["maxsteps"]);
  double hmin = as<double>(rxControl["hmin"]);
  double hini = as<double>(rxControl["hini"]);
  int maxordn = as<int>(rxControl["maxordn"]);
  int maxords = as<int>(rxControl["maxords"]);
  unsigned int cores = as<unsigned int>(rxControl["cores"]);
  int covs_interpolation = as<int>(rxControl["covsInterpolation"]);
  bool addCov = as<bool>(rxControl["addCov"]);
  int matrix = as<int>(rxControl["matrix"]);
  int nDisplayProgress = as<int>(rxControl["nDisplayProgress"]);
  NumericVector stateTrim = rxControl["stateTrim"];
  double stateTrimU= R_PosInf;
  double stateTrimL= R_NegInf;
  if (stateTrim.size() == 2){
    if (stateTrim[0] > stateTrim[1]){
      stateTrimU = stateTrim[0];
      stateTrimL = stateTrim[1];
    } else {
      stateTrimU = stateTrim[1];
      stateTrimL = stateTrim[0];
    }
  } else if (stateTrim.size() == 1){
    if (!ISNA(stateTrimU)){
      stateTrimU = fabs(stateTrim[0]);
      stateTrimL = -stateTrimU;
    }
  } else {
    stop("'stateTrim' must be a vector of 1-2 elements");
  }
  rxSolve_t rxSolveDat0;
  rxSolve_t* rxSolveDat = &rxSolveDat0;
  RObject object;
  rxSolveDat->updateObject = as<bool>(rxControl["updateObject"]);  
  rxSolveDat->isRxSolve = rxIs(obj, "rxSolve");
  rxSolveDat->isEnvironment = rxIs(obj, "environment");
  rxSolveDat->idFactor= as<bool>(rxControl["idFactor"]);
  rxSolveDat->warnIdSort = as<bool>(rxControl["warnIdSort"]);
  if (rxSolveDat->updateObject && !rxSolveDat->isRxSolve && !rxSolveDat->isEnvironment){
    if (rxIs(rxCurObj, "rxSolve")){
      object = rxCurObj;
      rxSolveDat->isRxSolve = true;
    } else {
      rxSolveFree();
      stop(_("cannot update this object"));
    }
  } else {
    if (method == 3){
      rxSolveDat->mv = rxModelVars(obj);
      List indLin = rxSolveDat->mv["indLin"];
      if (indLin.size() == 0){
	Function RxODE = getRxFn("RxODE");
	object = RxODE(obj, _["indLin"]=true);
	rxSolveDat->mv = rxModelVars(object);
      } else {
	object =obj;
      }
    } else {
      object =obj;
      rxSolveDat->mv = rxModelVars(object);
    }
  }
  if (rxSolveDat->isRxSolve || rxSolveDat->isEnvironment){
    return rxSolve_update(object, rxControl, specParams,
			  extraArgs, params, events, inits,
			  rxSolveDat);
  } else {
#ifdef rxSolveT
    REprintf("Time1: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif// rxSolveT
    // Load model
    if (!rxDynLoad(object)){
      rxSolveFree();
      stop(_("cannot load RxODE dlls for this model"));
    }
    // Get model 
    RObject trueParams;
    RObject trueEvents;
    bool swappedEvents=false;
    getRxModels();
    // Get the C solve object
    rx_solve* rx = getRxSolve_();
    rx_solving_options* op = rx->op;    
    if (rxIs(params, "rx.event")){
      trueEvents = params;
      trueParams = events;
      swappedEvents=true;
    } else if (rxIs(events, "rx.event")){
      trueEvents = events;
      trueParams = params;
    } else {
      stop(_("cannot solve without event information"));
    }
#ifdef rxSolveT
    REprintf("Time2: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif// rxSolveT
    if (setupOnly == 0 &&
	rxSolve_loaded(trueEvents, trueParams, rxControl, as<SEXP>(rxSolveDat->mv), inits)){
      rxSolveDat = &rxSolveDatLast;
      rxSolve_updateParams(trueParams,
			   obj, rxControl, specParams, extraArgs, params, events,
			   inits, rxSolveDat);
      if (_rxModels.exists(".ws") && !rxIs(_rxModels[".ws"], "NULL")){
	List ws = _rxModels[".ws"];
	for (int iii = 0; iii < ws.size(); ++iii){
	  warning(as<std::string>(ws[iii]));
	}
      } 
      return rxSolve_finalize(obj, rxControl, specParams, extraArgs, params, events,
			      inits, rxSolveDat);
      // par_solve(rx);
      // Free(op->indLin);
      // rxSolveFree();
    } else {
      _rxModels[".lastEvents"] = trueEvents;
      _rxModels[".lastParNames"] = trueParams.attr("names");
      _rxModels[".lastNrow"] = getNRows(trueParams);
      _rxModels[".lastMv"] = rxSolveDat->mv;
      _rxModels[".lastControl"] = rxControl;
      _rxModels[".lastInits"] = inits;
      Free(op->indLin);
    }
    rxSolveDat->addDosing = as<Nullable<LogicalVector>>(rxControl["addDosing"]);
#ifdef rxSolveT
    REprintf("Time3: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif// rxSolveT
    CharacterVector pars = rxSolveDat->mv["params"];
    rxSolveDat->npars = pars.size();
    rxSolveDat->hasCmt = false;
    for (int i = rxSolveDat->npars; i--;){
      if (as<std::string>(pars[i]) == "CMT"){
	rxSolveDat->hasCmt=true;
	break;
      }
    }
    // Assign Pointers
    rxAssignPtr(rxSolveDat->mv);
    rx->nKeep0 = 0;
    rx->nKeepF = 0;
    rx->stateTrimU = stateTrimU;
    rx->stateTrimL = stateTrimL;
    rx->matrix = matrix;
    rx->needSort = as<int>(rxSolveDat->mv["needSort"]);
    rx->nMtime = as<int>(rxSolveDat->mv["nMtime"]);
    rx->add_cov = (int)(addCov);
    rx->istateReset = as<int>(rxControl["istateReset"]);
    rx->safeZero = as<int>(rxControl["safeZero"]);
    op->stiff = method;
    op->linLog=as<int>(rxControl["linLog"]);
    op->advanLinCmt = as<int>(rxControl["advanLinCmt"]);
    if (method != 2){
      op->cores =1;
    } else {
      op->cores=cores;
    }
    seedEng(op->cores);
    // Now set up events and parameters
    RObject par0 = params;
    RObject ev0  = events;
    RObject ev1;
    rxSolveDat->swappedEvents = false;
    rxSolveDat->nsvar = 0;
    rxSolveDat->labelID=false;
    rxSolveDat->swappedEvents=swappedEvents;
    ev1 = trueEvents;
    rxSolveDat->par1 = trueParams;
    if (rxIs(rxSolveDat->par1, "NULL")){
      rxSolveDat->par1=rxInits(obj);
      rxSolveDat->fromIni=true;
    }
#ifdef rxSolveT
    REprintf("Time5: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif// rxSolveT

    // Update event table with observations if they are missing
    rxSolve_ev1Update(obj, rxControl, specParams, extraArgs, params,
		      ev1, inits, rxSolveDat);

#ifdef rxSolveT
    REprintf("Time6: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif// rxSolveT
    // Now get the parameters (and covariates)
    //
    // Unspecified parameters can be found in the modVars["ini"]
    rxSolveDat->mvIni = rxSolveDat->mv["ini"];
    // The event table can contain covariate information, if it is acutally a data frame or matrix.
    Nullable<CharacterVector> covnames0, simnames0;
    CharacterVector covnames, simnames;
    rxSolveDat->eGparPos = IntegerVector(rxSolveDat->npars);
    CharacterVector state = rxSolveDat->mv["state"];
    CharacterVector lhs = rxSolveDat->mv["lhs"];
    op->neq = state.size();
    op->badSolve = 0;
    op->abort = 0;
    op->ATOL = atolNV[0];          //absolute error
    op->RTOL = rtolNV[0];          //relative error

    op->minSS = as<int>(rxControl["minSS"]);
    op->maxSS = as<int>(rxControl["maxSS"]);
    op->infSSstep = as<double>(rxControl["infSSstep"]);
    if (op->infSSstep <= 0) stop(_("'infSSstep' needs to be positive"));
    op->indLinPhiTol=as<double>(rxControl["indLinPhiTol"]);
    op->indLinMatExpType=as<int>(rxControl["indLinMatExpType"]);
    op->indLinPhiM = as<int>(rxControl["indLinPhiM"]);
    op->indLinMatExpOrder=as<int>(rxControl["indLinMatExpOrder"]);
    List indLin = rxSolveDat->mv["indLin"];
    op->doIndLin=0;
    if (indLin.size() == 4){
      int me = rxIs(indLin[1], "NULL");
      if (as<bool>(indLin[2])){
	// Inductive linearization
	IntegerVector indLinItems = as<IntegerVector>(indLin[3]);
	op->indLinN = indLinItems.size();
	op->indLin = Calloc(op->indLinN,int);
	std::copy(indLinItems.begin(), indLinItems.end(), op->indLin);
	if (me){
	  // homogenous ME + IndLin
	  op->doIndLin=3;
	} else {
	  // inhomogenous ME + IndLin
	  op->doIndLin=4;
	}
      } else {
	op->indLinN = 0;
	op->indLin = NULL;
	if (me){
	  op->doIndLin=1;
	} else {
	  op->doIndLin=2;
	}
      }
    }
    
    op->H0 = hini;
    op->HMIN = hmin;
    op->mxstep = maxsteps;
    op->MXORDN = maxordn;
    op->MXORDS = maxords;
    // The initial conditions cannot be changed for each individual; If
    // they do they need to be a parameter.
    int transit = 0;
    if (transit_abs.isNull()){
      transit = rxSolveDat->mv["podo"];
      if (transit){
        warning(_("assumed transit compartment model since 'podo' is in the model"));
      }
    }  else {
      LogicalVector tr = LogicalVector(transit_abs);
      if (tr[0]){
        transit=  1;
      }
    }
    op->do_transit_abs = transit;
    op->nlhs = lhs.size();
    CharacterVector trans = rxSolveDat->mv["trans"];
    // Make sure the model variables are assigned...
    // This fixes random issues on windows where the solves are done and the data set cannot be solved.
    getRxModels();
    _rxModels[as<std::string>(trans["model_vars"])] = rxSolveDat->mv;
    sprintf(op->modNamePtr, "%s", (as<std::string>(trans["model_vars"])).c_str());
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
      rxSolveFree();
      stop(_("unknown covariate interpolation specified"));
    }
    op->extraCmt=op->neq+as<int>(rxSolveDat->mv["extraCmt"]);
    op->nDisplayProgress = nDisplayProgress;
    // Covariate options
    // Simulation variabiles
    // int *svar;
    rxSolveDat->usePar1 = false;
    rxSolveDat->addTimeUnits = false;
    if (rxIs(ev1, "rxEtTran")){
      CharacterVector cls = ev1.attr("class");
      List evT = cls.attr(".RxODE.lst");
      evT.attr("class") = R_NilValue;
      rxSolveDat->covUnits = evT["covUnits"];
    }
    rxSolveDat->par1ini = rxSolveDat->par1;
    // This will update par1 with simulated values
#ifdef rxSolveT
    REprintf("Time7: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif// rxSolveT
    rxSolve_simulate(obj, rxControl, specParams, extraArgs,
		     params, ev1, inits, rxSolveDat);
#ifdef rxSolveT
    REprintf("Time8: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif// rxSolveT
    // This will setup the parameters
    rxSolve_parSetup(obj, rxControl, specParams, extraArgs,
		     pars, ev1, inits, rxSolveDat);
#ifdef rxSolveT
    REprintf("Time9: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif// rxSolveT
    rxOptionsIniEnsure(rxSolveDat->nPopPar); // 1 simulation per parameter specification


    // Setup some data-based parameters like hmax
    rxSolve_datSetupHmax(obj, rxControl, specParams, extraArgs,
			 pars, ev1, inits, rxSolveDat);

#ifdef rxSolveT
    REprintf("Time10: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif // rxSolveT
    
    // Make sure the user input all the parameters.
    rxSolve_parOrder(obj, rxControl, specParams, extraArgs,
		     pars, ev1, inits, rxSolveDat);

#ifdef rxSolveT
    REprintf("Time11: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif // rxSolveT
    
    op->svar = &_globals.gsvar[0];
    op->nsvar = rxSolveDat->nsvar;
    if (op->nsvar == 0){
      getRxModels();
      if(_rxModels.exists(".sigma")){
	_rxModels.remove(".sigma");
      }
    }
    // Now setup the rest of the rx_solve object
    if (rxSolveDat->nPopPar != 1 && rxSolveDat->nPopPar % rx->nsub != 0){
      rxSolveFree();
      stop(_("number of parameters (%d) solved by RxODE for multi-subject data needs to be a multiple of the number of subjects (%d)"),rxSolveDat->nPopPar, rx->nsub);
    }
#ifdef rxSolveT
    REprintf("Time12: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif // rxSolveT
    rxSolveDat->nSize = rxSolveDat->nPopPar*rx->nsub;
    if (rxSolveDat->nPopPar % rx->nsub == 0) rx->nsim = rxSolveDat->nPopPar / rx->nsub;
    else rx->nsim=1;
    IntegerVector linCmtI = rxSolveDat->mv["flags"];
    int linNcmt = linCmtI[0];
    int linKa = linCmtI[1];
    op->nlin = linNcmt+linKa;
    int n0 = (rx->nall+3*rx->nsub)*(state.size())*rx->nsim;
    int nLin = linNcmt+linKa;
    if (nLin != 0) nLin = rx->nall*(nLin+2)*rx->nsim;
    if (op->advanLinCmt == 0) nLin = 0;
    int n2 = rx->nMtime*rx->nsub*rx->nsim;
    int n3 = op->neq*rxSolveDat->nSize;
#ifdef rxSolveT
    REprintf("Time12a: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif // rxSolveT

    rxSolveDat->initsC = rxInits(object, inits, state, 0.0);

#ifdef rxSolveT
    REprintf("Time12b: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif // rxSolveT


    int n4 = rxSolveDat->initsC.size();
    int n5 = lhs.size()*rxSolveDat->nSize;
    // The initial conditions cannot be changed for each individual; If
    // they do they need to be a parameter.
    NumericVector scaleC = rxSetupScale(object, scale, extraArgs);
    int n6 = scaleC.size();
    if (_globals.gsolve != NULL) free(_globals.gsolve);
    _globals.gsolve = (double*)calloc(n0+nLin+n2+n3+n4+n5+n6+ 4*op->neq, sizeof(double));// [n0]
    _globals.gadvan = _globals.gsolve+n0; // [nLin]
#ifdef rxSolveT
    REprintf("Time12c (double alloc %d): %f\n",n0+nLin+n2+n3+n4+n5+n6+ 4*op->neq,((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif // rxSolveT
    if (_globals.gsolve == NULL){
      rxSolveFree();
      stop(_("could not allocate enough memory for solving"));
    }
    _globals.gmtime = _globals.gadvan +nLin; // [n2]
    _globals.gInfusionRate = _globals.gmtime + n2; //[n3]
    _globals.ginits = _globals.gInfusionRate + n3; // [n4]
    std::copy(rxSolveDat->initsC.begin(), rxSolveDat->initsC.end(), &_globals.ginits[0]);
    op->inits = &_globals.ginits[0];
    _globals.glhs = _globals.ginits + n4; // [n5]
    // initially NA_REAL
    //std::fill_n(_globals.glhs,n5, NA_REAL); // TOO slow
    _globals.gscale = _globals.glhs + n5; //[n6]
    std::copy(scaleC.begin(),scaleC.end(),&_globals.gscale[0]);
    op->scale = &_globals.gscale[0];
    _globals.gatol2=_globals.gscale   + n6; //[op->neq]
    _globals.grtol2=_globals.gatol2   + op->neq;  //[op->neq]
    _globals.gssRtol=_globals.grtol2  + op->neq; //[op->neq]
    _globals.gssAtol=_globals.gssRtol + op->neq; //[op->neq]
    std::fill_n(&_globals.gatol2[0],op->neq, atolNV[0]);
    std::fill_n(&_globals.grtol2[0],op->neq, rtolNV[0]);
    std::copy(atolNV.begin(), atolNV.begin() + min2(op->neq, atolNV.size()), &_globals.gatol2[0]);
    std::copy(rtolNV.begin(), rtolNV.begin() + min2(op->neq, rtolNV.size()), &_globals.grtol2[0]);

    std::fill_n(&_globals.gssAtol[0],op->neq, atolNVss[0]);
    std::fill_n(&_globals.gssRtol[0],op->neq, rtolNVss[0]);
    std::copy(atolNVss.begin(), atolNVss.begin() + min2(op->neq, atolNVss.size()), &_globals.gssAtol[0]);
    std::copy(rtolNVss.begin(), rtolNVss.begin() + min2(op->neq, rtolNVss.size()), &_globals.gssRtol[0]);
    op->atol2 = &_globals.gatol2[0];
    op->rtol2 = &_globals.grtol2[0];
    op->ssAtol = _globals.gssAtol;
    op->ssRtol = _globals.gssRtol;
    // Not needed since we use Calloc.
    // std::fill_n(&_globals.gsolve[0], rx->nall*state.size()*rx->nsim, 0.0);
    int n1 = rx->nsub*rx->nsim*state.size();
#ifdef rxSolveT
    REprintf("Time12d (fill in!): %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif // rxSolveT
    if (_globals.gon != NULL) free(_globals.gon);
    _globals.gon = (int*)calloc(n1+n3 + 4*rxSolveDat->nSize + rx->nall*rx->nsim, sizeof(int)); // [n1]
#ifdef rxSolveT
    REprintf("Time12e (int alloc %d): %f\n", n1+n3 + 4*rxSolveDat->nSize + rx->nall*rx->nsim, ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif // rxSolveT
    std::fill_n(&_globals.gon[0], n1, 1);
    _globals.gBadDose = _globals.gon+n1; // [n3]
    _globals.grc = _globals.gBadDose + n3; //[nSize]
    _globals.slvr_counter = _globals.grc + rxSolveDat->nSize; //[nSize]
    _globals.dadt_counter = _globals.slvr_counter + rxSolveDat->nSize; // [nSize]
    _globals.jac_counter = _globals.dadt_counter + rxSolveDat->nSize; // [nSize]
    _globals.gix=_globals.jac_counter+rxSolveDat->nSize; // rx->nall*rx->nsim

#ifdef rxSolveT
    REprintf("Time13: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif // rxSolveT
    rxSolve_normalizeParms(obj, rxControl, specParams, extraArgs,
			   pars, ev1, inits, rxSolveDat);
    
    if (setupOnly){
      return as<SEXP>(LogicalVector::create(true));
    }
#ifdef rxSolveT
    REprintf("Time14: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif// rxSolveT
    SEXP ret = rxSolve_finalize(obj, rxControl, specParams, extraArgs, params, events,
				inits, rxSolveDat);
#ifdef rxSolveT
    REprintf("Time15: %f\n", ((double)(clock() - _lastT0))/CLOCKS_PER_SEC);
    _lastT0 = clock();
#endif // rxSolveT
    return ret;
  }
  return R_NilValue;
}


RObject rxSolveGet_rxSolve(RObject &obj, std::string &sarg, LogicalVector &exact,
			   List &lst){
  int i, j, n;
  rxCurObj = obj;
  CharacterVector cls = lst.attr("class");
  Environment e = as<Environment>(cls.attr(".RxODE.env"));
  if (sarg == "env"){
    return as<RObject>(e);
  }
  if (sarg == "model"){
    List mv = rxModelVars(obj);
    CharacterVector mods = mv["model"];
    CharacterVector retS = as<std::string>(mods["normModel"]);
    retS.attr("class") = "rxModelText";
    return(retS);
  }
  updateSolveEnvPost(e);
  if (e.exists(sarg)){
    return e[sarg];
  }
  if (sarg == "params" || sarg == "par" || sarg == "pars" || sarg == "param"){
    List ret = clone(as<List>(e[".params.dat"]));
    return ret;
  } else if (sarg == "inits" || sarg == "init"){
    NumericVector ret = clone(as<NumericVector>(e[".init.dat"]));
    return ret;
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
  List pars = clone(List(e[".params.dat"]));
  CharacterVector nmp = pars.names();
  n = pars.size();
  for (i = n; i--;){
    if (nmp[i] == sarg){
      return pars[sarg];
    }
  }
  // // Now inis.
  // Function sub("sub", R_BaseNamespace);
  NumericVector ini = clone(NumericVector(e[".init.dat"]));
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
  return R_NilValue;
}

//[[Rcpp::export]]
CharacterVector rxSolveDollarNames(RObject obj){
  CharacterVector names1 = obj.attr("names");
  List lst = as<List>(obj);
  CharacterVector cls = lst.attr("class");
  Environment e = as<Environment>(cls.attr(".RxODE.env"));
  updateSolveEnvPost(e);
  int nExtra = 6;
  if (e.exists(".theta")) nExtra++;
  if (e.exists(".sigmaL")) nExtra++;
  if (e.exists(".omegaL")) nExtra++;
  
  List pars = List(e[".params.dat"]);
  CharacterVector pn = pars.attr("names");
  
  NumericVector ini = NumericVector(e[".init.dat"]);
  CharacterVector nmi = ini.names();

  CharacterVector envl = e.ls(false);

  CharacterVector ret(nExtra + names1.size() +
		      pn.size() + nmi.size() + envl.size());
  
  int j = 0;
  for (int i = names1.size(); i--;){
    ret[j++] = names1[i];
  }

  for (int i = pn.size(); i--;){
    ret[j++] = pn[i];
  }

  for (int i = nmi.size(); i--;){
    ret[j++] = as<std::string>(nmi[i]) + "0";
  }

  for (int i = envl.size(); i--;){
    ret[j++] = envl[i];
  }
  
  ret[j++] = "env";
  ret[j++] = "model";
  ret[j++] = "params";
  ret[j++] = "inits";
  ret[j++] = "t";
  ret[j++] = "rxode";
  if (e.exists(".theta")) ret[j++] = "thetaMat";
  if (e.exists(".sigmaL")) ret[j++] = "sigmaList";
  if (e.exists(".omegaL")) ret[j++] = "omegaList";
  return ret;
}

//[[Rcpp::export]]
RObject rxSolveGet(RObject obj, RObject arg, LogicalVector exact = true){
  std::string sarg;
  int i, n;
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
	    warning(_("partial match of '%s' to '%s'"),sarg.c_str(), (as<std::string>(nm[i])).c_str());
	  }
	  return lst[i];
	}
      }
      if (rxIs(obj, "rxSolve")){
	RObject ret0 = rxSolveGet_rxSolve(obj, sarg, exact, lst);
	if (!rxIs(ret0, "NULL")){
	  return ret0;
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
      CharacterVector classattr = obj.attr("class");
      Environment e = as<Environment>(classattr.attr(".RxODE.env"));
      List rxControl = e[".args"];
      // updateSolveEnvPost(e);
      if (what.size() == 1){
	std::string sarg = as<std::string>(what[0]);
	// Now check to see if this is something that can be updated...
	if (sarg == "params"){
	  // rxControl["params"] = value;
	  return rxSolve_(obj,rxControl,
                          CharacterVector::create("params"),
			  R_NilValue,
                          value, //defrx_params,
                          defrx_events,
                          defrx_inits, 0);
	} else if (sarg == "events"){
	  return rxSolve_(obj,rxControl,
			  CharacterVector::create("events"),
			  R_NilValue,
			  defrx_params,
			  value, // defrx_events,
			  defrx_inits, 0);
	} else if (sarg == "inits"){
	  return rxSolve_(obj, rxControl,
                          CharacterVector::create("inits"),
			  R_NilValue,
                          defrx_params,
                          defrx_events,
                          as<RObject>(value), //defrx_inits,
			  0);
	} else if (sarg == "t" || sarg == "time"){
	  CharacterVector classattr = obj.attr("class");
          Environment e = as<Environment>(classattr.attr(".RxODE.env"));
          updateSolveEnvPost(e);
	  Function f = as<Function>(e[".replace.sampling"]);
	  return f(value);
        } else {
	  CharacterVector classattr = obj.attr("class");
	  Environment e = as<Environment>(classattr.attr(".RxODE.env"));
          updateSolveEnvPost(e);
          if (rxIs(e[".params.dat"], "NULL")) {
	    rxSolveFree();
	    stop(_("cannot update nonexistent parameters"));
	  }
          List pars = List(e[".params.dat"]);
	  CharacterVector nmp = pars.names();
	  int i, n, np, nc, j;
	  np = (as<NumericVector>(pars[0])).size();
	  List events = List(e[".args.events"]);
	  CharacterVector nmc;
	  nmc = events.names();
	  nc = (as<NumericVector>(events[0])).size();
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
		return rxSolve_(obj,rxControl,
				CharacterVector::create("params"),
				R_NilValue,
				pars, //defrx_params,
				defrx_events,
				defrx_inits, 0);
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
		List newEvents(events.size()+1);
		CharacterVector newEventsNames(events.size()+1);
		for (j = events.size(); j--;){
		  newEvents[j]      = events[j];
		  newEventsNames[j] = nmc[j];
		}
		j = events.size();
		newEvents[j]      = val;
		newEventsNames[j] = nmp[i];
		newEvents.attr("names") = newEventsNames;
		newEvents.attr("class") = "data.frame";
		newEvents.attr("row.names") = IntegerVector::create(NA_INTEGER,-nc);
		return rxSolve_(obj, rxControl,
				CharacterVector::create("params","events"),
				R_NilValue,
				newPars, //defrx_params,
				newEvents,
				defrx_inits, 0);
	      }
	      return R_NilValue;
	    }
	  }
	  ///////////////////////////////////////////////////////////////////////////////
	  // Update Covariates by covariate name
	  n = events.size();
	  for (i = n; i--;){
	    if (nmc[i] == sarg){
	      // Update solve.
	      NumericVector val = NumericVector(value);
	      if (val.size() == nc){
		// Update Covariate
		events[i]=val;
		return rxSolve_(obj, rxControl,
				CharacterVector::create("events"),
				R_NilValue,
				defrx_params,
				events,
				defrx_inits, 0);
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
		List newEvents(events.size()-1);
		CharacterVector newEventsNames(events.size()-1);
		for (j = 0; j < i; j++){
		  newEvents[j]      = events[j];
		  newEventsNames[j] = nmc[j];
		}
		for (j=i+1; j < events.size(); j++){
		  newEvents[j-1]      = events[j];
		  newEventsNames[j-1] = nmc[j];
		}
		newEvents.attr("names") = newEventsNames;
		newEvents.attr("class") = "data.frame";
		newEvents.attr("row.names") = IntegerVector::create(NA_INTEGER,-nc);
		return rxSolve_(obj,rxControl,
				CharacterVector::create("events", "params"),
				R_NilValue,
				newPars,//defrx_params,
				newEvents,
				defrx_inits, 0);
	      }
	    }
	  }
	  ////////////////////////////////////////////////////////////////////////////////
          // Update Initial Conditions
	  NumericVector ini = NumericVector(e[".init.dat"]);
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
              NumericVector ini = NumericVector(e[".init.dat"]);
	      double v = as<double>(value);
	      for (j = 0; j < n; j++){
		if (cur == as<std::string>(nmi[j])){
		  ini[j] = v;
		}
	      }
              return rxSolve_(obj, rxControl,
			      CharacterVector::create("inits"),
			      R_NilValue,
			      defrx_params,
			      defrx_events,
			      ini, 0);
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
    rxSolveFree();
    stop(_("Can not figure out the RxODE object"));
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
    Nullable<CharacterVector>pkgn = e["package"];
    if (!pkgn.isNull()){
      Function rxPkgDll = getRxFn(".rxPkgDll");
      return(as<std::string>(rxPkgDll(wrap(obj))));
    } 
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
      rxSolveFree();
      stop(_("Can not figure out the DLL for this object"));
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
      rxSolveFree();
      stop(_("can not figure out the DLL for this object"));
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
  std::string dydt = as<std::string>(trans["model_vars"]);
  bool ret = as<bool>(isLoaded(dydt));
  return ret;
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

//' Lock/unlocking of RxODE dll file
//'
//' @param obj A RxODE family of objects
//' @export
//[[Rcpp::export]]
RObject rxLock(RObject obj){
  getRxModels();
  std::string file = rxDll(obj);
  int ret = 1;
  if (_rxModels.exists(file)){
    ret = as<int>(_rxModels[file]);
    ret = ret+1;
    _rxModels[file] = ret;
  } else {
    _rxModels[file] = ret;
  }
  return R_NilValue;
}

//' @rdname rxLock
//' @export
//[[Rcpp::export]]
RObject rxUnlock(RObject obj){
  getRxModels();
  std::string file = rxDll(obj);
  int ret;
  if (_rxModels.exists(file)){
    ret = as<int>(_rxModels[file]);
    ret = ret - 1;
    if (ret > 0) _rxModels[file] = ret;
    else  _rxModels[file] = 0;
  }
  return R_NilValue;
}

bool rxCanUnload(RObject obj){
  getRxModels();
  std::string file = rxDll(obj);
  if(_rxModels.exists(file)){
    int ret = as<int>(_rxModels[file]);
    return (ret == 0L);
  }
  return true;
}

void rmRxModelsFromDll(std::string str){
  Function getInfo = getRxFn(".rxGetModelInfoFromDll");
  CharacterVector extra = getInfo(str);
  for (int j = extra.size(); j--;){
    if (_rxModels.exists(as<std::string>(extra[j]))){
      _rxModels.remove(as<std::string>(extra[j]));
    }
  }
  if (_rxModels.exists(str)){
    _rxModels.remove(str);
  }
}
bool rxUnload_ = true;
//' Allow unloading of dlls
//'
//' @param allow boolean indicating if garbage collection will unload of RxODE dlls.
//'
//' @examples
//'
//' # Garbage collection will not unload un-used RxODE dlls
//' rxAllowUnload(FALSE);
//'
//' # Garbage collection will unload unused RxODE dlls
//' rxAllowUnload(TRUE);
//' @export
//' @author Matthew Fidler
//[[Rcpp::export]]
bool rxAllowUnload(bool allow){
  rxUnload_=allow;
  return rxUnload_;
}


//' Unload all RxODE Dlls that are not locked for solving.
//' @return NULL
//' @export
//[[Rcpp::export]]
RObject rxUnloadAll(){
  getRxModels();
  Function dynUnload("dyn.unload", R_BaseNamespace);
  CharacterVector vars = _rxModels.ls(true);
  std::string exclude = ".rxSolveDat.";
  for (int i = vars.size(); i--;){
    std::string varC = as<std::string>(vars[i]);
    if (varC.find(exclude) == std::string::npos){
      if (rxIs(_rxModels[varC],"integer")){
	int val = as<int>(_rxModels[varC]);
	if (val > 1){
	} else if (val == 0 && rxUnload_){
	  dynUnload(varC);
	  rmRxModelsFromDll(varC);
	}
      }
    }
  }
  return R_NilValue;
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
  if (!rxUnload_) return false;
  if (rxIs(obj, "RxODE")){
    Environment e = as<Environment>(obj);
    Nullable<CharacterVector> pkg = e["package"];
    if (!pkg.isNull()){
      std::string sobj = as<std::string>(e["modName"]);
      if (sobj.find("_new")==std::string::npos){
	rxSolveFree();
	stop(_("package-based models cannot be unloaded"));
      }
    }
  }
  List mv = rxModelVars(obj);
  CharacterVector trans = mv["trans"];
  std::string ptr = as<std::string>(trans["model_vars"]);
  if (rxIsLoaded(obj)){
    Function dynUnload("dyn.unload", R_BaseNamespace);
    std::string file = rxDll(obj);
    rxUnlock(obj);
    if (rxCanUnload(obj)){
      dynUnload(file);
    } else {
      rxLock(obj);
      return false;
    }
  }
  rxRmModelLib_(ptr); // Clears all pointers
  std::string file = rxDll(obj);
  getRxModels();
  if(_rxModels.exists(file)){
    _rxModels.remove(file);
  }
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
  if (rxIs(obj, "RxODE")){
    Environment e = as<Environment>(obj);
    Nullable<CharacterVector> pkg = e["package"];
    if (!pkg.isNull()){
      std::string sobj = as<std::string>(e["modName"]);
      if (sobj.find("_new")==std::string::npos){
	rxSolveFree();
	stop(_("package-based models cannot be deleted"));
      }
    }
  }
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


int isProgSupportedI = 1;

//[[Rcpp::export]]
SEXP setProgSupported(int isSupported=1){
  isProgSupportedI=isSupported;
  return wrap(isProgSupportedI);
}

extern "C" int isProgSupported(){
  return isProgSupportedI;
}

//[[Rcpp::export]]
SEXP getProgSupported(){
  return wrap(isProgSupportedI);
}

//[[Rcpp::export]]
List rxUpdateTrans_(List ret, std::string prefix, std::string libName){
  CharacterVector oldTrans = as<CharacterVector>(ret["trans"]);
  CharacterVector newLst(22);
  CharacterVector newNme(22);
  newNme[0] ="lib.name";
  newLst[0] = libName;
  
  newNme[1] = "jac";
  newLst[1] = as<std::string>(oldTrans[1]);
  
  newNme[2] = "prefix";
  newLst[2] = prefix;
  
  newNme[3] = "dydt";
  newLst[3] = prefix + "dydt";
  
  newNme[4] = "calc_jac";
  newLst[4] = prefix + "calc_jac";
  
  newNme[5] = "calc_lhs";
  newLst[5] = prefix + "calc_lhs";
  
  newNme[6] = "model_vars";
  newLst[6] = prefix + "model_vars";
  
  newNme[7] = "theta";
  newLst[7] = prefix + "theta";
  
  newNme[8]="inis";
  newLst[8]= prefix + "inis";

  newNme[9]="dydt_lsoda";
  newLst[9]= prefix + "dydt_lsoda";
  
  newNme[10]="calc_jac_lsoda";
  newLst[10]= prefix + "calc_jac_lsoda";
  
  newNme[11]= "ode_solver_solvedata";
  newLst[11]= prefix + "ode_solver_solvedata";

  newNme[12] = "ode_solver_get_solvedata";
  newLst[12] = prefix+"ode_solver_get_solvedata";
  
  newNme[13]= "dydt_liblsoda";
  newLst[13]= prefix + "dydt_liblsoda";
  
  newNme[14]="F";
  newLst[14]= prefix + "F";

  newNme[15] ="Lag";
  newLst[15] = prefix + "Lag";

  newNme[16]= "Rate";
  newLst[16] = prefix + "Rate";
  
  newNme[17] ="Dur";
  newLst[17] = prefix + "Dur";
  
  newNme[18] = "mtime";
  newLst[18] = prefix + "mtime";
  
  newNme[19] = "assignFuns";
  newLst[19] = prefix + "assignFuns";

  newNme[20] = "ME";
  newLst[20] = prefix + "ME";
  
  newNme[21] = "IndF";
  newLst[21]= prefix + "IndF";

  newLst.attr("names") = newNme;
				   
  ret[3] = newLst;
  return(ret);
}

extern "C" {
  double getTime(int idx, rx_solving_options_ind *ind);
}

extern "C" void getWh(int evid, int *wh, int *cmt, int *wh100, int *whI, int *wh0);
extern "C" void doSort(rx_solving_options_ind *ind){
  // Reset indexes
  std::iota(&(ind->ix[0]),&(ind->ix[0])+ind->n_all_times, 0);
  // Reset times for infusion
  int wh, cmt, wh100, whI, wh0;
  for (int j = ind->n_all_times; j--;){
    getWh(ind->evid[j], &wh, &cmt, &wh100, &whI, &wh0);
    if (whI == 6 || whI == 7){
      ind->all_times[j] = ind->all_times[j-1];
    }
  }
  try {
    SORT(&(ind->ix[0]),&(ind->ix[0])+ind->n_all_times,
	      [&ind](int a, int b){
		double ta=getTime(a, ind);
		if (ind->err){
		  throw std::runtime_error("error");
		}
		double tb = getTime(b, ind);
		if (ind->err){
		  throw std::runtime_error("error");
		}
		if (ta == tb) return a < b;
		return ta < tb;
	      });
  } catch(...){
  }
}

//[[Rcpp::export]]
List dropUnitsRxSolve(List x){
  List ret;
  if (rxIs(x, "rxSolve")){
    ret = clone(x);
    for (int j = (int)ret.size(); j--;){
      if (rxIs(ret[j],"units")){
	RObject cur = ret[j];
	cur.attr("units") = R_NilValue;
	cur.attr("class") = R_NilValue;
      }
    }
  }
  return ret;
}

extern "C" void setSilentErr(int silent);
//[[Rcpp::export]]
bool rxSetSilentErr(int silent){
  setSilentErr(silent);
  return true;
}

#undef MV
