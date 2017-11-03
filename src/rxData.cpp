// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

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
    for (int i = 0; i < len; i++){
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
//' @param cls Type of class.  Only s3 classes and primitive classes are checked.
//'    For matrix types they are distinguished as \code{numeric.matrix}, \code{integer.matrix},
//'    \code{logical.matrix}, and \code{character.matrix} as well as the traditional \code{matrix}
//'    class. Additionally checks for \code{event.data.frame} which is an \code{data.frame} object
//'    with \code{time},  \code{evid} and \code{amt}. (UPPER, lower or Title cases accepted)
//'
//' @return A boolean indicating if the object is a member of the class.
//' @keywords internal
//' @export
// [[Rcpp::export]]
bool rxIs(const RObject &obj, std::string cls){
  if (cls == "rx.event"){
    return (rxIs(obj, "eventTable") || rxIs(obj, "event.data.frame") || rxIs(obj, "event.matrix"));
  } else if (cls == "event.data.frame"){
    if (rxIs(obj, "data.frame")){
      CharacterVector cv =as<CharacterVector>((as<DataFrame>(obj)).names());
      return rxHasEventNames(cv);
    } else {
      return false;
    }
  } else if (cls == "event.matrix"){
    if (rxIs(obj,"numeric.matrix") && obj.hasAttribute("dimnames")){
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
  } else if (obj.isObject()){
    CharacterVector classattr = obj.attr("class");
    for (int i = 0; i < classattr.size(); i++){
      if (as<std::string>(classattr[i]) == cls){
        return true;
      }
    }
  } else {
    int type = obj.sexp_type();
    bool hasDim = obj.hasAttribute("dim");
    if (type == REALSXP){
      if (hasDim){
	if (cls == "numeric.matrix" || cls == "matrix"){
	  return true;
	} else {
	  return false;
	}
      } else {
	if (cls == "numeric")
          return true;
        else 
          return false;
      }
    }
    if (type == INTSXP){
      if (hasDim){
	if (cls == "integer.matrix" || cls == "matrix"){
          return true;
        } else {
          return false;
        }
      } else {
	if (cls == "integer")
          return true;
        else
          return false;
      }
    }
    if (type == LGLSXP){
      if (hasDim){
        if (cls == "logical.matrix" || cls == "matrix"){
          return true;
        } else {
          return false;
        }
      } else {
	if (cls == "logical")
          return true;
        else
          return false;
      }
    }
    if (type == STRSXP){
      if (hasDim){
	if (cls == "character.matrix" || cls == "matrix"){
          return true;
        } else {
          return false;
        }
      } else {
	if (cls == "character")
          return true;
        else
          return false;
      }
    }
    if (type == VECSXP){
      if (cls == "list"){
        return true;
      } else {
        return false;
      }
    }
  }
  return false;
}

// List rxSetupParameters(RObject dllInfo, RObject param, RObject theta = R_NilValue, RObject eta = R_NilValue){
//   // Purpose: Sort parameters
// }

//' Setup a data frame for solving multiple subjects at once in RxODE.
//'
//' @param ro R object to setup; Must be in RxODE compatible format.
//' @param covNames Covariate names in dataset.
//' @param amountUnits Dosing amount units.
//' @param timeUnits Time units.
//'
//' @return A data structure to allow C-based for loop (ie solving each
//'       individual in C)
//'
//' @export
// [[Rcpp::export]]
List rxDataSetup(const RObject &ro, const RObject &covNames = R_NilValue,
		 const StringVector &amountUnits = NA_STRING, const StringVector &timeUnits = "hours"){
  // Purpose: get positions of each id and the length of each id's observations
  // Separate out dose vectors and observation vectors
  if (rxIs(ro,"EventTable")){
    List et = as<List>(ro);
    Function f = et["get.EventTable"];
    DataFrame df = f();
    f = et["get.units"];
    StringVector units = f();
    StringVector amt = as<StringVector>(units["dosing"]);
    StringVector time = as<StringVector>(units["time"]);
    return rxDataSetup(df, covNames, amt, time);
  } else if (rxIs(ro,"event.data.frame")||
      rxIs(ro,"event.matrix")){
    DataFrame df = as<DataFrame>(ro);
    
    int nSub = 0, nObs = 0, nDoses = 0, i = 0, j = 0, k=0;
    IntegerVector evid  = as<IntegerVector>(df[rxcEvid]);
    bool missingId = false;
    IntegerVector id(evid.size());
    if (rxcId > -1){
      id    = as<IntegerVector>(df[rxcId]);
    } else {
      for (i = 0; i < evid.size(); i++){
	id[i]=1;
      }
      missingId=true;
    }
    bool missingDv = false;
    NumericVector dv(evid.size());
    if (rxcDv > -1){
      dv = as<NumericVector>(df[rxcDv]);
    } else {
      for (i = 0; i < evid.size(); i++){
	dv[i] = NA_REAL;
      }
      missingDv = true;
    }
    NumericVector time0 = df[rxcTime];
    NumericVector amt   = df[rxcAmt];
    int ids = id.size();
    int lastId = id[0]-1;
    // Get the number of subjects
    // Get the number of observations
    // Get the number of doses
    for (i = 0; i < ids; i++){
      if (lastId != id[i]){
        nSub++;
        lastId=id[i];
      }
      if (evid[i]){
        nDoses++;
      } else {
        nObs++;
      }
    }
    // Now create data frames of observations and events
    NumericVector newDv(nObs);
    NumericVector newTimeO(nObs);
    int nCovs = 0;
    StringVector covN;
    bool dataCov = false;
    DataFrame covDf;
    if (rxIs(covNames, "character")){
      covN = (as<StringVector>(covNames));
      nCovs= covN.size();
    } else if (rxIs(covNames, "data.frame") || rxIs(covNames,"numeric.matrix")){
      covDf = as<DataFrame>(covNames);
      covN = (as<StringVector>(covDf.names()));
      nCovs = covN.size();
      dataCov = true;
    } else {
      nCovs=0;
    }
    NumericVector newCov(nObs*nCovs);
  
    IntegerVector newEvid(nDoses);
    NumericVector newAmt(nDoses);
    NumericVector newTimeA(nDoses);

    lastId = id[0]-1;
  
    IntegerVector newId(nSub);
    IntegerVector posDose(nSub);
    IntegerVector posObs(nSub);
    IntegerVector posCov(nSub);
    IntegerVector nCov(nSub);
    IntegerVector nDose(nSub);
    IntegerVector nObsN(nSub);
    IntegerVector posEt(nSub);
    IntegerVector nEtN(nSub);
    double minTime = NA_REAL;
    double maxTime = -1e10;
    double minIdTime = NA_REAL;
    int m = 0, nEt=0;
    for (i = 0; i < ids; i++){
      if (lastId != id[i]){
        lastId     = id[i];
        newId[m]   = id[i];
        posDose[m] = j;
        posObs[m]  = k;
        posEt[m] = i;
        if (m != 0){
          nDose[m-1] = nDoses;
          nObsN[m-1]  = nObs;
          nEtN[m-1] = nEt;
        }
        nDoses = 0;
        nObs = 0;
        nEt  = 0;
        m++;
	minIdTime = time0[i];
      }
      if (minIdTime > time0[i]){
	stop("Data need to be ordered by ID and TIME.");
      } else {
	minIdTime = time0[i];
      }
      if (evid[i]){
        // Dose
        newEvid[j]  = evid[i];
        newTimeA[j] = time0[i];
        newAmt[j]   = amt[i];
        nDoses++;
        j++;
        nEt++;
      } else {
        // Observation
        newDv[k]    = dv[i];
        newTimeO[k] = time0[i];
        if (nObs == 0){
          minTime = time0[i];
        } else if (time0[i] > maxTime) {
          maxTime = time0[i];
        }
        nObs++;
        k++;
        nEt++;
      }
    }
    nDose[m-1]=nDoses;
    nObsN[m-1]=nObs;
    nEtN[m-1] = nEt;
    if (dataCov && (as<NumericVector>(covDf[0])).size() != nObs){
      stop("Covariate data needs ");
    }
    // Covariates are stacked by id that is
    // id=cov1,cov1,cov1,cov2,cov2,cov2,...
    lastId = id[0]-1;
    int n0 = 0, n = 0, nc = 0;
    m = 0;
    for (i = 0; i < ids; i++){
      if (lastId != id[i]){
        lastId     = id[i];
        if (m != 0){
          n0 += nObsN[m-1]*nCovs;
        }
        posCov[m] = n0;
        nCov[m]   = nObsN[m]*nCovs;
        nc = 0;
        m++;
      }
      if (!evid[i]){
        // Observation
        for (n = 0; n < nCovs; n++){
	  if (dataCov){
	    newCov[n0 + nc + n*nObsN[m-1]] = (as<NumericVector>(covDf[as<std::string>(covN[n])]))[nc];
	  } else {
	    newCov[n0 + nc + n*nObsN[m-1]] = (as<NumericVector>(df[as<std::string>(covN[n])]))[i];
          }
        }
        nc++;
      }
    }
    // nCov[m-1] = nObs*nCovs;
    List ret = List::create(_["dose"] = DataFrame::create(_["evid"]    = newEvid,
                                                          _["time"]    = newTimeA,
                                                          _["amt"]     = newAmt),
                            _["obs"]  = DataFrame::create(_["dv"]      = newDv,
                                                          _["time"]    = newTimeO),
                            _["ids"]  = DataFrame::create(_["id"]      = newId,
                                                          _["posDose"] = posDose,
                                                          _["posObs"]  = posObs,
                                                          _["posCov"]  = posCov,
                                                          _["posEvent"]= posEt,
                                                          _["nDose"]   = nDose,
                                                          _["nObs"]    = nObsN,
                                                          _["nCov"]    = nCov,
                                                          _["nEvent"]   = nEtN),
                            _["et"] = DataFrame::create(_["evid"]=evid,
                                                        _["time"]=time0),
                            _["cov"]=newCov,
                            _["nSub"]=nSub,
                            _["nDoses"]=newEvid.size(),
                            _["nObs"]=newDv.size(),
                            _["min.time"] = minTime,
                            _["max.time"] = maxTime,
                            _["cov.names"]=(covN.size() == 0 ? R_NilValue : wrap(covN)),
                            _["amount.units"]=amountUnits,
                            _["time.units"]=timeUnits,
			    _["missing.id"]=missingId,
			    _["missing.dv"]=missingDv
                            );
    ret.attr("class") = "RxODE.multi.data";
    return ret;
  } else {
    stop("Data is not setup appropriately.");
  }
}
//' All model variables for a RxODE object
//'
//' Return all the known model variables for a specified RxODE object
//'
//' These items are only calculated after compilation; they are
//' built-into the RxODE compiled DLL.
//'
//' @param obj RxODE family of objects
//'
//' @return A list of RxODE model properties including:
//'
//' \item{params}{ a character vector of names of the model parameters}
//' \item{lhs}{ a character vector of the names of the model calculated parameters}
//' \item{state}{ a character vector of the compartments in RxODE object}
//' \item{trans}{ a named vector of translated model properties
//'       including what type of jacobian is specified, the \code{C} function prefixes,
//'       as well as the \code{C} functions names to be called through the compiled model.}
//' \item{md5}{a named vector that gives the digest of the model (\code{file_md5}) and the parsed model
//'      (\code{parsed_md5})}
//' \item{model}{ a named vector giving the input model (\code{model}),
//'    normalized model (no comments and standard syntax for parsing, \code{normModel}),
//'    and interim code that is used to generate the final C file \code{parseModel}}
//'
//' @keywords internal
//' @author Matthew L.Fidler
//' @export
// [[Rcpp::export]]
List rxModelVars(RObject obj = R_NilValue){
  if (rxIs(obj,"RxODE")) {
    Function f = as<Function>((as<List>((as<List>(obj))["cmpMgr"]))["rxDll"]);
    List lst = f();
    return lst["modVars"];
  } else if (rxIs(obj,"solveRxODE")){
    Environment e = as<Environment>((as<Environment>(obj.attr(".env")))["env"]);
    return  rxModelVars(e["out"]);
  } else if (rxIs(obj,"rxDll")){
    List lobj = (as<List>(obj))["modVars"];
    return lobj;
  } else if (rxIs(obj, "RxCompilationManager")){
    Function f =  as<Function>((as<List>(obj))["rxDll"]);
    List lst = f();
    return (lst["modVars"]);
  } else if (rxIs(obj, "character")){
    Environment RxODE("package:RxODE");
    Function f = RxODE["rxModelVars.character"];
    return f(obj);
  } else if (rxIs(obj,"list")){
    bool params=false, lhs=false, state=false, trans=false, ini=false, model=false, md5=false, podo=false, dfdy=false;
    List lobj  = as<List>(obj);
    CharacterVector nobj = lobj.names();
    for (int i = 0; i < nobj.size(); i++){
      if (!params && nobj[i]== "params"){
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
  } else {
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
RObject rxState(RObject obj = R_NilValue, RObject state = R_NilValue){
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

//' Parameters specified by the model
//'
//' This return the model's parameters that are required to solve the
//' ODE system.
//'
//' @inheritParams rxModelVars
//'
//' @return a character vector listing the parameters in the model.
//'
//' @author Matthew L.Fidler
//' @export
//[[Rcpp::export]]
CharacterVector rxParams(RObject obj = R_NilValue){
  List modVar = rxModelVars(obj);
  CharacterVector ret = modVar["params"];
  return ret;
}


//' Jacobain and parameter derivatives
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
CharacterVector rxDfdy(RObject obj = R_NilValue){
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
CharacterVector rxLhs(RObject obj = R_NilValue){
  List modVar = rxModelVars(obj);
  CharacterVector ret = modVar["lhs"];
  return ret;
}

RObject rxSolveCpp(List args, Environment e){
  List dll = as<List>(e["dll"]);
  List modVars = as<List>(dll["modVars"]);
  CharacterVector trans = modVars["trans"];
  CharacterVector state = modVars["state"];
  CharacterVector lhs = modVars["lhs"];
  CharacterVector pars = modVars["pars"];
  Nullable<IntegerVector> stateIgnore = modVars["state.ignore"];
  NumericVector params;
  RObject par0 = args["params"];
  RObject ev0  = args["events"];
  RObject covs = args["covs"];
  List events;
  // Setup events
  if (rxIs(par0, "rx.event")){
    events = rxDataSetup(par0, covs);
  } else if (rxIs(ev0,"rx.event")){
    events = rxDataSetup(ev0, covs);
  }
}
