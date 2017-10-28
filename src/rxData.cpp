// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Check the type of an object using Rcpp
//'
//' @param obj Object to check
//' @param cls Type of class.  Only s3 classes and primitive classes are checked.
//'    For matrix types they are distinguished as \code{numeric.matrix}, \code{integer.matrix},
//'    \code{logical.matrix}, and \code{character.matrix} as well as the traditional \code{matrix}
//'    class.
//'
//' @return A boolean indicating if the object is a member of the class.
//' @keywords internal
//' @export
// [[Rcpp::export]]
bool rxIs(RObject obj, std::string cls){
  if (obj.isObject()){
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
  }
  return false;
}

// List rxSetupParameters(RObject dllInfo, RObject param, RObject theta = R_NilValue, RObject eta = R_NilValue){
//   // Purpose: Sort parameters
// }

//' Setup a data frame for solving multiple subjects at once in RxODE.
//'
//' @param df dataframe to setup; Must be in RxODE compatible format.
//' @param covNames Covariate names in dataset.
//' @param amountUnits Dosing amount units.
//' @param timeUnits Time units.
//'
//' @return A data structure to allow C-based for loop (ie solving each
//'       individual in C)
//'
//' @export
// [[Rcpp::export]]
List rxDataSetup(const DataFrame &df, const Nullable<StringVector> &covNames = R_NilValue,
		 const std::string &amountUnits = "NA", const std::string &timeUnits = "hours"){
  // Purpose: get positions of each id and the length of each id's observations
  // Separate out dose vectors and observation vectors
  IntegerVector id    = df["ID"];
  IntegerVector evid  = df["EVID"];
  NumericVector dv    = df["DV"];
  NumericVector time0 = df["TIME"];
  NumericVector amt   = df["AMT"];
  int ids = id.size();
  int lastId = id[0]-1;
  // Get the number of subjects
  // Get the number of observations
  // Get the number of doses
  int nSub = 0, nObs = 0, nDoses = 0, i = 0, j = 0, k=0;
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
  int nCovs;
  if (covNames.isNotNull()){
    nCovs=(as<StringVector>(covNames)).size();
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
        newCov[n0 + nc + n*nObsN[m-1]] = (as<NumericVector>(df[as<std::string>((as<StringVector>(covNames))[n])]))[i];
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
                          _["cov.names"]=covNames,
			  _["amount.units"]=amountUnits,
			  _["time.units"]=timeUnits
			  );
  ret.attr("class") = "RxODE.multi.data";
  return ret;
}

//[[Rcpp::export]]
List rxEventTableExpand(const int &nsub,const DataFrame &df,
                        const std::string &amountUnits = "NA", const std::string &timeUnits = "hours",
			const LogicalVector &expandData = false){
  // Purpose: Expand current event table to have number of subjects,
  // and then return rxDataSetup list.
  //IntegerVector id    = df["id"];
  IntegerVector evid  = df["evid"];
  //NumericVector dv    = df["dv"];
  NumericVector time0 = df["time"];
  NumericVector amt   = df["amt"];
  int newdim = amt.size();
  int newsub = 1;
  if (expandData[0]){
    newdim *= nsub;
    newsub = nsub;
  }
  IntegerVector nId(newdim);
  IntegerVector nEvid(newdim);
  NumericVector nDv(newdim);
  NumericVector nTime(newdim);
  NumericVector nAmt(newdim);
  for (int i = 0; i < amt.size(); i++){
    for (int j = 0; j < newsub; j++){
      nId[i+j*amt.size()] = j+1;
      nEvid[i+j*amt.size()] = evid[i];
      nTime[i+j*amt.size()] = time0[i];
      nAmt[i+j*amt.size()] = amt[i];
      nDv[i+j*amt.size()] = NA_REAL;
    }
  }
  DataFrame ndf = DataFrame::create(_["ID"]=nId,
				    _["EVID"]=nEvid,
				    _["DV"]=nDv,
				    _["TIME"]=nTime,
				    _["AMT"]=nAmt);
  List ret = rxDataSetup(ndf, R_NilValue, amountUnits, timeUnits);
  if (expandData[0]){
    return ret;
  } else {
    ret["nSub"]=nsub;
    ret["nDoses"]=nsub*as<int>(ret["nDoses"]);
    ret["nObs"]=nsub*as<int>(ret["nObs"]);
    ret.attr("class") = "RxODE.multi.data.dup";
    return ret;
  }
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
  List events;
  // 
  if (rxIs(par0, "eventTable")){
    events = as<List>(par0);
  } else if (rxIs(ev0,"eventTable")){
    events = as<List>(ev0);
  }
  
  // if (!is.null(params)){
  //   if (is.null(events) && is(params,"EventTable")){
  //     events <- params;
  //     params <- c();
  //   }
  // }
}
