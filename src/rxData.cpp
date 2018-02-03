// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <thread>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
extern "C" {
#include "solve.h"
}

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
//' @author Matthew L. Fidler
//' @export
// [[Rcpp::export]]
bool rxIs(const RObject &obj, std::string cls){
  if (cls == "rx.event"){
    return (rxIs(obj, "EventTable") || rxIs(obj, "event.data.frame") || rxIs(obj, "event.matrix"));
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
	return false; // nocov
      }
    } else {
      return false;
    }
  } else if (obj.isObject()){
    CharacterVector classattr = obj.attr("class");
    for (int i = 0; i < classattr.size(); i++){
      if (as<std::string>(classattr[i]) == cls){
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
    if (type == ENVSXP){
      if (cls == "environment"){
	return true;
      } else {
	return false;
      }
    }
    if (type == EXTPTRSXP){
      if (cls == "externalptr" || cls == "refObject"){
	return true;
      } else {
	return false;
      }
    }
  }
  return false;
}

extern "C" int rxIsC(SEXP obj, const char *cls){
  std::string str(cls);
  if (rxIs(as<RObject>(obj),cls)){
    return 1;
  } else {
    return 0;
  }
}

RObject rxSimSigma(const RObject &sigma,
		   const RObject &df,
		   int ncores,
		   const bool &isChol,
		   int nObs){
  if (rxIs(sigma, "numeric.matrix")){
    // FIXME more distributions
    NumericMatrix sigmaM(sigma);
    if (sigmaM.nrow() != sigmaM.ncol()){
      stop("The matrix must be a square matrix.");
    }
    if (!sigmaM.hasAttribute("dimnames")){
      stop("The matrix must have named dimensions.");
    }
    List dimnames = sigmaM.attr("dimnames");
    StringVector simNames = as<StringVector>(dimnames[1]);
    Environment base("package:base");
    Function loadNamespace=base["loadNamespace"];
    Environment mvnfast = loadNamespace("mvnfast");
    NumericMatrix simMat(nObs,sigmaM.ncol());
    NumericVector m(sigmaM.ncol());
    // I'm unsure if this for loop is necessary.
    // for (int i = 0; i < m.size(); i++){
    //   m[i] = 0;
    // }
    // Ncores = 1?  Should it be parallelized when it can be...?
    // Note that if so, the number of cores also affects the output.
    if (df.isNULL()){
      Function rmvn = as<Function>(mvnfast["rmvn"]);
      rmvn(_["n"]=nObs, _["mu"]=m, _["sigma"]=sigmaM, _["ncores"]=ncores, _["isChol"]=isChol, _["A"] = simMat); // simMat is updated with the random deviates
    } else {
      double df2 = as<double>(df);
      if (R_FINITE(df2)){
	Function rmvt = as<Function>(mvnfast["rmvt"]);
        rmvt(_["n"]=nObs, _["mu"]=m, _["sigma"]=sigmaM, _["df"] = df, _["ncores"]=ncores, _["isChol"]=isChol, _["A"] = simMat);
      } else {
	Function rmvn = as<Function>(mvnfast["rmvn"]);
        rmvn(_["n"]=nObs, _["mu"]=m, _["sigma"]=sigmaM, _["ncores"]=ncores, _["isChol"]=isChol, _["A"] = simMat);
      }
    }
    simMat.attr("dimnames") = List::create(R_NilValue, simNames);
    return wrap(simMat);
  } else {
    return R_NilValue;
  }
}


bool foundEnv = false;
Environment _rxModels;
void getRxModels(){
  if (!foundEnv){ // minimize R call
    Environment RxODE("package:RxODE");
    Function f = as<Function>(RxODE["rxModels_"]);
    _rxModels = f();
    foundEnv = true;
  }
}

extern "C" SEXP rxSimSigmaC(rx_solving_options *op,
			    int nObs){
  bool isChol = false;
  if (op->isChol == 1){
    isChol=true;
  }
  RObject ret;
  int n = op->sigmaSize;
  if (n > 0){
    NumericMatrix sigma = NumericMatrix(n, n);
    n  = n*n;
    for (int i = 0; i < n; i++){
      sigma[i] = op->sigma[i];
    }
    getRxModels();
    sigma.attr("dimnames") = List::create(_rxModels[".simNames"], _rxModels[".simNames"]);
    Nullable<NumericVector> dfN(1);
    if (op->df < 0){
      dfN = R_NilValue;
    } else {
      NumericVector tmp = NumericVector(1);
      tmp[0] = op->df;
    }
    ret = rxSimSigma(as<RObject>(sigma), as<RObject>(dfN), op->ncoresRV, isChol, nObs);
  } else {
    ret = R_NilValue;
  }
  return wrap(ret);
}
// [[Rcpp::export]]
List rxDataSetup(const RObject &ro,
		 const RObject &covNames = R_NilValue,
		 const RObject &sigma = R_NilValue,
		 const RObject &df = R_NilValue,
		 const int &ncoresRV = 1,
		 const bool &isChol = false,
                 const int &nDisplayProgress = 10000,
		 const StringVector &amountUnits = NA_STRING,
		 const StringVector &timeUnits = "hours"){
  // Purpose: get positions of each id and the length of each id's observations
  // Separate out dose vectors and observation vectors
  if (rxIs(ro,"EventTable")){
    List et = List(ro);
    Function f = et["get.EventTable"];
    DataFrame dataf = f();
    f = et["get.units"];
    RObject unitsRO = f();
    CharacterVector units;
    int i, n;
    if (rxIs(unitsRO, "character")){
      units = as<CharacterVector>(unitsRO);
      n=units.size();
      for (i =0; i<n; i++){
	if (units[i] == "NA"){
	  units[i] = NA_STRING;
	}
      }
    } else {
      units = CharacterVector::create(_["dosing"]=NA_STRING,
				      _["time"]=NA_STRING);
    }
    // {
    //   units = StringVector(unitsRO);
    //   if (units[0] == "NA"){
    // 	units[0] = NA_STRING;
    //   }
    //   if (units[1] == "NA"){
    // 	units[1] = NA_STRING;
    //   }
    // } else {
    //   // Otherwise this is likely 2 NAs.
    //   units[0] = NA_STRING;
    //   units[1] = NA_STRING;
    //   StringVector units2(2);
    //   units2[0] = "dosing";
    //   units2[1] = "time";
    //   units.names() = units;
    // }
    CharacterVector amt = (units["dosing"] == NA_STRING) ? StringVector::create(NA_STRING) : as<StringVector>(units["dosing"]);
    CharacterVector time = (units["time"] == NA_STRING) ? StringVector::create(NA_STRING) : as<StringVector>(units["time"]);
    return rxDataSetup(dataf, covNames, sigma, df, ncoresRV, isChol, nDisplayProgress, amt, time);
  } else if (rxIs(ro,"event.data.frame")||
      rxIs(ro,"event.matrix")){
    DataFrame dataf = as<DataFrame>(ro);
    int nSub = 0, nObs = 0, nDoses = 0, i = 0, j = 0, k=0;
    // Since the event data frame can be "wild", these need to be
    // converted to integers.
    IntegerVector evid  = as<IntegerVector>(dataf[rxcEvid]);
    bool missingId = false;
    IntegerVector id(evid.size());
    if (rxcId > -1){
      id    = as<IntegerVector>(dataf[rxcId]);
    } else {
      for (i = 0; i < evid.size(); i++){
	id[i]=1;
      }
      missingId=true;
    }
    bool missingDv = false;
    NumericVector dv(evid.size());
    if (rxcDv > -1){
      dv = as<NumericVector>(dataf[rxcDv]);
    } else {
      for (i = 0; i < evid.size(); i++){
	dv[i] = NA_REAL;
      }
      missingDv = true;
    }
    NumericVector time0 = dataf[rxcTime];
    NumericVector amt   = dataf[rxcAmt];
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
    StringVector covN, simN;
    bool dataCov = false;
    DataFrame covDf;
    NumericMatrix simMat;
    StringVector simNames;
    bool simVals = false;
    RObject tmp_ro = rxSimSigma(sigma, df, ncoresRV, isChol, nObs);
    if (!tmp_ro.isNULL()){
      simMat = NumericMatrix(tmp_ro);
      List dimnames = simMat.attr("dimnames");
      simNames = StringVector(dimnames[1]);
      simVals = true;
    }
    int nCovObs = 0;
    if (rxIs(covNames, "character")){
      covN = StringVector(covNames);
      nCovObs = covN.size();
      if (simVals){
	nCovs= nCovObs + simNames.size();
      } else {
        nCovs = nCovObs;
      }
    } else if (rxIs(covNames, "data.frame") || rxIs(covNames,"numeric.matrix")){
      covDf = as<DataFrame>(covNames);
      covN = StringVector(covDf.names());
      nCovObs = covN.size();
      if (simVals){
	nCovs = nCovObs + simNames.size();
      } else {
	nCovs = nCovObs;
      }
      dataCov = true;
    } else if (simVals) {
      nCovs= simNames.size();
    }
    
    // Rprintf("nObs: %d; nCovs: %d\n", nObs, nCovs);
    NumericVector newCov(nObs*nCovs);
    
    IntegerVector newEvid(nDoses);
    IntegerVector idose(nDoses);
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
    NumericVector Hmax(nSub); // For LSODA default
    IntegerVector rc(nSub);
    
    double minTime = NA_REAL;
    double maxTime = -1e10;
    double minIdTime = NA_REAL;
    double lastTime = 0;
    //hmax <- max(abs(diff(event.table$time)))
    double mdiff = 0;
    double HmaxA = 0;
    double tmp;
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
	  Hmax[m-1] = mdiff;
	  rc[m-1] = 0;
        }
        nDoses = 0;
        nObs = 0;
        nEt  = 0;
        m++;
	minIdTime = time0[i];
	lastTime = time0[i];
	mdiff = 0;
      }
      if (minIdTime > time0[i]){
	stop("Data need to be ordered by ID and TIME.");
      } else {
	minIdTime = time0[i];
      }
      tmp = time0[i]-lastTime;
      if (tmp > mdiff){
	mdiff = tmp;
	if (tmp > HmaxA){
          HmaxA = tmp;
	}
      }
      lastTime = time0[i];
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
    Hmax[m-1] = mdiff;
    rc[m-1] = 0;
    k = 0;
    if (dataCov && covDf.nrow() != nObs){
      if (covDf.nrow() == ids){
	List covDf2(covDf.nrow());
	for (j = 0; j < covDf.ncol(); j++){
          covDf2[j] = NumericVector(nObs);
	}
	for (i = 0; i < ids; i++){
	  if (!evid[i]){
	    for (j = 0; j < covDf.ncol(); j++){
              NumericVector cur = covDf2[j];
	      cur[k] = (as<NumericVector>(covDf[j]))[i];
            }
	    k++;
	  }
	}
        covDf2.attr("names") = covDf.attr("names");
        covDf2.attr("class") = "data.frame";
        covDf2.attr("row.names") = IntegerVector::create(NA_INTEGER,-nObs);
	covDf = as<DataFrame>(covDf2);
      } else {
	stop("Covariate data needs to match the number of observations in the overall dataset.");
      }
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
	  k = n0 + nc + n*nObsN[m-1];
	  if (n < nCovObs){
	    if (dataCov){
              newCov[k] = (as<NumericVector>(covDf[as<std::string>(covN[n])]))[nc];
            } else {
              newCov[k] = (as<NumericVector>(dataf[as<std::string>(covN[n])]))[i];
            }
          } else {
	    newCov[k] = simMat(nc, n-nCovObs);
          }
        }
        nc++;
      }
    }
    // nCov[m-1] = nObs*nCovs;
    List ret = List::create(_["dose"] = DataFrame::create(_["evid"]         = newEvid,
                                                          _["time"]         = newTimeA,
                                                          _["amt"]          = newAmt),
                            _["obs"]  = DataFrame::create(_["dv"]           = newDv,
                                                          _["time"]         = newTimeO),
                            _["ids"]  = DataFrame::create(_["id"]           = newId,
                                                          _["posDose"]      = posDose,
                                                          _["posObs"]       = posObs,
                                                          _["posCov"]       = posCov,
                                                          _["posEvent"]     = posEt,
                                                          _["nDose"]        = nDose,
                                                          _["nObs"]         = nObsN,
                                                          _["nCov"]         = nCov,
                                                          _["nEvent"]       = nEtN,
							  _["HmaxDefault"]  = Hmax,
							  _["rc"]           = rc),
                            _["et"] = DataFrame::create(_["evid"]=evid,
                                                        _["time"]=time0),
                            _["cov"]=newCov,
                            _["nSub"]=nSub,
                            _["nDoses"]=newEvid.size(),
                            _["nObs"]=newDv.size(),
                            _["min.time"] = minTime,
                            _["max.time"] = maxTime,
                            _["cov.names"]=(covN.size() == 0 ? R_NilValue : wrap(covN)),
			    _["n.observed.covariates"] = nCovObs,
			    _["simulated.vars"] = (simNames.size()== 0 ? R_NilValue : wrap(simNames)),
			    _["sigma"] = wrap(sigma),
                            _["amount.units"]=(as<std::string>(amountUnits) == "NA") ? StringVector::create(NA_STRING) : amountUnits,
                            _["time.units"]=(as<std::string>(timeUnits)  == "NA") ? StringVector::create(NA_STRING) : timeUnits,
			    _["missing.id"]=missingId,
			    _["missing.dv"]=missingDv,
			    _["ncoresRV"] = wrap(ncoresRV),
			    _["isChol"] = wrap(isChol)
                            );
    // Not sure why, but putting this in above gives errors...
    ret["df"]= df;
    ret["idose"] = idose;
    ret["Hmax"] = HmaxA;
    ret["nDisplayProgress"] = nDisplayProgress;
    ret.attr("class") = "RxODE.multi.data";
    return ret;
  } else if (rxIs(ro,"list")){
    Function asDf("as.data.frame", R_BaseNamespace);
    return rxDataSetup(asDf(ro), covNames,sigma, df, ncoresRV, isChol, nDisplayProgress, amountUnits, timeUnits);
  } else {
    stop("Data is not setup appropriately.");
  }
}

extern "C" void set_solve(rx_solve *rx);

rx_solve *getRxSolve(SEXP ptr){
  if (rxIs(ptr,"RxODE.pointer.multi")){
    List lst = List(ptr);
    rxUpdateFuns(lst["trans"]);
    rx_solve *ret = getRxSolve_();
    // Also assign it.
    set_solve(ret);
    return ret;
  } else {
    stop("Cannot get the solving data (getRxSolve).");
  }
  rx_solve *o;
  return o;
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
List rxModelVars(const RObject &obj){
  if (rxIs(obj, "rxModelVars")){
    List ret(obj);
    return ret;
  } else if (rxIs(obj,"RxODE")) {
    Environment e = as<Environment>(obj);
    List rxDll = e["rxDll"];
    return rxDll["modVars"];
  } else if (rxIs(obj,"rxSolve")){
    CharacterVector cls = obj.attr("class");
    Environment e = as<Environment>(cls.attr(".RxODE.env"));
    return  rxModelVars(as<RObject>(e["args.object"]));
  } else if (rxIs(obj,"rxDll")){
    List lobj = (as<List>(obj))["modVars"];
    return lobj;
  } else if (rxIs(obj, "character")){
    getRxModels();
    std::string sobj =as<std::string>(obj);
    if (_rxModels.exists(sobj)){
      RObject obj1 = _rxModels.get(sobj);
      if (rxIs(obj1, "rxModelVars")){
	return as<List>(obj1);
      } else if (rxIs(obj1, "RxODE")){
	return rxModelVars(obj1);
      }
    }
    Environment RxODE("package:RxODE");
    Function f = as<Function>(RxODE["rxModelVars.character"]);
    return f(obj);
  } else if (rxIs(obj,"list")){
    bool params=false, lhs=false, state=false, trans=false, ini=false, model=false, md5=false, podo=false, dfdy=false;
    List lobj  = as<List>(obj);
    CharacterVector nobj = lobj.names();
    for (int i = 0; i < nobj.size(); i++){
      if (nobj[i] == "modVars"){
	return(rxModelVars(lobj["modVars"]));
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
CharacterVector rxParams(const RObject &obj){
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
NumericVector rxInits(const RObject &obj,
		      RObject vec = R_NilValue,
		      Nullable<CharacterVector> req = R_NilValue,
		      double defaultValue = 0,
		      bool noerror = false,
		      bool noini=false){
  if (vec.isNULL()){
    return rxInits0(obj, R_NilValue, req, defaultValue, noerror,noini);
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
    return rxInits0(obj, vec2, req, defaultValue, noerror,noini);
  } else if (rxIs(vec, "integer") || rxIs(vec, "numeric")){
    return rxInits0(obj, as<NumericVector>(vec), req, defaultValue, noerror,noini);
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

NumericMatrix rxSetupParamsThetaEta(const RObject &params = R_NilValue,
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
      // Create the matrix
      NumericVector thetaN;
      if (rxIs(theta,"numeric") || rxIs(theta,"integer")){
        thetaN = as<NumericVector>(theta);
      } else if (rxIs(theta, "matrix")){
        NumericMatrix thetaM = as<NumericMatrix>(theta);
        if (thetaM.nrow() == 1){
          thetaN = NumericVector(thetaM.ncol());
          for (i = 0 ; i < thetaM.ncol(); i++){
            thetaN[i] = thetaM(1,i);
          }
        } else if (thetaM.ncol() == 1){
          thetaN = NumericVector(thetaM.nrow());
          for (i = 0 ; i < thetaM.ncol(); i++){
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
          for (i = 0 ; i < etaM.ncol(); i++){
            etaN[i] = etaM(0, i);
          }
        } else if (etaM.ncol() == 1){
          etaN = NumericVector(etaM.nrow());
          for (i = 0 ; i < etaM.ncol(); i++){
            etaN[i] = etaM(i, 0);
          }
        } else {
          stop("'eta' is not compatible with params, check dimensions to make sure they are compatible.");
        }
      } else if (!eta.isNULL()){
        stop("'eta' is not compatible with params, check dimensions to make sure they are compatible.");
      }
      NumericMatrix tmp1(1, thetaN.size()+etaN.size());
      CharacterVector tmpN = CharacterVector(tmp1.size());
      for (i = 0; i < thetaN.size(); i++){
        tmp1(0, i) = thetaN[i];
        tmpN[i] = "THETA[" + std::to_string(i + 1) + "]";
      }
      for (; i < thetaN.size()+ etaN.size(); i++){
        tmp1(0, i) = etaN[i - thetaN.size()];
        tmpN[i] = "ETA[" + std::to_string(i - thetaN.size() + 1) + "]";
      }
      tmp1.attr("dimnames") = List::create(R_NilValue, tmpN);
      parMat = tmp1;
    }
  } else if (rxIs(params, "data.frame") || rxIs(params, "matrix")){
    if (rxIs(params,"data.frame")){
      DataFrame tmp = as<DataFrame>(params);
      int nr = tmp.nrows();
      NumericMatrix tmpM(nr,tmp.size());
      for (i = 0; i < tmp.size(); i++){
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
        for (i = 0 ; i < thetaM.ncol(); i++){
          thetaN[i] = thetaM(1,i);
        }
      } else if (thetaM.ncol() == 1){
        thetaN = NumericVector(thetaM.nrow());
        for (i = 0 ; i < thetaM.ncol(); i++){
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
        for (i = 0 ; i < etaM.ncol(); i++){
          etaN[i] = etaM(0, i);
        }
      } else if (etaM.ncol() == 1){
        etaN = NumericVector(etaM.nrow());
        for (i = 0 ; i < etaM.ncol(); i++){
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
  return parMat;
}

//' Setup Data and Parameters
//'
//' @inheritParams rxSolve
//' @param sigma Named sigma matrix.
//' @param sigmaDf The degrees of freedom of a t-distribution for
//'     simulation.  By default this is \code{NULL} which is
//'     equivalent to \code{Inf} degrees, or to simulate from a normal
//'     distribution instead of a t-distribution.
//' @param nCoresRV Number of cores for residual simulation.  This,
//'     along with the seed, affects both the outcome and speed of
//'     simulation. By default it is one.
//' @param sigmaIsChol Indicates if the \code{sigma} supplied is a
//'     Cholesky decomposed matrix instead of the traditional
//'     symmetric matrix.
//' @return Data setup for running C-based RxODE runs.
//' @author Matthew L. Fidler
//' @keywords internal
//' @export
//[[Rcpp::export]]
List rxDataParSetup(const RObject &object,
		    const RObject &params = R_NilValue,
		    const RObject &events = R_NilValue,
		    const RObject &inits = R_NilValue,
		    const RObject &covs  = R_NilValue,
		    const RObject &sigma= R_NilValue,
		    const RObject &sigmaDf= R_NilValue,
		    const int &nCoresRV= 1,
		    const bool &sigmaIsChol= false,
                    const int &nDisplayProgress = 10000,
		    const StringVector &amountUnits = NA_STRING,
		    const StringVector &timeUnits = "hours",
		    const RObject &theta = R_NilValue,
                    const RObject &eta = R_NilValue,
		    const RObject &scale = R_NilValue,
		    const Nullable<List> &extraArgs = R_NilValue){
  List modVars = rxModelVars(object);
  CharacterVector state = modVars["state"];
  // The initial conditions cannot be changed for each individual; If
  // they do they need to be a parameter.
  NumericVector initsC = rxInits(object, inits, state, 0.0);
  NumericVector scaleC = rxSetupScale(object, scale, extraArgs);
  // The parameter vector/matrix/data frame contains the parameters
  // that will be used.
  RObject par0 = params;
  RObject ev0  = events;
  RObject ev1;
  RObject par1;
  if (rxIs(par0, "rx.event")){
    // Swapped events and parameters
    ev1 = par0;
    par1 = ev0;
  } else if (rxIs(ev0, "rx.event")) {
    ev1 = ev0;
    par1 = par0;
  } else {
    stop("Need some event information (observation/dosing times) to solve.\nYou can use either 'eventTable' or an RxODE compatible data frame/matrix.");
  }
  // Now get the parameters (and covariates)
  //
  // Unspecified parameters can be found in the modVars["ini"]
  NumericVector modVarsIni = modVars["ini"];
  // The event table can contain covariate information, if it is acutally a data frame or matrix.
  Nullable<CharacterVector> covnames0, simnames0;
  CharacterVector covnames, simnames;
  CharacterVector pars = modVars["params"];
  int i, j, k = 0;
  CharacterVector tmpCv;
  List ret;
  if (!rxIs(ev1,"EventTable") &&  covs.isNULL()){
    // Now covnames is setup correctly, import into a setup data table.
    // In this case the events are a data frame or matrix
    CharacterVector tmpCv =as<CharacterVector>((as<DataFrame>(ev1)).names());
    for (i = 0; i < pars.size(); i++){
      for (j = 0; j < tmpCv.size(); j++){
	if (pars[i] == tmpCv[j]){
	  k++;
	  break;
	}
      }
    }
    covnames = CharacterVector(k);
    k = 0;
    for (i = 0; i < pars.size(); i++){
      for (j = 0; j < tmpCv.size(); j++){
	if (pars[i] == tmpCv[j]){
	  covnames[k] = pars[i];
	  k++;
	  break;
	}
      }
    }
    
    ret = rxDataSetup(ev1, (covnames.size() == 0 ? R_NilValue : wrap(covnames)),
		      sigma, sigmaDf, nCoresRV, sigmaIsChol, nDisplayProgress, 
		      amountUnits, timeUnits);
  } else {
    ret = rxDataSetup(ev1, covs, sigma, sigmaDf,
		      nCoresRV, sigmaIsChol, nDisplayProgress, amountUnits, timeUnits);
    covnames0 = as<Nullable<CharacterVector>>(ret["cov.names"]);
    if (!covnames0.isNull()){
      covnames = CharacterVector(covnames0);
    }
  }
  simnames0 = as<Nullable<CharacterVector>>(ret["simulated.vars"]);
  if (!simnames0.isNull()){
    simnames = CharacterVector(simnames0);
  }
  NumericMatrix parMat = rxSetupParamsThetaEta(par1, theta, eta);
  int nSub = as<int>(ret["nSub"]);
  if (parMat.nrow() % nSub != 0){
    stop("The Number of parameters must be a multiple of the number of subjects.");
  }
  k = 0;
  IntegerVector pcov(covnames.size());
  IntegerVector svar(simnames.size());
  for (i = 0; i < covnames.size(); i++){
    for (j = 0; j < pars.size(); j++){
      if (covnames[i] == pars[j]){
	pcov[i] = j + 1;
	break;
      }
    }
  }
  for (i = 0; i < simnames.size(); i++){
    for (j = 0; j < pars.size(); j++){
      if (simnames[i] == pars[j]){
        svar[i] = j;
        break;
      }
    }
  }
  // Now pcov gives the which for the covariate parameters.
  // Now check if we have all the parameters that are needed.
  std::string errStr = "";
  bool allPars = true;
  bool curPar = false;
  IntegerVector posPar(pars.size());
  CharacterVector nms = modVarsIni.names();
  Nullable<CharacterVector> nmP2 = (as<List>(parMat.attr("dimnames")))[1];
  CharacterVector nmP;
  if (!nmP2.isNull()){
    nmP = CharacterVector(nmP2);
  }
  for (i = 0; i < pars.size(); i++){
    curPar = false;
    // integers are faster to compare than strings.
    for (j = 0; j < pcov.size(); j++){
      if (pcov[j] == i + 1){
	posPar[i] = 0; // Covariates are zeroed out.
	curPar = true;
	break;
      }
    }
    // First Check to see if the user specified the parameter.
    if (!curPar){
      for (j = 0; j < nmP.size(); j++){
        if (nmP[j] == pars[i]){
          curPar = true;
          posPar[i] = j + 1;
          break;
        }
      }
    }
    // Now check $ini
    if (!curPar){
      for(j = 0; j < modVarsIni.size(); j++){
        if (nms[j] == pars[i]){
          curPar = true;
          posPar[i] = -j - 1;
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
    CharacterVector modSyntax = modVars["model"];
    Rcout << "Model:\n\n" + modSyntax[0] + "\n";
    stop(errStr);
  }
  // Now  the parameter names are setup.
  // The parameters are setup in a numeric vector in order of pars
  int nr = parMat.nrow();
  if (nr == 0) nr = 1;
  NumericVector parsVec(pars.size()*nr);
  j = 0;
  for (i = 0; i < parsVec.size(); i++){
    j = floor(i / pars.size());
    k = i % pars.size();
    if (posPar[k] == 0){
      parsVec[i] = 0;
    } else if (posPar[k] > 0){
      // posPar[i] = j + 1;
      parsVec[i] = parMat(j, posPar[k]-1);
    } else {
      // posPar[i] = -j - 1;
      parsVec[i] = modVarsIni[-(posPar[k]+1)];
    }
  }
  ret["pars"] = parsVec;
  nr = parMat.nrow() / nSub;
  if (nr == 0) nr = 1;
  ret["nsim"] = nr;
  // NumericVector initsS = NumericVector(initsC.size()*nSub*nr);
  // for (i = 0; i < initsS.size(); i++){
  //   j = i % initsS.size();
  //   initsS[i] =  initsC[j];
  // }
  ret["inits"] = initsC;
  ret["scale"] = scaleC;
  // ret["inits.full"] = initsS;
  ret["n.pars"] = (int)(pars.size());
  ret["pcov"] = pcov;
  ret["svar"] = svar;
  ret["neq"] = state.size();
  DataFrame et      = as<DataFrame>(ret["et"]);
  NumericVector solve(state.size()*et.nrow()*nr);
  ret["solve"] = solve;
  CharacterVector lhs = as<CharacterVector>(modVars["lhs"]);
  ret["lhs"] = NumericVector(lhs.size()*nSub*nr);
  ret["lhsSize"] = lhs.size();
  ret["InfusionRate"] =NumericVector(state.size()*nSub*nr);
  ret["BadDose"] =IntegerVector(state.size()*nSub*nr);
  ret["state.ignore"] = modVars["state.ignore"];
  ret["trans"] = modVars["trans"];
  NumericVector atol(state.size(), 1e-08);
  NumericVector rtol(state.size(), 1e-06);
  ret["atol"] = atol;
  ret["rtol"] = rtol;
  CharacterVector cls(2);
  cls(1) = "RxODE.par.data";
  cls(0) = "RxODE.multi.data";
  ret.attr("class") = cls;
  return ret;
}

SEXP rxSolvingOptions(const RObject &object,
                      const std::string &method = "liblsoda",
                      const Nullable<LogicalVector> &transit_abs = R_NilValue,
                      const double atol = 1.0e-8,
                      const double rtol = 1.0e-6,
                      const int maxsteps = 5000,
                      const double hmin = 0,
		      const double hini = 0,
		      const int maxordn = 12,
                      const int maxords = 5,
		      const int cores = 1,
		      const int ncov = 0,
		      int *par_cov = NULL,
		      int do_par_cov = 0,
		      double *inits = NULL,
		      double *scale = NULL,
		      std::string covs_interpolation = "linear",
		      double hmax2 = 0,
                      double *atol2 = NULL,
                      double *rtol2 = NULL,
                      int nDisplayProgress = 10000,
		      RObject sigma = R_NilValue,
                      RObject df = R_NilValue,
                      int ncoresRV = 1,
                      int isChol = 1,
                      int *svar =NULL){
  if (maxordn < 1 || maxordn > 12){
    stop("'maxordn' must be >1 and <= 12.");
  }
  if (maxords < 1 || maxords > 5){
    stop("'maxords' must be >1 and <= 5.");
  }
  if (hmin < 0){
    stop("'hmin' must be a non-negative value.");
  }
  // HMAX is determined by the problem since it can be thought of as the maximum difference of the event table's time
  if (hini < 0){
    stop("'hini' must be a non-negative value.");
  }
  List modVars = rxModelVars(object);
  int transit = 0;
  if (transit_abs.isNull()){
    transit = modVars["podo"];
    if (transit){
      warning("Assumed transit compartment model since 'podo' is in the model.");
    }
  } else {
    LogicalVector tr = LogicalVector(transit_abs);
    if (tr[0]){
      transit=  1;
    }
  }
  int is_locf = 0;
  double f1 = 1.0, f2 = 0.0;
  int kind=1;
  if (covs_interpolation == "linear"){
  } else if (covs_interpolation == "constant" || covs_interpolation == "locf" || covs_interpolation == "LOCF"){
    f2 = 0.0;
    f1 = 1.0;
    kind = 0;
    is_locf=1;
  } else if (covs_interpolation == "nocb" || covs_interpolation == "NOCB"){
    f2 = 1.0;
    f1 = 0.0;
    kind = 0;
    is_locf=2;
  }  else if (covs_interpolation == "midpoint"){
    f1 = f2 = 0.5;
    kind = 0;
    is_locf=3;
  } else {
    stop("Unknown covariate interpolation specified.");
  }
  int st=0;
  if (method == "liblsoda"){
    st = 2;
  } else if (method == "lsoda"){
    st = 1;
  } else if (method == "dop853"){
    st = 0;
  } else {
    stop("Unknown ODE solving method specified.");
  }
  CharacterVector lhs = as<CharacterVector>(modVars["lhs"]);
  CharacterVector state = as<CharacterVector>(modVars["state"]);
  CharacterVector params = as<CharacterVector>(modVars["params"]);
  CharacterVector trans = modVars["trans"];
  double dfN = -1;
  if (df.isNULL()){
    dfN = -1;
  } else if (rxIs(df,"integer") || rxIs(df,"numeric")){
    NumericVector df0 = as<NumericVector>(df);
    dfN = df0[0];
  }
  double *sigmaD = NULL;
  int sigmaSize = -1;
  Nullable<NumericMatrix> sigma1 = as<Nullable<NumericMatrix>>(sigma);
  getRxModels();
  if (!sigma1.isNull()){
    NumericMatrix sigma2 = NumericMatrix(sigma1);
    sigmaD = &sigma2[0];
    sigmaSize = sigma2.nrow();
    if (!sigma2.hasAttribute("dimnames")){
      stop("The sigma matrix must have named dimensions.");
    }
    List dimnames = sigma2.attr("dimnames");
    StringVector simNames = as<StringVector>(dimnames[1]);
    _rxModels[".simNames"] = simNames;
  }
  // Make sure the model variables are assigned...
  // This fixes random issues on windows where the solves are done and the data set cannot be solved.
  std::string ptrS = (as<std::string>(trans["ode_solver_ptr"]));
  _rxModels[ptrS] = modVars;
  return getSolvingOptionsPtr(atol,rtol,hini, hmin,
			      maxsteps, maxordn, maxords, transit,
			      lhs.size(), state.size(),
			      st, f1, f2, kind, is_locf, cores,
			      ncov,par_cov, do_par_cov, &inits[0], &scale[0],
			      ptrS.c_str(), hmax2, atol2, rtol2, 
			      nDisplayProgress, sigmaD, sigmaSize,
                              dfN, ncoresRV, isChol,svar);
}

SEXP rxSolvingData(const RObject &model,
                   const RObject &parData,
                   const std::string &method = "liblsoda",
                   const Nullable<LogicalVector> &transit_abs = R_NilValue,
		   const double atol = 1.0e-8,
		   const double rtol = 1.0e-6,
		   const int maxsteps = 5000,
		   const double hmin = 0,
		   const Nullable<NumericVector> &hmax = R_NilValue,
		   const double hini = 0,
		   const int maxordn = 12,
		   const int maxords = 5,
		   const int cores = 1,
		   std::string covs_interpolation = "linear",
		   bool addCov = false,
		   bool matrix = false) {
  if (rxIs(parData, "RxODE.par.data")){
    List opt = List(parData);
    DataFrame ids = as<DataFrame>(opt["ids"]);
    IntegerVector BadDose = as<IntegerVector>(opt["BadDose"]);
    NumericVector InfusionRate = as<NumericVector>(opt["InfusionRate"]);
    NumericVector par = as<NumericVector>(opt["pars"]);
    double hm;
    int nPar = as<int>(opt["n.pars"]);
    int nSub = as<int>(opt["nSub"]);
    NumericVector inits = as<NumericVector>(opt["inits"]);
    NumericVector scale = as<NumericVector>(opt["scale"]);
    DataFrame doseDf = as<DataFrame>(opt["dose"]);
    NumericVector amt = as<NumericVector>(doseDf["amt"]);
    IntegerVector idose = as<IntegerVector>(opt["idose"]);
    IntegerVector posDose = as<IntegerVector>(ids["posDose"]);
    IntegerVector posEvent = as<IntegerVector>(ids["posEvent"]);
    IntegerVector posCov = as<IntegerVector>(ids["posCov"]);
    IntegerVector nEvent = as<IntegerVector>(ids["nEvent"]);
    NumericVector solve  = as<NumericVector>(opt["solve"]);
    DataFrame et         = as<DataFrame>(opt["et"]);
    IntegerVector evid   = as<IntegerVector>(et["evid"]);
    NumericVector all_times = as<NumericVector>(et["time"]);
    int totSize = et.nrow();
    NumericVector lhs = as<NumericVector>(opt["lhs"]);
    int lhsSize = as<int>(opt["lhsSize"]);
    IntegerVector par_cov = as<IntegerVector>(opt["pcov"]);
    NumericVector cov = as<NumericVector>(opt["cov"]);
    IntegerVector rc=as<IntegerVector>(ids["rc"]);
    IntegerVector siV = as<IntegerVector>(opt["state.ignore"]);
    int do_par_cov = 0;
    if (par_cov.size() > 0){
      do_par_cov = 1;
    }
    rx_solving_options_ind *inds;
    int nsim = as<int>(opt["nsim"]);
    inds = rxOptionsIniEnsure(nSub*nsim);
    int neq = as<int>(opt["neq"]);
    int ncov =-1;
    int cid;
    for (int simNum = 0; simNum < nsim; simNum++){
      for (int id = 0; id < nSub; id++){
	cid = id+simNum*nSub;
	if (hmax.isNull()){
          // Get from data.
          NumericVector hmn = as<NumericVector>(ids["HmaxDefault"]);
          hm = hmn[id];
        } else {
          NumericVector hmn = NumericVector(hmax);
          if (R_FINITE(hmn[0])){
	    stop("'hmax' must be a non-negative value.");
            hm = hmn[0];
          } else {
            hm = 0.0;
          }
        }
	ncov = par_cov.size();
        getSolvingOptionsIndPtr(&InfusionRate[cid*neq],&BadDose[cid*neq], hm,
				&par[cid*nPar], &amt[posDose[id]],
				&idose[posDose[id]],
                                // Solve and lhs are written to in the solve...
                                &solve[cid*totSize*neq],
                                &lhs[cid*lhsSize],
                                // Doesn't change with the solve.
				&evid[posEvent[id]], &rc[id], &cov[posCov[id]],
                                nEvent[id], &all_times[posEvent[id]], id, simNum,
                                &inds[cid]);
      }
    }
    NumericVector atol2 = as<NumericVector>(opt["atol"]);
    NumericVector rtol2 = as<NumericVector>(opt["rtol"]);
    for (int i = 0; i < atol2.size(); i++){
      atol2[i]=atol;
      rtol2[i]=rtol;
    }
    double hmax2 = as<double>(opt["Hmax"]);
    IntegerVector svar = as<IntegerVector>(opt["svar"]);
    bool isCholB =  as<bool>(opt["isChol"]);
    int isChol = 0;
    if (isCholB) isChol = 1;
    SEXP op = rxSolvingOptions(model,method, transit_abs, atol, rtol, maxsteps, hmin, hini, maxordn,
			       maxords, cores, ncov, &par_cov[0], do_par_cov, &inits[0], &scale[0], covs_interpolation,
			       hmax2,&atol2[0],&rtol2[0], as<int>(opt["nDisplayProgress"]),
			       as<RObject>(opt["sigma"]), as<RObject>(opt["df"]),
			       as<int>(opt["ncoresRV"]),isChol, &svar[0]);

    int add_cov = 0;
    if (addCov) add_cov = 1;
    int nobs = as<int>(opt["nObs"]);
    int mat = 0;
    if (matrix) mat = 1;
    return rxSolveData(inds, nSub, nsim, &siV[0], nobs, add_cov, mat, op);
  } else {
    stop("This requires something setup by 'rxDataParSetup'.");
  }
  return R_NilValue;
}

extern "C" rx_solve *rxSingle(SEXP object, const int stiff,const int transit_abs,
			      const double atol, const double rtol, const int maxsteps,
			      const double hmin, const double hini, const int maxordn,
			      const int maxords, const int cores, const int ncov,
			      int *par_cov, int do_par_cov, 
			      int is_locf,
			      // Other single solve option
			      double hmax, double *par,
			      double *amt, double *solve, double *lhs,
			      int *evid, int *rc, double *cov,
			      int nTimes, double *all_times){
  List mv = rxModelVars(object);
  // Use the number of each element to speed calculation.
  List solveL = mv["solve"];
  NumericVector inits           = solveL["inits"];
  NumericVector scale           = solveL["scale"];
  NumericVector InfusionRate    = solveL["infusion"];
  IntegerVector BadDose         = solveL["badDose"];
  // Instead of having the correct length for idose, use idose length = length of ntime
  // Saves an additional for loop at the cost of a little memory.
  /* int *idose; */
  IntegerVector idose(nTimes);
  rx_solving_options_ind *inds;
  inds = rxOptionsIniEnsure(1);//(rx_solving_options_ind *)Calloc(1, rx_solving_options_ind);
  getSolvingOptionsIndPtr(&InfusionRate[0],&BadDose[0], hmax, par, amt, &idose[0], solve, 
			  lhs, evid, rc, cov, nTimes, all_times, 1, 1, &inds[0]);
  std::string method = "lsoda";
  if (stiff == 0) {
    method = "dop853";
  }
  std::string covs_interpolation = "linear";
  if (is_locf == 1){
    covs_interpolation="constant";
  } else if (is_locf==2){
    covs_interpolation="nocb";
  } else if (is_locf== 3){
    covs_interpolation="midpoint";
  }
  LogicalVector transit_absLV(1);
  if (transit_abs == 1) {
    transit_absLV[0] = true;
  }  else {
    transit_absLV[0] = false;
  }
  SEXP op = rxSolvingOptions(object,method, transit_absLV, atol, rtol, maxsteps, hmin, hini, maxordn,
			     maxords, 1, ncov, par_cov, do_par_cov, &inits[0], &scale[0], covs_interpolation);
  IntegerVector siV = mv["state.ignore"];
  rxSolveData(inds, 1, 1, &siV[0], -1, 0, 0, op);
  SEXP trans = mv["trans"];
  rxUpdateFuns(trans);
  rx_solve *ret = getRxSolve_();
  // Also assign it.
  set_solve(ret);
  return ret;
}

List rxData(const RObject &object,
            const RObject &params = R_NilValue,
            const RObject &events = R_NilValue,
            const RObject &inits = R_NilValue,
            const RObject &covs  = R_NilValue,
            const std::string &method = "liblsoda",
            const Nullable<LogicalVector> &transit_abs = R_NilValue,
            const double atol = 1.0e-8,
            const double rtol = 1.0e-6,
            const int maxsteps = 5000,
            const double hmin = 0,
            const Nullable<NumericVector> &hmax = R_NilValue,
            const double hini = 0,
            const int maxordn = 12,
            const int maxords = 5,
            const int cores = 1,
            std::string covs_interpolation = "linear",
	    bool addCov = false,
            bool matrix = false,
            const RObject &sigma= R_NilValue,
            const RObject &sigmaDf= R_NilValue,
            const int &nCoresRV= 1,
            const bool &sigmaIsChol= false,
            const int &nDisplayProgress = 10000,
            const StringVector &amountUnits = NA_STRING,
            const StringVector &timeUnits = "hours",
	    const RObject &theta = R_NilValue,
            const RObject &eta = R_NilValue,
	    const RObject &scale = R_NilValue,
	    const Nullable<List> &extraArgs = R_NilValue){
  List parData = rxDataParSetup(object,params, events, inits, covs, sigma, sigmaDf,
				nCoresRV, sigmaIsChol, nDisplayProgress, amountUnits, timeUnits,
				theta,eta, scale, extraArgs);
  parData["pointer"] = rxSolvingData(object, parData, method, transit_abs, atol,  rtol, maxsteps,
                                     hmin, hmax,  hini, maxordn, maxords, cores, covs_interpolation,
				     addCov, matrix);
  List modVars = rxModelVars(object);
  StringVector cls(3);
  cls(2) = "RxODE.par.data";
  cls(1) = "RxODE.multi.data";
  cls(0) = "RxODE.pointer.multi";
  parData.attr("class") = cls;
  return parData;
}

#define defrx_params R_NilValue
#define defrx_events R_NilValue
#define defrx_inits R_NilValue
#define defrx_covs R_NilValue
#define defrx_method "liblsoda"
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
#define defrx_covs_interpolation "linear"
#define defrx_addCov false
#define defrx_matrix false
#define defrx_sigma  R_NilValue
#define defrx_sigmaDf R_NilValue
#define defrx_nCoresRV 1
#define defrx_sigmaIsChol false
#define defrx_nDisplayProgress 10000
#define defrx_amountUnits NA_STRING
#define defrx_timeUnits "hours"
#define defrx_addDosing false


RObject rxCurObj;

Nullable<Environment> rxRxODEenv(RObject obj);

std::string rxDll(RObject obj);

bool rxDynLoad(RObject obj);

SEXP rxSolveC(const RObject &object,
              const Nullable<CharacterVector> &specParams = R_NilValue,
	      const Nullable<List> &extraArgs = R_NilValue,
	      const RObject &params = R_NilValue,
	      const RObject &events = R_NilValue,
	      const RObject &inits = R_NilValue,
	      const RObject &scale = R_NilValue,
	      const RObject &covs  = R_NilValue,
	      const CharacterVector &method = "liblsoda",
	      const Nullable<LogicalVector> &transit_abs = R_NilValue,
	      const double atol = 1.0e-8,
	      const double rtol = 1.0e-6,
	      const int maxsteps = 5000,
	      const double hmin = 0,
	      const Nullable<NumericVector> &hmax = R_NilValue,
	      const double hini = 0,
	      const int maxordn = 12,
	      const int maxords = 5,
	      const int cores = 1,
	      const CharacterVector &covs_interpolation = "linear",
	      bool addCov = false,
	      bool matrix = false,
	      const RObject &sigma= R_NilValue,
	      const RObject &sigmaDf= R_NilValue,
	      const int &nCoresRV= 1,
	      const bool &sigmaIsChol= false,
	      const int &nDisplayProgress = 10000,
	      const CharacterVector &amountUnits = NA_STRING,
	      const CharacterVector &timeUnits = "hours",
              const bool addDosing = false,
	      const RObject &theta = R_NilValue,
	      const RObject &eta = R_NilValue,
	      const bool updateObject = false,
	      const bool doSolve = true
              ){
  if (updateObject && !rxIs(object, "rxSolve")){
    return rxSolveC(rxCurObj, specParams, extraArgs, params, events, inits,
                    scale, covs, method, transit_abs, atol, rtol, maxsteps,
                    hmin, hmax, hini, maxordn, maxords, cores,covs_interpolation,
                    addCov, matrix, sigma, sigmaDf, nCoresRV, sigmaIsChol, nDisplayProgress,
                    amountUnits,timeUnits, addDosing, R_NilValue, R_NilValue, updateObject, false);
  } else if (rxIs(object, "rxSolve") || rxIs(object, "environment")){
    // Check to see what parameters were updated by specParams
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
      update_dosing = false;
    if (specParams.isNull()){
      warning("No additional parameters were specified; Returning fit.");
      return object;
    }
    CharacterVector specs = CharacterVector(specParams);
    int n = specs.size(), i;
    for (i = 0; i < n; i++){
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
    }
    // Now update
    Environment e;
    List obj;
    if (rxIs(object, "rxSolve")){
      obj = as<List>(obj);
      CharacterVector classattr = object.attr("class");
      e = as<Environment>(classattr.attr(".RxODE.env"));
    } else if (rxIs(object, "environment")) {
      e = as<Environment>(object);
      obj = as<List>(e["obj"]);
    }
    RObject new_params = update_params ? params : e["args.params"];
    RObject new_events = update_events ? events : e["args.events"];
    RObject new_inits = update_inits ? inits : e["args.inits"];
    RObject new_covs  = update_covs  ? covs  : e["args.covs"];
    CharacterVector new_method = update_method ? method : e["args.method"];
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
    CharacterVector new_covs_interpolation = update_covs_interpolation ? covs_interpolation : e["args.covs_interpolation"];
    bool new_addCov = update_addCov ? addCov : e["args.addCov"];
    bool new_matrix = update_matrix ? matrix : e["args.matrix"];
    RObject new_sigma = update_sigma ? sigma : e["args.sigma"];
    RObject new_sigmaDf = update_sigmaDf ? sigmaDf : e["args.sigmaDf"];
    int new_nCoresRV = update_nCoresRV ? nCoresRV : e["args.nCoresRV"];
    bool new_sigmaIsChol = update_sigmaIsChol ? sigmaIsChol : e["args.sigmaIsChol"];
    int new_nDisplayProgress = e["args.nDisplayProgress"];
    CharacterVector new_amountUnits = update_amountUnits ? amountUnits : e["args.amountUnits"];
    CharacterVector new_timeUnits = update_timeUnits ? timeUnits : e["args.timeUnits"];
    RObject new_scale = update_scale ? scale : e["args.scale"];
    bool new_addDosing = update_dosing ? addDosing : e["args.addDosing"];

    RObject new_object = as<RObject>(e["args.object"]);
    CharacterVector new_specParams(0);
    List dat = as<List>(rxSolveC(new_object, new_specParams, extraArgs, new_params, new_events, new_inits, new_scale, new_covs,
				 new_method, new_transit_abs, new_atol, new_rtol, new_maxsteps, new_hmin,
				 new_hmax, new_hini,new_maxordn, new_maxords, new_cores, new_covs_interpolation,
				 new_addCov, new_matrix, new_sigma, new_sigmaDf, new_nCoresRV, new_sigmaIsChol,
                                 new_nDisplayProgress, new_amountUnits, new_timeUnits, new_addDosing));
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
    if (!rxDynLoad(object)){
      stop("Cannot load RxODE dlls for this model.");
    }
    List parData = rxData(object, params, events, inits, covs, as<std::string>(method[0]), transit_abs, atol,
                          rtol, maxsteps, hmin,hmax, hini, maxordn, maxords, cores,
                          as<std::string>(covs_interpolation[0]), addCov, matrix, sigma, sigmaDf, nCoresRV, sigmaIsChol,
                          nDisplayProgress, amountUnits, timeUnits, theta, eta, scale, extraArgs);
    if (!doSolve){
      // Backwards Compatible; Create solving environment
      if (as<int>(parData["nsim"]) == 1 && as<int>(parData["nSub"]) == 1){
	int stiff = 0;
        if (as<std::string>(method[0]) == "liblsoda"){
	  stiff = 2;
	} else if (as<std::string>(method[0]) == "lsoda"){
	  stiff = 1;
	} else if (as<std::string>(method[0]) != "dop853") {
	  stop("Only lsoda or dop853 can be used with do.solve=FALSE");
	}
	List mv = rxModelVars(object);
	List et = parData["et"];
	List dose = parData["dose"];
        List ret;
	if (!extraArgs.isNull()){
	  ret = as<List>(extraArgs);
	}
        RObject par0 = params;
        RObject ev0  = events;
        RObject par1;
        if (rxIs(par0, "rx.event")){
          // Swapped events and parameters
          par1 = ev0;
        } else if (rxIs(ev0, "rx.event")) {
          par1 = par0;
        } else {
          stop("Need some event information (observation/dosing times) to solve.\nYou can use either 'eventTable' or an RxODE compatible data frame/matrix.");
        }
	NumericVector p = parData["pars"];
	p.attr("names") = mv["params"];
        ret["params"] = p;
	ret["inits"] = parData["inits"]; // named
	ret["scale"] = parData["scale"];
	// add.cov
	ret["state.ignore"] = parData["state.ignore"];
	ret["object"] = rxRxODEenv(mv);
	// event.table
	// events
	// extra.args
	ret["lhs_vars"] = mv["lhs"];
	ret["time"] = et["time"];
	ret["evid"] = et["evid"];
	ret["amt"] = dose["amt"];
	ret["pcov"] = parData["pcov"];
	ret["covs"] = parData["cov"];
	// FIXME: isLocf
        int is_locf = 0;
        if (as<std::string>(covs_interpolation[0]) == "linear"){
        } else if (as<std::string>(covs_interpolation[0]) == "constant" || as<std::string>(covs_interpolation[0]) == "locf" || as<std::string>(covs_interpolation[0]) == "LOCF"){
          is_locf=1;
        } else if (as<std::string>(covs_interpolation[0]) == "nocb" || as<std::string>(covs_interpolation[0]) == "NOCB"){
          is_locf=2;
        }  else if (as<std::string>(covs_interpolation[0]) == "midpoint"){
          is_locf=3;
        } else {
          stop("Unknown covariate interpolation specified.");
        }
	ret["isLocf"] = is_locf;
        ret["atol"] = atol;
        ret["rtol"] = rtol;
	ret["hmin"] = hmin;
        if (hmax.isNull()){
          // Get from data.
	  List ids = as<List>(parData["ids"]);
          NumericVector hmn = as<NumericVector>(ids["HmaxDefault"]);
          ret["hmax"] = hmn[0];
	} else {
	  NumericVector hmn = as<NumericVector>(hmax);
          ret["hmax"] = hmn[0];
	}
	ret["hini"] = hini;
	ret["maxordn"] = maxordn;
	ret["maxords"] = maxords;
	ret["maxsteps"] = maxsteps;
	// FIXME stiff
	ret["stiff"] = stiff;
        int transit = 0;
        if (transit_abs.isNull()){
          transit = mv["podo"];
          if (transit){
            warning("Assumed transit compartment model since 'podo' is in the model.");
          }
        } else {
          LogicalVector tr = LogicalVector(transit_abs);
          if (tr[0]){
            transit=  1;
          }
        }
	ret["transit_abs"] = transit;
	IntegerVector rc(1);
	ret["rc"] = rc;
        return as<SEXP>(ret);
      } else {
	stop("do.solve = TRUE only works with single subject data (currently).");
      }
    }
    DataFrame ret;
    
    rx_solve *rx;
    rx = getRxSolve(parData);
    par_solve(rx);
    rx_solving_options *op = (rx_solving_options*)R_ExternalPtrAddr(rx->op);
    if (op->abort){
      stop("Aborted solve.");
    }
    int doDose = 0;
    if (addDosing){
      doDose = 1;
    } else {
      doDose = 0;
    }
    List dat = RxODE_df(parData, doDose);
    List xtra;
    if (!rx->matrix) xtra = RxODE_par_df(parData);
    int nr = as<NumericVector>(dat[0]).size();
    int nc = dat.size();
    if (rx->matrix){
      dat.attr("class") = "data.frame";
      NumericMatrix tmpM(nr,nc);
      for (int i = 0; i < dat.size(); i++){
        tmpM(_,i) = as<NumericVector>(dat[i]);
      }
      tmpM.attr("dimnames") = List::create(R_NilValue,dat.names());
      return tmpM;
    } else {
      Function newEnv("new.env", R_BaseNamespace);
      Environment RxODE("package:RxODE");
      Environment e = newEnv(_["size"] = 29, _["parent"] = RxODE);
      e["check.nrow"] = nr;
      e["check.ncol"] = nc;
      e["check.names"] = dat.names();
      // Save information
      // Remove one final; Just for debug.
      // e["parData"] = parData;
      List pd = as<List>(xtra[0]);
      if (pd.size() == 0){
	e["params.dat"] = R_NilValue;
      } else {
	e["params.dat"] = pd;
      }
      if (as<int>(parData["nSub"]) == 1 && as<int>(parData["nsim"]) == 1){
        int n = pd.size();
        NumericVector par2(n);
        for (int i = 0; i <n; i++){
          par2[i] = (as<NumericVector>(pd[i]))[0];
        }
        par2.names() = pd.names();
	if (par2.size() == 0){
	  e["params.single"] = R_NilValue;
	} else {
	  e["params.single"] = par2;
        }
      } else {
        e["params.single"] = R_NilValue;
      }
      e["EventTable"] = xtra[1];
      e["dosing"] = xtra[3];
      e["sampling"] = xtra[2];
      e["obs.rec"] = xtra[4];
      e["covs"] = xtra[5];
      e["counts"] = xtra[6];
      e["inits.dat"] = parData["inits"];
      CharacterVector units(2);
      units[0] = as<std::string>(parData["amount.units"]);
      units[1] = as<std::string>(parData["time.units"]);
      CharacterVector unitsN(2);
      unitsN[0] = "dosing";
      unitsN[1] = "time";
      units.names() = unitsN;
      e["units"] = units;
      e["nobs"] = parData["nObs"];
    
      Function eventTable("eventTable",RxODE);
      List et = eventTable(_["amount.units"] = as<std::string>(parData["amount.units"]), _["time.units"] =as<std::string>(parData["time.units"]));
      Function importEt = as<Function>(et["import.EventTable"]);
      importEt(e["EventTable"]);
      e["events.EventTable"] = et;
      Function parse2("parse", R_BaseNamespace);
      Function eval2("eval", R_BaseNamespace);
      // eventTable style methods
      e["get.EventTable"] = eval2(_["expr"]   = parse2(_["text"]="function() EventTable"),
                                  _["envir"]  = e);
      e["get.obs.rec"] = eval2(_["expr"]   = parse2(_["text"]="function() obs.rec"),
                               _["envir"]  = e);
      e["get.nobs"] = eval2(_["expr"]   = parse2(_["text"]="function() nobs"),
                            _["envir"]  = e);
      e["add.dosing"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {et <- create.eventTable(); et$add.dosing(...); invisible(rxSolve(args.object,events=et,update.object=TRUE))}"),
                              _["envir"]  = e);
      e["clear.dosing"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {et <- create.eventTable(); et$clear.dosing(...); invisible(rxSolve(args.object,events=et,update.object=TRUE))}"),
                                _["envir"]  = e);
      e["get.dosing"] = eval2(_["expr"]   = parse2(_["text"]="function() dosing"),
                              _["envir"]  = e);

      e["add.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {et <- create.eventTable(); et$add.sampling(...); invisible(rxSolve(args.object,events=et,update.object=TRUE))}"),
                                _["envir"]  = e);
      
      e["clear.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {et <- create.eventTable(); et$clear.sampling(...); invisible(rxSolve(args.object,events=et,update.object=TRUE))}"),
                                  _["envir"]  = e);

      e["replace.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function(...) {et <- create.eventTable(); et$clear.sampling(); et$add.sampling(...); invisible(rxSolve(args.object,events=et,update.object=TRUE))}"),
                                _["envir"]  = e);

      e["get.sampling"] = eval2(_["expr"]   = parse2(_["text"]="function() sampling"),
				_["envir"]  = e);
      
      e["get.units"] = eval2(_["expr"]   = parse2(_["text"]="function() units"),
                             _["envir"]  = e);

      e["import.EventTable"] = eval2(_["expr"]   = parse2(_["text"]="function(imp) {et <- create.eventTable(imp); invisible(rxSolve(args.object,events=et,update.object=TRUE))}"),
				     _["envir"]  = e);
      
      e["create.eventTable"] = eval2(_["expr"]   = parse2(_["text"]="function(new.event) {et <- eventTable(amount.units=units[1],time.units=units[2]);if (missing(new.event)) {nev <- EventTable; } else {nev <- new.event;}; et$import.EventTable(nev); return(et);}"),
                                     _["envir"]  = e);
      // Note event.copy doesn't really make sense...?  The create.eventTable does basically the same thing.
      e["args.object"] = object;
      e["dll"] = rxDll(object);
      if (rxIs(events, "rx.event")){
	e["args.params"] = params;    
        e["args.events"] = events;
      } else {
	e["args.params"] = events;    
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
      e[".real.update"] = true;
      CharacterVector cls(2);
      cls(0) = "rxSolve";
      cls(1) = "data.frame";
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
                  opts[0], // const CharacterVector &method = "lsoda",
                  opts[1], // const Nullable<LogicalVector> &transit_abs = R_NilValue,
                  opts[2], //const double atol = 1.0e-8,
                  opts[3], // const double rtol = 1.0e-6,
                  opts[4], //const int maxsteps = 5000,
                  opts[5], //const double hmin = 0,
                  opts[6], //const Nullable<NumericVector> &hmax = R_NilValue,
                  opts[7], //const double hini = 0,
                  opts[8], //const int maxordn = 12,
                  opts[9], //const int maxords = 5,
                  opts[10], //const int cores = 1,
                  opts[11], //const CharacterVector &covs_interpolation = "linear",
                  opts[12], //bool addCov = false,
                  opts[13], //bool matrix = false,
                  opts[14], //const RObject &sigma= R_NilValue,
                  opts[15], //const RObject &sigmaDf= R_NilValue,
                  opts[16], //const int &nCoresRV= 1,
                  opts[17], //const bool &sigmaIsChol= false,
                  opts[18], // nDisplayProgress
                  opts[19], //const CharacterVector &amountUnits = NA_STRING,
                  opts[20], //const CharacterVector &timeUnits = "hours",
                  opts[21], //const RObject &theta = R_NilValue,
                  opts[22], //const RObject &eta = R_NilValue,
                  opts[23], //const bool addDosing = false
		  opts[24],
		  opts[25]);//const bool updateObject = false)
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
	  CharacterVector retS = CharacterVector::create(mods["model"]);
	  retS.attr("class") = "RxODE.modeltext";
	  return(retS);
	}
	if (e.exists(sarg)){
	  return e[sarg];
	}
	if (sarg == "params" || sarg == "par" || sarg == "pars" || sarg == "param"){
	  return e["params.dat"];
	} else if (sarg == "inits" || sarg == "init"){
	  return e["inits.dat"];
	} else if (sarg == "t"){
	  return lst["time"];
	}
	// Now parameters
	List pars = List(e["params.dat"]);
	CharacterVector nmp = pars.names();
	n = pars.size();
	for (i = 0; i < n; i++){
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
        for (i = 0; i < n; i++){
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
	// // Sensitivities -- last
	// // This is slower, defer to last.
	// for (i = 0; i < n; i++){
	//   // The regular expression came from rex;  It is a it long...
	//   if (as<std::string>(sub("rx__sens_((?:[a-zA-Z][_a-zA-Z0-9.]*|(?:\\.){1,}[_a-zA-Z][_a-zA-Z0-9.]*))_BY_((?:[a-zA-Z][_a-zA-Z0-9.]*|(?:\\.){1,}[_a-zA-Z][_a-zA-Z0-9.]*))__",
	//                              "_sens_\\1_\\2", nm[i])) == sarg){
	//        return lst[i];
	//   }
	// }
      }
    } else {
      if (rxIs(arg, "integer") || rxIs(arg, "numeric")){
	int iarg = as<int>(arg);
	if (iarg < lst.size()){
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
                          defrx_timeUnits, defrx_addDosing);
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
			  defrx_addDosing);
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
			  defrx_addDosing);
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
			  defrx_addDosing);
	} else if (sarg == "t" || sarg == "time"){
	  CharacterVector classattr = obj.attr("class");
          Environment e = as<Environment>(classattr.attr(".RxODE.env"));
	  Function f = as<Function>(e["replace.sampling"]);
	  return f(value);
        } else {
	  CharacterVector classattr = obj.attr("class");
	  Environment e = as<Environment>(classattr.attr(".RxODE.env"));
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
	  for (i = 0; i < n; i++){
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
				defrx_addDosing);
	      } else if (val.size() == nc){
		// Change parameter -> Covariate
		List newPars(pars.size()-1);
		CharacterVector newParNames(pars.size()-1);
		for (j = 0; j < i; j++){
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
		for (j = 0; j < covs.size(); j++){
		  newCovs[j]      = covs[j];
		  newCovsNames[j] = nmc[j];
		}
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
				defrx_addDosing);
	      }
	      return R_NilValue;
	    }
	  }
	  ///////////////////////////////////////////////////////////////////////////////
	  // Update Covariates by covariate name
	  n = covs.size();
	  for (i = 0; i < n; i++){
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
				defrx_addDosing);
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
				defrx_addDosing);
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
			      defrx_addDosing);
	    }
	  }
	  return R_NilValue;
	}
      }
    }
  }
  return R_NilValue;
}

extern "C" void rxAddModelLib(SEXP mv){
  getRxModels();
  CharacterVector trans = as<List>(mv)["trans"];
  std::string ptr =as<std::string>(trans["ode_solver_ptr"]);
  _rxModels[ptr]= mv;
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
extern "C" void RxODE_assign_fn_pointers_(SEXP mv, int addit);
//' Assign pointer based on model variables
//' @param object RxODE family of objects
//' @export
//[[Rcpp::export]]
void rxAssignPtr(SEXP object = R_NilValue){
  List mv=rxModelVars(as<RObject>(object));
  RxODE_assign_fn_pointers_(as<SEXP>(mv), 0);
  CharacterVector trans = mv["trans"];
  rxUpdateFuns(as<SEXP>(trans));
  rx_solve *ret = getRxSolve_();
  // Also assign it.
  set_solve(ret); 
  // Update rxModels environment.
  getRxModels();
  std::string ptr = as<std::string>(trans["ode_solver_ptr"]); 
  if (!_rxModels.exists(ptr)){
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
    Environment RxODE("package:RxODE");
    Function f = as<Function>(RxODE["rxCompile.character"]);
    RObject newO = f(as<std::string>(obj));
    return(rxDll(newO));
  } else {
    List mv = rxModelVars(obj);
    Nullable<Environment> en = rxRxODEenv(mv);
    if (en.isNull()){
      stop("Can't figure out the DLL for this object");
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
std::string rxC(RObject obj){
  if (rxIs(obj,"RxODE")){
    Environment e = as<Environment>(obj);
    return as<std::string>((as<List>(e["rxDll"]))["c"]);
  } else if (rxIs(obj,"rxSolve")) {
    CharacterVector cls = obj.attr("class");
    Environment e = as<Environment>(cls.attr(".RxODE.env"));
    return(as<std::string>(e["c"]));
  } else if (rxIs(obj, "rxDll")){
    return as<std::string>(as<List>(obj)["c"]);
  } else if (rxIs(obj, "character")){
    Environment RxODE("package:RxODE");
    Function f = as<Function>(RxODE["rxCompile.character"]);
    RObject newO = f(as<std::string>(obj));
    return(rxDll(newO));
  } else {
    List mv = rxModelVars(obj);
    Nullable<Environment> en = rxRxODEenv(mv);
    if (en.isNull()){
      stop("Can't figure out the DLL for this object");
    } else {
      Environment e = as<Environment>(en);
      return as<std::string>((as<List>(e["rxDll"]))["dll"]);
    }
  }
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
  Function isLoaded("is.loaded", R_BaseNamespace);
  List mv = rxModelVars(obj);
  CharacterVector trans = mv["trans"];
  std::string dydt = as<std::string>(trans["ode_solver"]);
  return as<bool>(isLoaded(dydt));
}

// 
inline bool fileExists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
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
      Function dynLoad("dyn.load", R_BaseNamespace);
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
  std::string ptr = as<std::string>(trans["ode_solver_ptr"]);
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
    std::string cfile = rxC(obj);
    if (fileExists(cfile)) remove(cfile.c_str());
    if (!fileExists(file)) return true;
    if (remove(file.c_str()) == 0) return true;
  }
  return false;
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
//' @param nSub Number between subject variabilities (ETAs) simulated.
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
//' @param nStud Number virtual studies to characterize uncertainty in fixed parameters.
//'
//' @param sigma Matrix for residual variation.  Adds a "NA" value for each of the 
//'     indivdual parameters, residuals are updated after solve is completed. 
//'
//' @inheritParams rxSolve
//'
//' @author Matthew L.Fidler
//'
//' @export
//[[Rcpp::export]]
List rxSimThetaOmega(const Nullable<NumericVector> &params    = R_NilValue,
		     const Nullable<NumericMatrix> &omega= R_NilValue,
		     const Nullable<NumericMatrix> &omegaDf= R_NilValue,
		     const bool &omegaIsChol = false,
		     int nSub = 1,
		     const Nullable<NumericMatrix> &thetaMat = R_NilValue,
		     const Nullable<NumericMatrix> &thetaDf  = R_NilValue,
		     const bool &thetaIsChol = false,
		     int nStud = 1,
                     const Nullable<NumericMatrix> &sigma = R_NilValue,
		     int nCoresRV = 1){
  NumericVector par;
  if (params.isNull()){
    stop("This function requires overall parameters.");
  } else {
    par = NumericVector(params);
    if (!par.hasAttribute("names")){
      stop("Parameters must be named.");
    }
  }
  
  NumericMatrix thetaM;
  CharacterVector thetaN;
  if (!thetaMat.isNull() && nStud > 0){
    thetaM = as<NumericMatrix>(rxSimSigma(as<RObject>(thetaMat), as<RObject>(thetaDf), nCoresRV, thetaIsChol, nStud));
    thetaN = as<CharacterVector>((as<List>(thetaM.attr("dimnames")))[1]);
  }

  NumericMatrix omegaM;
  CharacterVector omegaN;
  if (!omega.isNull() && nSub*nStud > 0){
    omegaM = as<NumericMatrix>(rxSimSigma(as<RObject>(omega), as<RObject>(omegaDf), nCoresRV, omegaIsChol, nSub*nStud));
    omegaN = as<CharacterVector>((as<List>(omegaM.attr("dimnames")))[1]);
  }
  // Now create data frame of parameter values
  List ret;
  int i, j, k;
  CharacterVector parN = CharacterVector(par.attr("names"));
  IntegerVector parI(parN.size());
  NumericVector tmpNV;
  int parNum = -1;
  for (i = 0; i < parN.size(); i++){
    parNum = -1;
    for (j = 0; j < thetaN.size(); j++){
      if (parN[i] == thetaN[j]){
	parNum = j;
	break;
      }
    }
    tmpNV = NumericVector(nSub*nStud);
    if (parNum == -1){
      for (j = 0; j < nSub*nStud; j++){
	tmpNV(j) = par(i);
      }
    } else {
      for (j = 0; j < nStud; j++){
	for (k = 0; k < nSub; k++){
	  tmpNV(j*nSub+k) = par(i) +  thetaM(j,parNum);
	}
      }
    }
    ret[as<std::string>(parN[i])] = tmpNV;
  }
  for (i = 0; i < omegaN.size(); i++){
    tmpNV = NumericVector(nSub*nStud);
    for (j = 0; j < nSub*nStud; j++){
      tmpNV[j] = omegaM(j,i);
    }
    ret[as<std::string>(omegaN[i])] = tmpNV;
  }
  CharacterVector sigmaN;
  NumericMatrix sigmaM;
  if (!sigma.isNull()){
    sigmaM = as<NumericMatrix>(sigma);
    sigmaN = as<CharacterVector>((as<List>(sigmaM.attr("dimnames")))[1]);
  }
  for (i = 0; i < sigmaN.size(); i++){
    tmpNV = NumericVector(nSub*nStud);
    ret[as<std::string>(sigmaN[i])] = tmpNV;
  }
  ret.attr("class") = "data.frame";
  ret.attr("row.names") = IntegerVector::create(NA_INTEGER,-nSub*nStud);
  return ret;
}

SEXP rxGetFromChar(const char *ptr, std::string var){
  std::string str(ptr);
  CharacterVector cv(1);
  cv[0] = str;
  List mv = rxModelVars(as<RObject>(cv));
  return wrap(mv[var]);
}

extern "C" SEXP rxStateNames(const char *ptr){
  return rxGetFromChar(ptr, "state");
}

extern "C" SEXP rxLhsNames(const char *ptr){
  return rxGetFromChar(ptr, "lhs");
}

extern "C" SEXP rxParamNames(const char *ptr){
  return rxGetFromChar(ptr, "params");
}
