// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Setup a data frame for solving multiple subjects at once in RxODE.
//'
//' @param df dataframe to setup; Must be in RxODE compatible format.
//' @param covNames Covariate names in dataset.
//'
//' @return A data structure to allow C-based for loop (ie solving each
//'       individual in C)
//'
//' @export
// [[Rcpp::export]]
List rxDataSetup(const DataFrame &df, const Nullable<StringVector> &covNames = R_NilValue){
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
    
  int m = 0;
  for (i = 0; i < ids; i++){
    if (lastId != id[i]){
      lastId     = id[i];
      newId[m]   = id[i];
      posDose[m] = j;
      posObs[m]  = k;
      if (m != 0){
        nDose[m-1] = nDoses;
        nObsN[m-1]  = nObs;
      }
      nDoses = 0;
      nObs = 0;
      m++;
    }
    if (evid[i]){
      // Dose
      newEvid[j]  = evid[i];
      newTimeA[j] = time0[i];
      newAmt[j]   = amt[i];
      nDoses++;
      j++;
    } else {
      // Observation
      newDv[k]    = dv[i];
      newTimeO[k] = time0[i];
      nObs++;
      k++;
    }
  }
  nDose[m-1]=nDoses;
  nObsN[m-1]=nObs;

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
                                                        _["nDose"]   = nDose,
                                                        _["nObs"]    = nObsN,
                                                        _["nCov"]    = nCov),
                          _["cov"]=newCov,
                          _["nSub"]=nSub,
                          _["nDoses"]=newEvid.size(),
                          _["nObs"]=newDv.size(),
                          _["cov.names"]=covNames);
  ret.attr("class") = "RxODE.multi.data";
  return ret;
}

//[[Rcpp::export]]
List rxEventTableExpand(const int &nsub,const DataFrame &df){
  // Purpose: Expand current event table to have number of subjects,
  // and then return rxDataSetup list.
  //IntegerVector id    = df["id"];
  IntegerVector evid  = df["evid"];
  //NumericVector dv    = df["dv"];
  NumericVector time0 = df["time"];
  NumericVector amt   = df["amt"];
  int newdim = nsub*amt.size();
  IntegerVector nId(newdim);
  IntegerVector nEvid(newdim);
  NumericVector nDv(newdim);
  NumericVector nTime(newdim);
  NumericVector nAmt(newdim);
  for (int i = 0; i < amt.size(); i++){
    for (int j = 0; j < nsub; j++){
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
  return rxDataSetup(ndf, R_NilValue);
}
