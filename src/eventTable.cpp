#include <Rcpp.h>
using namespace Rcpp;

bool rxIs(const RObject &obj, std::string cls);
void updateSolveEnvPost(Environment e);

// [[Rcpp::export]]
RObject add_dosing_ (RObject eventTable, double dose,
		     int nbr_doses = 1,
		     double dosing_interval = 24,
		     int dosing_to=1,
		     Nullable<NumericVector> rate=R_NilValue , 
		     CharacterVector amount_units = NA_STRING,
		     double start_time = 0.0,
		     bool do_sampling=false,
		     CharacterVector time_units = NA_STRING){
  if (rxIs(eventTable, "EventTable")){
    List et = as<List>(eventTable);
    Function f = as<Function>(et["add.dosing"]);
    f(_["dose"] = dose, 
      _["nbr.doses"] = nbr_doses, 
      _["dosing.interval"] = dosing_interval, 
      _["dosing.to"] = dosing_to, 
      _["rate"] = rate, 
      _["amount.units"] = amount_units, 
      _["start.time"]  = start_time, 
      _["do.sampling"] = do_sampling, 
      _["time.units"] =time_units);
    return eventTable;
  } else if (rxIs(eventTable, "rxSolve")){
    CharacterVector classattr = eventTable.attr("class");
    Environment e = as<Environment>(classattr.attr(".RxODE.env"));
    updateSolveEnvPost(e);
    e[".real.update"]=false;
    Function f = as<Function>(e["add.dosing"]);
    RObject ret = as<RObject>(f(_["dose"] = dose, 
      _["nbr.doses"] = nbr_doses, 
      _["dosing.interval"] = dosing_interval, 
      _["dosing.to"] = dosing_to, 
      _["rate"] = rate, 
      _["amount.units"] = amount_units, 
      _["start.time"]  = start_time, 
      _["do.sampling"] = do_sampling, 
      _["time.units"] =time_units));
    e[".real.update"]=true;
    return ret;
  }
  return R_NilValue;
}

// [[Rcpp::export]]
RObject add_sampling_ (RObject eventTable, 
		       NumericVector time, 
		       CharacterVector time_units = NA_STRING){
  if (rxIs(eventTable, "EventTable")){
    List et = as<List>(eventTable);
    Function f = as<Function>(et["add.sampling"]);
    f(_["time"] = time, 
      _["time.units"] =time_units);
    return eventTable;
  } else if (rxIs(eventTable, "rxSolve")){
    CharacterVector classattr = eventTable.attr("class");
    Environment e = as<Environment>(classattr.attr(".RxODE.env"));
    updateSolveEnvPost(e);
    e[".real.update"]=false;
    Function f = as<Function>(e["add.sampling"]);
    RObject ret = as<RObject>(f(_["time"] = time, 
				_["time.units"] =time_units));
    e[".real.update"]=true;
    return ret;
  }
  return R_NilValue;
}
