#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
List nmData(DataFrame inData){
  CharacterVector lName = as<CharacterVector>(inData.attr("names"));
  int i, idCol = -1, evidCol=-1, timeCol=-1, amtCol=-1, cmtCol=-1,
    dvCol=-1, ssCol=-1, mdvCol=-1, dateCol=-1, dat1Col=-1, dat2Col=-1,
    dat3Col=-1, rateCol=-1, addlCol=-1, iiCol=-1;
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
    else if (tmpS == "dv" || tmpS == "y") dvCol=-1;
    else if (tmpS == "ss")   ssCol=-1;
    else if (tmpS == "mdv")  mdvCol=-1;
    else if (tmpS == "date") dateCol=-1;
    else if (tmpS == "rate") rateCol=-1;
    else if (tmpS == "addl") addlCol=-1;
    else if (tmpS == "ii")   iiCol=-1;
  }
  // EVID = 0; Observations
  // EVID = 1 is illegal, but converted from NONMEM
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
  // Steady state events are not yet required, but would need a II data item.
  std::vector<int> id;
  std::vector<int> evid;
  std::vector<double> time;
  std::vector<double> amt;
  std::vector<double> ii;
  List lst(1);
  lst[0] = lName;
  return lst;
}
