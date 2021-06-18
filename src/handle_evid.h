#ifndef __HANDLE_EVID_H___
#define __HANDLE_EVID_H___

#include "../inst/include/RxODE.h"

#if defined(__cplusplus)
#define FLOOR(x) std::floor(x)
extern "C" {
#else
#define FLOOR(x) floor(x)
#endif

  int handle_evidL(int evid, double *yp, double xout, int id, rx_solving_options_ind *ind);
  void handleTlast(double *time, rx_solving_options_ind *ind);
  
#if defined(__cplusplus)
}
#endif  


static inline void getWh(int evid, int *wh, int *cmt, int *wh100, int *whI, int *wh0){
  *wh = evid;
  *cmt = 0;
  *wh100 = FLOOR(*wh/1e5L);
  *whI   = FLOOR(*wh/1e4L-*wh100*10);
  *wh    = *wh - *wh100*1e5 - (*whI-1)*1e4;
  *wh0 = FLOOR((*wh%10000)/100);
  *cmt = *wh0 - 1 + *wh100*100;
  *wh0 = evid - *wh100*1e5 - *whI*1e4 - *wh0*100;
  if (rx_global.linNcmt != 0) {
    if (rx_global.linKa) {
      switch (*cmt) {
      case 0:
	*cmt = op_global.neq;
	break;
      case 1:
	*cmt = op_global.neq+1;
	break;
      case 2:
	*cmt -= 2;
	break;
      }
    } else {
      if (*cmt == 0) {
	*cmt = op_global.neq;
      } else {
	*cmt -= 1;
      }
    }
  }
}

static inline void handleTlastInline(double *time, rx_solving_options_ind *ind) {
  rx_solving_options *op = &op_global;
  double _time = *time + ind->curShift;
  if (op->neq + op->extraCmt != 0 && ind->tlast != _time && isDose(ind->evid[ind->ix[ind->idx]]) &&
      ind->cmt < op->neq + op->extraCmt){
    ind->dosenum++;
    ind->tlast = _time;
    if (ISNA(ind->tfirst)) ind->tfirst = _time;
    ind->tlastS[ind->cmt] = _time;
    if (ISNA(ind->tfirstS[ind->cmt])) ind->tfirstS[ind->cmt] = _time;
  }
}



static inline int syncIdx(rx_solving_options_ind *ind) {
  if (ind->ix[ind->idx] != ind->idose[ind->ixds]) {
    // bisection https://en.wikipedia.org/wiki/Binary_search_algorithm
    int l = 0, r = ind->ndoses-1, m=0, successful = 0;
    while(l <= r){
      m = FLOOR((l+r)/2);
      if (ind->idose[m] < ind->ix[ind->idx]) l = m+1;
      else if (ind->idose[m] > ind->ix[ind->idx]) r = m-1;
      else {
	successful = 1;
	break;
      }
    }
    if (!successful) {
      //262144
      if (!(ind->err & 262144)){
	ind->err += 262144;
      }
      return 0;
    }
    if (ind->idose[m] == ind->ix[ind->idx]) {
      if (ind->evid[ind->ix[ind->idx]] != ind->evid[ind->ixds]) {
	//262144
	if (!(ind->err & 262144)){
	  ind->err += 262144;
	}
	return 0;
      } else {
	ind->ixds=m;
      }
    } else {
      //262144
      if (!(ind->err & 262144)){
	ind->err += 262144;
      }
      return 0;
      /* Rf_errorcall(R_NilValue, "Corrupted event table; EVID=%d: %d %d %d", evid, ind->idose[m], ind->ix[ind->idx], */
      /* 	ind->idx); */
    }
    // Need to adjust ixdsr
    for(int j = ind->ixds; j--;){
      if (ind->ix[ind->idx] == ind->idose[j]){
	ind->ixds = j;
	break;
      }
    }
    if (ind->ix[ind->idx] != ind->idose[ind->ixds]){
      for(int j = ind->ixds+1; j< ind->ndoses; j++){
	if (ind->ix[ind->idx] == ind->idose[j]){
	  ind->ixds = j;
	  break;
	}
      }
    }
    if (ind->ix[ind->idx] != ind->idose[ind->ixds]){
      //524288
      if (!(ind->err & 524288)){
	ind->err += 524288;
      }
      return 0;
    }
  }
  return 1;
}

extern t_F AMT;

static inline double getAmt(rx_solving_options_ind *ind, int id, int cmt, double dose, double t, double *y){
  double ret = AMT(id, cmt, dose, t, y);
  if (ISNA(ret)){
    rx_solving_options *op = &op_global;
    op->badSolve=1;
    op->naTime = 1;
  }
  return ret;
}
 
static inline int handle_evid(int evid, int neq, 
			      int *BadDose,
			      double *InfusionRate,
			      double *dose,
			      double *yp,
			      int do_transit_abs,
			      double xout, int id,
			      rx_solving_options_ind *ind) {
  if (isObs(evid)) return 0;
  int cmt, foundBad, j;
  double tmp;
  getWh(evid, &(ind->wh), &(ind->cmt), &(ind->wh100), &(ind->whI), &(ind->wh0));
  handleTlastInline(&xout, ind);
  if (ind->wh0 == 40){
    ind->ixds++;
    return 1;
  }
  /* wh100 = ind->wh100; */
  cmt = ind->cmt;
  if (cmt<0) {
    if (!(ind->err & 65536)){
      ind->err += 65536;
      /* Rprintf("Supplied an invalid EVID (EVID=%d; cmt %d)", evid, cmt); */
    }
    return 0;
  }
  if (cmt >= neq){
    foundBad = 0;
    for (j = 0; j < ind->nBadDose; j++){
      if (BadDose[j] == cmt+1){
	foundBad=1;
	break;
      }
    }
    if (!foundBad){
      BadDose[ind->nBadDose]=cmt+1;
      ind->nBadDose++;
    }
  } else {
    rx_solving_options *op = &op_global;
    
    REprintf("ind->ix[ind->idx]: %d ind->ixds: %d\n", ind->ix[ind->idx], ind->ixds);
    if (syncIdx(ind) == 0) return 0;
    REprintf("\tind->ix[ind->idx]: %d ind->ixds: %d\n", ind->ix[ind->idx], ind->ixds);
    if (ind->wh0 == 30){
      yp[cmt]=op_global.inits[cmt];
      InfusionRate[cmt] = 0;
      ind->cacheME=0;
      ind->on[cmt] = 0;
      return 1;
    }
    if (!ind->doSS && ind->wh0 == 20 && cmt < op->neq){
      // Save for adding at the end; Only for ODE systems
      memcpy(ind->solveSave, yp, op->neq*sizeof(double));
    }
    switch(ind->whI){
    case 9: // modeled rate.
    case 8: // modeled duration.
      // Rate already calculated and saved in the next dose record
      ind->on[cmt] = 1;
      ind->cacheME=0;
      REprintf("Turn on infusion at cmt: %d t: %f rate: %f; %d\n", cmt, xout, dose[ind->ixds+1], ind->ixds+1);
      InfusionRate[cmt] -= dose[ind->ixds+1];
      if (ind->wh0 == 20 && getAmt(ind, id, cmt, dose[ind->ixds], xout, yp) != dose[ind->ixds]){
	if (!(ind->err & 1048576)){
	  ind->err += 1048576;
	}
	return 0;
	/* Rf_errorcall(R_NilValue, "SS=2 & Modeled F does not work"); */
      }
      break;
    case 7: // End modeled rate
    case 6: // end modeled duration
      // In this case re-sort is not going to be assessed
      // If cmt is off, don't remove rate....
      // Probably should throw an error if the infusion rate is on still.
      REprintf("Turn OFF infusion at cmt: %d t: %f rate: %f; %d\n", cmt, xout, dose[ind->ixds], ind->ixds);
      InfusionRate[cmt] += dose[ind->ixds];
      ind->cacheME=0;
      if (ind->wh0 == 20 &&
	  getAmt(ind, id, cmt, dose[ind->ixds], xout, yp) !=
	  dose[ind->ixds]){
	if (!(ind->err & 2097152)){
	  ind->err += 2097152;
	}
	return 0;
      }
      break;
    case 2:
      // In this case bio-availability changes the rate, but the
      // duration remains constant.  rate = amt/dur
      ind->on[cmt] = 1;
      tmp = getAmt(ind, id, cmt, dose[ind->ixds], xout, yp);
      InfusionRate[cmt] += tmp;
      ind->cacheME=0;
      if (ind->wh0 == 20 && tmp != dose[ind->ixds]){
	if (!(ind->err & 4194304)){
	  ind->err += 4194304;
	}
	return 0;
      }
      break;
    case 1:
      ind->on[cmt] = 1;
      InfusionRate[cmt] += dose[ind->ixds];
      ind->cacheME=0;
      if (ind->wh0 == 20 && dose[ind->ixds] > 0 && getAmt(ind, id, cmt, dose[ind->ixds], xout, yp) != dose[ind->ixds]){
	if (!(ind->err & 4194304)){
	  ind->err += 4194304;
	}
      }
      break;
    case 4: // replace
      ind->on[cmt] = 1;
      ind->podo = 0;
      handleTlastInline(&xout, ind);
      yp[cmt] = getAmt(ind, id, cmt, dose[ind->ixds], xout, yp);     //dosing before obs
      break;
    case 5: //multiply
      ind->on[cmt] = 1;
      ind->podo = 0;
      handleTlastInline(&xout, ind);
      yp[cmt] *= getAmt(ind, id, cmt, dose[ind->ixds], xout, yp);     //dosing before obs
      break;
    case 0:
      if (do_transit_abs) {
	ind->on[cmt] = 1;
	if (ind->wh0 == 20){
	  tmp = getAmt(ind, id, cmt, dose[ind->ixds], xout, yp);
	  ind->podo = tmp;
	} else {
	  ind->podo = getAmt(ind, id, cmt, dose[ind->ixds], xout, yp);
	}
	handleTlastInline(&xout, ind);
      } else {
	ind->on[cmt] = 1;
	ind->podo = 0;
	handleTlastInline(&xout, ind);
	yp[cmt] += getAmt(ind, id, cmt, dose[ind->ixds], xout, yp);     //dosing before obs
      }
    }
    ind->ixds++;
    ind->solved = ind->idx;
    return 1;
  }
  return 0;
}

static inline int handleEvid1(int *i, rx_solve *rx, int *neq,
			      double *yp, double *xout) {
  rx_solving_options_ind *ind = &(rx->subjects[neq[1]]);
  rx_solving_options *op = rx->op;
  ind->idx = *i;
  return handle_evid(ind->evid[ind->ix[ind->idx]], neq[0] + op->extraCmt,
		     ind->BadDose, ind->InfusionRate, ind->dose, yp,
		     op->do_transit_abs, *xout, neq[1], ind);
}

#endif

