#include "getTime.h"

extern "C" double getTime(int idx, rx_solving_options_ind *ind) {
  rx_solving_options *op = &op_global;
  rx_solve *rx = &rx_global;
  int evid = ind->evid[idx];
  if (evid == 9) return 0.0;
  if (evid >= 10 && evid <= 99) return ind->mtime[evid-10];
  if (isObs(evid))  return ind->all_times[idx];
  double ret;
  getWh(evid, &(ind->wh), &(ind->cmt), &(ind->wh100), &(ind->whI), &(ind->wh0));
  double *yp;
  if (ind->wh0 == EVID0_SSINF){
  } else {
    // yp should be the current solve values
    //
    // Before solving the solve will be zero
    // After solving the yp will contain the solved values
    if (ind->idx < idx){
      yp = getSolve(ind->idx);
    } else {
      yp = getSolve(idx);
    }
    switch(ind->whI){
    case EVIDF_MODEL_DUR_OFF:
      if (idx > 0){
	yp = rx->ypNA;
	int wh, cmt, wh100, whI, wh0;
	getWh(ind->evid[idx-1], &wh, &cmt, &wh100, &whI, &wh0);
	if (whI != 8){
	  if (!(ind->err & 64)){
	    ind->err += 64;
	  }
	  return 0.0;
	  /* Rf_errorcall(R_NilValue, "Data error 686 (whI = %d; evid=%d)", whI, ind->evid[idx-1]); */
	}
	updateDur(idx-1, ind, yp);
      } else {
	if (!(ind->err & 128)){
	  ind->err += 128;
	}
	return 0.0;
	/* Rf_errorcall(R_NilValue, "Data Error -6\n"); */
      }
      break;
    case EVIDF_MODEL_DUR_ON:
      if (idx >= ind->n_all_times){
	// error: Last record, can't be used.
	if (!(ind->err & 256)){
	  ind->err += 256;
	}
	/* Rf_errorcall(R_NilValue, "Data Error 8\n"); */
	return 0.0;
      } else {
	int wh, cmt, wh100, whI, wh0;
	getWh(ind->evid[idx+1], &wh, &cmt, &wh100, &whI, &wh0);
	if (whI != 6){
	  if (!(ind->err & 512)){
	    ind->err += 512;
	  }
	  return 0.0;
	  /* Rf_errorcall(R_NilValue, "Data error 886 (whI=%d, evid=%d to %d)\n", whI, */
	  /*       ind->evid[idx], ind->evid[idx+1]); */
	}
	yp = rx->ypNA;
	updateDur(idx, ind, yp);
      }
      break;
    case EVIDF_MODEL_RATE_OFF:
      if (idx > 0){
	int wh, cmt, wh100, whI, wh0;
	getWh(ind->evid[idx-1], &wh, &cmt, &wh100, &whI, &wh0);
	if (whI != 9){
	  if (!(ind->err & 1024)){
	    ind->err += 1024;
	  }
	  /* Rf_errorcall(R_NilValue, "Data error 797 (whI = %d; evid=%d)", whI, ind->evid[idx-1]); */
	  return 0.0;
	}
	yp = rx->ypNA;
	updateRate(idx-1, ind, yp);
      } else {
	if (!(ind->err & 2048)){
	  ind->err += 2048;
	}
	/* Rf_errorcall(R_NilValue, "Data Error -7\n"); */
	return 0.0;
      }
      break;
    case EVIDF_MODEL_RATE_ON:
      // This calculates the rate and the duration and then assigns it to the next record
      if (idx >= ind->n_all_times){
	// error: Last record, can't be used.
	if (!(ind->err & 4096)){
	  ind->err += 4096;
	}
	/* Rf_errorcall(R_NilValue, "Data Error 9\n"); */
	return 0.0;
      } else {
	int wh, cmt, wh100, whI, wh0;
	getWh(ind->evid[idx+1], &wh, &cmt, &wh100, &whI, &wh0);
	if (whI != 7){
	  if (!(ind->err & 8192)){
	    ind->err += 8192;
	  }
	  return 0.0;
	}
	yp = rx->ypNA;
	updateRate(idx, ind, yp);
      }
      break;
    case EVIDF_INF_RATE:
      {
	double amt = ind->dose[idx];
	if (amt > 0){
	  ret = getLag(ind, ind->id, ind->cmt, ind->all_times[idx]);
	  return ret;
	} else if (amt < 0){
	  int j = getDoseNumberFromIndex(ind, idx);
	  if (j == -1){
	    if (!(ind->err & 16384)){
	      ind->err += 16384;
	    }
	    return 0.0;
	    /* Rf_errorcall(R_NilValue, "Corrupted event table during sort (1)."); */
	  }
	  int k;
	  for (k = j; k--;){
	    if (ind->evid[ind->idose[j]] == ind->evid[ind->idose[k]]) break;
	    if (k == 0) {
	      if (!(ind->err & 32768)){
		ind->err += 32768;
	      }
	      return 0.0;
	    }
	  }
	  rx_solve *rx = &rx_global;
	  double f = getAmt(ind, ind->id, ind->cmt, 1.0, ind->all_times[ind->idose[j-1]], rx->ypNA);
	  if (ISNA(f)){
	    rx_solving_options *op = &op_global;
	    op->badSolve=1;
	    op->naTime = 1;
	  }
	  double durOld = (ind->all_times[ind->idose[j]] -
			   ind->all_times[ind->idose[k]]); 
	  double dur = f*durOld;
	  double t = ind->all_times[ind->idose[k]]+dur;
	  ret = getLag(ind, ind->id, ind->cmt, t);
	  return ret;
	} else {
	  /* Rf_errorcall(R_NilValue, "Corrupted events."); */
	  if (!(ind->err & 131072)){
	    ind->err += 131072;
	  }
	  return 0.0;
	}
      }
      break;
    }
  }
  ret = getLag(ind, ind->id, ind->cmt, ind->all_times[idx]);
  return ret;
}


// Adapted from 
// https://github.com/Rdatatable/data.table/blob/588e0725320eacc5d8fc296ee9da4967cee198af/src/forder.c#L630-L649
extern "C" void sortRadix(rx_solving_options_ind *ind){
  rx_solve *rx = &rx_global;
  rx_solving_options *op = &op_global;
  uint8_t **key = rx->keys[0];
  // Reset times for infusion
  int wh, cmt, wh100, whI, wh0;
  int doSort = 1;
  double *time = new double[ind->n_all_times];
  uint64_t *all = new uint64_t[ind->n_all_times];
  uint64_t minD, maxD;
  ind->ixds = 0;
  ind->curShift = 0;
  for (int i = 0; i < ind->n_all_times; i++) {
    ind->ix[i] = i;
    ind->idx = i;
    if (!isObs(ind->evid[i])) {
      getWh(ind->evid[i], &wh, &cmt, &wh100, &whI, &wh0);
      if (whI == 6 || whI == 7) {
	// Reset on every sort (since sorted only once)
	ind->all_times[i] = ind->all_times[i-1];
      }
      time[i] = getTime(ind->ix[i], ind);
      ind->ixds++;
    } else {
      if (ind->evid[i] == 3) {
	ind->curShift -= rx->maxShift;
      }
      time[i] = getTime(ind->ix[i], ind);
    }
    all[i]  = dtwiddle(time, i);
    if (i == 0){
      minD = maxD = all[0];
    } else if (all[i] < minD){
      minD = all[i];
    } else if (all[i] > maxD) {
      maxD = all[i];
    }
    if (op->naTime == 1){
      doSort=0;
      break;
    }
  }
  if (doSort){
    int nradix=0, nbyte=0, spare=0;
    calcNradix(&nbyte, &nradix, &spare, &maxD, &minD);
    rx->nradix[0] = nradix;
    // Allocate more space if needed
    for (int b = 0; b < nbyte; b++){
      if (key[b] == NULL) {
	key[b] = (uint8_t *)calloc(rx->maxAllTimes+1, sizeof(uint8_t));
      }
    }
    for (int i = 0; i < ind->n_all_times; i++) {
      uint64_t elem = all[i] - minD;
      elem <<= spare;
      for (int b= nbyte-1; b>0; b--) {
	key[b][i] = (uint8_t)(elem & 0xff);
	elem >>= 8;
      }
      // RxODE uses key[0][i] = 0 | (uint8_t)(elem & 0xff) instead of
      //  key[0][i] |= (uint8_t)(elem & 0xff)
      // because unlike data.table, key[0][i] is not necessarily zero. 
      key[0][i] = 0 | (uint8_t)(elem & 0xff);
    }
    radix_r(0, ind->n_all_times-1, 0, ind, rx);
  }
  delete[] time;
  delete[] all;
}

