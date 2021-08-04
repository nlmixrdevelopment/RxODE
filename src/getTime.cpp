#include "getTime.h"

extern "C" double getTime(int idx, rx_solving_options_ind *ind) {
  return getTime__(idx, ind, 0);
}

// Adapted from 
// https://github.com/Rdatatable/data.table/blob/588e0725320eacc5d8fc296ee9da4967cee198af/src/forder.c#L630-L649
extern "C" void sortRadix(rx_solving_options_ind *ind){
  rx_solve *rx = &rx_global;
  rx_solving_options *op = &op_global;
  uint8_t **key = rx->keys[0];
  // Reset times for infusion
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
      time[i] = getTime__(ind->ix[i], ind, 1);
      ind->ixds++;
    } else {
      if (ind->evid[i] == 3) {
	ind->curShift -= rx->maxShift;
      }
      time[i] = getTime__(ind->ix[i], ind, 1);
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
