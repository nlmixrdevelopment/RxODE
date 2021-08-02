#ifndef __GETTIME_H___
#define __GETTIME_H___

#include "handle_evid.h"

#if defined(__cplusplus)
#define FLOOR(x) std::floor(x)
extern "C" {
#else
#define FLOOR(x) floor(x)
#endif

  double getTime(int idx, rx_solving_options_ind *ind);


  void radix_r(const int from, const int to, const int radix,
			rx_solving_options_ind *ind, rx_solve *rx);

  void calcNradix(int *nbyte, int *nradix, int *spare, uint64_t *maxD, uint64_t *minD);

  uint64_t dtwiddle(const void *p, int i);

  void sortRadix(rx_solving_options_ind *ind);

  extern t_dydt dydt;

  extern t_calc_jac calc_jac;

  extern t_calc_lhs calc_lhs;

  extern t_update_inis update_inis;

  extern t_dydt_lsoda_dum dydt_lsoda_dum;

  extern t_dydt_liblsoda dydt_liblsoda;

  extern t_jdum_lsoda jdum_lsoda;

  extern t_get_solve get_solve;

  extern t_assignFuns assignFuns;

  extern t_LAG LAG;
  extern t_RATE RATE;
  extern t_DUR DUR;
  extern t_calc_mtime calc_mtime;

  extern t_ME ME;
  extern t_IndF IndF;

  static inline void calcMtime(int solveid, double *mtime){
    calc_mtime(solveid,mtime);
  }

  static inline double getLag(rx_solving_options_ind *ind, int id, int cmt, double time){
    double ret = LAG(id, cmt, time);
    if (ISNA(ret)) {
      rx_solving_options *op = &op_global;
      op->badSolve=1;
      op->naTime = 1;
    }
    return ret;
  }

  static inline double getRate(rx_solving_options_ind *ind, int id, int cmt, double dose, double t){
    double ret = RATE(id, cmt, dose, t);
    if (ISNA(ret)){
      rx_solving_options *op = &op_global;
      op->badSolve=1;
      op->naTime = 1;
    }
    return ret;
  }

  static inline double getDur(rx_solving_options_ind *ind, int id, int cmt, double dose, double t){
    double ret = DUR(id, cmt, dose, t);
    if (ISNA(ret)){
      rx_solving_options *op = &op_global;
      op->badSolve=1;
      op->naTime = 1;
    }
    return ret;
  }

  static inline void updateDur(int idx, rx_solving_options_ind *ind, double *yp){
    double t = ind->all_times[idx];
    int oldIdx = ind->idx;
    ind->idx=idx;
    if (ind->all_times[idx+1] == t){
      double dur, rate, amt;
      // The duration and f cannot depend on state values
      amt  = getAmt(ind, ind->id, ind->cmt, ind->dose[idx], t, yp);
      dur  = getDur(ind, ind->id, ind->cmt, amt, t);
      if (dur > 0){
	rate = amt/dur;// mg/hr
	ind->dose[idx+1]      = -rate;
	ind->all_times[idx+1] = t + dur;
	ind->idx=oldIdx;
      } else {
	rx_solve *rx = &rx_global;
	rx_solving_options *op = &op_global;
	if (ind->cmt < op->neq){
	  if (rx->needSort & 4){
	    if (!(ind->err & 16)){
	      ind->err += 16;
	    }
	    return;
	    /* Rf_errorcall(R_NilValue, "Duration is zero/negative (dur=%f; cmt=%d; amt=%f)", dur, ind->cmt+1, amt); */
	  } else {
	    if (!(ind->err & 32)){
	      ind->err += 32;
	    }
	    return;
	    /* Rf_errorcall(R_NilValue, "Modeled duration requested in event table, but not in model; use 'dur(cmt) ='"); */
	  }
	}
      }
    }
    ind->idx=oldIdx;
  }

  static inline void updateRate(int idx, rx_solving_options_ind *ind, double *yp) {
    double t = ind->all_times[idx];
    int oldIdx = ind->idx;
    ind->idx=idx;
    if (ind->all_times[idx+1] == t){
      double dur, rate, amt;
      amt  = getAmt(ind, ind->id, ind->cmt, ind->dose[idx], t, yp);
      rate  = getRate(ind, ind->id, ind->cmt, amt, t);
      if (rate > 0){
	dur = amt/rate; // mg/hr
	ind->dose[idx+1]      = - rate;
	ind->all_times[idx+1] = t+dur;
	ind->idx=oldIdx;
      } else {
	rx_solve *rx;
	rx = &rx_global;
	rx_solving_options *op = &op_global;
	if (ind->cmt < op->neq){
	  if (rx->needSort & 8){
	    if (!(ind->err & 2)){
	      ind->err += 2;
	      /* Rf_errorcall(R_NilValue, "Rate is zero/negative"); */
	    }
	    return;
	  } else {
	    // FIXME don't error out with linear compartmental model
	    if (!(ind->err & 4)){
	      ind->err += 4;
	    }
	    return;
	    /* Rf_errorcall(R_NilValue, "Modeled rate requested in event table, but not in model; use 'rate(cmt) ='"); */
	  }
	}
	// error rate is zero/negative
      }
    }
    ind->idx=oldIdx;
  }

#if defined(__cplusplus)
}
#endif


#endif
