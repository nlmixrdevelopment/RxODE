#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>
#include "../inst/include/RxODE.h"
#include "dop853.h"
#include "common.h"
#include "lsoda.h"
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
// Yay easy parallel support
// For Mac, see: http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/ (as far as I can tell)
// and https://github.com/Rdatatable/data.table/wiki/Installation#openmp-enabled-compiler-for-mac
// It may have arrived, though I'm not sure...
// According to http://dirk.eddelbuettel.com/papers/rcpp_parallel_talk_jan2015.pdf
// OpenMP is excellent for parallelizing existing loops where the iterations are independent;
// OpenMP is used by part of the R core, therefore support will come for all platforms at some time in the future.
// Since these are independent, we will just use Open MP.
#ifdef _OPENMP
#include <omp.h>
#endif
#define RSprintf(fmt,...) if (_setSilentErr == 0) REprintf(fmt,__VA_ARGS__);
#define RSprintf0(fmt) if (_setSilentErr == 0) REprintf(fmt);
int _setSilentErr=0;
void setSilentErr(int silent){
  _setSilentErr = silent;
}

extern int getSilentErr(){return _setSilentErr;}

void printErr(int err, int id){
  RSprintf("Recovered solving errors for internal ID %d (%d):\n", id+1, err);
  if (err & 1){
    RSprintf0("  Corrupted event table during sort (1)\n");
  }
  if (err & 2){
    RSprintf0("  Rate is zero/negative\n");
  }
  if (err & 4){
    RSprintf0("  Modeled rate requested in event table, but not in model; use 'rate(cmt) ='\n");
  }
  if (err & 4){
    RSprintf0("  Modeled rate requested in event table, but not in model; use 'rate(cmt) ='\n");
  }
  if (err & 8){
    RSprintf0("  Corrupted event table during sort (2)\n");
  }
  if (err & 16){
    RSprintf0("  Duration is zero/negative\n");
  }
  if (err & 32){
    RSprintf0("  Modeled duration requested in event table, but not in model; use 'dur(cmt) ='\n");
  }
  if (err & 64){
    RSprintf0("  Data error 686\n");
  }
  if (err & 128){
    RSprintf0("  Data Error -6\n");
  }
  if (err & 256){
    RSprintf0("  Data Error 8\n");
  }
  if (err & 512){
    RSprintf0("  Data error 886\n");
  }
  if (err & 1024){
    RSprintf0("  Data error 797\n");
  }
  if (err & 2048){
    RSprintf0("  Data Error -7\n");
  }
  if (err & 4096){
    RSprintf0("  Data Error 9\n");
  }
  if (err & 8192){
    RSprintf0("  Data error 997\n");
  }
  if (err & 16384){
    RSprintf0("  Corrupted event table during sort (1)\n");
  }
  if (err & 32768){
    RSprintf0("  Corrupted event table\n");
  }
  if (err & 131072){
    RSprintf0("  Corrupted events\n");
  }
  if (err & 65536){
    RSprintf0("  Supplied an invalid EVID\n");
  }
  if (err & 262144){
    RSprintf0("  Corrupted event table\n");
  }
  if (err & 524288){
    RSprintf0("  The event table has been corrupted\n");
  }
  if (err & 1048576){
    RSprintf0("  SS=2 & Modeled F does not work\n");
  }
  if (err & 2097152){
    RSprintf0("  SS=2 & Modeled F does not work\n");
  }
  if (err & 4194304){
    RSprintf0("  SS=2 & Modeled F does not work\n");
  }
  if (err & 8388608){
    RSprintf0(" Rate is zero/negative\n");
  }
  
}

rx_solve rx_global;

rx_solving_options op_global;

rx_solving_options_ind *inds_global = NULL;
int gitol=0, gitask = 1, giopt = 0, gliw=0, glrw = 0;

void par_flush_console() {
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
  R_FlushConsole();
#endif
}

int isRstudio();
int isProgSupported();
int par_progress_0=0;
int par_progress(int c, int n, int d, int cores, clock_t t0, int stop){
  float progress = (float)(c)/((float)(n));
  if (progress < 0.) progress = 0.;
  if (progress > 1.) progress = 1.;
  if (progress == 0.) par_progress_0=0;
  if (c <= n){
    int nticks= (int)(progress * 50);
    int curTicks = d;
    if (nticks < 0) nticks=0;
    if (nticks > 50) nticks=50;
    if (curTicks < 0) curTicks=0;
    if (curTicks > 50) curTicks=50;
    int isSupported = isProgSupported();
    if (isSupported == -1){
    } else if (isSupported == 0){
      int i;
      for (i = curTicks; i < nticks; i++){
	if (i == 0) {
	  Rprintf("[");
	} else if (i % 5 == 0) {
	  Rprintf("|");
	} else {
	  Rprintf("=");
	}
      }
      if (nticks == 50){
	if (!par_progress_0){
	  par_progress_0 = 1;
	  Rprintf("] ");
	  clock_t t = clock() - t0;
	  double ts = ((double)t)/CLOCKS_PER_SEC;
	  if (ts < 60){
	    Rprintf("0:00:%02.f ", floor(ts));
	  } else {
	    double f = floor(ts/60);
	    double s = ts-f*60;
	    if (f >= 60){
	      double h = floor(f/60);
	      f = f-h*60;
	      Rprintf("%.0f:%02.f:%02.f ", h, f, floor(s));
	    } else {
	      Rprintf("0:%02.f:%02.f ", f, floor(s));
	    }
	  }
	}
      }
    } else if (isRstudio()){
      if (!par_progress_0){
	Rprintf("\r");
	int i;
	for (i = 0; i < nticks; i++){
	  if (i == 0) {
	    Rprintf("[");
	  } else if (i % 5 == 0) {
	    Rprintf("|");
	  } else {
	    Rprintf("=");
	  }
	}
	if (nticks < 50) {Rprintf(">");}
	else {par_progress_0 = 1;}
	for (i = nticks+1; i < 50; i++){
	  Rprintf("-");
	}
	Rprintf("] ");
	if (nticks < 50) Rprintf(" ");
	Rprintf("%02.f%%; ",100*progress,cores);
	clock_t t = clock() - t0;
	double ts = ((double)t)/CLOCKS_PER_SEC;
	if (ts < 60){
	  Rprintf("0:00:%02.f ", floor(ts));
	} else {
	  double f = floor(ts/60);
	  double s = ts-f*60;
	  if (f >= 60){
	    double h = floor(f/60);
	    f = f-h*60;
	    Rprintf("%.0f:%02.f:%02.f ", h, f, floor(s));
	  } else {
	    Rprintf("0:%02.f:%02.f ", f, floor(s));
	  }
	}
	if (stop){
	  Rprintf("Stopped Calculation!\n");
	}
      }
    } else {
      if (!par_progress_0){
	RSprintf0("\r");
	int i;
	for (i = 0; i < nticks; i++){
	  if (i == 0) {
	    RSprintf0("%%[");
	  } else if (i % 5 == 0) {
	    RSprintf0("|");
	  } else {
	    RSprintf0("=");
	  }
	}
	if (nticks < 50) { RSprintf0(">");}
	else {par_progress_0 = 1;}
	if (nticks + 1 < 50){
	  for (i = nticks+1; i < 50; i++){
	    RSprintf0("-");
	  }
	}
	RSprintf0("] ");
	if (nticks < 50) RSprintf0(" ");
	RSprintf("%02.f%%; ",100*progress,cores);
	clock_t t = clock() - t0;
	double ts = ((double)t)/CLOCKS_PER_SEC;
	if (ts < 60){
	  RSprintf("0:00:%02.f ", floor(ts));
	} else {
	  double f = floor(ts/60);
	  double s = ts-f*60;
	  if (f >= 60){
	    double h = floor(f/60);
	    f = f-h*60;
	    RSprintf("%.0f:%02.f:%02.f ", h, f, floor(s));
	  } else {
	    RSprintf("0:%02.f:%02.f", f, floor(s));
	  }
	}
	if (stop){
	  RSprintf0("Stopped Calculation!\n");
	}
      }
    }
    par_flush_console();
    return nticks;
  }
  return d;
}

typedef struct {
  int cur;
  int n;
  int d;
  int cores;
  clock_t t0;
} rx_tick;

rx_tick rxt;

SEXP _rxTick(){
  rxt.cur++;
  SEXP ret = PROTECT(allocVector(INTSXP, 1));
  rxt.d =par_progress(rxt.cur, rxt.n, rxt.d, rxt.cores, rxt.t0, 0);
  INTEGER(ret)[0]=rxt.d;
  UNPROTECT(1);
  return ret;
}

SEXP _rxProgress(SEXP num, SEXP core){
  rxt.t0 = clock();
  rxt.cores = INTEGER(core)[0];
  rxt.n = INTEGER(num)[0];
  rxt.d=0;
  rxt.cur = 0;
  return R_NilValue;
}

SEXP _rxProgressStop(SEXP clear){
  par_progress_0=0;
  int clearB = INTEGER(clear)[0];
  if (clearB){
    int doIt=isProgSupported();
    if (doIt == -1){
    } else if (isRstudio() || doIt==0){
      Rprintf("\n");
    } else {
      RSprintf0("\r                                                                                 \r");
    }
  } else {
    par_progress(rxt.n, rxt.n, rxt.d, rxt.cores, rxt.t0, 1);
    int doIt=isProgSupported();
    if (isRstudio() || doIt == 0){
      Rprintf("\n");
    }
  }
  rxt.d = rxt.n;
  rxt.cur = rxt.n;
  return R_NilValue;
}

SEXP _rxProgressAbort(SEXP str){
  par_progress_0=0;
  if (rxt.d != rxt.n || rxt.cur != rxt.n){
    par_progress(rxt.n, rxt.n, rxt.d, rxt.cores, rxt.t0, 0);
    error(CHAR(STRING_ELT(str,0)));
  }
  return R_NilValue;
}

t_set_solve set_solve = NULL;

void rxOptionsIniEnsure(int mx){
  Free(inds_global);
  inds_global = Calloc(mx, rx_solving_options_ind);
  rx_solve *rx=(&rx_global);
  rx->subjects = inds_global;  
}

t_dydt dydt = NULL;

t_calc_jac calc_jac = NULL;

t_calc_lhs calc_lhs = NULL;

t_update_inis update_inis = NULL;

t_dydt_lsoda_dum dydt_lsoda_dum = NULL;

t_dydt_liblsoda dydt_liblsoda = NULL;

t_jdum_lsoda jdum_lsoda = NULL;

t_get_solve get_solve = NULL;

t_assignFuns assignFuns=NULL;

t_get_theta get_theta = NULL;

t_F AMT = NULL;
t_LAG LAG = NULL;
t_RATE RATE = NULL;
t_DUR DUR = NULL;
t_calc_mtime calc_mtime = NULL;

int global_jt = 2;
int global_mf = 22;  
int global_debug = 0;

double *global_rworkp;
int *global_iworkp;

void rxUpdateFuns(SEXP trans){
  const char *lib, *s_dydt, *s_calc_jac, *s_calc_lhs, *s_inis, *s_dydt_lsoda_dum, *s_dydt_jdum_lsoda, 
    *s_ode_solver_solvedata, *s_ode_solver_get_solvedata, *s_dydt_liblsoda, *s_AMT, *s_LAG, *s_RATE,
    *s_DUR, *s_mtime, *s_theta, *s_assignFuns;
  lib = CHAR(STRING_ELT(trans, 0));
  s_dydt = CHAR(STRING_ELT(trans, 3));
  s_calc_jac = CHAR(STRING_ELT(trans, 4));
  s_calc_lhs = CHAR(STRING_ELT(trans, 5));
  s_inis = CHAR(STRING_ELT(trans, 8));
  s_dydt_lsoda_dum = CHAR(STRING_ELT(trans, 9));
  s_dydt_jdum_lsoda = CHAR(STRING_ELT(trans, 10));
  s_ode_solver_solvedata = CHAR(STRING_ELT(trans, 11));
  s_ode_solver_get_solvedata = CHAR(STRING_ELT(trans, 12));
  s_dydt_liblsoda = CHAR(STRING_ELT(trans, 13));
  s_AMT=CHAR(STRING_ELT(trans,14));
  s_LAG=CHAR(STRING_ELT(trans, 15));
  s_RATE=CHAR(STRING_ELT(trans, 16));
  s_DUR=CHAR(STRING_ELT(trans, 17));
  s_mtime=CHAR(STRING_ELT(trans, 18));
  s_assignFuns=CHAR(STRING_ELT(trans, 19));
  s_theta=CHAR(STRING_ELT(trans, 7));
  global_jt = 2;
  global_mf = 22;  
  global_debug = 0;
  if (strcmp(CHAR(STRING_ELT(trans, 1)),"fulluser") == 0){
    global_jt = 1;
    global_mf = 21;
  } else {
    global_jt = 2;
    global_mf = 22;
  }
  calc_lhs =(t_calc_lhs) R_GetCCallable(lib, s_calc_lhs);
  dydt =(t_dydt) R_GetCCallable(lib, s_dydt);
  calc_jac =(t_calc_jac) R_GetCCallable(lib, s_calc_jac);
  update_inis =(t_update_inis) R_GetCCallable(lib, s_inis);
  dydt_lsoda_dum =(t_dydt_lsoda_dum) R_GetCCallable(lib, s_dydt_lsoda_dum);
  jdum_lsoda =(t_jdum_lsoda) R_GetCCallable(lib, s_dydt_jdum_lsoda);
  set_solve = (t_set_solve)R_GetCCallable(lib, s_ode_solver_solvedata);
  get_solve = (t_get_solve)R_GetCCallable(lib, s_ode_solver_get_solvedata);
  dydt_liblsoda = (t_dydt_liblsoda)R_GetCCallable(lib, s_dydt_liblsoda);
  AMT = (t_F)R_GetCCallable(lib, s_AMT);
  LAG = (t_LAG) R_GetCCallable(lib, s_LAG);
  RATE = (t_RATE) R_GetCCallable(lib, s_RATE);
  DUR = (t_DUR) R_GetCCallable(lib, s_DUR);
  calc_mtime = (t_calc_mtime) R_GetCCallable(lib, s_mtime);
  get_theta = (t_get_theta) R_GetCCallable(lib, s_theta);
  assignFuns = R_GetCCallable(lib, s_assignFuns);
}

void rxClearFuns(){
  calc_lhs		= NULL;
  dydt			= NULL;
  calc_jac		= NULL;
  update_inis		= NULL;
  dydt_lsoda_dum	= NULL;
  jdum_lsoda		= NULL;
  set_solve		= NULL;
  get_solve		= NULL;
  dydt_liblsoda		= NULL;
}

void F77_NAME(dlsoda)(
                      void (*)(int *, double *, double *, double *),
                      int *, double *, double *, double *, int *, double *, double *,
                      int *, int *, int *, double *,int *,int *, int *,
                      void (*)(int *, double *, double *, int *, int *, double *, int *),
                      int *);

extern rx_solve *getRxSolve2_(){
  return &rx_global;
}
extern rx_solve *getRxSolve_(){
  /* if (set_solve == NULL) */
  /*   error("RxODE model function pointers are not setup."); */
  /* set_solve(&rx_global); */
  rx_solve *rx;
  rx = &rx_global;
  rx->subjects = inds_global;
  rx->op = &op_global;
  return &rx_global;
}

rx_solving_options_ind *getRxId(rx_solve *rx, unsigned int id){
  return &(rx->subjects[id]);
}

void doSort(rx_solving_options_ind *ind);

void getWh(int evid, int *wh, int *cmt, int *wh100, int *whI, int *wh0){
  *wh = evid;
  *cmt = 0;
  *wh100 = floor(*wh/1e5L);
  *whI   = floor(*wh/1e4L-*wh100*10);
  *wh    = *wh - *wh100*1e5 - (*whI-1)*1e4;
  *wh0 = floor((*wh%10000)/100);
  *cmt = *wh0 - 1 + *wh100*100;
  *wh0 = evid - *wh100*1e5 - *whI*1e4 - *wh0*100;
}

void updateRate(int idx, rx_solving_options_ind *ind){
  double t = ind->all_times[idx];
  int oldIdx = ind->idx;
  ind->idx=idx;
  if (ind->all_times[idx+1] == t){
    // Hasn't been calculated yet.
    int j;
    // Find the amount
    // bisection https://en.wikipedia.org/wiki/Binary_search_algorithm
    int l = 0, r = ind->ndoses-1, m=0;
    while(l <= r){
      m = floor((l+r)/2);
      if (ind->idose[m] < idx) l = m+1;
      else if (ind->idose[m] > idx) r = m-1;
      else break;
    }
    if (ind->idose[m] == idx){
      j=m;
    } else {
      if (!(ind->err & 1)){
	ind->err += 1;
      }
      return;
      /* error("Corrupted event table during sort (1)."); */
    }
    double dur, rate, amt;
    amt  = AMT(ind->id, ind->cmt, ind->dose[j], t);
    rate  = RATE(ind->id, ind->cmt, amt, t);
    if (rate > 0){
      dur = amt/rate;// mg/hr
      ind->dose[j+1] = -rate;
      ind->all_times[idx+1]=t+dur;
      ind->idx=oldIdx;
    } else {
      rx_solve *rx;
      rx = &rx_global;
      rx_solving_options *op = &op_global;
      if (ind->cmt < op->neq){
	if (rx->needSort & 8){
	  if (!(ind->err & 2)){
	    ind->err += 2;
	    /* error("Rate is zero/negative"); */
	  }
	  return;
	} else {
	  // FIXME don't error out with linear compartmental model
	  if (!(ind->err & 4)){
	    ind->err += 4;
	  }
	  return;
	  /* error("Modeled rate requested in event table, but not in model; use 'rate(cmt) ='"); */
	}
      }
      // error rate is zero/negative
    }
  }
}

void updateDur(int idx, rx_solving_options_ind *ind){
  double t = ind->all_times[idx];
  int oldIdx = ind->idx;
  ind->idx=idx;
  if (ind->all_times[idx+1] == t){
    // Hasn't been calculated yet.
    int j;
    // Find the amount
    // Find the amount
    // bisection https://en.wikipedia.org/wiki/Binary_search_algorithm
    int l = 0, r = ind->ndoses-1, m=0;
    while(l <= r){
      m = floor((l+r)/2);
      if (ind->idose[m] < idx) l = m+1;
      else if (ind->idose[m] > idx) r = m-1;
      else break;
    }
    if (ind->idose[m] == idx){
      j=m;
    } else {
      if (!(ind->err & 8)){
	ind->err += 8;
      }
      return;
      /* error("Corrupted event table during sort (2)."); */
    }
    double dur, rate, amt;
    amt  = AMT(ind->id, ind->cmt, ind->dose[j], t);
    dur  = DUR(ind->id, ind->cmt, amt, t);
    if (dur > 0){
      rate = amt/dur;// mg/hr
      ind->dose[j+1] = -rate;
      ind->all_times[idx+1]=t+dur;
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
	  /* error("Duration is zero/negative (dur=%f; cmt=%d; amt=%f)", dur, ind->cmt+1, amt); */
	} else {
	  if (!(ind->err & 32)){
	    ind->err += 32;
	  }
	  return;
	  /* error("Modeled duration requested in event table, but not in model; use 'dur(cmt) ='"); */
	}
      }
    }
  }
}

extern double getTime(int idx, rx_solving_options_ind *ind){
  int evid = ind->evid[idx];
  if (evid == 9) return 0;
  if (evid >= 10 && evid <= 99) return ind->mtime[evid-10];
  if (isObs(evid)) return ind->all_times[idx];
  getWh(evid, &(ind->wh), &(ind->cmt), &(ind->wh100), &(ind->whI), &(ind->wh0));
  if (ind->wh0 == 40){
  } else {
    switch(ind->whI){
    case 6:
      if (idx > 0){
	int wh, cmt, wh100, whI, wh0;
	getWh(ind->evid[idx-1], &wh, &cmt, &wh100, &whI, &wh0);
	if (whI != 8){
	  // FIXME can crash parallel runs and cause many issues.  Need to defer to end.
	  if (!(ind->err & 64)){
	    ind->err += 64;
	  }
	  return 0.0;
	  /* error("Data error 686 (whI = %d; evid=%d)", whI, ind->evid[idx-1]); */
	}
	updateDur(idx-1, ind);
      } else {
	if (!(ind->err & 128)){
	  ind->err += 128;
	}
	return 0.0;
	/* error("Data Error -6\n"); */
      }
      break;
    case 8:
      if (idx >= ind->n_all_times){
	// error: Last record, can't be used.
	if (!(ind->err & 256)){
	  ind->err += 256;
	}
	/* error("Data Error 8\n"); */
	return 0.0;
      } else {
	int wh, cmt, wh100, whI, wh0;
	getWh(ind->evid[idx+1], &wh, &cmt, &wh100, &whI, &wh0);
	if (whI != 6){
	  if (!(ind->err & 512)){
	    ind->err += 512;
	  }
	  return 0.0;
	  /* error("Data error 886 (whI=%d, evid=%d to %d)\n", whI, */
	  /*       ind->evid[idx], ind->evid[idx+1]); */
	}
	updateDur(idx, ind);
      }
      break;
    case 7:
      if (idx > 0){
	int wh, cmt, wh100, whI, wh0;
	getWh(ind->evid[idx-1], &wh, &cmt, &wh100, &whI, &wh0);
	if (whI != 9){
	  if (!(ind->err & 1024)){
	    ind->err += 1024;
	  }
	  /* error("Data error 797 (whI = %d; evid=%d)", whI, ind->evid[idx-1]); */
	  return 0.0;
	}
	updateRate(idx-1, ind);
      } else {
	if (!(ind->err & 2048)){
	  ind->err += 2048;
	}
	/* error("Data Error -7\n"); */
	return 0;
      }
      break;
    case 9:
      // This calculates the rate and the duration and then assigns it to the next record
      if (idx >= ind->n_all_times){
	// error: Last record, can't be used.
	if (!(ind->err & 4096)){
	  ind->err += 4096;
	}
	/* error("Data Error 9\n"); */
	return 0.0;
      } else {
	int wh, cmt, wh100, whI, wh0;
	getWh(ind->evid[idx+1], &wh, &cmt, &wh100, &whI, &wh0);
	if (whI != 7){
	  if (!(ind->err & 8192)){
	    ind->err += 8192;
	  }
	  return 0.0;
	  /* error("Data error 997 (whI=%d, evid=%d to %d)\n", whI, */
	  /*       ind->evid[idx], ind->evid[idx+1]); */
	}
	updateRate(idx, ind);
      }
      break;
    case 1:
      {
	int j;
	// Find the amount
	// bisection https://en.wikipedia.org/wiki/Binary_search_algorithm
	int l = 0, r = ind->ndoses-1, m=0;
	while(l <= r){
	  m = floor((l+r)/2);
	  if (ind->idose[m] < idx) l = m+1;
	  else if (ind->idose[m] > idx) r = m-1;
	  else break;
	}
	if (ind->idose[m] == idx){
	  j=m;
	} else {
	  if (!(ind->err & 16384)){
	    ind->err += 16384;
	  }
	  return 0.0;
	  /* error("Corrupted event table during sort (1)."); */
	}
	if (ind->dose[j] > 0){
	  return LAG(ind->id, ind->cmt, ind->all_times[idx]);
	} else if (ind->dose[j] < 0){
	  // f*amt/rate=dur
	  // amt/rate=durOld
	  // f = dur/durOld
	  // f*durOld = dur
	  int k;
	  for (k = j; k--;){
	    if (ind->evid[ind->idose[j]] == ind->evid[ind->idose[k]]) break;
	    if (k == 0) {
	      if (!(ind->err & 32768)){
		ind->err += 32768;
	      }
	      return 0.0;
	      /* error("corrupted event table"); */
	    }
	  }
	  double f = AMT(ind->id, ind->cmt, 1.0, ind->all_times[ind->idose[j-1]]);
	  double durOld = (ind->all_times[ind->idose[j]] - ind->all_times[ind->idose[k]]); 
	  double dur = f*durOld;
	  double t = ind->all_times[ind->idose[k]]+dur;
	  return LAG(ind->id, ind->cmt, t);
	} else {
	  /* error("Corrupted events."); */
	  if (!(ind->err & 131072)){
	    ind->err += 131072;
	  }
	  return 0.0;
	}
      }
      break;
    }
  }
  return LAG(ind->id, ind->cmt, ind->all_times[idx]);
}

inline static int syncIdx(rx_solving_options_ind *ind){
  if (ind->ix[ind->idx] != ind->idose[ind->ixds]){
    // bisection https://en.wikipedia.org/wiki/Binary_search_algorithm
    int l = 0, r = ind->ndoses-1, m=0;
    while(l <= r){
      m = floor((l+r)/2);
      if (ind->idose[m] < ind->ix[ind->idx]) l = m+1;
      else if (ind->idose[m] > ind->ix[ind->idx]) r = m-1;
      else break;
    }
    if (ind->idose[m] == ind->ix[ind->idx]){
      ind->ixds=m;
    } else {
      //262144
      if (!(ind->err & 262144)){
	ind->err += 262144;
      }
      return 0;
      /* error("Corrupted event table; EVID=%d: %d %d %d", evid, ind->idose[m], ind->ix[ind->idx], */
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
      /* error("The event table has been corrupted; ind->idx: %d ind->ixds: %d ind->idose: %d.", */
      /* 	ind->ix[ind->idx], ind->ixds, ind->idose[ind->ixds]); */
    }
  }
  return 1;
}

int handle_evid(int evid, int neq, 
		int *BadDose,
		double *InfusionRate,
		double *dose,
		double *yp,
		int do_transit_abs,
		double xout, int id,
		rx_solving_options_ind *ind){
  if (isObs(evid)) return 0;
  int wh = evid, cmt, foundBad, j;
  double tmp;
  if (wh) {
    getWh(evid, &(ind->wh), &(ind->cmt), &(ind->wh100), &(ind->whI), &(ind->wh0));
    if (ind->wh0 == 40){
      ind->ixds++;
      return 1;
    }
    /* wh100 = ind->wh100; */
    wh = ind->wh;
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
      if (syncIdx(ind) == 0) return 0;
      if (ind->wh0 == 30){
	yp[cmt]=op_global.inits[cmt];
	InfusionRate[cmt] = 0;
	ind->on[cmt] = 0;
	return 1;
      }
      if (!ind->doSS && ind->wh0 == 20){
	// Save for adding at the end
	memcpy(ind->solveSave, yp, neq*sizeof(double));
      }
      switch(ind->whI){
      case 9: // modeled rate.
      case 8: // modeled duration.
	// Rate already calculated and saved in the next dose record
	ind->on[cmt] = 1;
	InfusionRate[cmt] -= dose[ind->ixds+1];
	if (ind->wh0 == 20 && AMT(id, cmt, dose[ind->ixds], xout) != dose[ind->ixds]){
	  if (!(ind->err & 1048576)){
	    ind->err += 1048576;
	  }
	  return 0;
	  /* error("SS=2 & Modeled F does not work"); */
	}
	break;
      case 7: // End modeled rate
      case 6: // end modeled duration
	// If cmt is off, don't remove rate....
	// Probably should throw an error if the infusion rate is on still.
	InfusionRate[cmt] += dose[ind->ixds]*((double)(ind->on[cmt]));
	if (ind->wh0 == 20 && AMT(id, cmt, dose[ind->ixds], xout) != dose[ind->ixds]){
	  /* error("SS=2 & Modeled F does not work"); */
	  if (!(ind->err & 2097152)){
	    ind->err += 2097152;
	  }
	  return 0;
	}
	break;
      case 2:
	// In this case bio-availability changes the rate, but the duration remains constant.
	// rate = amt/dur
	ind->on[cmt] = 1;
	tmp = AMT(id, cmt, dose[ind->ixds], xout);
	InfusionRate[cmt] += tmp;
	if (ind->wh0 == 20 && tmp != dose[ind->ixds]){
	  /* error("SS=2 & Modeled F does not work"); */
	  if (!(ind->err & 4194304)){
	    ind->err += 4194304;
	  }
	  return 0;
	}
	break;
      case 1:
	ind->on[cmt] = 1;
	InfusionRate[cmt] += dose[ind->ixds];
	if (ind->wh0 == 20 && dose[ind->ixds] > 0 && AMT(id, cmt, dose[ind->ixds], xout) != dose[ind->ixds]){
	  /* error("SS=2 & Modeled F does not work"); */
	  if (!(ind->err & 4194304)){
	    ind->err += 4194304;
	  }
	}
	break;
      case 4: // replace
	ind->on[cmt] = 1;
	ind->podo = 0;
	ind->tlast = xout;
	yp[cmt] = AMT(id, cmt, dose[ind->ixds], xout);     //dosing before obs
	break;
      case 5: //multiply
	ind->on[cmt] = 1;
	ind->podo = 0;
	ind->tlast = xout;
	yp[cmt] *= AMT(id, cmt, dose[ind->ixds], xout);     //dosing before obs
	break;
      case 0:
	if (do_transit_abs) {
	  ind->on[cmt] = 1;
	  if (ind->wh0 == 20){
	    tmp = AMT(id, cmt, dose[ind->ixds], xout);
	    ind->podo = tmp;
	  } else {
	    ind->podo = AMT(id, cmt, dose[ind->ixds], xout);
	  }
	  ind->tlast = xout;
	} else {
	  ind->on[cmt] = 1;
	  ind->podo = 0;
	  ind->tlast = xout;
	  yp[cmt] += AMT(id, cmt, dose[ind->ixds], xout);     //dosing before obs
	}
      }
      /* istate = 1; */
      ind->ixds++;
      /* xp = xout; */
      return 1;
    }
  }
  return 0;
}

static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

int checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

static char *err_msg_ls[]=
    {
      "excess work done on this call (perhaps wrong jt).",
      "excess accuracy requested (tolerances too small).",
      "illegal input detected (see printed message).",
      "repeated error test failures (check all inputs).",
      "repeated convergence failures (perhaps bad jacobian supplied or wrong choice of jt or tolerances).",
      "error weight became zero during problem. (solution component i vanished, and atol or atol(i) = 0.)",
      "work space insufficient to finish (see messages)."
    };

//dummy solout fn
void solout(long int nr, double t_old, double t, double *y, int *nptr, int *irtrn){}

void solveSS_1(int *neq, 
	       int *BadDose,
	       double *InfusionRate,
	       double *dose,
	       double *yp,
	       int do_transit_abs,
	       double xout, double xp, int id,
	       int *i, int nx,
	       int *istate,
	       rx_solving_options *op,
	       rx_solving_options_ind *ind,
	       t_update_inis u_inis,
	       void *ctx){
  int idid;
  switch(op->stiff){
  case 2:
    lsoda(ctx, yp, &xp, xout);
    if (*istate <= 0) {
      RSprintf("IDID=%d, %s\n", *istate, err_msg_ls[-(*istate)-1]);
      /* ind->rc[0] = *istate; */
      ind->rc[0] = -2019;
      // Bad Solve => NA
      /* for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL; */
      /* op->badSolve = 1; */
      /* *i = ind->n_all_times-1; */ // Get out of here!
      /* j=op->maxSS; */
      break;
    } else if (ind->err){
      printErr(ind->err, ind->id);
      ind->rc[0] = -2019;
      /* for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL; */
      /* op->badSolve = 1; */
      *i = ind->n_all_times-1; // Get out of here!
      /* j=op->maxSS; */
      break;
    }
    break;
  case 1:
    F77_CALL(dlsoda)(dydt_lsoda_dum, neq, yp, &xp, &xout,
		     &gitol, &(op->RTOL), &(op->ATOL), &gitask,
		     istate, &giopt, global_rworkp,
		     &glrw, global_iworkp, &gliw, jdum_lsoda, &global_jt);
    if (*istate <= 0) {
      RSprintf("IDID=%d, %s\n", *istate, err_msg_ls[-(*istate)-1]);
      ind->rc[0] = -2019;/* *istate; */
      // Bad Solve => NA
      /* for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL; */
      /* op->badSolve = 1; */
      /* *i = ind->n_all_times-1; // Get out of here! */
      /* j=op->maxSS; */
      break;
    } else if (ind->err){
      printErr(ind->err, ind->id);
      ind->rc[0] = -2019;
      /* for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL; */
      /* op->badSolve = 1; */
      /* *i = ind->n_all_times-1; // Get out of here! */
      /* j=op->maxSS; */
      break;
    }
    break;
  case 0:
    idid = dop853(neq,       /* dimension of the system <= UINT_MAX-1*/
		  dydt,       /* function computing the value of f(x,y) */
		  xp,           /* initial x-value */
		  yp,           /* initial values for y */
		  xout,         /* final x-value (xend-x may be positive or negative) */
		  &(op->RTOL),          /* relative error tolerance */
		  &(op->ATOL),          /* absolute error tolerance */
		  gitol,         /* switch for rtoler and atoler */
		  solout,         /* function providing the numerical solution during integration */
		  0,         /* switch for calling solout */
		  NULL,           /* messages stream */
		  DBL_EPSILON,    /* rounding unit */
		  0,              /* safety factor */
		  0,              /* parameters for step size selection */
		  0,
		  0,              /* for stabilized step size control */
		  0,              /* maximal step size */
		  0,            /* initial step size */
		  op->mxstep,            /* maximal number of allowed steps */
		  1,            /* switch for the choice of the coefficients */
		  -1,                     /* test for stiffness */
		  0,                      /* number of components for which dense outpout is required */
		  NULL,           /* indexes of components for which dense output is required, >= nrdens */
		  0                       /* declared length of icon */
		  );
    if (idid < 0) {
      ind->rc[0] = -2019;
      // Bad Solve => NA
      /* for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL; */
      /* op->badSolve = 1; */
      /* *i = ind->n_all_times-1; // Get out of here! */
      /* j=op->maxSS; */
      break;
    } else if (ind->err){
      printErr(ind->err, ind->id);
      /* ind->rc[0] = -1000; */
      /* for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL; */
      /* op->badSolve = 1; */
      *i = ind->n_all_times-1; // Get out of here!
      /* j=op->maxSS; */
      break;
    }
    break;
  }
}

void handleSS(int *neq,
	      int *BadDose,
	      double *InfusionRate,
	      double *dose,
	      double *yp,
	      int do_transit_abs,
	      double xout, double xp, int id,
	      int *i, int nx,
	      int *istate,
	      rx_solving_options *op,
	      rx_solving_options_ind *ind,
	      t_update_inis u_inis,
	      void *ctx){
  rx_solve *rx = &rx_global;
  int j;
  int doSS2=0;
  int doSSinf=0;
  /* Rprintf("evid: %d\n", ind->evid[ind->ixds-1]); */
  if (((ind->wh0 == 20 || ind->wh0 == 10) &&
      ind->ii[ind->ixds-1] > 0) || ind->wh0 == 40){
    ind->doSS=1;
    ind->ixds--; // This dose stays in place; Reverse dose
    if (ind->wh0 == 20){
      doSS2=1;
    } else if (ind->wh0 == 40){
      doSSinf=1;
    }
    double dur = 0, dur2=0;
    int infBixds =0, infEixds = 0, ei=0, wh, cmt, wh100, whI, wh0, oldI;
    if (doSSinf){
    } else if (ind->whI == 1 || ind->whI == 2){
      oldI = ind->whI;
      infBixds = ind->ixds;
      // Find the next fixed length infusion that is turned off.
      for (j = ind->ixds+1; j < ind->ndoses; j++){
	if (ind->dose[j] == -ind->dose[ind->ixds]){
	  getWh(ind->evid[ind->idose[j]], &wh, &cmt, &wh100, &whI, &wh0);
	  if (whI == oldI && cmt == ind->cmt){
	    dur = getTime(ind->idose[j], ind) - getTime(ind->ix[*i], ind);
	    dur2 = ind->ii[ind->ixds] - dur;
	    /* Rprintf("000; dur: %f; dur2: %f; ii: %f;\n", dur, dur2, ind->ii[ind->ixds]); */
	    infEixds = j;
	    break;
	  }
	}
      }
    } else if (ind->whI == 9 || ind->whI == 8) {
      // These are right next to another.
      infBixds = ind->ixds;
      infEixds = ind->ixds+1;
      dur = getTime(ind->idose[infEixds], ind) - getTime(ind->idose[infBixds],ind);
      dur2 = ind->ii[ind->ixds] - dur;
    }
    /* bi = *i; */
    if (ind->wh0 == 40){
    } else if (ind->whI == 1 || ind->whI == 2 || ind->whI == 8 || ind->whI == 9){
      ei = *i;
      while(ind->ix[ei] != ind->idose[infEixds] && ei < ind->n_all_times){
	ei++;
      }
      if (ind->ix[ei] != ind->idose[infEixds]){
	/* error("Cannot figure out infusion end time."); */
	if (!(ind->err & 8388608)){
	  ind->err += 8388608;
	  /* error("Rate is zero/negative"); */
	}
	return;
      }
    }
    // First Reset
    for (j = neq[0]; j--;) {
      ind->InfusionRate[j] = 0;
      ind->on[j] = 1;
    }
    memcpy(yp,op->inits, neq[0]*sizeof(double));
    u_inis(neq[1], yp); // Update initial conditions @ current time
    if (rx->istateReset) *istate = 1;
    int k;
    double xp2, xout2;
    int canBreak=0;
    xp2 = xp;
    if (doSSinf){
      double rate;
      ind->idx=*i;
      // Rate is fixed, so modifying bio-availability doesn't change duration.
      if (ind->whI == 9){
	rate  = RATE(ind->id, ind->cmt, 0.0, ind->all_times[ind->idose[ind->ixds]]);
      } else {
	rate = ind->dose[ind->ixds];
      }
      ind->InfusionRate[ind->cmt] = rate;
      ind->on[ind->cmt] = 1;
      double infStep = op->infSSstep, a1=1.0, t1=xp2+1.0;
      // Based on http://www.rxkinetics.com/theo.html -- Chiou method
      for (j = 0; j < op->maxSS; j++){
	if (j == 0) xout2 = xp2+1.; // the first level drawn one hour after infusion
	else xout2 = xp2+infStep; 
	solveSS_1(neq, BadDose, InfusionRate, dose, yp, op->do_transit_abs,
		  xout2, xp2, id, i, nx, istate, op, ind, u_inis, ctx);
	canBreak=1;
	if (j <= op->minSS -1){
	  for (k = neq[0]; k--;) {
	    ind->solveLast[k] = yp[k];
	  }
	  if (j == 0) {
	    a1 = yp[ind->cmt];
	  }
	  canBreak=0;
	} else {
	  for (k = neq[0]; k--;){
	    if (op->ssRtol[k]*fabs(yp[k]) + op->ssAtol[k] <= fabs(yp[k]-ind->solveLast[k])){
	      canBreak=0;
	    }
	    ind->solveLast[k] = yp[k];
	  }
	  if (canBreak){
	    ind->InfusionRate[ind->cmt] = 0.0;
	    break;
	  } else {
	    // Assumes that this is at least one half life.
	    double a2 = yp[ind->cmt];	  
	    infStep = max2(infStep,M_LN2/(rate/(a1+a2) + 2*(a1-a2)/((a1+a2)*(xout-t1))));
	  }
	}
	xp2=xout;
	*istate=1;
      }
    } else if (dur == 0){
      // Oral or Steady State Infusion
      for (j = 0; j < op->maxSS; j++){
	xout2 = xp2+ind->ii[ind->ixds];
	// Use "real" xout for handle_evid functions.
	ind->idx=*i;
	handle_evid(ind->evid[ind->ix[*i]], neq[0], BadDose, InfusionRate, dose, yp,
		    op->do_transit_abs, xout, neq[1], ind);
	// yp is last solve or y0
	solveSS_1(neq, BadDose, InfusionRate, dose, yp, op->do_transit_abs,
		  xout2, xp2, id, i, nx, istate, op, ind, u_inis, ctx);
	ind->ixds--; // This dose stays in place
	canBreak=1;
	if (j <= op->minSS -1){
	  if (ind->rc[0]== -2019){
	    for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	    op->badSolve = 1;
	    *i = ind->n_all_times-1;
	    break;
	  }
 	  for (k = op->neq; k--;) {
	    ind->solveLast[k] = yp[k];
	  }
	  canBreak=0;
	} else if (j >= op->minSS){
	  if (ind->rc[0] == -2019){
	    for (k = neq[0]; k--;) {
	      yp[k] = ind->solveLast[k];
	    }
	    ind->rc[0] = 2019;
	    break;
	  }
	  for (k = neq[0]; k--;){
	    if (op->ssRtol[k]*fabs(yp[k]) + op->ssAtol[k] <= fabs(yp[k]-ind->solveLast[k])){
	      canBreak=0;
	    }
	    ind->solveLast[k] = yp[k];
	  }
	  if (canBreak){
	    break;
	  }
	}
	*istate=1;
	xp2 = xout2;
      }
    } else {
      if (dur >= ind->ii[ind->ixds]){
	ind->wrongSSDur=1;
	// Bad Solve => NA
	for (j = neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	op->badSolve = 1;
	*i = nx-1; // Get out of here!
      } else if (ind->err){
	printErr(ind->err, ind->id);
	for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	op->badSolve = 1;
	*i = nx-1; // Get out of here!
      } else {
	// Infusion
	for (j = 0; j < op->maxSS; j++){
	  // Turn on Infusion, solve (0-dur)
	  canBreak=1;
	  xout2 = xp2+dur;
	  ind->idx=*i;
	  ind->ixds = infBixds;
	  handle_evid(ind->evid[ind->idose[infBixds]], neq[0], BadDose, InfusionRate, dose, yp,
		      op->do_transit_abs, xout, neq[1], ind);
	  // yp is last solve or y0
	  *istate=1;
	  // yp is last solve or y0
	  solveSS_1(neq, BadDose, InfusionRate, dose, yp, op->do_transit_abs,
		    xout2, xp2, id, i, nx, istate, op, ind, u_inis, ctx);
	  xp2 = xout2;
	  // Turn off Infusion, solve (dur-ii)
	  xout2 = xp2+dur2;
	  ind->ixds = infEixds;
	  ind->idx=ei;
	  handle_evid(ind->evid[ind->idose[infEixds]], neq[0], BadDose, InfusionRate, dose, yp,
		      op->do_transit_abs, xout+dur, neq[1], ind);
	  if (j <= op->minSS -1){
	    if (ind->rc[0]== -2019){
	      for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	      op->badSolve = 1;
	      *i = ind->n_all_times-1;
	      break;
	    }
	    for (k = neq[0]; k--;) {
	      ind->solveLast[k] = yp[k];
	    }
	    canBreak=0;
	  } else if (j >= op->minSS){
	    if (ind->rc[0]== -2019){
	      if (op->strictSS){
                for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
                op->badSolve = 1;
                *i = ind->n_all_times-1;
              } else {
                for (k = neq[0]; k--;){
                  yp[k] = ind->solveLast[k];
                }
                ind->rc[0] = 2019;
              }
	    }
	    for (k = neq[0]; k--;) {
	      ind->solveLast[k] = yp[k];
	      if (op->ssRtol[k]*fabs(yp[k]) + op->ssAtol[k] <= fabs(yp[k]-ind->solveLast[k])){
		canBreak=0;
	      }
	    }
	  }
	  // yp is last solve or y0
	  *istate=1;
	  solveSS_1(neq, BadDose, InfusionRate, dose, yp, op->do_transit_abs,
		    xout2, xp2, id, i, nx, istate, op, ind, u_inis, ctx);
	  if (j <= op->minSS -1){
	    if (ind->rc[0]== -2019){
	      for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	      op->badSolve = 1;
	      *i = ind->n_all_times-1;
	      break;
	    }
	    for (k = neq[0]; k--;){
	      ind->solveLast2[k] = yp[k];
	    }
	    canBreak=0;
	  } else if (j >= op->minSS){
	    if (ind->rc[0]== -2019){
	      if (op->strictSS){
		for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
                op->badSolve = 1;
                *i = ind->n_all_times-1;
              } else {
		for (k = neq[0]; k--;){
                  yp[k] = ind->solveLast2[k];
                }
		ind->rc[0] = 2019;
              }
	      break;
	    }
	    for (k = neq[0]; k--;){
	      if (op->ssRtol[k]*fabs(yp[k]) + op->ssAtol[k] <= fabs(yp[k]-ind->solveLast2[k])){
		  canBreak=0;
	      }
	      ind->solveLast2[k] = yp[k];
	    }
	    if (canBreak){
	      break;
	    }
	  }
	  xp2 = xout2;
	}
	*istate=1;
	ind->idx=*i;
	ind->ixds = infBixds;
      }
    }
	  
    if (doSS2){
      // Add at the end
      for (j = neq[0];j--;) yp[j]+=ind->solveSave[j];
    }
    ind->idx=*i;
    if (!doSSinf){
      handle_evid(ind->evid[ind->ix[*i]], neq[0], BadDose, InfusionRate, dose, yp,
		  op->do_transit_abs, xout, neq[1], ind);
    }
    ind->doSS=0;
  }
}

extern void ind_liblsoda0(rx_solve *rx, rx_solving_options *op, struct lsoda_opt_t opt, int solveid, 
			  t_dydt_liblsoda dydt_liblsoda, t_update_inis u_inis){
  int i;
  int neq[2];
  neq[0] = op->neq;
  neq[1] = solveid;
  /* double *yp = &yp0[neq[1]*neq[0]]; */
  int nx;
  rx_solving_options_ind *ind;
  double *inits;
  int *evid;
  double *x;
  int *BadDose;
  double *InfusionRate;
  double *dose;
  double *ret;
  double xout;
  int *rc;
  double *yp;
  inits = op->inits;
  struct lsoda_context_t ctx = {
    .function = dydt_liblsoda,
    .neq = neq[0],
    .data = &neq,
    .state = 1
  };
  lsoda_prepare(&ctx, &opt);
  ind = &(rx->subjects[neq[1]]);
  ind->ixds = 0;
  ind->id = neq[1];
  nx = ind->n_all_times;
  evid = ind->evid;
  BadDose = ind->BadDose;
  InfusionRate = ind->InfusionRate;
  for (int j = neq[0]; j--;) {
    ind->InfusionRate[j] = 0;
    ind->on[j] = 1;
  }
  dose = ind->dose;
  ret = ind->solve;
  x = ind->all_times;
  rc= ind->rc;
  double xp = x[0];
  //--- inits the system
  memcpy(ret,inits, neq[0]*sizeof(double));
  u_inis(neq[1], ret); // Update initial conditions
  unsigned int j;
  if (rx->nMtime) calc_mtime(neq[1], ind->mtime);
  if (rx->needSort) doSort(ind);
  /* for(i=0; i<neq[0]; i++) yp[i] = inits[i]; */
  ind->_newind = 1;
  for(i=0; i<nx; i++) {
    ind->idx=i;
    xout = getTime(ind->ix[i], ind);
    yp = ret+neq[0]*i;
    if(ind->evid[ind->ix[i]] != 3 && xout-xp > DBL_EPSILON*max(fabs(xout),fabs(xp))){
      if (ind->err){
	*rc = -1000;
	// Bad Solve => NA
	for (j = neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	op->badSolve = 1;
	i = nx-1; // Get out of here!
      } else {
	lsoda(&ctx, yp, &xp, xout);
	if (ctx.state <= 0) {
	  /* RSprintf("IDID=%d, %s\n", istate, err_msg_ls[-*istate-1]); */
	  *rc = ctx.state;
	  // Bad Solve => NA
	  for (j = neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	  op->badSolve = 1;
	  i = nx-1; // Get out of here!
	} else if (ind->err){
	  /* RSprintf("IDID=%d, %s\n", istate, err_msg_ls[-*istate-1]); */
	  *rc = ctx.state;
	  // Bad Solve => NA
	  for (j = neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	  op->badSolve = 1;
	  i = nx-1; // Get out of here!
	} else {
	  if (R_FINITE(rx->stateTrim)){
	    double top=fabs(rx->stateTrim);
	    for (unsigned int j = neq[0]; j--;) yp[j]= max(-top, min(top,yp[j]));
	  }
	}
      }
    }
    ind->_newind = 2;
    if (!op->badSolve){
      ind->idx = i;
      if (ind->evid[ind->ix[i]] == 3){
	for (j = neq[0]; j--;) {
	  ind->InfusionRate[j] = 0;
	  ind->on[j] = 1;
	}
	memcpy(yp,inits, neq[0]*sizeof(double));
	u_inis(neq[1], yp); // Update initial conditions @ current time
	if (rx->istateReset) ctx.state = 1;
	xp=xout;
	ind->ixds++;
      } else if (handle_evid(evid[ind->ix[i]], neq[0], BadDose, InfusionRate, dose, yp,
			     op->do_transit_abs, xout, neq[1], ind)){
	handleSS(neq, BadDose, InfusionRate, dose, yp, op->do_transit_abs, xout,
		 xp, ind->id, &i, nx, &ctx.state, op, ind, u_inis, &ctx);
	if (ind->wh0 == 30){
	  ret[ind->cmt] = inits[ind->cmt];
	}
	if (rx->istateReset) ctx.state = 1;
	xp = xout;
      }
      if (i+1 != nx) memcpy(ret+neq[0]*(i+1), yp, neq[0]*sizeof(double));
      ind->slvr_counter[0]++; // doesn't need do be critical; one subject at a time.
      /* for(j=0; j<neq[0]; j++) ret[neq[0]*i+j] = yp[j]; */
    }
  }
  lsoda_free(&ctx);
}

extern void ind_liblsoda(rx_solve *rx, int solveid, 
			 t_dydt_liblsoda dydt, t_update_inis u_inis){
  rx_solving_options *op = &op_global;
  struct lsoda_opt_t opt = {0};
  opt.ixpr = 0; // No extra printing...
  // Unlike traditional lsoda, these are vectors.
  opt.rtol = op->rtol2;
  opt.atol = op->atol2;
  opt.itask = 1;
  opt.mxstep = op->mxstep;
  opt.mxhnil = op->mxhnil;
  opt.mxordn = op->MXORDN;
  opt.mxords = op->MXORDS;
  opt.h0 = op->H0;
  opt.hmax = op->hmax2;
  opt.hmin = op->HMIN;
  opt.hmxi = op->hmxi;
  /* ind_liblsoda0(rx, op, opt, solveid, dydt_liblsoda, update_inis); */
  ind_liblsoda0(rx, op, opt, solveid, dydt, u_inis);
}



extern void par_liblsoda(rx_solve *rx){
  rx_solving_options *op = &op_global;
#ifdef _OPENMP
  int cores = op->cores;
#else
  int cores = 1;
#endif
  int nsub = rx->nsub, nsim = rx->nsim;
  int displayProgress = (op->nDisplayProgress <= nsim*nsub);
  clock_t t0 = clock();
  /* double *yp0=(double*) malloc((op->neq)*nsim*nsub*sizeof(double)); */
  struct lsoda_opt_t opt = {0};
  opt.ixpr = 0; // No extra printing...
  // Unlike traditional lsoda, these are vectors.
  opt.rtol = op->rtol2;
  opt.atol = op->atol2;
  opt.itask = 1;
  opt.mxstep = op->mxstep;
  opt.mxhnil = op->mxhnil;
  opt.mxordn = op->MXORDN;
  opt.mxords = op->MXORDS;
  opt.h0 = op->H0;
  opt.hmax = op->hmax2;
  opt.hmin = op->HMIN;
  opt.hmxi = op->hmxi;
  int curTick=0;
  int cur=0;
  // Breaking of of loop ideas came from http://www.thinkingparallel.com/2007/06/29/breaking-out-of-loops-in-openmp/
  // http://permalink.gmane.org/gmane.comp.lang.r.devel/27627
  // It was buggy due to Rprint.  Use REprint instead since Rprint calls the interrupt every so often....
  int abort = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
  for (int solveid = 0; solveid < nsim*nsub; solveid++){
    if (abort == 0){
      ind_liblsoda0(rx, op, opt, solveid, dydt_liblsoda, update_inis);
      if (displayProgress){
#pragma omp critical
	  cur++;
#ifdef _OPENMP
	  if (omp_get_thread_num() == 0) // only in master thread!
#endif
	    {
	      curTick = par_progress(cur, nsim*nsub, curTick, cores, t0, 0);
	      if (abort == 0){
		if (checkInterrupt()) abort =1;
	      }
	    }
      }
    }
  }
  if (abort == 1){
    op->abort = 1;
    /* yp0 = NULL; */
    par_progress(cur, nsim*nsub, curTick, cores, t0, 1);
  } else {
    if (displayProgress && curTick < 50) par_progress(nsim*nsub, nsim*nsub, curTick, cores, t0, 0);
  }
  if (displayProgress) {
    int doIt = isProgSupported();
    if (doIt == -1){
    } else if (isRstudio() || doIt == 0){
      Rprintf("\n");
    } else {
      RSprintf0("\r                                                                                \r");
    }
  }
}

unsigned int global_rworki = 0;
double *global_rwork(unsigned int mx){ 
  if (mx >= global_rworki){
    global_rworki = mx+1024;
    global_rworkp = Realloc(global_rworkp, global_rworki, double);
  }
  return global_rworkp;
}



unsigned int global_iworki = 0;
int *global_iwork(unsigned int mx){
  if (mx >= global_iworki){
    global_iworki = mx+1024;
    global_iworkp = Realloc(global_iworkp, global_iworki, int);
  }
  return global_iworkp;
}

double *global_InfusionRatep;
unsigned int global_InfusionRatei = 0;
double *global_InfusionRate(unsigned int mx){
  if (mx >= global_InfusionRatei){
    global_InfusionRatei = mx+1024;
    global_InfusionRatep = Realloc(global_InfusionRatep, global_InfusionRatei, double);
  }
  return global_InfusionRatep;
}

double *global_scalep;
unsigned int global_scalei = 0;
double *global_scale(unsigned int mx){
  if (mx >= global_scalei){
    global_scalei = mx+1024;
    global_scalep = Realloc(global_scalep, global_scalei, double);
  }
  return global_scalep;
}


int *global_BadDosep;
unsigned int global_BadDosei = 0;
int *global_BadDose(unsigned int mx){
  if (mx >= global_BadDosei){
    global_BadDosei = mx+1024;
    global_BadDosep = Realloc(global_BadDosep, global_BadDosei, int);
  }
  return global_BadDosep;
}

void rxOptionsIni(){
  global_iworki = 1024*4;
  global_iworkp=Calloc(1024*4, int);
  
  global_rworki=4*1024;
  global_rworkp=Calloc(1024*4, double);
  
  global_InfusionRatei = 1024;
  global_InfusionRatep=Calloc(1024, double);

  global_BadDosei = 1024;
  global_BadDosep=Calloc(1024, int);

  global_scalei = 1024;
  global_scalep=Calloc(1024, double);

  rx_solve *rx=(&rx_global);

  rx->op = &op_global;
  rx->subjects = inds_global;
}

void rxOptionsFree(){
  global_iworki = 0;
  Free(global_iworkp);

  global_rworki = 0;
  Free(global_rworkp);

  global_InfusionRatei = 0;
  Free(global_InfusionRatep);

  global_BadDosei = 0;
  Free(global_BadDosep);

  global_scalei = 0;
  Free(global_scalep);
}

void rxFreeLast(){
  Free(inds_global);
  inds_global=NULL;
}

extern void ind_lsoda0(rx_solve *rx, rx_solving_options *op, int solveid, int *neq, double *rwork, int lrw, int *iwork, int liw, int jt,
                       t_dydt_lsoda_dum dydt_lsoda,
                       t_update_inis u_inis,
                       t_jdum_lsoda jdum){
  rx_solving_options_ind *ind;
  double *yp;
  void *ctx = NULL;
  
  
  int istate = 1, i = 0;
  gitol = 1; gitask = 1; giopt = 1;
  gliw = liw;
  glrw = lrw;

  /* memset(rwork,0.0,lrw+1); */ // Does not work since it is a double
  for (i = lrw+1; i--;) rwork[i]=0;
  memset(iwork,0,liw+1); // Works because it is a integer

  neq[1] = solveid;
  
  ind = &(rx->subjects[neq[1]]);
  ind->id = neq[1];

  rwork[4] = op->H0; // H0
  rwork[5] = ind->HMAX; // Hmax
  rwork[6] = op->HMIN; // Hmin
  
  iwork[4] = 0; // ixpr
  iwork[5] = op->mxstep; // mxstep 
  iwork[6] = op->mxhnil; // MXHNIL 
  iwork[7] = op->MXORDN; // MXORDN 
  iwork[8] = op->MXORDS;  // MXORDS
    
  ind->ixds = 0;
  double xp = ind->all_times[0];
  double xout;

  //--- inits the system
  for (int j = neq[0]; j--;) {
    ind->InfusionRate[j] = 0;
    ind->on[j] = 1;
  }
  memcpy(ind->solve, op->inits, neq[0]*sizeof(double));
  u_inis(neq[1], ind->solve); // Update initial conditions
  if (rx->nMtime) calc_mtime(neq[1], ind->mtime);
  if (rx->needSort) doSort(ind);
  unsigned int j;
  ind->_newind = 1;
  for(i=0; i < ind->n_all_times; i++) {
    ind->idx=i;
    xout = getTime(ind->ix[i], ind);
    yp   = ind->solve+neq[0]*i;
    if(ind->evid[ind->ix[i]] != 3 && xout - xp > DBL_EPSILON*max(fabs(xout),fabs(xp))) {
      if (ind->err){
	ind->rc[0] = -1000;
	// Bad Solve => NA
	for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	op->badSolve = 1;
	i = ind->n_all_times-1; // Get out of here!
      } else {
	F77_CALL(dlsoda)(dydt_lsoda, neq, yp, &xp, &xout, &gitol, &(op->RTOL), &(op->ATOL), &gitask,
			 &istate, &giopt, rwork, &lrw, iwork, &liw, jdum, &jt);
	if (istate <= 0) {
	  RSprintf("IDID=%d, %s\n", istate, err_msg_ls[-(istate)-1]);
	  ind->rc[0] = istate;
	  // Bad Solve => NA
	  for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	  op->badSolve = 1;
	  i = ind->n_all_times-1; // Get out of here!
	} else if (ind->err){
	  ind->rc[0] = -1000;
	  // Bad Solve => NA
	  for (j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	  op->badSolve = 1;
	  i = ind->n_all_times-1; // Get out of here!
	} else {
	  if (R_FINITE(rx->stateTrim)){
	    double top=fabs(rx->stateTrim);
	    for (j = neq[0]; j--;) yp[j]= max(-top, min(top,yp[j]));
	  }
	}
	ind->slvr_counter[0]++;
	//dadt_counter = 0;
      }
    }
    ind->_newind = 2;
    if (!op->badSolve){
      ind->idx = i;
      if (ind->evid[ind->ix[i]] == 3){
	for (j = neq[0]; j--;) {
	  ind->InfusionRate[j] = 0;
	  ind->on[j] = 1;
	}
	memcpy(yp, op->inits, neq[0]*sizeof(double));
	u_inis(neq[1], yp); // Update initial conditions @ current time
	if (rx->istateReset) istate = 1;
	ind->ixds++;
	xp = xout;
      } else if (handle_evid(ind->evid[ind->ix[i]], neq[0], ind->BadDose, ind->InfusionRate, ind->dose, yp,
			     op->do_transit_abs, xout, neq[1], ind)){
	handleSS(neq, ind->BadDose, ind->InfusionRate, ind->dose, yp, op->do_transit_abs, xout,
		 xp, ind->id, &i, ind->n_all_times, &istate, op, ind, u_inis, &ctx);
	if (ind->wh0 == 30){
	  ind->solve[ind->cmt] = op->inits[ind->cmt];
	}
	if (rx->istateReset) istate = 1;
	xp = xout;
      }
      // Copy to next solve so when assigned to yp=ind->solve[neq[0]*i]; it will be the prior values
      if (i+1 != ind->n_all_times) memcpy(ind->solve+neq[0]*(i+1), yp, neq[0]*sizeof(double));
    }
  }
}

extern void ind_lsoda(rx_solve *rx, int solveid,
                      t_dydt_lsoda_dum dydt_ls, t_update_inis u_inis, t_jdum_lsoda jdum,
		      int cjt){
  int neq[2];
  neq[0] = op_global.neq;
  neq[1] = 0;
  
  // Set jt to 1 if full is specified.
  int lrw=22+neq[0]*max(16, neq[0]+9), liw=20+neq[0];
  double *rwork;
  int *iwork;
  if (global_debug)
    RSprintf("JT: %d\n",cjt);
  rwork = global_rwork(lrw+1);
  iwork = global_iwork(liw+1);
  ind_lsoda0(rx, &op_global, solveid, neq, rwork, lrw, iwork, liw, cjt,
             dydt_ls, u_inis, jdum);
}

extern void par_lsoda(rx_solve *rx){
  int nsub = rx->nsub, nsim = rx->nsim;
  int displayProgress = (op_global.nDisplayProgress <= nsim*nsub);
  clock_t t0 = clock();
  int neq[2];
  neq[0] = op_global.neq;
  neq[1] = 0;
  /* yp = global_yp(neq[0]); */
  
  // Set jt to 1 if full is specified.
  int lrw=22+neq[0]*max(16, neq[0]+9), liw=20+neq[0], jt = global_jt;
  double *rwork;
  int *iwork;
  
  
  if (global_debug)
    RSprintf("JT: %d\n",jt);
  rwork = global_rwork(lrw+1);
  iwork = global_iwork(liw+1);
  
  int curTick = 0;
  int abort = 0;
  for (int solveid = 0; solveid < nsim*nsub; solveid++){
    ind_lsoda0(rx, &op_global, solveid, neq, rwork, lrw, iwork, liw, jt,
	       dydt_lsoda_dum, update_inis, jdum_lsoda);
    if (displayProgress){ // Can only abort if it is long enough to display progress.
      curTick = par_progress(solveid, nsim*nsub, curTick, 1, t0, 0);
      if (checkInterrupt()){
	abort =1;
	break;
      }
    }
  }
  if (abort == 1){
    op_global.abort = 1;
  } else {
    if (displayProgress && curTick < 50) par_progress(nsim*nsub, nsim*nsub, curTick, 1, t0, 0);
  }
}

extern void ind_dop0(rx_solve *rx, rx_solving_options *op, int solveid, int *neq, 
                     t_dydt c_dydt,
                     t_update_inis u_inis){
  double rtol=op->RTOL, atol=op->ATOL;
  int itol=0;           //0: rtol/atol scalars; 1: rtol/atol vectors
  int iout=0;           //iout=0: solout() NEVER called
  int idid=0;
  int i;
  double xout;
  double *yp;
  void *ctx = NULL;
  int istate = 0;
  static char *err_msg[]=
    {
      "input is not consistent",
      "larger nmax is needed",
      "step size becomes too small",
      "problem is probably stiff (interrupted)"
    };
  rx_solving_options_ind *ind;
  int *evid;
  double *x;
  int *BadDose;
  double *InfusionRate;
  double *dose;
  double *ret, *inits;
  int *rc;
  int nx;
  neq[1] = solveid;
  ind = &(rx->subjects[neq[1]]);
  ind->id = neq[1];
  ind->ixds = 0;
  nx = ind->n_all_times;
  inits = op->inits;
  evid = ind->evid;
  BadDose = ind->BadDose;
  InfusionRate = ind->InfusionRate;
  dose = ind->dose;
  ret = ind->solve;
  x = ind->all_times;
  rc= ind->rc;
  double xp = x[0];
  //--- inits the system
  for (int j = neq[0]; j--;) {
    ind->InfusionRate[j] = 0;
    ind->on[j] = 1;
  }
  memcpy(ret,inits, neq[0]*sizeof(double));
  u_inis(neq[1], ret); // Update initial conditions
  if (rx->nMtime) calc_mtime(neq[1], ind->mtime);
  if (rx->needSort) doSort(ind);
  //--- inits the system
  unsigned int j;
  ind->_newind = 1;
  for(i=0; i<nx; i++) {
    ind->idx=i;
    xout = getTime(ind->ix[i], ind);
    yp = &ret[neq[0]*i];
    if (global_debug){
      RSprintf("i=%d xp=%f xout=%f\n", i, xp, xout);
    }
    if(ind->evid[ind->ix[i]] != 3 && xout-xp>DBL_EPSILON*max(fabs(xout),fabs(xp)))
      {
	if (ind->err){
	  printErr(ind->err, ind->id);
	  *rc = idid;
	  // Bad Solve => NA
	  for (j = (ind->n_all_times)*neq[0];j--;) ret[i] = NA_REAL; 
	  op->badSolve = 1;
	  i = nx-1; // Get out of here!
	} else {
	  idid = dop853(neq,       /* dimension of the system <= UINT_MAX-1*/
			c_dydt,       /* function computing the value of f(x,y) */
			xp,           /* initial x-value */
			yp,           /* initial values for y */
			xout,         /* final x-value (xend-x may be positive or negative) */
			&rtol,          /* relative error tolerance */
			&atol,          /* absolute error tolerance */
			itol,         /* switch for rtoler and atoler */
			solout,         /* function providing the numerical solution during integration */
			iout,         /* switch for calling solout */
			NULL,           /* messages stream */
			DBL_EPSILON,    /* rounding unit */
			0,              /* safety factor */
			0,              /* parameters for step size selection */
			0,
			0,              /* for stabilized step size control */
			0,              /* maximal step size */
			0,            /* initial step size */
			op->mxstep, /* maximal number of allowed steps */
			1,            /* switch for the choice of the coefficients */
			-1,                     /* test for stiffness */
			0,                      /* number of components for which dense outpout is required */
			NULL,           /* indexes of components for which dense output is required, >= nrdens */
			0                       /* declared length of icon */
			);
	}
        if (idid<0) {
	  RSprintf("IDID=%d, %s\n", idid, err_msg[-idid-1]);
	  *rc = idid;
	  // Bad Solve => NA
	  for (j = (ind->n_all_times)*neq[0];j--;) ret[i] = NA_REAL; 
	  op->badSolve = 1;
	  i = nx-1; // Get out of here!
	} else if (ind->err){
	  printErr(ind->err, ind->id);
	  *rc = idid;
	  // Bad Solve => NA
	  for (j = (ind->n_all_times)*neq[0];j--;) ret[i] = NA_REAL; 
	  op->badSolve = 1;
	  i = nx-1; // Get out of here!
	} else {
	  if (R_FINITE(rx->stateTrim)){
	    double top=fabs(rx->stateTrim);
	    for (j = neq[0]; j--;) yp[j]= max(-top, min(top,yp[j]));
	  }
	}
        xp = xRead();
        ind->slvr_counter[0]++;
        //dadt_counter = 0;
      }
    ind->_newind = 1;
    if (!op->badSolve){
      ind->idx = i;
      if (ind->evid[ind->ix[i]] == 3){
	for (j = neq[0]; j--;) {
	  ind->InfusionRate[j] = 0;
	  ind->on[j] = 1;
	}
	memcpy(yp, op->inits, neq[0]*sizeof(double));
	u_inis(neq[1], yp); // Update initial conditions @ current time
	ind->ixds++;
	xp=xout;
      } else if (handle_evid(evid[ind->ix[i]], neq[0], BadDose, InfusionRate, dose, yp,
			     op->do_transit_abs, xout, neq[1], ind)){
	handleSS(neq, BadDose, InfusionRate, dose, yp, op->do_transit_abs, xout,
		 xp, ind->id, &i, nx, &istate, op, ind, u_inis, &ctx);
	if (ind->wh0 == 30){
	  ret[ind->cmt] = inits[ind->cmt];
	}
	xp = xout;
      }
      /* for(j=0; j<neq[0]; j++) ret[neq[0]*i+j] = yp[j]; */
      if (i+1 != nx) memcpy(ret+neq[0]*(i+1), ret + neq[0]*i, neq[0]*sizeof(double));
    }
  }
}

extern void ind_dop(rx_solve *rx, int solveid,
		    t_dydt c_dydt, t_update_inis u_inis){
  rx_solving_options *op = &op_global;
  int neq[2];
  neq[0] = op->neq;
  neq[1] = 0;
  ind_dop0(rx, &op_global, solveid, neq, c_dydt, u_inis);
}

void par_dop(rx_solve *rx){
  rx_solving_options *op = &op_global;
  int nsub = rx->nsub, nsim = rx->nsim;
  int displayProgress = (op->nDisplayProgress <= nsim*nsub);
  clock_t t0 = clock();
  int neq[2];
  neq[0] = op->neq;
  neq[1] = 0;
  
  //DE solver config vars
  // This part CAN be parallelized, if dop is thread safe...
  // Therefore you could use https://github.com/jacobwilliams/dop853, but I haven't yet
  
  int curTick = 0;
  int abort = 0;
  for (int solveid = 0; solveid < nsim*nsub; solveid++){
    if (abort == 0){
      ind_dop0(rx, &op_global, solveid, neq, dydt, update_inis);
      if (displayProgress && abort == 0){
        if (checkInterrupt()) abort =1;
      }
      if (displayProgress) curTick = par_progress(solveid, nsim*nsub, curTick, 1, t0, 0);
    }
  }
  if (abort == 1){
    op->abort = 1;
  } else {
    if (displayProgress && curTick < 50) par_progress(nsim*nsub, nsim*nsub, curTick, 1, t0, 0);
  }
  if (displayProgress){
    int doIt = isProgSupported();
    if (doIt == -1){
    } else if (isRstudio() || doIt == 0){
      Rprintf("\n");
    } else {
      RSprintf0("\r                                                                                \r");
    }
  }
}

void ind_solve(rx_solve *rx, unsigned int cid,
	       t_dydt_liblsoda dydt_lls,
	       t_dydt_lsoda_dum dydt_lsoda, t_jdum_lsoda jdum,
	       t_dydt c_dydt, t_update_inis u_inis,
	       int jt){
  assignFuns();
  rx_solving_options *op = &op_global;
  if (op->neq !=  0){
    switch (op->stiff){
    case 2: 
      ind_liblsoda(rx, cid, dydt_lls, u_inis);
      break;
    case 1:
      ind_lsoda(rx,cid, dydt_lsoda, u_inis, jdum, jt);
      break;
    case 0:
      ind_dop(rx, cid, c_dydt, u_inis);
      break;
    }
  } 
}

inline void par_solve(rx_solve *rx){
  assignFuns();
  rx_solving_options *op = &op_global;
  if (op->neq != 0){
    switch(op->stiff){
    case 2:
      par_liblsoda(rx);
      break;
    case 1:
      // lsoda
      par_lsoda(rx);
      break;
    case 0:
      // dop
      par_dop(rx);
      break;
    }
  } 
}

rx_solve *_globalRx = NULL;

extern void rxode_assign_rx(rx_solve *rx){
  _globalRx=rx;
}


extern double rxLhsP(int i, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  rx_solving_options *op = &op_global;
  if (i < op->nlhs){
    return(ind->lhs[i]);
  } else {
    error("Trying to access an equation that isn't calculated. lhs(%d/%d)\n",i, op->nlhs);
  }
  return 0;
}
extern void rxCalcLhsP(int i, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  rx_solving_options *op = &op_global;
  double *solve, *lhs;
  solve = ind->solve;
  lhs = ind->lhs;
  if (i < ind->n_all_times){
    ind->idx=i;
    if (ind->evid[ind->ix[i]]) ind->tlast = getTime(ind->ix[i], ind);
    calc_lhs((int)id, getTime(ind->ix[i], ind), solve+i*op->neq, lhs);
  } else {
    error("LHS cannot be calculated (%dth entry).",i);
  }
}
extern void setExtraCmtP(int xtra, rx_solve *rx){
  rx_solving_options *op = &op_global;
  if (xtra > op->extraCmt){
    op->extraCmt = xtra;
  }
}

void setExtraCmt(int xtra){
  setExtraCmtP(xtra, _globalRx);
}

SEXP rxStateNames(char *ptr);
SEXP rxLhsNames(char *ptr);
SEXP rxParamNames(char *ptr);

extern double *rxGetErrs();
extern int rxGetErrsNcol();

extern double get_ikeep(int col, int id);
extern const SEXP get_ikeepn();
extern double get_fkeep(int col, int id);
extern const SEXP get_fkeepn();

extern SEXP RxODE_df(int doDose0, int doTBS){
  rx_solve *rx;
  rx = &rx_global;
  rx_solving_options *op = &op_global;
  int add_cov = rx->add_cov;
  int ncov = op->ncov;
  int ncov0 = rx->nCov0;
  int nkeep0 = rx->nKeep0;
  int nkeep  = rx->nKeepF;
  int nlhs = op->nlhs;
  int nobs = rx->nobs - rx->nevid9;
  int nsim = rx->nsim;
  int nall = rx->nall - rx->nevid9;
  int errNcol = rxGetErrsNcol();
  if (op->nsvar != errNcol){
    error("The simulated residual errors do not match the model specification (%d=%d)",op->nsvar, errNcol);
  }
  int doDose;
  int evid0 = 0;
  int nmevid=0;
  int subsetEvid = 0;
  if (doDose0 == -1){
    nobs = rx->nobs2;
    doDose=0;
    evid0=1;
  } else if (doDose0 == 2 || doDose0 == 3){
    // rate dur ii ss
    doDose=1;
    nmevid=1;
    if (doDose0 == 3){
      subsetEvid=1;
      doDose0 = 2;
    }
  } else {
    doDose=doDose0;
  }
  int di = 0;
  double *dose;
  double *dfp;
  int *dfi;
  int ii=0, jj = 0, ntimes;
  double *solve;
  double *cov_ptr;
  int nBadDose;
  int *BadDose;
  int extraCmt = op->extraCmt;
  int *svar = op->svar;
  int kk = 0;
  int wh, cmt, wh100, whI, wh0;
  int //dullEvid = 1,
    dullRate=1, dullDur=1,
    dullSS=1, dullIi=1;
  int csub = 0, evid;
  int nsub = rx->nsub;
  int *rmState = rx->stateIgnore;
  int nPrnState =0;
  int i, j;
  int neq[2];
  double *scale;
  rx_solving_options_ind *ind;  
  if (subsetEvid == 1){
    rx->nr=0;
    for (int csim = 0; csim < nsim; csim++){
      for (csub = 0; csub < nsub; csub++){
	neq[1] = csub+csim*nsub;
	ind = &(rx->subjects[neq[1]]);
	ind->id = neq[1];
	ntimes = ind->n_all_times;
	dose = ind->dose;
	di = 0;
	for (i = 0; i < ntimes; i++){
	  evid = ind->evid[ind->ix[i]];
	  if (isDose(evid)){
	    getWh(evid, &wh, &cmt, &wh100, &whI, &wh0);
	    if (whI != 7  && whI != 6){
	      if (dose[di++] > 0){
		rx->nr++;
	      }
	    } else {
	      di++;
	    }
	  } else if (isObs(evid)){
	    if (evid < 9){
	      rx->nr++;
	    }
	  }
	}
      }
    }
    di = 0;
  } else {
      rx->nr = (doDose == 1 ? nall : nobs)*nsim;
  }
  scale = op->scale;
  neq[0] = op->neq;
  neq[1] = 0;
  for (i = 0; i < neq[0]; i++){
    nPrnState+= (1-rmState[i]);
  }
  // Mutiple ID data?
  int md = 0;
  if (rx->nsub > 1) md = 1;
  // Multiple simulation data?
  int sm = 0;
  if (rx->nsim > 1) sm = 1;
  int ncols =add_cov*(ncov+ncov0)+nkeep0+nkeep+1+nPrnState+nlhs;
  int doseCols = 0;
  if (doDose){
    doseCols = 2;
  }
  int nidCols = md + sm;
  int pro = 0;
  if (op->badSolve){
    if (nidCols == 0){
      error("Could not solve the system.");
    } else {
      warning("Some ID(s) could not solve the ODEs correctly; These values are replaced with NA.");
    }
  }  
  SEXP df = PROTECT(allocVector(VECSXP,ncols+nidCols+doseCols+doTBS*2+5*nmevid)); pro++;
  for (i = nidCols; i--;){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(INTSXP, rx->nr))); pro++;
  }
  i = nidCols;
  double *par_ptr;
  double *errs = rxGetErrs();
  int updateErr = 0;
  
  if (errNcol > 0){
    updateErr = 1;
  }
  if (doDose){
    //evid
    SET_VECTOR_ELT(df, i++, PROTECT(allocVector(INTSXP, rx->nr))); pro++;
    if (nmevid){
      // cmt
      SET_VECTOR_ELT(df, i++, PROTECT(allocVector(INTSXP, rx->nr))); pro++;
      // ss
      SET_VECTOR_ELT(df, i++, PROTECT(allocVector(INTSXP, rx->nr))); pro++;
    }
    // amt
    SET_VECTOR_ELT(df, i++, PROTECT(allocVector(REALSXP, rx->nr))); pro++;
  }
  for (i = md + sm + doseCols + 2*nmevid; i < ncols + doseCols + nidCols + 2*nmevid; i++){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(REALSXP, rx->nr))); pro++;
  }
  for (i = ncols + doseCols + nidCols + 2*nmevid; i < ncols + doseCols + nidCols + doTBS*2 + nmevid*5; i++){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(REALSXP, rx->nr))); pro++;
  }
  // Now create the data frame
  for (int csim = 0; csim < nsim; csim++){
    for (csub = 0; csub < nsub; csub++){
      neq[1] = csub+csim*nsub;
      ind = &(rx->subjects[neq[1]]);
      ind->id = neq[1];
      ind->_newind = 1;
      if (rx->nMtime) calc_mtime(neq[1], ind->mtime);
      if (rx->needSort) doSort(ind);
      nBadDose = ind->nBadDose;
      BadDose = ind->BadDose;
      ntimes = ind->n_all_times;
      solve =  ind->solve;
      cov_ptr = ind->cov_ptr;
      par_ptr = ind->par_ptr;
      dose = ind->dose;
      di = 0;
      if (nBadDose && csim == 0){
	for (i = 0; i < nBadDose; i++){
	  if (BadDose[i] > extraCmt){
	    warning("Dose to Compartment %d ignored (not in ODE; id=%d)", BadDose[i],csub+1);
	  }
	}
      }
      if (ind->allCovWarn && csim == 0){
	warning("One or more covariates were all NA for subject id=%d", csub+1);
      }	
      for (i = 0; i < ntimes; i++){
        evid = ind->evid[ind->ix[i]];
	if (evid == 9) continue;
	if (subsetEvid == 1){
	  if (isObs(evid) && evid >= 10) continue;
	  if (isDose(evid)){
	    getWh(evid, &wh, &cmt, &wh100, &whI, &wh0);
	    if (whI == 7 || whI == 6){
	      dullRate=0;
	      di++;
	      continue;
	    }
	    if (dose[di] <= 0){
	      di++;
	      continue;
	    }
	  }
	}
	if (isDose(evid)) ind->tlast = getTime(ind->ix[i], ind);
        if (updateErr){
          for (j=0; j < errNcol; j++){
	    par_ptr[svar[j]] = errs[rx->nr*j+kk];
          }
	  if ( (evid0 == 0 && isObs(evid)) || (evid0 == 1 && evid==0) || doDose){
	    // Only incerement if this is an observation or of this a
	    // simulation that requests dosing information too.
            kk++;
	  }
        }
        jj  = 0 ;
	if ((evid0 == 0 && isObs(evid)) || (evid0 == 1 && evid==0)  || doDose){
          // sim.id
          if (sm){
            dfi = INTEGER(VECTOR_ELT(df, jj));
            dfi[ii] = csim+1;
            jj++;
          }
	  // id
          if (md){
            dfi = INTEGER(VECTOR_ELT(df, jj));
            dfi[ii] = csub+1;
            jj++;
          }
	  if (doDose){
	    if (nmevid){
	      if (isObs(evid)){
		// evid
		dfi = INTEGER(VECTOR_ELT(df, jj++));
		if (evid >= 10){
		  dfi[ii] = evid+91; // mtime 101 102 103...
		  /* dullEvid=0; */
		} else {
		  /* if (evid == 2) dullEvid=0; */
		  dfi[ii] = evid;
		}
		// cmt
		dfi = INTEGER(VECTOR_ELT(df, jj++));
		dfi[ii] = NA_INTEGER; // Has all states, cmt makes no sense.
		// ss
		dfi = INTEGER(VECTOR_ELT(df, jj++));
		dfi[ii] = 0;
		// amt
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// rate
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// dur
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// ii
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
	      } else {
		getWh(evid, &wh, &cmt, &wh100, &whI, &wh0);
		dfi = INTEGER(VECTOR_ELT(df, jj++));
		double curAmt = dose[di];
		if (whI == 7){
		  dullRate=0;
		  dfi[ii] = -1;
		} else if (whI == 6){
		  dullRate=0;
		  dfi[ii] = -2; // evid
		} else {
		  if (curAmt > 0) {
		    if (whI == 4){
		      dfi[ii] = 5;
		    } else if (whI == 5){
		      dfi[ii] = 6;
		    } else {
		      dfi[ii] = 1; // evid
		    }
		  } else {
		    if (whI == 1){
		      dullRate=0;
		      dfi[ii] = -10; // evid
		    } else if (whI == 2) {
		      dullDur=0;
		      dfi[ii] = -20; // evid
		    } else if (whI == 4){
		      dfi[ii] = 5;
		    } else if (whI == 5){
		      dfi[ii] = 6;
		    } else {
		      dfi[ii] = 1;
		    }
		  }
		}
		// cmt
		dfi = INTEGER(VECTOR_ELT(df, jj++));
		if (evid == 2 || evid == 3){
		  dfi[ii] = NA_INTEGER;
		} else if (wh0 == 30){
		  dfi[ii] = -cmt-1;
		} else {
		  dfi[ii] = cmt+1;
		}
		// ss
		dfi = INTEGER(VECTOR_ELT(df, jj++));
		switch (wh0){
		/* case 30: */
		case 20:
		  dullSS=0;
		  dfi[ii] = 2;
		  break;
		case 10:
		  dullSS=0;
		  dfi[ii] = 1;
		  break;
		default:
		  dfi[ii] = 0;
		  break;
		}
	      }
	    } else {
	      // evid
	      dfi = INTEGER(VECTOR_ELT(df, jj++));
	      dfi[ii] = evid;
	      // amt
	      dfp = REAL(VECTOR_ELT(df, jj++));
	      dfp[ii] = isObs(evid) ? NA_REAL : dose[di++];
	    }
	    if (nmevid && isDose(evid)){
	      double curIi = ind->ii[di];
	      if (curIi != 0) dullIi=0;
	      double curAmt = dose[di++];
	      // rate dur ii ss
	      switch(ind->whI){
	      case 9: // modeled rate
		// amt
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = curAmt;
		// rate
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = -1.0;
		// dur
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// ii
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = curIi;
		break;
	      case 8: // modeled duration
		// amt
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = curAmt;
		// rate
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = -2.0;
		// dur
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// ii
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = curIi;
		break;
	      case 7: // End modeled rate
		// amt
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// rate
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// dur
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// ii
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = curIi;
		break;
	      case 6: // end modeled duration
		// amt
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// rate
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// dur
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// ii
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = curIi;
		break;
	      case 2: // Infusion specified by dur
		if (curAmt < 0){
		  // amt
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = NA_REAL;
		  // rate
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = NA_REAL;
		  // dur
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = NA_REAL;
		} else {
		  // Find the next fixed length infusion that is turned off.
		  double curDur=0.0;
		  for (int jjj = di; jjj < ind->ndoses; jjj++){
		    if (ind->dose[jjj] == -curAmt){
		      int nWh = 0, nCmt = 0, nWh100 = 0, nWhI = 0, nWh0 = 0;
		      getWh(ind->evid[ind->idose[jjj]], &nWh, &nCmt, &nWh100, &nWhI, &nWh0);
		      if (nWhI == whI && nCmt == cmt){
			curDur = getTime(ind->idose[jjj], ind) - getTime(ind->ix[i], ind);
			break;
		      }
		    }
		  }
		  // amt
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = curAmt*curDur;
		  // rate
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = NA_REAL;
		  // dur
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = curDur;
		}
		// ii
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = curIi;		  
		break;
	      case 1: // Infusion specified by rate
		if (curAmt < 0){
		  // amt
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = NA_REAL;
		  // rate
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = NA_REAL;
		  // dur
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = NA_REAL;
		} else {
		  double curDur=0.0;
		  for (int jjj = di; jjj < ind->ndoses; jjj++){
		    if (ind->dose[jjj] == -curAmt){
		      int nWh = 0, nCmt = 0, nWh100 = 0, nWhI = 0, nWh0 = 0;
		      getWh(ind->evid[ind->idose[jjj]], &nWh, &nCmt, &nWh100, &nWhI, &nWh0);
		      if (nWhI == whI && nCmt == cmt){
			curDur = getTime(ind->idose[jjj], ind) - getTime(ind->ix[i], ind);
			break;
		      }
		    }
		  }
		  // amt
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = curAmt*curDur;
		  // rate
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = curAmt;
		  // dur
		  dfp = REAL(VECTOR_ELT(df, jj++));
		  dfp[ii] = NA_REAL;
		}
		// ii
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = curIi;
		break;
	      default:
		// Non infusion dose.
		// amt
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = curAmt;
		// rate
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// dur
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = NA_REAL;
		// ii
		dfp = REAL(VECTOR_ELT(df, jj++));
		dfp[ii] = curIi;
	      }
	    }
	  }
          // time
          dfp = REAL(VECTOR_ELT(df, jj++));
          dfp[ii] = getTime(ind->ix[i], ind);
          // LHS
          if (nlhs){
	    rxCalcLhsP(i, rx, neq[1]);
	    for (j = 0; j < nlhs; j++){
	      dfp = REAL(VECTOR_ELT(df, jj));
               dfp[ii] =rxLhsP(j, rx, neq[1]);
	       jj++;
             }
          }
          // States
          if (nPrnState){
            for (j = 0; j < neq[0]; j++){
              if (!rmState[j]){
                dfp = REAL(VECTOR_ELT(df, jj));
                dfp[ii] = solve[j+i*neq[0]] / scale[j];
                jj++;
              }
            }
          }       
          // Cov
          if (add_cov*ncov > 0){
	    for (j = 0; j < add_cov*ncov; j++){
              dfp = REAL(VECTOR_ELT(df, jj));
	      // is this ntimes = nAllTimes or nObs time for this subject...?
	      dfp[ii] = isObs(evid)  ? cov_ptr[j*ntimes+i] : NA_REAL;
	      jj++;
	    }
          }
	  if (add_cov*ncov0 > 0){
	    for (j = 0; j < add_cov*ncov0; j++){
              dfp = REAL(VECTOR_ELT(df, jj));
	      // is this ntimes = nAllTimes or nObs time for this subject...?
	      dfp[ii] = isObs(evid) ? ind->par_ptr[rx->cov0[j]] : NA_REAL;
	      jj++;
	    }
	  }
	  for (j = 0; j < nkeep0; j++){
	    dfp = REAL(VECTOR_ELT(df, jj));
	    dfp[ii] = get_ikeep(j, neq[1]);
	  }
	  for (j = 0; j < nkeep; j++){
	    dfp = REAL(VECTOR_ELT(df, jj));
	    // is this ntimes = nAllTimes or nObs time for this subject...?
	    dfp[ii] = get_fkeep(j, i);
	    jj++;
	  }
	  // 
	  if (doTBS){
	    dfp = REAL(VECTOR_ELT(df, jj));
	    dfp[ii] = ind->lambda;
	    jj++;
            dfp = REAL(VECTOR_ELT(df, jj));
            dfp[ii] = ind->yj;
	    jj++;
	  }
          ii++;
        }
	ind->_newind = 2;
      }
      if (updateErr){
        for (j=0; j < errNcol; j++){
          par_ptr[svar[j]] = NA_REAL;
        }
      }
    }
  }
  SEXP sexp_rownames = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(sexp_rownames)[0] = NA_INTEGER;
  INTEGER(sexp_rownames)[1] = -rx->nr;
  setAttrib(df, R_RowNamesSymbol, sexp_rownames);
  SEXP sexp_colnames = PROTECT(allocVector(STRSXP,ncols+nidCols+doseCols+doTBS*2+5*nmevid)); pro++;
  jj = 0;
  if (sm){
    SET_STRING_ELT(sexp_colnames, jj, mkChar("sim.id"));
    jj++;
  }
  // id
  if (md){
    SET_STRING_ELT(sexp_colnames, jj, mkChar("id"));
    jj++;
  }
  if (doDose){
    SET_STRING_ELT(sexp_colnames, jj, mkChar("evid"));
    jj++;
    if (nmevid){
      SET_STRING_ELT(sexp_colnames, jj, mkChar("cmt"));
      jj++;
      SET_STRING_ELT(sexp_colnames, jj, mkChar("ss"));
      jj++;
    }
    SET_STRING_ELT(sexp_colnames, jj, mkChar("amt"));
    jj++;
    if (nmevid){
      SET_STRING_ELT(sexp_colnames, jj, mkChar("rate"));
      jj++;
      SET_STRING_ELT(sexp_colnames, jj, mkChar("dur"));
      jj++;
      SET_STRING_ELT(sexp_colnames, jj, mkChar("ii"));
      jj++;
    }
    
  }
  SET_STRING_ELT(sexp_colnames, jj, mkChar("time"));
  jj++;

  // Put in LHS names
  SEXP lhsNames = PROTECT(rxLhsNames(op->modNamePtr)); pro++;
  for (i = 0; i < nlhs; i++){
    SET_STRING_ELT(sexp_colnames, jj, STRING_ELT(lhsNames,i));
    jj++;
  }
  // Put in state names
  SEXP stateNames = PROTECT(rxStateNames(op->modNamePtr)); pro++;
  if (nPrnState){
    for (j = 0; j < neq[0]; j++){
      if (!rmState[j]){
        SET_STRING_ELT(sexp_colnames, jj, STRING_ELT(stateNames,j));
        jj++;
      }
    }
  }  
  // Put in Cov names
  SEXP paramNames = PROTECT(rxParamNames(op->modNamePtr)); pro++;
  int *par_cov = op->par_cov;
  for (i = 0; i < ncov*add_cov; i++){
    SET_STRING_ELT(sexp_colnames,jj, STRING_ELT(paramNames, par_cov[i]-1));
    jj++;
  }
  par_cov = rx->cov0;
  for (i = 0; i < ncov0*add_cov; i++){
    SET_STRING_ELT(sexp_colnames,jj, STRING_ELT(paramNames, par_cov[i]));
    jj++;
  }
  SEXP ikeepNames = PROTECT(get_ikeepn()); pro++;
  for (i = 0; i < nkeep0; i++){
    SET_STRING_ELT(sexp_colnames,jj, STRING_ELT(ikeepNames, i));
    jj++;
  }
  SEXP fkeepNames = PROTECT(get_fkeepn()); pro++;
  for (i = 0; i < nkeep; i++){
    SET_STRING_ELT(sexp_colnames,jj, STRING_ELT(fkeepNames, i));
    jj++;
  }
  if (doTBS){
    SET_STRING_ELT(sexp_colnames, jj, mkChar("rxLambda"));
    jj++;
    SET_STRING_ELT(sexp_colnames, jj, mkChar("rxYj"));
    jj++;
  }
  setAttrib(df, R_NamesSymbol, sexp_colnames);
  SEXP df2;
  if (nmevid){
    df2 = PROTECT(allocVector(VECSXP,ncols+nidCols+doseCols+doTBS*2+5*nmevid-
			      dullRate - dullDur-dullSS-dullIi)); pro++;
    SEXP sexp_colnames2 = PROTECT(allocVector(STRSXP,ncols+nidCols+doseCols+doTBS*2+5*nmevid-
					      dullRate - dullDur-dullSS-dullIi)); pro++;
    jj = 0;
    kk = 0;
    if (sm){
      SET_STRING_ELT(sexp_colnames2, jj, mkChar("sim.id"));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
    }
    // id
    if (md){
      SET_STRING_ELT(sexp_colnames2, jj, mkChar("id"));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
    }
    SET_STRING_ELT(sexp_colnames2, jj, mkChar("evid"));
    SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
    jj++;kk++;
    SET_STRING_ELT(sexp_colnames2, jj, mkChar("cmt"));
    SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
    jj++;kk++;
    if (dullSS){
      kk++;
    } else {
      SET_STRING_ELT(sexp_colnames2, jj, mkChar("ss"));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
    }
    SET_STRING_ELT(sexp_colnames2, jj, mkChar("amt"));
    SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
    jj++;kk++;
    if (dullRate){
      kk++;
    } else {
      SET_STRING_ELT(sexp_colnames2, jj, mkChar("rate"));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
    }
    if (dullDur){
      kk++;
    } else {
      SET_STRING_ELT(sexp_colnames2, jj, mkChar("dur"));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
    }
    if (dullIi){
      kk++;
    } else {
      SET_STRING_ELT(sexp_colnames2, jj, mkChar("ii"));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
    }
    SET_STRING_ELT(sexp_colnames2, jj, mkChar("time"));
    SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
    jj++;kk++;

    // Put in LHS names
    SEXP lhsNames2 = PROTECT(rxLhsNames(op->modNamePtr)); pro++;
    for (i = 0; i < nlhs; i++){
      SET_STRING_ELT(sexp_colnames2, jj, STRING_ELT(lhsNames2,i));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
    }
    // Put in state names
    SEXP stateNames2 = PROTECT(rxStateNames(op->modNamePtr)); pro++;
    if (nPrnState){
      for (j = 0; j < neq[0]; j++){
	if (!rmState[j]){
	  SET_STRING_ELT(sexp_colnames2, jj, STRING_ELT(stateNames2,j));
	  SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
	  jj++;kk++;
	}
      }
    }
    // Put in Cov names
    SEXP paramNames2 = PROTECT(rxParamNames(op->modNamePtr)); pro++;
    int *par_cov = op->par_cov;
    for (i = 0; i < ncov*add_cov; i++){
      SET_STRING_ELT(sexp_colnames2,jj, STRING_ELT(paramNames2, par_cov[i]-1));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
    }
    par_cov = rx->cov0;
    for (i = 0; i < ncov0*add_cov; i++){
      SET_STRING_ELT(sexp_colnames2,jj, STRING_ELT(paramNames2, par_cov[i]));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
    }
    for (i = 0; i < nkeep0; i++){
      SET_STRING_ELT(sexp_colnames2,jj, STRING_ELT(ikeepNames, i));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
    }
    for (i = 0; i < nkeep; i++){
      SET_STRING_ELT(sexp_colnames2,jj, STRING_ELT(fkeepNames, i));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++; kk++;
    }
    if (doTBS){
      SET_STRING_ELT(sexp_colnames2, jj, mkChar("rxLambda"));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
      SET_STRING_ELT(sexp_colnames2, jj, mkChar("rxYj"));
      SET_VECTOR_ELT(df2, jj, VECTOR_ELT(df, kk));
      jj++;kk++;
    }
    setAttrib(df2, R_NamesSymbol, sexp_colnames2);
    setAttrib(df2, R_RowNamesSymbol, sexp_rownames);
  } else {
    df2=df;
  }
  UNPROTECT(pro);
  return df2;
}


// rxSolveOldC
extern void rxSingleSolve(int subid, double *_theta, double *timep,
			  int *evidp, int *ntime,
			  double *initsp, double *dosep,
			  double *ii, double *retp,
			  double *lhsp, int *rc,
			  double *newTime, int *newEvid,
			  int *on, int *ix,
			  int *slvr_counter, int *dadt_counter, int *jac_counter,
			  double *InfusionRate, int *BadDose, int *idose,
			  double *scale, int *stateIgnore, double *mtime){
  double *theta = get_theta(_theta);
  rx_solve *rx = &rx_global;
  rx_solving_options *op = &op_global;
  rx_solving_options_ind *ind = &inds_global[subid];
  int i;
  ind->InfusionRate = InfusionRate;
  // Counters
  ind->slvr_counter = slvr_counter;
  ind->dadt_counter = dadt_counter;
  ind->jac_counter = jac_counter;

  ind->InfusionRate = InfusionRate;

  ind->BadDose = BadDose;
  ind->nBadDose = 0;

  ind->par_ptr = theta;
  ind->dose    = dosep;
  ind->ii      = ii;
  ind->solve   = retp;
  ind->lhs     = lhsp;
  ind->evid    = evidp;
  ind->rc      = rc;
  ind->n_all_times       = *ntime;
  ind->on = on;
  ind->ix = ix;
  ind->ixds = 0;
  ind->ndoses = -1;
  ind->all_times = timep;
  ind->idose = idose;
  ind->id = subid;
  ind->sim = 0;
  ind->ndoses=0;
  for (unsigned int i = 0; i < ind->n_all_times; i++){
    if (isDose(ind->evid[i])){
      ind->ndoses++;
      ind->idose[ind->ndoses-1] = i;
    }
  }
  op->badSolve=0;
  // No covariates not needed.
  // Linear is setup.
  op->ncov = 0;
  op->do_par_cov=0;
  //
  op->inits   = initsp;
  op->scale = scale;
  op->extraCmt = 0;
  rx->nsub =1;
  rx->nsim =1;
  rx->stateIgnore = stateIgnore;//gsiVSetup(op->neq);
  rx->nobs =-1;
  rx->add_cov =0;
  rx->matrix =0;
  ind->mtime = mtime;
  // Solve without the option of updating residuals.
  ind_solve(rx, subid, dydt_liblsoda, dydt_lsoda_dum, jdum_lsoda,
	      dydt, update_inis, global_jt);
  if (op->nlhs) {
    ind->_newind=1;
    for (i=0; i<*ntime; i++){
      ind->idx = i;
      newEvid[i] = ind->evid[ind->ix[i]];
      newTime[i] = getTime(ind->ix[i], ind);
      if (ind->evid[ind->ix[i]]) ind->tlast = newTime[i];
      // 0 = first subject; Calc lhs changed...
      calc_lhs(subid, newTime[i], retp+i*(op->neq), lhsp+i*(op->nlhs));
      ind->_newind=2;
    }
  }
}
