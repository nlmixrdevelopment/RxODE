#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>
#include "RxODE.h"
#include "dop853.h"
#include "common.h"
#include "lsoda.h"
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


rx_solve rx_global;

rx_solving_options op_global;

rx_solving_options_ind *inds_global;
int max_inds_global = 0;

void par_flush_console() {
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
  R_FlushConsole();
#endif
}

int isRstudio();

int par_progress(int c, int n, int d, int cores, clock_t t0, int stop){
  float progress = (float)(c)/((float)(n));
  if (c <= n){
    int nticks= (int)(progress * 50);
    int curTicks = d;
    if (isRstudio()){
      if (nticks > curTicks){
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
        for (i = nticks; i < 50; i++){
          Rprintf(" ");
        }
        Rprintf("] ");
        if (nticks < 50) Rprintf(" ");
        if (cores > 1){
          Rprintf("%02.f%%; ncores=%d; ",100*progress,cores);
        } else {
          Rprintf("%02.f%%; ",100*progress,cores);
        }
        clock_t t = clock() - t0;
        Rprintf(" %.3f sec ", ((double)t)/CLOCKS_PER_SEC);
        if (stop){
          Rprintf("Stopped Calculation!\n");
        }
      }
    } else {
      if (nticks > curTicks){
        REprintf("\r");
        int i;
        for (i = 0; i < nticks; i++){
          if (i == 0) {
            REprintf("%%[");
          } else if (i % 5 == 0) {
            REprintf("|");
          } else {
            REprintf("=");
          }
        }
        for (i = nticks; i < 50; i++){
          REprintf(" ");
        }
        REprintf("] ");
        if (nticks < 50) Rprintf(" ");
        if (cores > 1){
          REprintf("%02.f%%; ncores=%d; ",100*progress,cores);
        } else {
          REprintf("%02.f%%; ",100*progress,cores);
        }
        clock_t t = clock() - t0;
        REprintf(" %.3f sec ", ((double)t)/CLOCKS_PER_SEC);
        if (stop){
          REprintf("Stopped Calculation!\n");
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
  int clearB = INTEGER(clear)[0];
  if (clearB){
    if (isRstudio()){
      /* Rprintf("\r                                                                                \r"); */
      Rprintf("\n");
    } else {
      REprintf("\r                                                                                \r");
    }
  } else {
    par_progress(rxt.n, rxt.n, rxt.d, rxt.cores, rxt.t0, 1);
  }
  rxt.d = rxt.n;
  rxt.cur = rxt.n;
  return R_NilValue;
}

SEXP _rxProgressAbort(){
  if (rxt.d != rxt.n || rxt.cur != rxt.n){
    par_progress(rxt.n, rxt.n, rxt.d, rxt.cores, rxt.t0, 0);
    error("Aborted calculation");
  }
  return R_NilValue;
}

rx_solving_options_ind *rxOptionsIniEnsure(int mx){
  if (mx >= max_inds_global){
    Free(inds_global);
    inds_global =Calloc(mx+1024, rx_solving_options_ind);
    max_inds_global = mx+1024;
    rx_global.subjects = inds_global;
  }
  return inds_global;
}

t_dydt dydt = NULL;

t_calc_jac calc_jac = NULL;

t_calc_lhs calc_lhs = NULL;

t_update_inis update_inis = NULL;

t_dydt_lsoda_dum dydt_lsoda_dum = NULL;

t_dydt_liblsoda dydt_liblsoda = NULL;

t_jdum_lsoda jdum_lsoda = NULL;

t_set_solve set_solve = NULL;

t_get_solve get_solve = NULL;

int global_jt = 2;
int global_mf = 22;  
int global_debug = 0;

void rxUpdateFuns(SEXP trans){
  const char *lib, *s_dydt, *s_calc_jac, *s_calc_lhs, *s_inis, *s_dydt_lsoda_dum, *s_dydt_jdum_lsoda, 
    *s_ode_solver_solvedata, *s_ode_solver_get_solvedata, *s_dydt_liblsoda;
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

extern rx_solve *getRxSolve_(){
  set_solve(&rx_global);
  return &rx_global;
}

rx_solving_options *getRxOp(rx_solve *rx){
  /* if(!R_ExternalPtrAddr(rx->op)){ */
  /*   error("Cannot get global ode solver options."); */
  /* } */
  /* return (rx_solving_options*)(R_ExternalPtrAddr(rx->op)); */
  /* return &op_global; */
  return rx->op;
}

rx_solving_options_ind *getRxId(rx_solve *rx, unsigned int id){
  return &(rx->subjects[id]);
}

int handle_evid(int evid, int neq, 
		       int *BadDose,
		       double *InfusionRate,
		       double *dose,
		       double *yp,
		       int do_transit_abs,
		       double xout,
		       rx_solving_options_ind *ind){
  int wh = evid, wh100, cmt, foundBad, j;
  if (wh) {
    wh100 = floor(wh/1e5);
    wh = wh- wh100*1e5;
    cmt = (wh%10000)/100 - 1 + 100*wh100;
    if (cmt<0) {
      error("Supplied an invalid EVID (EVID=%d)", evid);
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
      if (wh>10000) {
	InfusionRate[cmt] += dose[ind->ixds];
      } else {
	if (do_transit_abs) {
	  ind->podo = dose[ind->ixds];
	  ind->tlast = xout;
	} else {
	  ind->podo = 0;
	  ind->tlast = xout;
	  yp[cmt] += dose[ind->ixds];     //dosing before obs
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
  nx = ind->n_all_times;
  evid = ind->evid;
  BadDose = ind->BadDose;
  InfusionRate = ind->InfusionRate;
  dose = ind->dose;
  ret = ind->solve;
  x = ind->all_times;
  rc= ind->rc;
  double xp = x[0];
  //--- inits the system
  /* Rprintf("ID: %d; nsim: %d; nsub: %d\n",neq[1], nsim, nsub); */
  /* Rprintf("inits[0]: %f\n",inits[0]); */
  /* Rprintf("ret[0]: %f\n",ret[0]); */
  memcpy(ret,inits, neq[0]*sizeof(double));
  u_inis(neq[1], ret); // Update initial conditions
  /* for(i=0; i<neq[0]; i++) yp[i] = inits[i]; */
  for(i=0; i<nx; i++) {
    xout = x[i];
    yp = ret+neq[0]*i;
    if(xout-xp > DBL_EPSILON*max(fabs(xout),fabs(xp))){
      ind->idx=i;
      lsoda(&ctx, yp, &xp, xout);
      if (ctx.state <= 0) {
        /* REprintf("IDID=%d, %s\n", istate, err_msg[-istate-1]); */
        *rc = ctx.state;
        // Bad Solve => NA
        for (unsigned int j = neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
        op->badSolve = 1;
        i = nx+42; // Get out of here!
      } else {
	if (R_FINITE(rx->stateTrim)){
	  double top=fabs(rx->stateTrim);
	  for (unsigned int j = neq[0]; j--;) yp[j]= max(-top, min(top,yp[j]));
	}
      }
    }
    if (!op->badSolve){
      if (handle_evid(evid[i], neq[0], BadDose, InfusionRate, dose, yp,
                      op->do_transit_abs, xout, ind)){
        ctx.state = 1;
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
  opt.mxhnil = 0;
  opt.mxordn = op->MXORDN;
  opt.mxords = op->MXORDS;
  opt.h0 = op->H0;
  opt.hmax = op->hmax2;
  opt.hmin = op->HMIN;
  opt.hmxi = 0.0;
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
  opt.mxhnil = 0;
  opt.mxordn = op->MXORDN;
  opt.mxords = op->MXORDS;
  opt.h0 = op->H0;
  opt.hmax = op->hmax2;
  opt.hmin = op->HMIN;
  opt.hmxi = 0.0;
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
    if (isRstudio()){
      /* Rprintf("\r                                                                                \r"); */
      Rprintf("\n");
    } else {
      REprintf("\r                                                                                \r");
    }
  }
}


double *global_rworkp;
unsigned int global_rworki = 0;
double *global_rwork(unsigned int mx){ 
  if (mx >= global_rworki){
    global_rworki = mx+1024;
    global_rworkp = Realloc(global_rworkp, global_rworki, double);
  }
  return global_rworkp;
}


int *global_iworkp;
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
  max_inds_global = 1024;
  inds_global =Calloc(1024, rx_solving_options_ind);

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

  max_inds_global = 0;
  Free(inds_global);

  global_InfusionRatei = 0;
  Free(global_InfusionRatep);

  global_BadDosei = 0;
  Free(global_BadDosep);

  global_scalei = 0;
  Free(global_scalep);
}

extern void ind_lsoda0(rx_solve *rx, rx_solving_options *op, int solveid, int *neq, double *rwork, int lrw, int *iwork, int liw, int jt,
                       t_dydt_lsoda_dum dydt_lsoda,
                       t_update_inis u_inis,
                       t_jdum_lsoda jdum){
  rx_solving_options_ind *ind;
  double *yp;
  
  static char *err_msg[]=
    {
      "excess work done on this call (perhaps wrong jt).",
      "excess accuracy requested (tolerances too small).",
      "illegal input detected (see printed message).",
      "repeated error test failures (check all inputs).",
      "repeated convergence failures (perhaps bad jacobian supplied or wrong choice of jt or tolerances).",
      "error weight became zero during problem. (solution component i vanished, and atol or atol(i) = 0.)",
      "work space insufficient to finish (see messages)."
    };
  int itask = 1, istate = 1, iopt = 1, i = 0, itol = 1;

  /* memset(rwork,0.0,lrw+1); */ // Does not work since it is a double
  for (i = lrw+1; i--;) rwork[i]=0;
  memset(iwork,0,liw+1); // Works because it is a integer

  neq[1] = solveid;
  
  ind = &(rx->subjects[neq[1]]);

  rwork[4] = op->H0; // H0
  rwork[5] = ind->HMAX; // Hmax
  rwork[6] = op->HMIN; // Hmin
  
  iwork[4] = 0; // ixpr
  iwork[5] = op->mxstep; // mxstep 
  iwork[6] = 0; // MXHNIL 
  iwork[7] = op->MXORDN; // MXORDN 
  iwork[8] = op->MXORDS;  // MXORDS
    
  ind->ixds = 0;
  double xp = ind->all_times[0];
  double xout;

  //--- inits the system
  memcpy(ind->solve, op->inits, neq[0]*sizeof(double));
  u_inis(neq[1], ind->solve); // Update initial conditions
  
  for(i=0; i < ind->n_all_times; i++) {
    xout = ind->all_times[i];
    yp   = ind->solve+neq[0]*i;
    if(xout - xp > DBL_EPSILON*max(fabs(xout),fabs(xp)))
      {
        ind->idx=i;
        F77_CALL(dlsoda)(dydt_lsoda, neq, yp, &xp, &xout, &itol, &(op->RTOL), &(op->ATOL), &itask,
                         &istate, &iopt, rwork, &lrw, iwork, &liw, jdum, &jt);

        if (istate <= 0) {
	  REprintf("IDID=%d, %s\n", istate, err_msg[-istate-1]);
	  ind->rc[0] = istate;
	  // Bad Solve => NA
	  for (unsigned int j=neq[0]*(ind->n_all_times); j--;) ind->solve[j] = NA_REAL;
	  op->badSolve = 1;
	  i = ind->n_all_times+42; // Get out of here!
	} else {
	  if (R_FINITE(rx->stateTrim)){
	    double top=fabs(rx->stateTrim);
	    for (unsigned int j = neq[0]; j--;) yp[j]= max(-top, min(top,yp[j]));
	  }
	}
        ind->slvr_counter[0]++;
        //dadt_counter = 0;
      }
    if (!op->badSolve){
      if (handle_evid(ind->evid[i], neq[0], ind->BadDose, ind->InfusionRate, ind->dose, yp,
                      op->do_transit_abs, xout, ind)){
        istate = 1;
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
    REprintf("JT: %d\n",cjt);
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
    REprintf("JT: %d\n",jt);
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

//dummy solout fn
void solout(long int nr, double t_old, double t, double *y, int *nptr, int *irtrn){}

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
  memcpy(ret,inits, neq[0]*sizeof(double));
  u_inis(neq[1], ret); // Update initial conditions
  //--- inits the system
  for(i=0; i<nx; i++) {
    xout = x[i];
    yp = &ret[neq[0]*i];
    if (global_debug){
      REprintf("i=%d xp=%f xout=%f\n", i, xp, xout);
    }
    if(xout-xp>DBL_EPSILON*max(fabs(xout),fabs(xp)))
      {
        ind->idx=i;
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
                      0,            /* maximal number of allowed steps */
                      1,            /* switch for the choice of the coefficients */
                      -1,                     /* test for stiffness */
                      0,                      /* number of components for which dense outpout is required */
                      NULL,           /* indexes of components for which dense output is required, >= nrdens */
                      0                       /* declared length of icon */
                      );
        if (idid<0) {
            REprintf("IDID=%d, %s\n", idid, err_msg[-idid-1]);
            *rc = idid;
            // Bad Solve => NA
            for (unsigned int j = (ind->n_all_times)*neq[0];j--;) ret[i] = NA_REAL; 
            op->badSolve = 1;
            i = nx+42; // Get out of here!
	} else {
	  if (R_FINITE(rx->stateTrim)){
	    double top=fabs(rx->stateTrim);
	    for (unsigned int j = neq[0]; j--;) yp[j]= max(-top, min(top,yp[j]));
	  }
	}
        xp = xRead();
        ind->slvr_counter[0]++;
        //dadt_counter = 0;
      }
    if (!op->badSolve){
      if (handle_evid(evid[i], neq[0], BadDose, InfusionRate, dose, yp,
                      op->do_transit_abs, xout, ind)){
        xp = xout;
      }
      /* for(j=0; j<neq[0]; j++) ret[neq[0]*i+j] = yp[j]; */
      if (i+1 != nx) memcpy(ret+neq[0]*(i+1), ret + neq[0]*i, neq[0]*sizeof(double));
    }
    //REprintf("wh=%d cmt=%d tm=%g rate=%g\n", wh, cmt, xp, InfusionRate[cmt]);

    /* if (global_debug){ */
    /*   REprintf("IDID=%d, ", idid); */
    /*   for(j=0; j<neq[0]; j++) */
    /*     { */
    /*       REprintf("%f ", yp[j]); */
    /*     } */
    /*   REprintf("\n"); */
    /* } */
    /* if (rc[0]){ */
    /*   REprintf("Error sovling using dop853\n"); */
    /*   return; */
    /* } */
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
    if (isRstudio()){
      /* Rprintf("\r                                                                                \r"); */
      Rprintf("\n");
    } else {
      REprintf("\r                                                                                \r");
    }
  }
}

void ind_solve(rx_solve *rx, unsigned int cid,
	       t_dydt_liblsoda dydt_lls,
	       t_dydt_lsoda_dum dydt_lsoda, t_jdum_lsoda jdum,
	       t_dydt c_dydt, t_update_inis u_inis,
	       int jt){
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
  rx_solving_options *op = &op_global;
  if (op->neq != 0){
    if (op->stiff == 2){
      par_liblsoda(rx);
    } else if (op->stiff == 1){
      // lsoda
      par_lsoda(rx);
    } else if (op->stiff == 0){
      // dop
      par_dop(rx);
    }
  }
}

rx_solve *_globalRx = NULL;

extern void rxode_assign_rx(rx_solve *rx){
  _globalRx=rx;
}

extern int rxEvidP(int i, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  if (i < ind->n_all_times){
    return(ind->evid[i]); // invalid read of size 4
  } else {
    error("Trying to access EVID outside of defined events.\n");
  }
}
extern unsigned int nDosesP(rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  if (ind->ndoses < 0){
    ind->ndoses=0;
    for (int i = 0; i < ind->n_all_times; i++){
      if (rxEvidP(i, rx, id)){
        ind->ndoses++;
        ind->idose[ind->ndoses-1] = i;
      }
    }
    return ind->ndoses;
  } else {
    return ind->ndoses;
  }
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
    calc_lhs((int)id, ind->all_times[i], solve+i*op->neq, lhs);
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

double rxDoseP(int i, rx_solve *rx, unsigned int id){
  if ((unsigned int)i < nDosesP(rx, id)){
    rx_solving_options_ind *ind;
    ind = getRxId(rx, id);
    return(ind->dose[i]);
  } else {
    error("Dose cannot be retrived (%dth entry).",i);
  }
  return 0;
}
double rxDose(int i){
  return(rxDoseP(i, _globalRx, 0));
}


SEXP rxStateNames(char *ptr);
SEXP rxLhsNames(char *ptr);
SEXP rxParamNames(char *ptr);

extern double *rxGetErrs();
extern int rxGetErrsNcol();

extern SEXP RxODE_df(int doDose, int doTBS){
  rx_solve *rx;
  rx = &rx_global;
  rx_solving_options *op = &op_global;
  int add_cov = rx->add_cov;
  int ncov = op->ncov;
  int nlhs = op->nlhs;
  int nobs = rx->nobs;
  int nsim = rx->nsim;
  int nall = rx->nall;
  rx->nr = (doDose == 1 ? nall : nobs)*nsim;
  int *rmState = rx->stateIgnore;
  int nPrnState =0;
  int i, j;
  int neq[2];
  double *scale;
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
  int ncols =add_cov*ncov+1+nPrnState+nlhs;
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
  int csub = 0, evid;
  int nsub = rx->nsub;
  rx_solving_options_ind *ind;
  SEXP df = PROTECT(allocVector(VECSXP,ncols+nidCols+doseCols+doTBS*2)); pro++;
  for (i = nidCols; i--;){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(INTSXP, rx->nr))); pro++;
  }
  i = nidCols;
  double *par_ptr;
  double *errs = rxGetErrs();
  
  int updateErr = 0;
  int errNcol = rxGetErrsNcol();
  if (errNcol > 0){
    updateErr = 1;
  }
  if (doDose){
    SET_VECTOR_ELT(df, i++, PROTECT(allocVector(INTSXP, rx->nr))); pro++;
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(REALSXP, rx->nr))); pro++;
  }
  for (i = md + sm + doseCols; i < ncols + doseCols+nidCols; i++){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(REALSXP, rx->nr))); pro++;
  }
  for (i = ncols + doseCols+nidCols; i < ncols + doseCols+nidCols + doTBS*2; i++){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(REALSXP, rx->nr))); pro++;
  }
  
  // Now create the data frame
  double *dfp;
  int *dfi;
  int ii=0, jj = 0, ntimes;
  double *solve;
  double *cov_ptr;
  int nBadDose;
  int *BadDose;
  int extraCmt = op->extraCmt;
  double *dose;
  int *svar = op->svar;
  int di = 0;
  int kk = 0;
  for (int csim = 0; csim < nsim; csim++){
    for (csub = 0; csub < nsub; csub++){
      neq[1] = csub+csim*nsub;
      ind = &(rx->subjects[neq[1]]);
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
      for (i = 0; i < ntimes; i++){
        evid = ind->evid[i];
        if (updateErr){
          for (j=0; j < errNcol; j++){
	    par_ptr[svar[j]] = errs[rx->nr*j+kk];
          }
	  if (evid == 0 || doDose){
	    // Only incerement if this is an observation or of this a
	    // simulation that requests dosing infomration too.
            kk++;
	  }
        }
        jj  = 0 ;
	if (evid==0 || doDose){
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
	    // evid
            dfi = INTEGER(VECTOR_ELT(df, jj++));
            dfi[ii] = evid;
            // amt
            dfp = REAL(VECTOR_ELT(df, jj++));
            dfp[ii] = (evid == 0 ? NA_REAL : dose[di++]);
	  }
          // time
          dfp = REAL(VECTOR_ELT(df, jj++));
          dfp[ii] = ind->all_times[i];
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
	      dfp[ii] = (evid == 0 ? cov_ptr[j*ntimes+i] : NA_REAL);
	      jj++;
	    }
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
  SEXP sexp_colnames = PROTECT(allocVector(STRSXP,ncols+nidCols+doseCols+doTBS*2)); pro++;
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
    SET_STRING_ELT(sexp_colnames, jj, mkChar("amt"));
    jj++;
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
  if (doTBS){
    SET_STRING_ELT(sexp_colnames, jj, mkChar("rxLambda"));
    jj++;
    SET_STRING_ELT(sexp_colnames, jj, mkChar("rxYj"));
    jj++;
  }
  setAttrib(df, R_NamesSymbol, sexp_colnames);
  UNPROTECT(pro);
  return df;
}

int *gidoseSetup(int n);
int *gsiVSetup(int n);
int *gslvr_counterSetup(int n);
int *gdadt_counterSetup(int n);
int *gjac_counterSetup(int n);
// rxSolveOldC
void protectOld();
extern void rxSolveOldC(int *neqa,
                        double *theta,  //order:
                        double *timep,
                        int *evidp,
                        int *ntime,
                        double *initsp,
                        double *dosep,
                        double *retp,
                        double *atol,
                        double *rtol,
                        int *stiffa,
                        int *transit_abs,
                        int *nlhsa,
                        double *lhsp,
                        int *rc){
  protectOld();
  rx_solve *rx = &rx_global;
  rx_solving_options *op = &op_global;
  rx_solving_options_ind *ind = &inds_global[0];
  int i;
  rx->stateTrim=R_PosInf;
  // Counters
  ind->slvr_counter = gslvr_counterSetup(1);
  ind->dadt_counter = gdadt_counterSetup(1);
  ind->jac_counter = gjac_counterSetup(1);
  ind->slvr_counter[0]   = 0;
  ind->dadt_counter[0]   = 0;
  ind->jac_counter[0]   = 0;

  ind->InfusionRate = global_InfusionRate(*neqa);
  /* memset(ind->InfusionRate, 0.0, *neqa);  not for doubles*/
  for (unsigned int j = *neqa; j--;) ind->InfusionRate[j]=0.0;
  
  ind->BadDose = global_BadDose(*neqa);
  memset(ind->BadDose, 0, *neqa); // int ok
  ind->nBadDose = 0;

  ind->HMAX = 0;
  ind->tlast = 0.0;
  ind->podo = 0;
  ind->par_ptr = theta;
  ind->dose    = dosep;
  ind->solve   = retp;
  ind->lhs     = lhsp;
  ind->evid    = evidp;
  ind->rc      = rc;
  /* double *cov_ptr; */
  /* ind->cov_ptr = cov_ptr; */
  ind->n_all_times       = *ntime;
  ind->ixds = 0;
  ind->ndoses = -1;
  ind->all_times = timep;
  ind->idose = gidoseSetup(*ntime);
  ind->id = -1;
  ind->sim = -1;
  
  op->badSolve=0;
  op->ATOL = *atol;
  op->RTOL = *rtol;
  op->H0 = 0;
  op->HMIN = 0;
  op->mxstep = 5000; // Not LSODA default but RxODE default
  op->MXORDN         = 0;
  op->MXORDS         = 0;
  op->do_transit_abs = *transit_abs;
  op->nlhs           = *nlhsa;
  op->neq            = *neqa;
  op->stiff          = *stiffa;
  // No covariates not needed.
  // Linear is setup.
  op->f1 = 1.0;
  op->f2 = 0.0;
  op->kind = 1;
  op->is_locf = 0;
  op->ncov = 0;
  op->do_par_cov=0;
  //
  op->inits   = initsp;
  op->scale = global_scale(*neqa);
  /* memset(op->scale, 1.0, *neqa); */
  for (unsigned int j = *neqa; j--;) op->scale[j] = 1.0;
  op->extraCmt = 0;
  op->hmax2=0;
  /* double *rtol2, *atol2; */
  /* op->rtol2 = rtol2; */
  /* op->atol2 = atol2; */
  op->cores = 1;
  op->nDisplayProgress = 100;
  op->ncoresRV = 1;
  op->isChol = 0;
  /* int *svar; */
  /* op->svar = svar; */
  op->abort = 0;  
  // FIXME? modNamePtr?
  /* op->modNamePtr */
  rx->subjects = ind;
  rx->nsub =1;
  rx->nsim =1;
  rx->stateIgnore = gsiVSetup(*neqa);
  memset(rx->stateIgnore, 0, *neqa); // int OK
  rx->nobs =-1;
  rx->add_cov =0;
  rx->matrix =0;
  /* int i =0; */
  _globalRx=rx;
  rx->op = &op_global;
  /* rxode_assign_rx(rx); */
  set_solve(rx);
  par_solve(rx); // Solve without the option of updating residuals.
  if (*nlhsa) {
    for (i=0; i<*ntime; i++){
      // 0 = first subject; Calc lhs changed...
      ind->idx = i;
      calc_lhs(0, timep[i], retp+i*(*neqa), lhsp+i*(*nlhsa));
    }
  }
}
