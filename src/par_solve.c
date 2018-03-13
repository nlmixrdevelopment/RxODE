#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>
#include "solve.h"
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

int par_progress(int c, int n, int d, int cores, clock_t t0, int stop){
  float progress = (float)(c)/((float)(n));
  if (c <= n){
    int nticks= (int)(progress * 50);
    int curTicks = d;
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
      if (nticks < 50) REprintf(" ");
      REprintf("%02.f%%; ncores=%d; ",100*progress,cores);
      clock_t t = clock() - t0;
      REprintf(" %.3f sec ", ((double)t)/CLOCKS_PER_SEC);
      if (stop){
	REprintf("Stopped Calculation!\n");
      }
      if (nticks >= 50){
	REprintf("\n");
      }
    }
    par_flush_console();
    return nticks;
  }
  return d;
}

rx_solving_options_ind *rxOptionsIniEnsure(int mx){
  if (mx >= max_inds_global){
    max_inds_global = mx+1024;
    inds_global = Realloc(inds_global, max_inds_global, rx_solving_options_ind);
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
  global_jt = 2;
  global_mf = 22;  
  global_debug = 0;
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

inline rx_solving_options_ind *getRxId(rx_solve *rx, unsigned int id){
  return &(rx->subjects[id]);
}

inline int handle_evid(int evid, int neq, 
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


extern void par_liblsoda(rx_solve *rx){
  rx_solving_options *op = &op_global;
#ifdef _OPENMP
  int cores = op->cores;
#else
  int cores = 1;
#endif
  int nsub = rx->nsub, nsim = rx->nsim;
  int displayProgress = (op->nDisplayProgress <= nsim*nsub);
  clock_t t0;
  if (displayProgress) t0 = clock();
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
      /* memset(ret + neq[0],0.0, (nx-1)*neq[0]); */
      memcpy(ret,inits, neq[0]*sizeof(double));
      update_inis(neq[1], ret); // Update initial conditions
      /* for(i=0; i<neq[0]; i++) yp[i] = inits[i]; */
      for(i=0; i<nx; i++) {
	xout = x[i];
        yp = ret+neq[0]*i;
	if(xout-xp > DBL_EPSILON*max(fabs(xout),fabs(xp))){
	  lsoda(&ctx, yp, &xp, xout);
	  if (ctx.state <= 0) {
	    /* REprintf("IDID=%d, %s\n", istate, err_msg[-istate-1]); */
	    *rc = ctx.state;
	    // Bad Solve => NA
            memset(ret,NA_REAL, nx*neq[0]);
	    /* for (i = 0; i < nx*neq[0]; i++) ret[i] = NA_REAL; */
	    op->badSolve = 1;
	    i = nx+42; // Get out of here!
	  }
	}
	if (handle_evid(evid[i], neq[0], BadDose, InfusionRate, dose, yp,
			op->do_transit_abs, xout, ind)){
	  ctx.state = 1;
	  xp = xout;
	}
	if (i+1 != nx) memcpy(ret+neq[0]*(i+1), ret + neq[0]*i, neq[0]*sizeof(double));
	/* for(j=0; j<neq[0]; j++) ret[neq[0]*i+j] = yp[j]; */
      }
      lsoda_free(&ctx);
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
}


double *global_rworkp;
unsigned int global_rworki = 0;
inline double *global_rwork(unsigned int mx){ 
  if (mx >= global_rworki){
    global_rworki = mx+1024;
    global_rworkp = Realloc(global_rworkp, global_rworki, double);
  }
  return global_rworkp;
}


int *global_iworkp;
unsigned int global_iworki = 0;
inline int *global_iwork(unsigned int mx){
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
inline double *global_scale(unsigned int mx){
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
  inds_global =Calloc(1024, rx_solving_options_ind);
  global_iworkp=Calloc(1024*4, int);
  global_rworkp=Calloc(1024*4, double);
  global_iworki=4*1024;
  global_rworki=4*1024;
  max_inds_global = 1024;
  global_InfusionRatep=Calloc(1024, double);
  global_InfusionRatei = 1024;
  global_BadDosep=Calloc(1024, int);
  global_BadDosei = 1024;
  global_scalep=Calloc(1024, double);
  global_scalei = 1024;
  rx_solve *rx=(&rx_global);
  rx->op = &op_global;
  rx->subjects = inds_global;
}

void rxOptionsFree(){
  Free(global_iworkp);
  Free(global_rworkp);
  Free(inds_global);
  Free(global_InfusionRatep);
  Free(global_BadDosep);
}

extern void par_lsoda(rx_solve *rx){
  rx_solving_options *op = &op_global;
  int nsub = rx->nsub, nsim = rx->nsim;
  int displayProgress = (op->nDisplayProgress <= nsim*nsub);
  clock_t t0 = NULL;
  if (displayProgress)
    t0 = clock();
  int i;
  double xout;
  double *yp;
  int neq[2];
  neq[0] = op->neq;
  neq[1] = 0;
  /* yp = global_yp(neq[0]); */
  int itol = 1;
  double  rtol = op->RTOL, atol = op->ATOL;
  // Set jt to 1 if full is specified.
  int itask = 1, istate = 1, iopt = 0, lrw=22+neq[0]*max(16, neq[0]+9), liw=20+neq[0], jt = global_jt;
  double *rwork;
  int *iwork;
  
  char *err_msg[]=
    {
      "excess work done on this call (perhaps wrong jt).",
      "excess accuracy requested (tolerances too small).",
      "illegal input detected (see printed message).",
      "repeated error test failures (check all inputs).",
      "repeated convergence failures (perhaps bad jacobian supplied or wrong choice of jt or tolerances).",
      "error weight became zero during problem. (solution component i vanished, and atol or atol(i) = 0.)",
      "work space insufficient to finish (see messages)."
    };
  if (global_debug)
    REprintf("JT: %d\n",jt);
  rwork = global_rwork(lrw+1);
  iwork = global_iwork(liw+1);
  
  /* iopt = 1; */
  
  rx_solving_options_ind *ind;
  /* int cores = op->cores; */
  
  int curTick = 0;
  int abort = 0;
  for (int solveid = 0; solveid < nsim*nsub; solveid++){
    itask = 1; 
    istate = 1;
    iopt = 1;
    memset(rwork,0.0,lrw+1);
    /* for (i = 0; i < lrw+1; i++) rwork[i]=0; */
    memset(iwork,0,liw+1);
    /* for (i = 0; i < liw+1; i++) iwork[i]=0; */
    /* for (i = 0; i < neq[0]; i++) yp[i]=0; */
    rwork[4] = op->H0; // H0 -- determined by solver
    rwork[6] = op->HMIN; // Hmin -- 0
  
    iwork[4] = 0; // ixpr  -- No extra printing.
    iwork[5] = op->mxstep; // mxstep 
    iwork[6] = 0; // MXHNIL 
    iwork[7] = op->MXORDN; // MXORDN 
    iwork[8] = op->MXORDS;  // MXORDS
    neq[1] = solveid;
    ind = &(rx->subjects[neq[1]]);
    ind->ixds = 0;
    rwork[5] = ind->HMAX; // Hmax -- Infinite
    double xp = ind->all_times[0];
    //--- inits the system
    memcpy(ind->solve,op->inits, neq[0]*sizeof(double));
    update_inis(neq[1], ind->solve); // Update initial conditions
    /* for(i=0; i<neq[0]; i++) yp[i] = inits[i]; */
    /* memcpy(yp,inits, neq[0]*sizeof(double)); */
    for(i=0; i<ind->n_all_times; i++) {
      xout = ind->all_times[i];
      yp = &ind->solve[neq[0]*i];
      memset(yp,0.0, neq[0]);
      /* if (global_debug){ */
      /*   REprintf("i=%d xp=%f xout=%f\n", i, xp, xout); */
      /* } */
      if(xout-xp > DBL_EPSILON*max(fabs(xout),fabs(xp)))
	{
	  F77_CALL(dlsoda)(dydt_lsoda_dum, neq, yp, &xp, &xout, &itol, &rtol, &atol, &itask,
			   &istate, &iopt, rwork, &lrw, iwork, &liw, jdum_lsoda, &jt);

	  if (istate<0)
	    {
	      REprintf("IDID=%d, %s\n", istate, err_msg[-istate-1]);
	      ind->rc[0] = istate;
	      // Bad Solve => NA
	      memset(ind->solve,NA_REAL, (ind->n_all_times)*neq[0]);
	      /* for (i = 0; i < nx*neq[0]; i++) ret[i] = NA_REAL; */
	      op->badSolve = 1;
	      i = ind->n_all_times+42; // Get out of here!
	    }
	  ind->slvr_counter++;
	  //dadt_counter = 0;
	}
      if (handle_evid(ind->evid[i], neq[0], ind->BadDose, ind->InfusionRate, ind->dose, yp,
		      op->do_transit_abs, xout, ind)){
	istate = 1;
	xp = xout;
      }
      if (i+1 != ind->n_all_times) memcpy(ind->solve+neq[0]*(i+1), yp, neq[0]*sizeof(double));
      /* for(j=0; j<neq[0]; j++) ret[neq[0]*i+j] = yp[j]; */
      /* memcpy(&ret[neq[0]*i],yp, neq[0]*sizeof(double)); */
      //REprintf("wh=%d cmt=%d tm=%g rate=%g\n", wh, cmt, xp, InfusionRate[cmt]);

      /* if (global_debug){ */
      /*   REprintf("ISTATE=%d, ", istate); */
      /*   for(j=0; j<neq[0]; j++) */
      /*     { */
      /*       REprintf("%f ", yp[j]); */
      /*     } */
      /*   REprintf("\n"); */
      /* } */
    }
    if (displayProgress){ // Can only abort if it is long enough to display progress.
      curTick = par_progress(solveid, nsim*nsub, curTick, 1, t0, 0);
      if (checkInterrupt()){
	abort =1;
	break;
      }
    }
  }
  if (abort == 1){
    op->abort = 1;
  } else {
    if (displayProgress && curTick < 50) par_progress(nsim*nsub, nsim*nsub, curTick, 1, t0, 0);
  }
  /* if (rc[0]){ */
  /*   /\* REprintf("Error solving using LSODA\n"); *\/ */
  /*   /\* Free(rwork); *\/ */
  /*   /\* Free(iwork); *\/ */
  /*   /\* Free(yp); *\/ */
  /*   /\* return; *\/ */
  /* } */
}

//dummy solout fn
void solout(long int nr, double t_old, double t, double *y, int *nptr, int *irtrn){}

void par_dop(rx_solve *rx){
  rx_solving_options *op = &op_global;
  int nsub = rx->nsub, nsim = rx->nsim;
  int displayProgress = (op->nDisplayProgress <= nsim*nsub);
  clock_t t0;
  if (displayProgress)
    t0 = clock();
  int i, j;
  double xout;
  double *yp;
  int neq[2];
  neq[0] = op->neq;
  neq[1] = 0;
  
  //DE solver config vars
  double rtol=op->RTOL, atol=op->ATOL;
  int itol=0;           //0: rtol/atol scalars; 1: rtol/atol vectors
  int iout=0;           //iout=0: solout() NEVER called
  int idid=0;
  char *err_msg[]=
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
  // This part CAN be parallelized, if dop is thread safe...
  // Therefore you could use https://github.com/jacobwilliams/dop853, but I haven't yet
  
  int curTick = 0;
  int abort = 0;
  for (int solveid = 0; solveid < nsim*nsub; solveid++){
    if (abort == 0){
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
      update_inis(neq[1], ret); // Update initial conditions
      //--- inits the system
      /* for(i=0; i<neq[0]; i++) yp[i] = inits[i]; */

      for(i=0; i<nx; i++) {
	xout = x[i];
        yp = &ret[neq[0]*i];
	if (global_debug){
	  REprintf("i=%d xp=%f xout=%f\n", i, xp, xout);
	}
	if(xout-xp>DBL_EPSILON*max(fabs(xout),fabs(xp)))
	  {
	    idid = dop853(neq,       /* dimension of the system <= UINT_MAX-1*/
			  dydt,       /* function computing the value of f(x,y) */
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
	    if (idid<0)
	      {
		REprintf("IDID=%d, %s\n", idid, err_msg[-idid-1]);
		*rc = idid;
		// Bad Solve => NA
		/* for (i = 0; i < nx*neq[0]; i++) ret[i] = NA_REAL; */
                memset(ret,NA_REAL, nx*neq[0]);
		op->badSolve = 1;
		i = nx+42; // Get out of here!
	      }
	    xp = xRead();
	    ind->slvr_counter++;
	    //dadt_counter = 0;
	  }
	if (handle_evid(evid[i], neq[0], BadDose, InfusionRate, dose, yp,
			op->do_transit_abs, xout, ind)){
	  xp = xout;
	}
	/* for(j=0; j<neq[0]; j++) ret[neq[0]*i+j] = yp[j]; */
        if (i+1 != nx) memcpy(ret+neq[0]*(i+1), ret + neq[0]*i, neq[0]*sizeof(double));
	//REprintf("wh=%d cmt=%d tm=%g rate=%g\n", wh, cmt, xp, InfusionRate[cmt]);

	if (global_debug){
	  REprintf("IDID=%d, ", idid);
	  for(j=0; j<neq[0]; j++)
	    {
	      REprintf("%f ", yp[j]);
	    }
	  REprintf("\n");
	}
	/* if (rc[0]){ */
	/*   REprintf("Error sovling using dop853\n"); */
	/*   return; */
	/* } */
      }
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
}

inline void par_solve(rx_solve *rx){
  rx_solving_options *op = &op_global;
  if (op->neq > 0){
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

extern int nEqP (rx_solve *rx){
  rx_solving_options *op = &op_global;
  return op->neq;
}

extern int nEq (rx_solve *rx){
  return nEqP(_globalRx);
}

extern int rxEvidP(int i, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  if (i < ind->n_all_times){
    return(ind->evid[i]);
  } else {
    error("Trying to access EVID outside of defined events.\n");
  }
}

extern int rxEvid(int i, rx_solve *rx, unsigned int id){
  return rxEvidP(i, _globalRx, 0);
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

extern unsigned int nDoses(){
  return nDosesP(_globalRx, 0);
}

extern unsigned int nObsP (rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(_globalRx, id);
  return (unsigned int)(ind->n_all_times - nDosesP(rx, id));
}

extern unsigned int nObs (){
  return nObsP(_globalRx, 0);
}

extern unsigned int nLhsP(rx_solve *rx){
  rx_solving_options *op = &op_global;
  return (unsigned int)(op->nlhs);
}
extern unsigned int nLhs(){
  return nLhsP(_globalRx);
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

extern double rxLhs(int i, rx_solve *rx, unsigned int id){
  return rxLhsP(i, _globalRx, 0);
}


extern void rxCalcLhsP(int i, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  rx_solving_options *op = &op_global;
  double *solve, *lhs;
  solve = ind->solve;
  lhs = ind->lhs;
  if (i < ind->n_all_times){
    calc_lhs((int)id, ind->all_times[i], solve+i*op->neq, lhs);
  } else {
    error("LHS cannot be calculated (%dth entry).",i);
  }
}

void rxCalcLhs(int i, rx_solve *rx, unsigned int id){
  rxCalcLhsP(i, _globalRx, 0);
}
extern unsigned int nAllTimesP(rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  return (unsigned int)(ind->n_all_times);
}

extern unsigned int nAllTimes(){
  return nAllTimesP(_globalRx, 0);
}


extern double RxODE_InfusionRateP(int val, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  return ind->InfusionRate[val];
}

extern double RxODE_InfusionRate(int val){
  return RxODE_InfusionRateP(val, _globalRx, 0);
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

double rxDosingTimeP(int i, rx_solve *rx, unsigned int id){
  if (i < nDosesP(rx,id)){
    rx_solving_options_ind *ind;
    ind = getRxId(rx, id);
    return ind->all_times[ind->idose[i]];
  } else {
    error("Dosing cannot retreived (%dth dose).", i);
  }
  return 0;
}
double rxDosingTime(int i){
  return rxDosingTimeP(i, _globalRx, 0);
}

int rxDosingEvidP(int i, rx_solve *rx, unsigned int id){
  if (i < nDosesP(rx, id)){
    rx_solving_options_ind *ind;
    ind = getRxId(rx, id);
    return (ind->evid[ind->idose[i]]);
  } else {
    error("Dosing cannot retreived (%dth dose).", i);
  }
  return 0;
}

int rxDosingEvid(int i){
  return rxDosingEvidP(i, _globalRx, 0);
}

double rxDoseP(int i, rx_solve *rx, unsigned int id){
  if (i < nDosesP(rx, id)){
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

extern SEXP RxODE_par_df(SEXP sd){
  rx_solve *rx;
  rx = &rx_global;
  rx_solving_options *op = &op_global;
  // Mutiple ID data?
  int md = 0;
  if (rx->nsub > 1) md = 1;
  // Multiple simulation data?
  int sm = 0;
  if (rx->nsim > 1) sm = 1;
  int nsub = rx->nsub;
  int nsim = rx->nsim;
  int pro=0, i;
  // paramNames
  SEXP paramNames = PROTECT(rxParamNames(op->modNamePtr)); pro++;
  int *par_cov = op->par_cov;
  int ncov = op->ncov;
  int npar = length(paramNames);
  int n = npar - ncov;
  SEXP df = PROTECT(allocVector(VECSXP,n+md+sm)); pro++;
  SEXP countsDf = PROTECT(allocVector(VECSXP,md+sm+3)); pro++;
  for (i = 0; i < md+sm; i++){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(INTSXP, nsim*nsub))); pro++;
    SET_VECTOR_ELT(countsDf, i, PROTECT(allocVector(INTSXP, nsim*nsub))); pro++;
  }
  for (i = 0; i < n; i++){
    SET_VECTOR_ELT(df, i+md+sm, PROTECT(allocVector(REALSXP, nsim*nsub))); pro++;
  }
  SET_VECTOR_ELT(countsDf, md+sm, PROTECT(allocVector(INTSXP, nsim*nsub))); pro++;
  SET_VECTOR_ELT(countsDf, md+sm+1, PROTECT(allocVector(INTSXP, nsim*nsub))); pro++;
  SET_VECTOR_ELT(countsDf, md+sm+2, PROTECT(allocVector(INTSXP, nsim*nsub))); pro++;
  int jj = 0, ii = 0, j = 0, nall = 0;
  double *par_ptr;
  rx_solving_options_ind *ind;
  int *dfi;
  double *dfp;
  int nobs = rx->nobs;
  int csub;
  int is_cov = 0;
  for (csub = 0; csub < nsub; csub++){
    ind = &(rx->subjects[csub]);
    nall+=ind->n_all_times;
  }
  // Event table information.
  SEXP dfe = PROTECT(allocVector(VECSXP,md+3)); pro++; // Events
  SEXP dfd = PROTECT(allocVector(VECSXP,md+3)); pro++; // Dosing
  SEXP dfs = PROTECT(allocVector(VECSXP,md+3)); pro++; // Sampling
  SEXP dfn = PROTECT(allocVector(STRSXP,md+3)); pro++;
  SEXP dfn1 = PROTECT(allocVector(STRSXP,md+3)); pro++;
  SEXP dfn2 = PROTECT(allocVector(STRSXP,md+3)); pro++;
  // Covariate Table
  SEXP covs = PROTECT(allocVector(VECSXP,md+ncov)); pro++;
  SEXP covsn = PROTECT(allocVector(STRSXP,md+ncov)); pro++;
  if (md){
    SET_VECTOR_ELT(dfe, 0, PROTECT(allocVector(INTSXP, nall))); pro++;
    SET_VECTOR_ELT(dfd, 0, PROTECT(allocVector(INTSXP, nall-nobs))); pro++;
    SET_VECTOR_ELT(dfs, 0, PROTECT(allocVector(INTSXP, nobs))); pro++;
    SET_VECTOR_ELT(covs, 0, PROTECT(allocVector(INTSXP, nobs))); pro++;
    SET_STRING_ELT(dfn, 0, mkChar("id"));
    SET_STRING_ELT(dfn1, 0, mkChar("id"));
    SET_STRING_ELT(dfn2, 0, mkChar("id"));
    SET_STRING_ELT(covsn, 0, mkChar("id"));
    i++;
  }
  for (i = 0; i < ncov; i++){
    SET_VECTOR_ELT(covs, md+i, PROTECT(allocVector(REALSXP, nobs))); pro++;
  }
  SEXP isEt = PROTECT(allocVector(LGLSXP,nall)); pro++;
  
  SET_VECTOR_ELT(dfe, md, PROTECT(allocVector(REALSXP, nall))); pro++;
  SET_VECTOR_ELT(dfd, md, PROTECT(allocVector(REALSXP, nall-nobs))); pro++;
  SET_VECTOR_ELT(dfs, md, PROTECT(allocVector(REALSXP, nobs))); pro++;
  SET_STRING_ELT(dfn, md, mkChar("time"));
  SET_STRING_ELT(dfn1, md, mkChar("time"));
  SET_STRING_ELT(dfn2, md, mkChar("time"));


  SET_VECTOR_ELT(dfe, md+1, PROTECT(allocVector(INTSXP, nall))); pro++;
  SET_VECTOR_ELT(dfd, md+1, PROTECT(allocVector(INTSXP, nall-nobs))); pro++;
  SET_VECTOR_ELT(dfs, md+1, PROTECT(allocVector(INTSXP, nobs))); pro++;
  SET_STRING_ELT(dfn, md+1, mkChar("evid"));
  SET_STRING_ELT(dfn1, md+1, mkChar("evid"));
  SET_STRING_ELT(dfn2, md+1, mkChar("evid"));

  SET_VECTOR_ELT(dfe, md+2, PROTECT(allocVector(REALSXP, nall))); pro++;
  SET_VECTOR_ELT(dfd, md+2, PROTECT(allocVector(REALSXP, nall-nobs))); pro++;
  SET_VECTOR_ELT(dfs, md+2, PROTECT(allocVector(REALSXP, nobs))); pro++;
  SET_STRING_ELT(dfn, md+2, mkChar("amt"));
  SET_STRING_ELT(dfn1, md+2, mkChar("amt"));
  SET_STRING_ELT(dfn2, md+2, mkChar("amt"));

  i++;

  SEXP dfre = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(dfre)[0] = NA_INTEGER;
  INTEGER(dfre)[1] = -nall;
  setAttrib(dfe, R_RowNamesSymbol, dfre);

  SEXP dfrd = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(dfrd)[0] = NA_INTEGER;
  INTEGER(dfrd)[1] = -(nall-nobs);
  setAttrib(dfd, R_RowNamesSymbol, dfrd);

  SEXP dfrs = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(dfrs)[0] = NA_INTEGER;
  INTEGER(dfrs)[1] = -nobs;
  
  SEXP dfrs1 = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(dfrs1)[0] = NA_INTEGER;
  INTEGER(dfrs1)[1] = -nobs;
  
  setAttrib(dfs, R_RowNamesSymbol, dfrs);
  setAttrib(covs, R_RowNamesSymbol, dfrs1);


  setAttrib(dfs, R_NamesSymbol, dfn);
  setAttrib(dfd, R_NamesSymbol, dfn1);
  setAttrib(dfe, R_NamesSymbol, dfn2);

  SEXP clse = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(clse, 0, mkChar("data.frame"));

  SEXP clse1 = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(clse1, 0, mkChar("data.frame"));
  
  SEXP clse2 = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(clse2, 0, mkChar("data.frame"));

  classgets(dfs, clse);
  classgets(dfd, clse1);
  classgets(dfe, clse2);
  
  double *times, *doses, **cov_ptrA, *cov_ptr;
  int evid, iie = 0, iis = 0, iid = 0, curdose, ntimes;
  int k, kk;
  for (csub = 0; csub < nsub; csub++){
    ind = &(rx->subjects[csub]);
    ntimes = ind->n_all_times;
    times = ind->all_times;
    doses = ind->dose;
    cov_ptrA = ind->cov_ptr;
    curdose=0;
    for (i = 0; i < ntimes; i++){
      evid = rxEvidP(i,rx,csub);
      if (evid ==0){
	// Sampling Record
	// id
	if (md){
	  INTEGER(VECTOR_ELT(dfe, 0))[iie] = csub+1;
	  INTEGER(VECTOR_ELT(dfs, 0))[iis] = csub+1;
	}
	// time
	REAL(VECTOR_ELT(dfe, md))[iie] = times[i];
        REAL(VECTOR_ELT(dfs, md))[iis] = times[i];
	// evid
	INTEGER(VECTOR_ELT(dfe, md+1))[iie] = evid;
        INTEGER(VECTOR_ELT(dfs, md+1))[iis] = evid;
	// amt
	REAL(VECTOR_ELT(dfe, md+2))[iie] = NA_REAL;
        REAL(VECTOR_ELT(dfs, md+2))[iis] = NA_REAL;
	LOGICAL(isEt)[iie] = 1;
	kk = 0;
	for (k = 0; k < npar; k++){
          is_cov=0;
          for (j = 0; j < ncov; j++){
            if (par_cov[j]-1 == k){
              is_cov=1;
              break;
            }
          }
          if (is_cov){
	    /* REprintf("covs[%d, %d]\n", kk,iis); */
	    dfp = REAL(VECTOR_ELT(covs, kk+md));
	    cov_ptr = cov_ptrA[kk];
            dfp[iis] = cov_ptr[iis];
            kk++;
          }
	}
	iie++;
        iis++;        
      } else {
	// Dosing Record
	if (md){
          INTEGER(VECTOR_ELT(dfe, 0))[iie] = csub+1;
          INTEGER(VECTOR_ELT(dfd, 0))[iid] = csub+1;
        }
        // time
        REAL(VECTOR_ELT(dfe, md))[iie] = times[i];
        REAL(VECTOR_ELT(dfd, md))[iid] = times[i];
        // evid
        INTEGER(VECTOR_ELT(dfe, md+1))[iie] = evid;
        INTEGER(VECTOR_ELT(dfd, md+1))[iid] = evid;
        // amt
        REAL(VECTOR_ELT(dfe, md+2))[iie] = doses[curdose];
        REAL(VECTOR_ELT(dfd, md+2))[iid] = doses[curdose];
	
	LOGICAL(isEt)[iie] = 0;
	
	curdose++;
	iie++;
        iid++;
      }
    }
  }
  for (int csim = 0; csim < nsim; csim++){
    for (csub = 0; csub < nsub; csub++){
      j = csub+csim*nsub;
      ind = &(rx->subjects[j]);
      par_ptr = ind->par_ptr;
      jj = 0;
      // sim.id
      if (sm){
        dfi = INTEGER(VECTOR_ELT(df, jj));
        dfi[ii] = csim+1;
	dfi = INTEGER(VECTOR_ELT(countsDf, jj));
        dfi[ii] = csim+1;
        jj++;
      }
      // id
      if (md){
        dfi = INTEGER(VECTOR_ELT(df, jj));
        dfi[ii] = csub+1;
	dfi = INTEGER(VECTOR_ELT(countsDf, jj));
        dfi[ii] = csub+1;
        jj++;
      }
      INTEGER(VECTOR_ELT(countsDf, sm+md))[ii] = (int)(ind->slvr_counter);
      INTEGER(VECTOR_ELT(countsDf, sm+md+1))[ii] = (int)(ind->dadt_counter);
      INTEGER(VECTOR_ELT(countsDf, sm+md+2))[ii] = (int)(ind->jac_counter);
      for (i = 0; i < npar; i++){
	is_cov=0;
	for (j = 0; j < ncov; j++){
	  if (par_cov[j]-1 == i){
	    is_cov=1;
	    break;
	  }
	}
	if (!is_cov){
	  dfp=REAL(VECTOR_ELT(df, jj));
          dfp[ii] = par_ptr[i];
          jj++;
	}
      }
      ii++;
    }
  }
  SEXP sexp_rownames = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(sexp_rownames)[0] = NA_INTEGER;
  INTEGER(sexp_rownames)[1] = -nsub*nsim;
  setAttrib(df, R_RowNamesSymbol, sexp_rownames);

  SEXP sexp_rownamesCount = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(sexp_rownamesCount)[0] = NA_INTEGER;
  INTEGER(sexp_rownamesCount)[1] = -nsub*nsim;
  setAttrib(countsDf , R_RowNamesSymbol, sexp_rownamesCount);

  SEXP sexp_colnames = PROTECT(allocVector(STRSXP,n+sm+md)); pro++;
  SEXP sexp_colnamesCount = PROTECT(allocVector(STRSXP,3+sm+md)); pro++;
  jj = 0;
  if (sm){
    SET_STRING_ELT(sexp_colnames, jj, mkChar("sim.id"));
    SET_STRING_ELT(sexp_colnamesCount, jj, mkChar("sim.id"));
    jj++;
  }
  // id
  kk = 0;
  if (md){
    SET_STRING_ELT(sexp_colnames, jj, mkChar("id"));
    SET_STRING_ELT(sexp_colnamesCount, jj, mkChar("id"));
    jj++; kk++;
  }

  SET_STRING_ELT(sexp_colnamesCount, sm+md, mkChar("slvr"));
  SET_STRING_ELT(sexp_colnamesCount, sm+md+1, mkChar("dadt"));
  SET_STRING_ELT(sexp_colnamesCount, sm+md+2, mkChar("jac"));
  
  kk=0;
  for (i = 0; i < ncov; i++){
    SET_STRING_ELT(covsn,kk, STRING_ELT(paramNames, par_cov[i]-1));
    kk++;
  }
  
  for (i = 0; i < npar; i++){
    is_cov=0;
    for (j = 0; j < ncov; j++){
      if (par_cov[j]-1 == i){
        is_cov=1;
        break;
      }
    }
    if (!is_cov){
      SET_STRING_ELT(sexp_colnames, jj, STRING_ELT(paramNames,i));
      jj++;
    }
  }
  setAttrib(df, R_NamesSymbol, sexp_colnames);
  setAttrib(countsDf , R_NamesSymbol, sexp_colnamesCount);
  setAttrib(covs, R_NamesSymbol, covsn);
  
  SEXP cls = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(cls, 0, mkChar("data.frame"));
  
  SEXP cls1 = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(cls1, 0, mkChar("data.frame"));
  classgets(df, cls);
  
  SEXP cls2 = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(cls2, 0, mkChar("data.frame"));
  classgets(covs, cls1);

  SEXP cls3 = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(cls3, 0, mkChar("data.frame"));
  classgets(countsDf, cls3);
  
  SEXP ret = PROTECT(allocVector(VECSXP,7)); pro++;
  SET_VECTOR_ELT(ret, 0, df);
  SET_VECTOR_ELT(ret, 1, dfe);
  SET_VECTOR_ELT(ret, 2, dfs);
  SET_VECTOR_ELT(ret, 3, dfd);
  SET_VECTOR_ELT(ret, 4, isEt);
  if (ncov == 0){
    SEXP covsn2 = PROTECT(R_NilValue);pro++;
    SET_VECTOR_ELT(ret, 5, covsn2);
  } else {
    SET_VECTOR_ELT(ret, 5, covs);
  }
  SET_VECTOR_ELT(ret, 6, countsDf);
  UNPROTECT(pro);
  return ret;
}

extern double *rxGetErrs();
extern int rxGetErrsNcol();

extern SEXP RxODE_df(int doDose){
  rx_solve *rx;
  rx = &rx_global;
  rx_solving_options *op = &op_global;
  int add_cov = rx->add_cov;
  int ncov = op->ncov;
  int nlhs = op->nlhs;
  int nobs = rx->nobs;
  int nsim = rx->nsim;
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
  int nall = 0;
  for (csub = 0; csub < nsub; csub++){
    ind = &(rx->subjects[csub]);
    nall+=ind->n_all_times;
  }
  SEXP df = PROTECT(allocVector(VECSXP,ncols+nidCols+doseCols)); pro++;
  for (i = 0; i < nidCols; i++){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(INTSXP, (doDose == 1 ? nall : nobs)*nsim))); pro++;
  }
  double *par_ptr;
  int nrow = (doDose == 1 ? nall : nobs)*nsim;
  
  double *errs = rxGetErrs();
  int updateErr = 0;
  int errNcol = rxGetErrsNcol();
  if (errNcol > 0){
    updateErr = 1;
  }
  if (doDose){
    SET_VECTOR_ELT(df, i++, PROTECT(allocVector(INTSXP, nall*nsim))); pro++;
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(REALSXP, nall*nsim))); pro++;
  }
  for (i = md + sm + doseCols; i < ncols + doseCols+nidCols; i++){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(REALSXP, (doDose == 1 ? nall : nobs)*nsim))); pro++;
  }
  
  // Now create the data frame
  double *dfp;
  int *dfi;
  int ii=0, jj = 0, ntimes;
  double *solve;
  double **cov_ptrA;
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
      cov_ptrA = ind->cov_ptr;
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
        if (updateErr){
          for (j=0; j < errNcol; j++){
	    par_ptr[svar[j]] = errs[nrow*j+kk];
          }
	  kk++;
        }
        jj  = 0 ;
	evid = rxEvidP(i,rx,neq[1]);
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
          // LHS
          if (nlhs){
	    rxCalcLhsP(i, rx, neq[1]);
	    for (j = 0; j < nlhs; j++){
	      dfp = REAL(VECTOR_ELT(df, jj));
               dfp[ii] =rxLhsP(j, rx, neq[1]);
	       jj++;
             }
          }
          // Cov
          if (add_cov*ncov > 0){
	    for (j = 0; j < add_cov*ncov; j++){
              dfp = REAL(VECTOR_ELT(df, jj));
	      // is this ntimes = nAllTimes or nObs time for this subject...?
	      cov_ptr = cov_ptrA[j];
	      dfp[ii] = (evid == 0 ? cov_ptr[i] : NA_REAL);
	      jj++;
	    }
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
  INTEGER(sexp_rownames)[1] = -(doDose == 1 ? nall : nobs)*nsim;
  setAttrib(df, R_RowNamesSymbol, sexp_rownames);
  SEXP sexp_colnames = PROTECT(allocVector(STRSXP,ncols+nidCols+doseCols)); pro++;
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
  // Put in LHS names
  SEXP lhsNames = PROTECT(rxLhsNames(op->modNamePtr)); pro++;
  for (i = 0; i < nlhs; i++){
    SET_STRING_ELT(sexp_colnames, jj, STRING_ELT(lhsNames,i));
    jj++;
  }
  // Put in Cov names
  SEXP paramNames = PROTECT(rxParamNames(op->modNamePtr)); pro++;
  int *par_cov = op->par_cov;
  for (i = 0; i < ncov*add_cov; i++){
    SET_STRING_ELT(sexp_colnames,jj, STRING_ELT(paramNames, par_cov[i]-1));
    jj++;
  }
  setAttrib(df, R_NamesSymbol, sexp_colnames);
  UNPROTECT(pro);
  return df;
}

int *gidoseSetup(int n);
int *gsiVSetup(int n);
// rxSolveOldC
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
  rx_solve *rx = &rx_global;
  rx_solving_options *op = &op_global;
  rx_solving_options_ind *ind = &inds_global[0];
  int i;
  // Counters
  ind->slvr_counter   = 0;
  ind->dadt_counter   = 0;
  ind->jac_counter    = 0;

  ind->InfusionRate = global_InfusionRate(*neqa);
  memset(ind->InfusionRate, 0.0, *neqa);

  ind->BadDose = global_BadDose(*neqa);
  memset(ind->BadDose, 0, *neqa);
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
  memset(op->scale, 1.0, *neqa);
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
  memset(rx->stateIgnore, 0, *neqa);
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
      calc_lhs(0, timep[i], retp+i*(*neqa), lhsp+i*(*nlhsa));
    }
  }
}


#define aexists(a, env) if (Rf_findVarInFrame(env, Rf_install(a)) == R_UnboundValue){ error("need '%s' in environment for solving.",a);}
void RxODE_ode_solve_env(SEXP sexp_rho){
  if(!isEnvironment(sexp_rho)){
    error("Calling RxODE_ode_solve_env without an environment...");
  }
  int pro = 0, i = 0;
  aexists("params",sexp_rho);
  SEXP sexp_theta = PROTECT(findVar(installChar(mkChar("params")),sexp_rho));pro++;
  aexists("inits",sexp_rho);
  SEXP sexp_inits = PROTECT(findVar(installChar(mkChar("inits")),sexp_rho)); pro++;
  aexists("lhs_vars",sexp_rho);
  SEXP sexp_lhs   = PROTECT(findVar(installChar(mkChar("lhs_vars")),sexp_rho)); pro++;
  // Events
  aexists("time",sexp_rho);
  SEXP sexp_time = PROTECT(findVar(installChar(mkChar("time")),sexp_rho)); pro++;
  aexists("evid",sexp_rho);
  SEXP sexp_evid = PROTECT(findVar(installChar(mkChar("evid")),sexp_rho)); pro++;
  aexists("amt",sexp_rho);
  SEXP sexp_dose = PROTECT(findVar(installChar(mkChar("amt")),sexp_rho)); pro++;
  // Covariates
  aexists("pcov",sexp_rho);
  SEXP sexp_pcov = PROTECT(findVar(installChar(mkChar("pcov")),sexp_rho)); pro++;
  aexists("covs",sexp_rho);
  SEXP sexp_cov = PROTECT(findVar(installChar(mkChar("covs")),sexp_rho)); pro++;
  aexists("isLocf",sexp_rho);
  SEXP sexp_locf = PROTECT(findVar(installChar(mkChar("isLocf")),sexp_rho)); pro++;
  // Solver Options
  aexists("atol",sexp_rho);
  SEXP sexp_atol = PROTECT(findVar(installChar(mkChar("atol")),sexp_rho)); pro++;
  aexists("rtol",sexp_rho);
  SEXP sexp_rtol = PROTECT(findVar(installChar(mkChar("rtol")),sexp_rho)); pro++;
  aexists("hmin",sexp_rho);
  SEXP sexp_hmin = PROTECT(findVar(installChar(mkChar("hmin")),sexp_rho)); pro++;
  aexists("hmax",sexp_rho);
  SEXP sexp_hmax = PROTECT(findVar(installChar(mkChar("hmax")),sexp_rho)); pro++;
  aexists("hini",sexp_rho);
  SEXP sexp_h0 = PROTECT(findVar(installChar(mkChar("hini")),sexp_rho)); pro++;
  aexists("maxordn",sexp_rho);
  SEXP sexp_mxordn = PROTECT(findVar(installChar(mkChar("maxordn")),sexp_rho)); pro++;
  aexists("maxords",sexp_rho);
  SEXP sexp_mxords = PROTECT(findVar(installChar(mkChar("maxords")),sexp_rho)); pro++;
  aexists("maxsteps",sexp_rho);
  SEXP sexp_mx = PROTECT(findVar(installChar(mkChar("maxsteps")),sexp_rho)); pro++;
  aexists("stiff",sexp_rho);
  SEXP sexp_stiff = PROTECT(findVar(installChar(mkChar("stiff")),sexp_rho)); pro++;
  aexists("transit_abs",sexp_rho);
  SEXP sexp_transit_abs = PROTECT(findVar(installChar(mkChar("transit_abs")),sexp_rho)); pro++;
  aexists("rc",sexp_rho);
  SEXP sexp_rc = PROTECT(findVar(installChar(mkChar("rc")),sexp_rho)); pro++;
  rx_solve *rx = &rx_global;
  rx_solving_options *op = &op_global;
  rx_solving_options_ind *ind = &inds_global[0];
  ind->rc=INTEGER(sexp_rc);

  // Let R handle deallocating the solve and lhs expressions; Should disappear with evironment
  SEXP sexp_solve = PROTECT(allocVector(REALSXP,length(sexp_time)*length(sexp_inits))); pro++;
  ind->solve = REAL(sexp_solve);
  memset(ind->solve,0,length(sexp_solve));
  SEXP sexp_lhsV = PROTECT(allocVector(REALSXP,length(sexp_time)*length(sexp_lhs))); pro++;
  ind->lhs = REAL(sexp_lhsV);
  memset(ind->lhs,0,length(sexp_time)*length(sexp_lhs));
  op->stiff = INTEGER(sexp_stiff)[0];
  op->do_transit_abs = INTEGER(sexp_transit_abs)[0];
  op->ATOL = REAL(sexp_atol)[0];
  op->RTOL= REAL(sexp_rtol)[0];
  op->mxstep=  INTEGER(sexp_mx)[0];
  op->HMIN = REAL(sexp_hmin)[0];
  op->H0 = REAL(sexp_h0)[0];
  op->MXORDN= INTEGER(sexp_mxordn)[0];
  op->MXORDS = INTEGER(sexp_mxords)[0];
  op->par_cov = INTEGER(sexp_pcov);
  op->is_locf = INTEGER(sexp_locf)[0];
  ind->HMAX = REAL(sexp_hmax)[0];
  ind->par_ptr = REAL(sexp_theta);
  op->inits   = REAL(sexp_inits);
  ind->dose    = REAL(sexp_dose);
  /* ind->solve   = retp; */
  ind->evid    = INTEGER(sexp_evid);
  ind->cov_ptr = REAL(sexp_cov);
  
  ind->ndoses        = -1;
  ind->all_times     = REAL(sexp_time);
  ind->n_all_times   = length(sexp_time);
  // Covariates
  op->do_par_cov    = 1;
  op->ncov  = length(sexp_pcov);
  // Solver options
  op->do_transit_abs = INTEGER(sexp_transit_abs)[0];
  op->stiff          = INTEGER(sexp_stiff)[0];
  ind->slvr_counter   = 0;
  ind->dadt_counter   = 0;
  ind->jac_counter    = 0;
  // LOCF
  if (op->is_locf == 1){
    op->f2 = 0.0; //= f=0 
    op->f1 = 1.0; // = 1-f = 1;
    op->kind = 0;
  } else if (op->is_locf == 2) {
    // NOCB
    op->f2 = 1.0; //= f=1
    op->f1 = 0.0;
    op->kind = 0;
  } else if (op->is_locf == 3){
    op->f2 = 0.5; //= f=0.5
    op->f1 = 0.5;
    op->kind = 0;
  } else {
    // Linear
    op->f2 = 1.0; //= f=0
    op->f1 = 0.0;
    op->kind = 1;
  }
  op->nlhs          = length(sexp_lhs);
  op->neq           = length(sexp_inits);
  ind->nBadDose = 0;
  ind->InfusionRate = global_InfusionRate(op->neq);
  memset(ind->InfusionRate, 0.0, op->neq);
  ind->BadDose = global_BadDose(op->neq);
  memset(ind->BadDose, 0, op->neq);
  /* rx_solve *rx = rxSingle(mv, stiff, transit_abs, atol, rtol, mx, hmin, h0,  mxordn, */
  /*                         mxords, 1, length(sexp_pcov), pcov, 1,  locf, */
  /*                         hmax, theta, dose, solve, lhs, evid, rce, cov, */
  /*                         length(sexp_time), time); */
  /* rxode_assign_rx(rx); */
  _globalRx=rx;
  rx->op = &op_global;
  par_solve(rx); // Solve without the option of updating residuals.
  if (op->nlhs) {
    for (i=0; i<ind->n_all_times; i++){
      // 0 = first subject; Calc lhs changed...
      calc_lhs(0, ind->all_times[i], ind->solve+i*(op->neq), ind->lhs+i*(op->nlhs));
    }
  }
  defineVar(install(".lhs"), sexp_lhsV, sexp_rho);
  defineVar(install(".solve"), sexp_solve, sexp_rho);
  UNPROTECT(pro);
}
