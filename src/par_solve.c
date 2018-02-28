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

inline rx_solving_options_ind *rxOptionsIniEnsure(int mx){
  if (mx >= max_inds_global){
    max_inds_global = mx+1024;
    inds_global = Realloc(inds_global, max_inds_global, rx_solving_options_ind);
  }
  return inds_global;
}

t_dydt dydt = NULL;

t_calc_jac calc_jac = NULL;

t_calc_lhs g_calc_lhs = NULL;
extern void calc_lhs(int cSub, double t, double *A, double *lhs){
  if (g_calc_lhs == NULL) error("RxODE library not setup correctly.");
  g_calc_lhs(cSub, t, A, lhs);
}

t_update_inis update_inis = NULL;

t_dydt_lsoda_dum dydt_lsoda_dum = NULL;

t_dydt_liblsoda dydt_liblsoda = NULL;

t_jdum_lsoda jdum_lsoda = NULL;

t_set_solve g_set_solve = NULL;
extern void set_solve(rx_solve *rx){
  if (g_set_solve == NULL) error("RxODE library not setup correctly.");
  g_set_solve(rx);
}

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
  g_calc_lhs =(t_calc_lhs) R_GetCCallable(lib, s_calc_lhs);
  dydt =(t_dydt) R_GetCCallable(lib, s_dydt);
  calc_jac =(t_calc_jac) R_GetCCallable(lib, s_calc_jac);
  update_inis =(t_update_inis) R_GetCCallable(lib, s_inis);
  dydt_lsoda_dum =(t_dydt_lsoda_dum) R_GetCCallable(lib, s_dydt_lsoda_dum);
  jdum_lsoda =(t_jdum_lsoda) R_GetCCallable(lib, s_dydt_jdum_lsoda);
  g_set_solve = (t_set_solve)R_GetCCallable(lib, s_ode_solver_solvedata);
  get_solve = (t_get_solve)R_GetCCallable(lib, s_ode_solver_get_solvedata);
  dydt_liblsoda = (t_dydt_liblsoda)R_GetCCallable(lib, s_dydt_liblsoda);
  global_jt = 2;
  global_mf = 22;  
  global_debug = 0;
}

void rxClearFuns(){
  g_calc_lhs		= NULL;
  dydt			= NULL;
  calc_jac		= NULL;
  update_inis		= NULL;
  dydt_lsoda_dum	= NULL;
  jdum_lsoda		= NULL;
  g_set_solve		= NULL;
  get_solve		= NULL;
  dydt_liblsoda		= NULL;
}

void getSolvingOptionsIndPtr(double *InfusionRate,
                             int *BadDose,
                             double HMAX, // Determined by diff
                             double *par_ptr,
                             double *dose,
			     int *idose,
                             double *solve,
                             double *lhs,
                             int *evid,
                             int *rc,
                             double *cov_ptr,
                             int n_all_times,
                             double *all_times,
                             int id,
                             int sim,
                             rx_solving_options_ind *o){
  o->slvr_counter = 0;
  o->dadt_counter = 0;
  o->jac_counter  = 0;
  o->InfusionRate = InfusionRate;
  o->BadDose = BadDose;
  o->nBadDose = 0;
  o->HMAX = HMAX; // Determined by diff
  o->tlast = 0.0;
  o->podo = 0.0;
  o->par_ptr = par_ptr;
  o->dose = dose;
  o->solve = solve;
  o->lhs = lhs;
  o->evid = evid;
  o->rc = rc;
  o->cov_ptr = cov_ptr;
  o->n_all_times = n_all_times;
  o->ixds = 0;
  o->ndoses = -1;
  o->all_times = all_times;
  o->idose = idose;
  o->id = id;
  o->sim = sim;
}

void getSolvingOptionsPtr(double ATOL,          //absolute error
                          double RTOL,          //relative error
                          double H0,
                          double HMIN,
                          int mxstep,
                          int MXORDN,
                          int MXORDS,
                          // Approx options
                          int do_transit_abs,
                          int nlhs,
                          int neq,
                          int stiff,
                          double f1,
                          double f2,
                          int kind,
                          int is_locf,
			  int cores,
			  int ncov,
			  int *par_cov,
                          int do_par_cov,
                          double *inits,
			  double *scale,
                          const char *modNamePtr,
			  double hmax2,
			  double *atol2,
			  double *rtol2,
                          int nDisplayProgress,
                          int ncoresRV,
                          int isChol,
			  int *svar){
  // This really should not be called very often, so just allocate one for now.
  rx_solving_options *o;
  o = &op_global;
  o->badSolve = 0;
  o->ATOL = ATOL;          //absolute error
  o->RTOL = RTOL;          //relative error
  o->H0 = H0;
  o->HMIN = HMIN;
  o->mxstep = mxstep;
  o->MXORDN = MXORDN;
  o->MXORDS = MXORDS;
  o->do_transit_abs = do_transit_abs;
  o->nlhs = nlhs;
  o->neq = neq;
  o->stiff = stiff;
  o->f1 = f1;
  o->f2 = f2;
  o->kind = kind;
  o->is_locf = is_locf;
  o->ncov=ncov;
  o->par_cov = par_cov;
  o->do_par_cov = do_par_cov;
  o->inits = inits;
  o->scale = scale;
  o->extraCmt = 0;
  o->hmax2 = hmax2; // Determined by diff
  o->rtol2 = rtol2;
  o->atol2 = atol2;
  o->cores = cores;
  o->nDisplayProgress = nDisplayProgress;
  o->ncoresRV = ncoresRV;
  o->isChol = isChol;
  o->svar = svar;
  o->abort = 0;
  sprintf(o->modNamePtr, "%s", modNamePtr);
}

void rxSolveData(rx_solving_options_ind *subjects,
                 int nsub,
                 int nsim,
		 int *stateIgnore,
		 int nobs,
                 int add_cov,
                 int matrix){
  rx_solve *o;
  o = &rx_global;//(rx_solve*)R_chk_calloc(1,sizeof(*o));
  o->subjects = subjects;
  o->nsub = nsub;
  o->nsim = nsim;
  o->stateIgnore = stateIgnore;
  o->nobs = nobs;
  o->add_cov = add_cov;
  o->matrix = matrix;
}

void F77_NAME(dlsoda)(
                      void (*)(int *, double *, double *, double *),
                      int *, double *, double *, double *, int *, double *, double *,
                      int *, int *, int *, double *,int *,int *, int *,
                      void (*)(int *, double *, double *, int *, int *, double *, int *),
                      int *);

rx_solve *getRxSolve(SEXP ptr);
extern rx_solve *getRxSolve_(){
  return &rx_global;
}

rx_solving_options *getRxOp(rx_solve *rx){
  /* if(!R_ExternalPtrAddr(rx->op)){ */
  /*   error("Cannot get global ode solver options."); */
  /* } */
  /* return (rx_solving_options*)(R_ExternalPtrAddr(rx->op)); */
  return &op_global;
}

rx_solving_options_ind *getRxId(rx_solve *rx, unsigned int id){
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
inline double *global_InfusionRate(unsigned int mx){
  if (mx >= global_InfusionRatei){
    global_InfusionRatei = mx+1024;
    global_InfusionRatep = Realloc(global_InfusionRatep, global_InfusionRatei, double);
  }
  return global_InfusionRatep;
}


int *global_BadDosep;
unsigned int global_BadDosei = 0;
inline int *global_BadDose(unsigned int mx){
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


/* Authors: Robert Gentleman and Ross Ihaka and The R Core Team */
/* Taken directly from https://github.com/wch/r-source/blob/922777f2a0363fd6fe07e926971547dd8315fc24/src/library/stats/src/approx.c*/
/* Changed as follows:
   - Different Name
   - Use RxODE structure
*/

static double rx_approxP(double v, double *x, double *y, int n,
                         rx_solving_options *Meth, rx_solving_options_ind *id)
{
  /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
  int i, j, ij;

  if(!n) return R_NaN;

  i = 0;
  j = n - 1;

  /* handle out-of-domain points */
  if(v < x[i]) return id->ylow;
  if(v > x[j]) return id->yhigh;

  /* find the correct interval by bisection */
  while(i < j - 1) { /* x[i] <= v <= x[j] */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    if(v < x[ij]) j = ij; else i = ij;
    /* still i < j */
  }
  /* provably have i == j-1 */

  /* interpolation */

  if(v == x[j]) return y[j];
  if(v == x[i]) return y[i];
  /* impossible: if(x[j] == x[i]) return y[i]; */

  if(Meth->kind == 1) /* linear */
    return y[i] + (y[j] - y[i]) * ((v - x[i])/(x[j] - x[i]));
  else /* 2 : constant */
    return (Meth->f1 != 0.0 ? y[i] * Meth->f1 : 0.0)
      + (Meth->f2 != 0.0 ? y[j] * Meth->f2 : 0.0);
}/* approx1() */

/* End approx from R */

rx_solve *_globalRx = NULL;

extern void rxode_assign_rx(rx_solve *rx){
  _globalRx=rx;
}

extern void update_par_ptrP(double t, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  rx_solving_options *op = &op_global;
  if (op->neq > 0){
    // Update all covariate parameters
    int k;
    int *par_cov = op->par_cov;
    double *par_ptr = ind->par_ptr;
    double *all_times = ind->all_times;
    double *cov_ptr = ind->cov_ptr;
    int ncov = op->ncov;
    if (op->do_par_cov){
      for (k = 0; k < ncov; k++){
        if (par_cov[k]){
          // Use the same methodology as approxfun.
          // There is some rumor the C function may go away...
          ind->ylow = cov_ptr[ind->n_all_times*k];
          ind->yhigh = cov_ptr[ind->n_all_times*k+ind->n_all_times-1];
          par_ptr[par_cov[k]-1] = rx_approxP(t, all_times, cov_ptr+ind->n_all_times*k, ind->n_all_times, op, ind);
        }
        if (global_debug){
          REprintf("par_ptr[%d] (cov %d/%d) = %f\n",par_cov[k]-1, k,ncov,cov_ptr[par_cov[k]-1]);
        }
      }
    }
  }
}
extern void update_par_ptr(double t){
  update_par_ptrP(t, _globalRx, 0);
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

extern double RxODE_par_ptrP(int val, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  double ret =ind->par_ptr[val];
  return ret;
}

extern double RxODE_par_ptr(int val){
  return RxODE_par_ptrP(val, _globalRx, 0);
}

extern long RxODE_jac_counter_valP(rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  return ind->jac_counter;
}
extern long RxODE_jac_counter_val(){
  return RxODE_jac_counter_valP(_globalRx, 0); // Not sure this function is needed...
}

extern long RxODE_dadt_counter_valP(rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  return ind->dadt_counter;
}

extern long RxODE_dadt_counter_val(){
  return RxODE_dadt_counter_val(_globalRx, 0);
}

extern void RxODE_jac_counter_incP(rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  ind->jac_counter++;
}

extern void RxODE_jac_counter_inc(){
  RxODE_jac_counter_incP(_globalRx, 0);
}

extern void RxODE_dadt_counter_incP(rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  ind->dadt_counter++;
}

extern void RxODE_dadt_counter_inc(){
  RxODE_dadt_counter_incP(_globalRx, 0);
}

extern double RxODE_podoP(rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  return ind->podo;
}

extern double RxODE_podo(){
  return RxODE_podo(_globalRx, 0);
}

extern double RxODE_tlastP(rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  return ind->tlast;
}

extern double RxODE_tlast(){
  return RxODE_tlastP(_globalRx, 0);
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

extern double RxODE_transit4P(double t, rx_solve *rx, unsigned int id, double n, double mtt, double bio){
  double ktr = (n+1)/mtt;
  double lktr = log(n+1)-log(mtt);
  return exp(log(bio*RxODE_podoP(rx, id))+lktr+n*(lktr+log(t))-ktr*t-lgamma1p(n));
}

extern double RxODE_transit4(double t, double n, double mtt, double bio){
  return RxODE_transit4P(t, _globalRx, 0, n, mtt, bio);
}

extern double RxODE_transit3P(double t, rx_solve *rx, unsigned int id, double n, double mtt){
  return RxODE_transit4P(t, rx, id, n,mtt, 1.0);
}

extern double RxODE_transit3(double t, double n, double mtt, rx_solve *rx, unsigned int id){
  return RxODE_transit4P(t, _globalRx, 0, n,mtt, 1.0);
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
  rx = getRxSolve(sd);
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
  
  double *times, *doses, *cov_ptr;
  int evid, iie = 0, iis = 0, iid = 0, curdose, ntimes;
  int k, kk;
  for (csub = 0; csub < nsub; csub++){
    ind = &(rx->subjects[csub]);
    ntimes = ind->n_all_times;
    times = ind->all_times;
    doses = ind->dose;
    cov_ptr = ind->cov_ptr;
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
            dfp[iis] = cov_ptr[kk*ntimes+iis];
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

extern SEXP rxSimSigmaC(rx_solving_options *op, int nObs);
extern double *rxGetErrs();
extern int rxGetErrsNcol();

extern SEXP RxODE_df(SEXP sd, int doDose){
  rx_solve *rx;
  rx = getRxSolve(sd);
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
	      dfp[ii] = (evid == 0 ? cov_ptr[j*ntimes+i] : NA_REAL);
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
  ind->InfusionRate = global_InfusionRate(*neqa);
  memset(ind->InfusionRate, 0.0, *neqa);
  ind->BadDose = global_BadDose(*neqa);
  memset(ind->BadDose, 0, *neqa);
  ind->ndoses = -1;
  ind->all_times         = timep;
  ind->n_all_times       = *ntime;
  //RxODE_ode_dosing_ = RxODE_ode_get_dosing();
  // par_cov
  op->do_par_cov        = 0;
  // cov_ptr
  op->ncov              = 0;
  op->is_locf           = 0;
  // Solver Options
  op->ATOL = *atol;
  op->RTOL = *rtol;
  // Assign to default LSODA behvior, or 0
  op->HMIN           = 0;
  ind->HMAX           = 0;
  op->H0             = 0;
  op->MXORDN         = 0;
  op->MXORDS         = 0;
  op->mxstep         = 5000; // Not LSODA default but RxODE default
  // Counters
  ind->slvr_counter   = 0;
  ind->dadt_counter   = 0;
  ind->jac_counter    = 0;

  op->nlhs           = *nlhsa;
  op->neq            = *neqa;
  op->stiff          = *stiffa;
  
  

  ind->nBadDose = 0;
  /* int *BadDose; */
  op->do_transit_abs = *transit_abs;
  
  ind->par_ptr = theta;
  op->inits   = initsp;
  ind->dose    = dosep;
  ind->solve   = retp;
  ind->lhs     = lhsp;
  ind->evid    = evidp;

  /* int i =0; */
  /* rx = rxSingle(mv, *stiffa,*transit_abs, *atol, *rtol, 5000,//maxsteps */
  /*               0, 0, 12, 5, 1, 0, par_cov, 0, 0, 0, theta, */
  /*               dosep, retp, lhsp, evidp, rc, cov, *ntime,timep); */
  _globalRx=rx;
  par_solve(rx); // Solve without the option of updating residuals.
  if (*nlhsa) {
    for (i=0; i<*ntime; i++){
      // 0 = first subject; Calc lhs changed...
      calc_lhs(0, timep[i], retp+i*(*neqa), lhsp+i*(*nlhsa));
    }
  }
}
