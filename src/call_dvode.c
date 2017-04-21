#include <stdio.h>
#include <stdlib.h>
#include "dop853.h"
#define max(a, b) ((a) > (b) ? (a) : (b))
#define RLEN (neq+nlhs+ncov+1)
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>


void F77_NAME(dlsoda)(
		      void (*)(int *, double *, double *, double *),
		      int *, double *, double *, double *, int *, double *, double *,
		      int *, int *, int *, double *,int *,int *, int *,
		      void (*)(int *, double *, double *, int *, int *, double *, int *),
		      int *);

long slvr_counter, dadt_counter, jac_counter;
double InfusionRate[100];
double ATOL;		//absolute error
double RTOL;		//relative error
double HMAX;
double H0;
double HMIN;

int    do_transit_abs=0;
double tlast=0;
double podo=0;
double *par_ptr, *inits, *dose, *solve, *lhs;
int    *par_cov, *evid;
int rc[1];
double *cov_ptr;
int    ncov, nlhs, neq, stiff;
int    is_locf;
int    n_all_times;
int    mxstep;
int    MXORDN;
int    MXORDS;
int    global_jt, global_mf, global_debug, ixds;
double *all_times;
FILE *fp;

/* void __DYDT__(unsigned int neq, double t, double *A, double *DADT); */
/* void __CALC_LHS__(double t, double *A, double *lhs); */
/* void __CALC_JAC__(unsigned int neq, double t, double *A, double *JAC, unsigned int __NROWPD__); */

void (*dydt)(unsigned int neq, double t, double *A, double *DADT);
void (*calc_jac)(unsigned int neq, double t, double *A, double *JAC, unsigned int __NROWPD__);
void (*calc_lhs)(double t, double *A, double *lhs);
void (*update_inis)(SEXP _ini_sexp);

void rxInner(SEXP rho);


void rxCalcLhs(int i){
  if (i < n_all_times){
    calc_lhs(all_times[i], solve+i*neq, lhs);
  } else {
    error("LHS cannot be calculated (%dth entry).",i);
  }
}

double rxLhs(int i){
  if (i < nlhs){
    return(lhs[i]);
  } else {
    error("Trying to access an equation that isn't calculated. lhs(%d)\n",i);
  }
}

int rxEvid(int i){
  if (i < n_all_times){
    return(evid[i]);
  } else {
    error("Trying to access EVID outside of defined events.\n");
  }
}

unsigned int nAllTimes (){
  return (unsigned int)(n_all_times);
}

unsigned int nObs(){
  return (unsigned int)(n_all_times - ixds);
}

unsigned int nLhs (){
  return (unsigned int)(nlhs);
}

int nEq (){
  return neq;
}

double RxODE_as_zero(double x){
  if (abs(x) < sqrt(DOUBLE_EPS)){
    return(0.0);
  } else {
    return(x);
  }
}

double RxODE_safe_log(double x){
  if (x == 0.0){
    return log(sqrt(DOUBLE_EPS));
  } else {
    return log(x);
  }
}

double RxODE_safe_zero(double x){
  if (x == 0){
    return sqrt(DOUBLE_EPS);
  } else {
    return(x);
  }
}


/* Taken directly from https://github.com/wch/r-source/blob/922777f2a0363fd6fe07e926971547dd8315fc24/src/library/stats/src/approx.c*/

typedef struct {
  double ylow;
  double yhigh;
  double f1;
  double f2;
  int kind;
} rx_appr_meth;

rx_appr_meth rx_aprox_M = {0.0, 0.0, 0.0, 0.0, 0};

static double rx_approx1(double v, double *x, double *y, int n,
			 rx_appr_meth *Meth)
{
  /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
  int i, j, ij;

  if(!n) return R_NaN;

  i = 0;
  j = n - 1;

  /* handle out-of-domain points */
  if(v < x[i]) return Meth->ylow;
  if(v > x[j]) return Meth->yhigh;

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

void update_par_ptr(double t){
  // Update all covariate parameters
  int k;
  for (k = 0; k < ncov; k++){
    if (par_cov[k]){
      // Use the same methodology as approxfun.
      // There is some rumor the C function may go away...
      rx_aprox_M.ylow = cov_ptr[nAllTimes()*k];
      rx_aprox_M.yhigh = cov_ptr[nAllTimes()*k+nAllTimes()-1];
      par_ptr[par_cov[k]-1] = rx_approx1(t, all_times, cov_ptr+nAllTimes()*k, nAllTimes(), &rx_aprox_M);
    }
    if (global_debug){
      Rprintf("par_ptr[%d] (cov %d/%d) = %f\n",par_cov[k]-1, k,ncov,cov_ptr[par_cov[k]-1]);
    }
  }
}

//--------------------------------------------------------------------------
void dydt_lsoda_dum(int *neq, double *t, double *A, double *DADT)
{
  dydt(*neq, *t, A, DADT);
}
void jdum_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){
  // Update all covariate parameters
  calc_jac(*neq, *t, A, JAC, *nrowpd);
}
void call_lsoda(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc)
{
  int i, j;
  double xout, xp=x[0], yp[99];
  int itol = 1;
  double  rtol = RTOL, atol = ATOL;
  // Set jt to 1 if full is specified.
  int itask = 1, istate = 1, iopt = 0, lrw=22+neq*max(16, neq+9), liw=20+neq, jt = global_jt;
  double *rwork;
  int *iwork;
  int wh, cmt;

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
    Rprintf("JT: %d\n",jt);
  rwork = (double*)R_alloc(lrw+1, sizeof(double));
  iwork = (int*)R_alloc(liw+1, sizeof(int));

  iopt = 1;
  
  rwork[4] = H0; // H0 -- determined by solver
  rwork[5] = HMAX; // Hmax -- Infinite
  rwork[6] = HMIN; // Hmin -- 0
  
  iwork[4] = 0; // ixpr  -- No extra printing.
  iwork[5] = mxstep; // mxstep 
  iwork[6] = 0; // MXHNIL 
  iwork[7] = MXORDN; // MXORDN 
  iwork[8] = MXORDS;  // MXORDS
  
  //--- inits the system
  for(i=0; i<neq; i++) yp[i] = inits[i];

  for(i=0; i<nx; i++)
    {
      wh = rxEvid(i);
      xout = x[i];
      if (global_debug){
	Rprintf("i=%d xp=%f xout=%f\n", i, xp, xout);
        fprintf(fp, "i=%d xp=%f xout=%f\n", i, xp, xout);
      }
      if(xout-xp> DBL_EPSILON*max(fabs(xout),fabs(xp)))
	{
	  F77_CALL(dlsoda)(dydt_lsoda_dum, &neq, yp, &xp, &xout, &itol, &rtol, &atol, &itask,
			   &istate, &iopt, rwork, &lrw, iwork, &liw, &jdum_lsoda, &jt);

	  if (istate<0)
	    {
	      Rprintf("IDID=%d, %s\n", istate, err_msg[-istate-1]);
	      *rc = istate;
	      error("Error solving using LSODA");
	      return;
	    }
	  slvr_counter++;
	  //dadt_counter = 0;
	}
      if (wh)
	{
	  cmt = (wh%10000)/100 - 1;
	  if (wh>10000)
	    {
	      InfusionRate[cmt] += dose[ixds];
	    }
	  else
	    {
	      if (do_transit_abs)
		{
		  podo = dose[ixds];
		  tlast = xout;
		}
	      else yp[cmt] += dose[ixds];	//dosing before obs
	    }
	  istate = 1;

	  ixds++;
	  xp = xout;
	}
      for(j=0; j<neq; j++) ret[neq*i+j] = yp[j];
      //Rprintf("wh=%d cmt=%d tm=%g rate=%g\n", wh, cmt, xp, InfusionRate[cmt]);

      if (global_debug){
	Rprintf("ISTATE=%d, ", istate);
	fprintf(fp, "ISTATE=%d, ", istate);
	for(j=0; j<neq; j++)
	  {
	    Rprintf("%f ", yp[j]);
	    fprintf(fp, "%f ", yp[j]);
	  }
	Rprintf("\n");
	fprintf(fp, "\n");
      }
    }
}

//dummy solout fn
void solout(long int nr, double t_old, double t,
	    double *y, unsigned int n, int *irtrn){}
void call_dop(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc)
{
  int i, j;
  //DE solver config vars
  double xout, xp=x[0], yp[99];
  double rtol=RTOL, atol=ATOL;
  int itol=0;		//0: rtol/atol scalars; 1: rtol/atol vectors
  int iout=0;		//iout=0: solout() NEVER called
  int idid=0;
  int wh, cmt;
  char *err_msg[]=
    {
      "input is not consistent",
      "larger nmax is needed",
      "step size becomes too small",
      "problem is probably stiff (interrupted)"
    };

  //--- inits the system
  for(i=0; i<neq; i++) yp[i] = inits[i];

  for(i=0; i<nx; i++)
    {
      wh = rxEvid(i);
      xout = x[i];
      if (global_debug){
	Rprintf("i=%d xp=%f xout=%f\n", i, xp, xout);
	fprintf(fp, "i=%d xp=%f xout=%f\n", i, xp, xout);
      }
		
      if(xout-xp>DBL_EPSILON*max(fabs(xout),fabs(xp)))
	{
	  idid = dop853(
			neq,      	/* dimension of the system <= UINT_MAX-1*/
			dydt,    	/* function computing the value of f(x,y) */
			xp,           /* initial x-value */
			yp,           /* initial values for y */
			xout,         /* final x-value (xend-x may be positive or negative) */
			&rtol,      	/* relative error tolerance */
			&atol,      	/* absolute error tolerance */
			itol,         /* switch for rtoler and atoler */
			solout,     	/* function providing the numerical solution during integration */
			iout,         /* switch for calling solout */
			NULL,       	/* messages stream */
			DBL_EPSILON, 	/* rounding unit */
			0,         	/* safety factor */
			0,         	/* parameters for step size selection */
			0,
			0,         	/* for stabilized step size control */
			0,         	/* maximal step size */
			0,            /* initial step size */
			0,            /* maximal number of allowed steps */
			1,            /* switch for the choice of the coefficients */
			-1,     		/* test for stiffness */
			0, 			/* number of components for which dense outpout is required */
			NULL, 		/* indexes of components for which dense output is required, >= nrdens */
			0  			/* declared length of icon */
			);
	  if (idid<0)
	    {
	      Rprintf("IDID=%d, %s\n", idid, err_msg[-idid-1]);
	      *rc = idid;
	      error("Error sovling using dop853");
	      return;  //exit(1);  // dj: should not abort R
	    }

	  xp = xRead();
	  slvr_counter++;
	  //dadt_counter = 0;
	}
      if (wh)
	{
	  cmt = (wh%10000)/100 - 1;
	  if (wh>10000)
	    {
	      InfusionRate[cmt] += dose[ixds];
	    }
	  else
	    {
	      if (do_transit_abs)
		{
		  podo = dose[ixds];
		  tlast = xout;
		}
	      else yp[cmt] += dose[ixds];	//dosing before obs
	    }
	  ixds++;
	  xp = xout;
	}
      for(j=0; j<neq; j++) ret[neq*i+j] = yp[j];
      //Rprintf("wh=%d cmt=%d tm=%g rate=%g\n", wh, cmt, xp, InfusionRate[cmt]);

      if (global_debug){
	Rprintf("IDID=%d, ", idid);
	fprintf(fp, "IDID=%d, ", idid);
	for(j=0; j<neq; j++)
	  {
	    Rprintf("%f ", yp[j]);
	    fprintf(fp, "%f ", yp[j]);
	  }
	Rprintf("\n");
	fprintf(fp, "\n");
      }
    }
}

void RxODE_ode_solver_c(int neq, int stiff, int *evid, double *inits, double *dose, double *solve, int *rc){
  ixds = 0;
  if (neq) {
    if (stiff==0){
      call_dop(neq, all_times, evid, nAllTimes(), inits, dose, solve, rc);
    } else{
      call_lsoda(neq, all_times, evid, nAllTimes(), inits, dose, solve, rc);
    }
  }
}

void RxODE_ode_solver_old_c(int *neq,
                            double *theta,  //order:
                            double *time,
                            int *evid,
                            int *ntime,
                            double *inits,
                            double *dose,
                            double *ret,
                            double *atol,
                            double *rtol,
                            int *stiff,
                            int *transit_abs,
                            int *nlhs,
                            double *lhs,
                            int *rc){
  int i;
  /* for (i=0; i< *neq; i++) InfusionRate[i] = 0.0; */
  ATOL = *atol;
  RTOL = *rtol;
  do_transit_abs = *transit_abs;
  par_ptr = theta;
  // Assign to default LSODA behvior, or 0
  HMIN           = 0;
  HMAX           = 0;
  H0             = 0;
  MXORDN         = 0;
  MXORDS         = 0;
  mxstep         = 0;
  slvr_counter   = 0;
  dadt_counter   = 0;
  jac_counter    = 0;
  // Assign global time information
  all_times     = time; 
  n_all_times   = *ntime;
  // Call solver
  RxODE_ode_solver_c(*neq, *stiff, evid, inits, dose, ret, rc);
  // Update LHS
  if (*nlhs) {
    for (i=0; i<*ntime; i++){
      calc_lhs(time[i], ret+i*(*neq), lhs+i*(*nlhs));
    }
  }
}

SEXP RxODE_ode_get_dosing(){
  SEXP sexp_ret;
  int ndoses = 0, i, j = 0;
  for (i = 0; i < n_all_times; i++){
    if (rxEvid(i)){
      ndoses++;
    }
  }
  sexp_ret  = PROTECT(allocMatrix(REALSXP, ndoses, 3));
  double *ret   = REAL(sexp_ret);
  for (i = 0; i < n_all_times; i++){
    if (rxEvid(i)){
      ret[j] = all_times[i];
      ret[ndoses+j] =(double)(evid[i]);
      ret[ndoses*2+j] = dose[j];
      /* Rprintf("\tj:%d(%d)\tt: %f\tevid:%f\tdose:%f\n", */
      /* 	      j, ndoses, */
      /* 	      ret[j],ret[ndoses+j],ret[ndoses*2+j]); */
      j++;
    }
  }
  // Unprotect after solving.
  UNPROTECT(1);
  return sexp_ret;
}

SEXP RxODE_ode_dosing_;

SEXP RxODE_ode_dosing(){
  return RxODE_ode_dosing_;
}

void RxODE_ode_free(){
  /* Free(InfusionRate); */
  /* Free(solve); */
  /* Free(lhs); */
  /* Free(rc); */
}

void RxODE_ode_alloc(){
  solve = (double*) R_alloc(neq*n_all_times+1, sizeof(double));
  lhs = (double *) R_alloc(nlhs,sizeof(double));
  /* InfusionRate = (double *)R_alloc(neq+2,sizeof(double)); */
  /* rc = (int *) R_alloc(1,sizeof(int)); */
  /* solve         = (double *) Calloc(neq*nAllTimes()+1,double); */
  /* lhs           = (double *) Calloc(nlhs,double); */
  /* InfusionRate  = (double *) Calloc(neq+2,double); */
  /* rc = (int *) Calloc(1,int); */
  rc[0] = 0;
}

void RxODE_ode_solver_0_6_c(int *neq,
                            double *theta,  //order:
                            double *time,
                            int *evid,
                            int *ntime,
                            double *inits,
                            double *dose,
                            double *ret,
                            double *atol,
                            double *rtol,
                            int *stiff,
                            int *transit_abs,
                            int *nlhs,
                            double *lhs,
                            int *rc,
			    double hmin,
			    double hmax,
			    double h0,
			    int mxordn,
			    int mxords,
			    int mxstep){
  int i;
  /* for (i=0; i<*neq; i++) InfusionRate[i] = 0.0; */
  ATOL = *atol;
  RTOL = *rtol;
  do_transit_abs = *transit_abs;
  par_ptr = theta;
  HMIN           = hmin;
  HMAX           = hmax;
  H0             = h0;
  MXORDN         = mxordn;
  MXORDS         = mxords;
  mxstep         = mxstep;
  // Counters
  slvr_counter   = 0;
  dadt_counter   = 0;
  jac_counter    = 0;
  // Assign global time information
  all_times     = time; 
  n_all_times   = *ntime;
  RxODE_ode_alloc();
  // Call solver
  RxODE_ode_solver_c(*neq, *stiff, evid, inits, dose, ret, rc);
  if (*nlhs) {
    for (i=0; i<*ntime; i++){
      calc_lhs(time[i], ret+i*(*neq), lhs+i*(*nlhs));
    }
  }
  RxODE_ode_free();
}

void RxODE_assign_fn_pointers(void (*fun_dydt)(unsigned int, double, double *, double *),
                              void (*fun_calc_lhs)(double, double *, double *),
                              void (*fun_calc_jac)(unsigned int, double, double *, double *, unsigned int),
			      void (*fun_update_inis)(SEXP _ini_sexp),
                              int fun_jt,
                              int fun_mf,
                              int fun_debug){
  // Assign functions pointers
  dydt     = fun_dydt;
  calc_jac = fun_calc_jac;
  calc_lhs = fun_calc_lhs;
  update_inis = fun_update_inis;
  // Assign solver options
  global_jt     = fun_jt;
  global_mf     = fun_mf;
  global_debug  = fun_debug;
}

void RxODE_ode_setup(SEXP sexp_inits,
                     SEXP sexp_lhs,
        	     // Events
                     SEXP sexp_time,
                     SEXP sexp_evid,
                     SEXP sexp_dose,
                     // Covariates

                     SEXP sexp_pcov,
                     SEXP sexp_cov,
                     SEXP sexp_locf,
                     // Solver Options
                     SEXP sexp_atol,
                     SEXP sexp_rtol,
                     SEXP sexp_hmin,
                     SEXP sexp_hmax,
                     SEXP sexp_h0,
                     SEXP sexp_mxordn,
                     SEXP sexp_mxords,
                     SEXP sexp_mx,
                     SEXP sexp_stiff,
                     SEXP sexp_transit_abs){
  /* RxODE_ode_free(); */
  // Events
  all_times     = REAL(sexp_time);
  n_all_times   = length(sexp_time);
  evid          = INTEGER(sexp_evid);
  dose          = REAL(sexp_dose);
  RxODE_ode_dosing_ = RxODE_ode_get_dosing();
  // Covariates
  par_cov       = INTEGER(sexp_pcov);
  cov_ptr       = REAL(sexp_cov);
  ncov          = length(sexp_pcov);
  is_locf       = INTEGER(sexp_locf)[0];
  // Solver options
  ATOL           = REAL(sexp_atol)[0];
  RTOL           = REAL(sexp_rtol)[0];
  HMAX           = REAL(sexp_hmax)[0];
  H0             = REAL(sexp_h0)[0];
  MXORDN         = INTEGER(sexp_mxordn)[0];
  MXORDS         = INTEGER(sexp_mxords)[0];
  mxstep         = INTEGER(sexp_mx)[0];
  do_transit_abs = INTEGER(sexp_transit_abs)[0];
  stiff          = INTEGER(sexp_stiff)[0];
  slvr_counter   = 0;
  dadt_counter   = 0;
  jac_counter    = 0;
  rx_aprox_M.f2 = 0.0; //= f=0 
  rx_aprox_M.f1 = 1.0; // = 1-f = 1;
  rx_aprox_M.kind = !is_locf;
  nlhs          = length(sexp_lhs);
  neq           = length(sexp_inits);
  update_inis(sexp_inits); // Update any run-time initial conditions.
}

void RxODE_ode_solve_env(SEXP sexp_rho){
  SEXP sexp_theta = findVar(installChar(mkChar("params")),sexp_rho);
  SEXP sexp_inits = findVar(installChar(mkChar("inits")),sexp_rho);
  SEXP sexp_lhs   = findVar(installChar(mkChar("lhs_vars")),sexp_rho);
  // Events
  SEXP sexp_time = findVar(installChar(mkChar("time")),sexp_rho);
  SEXP sexp_evid = findVar(installChar(mkChar("evid")),sexp_rho);
  SEXP sexp_dose = findVar(installChar(mkChar("amt")),sexp_rho);
  // Covariates
  SEXP sexp_pcov = findVar(installChar(mkChar("pcov")),sexp_rho);
  SEXP sexp_cov = findVar(installChar(mkChar("cov")),sexp_rho);
  SEXP sexp_locf = findVar(installChar(mkChar("isLocf")),sexp_rho);
  // Solver Options
  SEXP sexp_atol = findVar(installChar(mkChar("atol")),sexp_rho);
  SEXP sexp_rtol = findVar(installChar(mkChar("rtol")),sexp_rho);
  SEXP sexp_hmin = findVar(installChar(mkChar("hmin")),sexp_rho);
  SEXP sexp_hmax = findVar(installChar(mkChar("hmax")),sexp_rho);
  SEXP sexp_h0 = findVar(installChar(mkChar("hini")),sexp_rho);
  SEXP sexp_mxordn = findVar(installChar(mkChar("maxordn")),sexp_rho);
  SEXP sexp_mxords = findVar(installChar(mkChar("maxords")),sexp_rho);
  SEXP sexp_mx = findVar(installChar(mkChar("maxsteps")),sexp_rho);
  SEXP sexp_stiff = findVar(installChar(mkChar("stiff")),sexp_rho);
  SEXP sexp_transit_abs = findVar(installChar(mkChar("transit_abs")),sexp_rho);

  par_ptr       = REAL(sexp_theta);
  inits         = REAL(sexp_inits);

  RxODE_ode_setup(sexp_inits, sexp_lhs, sexp_time, sexp_evid, sexp_dose, sexp_pcov, sexp_cov, sexp_locf, sexp_atol, sexp_rtol, sexp_hmin, sexp_hmax, sexp_h0, sexp_mxordn, sexp_mxords, sexp_mx, sexp_stiff, sexp_transit_abs);
  RxODE_ode_alloc();
  
  RxODE_ode_solver_c(neq, stiff, evid, inits, dose, solve, rc);

}

SEXP RxODE_ode_solver (// Parameters
		       SEXP sexp_theta,
		       SEXP sexp_inits,
		       SEXP sexp_lhs,
		       // Events
		       SEXP sexp_time,
		       SEXP sexp_evid,
		       SEXP sexp_dose,
		       // Covariates
		       SEXP sexp_pcov,
		       SEXP sexp_cov,
		       SEXP sexp_locf,
		       // Solver Options
		       SEXP sexp_atol,
		       SEXP sexp_rtol,
		       SEXP sexp_hmin,
		       SEXP sexp_hmax,
		       SEXP sexp_h0,
		       SEXP sexp_mxordn,
		       SEXP sexp_mxords,
		       SEXP sexp_mx,
		       SEXP sexp_stiff,
		       SEXP sexp_transit_abs,
		       // Object Creation
		       SEXP sexp_object,
		       SEXP sexp_extra_args){
  // TODO: Absorption lag?
  // TODO: Annotation?
  // TODO: Units
  // TODO: Bioavailiability
  // Parameters
  par_ptr       = REAL(sexp_theta);
  inits         = REAL(sexp_inits);
  // Events
  RxODE_ode_setup(sexp_inits, sexp_lhs, sexp_time, sexp_evid, sexp_dose, sexp_pcov, sexp_cov, sexp_locf, sexp_atol, sexp_rtol, sexp_hmin, sexp_hmax, sexp_h0, sexp_mxordn, sexp_mxords, sexp_mx, sexp_stiff, sexp_transit_abs);
  RxODE_ode_alloc();
  
  int i = 0, j = 0;
  SEXP sexp_ret     = PROTECT(allocMatrix(REALSXP, nAllTimes(), ncov+1+neq+nlhs));
  SEXP sexp_counter = PROTECT(allocVector(INTSXP,4));
  int    *counts    = INTEGER(sexp_counter);
  
  double *ret   = REAL(sexp_ret);
  
  /* for (i=0; i<neq; i++) InfusionRate[i] = 0.0; */
  
  RxODE_ode_solver_c(neq, stiff, evid, inits, dose, solve, rc);

  // Now create the matrix.
  for (i = 0; i < nAllTimes(); i++){
    // Time
    ret[i] = all_times[i];
    // State
    if (neq){
      for (j = 0; j < neq; j++){
      	ret[nAllTimes()*(j+1)+i] = solve[j+i*neq];
      }
    }
    // LHS
    if (nlhs){
      rxCalcLhs(i);
      for (j = 0; j < nlhs; j++){
	ret[nAllTimes()*(j+1+neq)+i] = rxLhs(j);
      }
    }
    // Cov
    if (ncov > 0){
      for (j = 0; j < ncov; j++){
      	ret[nAllTimes()*(j+1+neq+nlhs)+i] = cov_ptr[j*nAllTimes()+i];
      }
    }
  }
  SEXP sexp_dimnames = PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(sexp_dimnames, 0, R_NilValue);
  SEXP sexp_colnames = PROTECT(allocVector(STRSXP,1+neq+nlhs+ncov));
  SEXP sexp_solve    = PROTECT(allocVector(REALSXP,1+neq+nlhs+ncov));
  SET_STRING_ELT(sexp_colnames, 0, mkChar("time"));
  double *solver = REAL(sexp_solve);
  solver[0] = all_times[nAllTimes()-1];
  SEXP temp = getAttrib(sexp_inits, R_NamesSymbol);
  for (i = 0; i < neq; i++){
  SET_STRING_ELT(sexp_colnames, 1+i, STRING_ELT(temp,i));
    solver[1+i] = inits[i];
  }
  for (i = 0; i < nlhs; i++){
    SET_STRING_ELT(sexp_colnames,1+neq+i, STRING_ELT(sexp_lhs,i));
    solver[1+neq+i] = NA_REAL;
  }
  temp = getAttrib(sexp_theta,R_NamesSymbol);
  for (i = 0; i < ncov; i++){
    SET_STRING_ELT(sexp_colnames,1+neq+nlhs+i, STRING_ELT(temp, par_cov[i]-1));
    solver[1+neq+nlhs+i] = NA_REAL;
  }
  SET_VECTOR_ELT(sexp_dimnames,1,sexp_colnames);
  setAttrib(sexp_ret, R_DimNamesSymbol, sexp_dimnames);
  SEXP sexp_ncounter = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(sexp_ncounter, 0, mkChar("solver"));
  counts[0] = slvr_counter;
  SET_STRING_ELT(sexp_ncounter, 1, mkChar("dadt"));
  counts[1] = dadt_counter;
  SET_STRING_ELT(sexp_ncounter, 2, mkChar("user_jac"));
  counts[2] = jac_counter;
  SET_STRING_ELT(sexp_ncounter, 3, mkChar("rc"));
  counts[3] = rc[0];
  setAttrib(sexp_counter, R_NamesSymbol, sexp_ncounter);
  
  SEXP sexp_lst   = PROTECT(allocVector(VECSXP, 5 + length(sexp_extra_args)));
  SEXP sexp_lstn  = PROTECT(allocVector(STRSXP, 5 + length(sexp_extra_args)));

  SET_STRING_ELT(sexp_lstn,0,mkChar("matrix"));
  SET_VECTOR_ELT(sexp_lst, 0,sexp_ret);
  
  SET_STRING_ELT(sexp_lstn,1,mkChar("counts"));
  SET_VECTOR_ELT(sexp_lst, 1,sexp_counter);
  
  SET_STRING_ELT(sexp_lstn,2,mkChar("inits"));
  SET_VECTOR_ELT(sexp_lst, 2,sexp_inits);
  
  SET_STRING_ELT(sexp_lstn,3,mkChar("params"));
  SET_VECTOR_ELT(sexp_lst, 3,sexp_theta);
  
  SET_STRING_ELT(sexp_lstn,4,mkChar("object"));
  SET_VECTOR_ELT(sexp_lst, 4,sexp_object);

  temp = getAttrib(sexp_extra_args,R_NamesSymbol);
  for (i = 0; i < length(sexp_extra_args); i++){
    SET_STRING_ELT(sexp_lstn,5+i,STRING_ELT(temp,i));
    SET_VECTOR_ELT(sexp_lst, 5+i,VECTOR_ELT(sexp_extra_args,i));
  }
  
  setAttrib(sexp_lst, R_NamesSymbol, sexp_lstn);
  setAttrib(sexp_solve, install("solveRxDll"), sexp_lst);
  SEXP sexp_class = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(sexp_class, 0, mkChar("solveRxDll"));
  setAttrib(sexp_solve, R_ClassSymbol, sexp_class);
  setAttrib(sexp_solve, R_NamesSymbol, sexp_colnames);

  UNPROTECT(9);
  if (fp) fclose(fp);
  RxODE_ode_free();
  return sexp_solve;
}

double RxODE_InfusionRate(int val){
  return InfusionRate[val];
}

double RxODE_par_ptr(int val){
  return par_ptr[val];
}

long RxODE_jac_counter_val(){
  return jac_counter;
}

long RxODE_dadt_counter_val(){
  return dadt_counter;
}

void RxODE_jac_counter_inc(){
  jac_counter++;
}

void RxODE_dadt_counter_inc(){
  dadt_counter++;
}

double RxODE_podo(){
  return podo;
}

double RxODE_tlast(){
  return tlast;
}

double RxODE_transit4(double t, double n, double mtt, double bio){
  double ktr = (n+1)/mtt;
  double lktr = log(n+1)-log(mtt);
  return exp(log(bio*podo)+lktr+n*(lktr+log(t))-ktr*t-lgamma1p(n));
}

double RxODE_transit3(double t, double n, double mtt){
  return RxODE_transit4(t, n,mtt, 1.0);
}

double RxODE_factorial(double x){
  return exp(lgamma1p(x));
}

SEXP RxODE_rxInner(SEXP etanewsSEXP, SEXP rhoSEXP);
SEXP RxODE_rxHessian(SEXP rhoSEXP);
SEXP RxODE_RxODE_focei_eta_lik(SEXP sexp_etaSEXP, SEXP sexp_rhoSEXP);
SEXP RxODE_RxODE_focei_eta_lp(SEXP sexp_etaSEXP, SEXP sexp_rhoSEXP);
SEXP RxODE_RxODE_focei_eta(SEXP fstrSEXP);
SEXP RxODE_RxODE_focei_finalize_llik(SEXP rhoSEXP);
SEXP RxODE_RxODE_finalize_log_det_OMGAinv_5(SEXP rhoSEXP);
SEXP RxODE_rxDetaDomega(SEXP rhoSEXP);
SEXP RxODE_rxOuter_(SEXP rhoSEXP);
SEXP RxODE_rxDetaDtheta(SEXP rhoSEXP);
SEXP RxODE_rxOuter(SEXP rhoSEXP);
SEXP RxODE_rxUpdateEtas(SEXP DnDhSSEXP, SEXP DhSSEXP, SEXP initSSEXP, SEXP acceptNSSEXP);
SEXP RxODE_RxODE_finalize_focei_omega(SEXP rho);
SEXP trans(SEXP orig_file, SEXP parse_file, SEXP c_file, SEXP extra_c, SEXP prefix, SEXP model_md5, SEXP parse_model,SEXP parse_model3);
SEXP RxODE_getMacroConstants(SEXP rho);
SEXP RxODE_getLinDerivs(SEXP rho);
SEXP RxODE_linCmtEnv(SEXP rho);
SEXP RxODE_rxInv(SEXP matrix);

double solvedC(double t, int parameterization, unsigned int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);

void R_init_RxODE(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"RxODE_ode_solver", (DL_FUNC) &RxODE_ode_solver, 21},
    {"RxODE_rxInner", (DL_FUNC) &RxODE_rxInner, 2},
    {"RxODE_rxHessian", (DL_FUNC) &RxODE_rxHessian, 1},
    {"RxODE_RxODE_focei_eta_lik", (DL_FUNC) &RxODE_RxODE_focei_eta_lik, 2},
    {"RxODE_RxODE_focei_eta_lp", (DL_FUNC) &RxODE_RxODE_focei_eta_lp, 2},
    {"RxODE_RxODE_focei_eta", (DL_FUNC) &RxODE_RxODE_focei_eta, 1},
    {"RxODE_RxODE_focei_finalize_llik", (DL_FUNC) &RxODE_RxODE_focei_finalize_llik, 1},
    {"RxODE_RxODE_finalize_log_det_OMGAinv_5", (DL_FUNC) &RxODE_RxODE_finalize_log_det_OMGAinv_5, 1},
    {"RxODE_RxODE_finalize_focei_omega", (DL_FUNC) &RxODE_RxODE_finalize_focei_omega, 1},
    {"RxODE_rxDetaDomega", (DL_FUNC) &RxODE_rxDetaDomega, 1},
    {"RxODE_rxOuter_", (DL_FUNC) &RxODE_rxOuter_, 1},
    {"RxODE_rxDetaDtheta", (DL_FUNC) &RxODE_rxDetaDtheta, 1},
    {"RxODE_rxOuter", (DL_FUNC) &RxODE_rxOuter, 1},
    {"RxODE_rxUpdateEtas", (DL_FUNC) &RxODE_rxUpdateEtas, 4},
    {"trans", (DL_FUNC) &trans, 8},
    {"RxODE_getMacroConstants", (DL_FUNC) &RxODE_getMacroConstants, 1},
    {"RxODE_getLinDerivs", (DL_FUNC) &RxODE_getLinDerivs, 1},
    {"RxODE_linCmtEnv", (DL_FUNC) &RxODE_linCmtEnv, 1},
    {"RxODE_rxInv", (DL_FUNC) &RxODE_rxInv, 1},
    {NULL, NULL, 0}
  };
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  
  //Functions
  R_RegisterCCallable("RxODE","RxODE_ode_solver",       (DL_FUNC) RxODE_ode_solver);
  R_RegisterCCallable("RxODE","RxODE_assign_fn_pointers", (DL_FUNC) RxODE_assign_fn_pointers);
  R_RegisterCCallable("RxODE","RxODE_ode_solver_old_c", (DL_FUNC) RxODE_ode_solver_old_c);
  R_RegisterCCallable("RxODE","RxODE_ode_solver_0_6_c", (DL_FUNC) RxODE_ode_solver_0_6_c);
  R_RegisterCCallable("RxODE","RxODE_ode_setup",         (DL_FUNC) RxODE_ode_setup);
  R_RegisterCCallable("RxODE","RxODE_ode_free", (DL_FUNC) RxODE_ode_free);
  
  //Infusion
  R_RegisterCCallable("RxODE","RxODE_InfusionRate",     (DL_FUNC) RxODE_InfusionRate);
  // Parameters
  R_RegisterCCallable("RxODE","RxODE_par_ptr",          (DL_FUNC) RxODE_par_ptr);
  R_RegisterCCallable("RxODE","RxODE_update_par_ptr",   (DL_FUNC) update_par_ptr);
  // Counters
  R_RegisterCCallable("RxODE","RxODE_dadt_counter_val", (DL_FUNC) RxODE_dadt_counter_val);
  R_RegisterCCallable("RxODE","RxODE_jac_counter_val",  (DL_FUNC) RxODE_jac_counter_val);
  R_RegisterCCallable("RxODE","RxODE_dadt_counter_inc", (DL_FUNC) RxODE_dadt_counter_inc);
  R_RegisterCCallable("RxODE","RxODE_jac_counter_inc",  (DL_FUNC) RxODE_jac_counter_inc);
  // podo or tlast
  R_RegisterCCallable("RxODE","RxODE_podo",             (DL_FUNC) RxODE_podo);
  R_RegisterCCallable("RxODE","RxODE_tlast",            (DL_FUNC) RxODE_tlast);
  // tranit compartment models
  R_RegisterCCallable("RxODE","RxODE_transit4",         (DL_FUNC) RxODE_transit4);
  R_RegisterCCallable("RxODE","RxODE_transit3",         (DL_FUNC) RxODE_transit3);
  R_RegisterCCallable("RxODE","RxODE_factorial",        (DL_FUNC) RxODE_factorial);
  R_RegisterCCallable("RxODE","RxODE_safe_log",         (DL_FUNC) RxODE_safe_log);
  R_RegisterCCallable("RxODE","RxODE_safe_zero",        (DL_FUNC) RxODE_safe_zero);
  R_RegisterCCallable("RxODE","RxODE_as_zero",          (DL_FUNC) RxODE_as_zero);
  R_RegisterCCallable("RxODE","solvedC",                (DL_FUNC) solvedC);
}
