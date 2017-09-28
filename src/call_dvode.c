#include <stdio.h>
#include <stdlib.h>
#include "dop853.h"
#define NCMT 100
#define max(a, b) ((a) > (b) ? (a) : (b))
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
double InfusionRate[NCMT];
int BadDose[NCMT];
int nBadDose=0;
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
int   do_par_cov=1;
int rc[1];
double *cov_ptr;
int    ncov, nlhs, neq, stiff;
int    is_locf;
int    n_all_times;
int    mxstep;
int    MXORDN;
int    MXORDS;
int    global_jt, global_mf, global_debug, ixds,ndoses = -1;
double *all_times;
int *idose;
int idosen = 0;
int extraCmt = 0;
FILE *fp;

/* void __DYDT__(unsigned int neq, double t, double *A, double *DADT); */
/* void __CALC_LHS__(double t, double *A, double *lhs); */
/* void __CALC_JAC__(unsigned int neq, double t, double *A, double *JAC, unsigned int __NROWPD__); */

void (*dydt)(unsigned int neq, double t, double *A, double *DADT);
void (*calc_jac)(unsigned int neq, double t, double *A, double *JAC, unsigned int __NROWPD__);
void (*calc_lhs)(double t, double *A, double *lhs);
void (*update_inis)(SEXP _ini_sexp);

void setExtraCmt(int xtra){
  if (xtra > extraCmt){
    extraCmt = xtra;
  }
}

double rxTime(int i){
  if (i < n_all_times){
    return(all_times[i]);
  } else {
    error("Time cannot be retrived (%dth entry).",i);
  }
  return 0;
}

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

unsigned int nDoses(){
  if (ndoses < 0){
    ndoses=0;
    for (int i = 0; i < n_all_times; i++){
      if (rxEvid(i)){
        ndoses++;
	if (ndoses >= idosen){
	  if (idosen == 0){
	    idose = Calloc(32,int);
	    idosen = 32;
	  } else {
	    idosen *= 2;
            Rprintf("Reallocating to %d\n", idosen);
            idose = Realloc(idose, idosen, int);
	  }
	}
	idose[ndoses-1] = i;
      }
    }
    return ndoses;
  } else {
    return ndoses;
  }
}

double rxDosingTime(int i){
  if (i < nDoses()){
    return (all_times[idose[i]]);
  } else {
    error("Dosing cannot retreived (%dth dose).", i);
  }
  return 0;
}

int rxDosingEvid(int i){
  if (i < nDoses()){
    return (evid[idose[i]]);
  } else {
    error("Dosing cannot retreived (%dth dose).", i);
  }
  return 0;
}

double rxDose(int i){
  if (i < nDoses()){
    return(dose[i]);
  } else {
    error("Dose cannot be retrived (%dth entry).",i);
  }
  return 0;
}



unsigned int nObs(){
  return (unsigned int)(n_all_times - nDoses());
}

unsigned int nLhs (){
  return (unsigned int)(nlhs);
}

int nEq (){
  return neq;
}

double RxODE_as_zero(double x){
  if (fabs(x) < sqrt(DOUBLE_EPS)){
    return(0.0);
  } else {
    return(x);
  }
}

extern double RxODE_safe_log(double x){
  if (x <= 0){
    // Warning?
    return log(DOUBLE_EPS);
  } else {
    return log(x);
  }
}

double RxODE_safe_zero(double x){
  if (x == 0){
    // Warning?
    return DOUBLE_EPS;
  } else {
    return(x);
  }
}

double RxODE_pow(double x, double y){
  if (x == 0 && y <= 0){
    return R_pow(DOUBLE_EPS, y);
  } else {
    return R_pow(x, y);
  }
}
double RxODE_pow_di(double x, int i){
  if (x == 0 && i <= 0){
    return R_pow_di(DOUBLE_EPS, i);
  } else {
    return R_pow_di(x, i);
  }
}

double RxODE_sign_exp(double sgn, double x){
  if (sgn > 0.0){
    return(exp(x));
  } else if (sgn < 0.0){
    return(-exp(x));
  } else {
    return(0.0);
  }
}

static R_NativePrimitiveArgType RxODE_sign_exp_t[] = {
  REALSXP, REALSXP
};

double RxODE_abs_log(double x){
  if  (fabs(x) <= sqrt(DOUBLE_EPS)){
    return log(sqrt(DOUBLE_EPS));
  } else if (x > 0.0){
    return log(x);
  } else if (x < 0.0){
    return log(-x);
  } else {
    return 0.0;
  }
}

double RxODE_abs_log1p(double x){
  if (x + 1.0 > 0.0){
    return(log1p(x));
  } else if (x + 1.0 > 0.0){
    return(log1p(-x));
  } else {
    return 0.0;
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
  if (do_par_cov){
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
  int i, j, foundBad;
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
  rwork = (double*)Calloc(lrw+1, double);
  iwork = (int*)Calloc(liw+1, int);
  
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
	      Rprintf("Error solving using LSODA\n");
	      return;
	    }
	  slvr_counter++;
	  //dadt_counter = 0;
	}
      if (wh)
	{
	  cmt = (wh%10000)/100 - 1;
	  if (cmt >= nEq()){
	    foundBad = 0;
            for (j = 0; j <nBadDose; j++){
	      if (BadDose[j] == cmt+1){
		foundBad=1;
		break;
	      }
	    }
	    if (!foundBad){
	      BadDose[nBadDose]=cmt+1;
	      nBadDose++;
	    }
	  } else {
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
                else yp[cmt] += dose[ixds];     //dosing before obs
              }
	    
	    istate = 1;

	    ixds++;
	    xp = xout;
	  }
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
  Free(rwork);
  Free(iwork);
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
  int wh, cmt, foundBad=0;
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
	      Rprintf("Error sovling using dop853");
	      return;  //exit(1);  // dj: should not abort R
	    }

	  xp = xRead();
	  slvr_counter++;
	  //dadt_counter = 0;
	}
      if (wh)
	{
	  cmt = (wh%10000)/100 - 1;
	  if (cmt >= nEq()){
            foundBad = 0;
            for (j = 0; j <nBadDose; j++){
              if (BadDose[j] == cmt+1){
                foundBad=1;
                break;
              }
            }
            if (!foundBad){
              BadDose[nBadDose]=cmt+1;
              nBadDose++;
            }
          } else {
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
    if (neq > NCMT){
      error("RxODE does not support %d compartments (Currently only %d compartments)", neq, NCMT);
    }
    if (stiff==0){
      call_dop(neq, all_times, evid, nAllTimes(), inits, dose, solve, rc);
    } else{
      call_lsoda(neq, all_times, evid, nAllTimes(), inits, dose, solve, rc);
    }
  }
}

void RxODE_ode_solver_old_c(int *neqa,
                            double *theta,  //order:
                            double *time,
                            int *evidp,
                            int *ntime,
                            double *initsp,
                            double *dosep,
                            double *ret,
                            double *atol,
                            double *rtol,
                            int *stiffa,
                            int *transit_abs,
                            int *nlhsa,
                            double *lhsp,
                            int *rc){
  if (*neqa > NCMT){
    error("RxODE does not support %d compartments (Currently only %d compartments)", neq, NCMT);
  }
  int i;
  for (i=0; i< *neqa; i++) InfusionRate[i] = 0.0;
  ndoses = -1;
  all_times         = time;
  n_all_times       = *ntime;
  //RxODE_ode_dosing_ = RxODE_ode_get_dosing();
  // par_cov
  do_par_cov        = 0;
  // cov_ptr
  ncov              = 0;
  is_locf           = 0;
  // Solver Options
  ATOL = *atol;
  RTOL = *rtol;
  // Assign to default LSODA behvior, or 0
  HMIN           = 0;
  HMAX           = 0;
  H0             = 0;
  MXORDN         = 0;
  MXORDS         = 0;
  mxstep         = 5000; // Not LSODA default but RxODE default
  // Counters
  slvr_counter   = 0;
  dadt_counter   = 0;
  jac_counter    = 0;

  nlhs           = *nlhsa;
  neq            = *neqa;
  stiff          = *stiffa;
  
  

  nBadDose = 0;
  do_transit_abs = *transit_abs;
  
  par_ptr = theta;
  inits   = initsp;
  dose    = dosep;
  solve   = ret;
  lhs     = lhsp;
  evid    = evidp;

  // Assign global time information
  // Call solver
  
  
  /* Rprintf("Call Solver; par_ptr[0] = %f; evid[0]=%d; inits[0]=%d\n",par_ptr[0],evid[0],inits[0]); */
  RxODE_ode_solver_c(*neqa, *stiffa, evidp, initsp, dosep, ret, rc);
  /* Rprintf("Update LHS\n"); */
  // Update LHS
  if (*nlhsa) {
    for (i=0; i<*ntime; i++){
      calc_lhs(time[i], ret+i*(*neqa), lhsp+i*(*nlhsa));
    }
  }
}

/* SEXP RxODE_ode_get_dosing(){ */
/*   SEXP sexp_ret; */
/*   int  i, j = 0; */
/*   sexp_ret  = PROTECT(allocMatrix(REALSXP, nDoses(), 3)); */
/*   double *ret   = REAL(sexp_ret); */
/*   for (i = 0; i < n_all_times; i++){ */
/*     if (rxEvid(i)){ */
/*       ret[j] = all_times[i]; */
/*       ret[nDoses()+j] =(double)(evid[i]); */
/*       ret[nDoses()*2+j] = dose[j]; */
/*       /\* Rprintf("\tj:%d(%d)\tt: %f\tevid:%f\tdose:%f\n", *\/ */
/*       /\* 	      j, ndoses, *\/ */
/*       /\* 	      ret[j],ret[ndoses+j],ret[ndoses*2+j]); *\/ */
/*       j++; */
/*     } */
/*   } */
/*   // Unprotect after solving. */
/*   UNPROTECT(1); */
/*   return sexp_ret; */
/* } */

void RxODE_ode_free(){
  /* Free(InfusionRate); */
  int j;
  if (nBadDose){
    for (j=0; j < nBadDose; j++){
      if (BadDose[j] > extraCmt){
	warning("Dose to Compartment %d ignored (not in ODE)",BadDose[j]);
      }
    }
  }
  Free(solve);
  Free(lhs);
  Free(idose);
  idosen       = 0;
  extraCmt     = 0;
  /* Free(rc); */
}

void RxODE_ode_alloc(){
  if (neq > NCMT){
    error("RxODE does not support %d compartments (Currently only %d compartments)", neq, NCMT);
  }
  /* RxODE_ode_dosing_calc = 0; */
  /* solve = (double*)  R_alloc(neq*n_all_times+1, sizeof(double)); */
  /* lhs   = (double*)  R_alloc(nlhs,sizeof(double)); */
  /* InfusionRate = (double *) R_alloc(neq+2,sizeof(double)); */
  /* rc           = (int *)    R_alloc(1,sizeof(int)); */
  solve        = Calloc(neq*nAllTimes()+1,double);
  lhs          = Calloc(nlhs,double);
  /* InfusionRate = (double *) Calloc(neq+2,double); */
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
  if (*neq > NCMT){
    error("RxODE does not support %d compartments (Currently only %d compartments)", neq, NCMT);
  }
  int i;
  for (i=0; i<*neq; i++) InfusionRate[i] = 0.0;
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
  ndoses        = -1;
  all_times     = REAL(sexp_time);
  n_all_times   = length(sexp_time);
  evid          = INTEGER(sexp_evid);
  dose          = REAL(sexp_dose);
  // Covariates
  par_cov       = INTEGER(sexp_pcov);
  do_par_cov    = 1;
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
  nBadDose = 0;
  if (neq > NCMT){
    error("RxODE does not support %d compartments (Currently only %d compartments)", neq, NCMT);
  }
  if (length(sexp_inits) > 0){
    update_inis(sexp_inits); // Update any run-time initial conditions.
  }
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
  SEXP sexp_rc = findVar(installChar(mkChar("rc")),sexp_rho);
  int *rce    = INTEGER(sexp_rc);

  par_ptr       = REAL(sexp_theta);
  inits         = REAL(sexp_inits);

  RxODE_ode_setup(sexp_inits, sexp_lhs, sexp_time, sexp_evid, sexp_dose, sexp_pcov, sexp_cov, sexp_locf, sexp_atol, sexp_rtol, sexp_hmin, sexp_hmax, sexp_h0, sexp_mxordn, sexp_mxords, sexp_mx, sexp_stiff, sexp_transit_abs);
  RxODE_ode_alloc();
  RxODE_ode_solver_c(neq, stiff, evid, inits, dose, solve, rc);
  // Send rc to environment
  rce[0] = rc[0];
}

SEXP RxODE_ode_solver (// Parameters
		       SEXP sexp_theta,
		       SEXP sexp_inits,
		       SEXP sexp_scale,
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
		       SEXP sexp_extra_args,
		       SEXP sexp_matrix,
		       SEXP sexp_add_cov,
		       SEXP sexp_ignore_state){
  // TODO: Absorption lag?
  // TODO: Annotation? -- in nlmixr
  // TODO: Units -- Can be in nlmixr
  // TODO: Bioavailiability?
  // Parameters
  int pro=0;
  par_ptr       = REAL(sexp_theta);
  inits         = REAL(sexp_inits);
  // Events
  RxODE_ode_setup(sexp_inits, sexp_lhs, sexp_time, sexp_evid, sexp_dose, sexp_pcov, sexp_cov, sexp_locf, sexp_atol, sexp_rtol, sexp_hmin, sexp_hmax, sexp_h0, sexp_mxordn, sexp_mxords, sexp_mx, sexp_stiff, sexp_transit_abs);
  RxODE_ode_alloc();
  
  
  int i = 0, j = 0;
  int *rmState = INTEGER(sexp_ignore_state);
  int nPrnState =0;
  for (i = 0; i < length(sexp_ignore_state); i++){
    nPrnState+= (1-rmState[i]);
  }
  SEXP sexp_counter = PROTECT(allocVector(INTSXP,4));pro++;
  int    *counts    = INTEGER(sexp_counter);

  int *matrix = INTEGER(sexp_matrix);
  
  for (i=0; i<neq; i++) InfusionRate[i] = 0.0;
  if (neq > 0){
    RxODE_ode_solver_c(neq, stiff, evid, inits, dose, solve, rc);
  }
  double *scale = REAL(sexp_scale);
  int add_cov = INTEGER(sexp_add_cov)[0];
  if (matrix[0]){
    SEXP sexp_ret     = PROTECT(allocMatrix(REALSXP, nObs(), add_cov*ncov+1+nPrnState+nlhs)); pro++;
    double *ret   = REAL(sexp_ret);
  
    // Now create the matrix.
    int ii=0, jj;
    for (i = 0; i < nAllTimes(); i++){
      // Time
      if (rxEvid(i)==0){
        ret[ii] = all_times[i];
        // State
        if (nPrnState){
	  jj = 0;
          for (j = 0; j < neq; j++){
	    if (rmState[j] == 0){
	      ret[nObs()*(jj+1)+ii] = solve[j+i*neq]/scale[j];
	      jj++;
            }
          }
        }
        // LHS
        if (nlhs){
          rxCalcLhs(i);
          for (j = 0; j < nlhs; j++){
            ret[nObs()*(j+1+nPrnState)+ii] = rxLhs(j);
          }
        }
        // Cov
        if (add_cov*ncov > 0){
          for (j = 0; j < ncov; j++){
            ret[nObs()*(j+1+nPrnState+nlhs)+ii] = cov_ptr[j*nAllTimes()+i];
          }
        }
        ii++;
      }
    }
    SEXP sexp_dimnames = PROTECT(allocVector(VECSXP,2));pro++;
    SET_VECTOR_ELT(sexp_dimnames, 0, R_NilValue);
    SEXP sexp_colnames = PROTECT(allocVector(STRSXP,1+nPrnState+nlhs+add_cov*ncov)); pro++;
    SET_STRING_ELT(sexp_colnames, 0, mkChar("time"));
    SEXP temp = getAttrib(sexp_inits, R_NamesSymbol);
    ii = 0;
    for (i = 0; i < neq; i++){
      if (!rmState[i]){
	SET_STRING_ELT(sexp_colnames, 1+ii, STRING_ELT(temp,i));
	ii++;
      }
    }
    for (i = 0; i < nlhs; i++){
      SET_STRING_ELT(sexp_colnames,1+nPrnState+i, STRING_ELT(sexp_lhs,i));
    }
    temp = getAttrib(sexp_theta,R_NamesSymbol);
    for (i = 0; i < add_cov*ncov; i++){
      SET_STRING_ELT(sexp_colnames,1+nPrnState+nlhs+i, STRING_ELT(temp, par_cov[i]-1));
    }
    SET_VECTOR_ELT(sexp_dimnames,1,sexp_colnames);
    setAttrib(sexp_ret, R_DimNamesSymbol, sexp_dimnames);
    SEXP sexp_solve2   = PROTECT(allocVector(VECSXP, 2)); pro++;
    SEXP sexp_rc = PROTECT(allocVector(INTSXP,1)); pro++;
    int *rc2 = INTEGER(sexp_rc);
    rc2[0] = rc[0];
    SET_VECTOR_ELT(sexp_solve2, 0,sexp_ret);
    SET_VECTOR_ELT(sexp_solve2, 1,sexp_rc);
    UNPROTECT(pro);
    if (fp) fclose(fp);
    RxODE_ode_free();
    return sexp_solve2;
  } else {
    int ncols =add_cov*ncov+1+nPrnState+nlhs,
      nobs =nObs(),
      ntimes = nAllTimes();
    SEXP df = PROTECT(allocVector(VECSXP,ncols)); pro++;
    for (i = 0; i < ncols; i++){
      SET_VECTOR_ELT(df, i, PROTECT(allocVector(REALSXP, nobs))); pro++;
    }
    // Now create the matrix.
    double *dfp;
    int ii=0, jj = 0;
    for (i = 0; i < ntimes; i++){
      // Time
      if (rxEvid(i)==0){
	dfp = REAL(VECTOR_ELT(df, 0));
        dfp[ii] = all_times[i];
        // State
        if (nPrnState){
	  jj = 0;
          for (j = 0; j < neq; j++){
	    if (!rmState[j]){
	      dfp = REAL(VECTOR_ELT(df, jj+1));
              dfp[ii] = solve[j+i*neq]/scale[j];
	      jj++;
            }
          }
        }
        // LHS
        if (nlhs){
          rxCalcLhs(i);
          for (j = 0; j < nlhs; j++){
	    dfp = REAL(VECTOR_ELT(df, j+1+nPrnState));
	    dfp[ii] =rxLhs(j);
          }
        }
        // Cov
        if (add_cov*ncov > 0){
          for (j = 0; j < add_cov*ncov; j++){
	    dfp = REAL(VECTOR_ELT(df, j+1+nPrnState+nlhs));
            dfp[ii] = cov_ptr[j*nAllTimes()+i];
          }
        }
        ii++;
      }
    }
    SEXP sexp_solve    = PROTECT(allocVector(REALSXP,1+nPrnState+nlhs+ncov*add_cov)); pro++;
    SEXP sexp_colnames = PROTECT(allocVector(STRSXP,1+nPrnState+nlhs+ncov*add_cov)); pro++;
    SEXP sexp_rownames = PROTECT(allocVector(INTSXP,2)); pro++;
    SEXP temp = getAttrib(sexp_inits, R_NamesSymbol);
    SET_STRING_ELT(sexp_colnames, 0, mkChar("time"));
    double *solver = REAL(sexp_solve);
    solver[0] = all_times[ntimes-1];
    ii = 0;
    for (i = 0; i < neq; i++){
      if (!rmState[i]){
	SET_STRING_ELT(sexp_colnames, 1+ii, STRING_ELT(temp,i));
        solver[1+ii] = inits[i];
	ii++;
      }
    }
    for (i = 0; i < nlhs; i++){
      SET_STRING_ELT(sexp_colnames,1+nPrnState+i, STRING_ELT(sexp_lhs,i));
      solver[1+nPrnState+i] = NA_REAL;
    }
    temp = getAttrib(sexp_theta,R_NamesSymbol);
    for (i = 0; i < ncov*add_cov; i++){
      SET_STRING_ELT(sexp_colnames,1+nPrnState+nlhs+i, STRING_ELT(temp, par_cov[i]-1));
      solver[1+nPrnState+nlhs+i] = NA_REAL;
    }
    /* SET_VECTOR_ELT(df,1,sexp_colnames); */
    INTEGER(sexp_rownames)[0] = NA_INTEGER;
    INTEGER(sexp_rownames)[1] = -nobs;
    
    setAttrib(df, R_RowNamesSymbol, sexp_rownames);
    setAttrib(df, R_NamesSymbol, sexp_colnames);
    SEXP sexp_ncounter = PROTECT(allocVector(STRSXP, 4)); pro++;
    SET_STRING_ELT(sexp_ncounter, 0, mkChar("solver"));
    counts[0] = slvr_counter;
    SET_STRING_ELT(sexp_ncounter, 1, mkChar("dadt"));
    counts[1] = dadt_counter;
    SET_STRING_ELT(sexp_ncounter, 2, mkChar("user_jac"));
    counts[2] = jac_counter;
    SET_STRING_ELT(sexp_ncounter, 3, mkChar("rc"));
    counts[3] = rc[0];
    setAttrib(sexp_counter, R_NamesSymbol, sexp_ncounter);

    SEXP env = eval(lang1(install("new.env")),R_GlobalEnv);

    defineVar(install("counts"), sexp_counter, env);
    defineVar(install("inits"), sexp_inits, env);
    defineVar(install("params"), sexp_theta, env);

    defineVar(install("lhs_vars"), sexp_lhs, env);
    // Events
    defineVar(install("time"), sexp_time, env);
    defineVar(install("evid"), sexp_evid, env);
    defineVar(install("amt"), sexp_dose, env);
    // Covariates
    defineVar(install("pcov"), sexp_pcov, env);
    defineVar(install("cov"), sexp_cov, env);
    defineVar(install("isLocf"), sexp_locf, env);
    // Solver Options
    defineVar(install("atol"), sexp_atol, env);
    defineVar(install("rtol"), sexp_rtol, env);
    defineVar(install("hmin"), sexp_hmin, env);
    defineVar(install("hmax"), sexp_hmax, env);
    defineVar(install("hini"), sexp_h0, env);
    defineVar(install("maxordn"), sexp_mxordn, env);
    defineVar(install("maxords"), sexp_mxords, env);
    defineVar(install("maxsteps"), sexp_mx, env);
    defineVar(install("stiff"), sexp_stiff, env);
    defineVar(install("transit_abs"), sexp_transit_abs, env);
    defineVar(install("env"), sexp_object, env);
    defineVar(install("extra.args"), sexp_extra_args, env);
    /* defineVar(install("rc"), sexp_rc, env); */
    
    setAttrib(df,install(".env"), env);

    SEXP cls = PROTECT(allocVector(STRSXP, 2)); pro++;
    SET_STRING_ELT(cls, 0, mkChar("solveRxODE"));
    SET_STRING_ELT(cls, 1, mkChar("data.frame"));

    classgets(df, cls);
    
    SEXP sexp_solve2   = PROTECT(allocVector(VECSXP, 2)); pro++;
    SEXP sexp_rc = PROTECT(allocVector(INTSXP,1)); pro++;
    int *rc2 = INTEGER(sexp_rc);
    rc2[0] = rc[0];
    SET_VECTOR_ELT(sexp_solve2, 0,df);
    SET_VECTOR_ELT(sexp_solve2, 1,sexp_rc);

    UNPROTECT(pro);
    if (fp) fclose(fp);
    RxODE_ode_free();
    return sexp_solve2;
  }
}

double RxODE_InfusionRate(int val){
  return InfusionRate[val];
}

static R_NativePrimitiveArgType RxODE_one_int_t[] = {
  INTSXP
};

double RxODE_par_ptr(int val){
  double ret =par_ptr[val];
  return ret;
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

static R_NativePrimitiveArgType RxODE_transit4_t[] = {
  REALSXP, REALSXP, REALSXP, REALSXP
};

double RxODE_transit3(double t, double n, double mtt){
  return RxODE_transit4(t, n,mtt, 1.0);
}

static R_NativePrimitiveArgType RxODE_transit3_t[] = {
  REALSXP, REALSXP, REALSXP
};

double RxODE_factorial(double x){
  return exp(lgamma1p(x));
}

static R_NativePrimitiveArgType RxODE_one_dbl_t[] = {
  REALSXP
};

SEXP trans(SEXP orig_file, SEXP parse_file, SEXP c_file, SEXP extra_c, SEXP prefix, SEXP model_md5, SEXP parse_model,SEXP parse_model3);
SEXP _RxODE_linCmtEnv(SEXP rho);
SEXP _RxODE_rxInv(SEXP matrix);
SEXP _RxODE_removableDrive(SEXP letter);
SEXP _RxODE_rxCoutEcho(SEXP number);
SEXP _RxODE_W_Cpp(SEXP zSEXP, SEXP branchSEXP);
SEXP _RxODE_RxODE_finalize_focei_omega(SEXP);
SEXP _RxODE_RxODE_finalize_log_det_OMGAinv_5(SEXP);

double RxODE_solveLinB(double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);
static R_NativePrimitiveArgType RxODE_solveLinB_t[] = {
  //t,    linCmt,  diff1,  diff2,  A,       alpha,  B,       beta,     C,       gamma, double ka, double tlag)
  REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

SEXP _rxKahanSum(SEXP input);
SEXP _rxNeumaierSum(SEXP input);

SEXP _rxPythonSum(SEXP input);
SEXP _rxProd(SEXP input);
SEXP _rxSum(SEXP input);
SEXP _rxSetSum(SEXP input);
SEXP _rxSetProd(SEXP input);
SEXP _rxPairwiseSum(SEXP input);

double RxODE_sum(double *input, int len);
double RxODE_prod(double *input, int len);

static R_NativePrimitiveArgType RxODE_Sum_t[] = {
  REALSXP, INTSXP
};


void R_init_RxODE(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"RxODE_ode_solver", (DL_FUNC) &RxODE_ode_solver, 25},
    {"trans", (DL_FUNC) &trans, 8},
    {"_RxODE_rxInv", (DL_FUNC) &_RxODE_rxInv, 1},
    {"_RxODE_RxODE_finalize_focei_omega",(DL_FUNC) &_RxODE_RxODE_finalize_focei_omega, 1},
    {"_RxODE_RxODE_finalize_log_det_OMGAinv_5",(DL_FUNC) &_RxODE_RxODE_finalize_log_det_OMGAinv_5, 1},
    {"_RxODE_rxCoutEcho", (DL_FUNC) &_RxODE_rxCoutEcho, 1},
    {"_RxODE_W_Cpp", (DL_FUNC) &_RxODE_W_Cpp,2},
    {"_rxKahanSum", (DL_FUNC) &_rxKahanSum,1},
    {"_rxNeumaierSum", (DL_FUNC) &_rxNeumaierSum,1},
    {"_rxPythonSum", (DL_FUNC) &_rxPythonSum, 1},
    {"_rxPairwiseSum", (DL_FUNC) &_rxPairwiseSum, 1},
    {"_rxSum", (DL_FUNC) &_rxSum, 1},
    {"_rxProd", (DL_FUNC) &_rxProd, 1},
    {"_rxSetSum",(DL_FUNC) &_rxSetSum, 1},
    {"_rxSetProd",(DL_FUNC) &_rxSetProd, 1},
    {"_RxODE_removableDrive", (DL_FUNC) &_RxODE_removableDrive, 1},
    {NULL, NULL, 0}
  };

  // C callables needed in FOCEi
  R_RegisterCCallable("RxODE","nEq",                 (DL_FUNC) nEq);
  R_RegisterCCallable("RxODE","nLhs",                (DL_FUNC) nLhs);
  R_RegisterCCallable("RxODE","rxLhs",               (DL_FUNC) rxLhs);
  R_RegisterCCallable("RxODE","nAllTimes",           (DL_FUNC) nAllTimes);
  R_RegisterCCallable("RxODE","rxEvid",              (DL_FUNC) rxEvid);
  R_RegisterCCallable("RxODE","rxCalcLhs",           (DL_FUNC) rxCalcLhs);
  R_RegisterCCallable("RxODE","nObs",                (DL_FUNC) nObs);
  R_RegisterCCallable("RxODE","RxODE_ode_solve_env", (DL_FUNC) RxODE_ode_solve_env);
  R_RegisterCCallable("RxODE","RxODE_ode_free",      (DL_FUNC) RxODE_ode_free);
  R_RegisterCCallable("RxODE","RxODE_safe_zero",     (DL_FUNC) RxODE_safe_zero);
  R_RegisterCCallable("RxODE","RxODE_safe_log",      (DL_FUNC) RxODE_safe_log);
  R_RegisterCCallable("RxODE","RxODE_sign_exp",      (DL_FUNC) RxODE_sign_exp);
  R_RegisterCCallable("RxODE","RxODE_abs_log",       (DL_FUNC) RxODE_abs_log);

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
  R_RegisterCCallable("RxODE","RxODE_sign_exp",         (DL_FUNC) RxODE_sign_exp);
  R_RegisterCCallable("RxODE","RxODE_abs_log",          (DL_FUNC) RxODE_abs_log);
  R_RegisterCCallable("RxODE","RxODE_abs_log1p",        (DL_FUNC) RxODE_abs_log1p);
  R_RegisterCCallable("RxODE","RxODE_solveLinB",        (DL_FUNC) RxODE_solveLinB);

  R_RegisterCCallable("RxODE","RxODE_sum",              (DL_FUNC) RxODE_sum);
  R_RegisterCCallable("RxODE","RxODE_prod",             (DL_FUNC) RxODE_prod);

  R_RegisterCCallable("RxODE","RxODE_pow",              (DL_FUNC) RxODE_pow);
  R_RegisterCCallable("RxODE","RxODE_pow_di",           (DL_FUNC) RxODE_pow_di);


  static const R_CMethodDef cMethods[] = {
    {"RxODE_InfusionRate", (DL_FUNC) &RxODE_InfusionRate, 1, RxODE_one_int_t},
    {"RxODE_par_ptr", (DL_FUNC) &RxODE_par_ptr, 1, RxODE_one_int_t},
    {"RxODE_jac_counter_val", (DL_FUNC) &RxODE_jac_counter_val, 0},
    {"RxODE_dadt_counter_val",(DL_FUNC) &RxODE_dadt_counter_val, 0},
    {"RxODE_jac_counter_inc", (DL_FUNC) &RxODE_jac_counter_inc, 0},
    {"RxODE_dadt_counter_inc",(DL_FUNC) &RxODE_dadt_counter_inc, 0},
    {"RxODE_podo",(DL_FUNC) &RxODE_podo, 0},
    {"RxODE_tlast",(DL_FUNC) &RxODE_tlast, 0},
    {"RxODE_transit4",(DL_FUNC) &RxODE_transit4, 4, RxODE_transit4_t},
    {"RxODE_transit3", (DL_FUNC) &RxODE_transit3, 4, RxODE_transit3_t},
    {"RxODE_factorial", (DL_FUNC) &RxODE_factorial, 1, RxODE_one_dbl_t},
    {"RxODE_safe_log", (DL_FUNC) &RxODE_safe_log, 1, RxODE_one_dbl_t},
    {"RxODE_safe_zero", (DL_FUNC) &RxODE_safe_zero, 1, RxODE_one_dbl_t},
    {"RxODE_as_zero", (DL_FUNC) &RxODE_as_zero, 1, RxODE_one_dbl_t},
    {"RxODE_sign_exp", (DL_FUNC) &RxODE_sign_exp, 2, RxODE_sign_exp_t},
    {"RxODE_abs_log", (DL_FUNC) &RxODE_abs_log, 1, RxODE_one_dbl_t},
    {"RxODE_abs_log1p", (DL_FUNC) &RxODE_abs_log1p, 1, RxODE_one_dbl_t},
    {"RxODE_solveLinB", (DL_FUNC) &RxODE_solveLinB, 12, RxODE_solveLinB_t},
    {"RxODE_sum", (DL_FUNC) &RxODE_sum, 2, RxODE_Sum_t},
    {"RxODE_prod", (DL_FUNC) &RxODE_prod, 2, RxODE_Sum_t},
    {NULL, NULL, 0, NULL}
  };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);

}

