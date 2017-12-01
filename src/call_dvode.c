#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "dop853.h"
#define NCMT 100
#define max(a, b) ((a) > (b) ? (a) : (b))
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>
#include <PreciseSums.h>
#include "solve.h"

void getSolvingOptionsIndPtr(double *InfusionRate,
			     int *BadDose,
			     double HMAX, // Determined by diff
			     double *par_ptr,
			     double *inits,
			     double *dose,
			     double *solve,
			     double *lhs,
			     int *par_cov,
			     int *evid,
			     int do_par_cov,
			     int *rc,
			     double *cov_ptr,
			     int ncov,
			     int n_all_times,
			     double *all_times,
			     int id,
			     int sim,
			     rx_solving_options_ind *o){
  o->slvr_counter = 0;
  o->dadt_counter = 0;
  o->jac_counter = 0;
  o->InfusionRate = InfusionRate;
  o->BadDose = BadDose;
  o->nBadDose = 0;
  o->HMAX = HMAX; // Determined by diff
  o->tlast = 0.0;
  o->podo = 0.0;
  o->par_ptr = par_ptr;
  o->inits = inits;
  o->dose = dose;
  o->solve = solve;
  o->lhs = lhs;
  o->par_cov = par_cov;
  o->evid = evid;
  o->do_par_cov = do_par_cov;
  o->rc = rc;
  o->cov_ptr = cov_ptr;
  o->ncov = ncov;
  o->n_all_times = n_all_times;
  o->ixds = 0;
  o->ndoses = -1;
  o->all_times = all_times;
  /* o.idose = idose;  allocated at run-time*/
  o->idosen = 0;
  o->id = id;
  o->sim = sim;
  o->extraCmt = 0;
}


static void getSolvingOptionsPtrFree(SEXP ptr)
{
  if(!R_ExternalPtrAddr(ptr)) return;
  rx_solving_options *o;
  o  = R_ExternalPtrAddr(ptr);
  Free(o);
  R_ClearExternalPtr(ptr);
}


SEXP getSolvingOptionsPtr(double ATOL,          //absolute error
			  double RTOL,          //relative error
			  double H0,
			  double HMIN,
			  int global_jt,
			  int global_mf,
			  int global_debug,
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
			  SEXP dydt,
			  SEXP calc_jac,
			  SEXP calc_lhs,
			  SEXP update_inis,
			  SEXP dydt_lsoda_dum,
			  SEXP jdum_lsoda){
  // This really should not be called very often, so just allocate one for now.
  rx_solving_options *o;
  o = Calloc(1,rx_solving_options);
  o->ATOL = ATOL;          //absolute error
  o->RTOL = RTOL;          //relative error
  o->H0 = H0;
  o->HMIN = HMIN;
  o->global_jt = global_jt;
  o->global_mf = global_mf;
  o->global_debug = global_debug;
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
  o->dydt = (t_dydt)(R_ExternalPtrAddr(dydt));
  o->calc_jac = (t_calc_jac)(R_ExternalPtrAddr(calc_jac));
  o->calc_lhs = (t_calc_lhs)(R_ExternalPtrAddr(calc_lhs));
  o->update_inis = (t_update_inis)(R_ExternalPtrAddr(update_inis));
  o->dydt_lsoda_dum = (t_dydt_lsoda_dum)(R_ExternalPtrAddr(dydt_lsoda_dum));
  o->jdum_lsoda = (t_jdum_lsoda)(R_ExternalPtrAddr(jdum_lsoda));
  SEXP ret = PROTECT(R_MakeExternalPtr(o, install("rx_solving_options"), R_NilValue));
  R_RegisterCFinalizerEx(ret, getSolvingOptionsPtrFree, TRUE);
  UNPROTECT(1);
  return(ret);
}

static void rxSolveDataFree(SEXP ptr) {
  if(!R_ExternalPtrAddr(ptr)) return;
  rx_solve *o;
  o  = R_ExternalPtrAddr(ptr);
  rx_solving_options_ind *inds;
  int n = (o->nsub)*(o->nsim);
  int *idose;
  inds = o->subjects;
  // Free all the idoses.
  for (int i = 0; i < n; i++){
    idose = (&inds[i])->idose;
    Free(idose);
  }
  // Free individuals;
  Free(inds);
  // Now free global options
  SEXP op = o->op;
  getSolvingOptionsPtrFree(op);
  // Now free object
  Free(o);
  R_ClearExternalPtr(ptr);
}

SEXP rxSolveData(rx_solving_options_ind *subjects,
                 int nsub,
                 int nsim,
                 SEXP op){
  rx_solve *o;
  o = Calloc(1,rx_solve);
  o->subjects = subjects;
  o->nsub = nsub;
  o->nsim = nsim;
  o->op = op;
  SEXP ret = PROTECT(R_MakeExternalPtr(o, install("rx_solve"), R_NilValue));
  R_RegisterCFinalizerEx(ret, rxSolveDataFree, TRUE);
  UNPROTECT(1);
  return(ret);
}

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

/* Authors: Robert Gentleman and Ross Ihaka and The R Core Team */
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

// Allow pointers to be called directly
void call_lsoda0(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc,
		 void (*fun_dydt_lsoda_dum)(int *neq, double *t, double *A, double *DADT),
		 void (*fun_jdum_lsoda)(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd),
		 double *ropts, int *iopts)
{
  int i, j, foundBad;
  double xout, xp=x[0], yp[99];
  int itol = 1;
  double  rtol = ropts[0], atol = ropts[1];
  // Set jt to 1 if full is specified.
  int itask = 1, istate = 1, iopt = 0, lrw=22+neq*max(16, neq+9), liw=20+neq, jt = iopts[0];
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

  rwork[4] = ropts[2]; // H0 -- determined by solver
  rwork[5] = ropts[3]; // Hmax -- Infinite
  rwork[6] = ropts[4]; // Hmin -- 0
  
  iwork[4] = 0; // ixpr  -- No extra printing.
  iwork[5] = iopts[1]; // mxstep 
  iwork[6] = 0; // MXHNIL 
  iwork[7] = iopts[2]; // MXORDN 
  iwork[8] = iopts[3];  // MXORDS
  
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
	  F77_CALL(dlsoda)(fun_dydt_lsoda_dum, &neq, yp, &xp, &xout, &itol, &rtol, &atol, &itask,
			   &istate, &iopt, rwork, &lrw, iwork, &liw, fun_jdum_lsoda, &jt);

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

void call_lsoda(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc){
  void (*fun_dydt_lsoda_dum)(int *neq, double *t, double *A, double *DADT);
  void (*fun_jdum_lsoda)(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd);
  fun_dydt_lsoda_dum = dydt_lsoda_dum;
  fun_jdum_lsoda = jdum_lsoda;
  double ropts[5];
  ropts[0] = RTOL;
  ropts[1] = ATOL;
  ropts[2] = H0;
  ropts[3] = HMAX;
  ropts[4] = HMIN;
  int iopts[4];
  iopts[0] = global_jt;
  iopts[1] = mxstep;
  iopts[2] = MXORDN;
  iopts[3] = MXORDS;
  call_lsoda0(neq, x, evid, nx, inits, dose, ret, rc, fun_dydt_lsoda_dum, fun_jdum_lsoda, ropts, iopts);
}


//dummy solout fn
void solout(long int nr, double t_old, double t,
	    double *y, unsigned int n, int *irtrn){}
void call_dop0(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc,
	       void (*fun_dydt)(unsigned int neq, double t, double *A, double *DADT),
	       double *ropts, int *iopts)
{
  int i, j;
  //DE solver config vars
  double xout, xp=x[0], yp[99];
  double rtol=ropts[0], atol=ropts[1];
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
			fun_dydt,    	/* function computing the value of f(x,y) */
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

void call_dop(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc){
  void (*fun_dydt)(unsigned int neq, double t, double *A, double *DADT);
  fun_dydt=dydt;
  double ropts[5];
  ropts[0] = RTOL;
  ropts[1] = ATOL;
  ropts[2] = H0;
  ropts[3] = HMAX;
  ropts[4] = HMIN;
  int iopts[4];
  iopts[0] = global_jt;
  iopts[1] = mxstep;
  iopts[2] = MXORDN;
  iopts[3] = MXORDS;
  call_dop0(neq, x, evid, nx, inits, dose, ret, rc,fun_dydt,ropts, iopts);
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
			    int mxstepA){
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
  mxstep         = mxstepA;
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

SEXP RxODE_get_fn_pointers(void (*fun_dydt)(unsigned int, double, double *, double *),
			   void (*fun_calc_lhs)(double, double *, double *),
			   void (*fun_calc_jac)(unsigned int, double, double *, double *, unsigned int),
			   void (*fun_update_inis)(SEXP _ini_sexp),
			   void (*fun_dydt_lsoda_dum)(int *, double *, double *, double *),
                           void (*fun_jdum_lsoda)(int *, double *, double *,int *, int *, double *, int *),
			   int fun_jt,
                           int fun_mf,
                           int fun_debug){
  SEXP dydt, lhs, jac, inis, dydt_lsoda, jdum;
  int pro=0;
  SEXP lst      = PROTECT(allocVector(VECSXP, 9)); pro++;
  SEXP names    = PROTECT(allocVector(STRSXP, 9)); pro++;

  void (*dydtf)(unsigned int neq, double t, double *A, double *DADT);
  void (*calc_jac)(unsigned int neq, double t, double *A, double *JAC, unsigned int __NROWPD__);
  void (*calc_lhs)(double t, double *A, double *lhs);
  void (*update_inis)(SEXP _ini_sexp);
  void (*dydt_lsoda_dum)(int *neq, double *t, double *A, double *DADT);
  void (*jdum_lsoda)(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd);
  
  dydtf		 = fun_dydt;
  calc_jac	 = fun_calc_jac;
  calc_lhs	 = fun_calc_lhs;
  update_inis	 = fun_update_inis;
  dydt_lsoda_dum = fun_dydt_lsoda_dum;
  jdum_lsoda     = fun_jdum_lsoda;
  
  SET_STRING_ELT(names,0,mkChar("dydt"));
  dydt=R_MakeExternalPtr(dydtf, install("RxODE_dydt"), R_NilValue);
  PROTECT(dydt); pro++;
  SET_VECTOR_ELT(lst,  0, dydt);

  SET_STRING_ELT(names,1,mkChar("lhs"));
  lhs=R_MakeExternalPtr(calc_lhs, install("RxODE_lhs"), R_NilValue);
  PROTECT(lhs); pro++;
  SET_VECTOR_ELT(lst,  1, lhs);

  SET_STRING_ELT(names,2,mkChar("jac"));
  jac=R_MakeExternalPtr(calc_jac, install("RxODE_jac"), R_NilValue);
  PROTECT(jac); pro++;
  SET_VECTOR_ELT(lst,  2, jac);

  SET_STRING_ELT(names,3,mkChar("inis"));
  inis=R_MakeExternalPtr(update_inis, install("RxODE_inis"), R_NilValue);
  PROTECT(inis); pro++;
  SET_VECTOR_ELT(lst,  3, inis);

  SET_STRING_ELT(names,4,mkChar("dydt_lsoda"));
  dydt_lsoda=R_MakeExternalPtr(dydt_lsoda_dum, install("RxODE_dydt_lsoda"), R_NilValue);
  PROTECT(dydt_lsoda); pro++;
  SET_VECTOR_ELT(lst,  4, dydt_lsoda);
  
  SET_STRING_ELT(names,5,mkChar("jdum"));
  jdum=R_MakeExternalPtr(jdum_lsoda, install("RxODE_jdum"), R_NilValue);
  PROTECT(jdum); pro++;
  SET_VECTOR_ELT(lst,  5, jdum);
  
  SET_STRING_ELT(names,6,mkChar("jt"));
  SEXP jt = PROTECT(allocVector(INTSXP, 1)); pro++;
  INTEGER(jt)[0] = fun_jt;
  SET_VECTOR_ELT(lst,  6, jt);

  SET_STRING_ELT(names,7,mkChar("mf"));
  SEXP mf = PROTECT(allocVector(INTSXP, 1)); pro++;
  INTEGER(mf)[0] = fun_mf;
  SET_VECTOR_ELT(lst,  7, mf);

  SET_STRING_ELT(names,8,mkChar("debug"));
  SEXP debug = PROTECT(allocVector(INTSXP, 1)); pro++;
  INTEGER(debug)[0] = fun_debug;
  SET_VECTOR_ELT(lst,  8, debug);
  setAttrib(lst, R_NamesSymbol, names);

  UNPROTECT(pro);
  return(lst);
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
  // LOCF
  if (is_locf == 1){
    rx_aprox_M.f2 = 0.0; //= f=0 
    rx_aprox_M.f1 = 1.0; // = 1-f = 1;
    rx_aprox_M.kind = 0;
  } else if (is_locf == 2) {
    // NOCB
    rx_aprox_M.f2 = 1.0; //= f=1
    rx_aprox_M.f1 = 0.0;
    rx_aprox_M.kind = 0;
  } else if (is_locf == 3){
    rx_aprox_M.f2 = 0.5; //= f=0.5
    rx_aprox_M.f1 = 0.5;
    rx_aprox_M.kind = 0;
  } else {
    // Linear
    rx_aprox_M.f2 = 1.0; //= f=0
    rx_aprox_M.f1 = 0.0;
    rx_aprox_M.kind = 1;
  }
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
  int pro = 0;
  SEXP sexp_theta = PROTECT(findVar(installChar(mkChar("params")),sexp_rho));pro++;
  SEXP sexp_inits = PROTECT(findVar(installChar(mkChar("inits")),sexp_rho)); pro++;
  SEXP sexp_lhs   = PROTECT(findVar(installChar(mkChar("lhs_vars")),sexp_rho)); pro++;
  // Events
  SEXP sexp_time = PROTECT(findVar(installChar(mkChar("time")),sexp_rho)); pro++;
  SEXP sexp_evid = PROTECT(findVar(installChar(mkChar("evid")),sexp_rho)); pro++;
  SEXP sexp_dose = PROTECT(findVar(installChar(mkChar("amt")),sexp_rho)); pro++;
  // Covariates
  SEXP sexp_pcov = PROTECT(findVar(installChar(mkChar("pcov")),sexp_rho)); pro++;
  SEXP sexp_cov = PROTECT(findVar(installChar(mkChar("cov")),sexp_rho)); pro++;
  SEXP sexp_locf = PROTECT(findVar(installChar(mkChar("isLocf")),sexp_rho)); pro++;
  // Solver Options
  SEXP sexp_atol = PROTECT(findVar(installChar(mkChar("atol")),sexp_rho)); pro++;
  SEXP sexp_rtol = PROTECT(findVar(installChar(mkChar("rtol")),sexp_rho)); pro++;
  SEXP sexp_hmin = PROTECT(findVar(installChar(mkChar("hmin")),sexp_rho)); pro++;
  SEXP sexp_hmax = PROTECT(findVar(installChar(mkChar("hmax")),sexp_rho)); pro++;
  SEXP sexp_h0 = PROTECT(findVar(installChar(mkChar("hini")),sexp_rho)); pro++;
  SEXP sexp_mxordn = PROTECT(findVar(installChar(mkChar("maxordn")),sexp_rho)); pro++;
  SEXP sexp_mxords = PROTECT(findVar(installChar(mkChar("maxords")),sexp_rho)); pro++;
  SEXP sexp_mx = PROTECT(findVar(installChar(mkChar("maxsteps")),sexp_rho)); pro++;
  SEXP sexp_stiff = PROTECT(findVar(installChar(mkChar("stiff")),sexp_rho)); pro++;
  SEXP sexp_transit_abs = PROTECT(findVar(installChar(mkChar("transit_abs")),sexp_rho)); pro++;
  SEXP sexp_rc = PROTECT(findVar(installChar(mkChar("rc")),sexp_rho)); pro++;
  int *rce    = INTEGER(sexp_rc);

  par_ptr       = REAL(sexp_theta);
  inits         = REAL(sexp_inits);

  RxODE_ode_setup(sexp_inits, sexp_lhs, sexp_time, sexp_evid, sexp_dose, sexp_pcov, sexp_cov,
		  sexp_locf, sexp_atol, sexp_rtol, sexp_hmin, sexp_hmax, sexp_h0, sexp_mxordn,
		  sexp_mxords, sexp_mx, sexp_stiff, sexp_transit_abs);
  RxODE_ode_alloc();
  RxODE_ode_solver_c(neq, stiff, evid, inits, dose, solve, rc);
  // Send rc to environment
  rce[0] = rc[0];
  UNPROTECT(pro);
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
		       SEXP sexp_opts,
                       SEXP sexp_extra_args){
  // TODO: Absorption lag?
  // TODO: Annotation? -- in nlmixr
  // TODO: Units -- Can be in nlmixr
  // TODO: Bioavailiability?
  // Parameters
  int pro=0;
  par_ptr       = REAL(sexp_theta);
  inits         = REAL(sexp_inits);
  // Events
  RxODE_ode_setup(sexp_inits, sexp_lhs, sexp_time, sexp_evid, sexp_dose, sexp_pcov, sexp_cov,
		  sexp_locf, sexp_atol, sexp_rtol, sexp_hmin, sexp_hmax, sexp_h0, sexp_mxordn,
		  sexp_mxords, sexp_mx, sexp_stiff, sexp_transit_abs);
  RxODE_ode_alloc();
  
  
  int i = 0, j = 0;
  int *rmState = INTEGER(sexp_opts);
  int nPrnState =0;
  for (i = 0; i < length(sexp_opts)-2; i++){
    nPrnState+= (1-rmState[i]);
  }
  SEXP sexp_counter = PROTECT(allocVector(INTSXP,4));pro++;
  int    *counts    = INTEGER(sexp_counter);

  int matrix = INTEGER(sexp_opts)[length(sexp_opts)-1];
  
  for (i=0; i<neq; i++) InfusionRate[i] = 0.0;
  if (neq > 0){
    RxODE_ode_solver_c(neq, stiff, evid, inits, dose, solve, rc);
  }
  double *scale = REAL(sexp_scale);
  int add_cov = INTEGER(sexp_opts)[length(sexp_opts)-2];
  if (matrix){
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
    SEXP temp = PROTECT(getAttrib(sexp_inits, R_NamesSymbol)); pro++;
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

    SEXP env = PROTECT(eval(lang1(install("new.env")),R_GlobalEnv));pro++;

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
SEXP _RxODE_RxODE_finalize_focei_omega(SEXP);
SEXP _RxODE_RxODE_finalize_log_det_OMGAinv_5(SEXP);
SEXP _rxCholInv(SEXP dms, SEXP theta, SEXP tn);
SEXP _RxODE_rxSymInvCholEnvCalculate(SEXP, SEXP, SEXP);
SEXP _RxODE_rxInvWishartVar(SEXP, SEXP);
SEXP _RxODE_rxSymInvChol(SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxDataSetup(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxIs(SEXP,SEXP);
SEXP _RxODE_rxModelVars(SEXP);
SEXP _RxODE_rxState(SEXP, SEXP);
SEXP _RxODE_rxParams(SEXP);
SEXP _RxODE_rxDfdy(SEXP);
SEXP _RxODE_rxLhs(SEXP);
SEXP _RxODE_rxInits(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxUpdateResiduals(SEXP);
SEXP _RxODE_rxSetupIni(SEXP, SEXP);
SEXP _RxODE_rxDataParSetup(SEXP, SEXP, SEXP, SEXP, SEXP,
			   SEXP, SEXP, SEXP, SEXP, SEXP,
			   SEXP);
SEXP _RxODE_rxSolvingOptions(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _RxODE_rxSolvingData(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

double RxODE_solveLinB(double t, int linCmt, int diff1, int diff2, double A, double alpha, double B, double beta, double C, double gamma, double ka, double tlag);
static R_NativePrimitiveArgType RxODE_solveLinB_t[] = {
  //t,    linCmt,  diff1,  diff2,  A,       alpha,  B,       beta,     C,       gamma, double ka, double tlag)
  REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

SEXP _RxODE_rxToOmega(SEXP cholInv);

double RxODE_sum(double *input, int len){
  return PreciseSums_sum(input, len);
}

extern double RxODE_sumV(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = PreciseSums_sum(p, n);
  Free(p);
  return s;
}

double RxODE_prod(double *input, int len){
  return PreciseSums_prod(input, len);
}

extern double RxODE_prodV(int n, ...){
  va_list valist;
  va_start(valist, n);
  double *p = Calloc(n, double);
  for (unsigned int i = 0; i < n; i++){
    p[i] = va_arg(valist, double);
  }
  va_end(valist);
  double s = PreciseSums_prod(p, n);
  Free(p);
  return s;
}

static R_NativePrimitiveArgType RxODE_Sum_t[] = {
  REALSXP, INTSXP
};


void R_init_RxODE(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"RxODE_ode_solver", (DL_FUNC) &RxODE_ode_solver, 23},
    {"trans", (DL_FUNC) &trans, 8},
    {"_RxODE_rxInv", (DL_FUNC) &_RxODE_rxInv, 1},
    {"_RxODE_RxODE_finalize_focei_omega",(DL_FUNC) &_RxODE_RxODE_finalize_focei_omega, 1},
    {"_RxODE_RxODE_finalize_log_det_OMGAinv_5",(DL_FUNC) &_RxODE_RxODE_finalize_log_det_OMGAinv_5, 1},
    {"_RxODE_rxCoutEcho", (DL_FUNC) &_RxODE_rxCoutEcho, 1},
    {"_RxODE_removableDrive", (DL_FUNC) &_RxODE_removableDrive, 1},
    {"_rxCholInv", (DL_FUNC) &_rxCholInv, 3},
    {"_RxODE_rxToOmega", (DL_FUNC) &_RxODE_rxToOmega, 1},
    {"_RxODE_rxSymInvCholEnvCalculate", (DL_FUNC) &_RxODE_rxSymInvCholEnvCalculate, 3},
    {"_RxODE_rxInvWishartVar", (DL_FUNC) &_RxODE_rxInvWishartVar, 2},
    {"_RxODE_rxSymInvChol", (DL_FUNC) &_RxODE_rxSymInvChol, 4},
    {"_RxODE_rxDataSetup", (DL_FUNC) &_RxODE_rxDataSetup, 8},
    {"_RxODE_rxIs", (DL_FUNC) &_RxODE_rxIs, 2},
    {"_RxODE_rxModelVars", (DL_FUNC) &_RxODE_rxModelVars, 1},
    {"_RxODE_rxState", (DL_FUNC) &_RxODE_rxState, 2},
    {"_RxODE_rxParams", (DL_FUNC) &_RxODE_rxParams, 1},
    {"_RxODE_rxDfdy", (DL_FUNC) &_RxODE_rxDfdy, 1},
    {"_RxODE_rxLhs", (DL_FUNC) &_RxODE_rxLhs, 1},
    {"_RxODE_rxInits", (DL_FUNC) &_RxODE_rxInits, 6},
    {"_RxODE_rxUpdateResiduals", (DL_FUNC) &_RxODE_rxUpdateResiduals, 1},
    {"_RxODE_rxSetupIni", (DL_FUNC) &_RxODE_rxSetupIni, 2},
    {"_RxODE_rxDataParSetup", (DL_FUNC) &_RxODE_rxDataParSetup, 11},
    {"_RxODE_rxSolvingOptions",(DL_FUNC) &_RxODE_rxSolvingOptions, 11},
    {"_RxODE_rxSolvingData", (DL_FUNC) &_RxODE_rxSolvingData, 13},
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
  R_RegisterCCallable("RxODE","RxODE_get_fn_pointers", (DL_FUNC) RxODE_get_fn_pointers);
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
    {"RxODE_InfusionRate",	(DL_FUNC) &RxODE_InfusionRate, 1, RxODE_one_int_t},
    {"RxODE_par_ptr",		(DL_FUNC) &RxODE_par_ptr, 1, RxODE_one_int_t},
    {"RxODE_jac_counter_val",	(DL_FUNC) &RxODE_jac_counter_val, 0},
    {"RxODE_dadt_counter_val",	(DL_FUNC) &RxODE_dadt_counter_val, 0},
    {"RxODE_jac_counter_inc",	(DL_FUNC) &RxODE_jac_counter_inc, 0},
    {"RxODE_dadt_counter_inc",	(DL_FUNC) &RxODE_dadt_counter_inc, 0},
    {"RxODE_podo",		(DL_FUNC) &RxODE_podo, 0},
    {"RxODE_tlast",		(DL_FUNC) &RxODE_tlast, 0},
    {"RxODE_transit4",		(DL_FUNC) &RxODE_transit4, 4, RxODE_transit4_t},
    {"RxODE_transit3",		(DL_FUNC) &RxODE_transit3, 4, RxODE_transit3_t},
    {"RxODE_factorial",		(DL_FUNC) &RxODE_factorial, 1, RxODE_one_dbl_t},
    {"RxODE_safe_log",		(DL_FUNC) &RxODE_safe_log, 1, RxODE_one_dbl_t},
    {"RxODE_safe_zero",		(DL_FUNC) &RxODE_safe_zero, 1, RxODE_one_dbl_t},
    {"RxODE_as_zero",		(DL_FUNC) &RxODE_as_zero, 1, RxODE_one_dbl_t},
    {"RxODE_sign_exp",		(DL_FUNC) &RxODE_sign_exp, 2, RxODE_sign_exp_t},
    {"RxODE_abs_log",		(DL_FUNC) &RxODE_abs_log, 1, RxODE_one_dbl_t},
    {"RxODE_abs_log1p",		(DL_FUNC) &RxODE_abs_log1p, 1, RxODE_one_dbl_t},
    {"RxODE_solveLinB",		(DL_FUNC) &RxODE_solveLinB, 12, RxODE_solveLinB_t},
    {"RxODE_sum",		(DL_FUNC) &RxODE_sum, 2, RxODE_Sum_t},
    {"RxODE_prod",		(DL_FUNC) &RxODE_prod, 2, RxODE_Sum_t},
    {NULL, NULL, 0, NULL}
  };

  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);

}

