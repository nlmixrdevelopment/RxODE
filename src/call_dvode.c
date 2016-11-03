#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dop853.h"
#define max(a, b) ((a) > (b) ? (a) : (b))
#define RLEN (neq+nlhs+ncov+1)
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>


void F77_NAME(dlsoda)(
     void (*)(int *, double *, double *, double *),
     int *, double *, double *, double *, int *, double *, double *,
     int *, int *, int *, double *,int *,int *, int *,
     void (*)(int *, double *, double *, int *, int *, double *, int *),
     int *);

void F77_NAME(dvode)(
     void (*)(int *, double *, double *, double *, double *,int *),
     int *, double *, double *, double *, int *, double *, double *,
     int *, int *, int *, double *,int *,int *, int *,
     void (*)(int *, double *, double *, int *, int *, double *, int *, double *, int *),
     int *, double *, int *);


long slvr_counter, dadt_counter, jac_counter;
double InfusionRate[99];
double ATOL;		//absolute error
double RTOL;		//relative error
double HMAX;
double H0;
double HMIN;

int    do_transit_abs=0;
double tlast=0;
double podo=0;
double *par_ptr;
int    *par_cov;
double *cov_ptr;
int    ncov;
int    is_locf;
int    n_all_times;
int    mxstep;
int    MXORDN;
int    MXORDS;
int    global_jt, global_mf, global_debug;
double *all_times;
FILE *fp;

/* void __DYDT__(unsigned int neq, double t, double *A, double *DADT); */
/* void __CALC_LHS__(double t, double *A, double *lhs); */
/* void __CALC_JAC__(unsigned int neq, double t, double *A, double *JAC, unsigned int __NROWPD__); */

void (*dydt)(unsigned int neq, double t, double *A, double *DADT);
void (*calc_jac)(unsigned int neq, double t, double *A, double *JAC, unsigned int __NROWPD__);
void (*calc_lhs)(double t, double *A, double *lhs);


void update_par_ptr(double t){
  // Update all covariate parameters
  int i;
  int j = -1;
  int k = -1;
  for (i = 0; i < ncov; i++){
    if (par_cov[i]){
      if (j == -1){
	if (t <= all_times[0]){
	  j = 0;
	} else if (t >= all_times[n_all_times-1]){
	  j = n_all_times-1;
	} else {
	  for (j = 0; j < n_all_times; j++){
	    if (all_times[j] > t){
              k = j-1;
              break;
	    } else if (all_times[j] == t){
	      break;
	    }
	  }
	}
      }
      if (k == -1){
	par_ptr[par_cov[i]-1] = cov_ptr[i*n_all_times+j];
      } else {
	if (!is_locf){
	  // Linear Interpolation
	  par_ptr[par_cov[i]-1] = (cov_ptr[i*n_all_times + k] - cov_ptr[i*n_all_times + j])/(all_times[k] - all_times[j])*(t-all_times[j])+cov_ptr[i*n_all_times+j];
	} else {
	  // LOCF
          par_ptr[par_cov[i]-1] = cov_ptr[i*n_all_times + k];
	}	
	// Spline?  I don't think it has much return on investment
      }
      if (global_debug){
	Rprintf("par_ptr[%d] (cov %d/%d) = %f\n",par_cov[i]-1, i,ncov,cov_ptr[par_cov[i]-1]);
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
  int ixds=0, i, j;
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
      wh = evid[i];
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


void dydt_dvode_dum(int *neq, double *t, double *A, double *DADT, double *RPAR, int *IPAR)
{
  dydt(*neq, *t, A, DADT);
}

void jdum_dvode(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd, double *RPAR, int *IPAR) {
  calc_jac(*neq, *t, A, JAC, *nrowpd);
}

void call_dvode(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc)
{
	int ixds=0, i, j;
	//DE solver config vars
	double xout, xp=x[0], yp[99];
    int itol = 1;
    double  rtol = RTOL, atol = ATOL;
    int itask = 1, istate = 1, iopt = 0, mf=global_mf;
    int lrw = 22+9*neq+2*neq*neq, liw = 30+neq;
    double *rwork;
    int *iwork;
    double *rpar=x;
    int *ipar=&neq;
	int wh, cmt;


	char *err_msg[]=
		{
			"excess work done on this call",
			"excess accuracy requested",
			"illegal input detected",
			"repeated error test failures",
			"repeated convergence failures",
			"error weight became zero during problem"
		};
	if (global_debug)
	  Rprintf("MF: %d\n",mf);
        rwork = (double*)R_alloc(lrw+1, sizeof(double));
	iwork = (int*)R_alloc(liw+1, sizeof(int));
	
	//--- inits the system
	for(i=0; i<neq; i++) yp[i] = inits[i];

	for(i=0; i<nx; i++)
	{
		wh = evid[i];
		xout = x[i];
		if (global_debug){
		  Rprintf("i=%d xp=%f xout=%f\n", i, xp, xout);
		  fprintf(fp, "i=%d xp=%f xout=%f\n", i, xp, xout);
		}
		
		if(xout-xp> DBL_EPSILON*max(fabs(xout), fabs(xp)))
		{
		  F77_CALL(dvode)(dydt_dvode_dum, &neq, yp, &xp, &xout, &itol, &rtol, &atol, &itask,
				&istate, &iopt, rwork, &lrw, iwork, &liw, &jdum_dvode, &mf, rpar, ipar);

			if (istate<0)
			{
			  Rprintf("IDID=%d, %s\n", istate, err_msg[-istate-1]);
			  error("Error solving using DVODE.");
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
	int ixds=0, i, j;
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
		wh = evid[i];
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
		       SEXP sexp_extra_args,
		       void (*fun_dydt)(unsigned int, double, double *, double *),
		       void (*fun_calc_lhs)(double, double *, double *),
		       void (*fun_calc_jac)(unsigned int, double, double *, double *, unsigned int),
		       int fun_jt,
		       int fun_mf,
		       int fun_debug){
  // TODO: Absorption lag?
  // TODO: Annotation?
  // TODO: Units
  // TODO: Bioavailiability
  // Assign functions pointers
  dydt     = fun_dydt;
  calc_jac = fun_calc_jac;
  calc_lhs = fun_calc_lhs;
  // Assign solver options
  global_jt	= fun_jt;
  global_mf	= fun_mf;
  global_debug	= fun_debug;
  // Parameters
  par_ptr       = REAL(sexp_theta);
  double *inits = REAL(sexp_inits);
  int nlhs      = length(sexp_lhs);
  // Events
  all_times     = REAL(sexp_time); 
  n_all_times   = length(sexp_time);
  int    *evid  = INTEGER(sexp_evid);
  double *dose  = REAL(sexp_dose);
  int    neq    = length(sexp_inits);
  // Covariates
  par_cov       = INTEGER(sexp_pcov);
  cov_ptr       = REAL(sexp_cov);
  ncov          = length(sexp_pcov);
  is_locf       = INTEGER(sexp_locf)[0];
  // Solver options
  ATOL           = REAL(sexp_atol)[0];
  RTOL           = REAL(sexp_rtol)[0];
  HMIN           = REAL(sexp_hmin)[0];
  HMAX           = REAL(sexp_hmax)[0];
  H0             = REAL(sexp_h0)[0];
  MXORDN         = INTEGER(sexp_mxordn)[0];
  MXORDS         = INTEGER(sexp_mxords)[0];
  mxstep         = INTEGER(sexp_mx)[0];
  do_transit_abs = INTEGER(sexp_transit_abs)[0];
  int stiff      = INTEGER(sexp_stiff)[0];
  slvr_counter   = 0;
  dadt_counter   = 0;
  jac_counter    = 0;
  int i = 0, j = 0;
  SEXP sexp_ret     = PROTECT(allocMatrix(REALSXP, n_all_times, ncov+1+neq+nlhs));
  SEXP sexp_counter = PROTECT(allocVector(INTSXP,4));
  int    *counts    = INTEGER(sexp_counter);

  double *solve, *lhs;
  double *ret   = REAL(sexp_ret);
  int *rc;
  
  rc = (int *) R_alloc(1,sizeof(int));
  rc[0] = 0;
  
  solve         = (double *) R_alloc(neq*n_all_times+1,sizeof(double));
  lhs           = (double *) R_alloc(nlhs,sizeof(double));

  for (i=0; i< 99; i++) InfusionRate[i] = 0.0;

  if (neq) {
    if (stiff==0){
      call_dop(neq, all_times, evid, n_all_times, inits, dose, solve, rc);
    } else{
      call_lsoda(neq, all_times, evid, n_all_times, inits, dose, solve, rc);
    }
  }
  // Now create the matrix.
  for (i = 0; i < n_all_times; i++){
    // Time
    ret[i] = all_times[i];
    // State
    if (neq){
      for (j = 0; j < neq; j++){
      	ret[n_all_times*(j+1)+i] = solve[j+i*neq];
      }
    }
    // LHS
    if (nlhs){
      calc_lhs(all_times[i], solve+i*neq, lhs);
      for (j = 0; j < nlhs; j++){
	ret[n_all_times*(j+1+neq)+i] = lhs[j];
      }
    }
    // Cov
    if (ncov > 0){
      for (j = 0; j < ncov; j++){
      	ret[n_all_times*(j+1+neq+nlhs)+i] = cov_ptr[j*n_all_times+i];
      }
    }
  }
  SEXP sexp_dimnames = PROTECT(allocVector(VECSXP,2));
  SET_VECTOR_ELT(sexp_dimnames, 0, R_NilValue);
  SEXP sexp_colnames = PROTECT(allocVector(STRSXP,1+neq+nlhs+ncov));
  SEXP sexp_solve    = PROTECT(allocVector(REALSXP,1+neq+nlhs+ncov));
  SET_STRING_ELT(sexp_colnames, 0, mkChar("time"));
  double *solver = REAL(sexp_solve);
  solver[0] = all_times[n_all_times-1];
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

void R_init_RxODE(DllInfo *info){
  //Function
  R_RegisterCCallable("RxODE","RxODE_ode_solver",       (DL_FUNC) RxODE_ode_solver);
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
}
