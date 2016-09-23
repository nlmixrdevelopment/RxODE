#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dop853.h"
#define max(a, b) ((a) > (b) ? (a) : (b))
#ifdef __STANDALONE__
#define Rprintf printf
#define R_alloc calloc
#else
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#endif


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
double *all_times;
FILE *fp;


void __DYDT__(unsigned int neq, double t, double *A, double *DADT);
void __CALC_LHS__(double t, double *A, double *lhs);
void __CALC_JAC__(unsigned int neq, double t, double *A, double *JAC, unsigned int __NROWPD__);

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
    }
  }
}

//--------------------------------------------------------------------------
void dydt_lsoda_dum(int *neq, double *t, double *A, double *DADT)
{
  __DYDT__(*neq, *t, A, DADT);
}
void jdum_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){
  // Update all covariate parameters
  __CALC_JAC__(*neq, *t, A, JAC, *nrowpd);
}
void call_lsoda(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc)
{
  int ixds=0, i, j;
  double xout, xp=x[0], yp[99];
  int itol = 1;
  double  rtol = RTOL, atol = ATOL;
  // Set jt to 1 if full is specified.
  int itask = 1, istate = 1, iopt = 0, lrw=22+neq*max(16, neq+9), liw=20+neq, jt = __JT__;
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
#ifdef __DEBUG__
  Rprintf("JT: %d\n",jt);
#endif
  rwork = (double*)R_alloc(lrw, sizeof(double));
  iwork = (int*)R_alloc(liw, sizeof(int));

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
#ifdef __DEBUG__
      fprintf(fp, "i=%d xp=%f xout=%f\n", i, xp, xout);
#endif

      if(xout>xp)
	{
	  F77_CALL(dlsoda)(dydt_lsoda_dum, &neq, yp, &xp, &xout, &itol, &rtol, &atol, &itask,
			   &istate, &iopt, rwork, &lrw, iwork, &liw, &jdum_lsoda, &jt);

	  if (istate<0)
	    {
	      Rprintf("IDID=%d, %s\n", istate, err_msg[-istate-1]);
#ifdef __STANDALONE__
	      exit(1);
#else
	      *rc = istate;
	      return;  // exit(1);   // dj: should not abort R
#endif
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

#ifdef __DEBUG__
      Rprintf("ISTATE=%d, ", istate);
      fprintf(fp, "ISTATE=%d, ", istate);
      for(j=0; j<neq; j++)
	{
	  Rprintf("%f ", yp[j]);
	  fprintf(fp, "%f ", yp[j]);
	}
      Rprintf("\n");
      fprintf(fp, "\n");
#endif
    }

#ifdef __STANDALONE__
  free(rwork);
  free(iwork);
#endif

}


void dydt_dvode_dum(int *neq, double *t, double *A, double *DADT, double *RPAR, int *IPAR)
{
  __DYDT__(*neq, *t, A, DADT);
}

void jdum_dvode(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd, double *RPAR, int *IPAR) {
  __CALC_JAC__(*neq, *t, A, JAC, *nrowpd);
}

void call_dvode(int neq, double *x, int *evid, int nx, double *inits, double *dose, double *ret, int *rc)
{
	int ixds=0, i, j;
	//DE solver config vars
	double xout, xp=x[0], yp[99];
    int itol = 1;
    double  rtol = RTOL, atol = ATOL;
    int itask = 1, istate = 1, iopt = 0, mf=__MF__;
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
#ifdef __DEBUG__
        Rprintf("MF: %d\n",mf);
#endif
        rwork = (double*)R_alloc(lrw, sizeof(double));
	iwork = (int*)R_alloc(liw, sizeof(int));
#ifdef __STANDALONE__
	if (!(rwork && iwork))
	{
		Rprintf("failed to alloc memory\n");
		exit(1);  
	}
#endif

	//--- inits the system
	for(i=0; i<neq; i++) yp[i] = inits[i];

	for(i=0; i<nx; i++)
	{
		wh = evid[i];
		xout = x[i];
#ifdef __DEBUG__
		fprintf(fp, "i=%d xp=%f xout=%f\n", i, xp, xout);
#endif

		if(xout>xp)
		{
	        F77_CALL(dvode)(dydt_dvode_dum, &neq, yp, &xp, &xout, &itol, &rtol, &atol, &itask,
				&istate, &iopt, rwork, &lrw, iwork, &liw, &jdum_dvode, &mf, rpar, ipar);

			if (istate<0)
			{
				Rprintf("IDID=%d, %s\n", istate, err_msg[-istate-1]);
#ifdef __STANDALONE__
            exit(1);
#else
				*rc = istate;
				return;    //exit(1);  // dj: should not abort R
#endif
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

#ifdef __DEBUG__
		Rprintf("ISTATE=%d, ", istate);
		fprintf(fp, "ISTATE=%d, ", istate);
		for(j=0; j<neq; j++)
		{
			Rprintf("%f ", yp[j]);
			fprintf(fp, "%f ", yp[j]);
		}
		Rprintf("\n");
		fprintf(fp, "\n");
#endif
	}

#ifdef __STANDALONE__
	free(rwork);
	free(iwork);
#endif
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
#ifdef __DEBUG__
		fprintf(fp, "i=%d xp=%f xout=%f\n", i, xp, xout);
#endif

		if(xout>xp+DBL_EPSILON)
		{
			idid = dop853(
							  neq,      	/* dimension of the system <= UINT_MAX-1*/
							  __DYDT__,    	/* function computing the value of f(x,y) */
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
#ifdef __STANDALONE__
            exit(1);
#else
				*rc = idid;
				return;  //exit(1);  // dj: should not abort R
#endif
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

#ifdef __DEBUG__
		Rprintf("IDID=%d, ", idid);
		fprintf(fp, "IDID=%d, ", idid);
		for(j=0; j<neq; j++)
		{
			Rprintf("%f ", yp[j]);
			fprintf(fp, "%f ", yp[j]);
		}
		Rprintf("\n");
		fprintf(fp, "\n");
#endif
	}
}

//wrapper
void __ODE_SOLVER__(
	int *neq,
	double *theta,	//order:
	double *time,
	int *evid,
	int *ntime,
	double *inits,
	double *dose,
	double *ret,
	double *atol,
	double *rtol,
	int *mx,
	int *stiff,
	int *transit_abs,
	int *nlhs,
	double *lhs,
	// Covariance terms
	int *pcov,
	double *cov,
	int *n_cov,
	int *locf,
	// Solver options
	double *h0,
	double *hmin,
	double *hmax,
	int *mxordn,
	int *mxords,
	// Return code
	int *counts,
 	int *rc){
  
        int i;
	for (i=0; i< 99; i++) InfusionRate[i] = 0.0;
	ATOL = *atol;
	RTOL = *rtol;
	HMIN = *hmin;
	HMAX = *hmax;
	H0   = *h0;
	MXORDN = *mxordn;
	MXORDS = *mxords;
	do_transit_abs = *transit_abs;
	par_ptr = theta;
	all_times = time;
	n_all_times = *ntime;
	is_locf = *locf;
	mxstep = *mx;
	
	par_cov = pcov;
	cov_ptr = cov;
	ncov    = *n_cov;

	slvr_counter = 0;
	dadt_counter = 0;
	jac_counter  = 0;
        if (*neq) {
	  if (*stiff==0)
            call_dop(*neq, time, evid, *ntime, inits, dose, ret, rc);
          else
            call_lsoda(*neq, time, evid, *ntime, inits, dose, ret, rc);
	}
	if (*nlhs) for (i=0; i<*ntime; i++){
	    // Update covariate parameters
	    __CALC_LHS__(time[i], ret+i*(*neq), lhs+i*(*nlhs));
        }
	counts[0] = slvr_counter;
	counts[1] = dadt_counter;
	counts[2] = jac_counter;
	if (fp) fclose(fp);
}
