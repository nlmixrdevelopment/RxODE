#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>
#include "solve.h"
#include "dop853.h"
// Yay easy parallel support
// For Mac, see: http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/ (as far as I can tell)
// It may have arrived, though I'm not sure...
// According to http://dirk.eddelbuettel.com/papers/rcpp_parallel_talk_jan2015.pdf
// OpenMP is excellent for parallelizing existing loops where the iterations are independent;
// OpenMP is used by part of the R core, therefore support will come for all platforms at some time in the future.
// Since these are independent, we will just use Open MP.
#ifdef _OPENMP
#include <omp.h>
#endif
#define max(a, b) ((a) > (b) ? (a) : (b))

int rxUpdateResiduals_(SEXP md);

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
			  int cores,
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
rx_solve *getRxSolve(SEXP ptr);

extern rx_solve *getRxSolve_(SEXP ptr){
  if(!R_ExternalPtrAddr(ptr)){
    error("Cannot get the solving data.");
  }
  rx_solve *o;
  o  = R_ExternalPtrAddr(ptr);
  return o;
}

rx_solving_options *getRxOp(rx_solve *rx){
  if(!R_ExternalPtrAddr(rx->op)){
    error("Cannot get global ode solver options.");
  }
  return (R_ExternalPtrAddr(rx->op));
}

rx_solving_options_ind *getRxId(rx_solve *rx, unsigned int id){
  return &(rx->subjects[id]);
}

extern int nEqP (rx_solve *rx){
  rx_solving_options *op;
  op =getRxOp(rx);
  return op->neq;
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

extern unsigned int nDosesP(rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  if (ind->ndoses < 0){
    ind->ndoses=0;
    for (int i = 0; i < ind->n_all_times; i++){
      if (rxEvidP(i, rx, id)){
        ind->ndoses++;
        if (ind->ndoses >= ind->idosen){
          if (ind->idosen == 0){
            ind->idose = Calloc(32,int);
            ind->idosen = 32;
          } else {
            ind->idosen *= 2;
            /* Rprintf("Reallocating to %d\n", ind->idosen); */
            ind->idose = Realloc(ind->idose, ind->idosen, int);
          }
        }
        ind->idose[ind->ndoses-1] = i;
      }
    }
    return ind->ndoses;
  } else {
    return ind->ndoses;
  }
}

extern unsigned int nObsP (rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  return (unsigned int)(ind->n_all_times - nDosesP(rx, id));
}

extern unsigned int nLhsP(rx_solve *rx){
  rx_solving_options *op;
  op =getRxOp(rx);
  return (unsigned int)(op->nlhs);
}
extern double rxLhsP(int i, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  rx_solving_options *op;
  op =getRxOp(rx);
  if (i < op->nlhs){
    return(ind->lhs[i]);
  } else {
    error("Trying to access an equation that isn't calculated. lhs(%d)\n",i);
  }
  return 0;
}
extern void rxCalcLhsP(int i, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  rx_solving_options *op;
  op =getRxOp(rx);
  t_calc_lhs calc_lhs = op->calc_lhs;
  double *solve, *lhs;
  solve = ind->solve;
  lhs = ind->lhs;
  if (i < ind->n_all_times){
    calc_lhs((int)id, ind->all_times[i], solve+i*op->neq, lhs);
  } else {
    error("LHS cannot be calculated (%dth entry).",i);
  }
}

extern unsigned int nAllTimesP(rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  return (unsigned int)(ind->n_all_times);
}


void F77_NAME(dlsoda)(
                      void (*)(int *, double *, double *, double *),
                      int *, double *, double *, double *, int *, double *, double *,
                      int *, int *, int *, double *,int *,int *, int *,
                      void (*)(int *, double *, double *, int *, int *, double *, int *),
                      int *);

void par_lsoda(SEXP sd){
  rx_solve *rx;
  rx = getRxSolve(sd);
  int i, j, foundBad;
  double xout;
  double *yp;
  rx_solving_options *op;
  if(!R_ExternalPtrAddr(rx->op)){
    error("Cannot get global ode solver options.");
  }
  op = R_ExternalPtrAddr(rx->op);
  int neq[2];
  neq[0] = op->neq;
  neq[1] = 0;
  yp = Calloc(neq[0], double);
  int itol = 1;
  double  rtol = op->RTOL, atol = op->ATOL;
  // Set jt to 1 if full is specified.
  int itask = 1, istate = 1, iopt = 0, lrw=22+neq[0]*max(16, neq[0]+9), liw=20+neq[0], jt = op->global_jt;
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
  int global_debug = op->global_debug;
  if (global_debug)
    Rprintf("JT: %d\n",jt);
  rwork = (double*)Calloc(lrw+1, double);
  iwork = (int*)Calloc(liw+1, int);
  
  iopt = 1;

  rwork[4] = op->H0; // H0 -- determined by solver
  rwork[6] = op->HMIN; // Hmin -- 0
  
  iwork[4] = 0; // ixpr  -- No extra printing.
  iwork[5] = op->mxstep; // mxstep 
  iwork[6] = 0; // MXHNIL 
  iwork[7] = op->MXORDN; // MXORDN 
  iwork[8] = op->MXORDS;  // MXORDS

  t_dydt_lsoda_dum dydt = (t_dydt_lsoda_dum)(op->dydt_lsoda_dum);
  t_jdum_lsoda jac = (t_jdum_lsoda)(op->jdum_lsoda);
  t_update_inis uini = (t_update_inis)(op->update_inis);
  int nx;
  rx_solving_options_ind *ind;
  double *inits;
  int *evid;
  double *x;
  int *BadDose;
  double *InfusionRate;
  double *dose;
  double *ret;
  int *rc;
  int nsub = rx->nsub;
  int nsim = rx->nsim;
  int cores = op->cores;
  int updateR = 1;
  for (int csim = 0; csim < nsim; csim++){
    // This part CAN be parallelized.
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int csub = 0; csub < nsub; csub++){
      neq[1] = csub+nsub*nsim;
      ind = &(rx->subjects[neq[1]]);
      ind->ixds = 0;
      nx = ind->n_all_times;
      inits = ind->inits;
      evid = ind->evid;
      BadDose = ind->BadDose;
      InfusionRate = ind->InfusionRate;
      dose = ind->dose;
      ret = ind->solve;
      x = ind->all_times;
      rc= ind->rc;
      rwork[5] = ind->HMAX; // Hmax -- Infinite
      double xp = x[0];
      //--- inits the system
      uini(neq[1], inits); // Update initial conditions
      for(i=0; i<neq[0]; i++) yp[i] = inits[i];
      for(i=0; i<nx; i++) {
        wh = evid[i];
        xout = x[i];
        if (global_debug){
          Rprintf("i=%d xp=%f xout=%f\n", i, xp, xout);
        }
        if(xout-xp> DBL_EPSILON*max(fabs(xout),fabs(xp)))
          {
            F77_CALL(dlsoda)(dydt, neq, yp, &xp, &xout, &itol, &rtol, &atol, &itask,
                             &istate, &iopt, rwork, &lrw, iwork, &liw, jac, &jt);

            if (istate<0)
              {
                Rprintf("IDID=%d, %s\n", istate, err_msg[-istate-1]);
                *rc = istate;
		i = nx+42; // Get out of here!
              }
            ind->slvr_counter++;
            //dadt_counter = 0;
          }
        if (wh)
          {
            cmt = (wh%10000)/100 - 1;
            if (cmt >= neq[0]){
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
              if (wh>10000)
                {
                  InfusionRate[cmt] += dose[ind->ixds];
                }
              else
                {
                  if (op->do_transit_abs)
                    {
                      ind->podo = dose[ind->ixds];
                      ind->tlast = xout;
                    }
                  else yp[cmt] += dose[ind->ixds];     //dosing before obs
                }
            
              istate = 1;

              ind->ixds++;
              xp = xout;
            }
          }
        for(j=0; j<neq[0]; j++) ret[neq[0]*i+j] = yp[j];
        //Rprintf("wh=%d cmt=%d tm=%g rate=%g\n", wh, cmt, xp, InfusionRate[cmt]);

        if (global_debug){
          Rprintf("ISTATE=%d, ", istate);
          for(j=0; j<neq[0]; j++)
            {
              Rprintf("%f ", yp[j]);
            }
          Rprintf("\n");
        }
      }
      
    }
    if (rc[0]){
      Rprintf("Error solving using LSODA\n");
      Free(rwork);
      Free(iwork);
      Free(yp);
      return;
    }
    if (updateR)
      updateR=rxUpdateResiduals_(sd);
  }
  Free(rwork);
  Free(iwork);
  Free(yp);
}

//dummy solout fn
void solout(long int nr, double t_old, double t,
            double *y, int *nptr, int *irtrn);//{}

void par_dop(SEXP sd){
  rx_solve *rx;
  rx = getRxSolve(sd);
  int i, j, foundBad;
  double xout;
  double *yp;
  rx_solving_options *op;
  if(!R_ExternalPtrAddr(rx->op)){
    error("Cannot get global ode solver options.");
  }
  op = R_ExternalPtrAddr(rx->op);
  int neq[2];
  neq[0] = op->neq;
  neq[1] = 0;
  yp = Calloc(neq, double);
  
  //DE solver config vars
  double rtol=op->RTOL, atol=op->ATOL;
  int itol=0;           //0: rtol/atol scalars; 1: rtol/atol vectors
  int iout=0;           //iout=0: solout() NEVER called
  int idid=0;
  int wh, cmt;
  char *err_msg[]=
    {
      "input is not consistent",
      "larger nmax is needed",
      "step size becomes too small",
      "problem is probably stiff (interrupted)"
    };
  t_dydt dydt = (t_dydt)(op->dydt);
  t_update_inis uini = (t_update_inis)(op->update_inis);
  rx_solving_options_ind *ind;
  double *inits;
  int *evid;
  double *x;
  int *BadDose;
  double *InfusionRate;
  double *dose;
  double *ret;
  int *rc;
  int nsub = rx->nsub;
  int nsim = rx->nsim;
  int cores = op->cores;
  int global_debug = op->global_debug;
  int nx;
  for (int csim = 0; csim < nsim; csim++){
    // This part CAN be parallelized.
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int csub = 0; csub < nsub; csub++){
      neq[1] = csub+nsub*nsim;
      ind = &(rx->subjects[neq[1]]);
      ind->ixds = 0;
      nx = ind->n_all_times;
      inits = ind->inits;
      evid = ind->evid;
      BadDose = ind->BadDose;
      InfusionRate = ind->InfusionRate;
      dose = ind->dose;
      ret = ind->solve;
      x = ind->all_times;
      rc= ind->rc;
      double xp = x[0];
      //--- inits the system
      uini(csub+nsub*nsim, inits); // Update initial conditions
      
      //--- inits the system
      for(i=0; i<neq[0]; i++) yp[i] = inits[i];

      for(i=0; i<nx; i++)
	{
	  wh = evid[i];
	  xout = x[i];
	  if (global_debug){
	    Rprintf("i=%d xp=%f xout=%f\n", i, xp, xout);
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
		  Rprintf("IDID=%d, %s\n", idid, err_msg[-idid-1]);
		  *rc = idid;
		  i = nx+42; // Get out of here!
		}
	      xp = xRead();
	      ind->slvr_counter++;
	      //dadt_counter = 0;
	    }
	  if (wh)
	    {
	      cmt = (wh%10000)/100 - 1;
	      if (cmt >= neq[0]){
		foundBad = 0;
		for (j = 0; j <ind->nBadDose; j++){
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
		if (wh>10000)
		  {
		    InfusionRate[cmt] += dose[ind->ixds];
		  }
		else
		  {
		    if (op->do_transit_abs)
		      {
			ind->podo = dose[ind->ixds];
                        ind->tlast = xout;
		      }
		    else yp[cmt] += dose[ind->ixds];     //dosing before obs
		  }
	      }
	      ind->ixds++;
	      xp = xout;
	    }
	  for(j=0; j<neq[0]; j++) ret[neq[0]*i+j] = yp[j];
	  //Rprintf("wh=%d cmt=%d tm=%g rate=%g\n", wh, cmt, xp, InfusionRate[cmt]);

	  if (global_debug){
	    Rprintf("IDID=%d, ", idid);
	    for(j=0; j<neq[0]; j++)
	      {
		Rprintf("%f ", yp[j]);
	      }
	    Rprintf("\n");
	  }
	}
    }
    if (rc[0]){
      Rprintf("Error sovling using dop853\n");
      return;
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

extern void update_par_ptrP(double t, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  rx_solving_options *op;
  op =getRxOp(rx);
  // Update all covariate parameters
  int k;
  int *par_cov = ind->par_cov;
  double *par_ptr = ind->par_ptr;
  double *all_times = ind->all_times;
  double *cov_ptr = ind->cov_ptr;
  int ncov = ind->ncov;
  if (ind->do_par_cov){
    for (k = 0; k < ind->ncov; k++){
      if (par_cov[k]){
        // Use the same methodology as approxfun.
        // There is some rumor the C function may go away...
        ind->ylow = cov_ptr[ind->n_all_times*k];
        ind->yhigh = cov_ptr[ind->n_all_times*k+ind->n_all_times-1];
        par_ptr[par_cov[k]-1] = rx_approxP(t, all_times, cov_ptr+ind->n_all_times*k, ind->n_all_times, op, ind);
      }
      if (op->global_debug){
        Rprintf("par_ptr[%d] (cov %d/%d) = %f\n",par_cov[k]-1, k,ncov,cov_ptr[par_cov[k]-1]);
      }
    }
  }
}
