#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> //Rmath includes math.
#include <R_ext/Rdynload.h>
#include "solve.h"
#include "dop853.h"
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
#define max(a, b) ((a) > (b) ? (a) : (b))

extern SEXP RxODE_get_fn_pointers(void (*fun_dydt)(int *neq, double t, double *A, double *DADT),
                                  void (*fun_calc_lhs)(int, double t, double *A, double *lhs),
                                  void (*fun_calc_jac)(int neq, double t, double *A, double *JAC, unsigned int __NROWPD__),
                                  void (*fun_update_inis)(int, double *),
                                  void (*fun_dydt_lsoda_dum)(int *, double *, double *, double *),
                                  void (*fun_jdum_lsoda)(int *, double *, double *,int *, int *, double *, int *),
                                  void (*fun_set_solve)(rx_solve *),
                                  rx_solve *(*fun_get_solve)(),
                                  int fun_jt,
                                  int fun_mf,
                                  int fun_debug){
  SEXP dydt, lhs, jac, inis, dydt_lsoda, jdum, get_solveS, set_solveS;
  int pro=0;
  SEXP lst      = PROTECT(allocVector(VECSXP, 11)); pro++;
  SEXP names    = PROTECT(allocVector(STRSXP, 11)); pro++;

  void (*dydtf)(int *neq, double t, double *A, double *DADT);
  void (*calc_jac)(int neq, double t, double *A, double *JAC, unsigned int __NROWPD__);
  void (*calc_lhs)(int, double t, double *A, double *lhs);
  void (*update_inis)(int neq, double *);
  void (*dydt_lsoda_dum)(int *neq, double *t, double *A, double *DADT);
  void (*jdum_lsoda)(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd);

  void (*set_solve)(rx_solve *);
  rx_solve *(*get_solve)();
  
  dydtf          = fun_dydt;
  calc_jac       = fun_calc_jac;
  calc_lhs       = fun_calc_lhs;
  update_inis    = fun_update_inis;
  dydt_lsoda_dum = fun_dydt_lsoda_dum;
  jdum_lsoda     = fun_jdum_lsoda;
  get_solve      = fun_get_solve;
  set_solve      = fun_set_solve;
  
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

  SET_STRING_ELT(names,9,mkChar("get_solve"));
  get_solveS=R_MakeExternalPtr(get_solve, install("RxODE_get_solve"), R_NilValue);
  PROTECT(get_solveS); pro++;
  SET_VECTOR_ELT(lst,  9, get_solveS);

  SET_STRING_ELT(names,10,mkChar("set_solve"));
  set_solveS=R_MakeExternalPtr(set_solve, install("RxODE_set_solve"), R_NilValue);
  PROTECT(set_solveS); pro++;
  SET_VECTOR_ELT(lst,  10, set_solveS);
  setAttrib(lst, R_NamesSymbol, names);

  UNPROTECT(pro);
  return(lst);
}

int rxUpdateResiduals_(SEXP md);

void getSolvingOptionsIndPtr(double *InfusionRate,
                             int *BadDose,
                             double HMAX, // Determined by diff
                             double *par_ptr,
                             double *dose,
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
  o->jac_counter = 0;
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
  /* o.idose = idose;  allocated at run-time*/
  o->idosen = 0;
  o->id = id;
  o->sim = sim;
}

static void getSolvingOptionsPtrFree(SEXP ptr)
{
  if(!R_ExternalPtrAddr(ptr)) return;
  rx_solving_options *o;
  o  = (rx_solving_options*)R_ExternalPtrAddr(ptr);
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
			  int ncov,
			  int *par_cov,
                          int do_par_cov,
                          double *inits,
                          SEXP stateNames,
			  SEXP lhsNames,
			  SEXP paramNames,
                          SEXP dydt,
                          SEXP calc_jac,
                          SEXP calc_lhs,
                          SEXP update_inis,
                          SEXP dydt_lsoda_dum,
                          SEXP jdum_lsoda,
			  SEXP set_solve,
			  SEXP get_solve){
  // This really should not be called very often, so just allocate one for now.
  rx_solving_options *o;
  o = (rx_solving_options*)R_chk_calloc(1,sizeof(*o));
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
  o->par_cov = par_cov;
  o->do_par_cov = do_par_cov;
  o->stateNames = stateNames;
  o->lhsNames = lhsNames;
  o->paramNames = paramNames;
  o->inits = inits;
  o->extraCmt = 0;
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

extern void rxSolveDataFree(SEXP ptr) {
  if(!R_ExternalPtrAddr(ptr)) return;
  rx_solve *o;
  o  = (rx_solve*)R_ExternalPtrAddr(ptr);
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
		 int *stateIgnore,
		 int nobs,
                 int add_cov,
                 int matrix,
                 SEXP op){
  rx_solve *o;
  o = (rx_solve*)R_chk_calloc(1,sizeof(*o));
  o->subjects = subjects;
  o->nsub = nsub;
  o->nsim = nsim;
  o->op = op;
  o->stateIgnore = stateIgnore;
  o->nobs = nobs;
  o->add_cov = add_cov;
  o->matrix = matrix;
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

rx_solve *getRxSolve(SEXP ptr);
extern rx_solve *getRxSolve_(SEXP ptr){
  if(!R_ExternalPtrAddr(ptr)){
    error("Cannot get the solving data (getRxSolve_).");
  }
  rx_solve *o;
  o  = (rx_solve*)R_ExternalPtrAddr(ptr);
  return o;
}

rx_solving_options *getRxOp(rx_solve *rx){
  if(!R_ExternalPtrAddr(rx->op)){
    error("Cannot get global ode solver options.");
  }
  return (rx_solving_options*)(R_ExternalPtrAddr(rx->op));
}

rx_solving_options_ind *getRxId(rx_solve *rx, unsigned int id){
  return &(rx->subjects[id]);
}

extern void par_lsoda(SEXP sd){
  rx_solve *rx;
  rx = getRxSolve(sd);
  int i, j, foundBad;
  double xout;
  double *yp;
  rx_solving_options *op;
  if(!R_ExternalPtrAddr(rx->op)){
    error("Cannot get global ode solver options.");
  }
  op = (rx_solving_options*)R_ExternalPtrAddr(rx->op);
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
  inits = op->inits;

  for (int csim = 0; csim < nsim; csim++){
    // This part CAN be parallelized.
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int csub = 0; csub < nsub; csub++){
      neq[1] = csub+csim*nsub;
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
		// Bad Solve => NA
		for (i = 0; i < nx*neq[0]; i++) ret[i] = NA_REAL;
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
    /* if (rc[0]){ */
    /*   /\* Rprintf("Error solving using LSODA\n"); *\/ */
    /*   /\* Free(rwork); *\/ */
    /*   /\* Free(iwork); *\/ */
    /*   /\* Free(yp); *\/ */
    /*   /\* return; *\/ */
    /* } */
    if (updateR)
      updateR=rxUpdateResiduals_(sd);
  }
  Free(rwork);
  Free(iwork);
  Free(yp);
}

//dummy solout fn
void solout(long int nr, double t_old, double t, double *y, int *nptr, int *irtrn){}

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
  op = (rx_solving_options*)R_ExternalPtrAddr(rx->op);
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
      neq[1] = csub+csim*nsub;
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
		  // Bad Solve => NA
                  for (i = 0; i < nx*neq[0]; i++) ret[i] = NA_REAL;
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

rx_solve *_globalRx = NULL;

extern void rxode_assign_rx(rx_solve *rx){
  _globalRx=rx;
}

extern void update_par_ptrP(double t, rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(rx, id);
  rx_solving_options *op;
  op =getRxOp(rx);
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
      if (op->global_debug){
        Rprintf("par_ptr[%d] (cov %d/%d) = %f\n",par_cov[k]-1, k,ncov,cov_ptr[par_cov[k]-1]);
      }
    }
  }
}
extern void update_par_ptr(double t){
  update_par_ptrP(t, _globalRx, 0);
}

extern int nEqP (rx_solve *rx){
  rx_solving_options *op;
  op =getRxOp(rx);
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

extern unsigned int nDoses(){
  return nDosesP(_globalRx, 0);
}

extern unsigned int nObsP (rx_solve *rx, unsigned int id){
  rx_solving_options_ind *ind;
  ind = getRxId(_globalRx, id);
  return (unsigned int)(ind->n_all_times - nDosesP(rx, id));
}

extern unsigned int nObs (){
  return nObsP (_globalRx, 0);
}

extern unsigned int nLhsP(rx_solve *rx){
  rx_solving_options *op;
  op =getRxOp(rx);
  return (unsigned int)(op->nlhs);
}
extern unsigned int nLhs(){
  return nObsP (_globalRx, 0);
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

extern double rxLhs(int i, rx_solve *rx, unsigned int id){
  return rxLhsP(i, _globalRx, 0);
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
  rx_solving_options *op;
  op =getRxOp(rx);
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
  if (i < nDoses(rx, id)){
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

extern SEXP RxODE_par_df(SEXP sd){
  rx_solve *rx;
  rx = getRxSolve(sd);
  rx_solving_options *op;
  if(!R_ExternalPtrAddr(rx->op)){
    error("Cannot get global ode solver options.");
  }
  op = (rx_solving_options*)R_ExternalPtrAddr(rx->op);
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
  SEXP paramNames = op->paramNames;
  int n =length(paramNames);
  SEXP df = PROTECT(allocVector(VECSXP,n+md+sm)); pro++;
  for (i = 0; i < md+sm; i++){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(INTSXP, nsim*nsub))); pro++;
  }
  for (i = 0; i < n; i++){
    SET_VECTOR_ELT(df, i+md+sm, PROTECT(allocVector(REALSXP, nsim*nsub))); pro++;
  }
  int jj = 0, ii = 0, j = 0, nall = 0;
  double *par_ptr;
  rx_solving_options_ind *ind;
  int *dfi;
  double *dfp;
  int nobs = rx->nobs;
  int csub;
  for (csub = 0; csub < nsub; csub++){
    ind = &(rx->subjects[csub]);
    nall+=ind->n_all_times;
  }
  // Event table information.
  SEXP dfe = PROTECT(allocVector(VECSXP,md+3)); pro++; // Events
  SEXP dfd = PROTECT(allocVector(VECSXP,md+3)); pro++; // Dosing
  SEXP dfs = PROTECT(allocVector(VECSXP,md+3)); pro++; // Sampling
  SEXP dfn = PROTECT(allocVector(STRSXP,md+3)); pro++;
  if (md){
    SET_VECTOR_ELT(dfe, 0, PROTECT(allocVector(INTSXP, nall))); pro++;
    SET_VECTOR_ELT(dfd, 0, PROTECT(allocVector(INTSXP, nall-nobs))); pro++;
    SET_VECTOR_ELT(dfs, 0, PROTECT(allocVector(INTSXP, nobs))); pro++;
    SET_STRING_ELT(dfn, 0, mkChar("id"));
    i++;
  }
  SEXP isEt = PROTECT(allocVector(LGLSXP,nall)); pro++;
  
  SET_VECTOR_ELT(dfe, md, PROTECT(allocVector(REALSXP, nall))); pro++;
  SET_VECTOR_ELT(dfd, md, PROTECT(allocVector(REALSXP, nall-nobs))); pro++;
  SET_VECTOR_ELT(dfs, md, PROTECT(allocVector(REALSXP, nobs))); pro++;
  SET_STRING_ELT(dfn, md, mkChar("time"));

  SET_VECTOR_ELT(dfe, md+1, PROTECT(allocVector(INTSXP, nall))); pro++;
  SET_VECTOR_ELT(dfd, md+1, PROTECT(allocVector(INTSXP, nall-nobs))); pro++;
  SET_VECTOR_ELT(dfs, md+1, PROTECT(allocVector(INTSXP, nobs))); pro++;
  SET_STRING_ELT(dfn, md+1, mkChar("evid"));

  SET_VECTOR_ELT(dfe, md+2, PROTECT(allocVector(REALSXP, nall))); pro++;
  SET_VECTOR_ELT(dfd, md+2, PROTECT(allocVector(REALSXP, nall-nobs))); pro++;
  SET_VECTOR_ELT(dfs, md+2, PROTECT(allocVector(REALSXP, nobs))); pro++;
  SET_STRING_ELT(dfn, md+2, mkChar("amt"));
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
  setAttrib(dfs, R_RowNamesSymbol, dfrs);

  setAttrib(dfs, R_NamesSymbol, dfn);
  setAttrib(dfd, R_NamesSymbol, dfn);
  setAttrib(dfe, R_NamesSymbol, dfn);

  SEXP clse = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(clse, 0, mkChar("data.frame"));

  classgets(dfs, clse);
  classgets(dfd, clse);
  classgets(dfe, clse);
  double *times, *doses;
  int evid, iie = 0, iis = 0, iid = 0, curdose, ntimes;
  for (csub = 0; csub < nsub; csub++){
    ind = &(rx->subjects[csub]);
    ntimes = ind->n_all_times;
    times = ind->all_times;
    doses = ind->dose;
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
        jj++;
      }
      // id
      if (md){
        dfi = INTEGER(VECTOR_ELT(df, jj));
        dfi[ii] = csub+1;
        jj++;
      }
      for (i = 0; i < n; i++){
	dfp=REAL(VECTOR_ELT(df, jj));
	dfp[ii] = par_ptr[i];
	jj++;
      }
      ii++;
    }
  }
  SEXP sexp_rownames = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(sexp_rownames)[0] = NA_INTEGER;
  INTEGER(sexp_rownames)[1] = -nsub*nsim;
  setAttrib(df, R_RowNamesSymbol, sexp_rownames);
  SEXP sexp_colnames = PROTECT(allocVector(STRSXP,n+sm+md)); pro++;
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
  for (i = 0; i < n; i++){
    SET_STRING_ELT(sexp_colnames, jj, STRING_ELT(paramNames,i));
    jj++;
  }
  setAttrib(df, R_NamesSymbol, sexp_colnames);
  SEXP cls = PROTECT(allocVector(STRSXP, 1)); pro++;
  SET_STRING_ELT(cls, 0, mkChar("data.frame"));
  classgets(df, cls);
  SEXP ret = PROTECT(allocVector(VECSXP,5)); pro++;
  SET_VECTOR_ELT(ret, 0, df);
  SET_VECTOR_ELT(ret, 1, dfe);
  SET_VECTOR_ELT(ret, 2, dfs);
  SET_VECTOR_ELT(ret, 3, dfd);
  SET_VECTOR_ELT(ret, 4, isEt);
  UNPROTECT(pro);
  return ret;
}

extern SEXP RxODE_df(SEXP sd){
  rx_solve *rx;
  rx = getRxSolve(sd);
  rx_solving_options *op;
  if(!R_ExternalPtrAddr(rx->op)){
    error("Cannot get global ode solver options.");
  }
  op = (rx_solving_options*)R_ExternalPtrAddr(rx->op);
  int add_cov = rx->add_cov;
  int ncov = op->ncov;
  int nlhs = op->nlhs;
  int nobs = rx->nobs;
  int nsim = rx->nsim;
  int *rmState = rx->stateIgnore;
  int nPrnState =0;
  int i, j;
  int neq[2];
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
  int nidCols = md + sm;
  int pro = 0;
  SEXP df = PROTECT(allocVector(VECSXP,ncols+nidCols)); pro++;
  for (i = 0; i < nidCols; i++){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(INTSXP, nobs*nsim))); pro++;
  }
  for (i = md + sm; i < ncols + nidCols; i++){
    SET_VECTOR_ELT(df, i, PROTECT(allocVector(REALSXP, nobs*nsim))); pro++;
  }
  int csub = 0, evid;
  int nsub = rx->nsub;
  rx_solving_options_ind *ind;
  // Now create the data frame
  double *dfp;
  int *dfi;
  int ii=0, jj = 0, ntimes;
  double *solve;
  double *cov_ptr;
  
  for (int csim = 0; csim < nsim; csim++){
    for (csub = 0; csub < nsub; csub++){
      neq[1] = csub+csim*nsub;
      ind = &(rx->subjects[neq[1]]);
      ntimes = ind->n_all_times;
      solve =  ind->solve;
      cov_ptr = ind->cov_ptr;
      for (i = 0; i < ntimes; i++){
	jj  = 0 ;
	evid = rxEvidP(i,rx,neq[1]);
	if (evid != 0 && csub == 0){
        } else if (evid==0){
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
          // time
          dfp = REAL(VECTOR_ELT(df, jj));
          dfp[ii] = ind->all_times[i];
	  jj++;
	  
          // States
          if (nPrnState){
            for (j = 0; j < neq[0]; j++){
               if (!rmState[j]){
                 dfp = REAL(VECTOR_ELT(df, jj));
                 dfp[ii] = solve[j+i*neq[0]];//scale[j];
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
	      dfp[ii] = cov_ptr[j*ntimes+i];
	      jj++;
	    }
          }
          ii++;
        }
      }
    }
  }
  SEXP sexp_rownames = PROTECT(allocVector(INTSXP,2)); pro++;
  INTEGER(sexp_rownames)[0] = NA_INTEGER;
  INTEGER(sexp_rownames)[1] = -nobs*nsim;
  setAttrib(df, R_RowNamesSymbol, sexp_rownames);
  SEXP sexp_colnames = PROTECT(allocVector(STRSXP,ncols+nidCols)); pro++;
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
  SET_STRING_ELT(sexp_colnames, jj, mkChar("time"));
  jj++;

  // Put in state names
  SEXP stateNames = op->stateNames;
  if (nPrnState){
    for (j = 0; j < neq[0]; j++){
      if (!rmState[j]){
	SET_STRING_ELT(sexp_colnames, jj, STRING_ELT(stateNames,j));
        jj++;
      }
    }
  }
  // Put in LHS names
  SEXP lhsNames = op->lhsNames;
  for (i = 0; i < nlhs; i++){
    SET_STRING_ELT(sexp_colnames, jj, STRING_ELT(lhsNames,i));
    jj++;
  }
  // Put in Cov names
  SEXP paramNames = op->paramNames;
  int *par_cov = op->par_cov;
  for (i = 0; i < ncov*add_cov; i++){
    SET_STRING_ELT(sexp_colnames,jj, STRING_ELT(paramNames, par_cov[i]-1));
    jj++;
  }
  setAttrib(df, R_NamesSymbol, sexp_colnames);
  UNPROTECT(pro);
  return df;
}

