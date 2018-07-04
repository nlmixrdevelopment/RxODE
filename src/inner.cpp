#include <RcppArmadillo.h>
#define NETAs 20
#define NTHETAs 20
#define NSUBs 100
#define innerOde(id) ind_solve(rx, id, inner_dydt_liblsoda, inner_dydt_lsoda_dum, inner_jdum_lsoda, inner_dydt, inner_update_inis, inner_global_jt)
#define getOmegaInv() (as<arma::mat>(rxSymInvCholEnvCalculate(_rxInv, "omegaInv", R_NilValue)))
#define getOmegaDet() trace(as<arma::mat>(rxSymInvCholEnvCalculate(_rxInv, "log.det.OMGAinv.5", R_NilValue)))
#define getOmegaN() as<int>(rxSymInvCholEnvCalculate(_rxInv, "ntheta", R_NilValue))
#define getOmegaTheta() as<NumericVector>(rxSymInvCholEnvCalculate(_rxInv, "theta", R_NilValue));
#define setOmegaTheta(x) rxSymInvCholEnvCalculate(_rxInv, "theta", x)
#define tbs(x) powerD(x, op_focei.lambda, op_focei.yj)
#define tbsL(x) powerL(x, op_focei.lambda, op_focei.yj)
#define tbsDL(x) powerDL(x, op_focei.lambda, op_focei.yj)

using namespace Rcpp;
using namespace arma;
extern "C"{
#include "solve.h"
  typedef void (*S2_fp) (int *, int *, double *, double *, double *, int *, float *, double *);
  typedef void (*n1qn1_fp)(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
			   int mode[], int niter[], int nsim[], int imp[], int lp[], double zm[], int izs[], 
			   float rzs[], double dzs[]);

  void n1qn1_(S2_fp simul, int n[], double x[], double f[], double g[], double var[], double eps[],
              int mode[], int niter[], int nsim[], int imp[], int lp[], double zm[], int izs[], 
              float rzs[], double dzs[]) {
    static n1qn1_fp fun=NULL;
    if (fun == NULL) fun = (n1qn1_fp) R_GetCCallable("n1qn1","n1qn1_");
    fun(simul, n, x, f, g, var, eps, mode, niter, nsim, imp, lp, zm, izs, rzs, dzs);
  }

  void ind_solve(rx_solve *rx, unsigned int cid, t_dydt_liblsoda dydt_lls, 
		 t_dydt_lsoda_dum dydt_lsoda, t_jdum_lsoda jdum,
                 t_dydt c_dydt, t_update_inis u_inis, int jt);

}

RObject rxSymInvCholEnvCalculate(List obj, std::string what, Nullable<NumericVector> theta = R_NilValue);
bool rxIs(const RObject &obj, std::string cls);
double powerD(double x, double lambda, int yj);
double powerL(double x, double lambda, int yj);
double powerDL(double x, double lambda, int yj);

List _rxInv;

// These are focei inner options
typedef struct {
  // 
  double *geta;
  double *goldEta;
  double *gsaveEta;
  double *gthetaGrad;
  // n1qn1 specific vectors
  double *gZm;
  double *gG;
  double *gVar;

  // Integer of ETAs
  unsigned int etaTransN;
  unsigned int gEtaTransN;
  unsigned int gEtaGTransN;
  unsigned int gThetaGTransN;
  unsigned int gZmN;
  // Where likelihood is saved.
  
  double cEps; // Central difference EPS
  
  int *etaTrans;

  int neta;
  unsigned int ntheta;
  int npars;
  int thetan;
  int omegan;

  double *fullTheta;
  double *theta;
  double *thetaGrad;
  double *initPar;

  unsigned int thetaTransN;

  int *fixedTrans;
  int *thetaTrans;

  double scaleTo;
  double epsilon;

  unsigned int maxOuterIterations;
  int maxInnerIterations;

  int nsim;
  int nzm;

  int imp;

  int yj;
  double lambda;
  int estLambda;

  mat omegaInv;

} focei_options;

focei_options op_focei;

extern "C" void rxOptionsIniFocei(){
  op_focei.etaTransN = 0;
  op_focei.thetaTransN = 0;
  op_focei.gEtaGTransN=0;
}

void foceiThetaN(unsigned int n){
  if (op_focei.thetaTransN < n){
    unsigned int cur = op_focei.thetaTransN;
    while (cur < n){
      cur += NTHETAs;
    }
    Free(op_focei.thetaTrans);
    Free(op_focei.theta);
    Free(op_focei.fullTheta);
    Free(op_focei.initPar);
    Free(op_focei.fixedTrans);
    op_focei.fullTheta   = Calloc(cur, double);
    op_focei.thetaTrans  = Calloc(cur, int);
    op_focei.fixedTrans  = Calloc(cur, int);
    op_focei.theta       = Calloc(cur, double);
    op_focei.initPar     = Calloc(cur, double);
    op_focei.thetaTransN = cur;
  }
}

void foceiEtaN(unsigned int n){
  if (op_focei.etaTransN < n){
    unsigned int cur = op_focei.etaTransN;
    while (cur < n){
      cur += NETAs;
    }
    Free(op_focei.etaTrans);
    op_focei.etaTrans = Calloc(cur, int);
    op_focei.etaTransN=cur;
  }
}

void foceiGThetaN(unsigned int n){
  if (op_focei.gThetaGTransN < n){
    unsigned int cur = op_focei.gThetaGTransN;
    while (cur < n){
      cur += NSUBs*NTHETAs;
    }
    Free(op_focei.gthetaGrad);
    op_focei.gthetaGrad = Calloc(cur, double);
    op_focei.gThetaGTransN=cur;
  }
}

void foceiGEtaN(unsigned int n){
  if (op_focei.gEtaGTransN < n){
    unsigned int cur = op_focei.gEtaGTransN;
    while (cur < n){
      cur += NETAs*NSUBs;
    }
    Free(op_focei.geta);
    Free(op_focei.goldEta);
    Free(op_focei.gsaveEta);
    Free(op_focei.gG);
    Free(op_focei.gVar);
    op_focei.geta = Calloc(cur, double);
    op_focei.goldEta = Calloc(cur, double);
    op_focei.gsaveEta = Calloc(cur, double);
    op_focei.gG = Calloc(cur, double);
    op_focei.gVar = Calloc(cur, double);
    // Prefill to 0.1 or 10%
    std::fill_n(&op_focei.gVar[0], cur, 0.1);
    op_focei.gEtaGTransN = cur;
  }
}

void foceiGgZm(unsigned int n){
  if (op_focei.gZmN < n){
    unsigned int cur = op_focei.gZmN;
    while (cur < n){
      cur += NETAs*(NETAs+13)/2*NSUBs;
    }
    Free(op_focei.gZm);
    op_focei.gZm = Calloc(cur, double);
    op_focei.gZmN = cur;
  }
}

extern "C" void rxOptionsFreeFocei(){
  Free(op_focei.thetaTrans);
  Free(op_focei.theta);
  Free(op_focei.fullTheta);
  Free(op_focei.initPar);
  Free(op_focei.fixedTrans);
  op_focei.thetaTransN=0;
  Free(op_focei.etaTrans);
  op_focei.etaTransN=0;
  Free(op_focei.geta);
  Free(op_focei.goldEta);
  Free(op_focei.gsaveEta);
  op_focei.gEtaGTransN=0;
  Free(op_focei.gthetaGrad);
  op_focei.gThetaGTransN=0;
}

typedef struct {
  double lik[3]; // lik[0] = liklihood; For central difference: lik[1] = lower lik[2] = upper
  double *eta; // Eta includes the ID number for the patient
  //
  double *thetaGrad; // Theta gradient; Calculated on the individual level for S matrix calculation
  double thVal[2]; // thVal[0] = lower; thVal[2] = upper
  //
  // F and varaibility
  unsigned int nobs;
  unsigned int setup;
  
  double *saveEta; // Saved when lik[0] is saved.
  double *oldEta;
  
  // Likilihood gradient
  double llik;
  mat a;
  mat B;
  mat c;
  mat lp;// = mat(neta,1);

  double *g;
  double *var;

  int izs;
  float rzs; 
  double dzs;

  int mode; // 1 = dont use zm, 2 = use zm.
  double *zm;
  unsigned int uzm;
} focei_ind;

focei_ind *inds_focei = NULL;
int max_inds_focei = 0;

rx_solve *getRxSolve_();

focei_ind *rxFoceiEnsure(int mx){
  if (mx >= max_inds_focei){
    Free(inds_focei);
    inds_focei =Calloc(mx+1024, focei_ind);
    max_inds_focei = mx+1024;
  }
  return inds_focei;
}


t_dydt inner_dydt = NULL;

t_calc_jac inner_calc_jac = NULL;

t_calc_lhs inner_calc_lhs = NULL;

t_update_inis inner_update_inis = NULL;

t_dydt_lsoda_dum inner_dydt_lsoda_dum = NULL;

t_dydt_liblsoda inner_dydt_liblsoda = NULL;

t_jdum_lsoda inner_jdum_lsoda = NULL;

t_set_solve inner_set_solve = NULL;

t_get_solve inner_get_solve = NULL;

int inner_global_jt = 2;
int inner_global_mf = 22;  
int inner_global_debug = 0;

void rxUpdateInnerFuns(SEXP trans){
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
  inner_global_jt = 2;
  inner_global_mf = 22;  
  inner_global_debug = 0;
  if (strcmp(CHAR(STRING_ELT(trans, 1)),"fulluser") == 0){
    inner_global_jt = 1;
    inner_global_mf = 21;
  } else {
    inner_global_jt = 2;
    inner_global_mf = 22;
  }
  inner_calc_lhs =(t_calc_lhs) R_GetCCallable(lib, s_calc_lhs);
  inner_dydt =(t_dydt) R_GetCCallable(lib, s_dydt);
  inner_calc_jac =(t_calc_jac) R_GetCCallable(lib, s_calc_jac);
  inner_update_inis =(t_update_inis) R_GetCCallable(lib, s_inis);
  inner_dydt_lsoda_dum =(t_dydt_lsoda_dum) R_GetCCallable(lib, s_dydt_lsoda_dum);
  inner_jdum_lsoda =(t_jdum_lsoda) R_GetCCallable(lib, s_dydt_jdum_lsoda);
  inner_set_solve = (t_set_solve)R_GetCCallable(lib, s_ode_solver_solvedata);
  inner_get_solve = (t_get_solve)R_GetCCallable(lib, s_ode_solver_get_solvedata);
  inner_dydt_liblsoda = (t_dydt_liblsoda)R_GetCCallable(lib, s_dydt_liblsoda);
}

void rxClearInnerFuns(){
  inner_calc_lhs              = NULL;
  inner_dydt                  = NULL;
  inner_calc_jac              = NULL;
  inner_update_inis           = NULL;
  inner_dydt_lsoda_dum        = NULL;
  inner_jdum_lsoda            = NULL;
  inner_set_solve             = NULL;
  inner_get_solve             = NULL;
  inner_dydt_liblsoda         = NULL;
}

rx_solve* rx;

////////////////////////////////////////////////////////////////////////////////
// n1qn1 functions
uvec lowerTri(mat H, bool diag = false){
  unsigned int d = H.n_rows;
  mat o(d, d, fill::ones);
  if (!diag){
    return find(trimatl(o,-1));
  } else {
    return find(trimatl(o));
  }
}

void updateZm(focei_ind *indF){
  if (!indF->uzm){
    // Udate the curvature to Hessian to restart n1qn1
    int n = op_focei.neta;
    mat L = eye(n, n);
    mat D = mat(n, n, fill::zeros);
    mat H = mat(n, n);
    unsigned int l_n = n * (n + 1)/2;
    vec zmV(l_n);
    std::copy(&indF->zm[0], &indF->zm[0]+l_n, zmV.begin());
    H.elem(lowerTri(H,true)) = zmV;
    L.elem(lowerTri(H,false)) = H.elem(lowerTri(H,0));
    D.diag() = H.diag();
    H = L*D*L.t();
    // Hessian -> c.hess
    std::fill(&indF->zm[0], &indF->zm[0]+op_focei.nzm,0.0);
    vec hessV = H.elem(lowerTri(H,1));
    std::copy(hessV.begin(),hessV.end(),&indF->zm[0]);
    indF->uzm = 1;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Likelihood for inner functions

void updateTheta(double *theta){
  // Theta is the acutal theta
  unsigned int j, k;
  if (op_focei.scaleTo > 0){ // Scaling
    for (k = op_focei.npars; k--;){
      j=op_focei.fixedTrans[k];
      op_focei.fullTheta[j] = theta[k] * op_focei.initPar[j] / op_focei.scaleTo; //pars <- pars * inits.vec / con$scale.to
    }
  } else { // No scaling.
    for (k = op_focei.npars; k--;){
      j=op_focei.fixedTrans[k];
      op_focei.fullTheta[j] = theta[k]; //pars <- pars * inits.vec / con$scale.to
    }
  }
  // Update theta parameters in each individual
  rx = getRxSolve_();
  for (int id = rx->nsub; id--;){
    rx_solving_options_ind ind = rx->subjects[id];
    double *par_ptr = ind.par_ptr;
    for (j = op_focei.ntheta; j--;){
      par_ptr[op_focei.thetaTrans[j]] = op_focei.fullTheta[j];
    }
  }
  // Update setOmegaTheta
  NumericVector omegaTheta(op_focei.omegan);
  std::copy(&op_focei.fullTheta[0] + op_focei.thetan, 
	    &op_focei.fullTheta[0] + op_focei.thetan + op_focei.omegan, 
	    omegaTheta.begin());
  setOmegaTheta(omegaTheta);
  // Update Lambda, if needed.
  if (op_focei.estLambda){
    op_focei.lambda = op_focei.fullTheta[op_focei.npars-1];
  }
}

void updateTheta1(double newTheta0, int k){
  // Theta is the acutal theta
  unsigned int j = op_focei.fixedTrans[k];
  double newTheta;
  if (op_focei.scaleTo > 0){ // Scaling
    newTheta = newTheta0 * op_focei.initPar[j] / op_focei.scaleTo; //pars <- pars * inits.vec / con$scale.to
  } else { // No scaling.
    newTheta = newTheta0; 
  }
  // Update theta parameters in each individual
  rx = getRxSolve_();
  if (j < op_focei.ntheta){
    for (int id = rx->nsub; id--;){
      rx_solving_options_ind ind = rx->subjects[id];
      double *par_ptr = ind.par_ptr;
      par_ptr[op_focei.thetaTrans[j]] = newTheta;
    }
  } else {
    // Update setOmegaTheta
    NumericVector omegaTheta(op_focei.omegan);
    std::copy(&op_focei.fullTheta[0] + op_focei.thetan, 
              &op_focei.fullTheta[0] + op_focei.thetan + op_focei.omegan, 
              omegaTheta.begin());
    omegaTheta[j-op_focei.ntheta] = newTheta;
    setOmegaTheta(omegaTheta);
  }
  // Lambda update is not handled by finite difference, so not handled here.
}

double likInner(double *eta){
  // id = eta[-1]
  // eta = eta
  double *eprev = &eta[0]-1;
  unsigned int id = (unsigned int)(eprev[0]);
  rx_solving_options_ind ind = rx->subjects[id];
  rx_solving_options *op = rx->op;
  focei_ind fInd = (inds_focei[id]);
  double *par_ptr = ind.par_ptr;
  int i, j;
  bool recalc = false;
  if (!fInd.setup){
    recalc=true;
    fInd.nobs = ind.n_all_times - ind.ndoses;
    fInd.lp  = arma::mat(op_focei.neta, 1);
    fInd.a   = arma::mat(fInd.nobs, op_focei.neta);
    fInd.B   = arma::mat(fInd.nobs, 1);
    fInd.c   = arma::mat(fInd.nobs, op_focei.neta);
    fInd.setup = 1;
  } else {
    // Check to see if old ETA matches.
    for (j = op_focei.neta; j--;){
      if (fInd.oldEta[j] != eta[j]){
	recalc=true;
	break;
      }
    }
  }
  if (recalc){
    // Update eta.
    for (j = op_focei.neta; j--;){
      par_ptr[op_focei.etaTrans[j]] = eta[j];
    }
    // Solve ODE
    innerOde(id);
    // Calculate matricies
    unsigned int k = fInd.nobs - 1;
    fInd.lp.fill(0.0);
    fInd.llik=0.0;
    double f, err, r, fpm, rp;
    for (j = ind.n_all_times; j--;){
      if (!ind.evid[j]){
	// Observation; Calc LHS.
	inner_calc_lhs((int)id, ind.all_times[j], ind.solve + j * op->neq, ind.lhs);
        f = tbs(ind.lhs[0]);
	// fInd.f(k, 0) = ind.lhs[0];
	err = f - tbs(ind.dv[j]);
	// fInd.err(k, 0) = ind.lhs[0] - ind.dv[k]; // pred-dv
        r = ind.lhs[op_focei.neta + 1];
	// fInd.r(k, 0) = ind.lhs[op_focei.neta+1];
        fInd.B(k, 0) = 2.0/r;
	// fInd.B(k, 0) = 2.0/ind.lhs[op_focei.neta+1];
	// lhs 0 = F
	// lhs 1-eta = df/deta
	// FIXME faster initiliaitzation via copy or elim
	for (i = op_focei.neta; i--; ){
	  fpm = fInd.a(k, i) = ind.lhs[i + 1]; // Almquist uses different a (see eq #15)
	  rp  = ind.lhs[i + op_focei.neta + 2];
	  fInd.c(k, i) = rp/r;
	  // lp is eq 12 in Almquist 2015
	  // // .5*apply(eps*fp*B + .5*eps^2*B*c - c, 2, sum) - OMGAinv %*% ETA
	  fInd.lp(i, 0)  += 0.25 * err * err * fInd.B(k, 0) * fInd.c(k, i) - 0.5 * fInd.c(k, i) - 
	    0.5 * err * fpm * fInd.B(k, 0);
	}
	// Eq #10
	fInd.llik += -0.5 * err * err/r - 0.5*log(r);
        k--;
      }
    }
    // Now finalize lp
    mat etam = arma::mat(op_focei.neta, 1);
    std::copy(&eta[0], &eta[0] + op_focei.neta, etam.begin()); // fill in etam
    // Finalize eq. #12
    fInd.lp = -(fInd.lp - op_focei.omegaInv * etam);
    // Partially finalize #10
    fInd.llik = trace(-(fInd.llik - 0.5*(etam.t() * op_focei.omegaInv * etam)));
    std::copy(&eta[0], &eta[0] + op_focei.neta, &fInd.oldEta[0]);
  }
  return fInd.llik;
}

double *lpInner(double *eta){
  double *eprev = &eta[0]-1;
  unsigned int id = (unsigned int)(eprev[0]);
  focei_ind fInd = (inds_focei[id]);
  likInner(eta);
  std::copy(&fInd.g[0], &fInd.g[0] + op_focei.neta, fInd.lp.begin());
  return &fInd.oldEta[0];
}

mat HInner(double *eta){
  double *eprev = &eta[0]-1;
  unsigned int id = (unsigned int)(eprev[0]);
  focei_ind fInd = (inds_focei[id]);
  // Hessian 
  mat H(op_focei.neta, op_focei.neta, fill::zeros);
  int k, l;
  for (k = op_focei.neta; k--;){
    for (l = k; l--;){
      H(k, l) = trace(fInd.a.col(k) * fInd.B * fInd.a.col(l) - 
		      trace(fInd.c.col(l)*fInd.c.col(k)) - 
		      op_focei.omegaInv(k, l));
      H(k, l) = H(l, k);
    }
  }
  return H;
}

double LikInner2(double *eta, int likId){
  double *eprev = &eta[0]-1;
  unsigned int id = (unsigned int)(eprev[0]);
  focei_ind fInd = (inds_focei[id]);
  double lik = -likInner(eta) + getOmegaDet();
  rx = getRxSolve_();
  rx_solving_options_ind ind = rx->subjects[id];
  // Calclaute lik first to calculate components for Hessian
  mat H = -HInner(eta);
  H = chol(H);
  lik += trace(log(H.diag()));

  // Add likelihood contribution based on transform both sides.
  if (op_focei.lambda != 1.0){
    for (unsigned int j = ind.n_all_times; j--;){
      if (!ind.evid[j]){
	lik +=tbsL(ind.dv[j]);
      }
    }
  }
  fInd.lik[likId] = lik;
  if (likId == 0) std::copy(&fInd.eta[0], &fInd.eta[0] + op_focei.neta, &fInd.saveEta[0]);
  return lik;
}

// Scli-lab style cost function for inner
void innerCost(int *ind, int *n, double *x, double *f, double *g, int *ti, float *tr, double *td){
  if (*ind==2 || *ind==4) {
    // Function
    f[0] = likInner(x);
  }
  if (*ind==3 || *ind==4) {
    // Gradient
    g = lpInner(x);
  }
}

static inline void innerEval(int id){
  focei_ind fInd = (inds_focei[id]);
  // Use eta
  double *eta = &(fInd.eta[0])+1; // id#, eta
  likInner(eta);
  LikInner2(eta, 0);
}

static inline void innerOpt1(int id, int likId){
  focei_ind fInd = (inds_focei[id]);
  // Use eta
  double *eta = &(fInd.eta[0])+1; // id#, eta
  // Convert Zm to Hessian, if applicable.
  updateZm(&fInd);
  int lp = 6;
  n1qn1_(innerCost, &(op_focei.neta), eta, &(fInd.llik), fInd.g,  fInd.var,
	 &(op_focei.epsilon), &(fInd.mode), &(op_focei.maxInnerIterations),
	 &(op_focei.nsim), &(op_focei.imp), &lp, fInd.zm, &fInd.izs,
	 &fInd.rzs, &fInd.dzs);
  // Use saved Hessian on next opimization.
  fInd.mode=1;
  fInd.uzm =0;
  LikInner2(eta, likId);
}

void innerOpt(){
#ifdef _OPENMP
  int cores = rx->op->cores;
#else
  int cores = 1;
#endif
  rx = getRxSolve_();
  if (op_focei.maxInnerIterations <= 0){
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int id = 0; id < rx->nsub; id++){
      innerEval(id);
    }
  } else {
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif    
    for (int id = 0; id < rx->nsub; id++){
      innerOpt1(id, 0);
    }
  }
}

double foceiLik(double *theta){
  updateTheta(theta);
  innerOpt();
  double lik = 0.0;
  for (int id=rx->nsub; id--;){
    focei_ind fInd = (inds_focei[id]);
    lik += fInd.lik[0];
  }
  return lik;
}

void numericGrad(double *theta){
  rx = getRxSolve_();
#ifdef _OPENMP
  int cores = rx->op->cores;
#else
  int cores = 1;
#endif
  int npars = op_focei.npars;
  if (op_focei.estLambda) npars--;
  int cpar;
  for (cpar = npars; cpar--;){
    // Gradient needs can be parallelized for each parameter then gradient calculated.
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
    for (int gid=0; gid < rx->nsub*2; gid++){
      int likId, id;
      focei_ind fInd;
      likId = (gid % 2);
      if (likId){
        id = (gid-1)/2;
        fInd = (inds_focei[id]);
        fInd.thVal[likId]= theta[cpar]+op_focei.cEps*theta[cpar];
	updateTheta1(fInd.thVal[likId], cpar);
        // Upper
        innerOpt1(id, 2);
      } else {
        id = gid/2;
        fInd = (inds_focei[id]);
        fInd.thVal[likId]= theta[cpar]-op_focei.cEps*theta[cpar];
        updateTheta1(fInd.thVal[likId], cpar);
        // Lower
        innerOpt1(id, 1);
      }
    }
    // Now calculate individual gradient components
    for (int gid=rx->nsub; gid--;){
      focei_ind fInd = (inds_focei[gid]);
      fInd.thetaGrad[cpar] = (fInd.lik[2]-fInd.lik[1])/(fInd.thVal[1]-fInd.thVal[0]);
      // thetaGrad[cpar] += fInd.thetaGrad[cpar];
    }
    // Reset theta
    updateTheta1(theta[cpar], cpar);
  }
  // Calculate Overall gradient based on central differences
  std::fill_n(&op_focei.thetaGrad[0], npars, 0.0);
  for (cpar = npars; cpar--;){
    for (int gid=0; gid < rx->nsub; gid++){
      focei_ind fInd = (inds_focei[gid]);
      op_focei.thetaGrad[cpar]+=fInd.thetaGrad[cpar];
    }
  }
  // Calculate exact gradient for Cox-Box Yeo-Johnson
  if (op_focei.estLambda){
    for (int id=0; id < rx->nsub; id++){
      rx_solving_options_ind ind = rx->subjects[id];
      for (unsigned int j = ind.n_all_times; j--;){
	if (!ind.evid[j]){
	  op_focei.thetaGrad[npars] += tbsDL(ind.dv[j]);
	}
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Setup FOCEi functions
CharacterVector rxParams_(const RObject &obj);
List rxModelVars_(const RObject &obj);

static inline void foceiSetupTrans_(CharacterVector pars){
  unsigned int k, j,  ps = pars.size();
  k=ps;
  std::string thetaS;
  std::string etaS;
  std::string cur;
  foceiEtaN(ps+1);
  foceiThetaN(ps+1);
  op_focei.neta = 0;
  op_focei.ntheta = 0;
  for (;k--;){
    for (j = ps; j--;){
      // Compare ETAS first since they are smaller strings.
      cur = as<std::string>(pars[k]);
      etaS = "ETA[" + std::to_string(j+1) + "]";
      if (cur == etaS){
        op_focei.etaTrans[j] = k;
        op_focei.neta++;
        break;
      } else {
        thetaS = "THETA[" + std::to_string(j+1) + "]";
        if (cur == thetaS){
          op_focei.thetaTrans[j] = k;
          op_focei.ntheta++;
          break;
        }
      }
    }
  }
  op_focei.nzm = op_focei.neta * (op_focei.neta + 13) / 2;
}

static inline void foceiSetupTheta_(RObject &obj,
		      NumericVector theta,
                      Nullable<LogicalVector> thetaFixed, 
                      double lambda,
		      bool estLambda,
                      double scaleTo){
  // Get the fixed thetas
  // fixedTrans gives the theta->full theta translation
  // initPar is the initial parameters used for parameter scaling.
  op_focei.scaleTo=scaleTo;
  int thetan = theta.size();
  int omegan = getOmegaN();
  NumericVector omegaTheta = getOmegaTheta();
  int fixedn = 0;
  int j;
  LogicalVector thetaFixed2;
  if (!thetaFixed.isNull()){
    thetaFixed2 = LogicalVector(thetaFixed);
    for (j = thetaFixed2.size(); j--;){
      if (thetaFixed2[j]) fixedn++;
    }
  } else {
    thetaFixed2 =LogicalVector(0);
  }
  int npars = thetan+omegan-fixedn;
  List mvi = rxModelVars_(obj);
  rxUpdateInnerFuns(as<SEXP>(mvi["trans"]));
  
  if (estLambda){
    npars++;
    foceiThetaN(npars);
    foceiSetupTrans_(as<CharacterVector>(mvi["params"]));
    op_focei.fullTheta[npars-1]=lambda;
    op_focei.estLambda = 1;
  } else {
    foceiSetupTrans_(as<CharacterVector>(mvi["params"]));
    foceiThetaN(npars);
    op_focei.estLambda = 0;
  }
  std::copy(theta.begin(), theta.end(), &op_focei.fullTheta[0]);  
  std::copy(omegaTheta.begin(), omegaTheta.end(), &op_focei.fullTheta[0]+thetan);
  op_focei.npars = npars;
  op_focei.thetan = thetan;
  op_focei.omegan = omegan;
  int k = 0;
  for (j = 0; j < npars+fixedn; j++){
    if (j < thetaFixed2.size() && !thetaFixed2[j]){
      if (j < theta.size()){
        op_focei.initPar[k] = theta[j];
      } else if (j < theta.size() + omegan){
        op_focei.initPar[k] = omegaTheta[j-theta.size()];
      } else {
        op_focei.initPar[k] = lambda;
      }
      op_focei.fixedTrans[k++] = j;
    } else if (j >= thetaFixed2.size()){
      if (j < theta.size()){
        op_focei.initPar[k] = theta[j];
      } else if (j < theta.size() + omegan){
        op_focei.initPar[k] = omegaTheta[j-theta.size()];
      } else {
        op_focei.initPar[k] = lambda;
      }
      op_focei.fixedTrans[k++] = j;
    }
  }
  std::fill(&op_focei.theta[0], &op_focei.theta[0] + op_focei.ntheta, op_focei.scaleTo);
}

static inline void foceiSetupEta_(Nullable<NumericMatrix> etaMat = R_NilValue){
  rx = getRxSolve_();
  NumericMatrix etaMat0;
  if (etaMat.isNull()){
    etaMat0 = NumericMatrix(rx->nsub, op_focei.neta, 0);
  } else {
    etaMat0 = as<NumericMatrix>(etaMat);
    if ((int)(etaMat0.nrow()) != rx->nsub || (int)(etaMat0.ncol()) != op_focei.neta){
      stop("The Eta matrix needs %d rows x %d columns.", 
	   rx->nsub, op_focei.neta);
    }
  }
  etaMat0 = transpose(etaMat0);
  foceiGEtaN((op_focei.neta+1)*rx->nsub);
  foceiGThetaN(op_focei.npars*(rx->nsub + 1));
  foceiGgZm(op_focei.neta*(op_focei.neta+13)/2*rx->nsub);
  unsigned int i, j = 0, k = 0, ii=0, jj = 0;
  rxFoceiEnsure(rx->nsub);
  focei_ind fInd;
  for (i = rx->nsub; i--;){
    fInd = (inds_focei[i]);
    fInd.eta = &op_focei.geta[0]+j;
    fInd.eta[j] = i;
    // Copy in etaMat0 to the inital eta stored (0 if unspecified)
    std::copy(&fInd.eta[0]+1, &fInd.eta[0]+1+op_focei.neta, 
	      &etaMat0[0]+i*op_focei.neta);
    fInd.oldEta = &op_focei.goldEta[0]+k;
    std::copy(&fInd.oldEta[0], &fInd.oldEta[0]+op_focei.neta, 
              &etaMat0[0]+i*op_focei.neta);
    fInd.saveEta = &op_focei.gsaveEta[0]+k;
    std::copy(&fInd.saveEta[0], &fInd.saveEta[0]+op_focei.neta, 
              &etaMat0[0]+i*op_focei.neta);
    fInd.g = &op_focei.gG[0]+k;
    fInd.var = &op_focei.gVar[0]+k;
    fInd.zm = &op_focei.gZm[0] + ii;
    j+=op_focei.neta+1;
    k+=op_focei.neta;
    ii+=op_focei.neta * (op_focei.neta + 13) / 2;
    fInd.thetaGrad = &op_focei.gthetaGrad[0]+jj;
    jj+= op_focei.npars;
    fInd.mode = 1;
    fInd.uzm = 1;
  }
  op_focei.thetaGrad = &op_focei.gthetaGrad[0]+jj;
}

// [[Rcpp::export]]
RObject foceiSetup_(RObject &obj,
                    NumericVector theta,
		    Nullable<LogicalVector> thetaFixed = R_NilValue,
		    RObject rxInv = R_NilValue,
		    Nullable<NumericVector> lower  = R_NilValue,
                    Nullable<NumericVector> upper  = R_NilValue,
                    Nullable<NumericMatrix> etaMat = R_NilValue,
                    double scaleTo = 1,
		    double lambda = 1.0,
                    bool yjTrans = false,
                    bool estLambda = false,
		    Nullable<NumericVector> epsilon = R_NilValue,
		    unsigned int maxInnerEvals = 100,
		    unsigned int maxOuterEvals = 10000,
		    Nullable<IntegerVector> nsim = R_NilValue,
		    bool printInner = false){
  if (!rxIs(rxInv, "rxSymInvCholEnv")){
    stop("Omega isn't in the proper format.");
  } else {
    _rxInv = as<List>(_rxInv);
  }
  if (yjTrans){
    op_focei.yj=1;
  } else {
    op_focei.yj=0;
  }
  rx = getRxSolve_();
  foceiSetupTheta_(obj, theta, thetaFixed,  lambda, estLambda, scaleTo);
  foceiSetupEta_(etaMat);
  
  if (epsilon.isNull()){
    op_focei.epsilon=DOUBLE_EPS;
  } else {
    NumericVector epsilon2 = NumericVector(epsilon);
    if (epsilon2.size() == 1 && epsilon2[0] > 0){
      op_focei.epsilon=epsilon2[0];
    } else {
      stop("epsilon has to be a positive number.");
    }
  }
  op_focei.maxOuterIterations = maxOuterEvals;
  op_focei.maxInnerIterations = maxInnerEvals;
  if (nsim.isNull()){
    op_focei.nsim=maxInnerEvals*10;
  } else {
    IntegerVector nsim2 = IntegerVector(nsim);
    if (nsim2.size() == 1 && nsim2[0] > 1){
      op_focei.nsim=nsim2[0];
    } else {
      stop("nsim has to be a positive integer.");
    }
  }
  if (printInner){
    op_focei.imp=1;
  } else {
    op_focei.imp=0;
  }
  op_focei.cEps=sqrt(DOUBLE_EPS);
  return as<RObject>(R_NilValue);
}
