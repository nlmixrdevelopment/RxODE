#include <RcppArmadillo.h>
#define NETAs 20
#define NTHETAs 20
#define innerOde(id) ind_solve(rx, id, inner_dydt_liblsoda, inner_dydt_lsoda_dum, inner_jdum_lsoda, inner_dydt, inner_update_inis, inner_global_jt)
#define getOmegaInv() (as<arma::mat>(rxSymInvCholEnvCalculate(_rxInv, "omegaInv", R_NilValue)))
#define getOmegaDet() (as<arma::mat>(rxSymInvCholEnvCalculate(_rxInv, "log.det.OMGAinv.5", R_NilValue)))
#define setOmegaTheta(x) rxSymInvCholEnvCalculate(_rxInv, "theta", x)

using namespace Rcpp;
using namespace arma;
extern "C"{
#include "solve.h"
  void ind_solve(rx_solve *rx, unsigned int cid, t_dydt_liblsoda dydt_lls, 
		 t_dydt_lsoda_dum dydt_lsoda, t_jdum_lsoda jdum,
                 t_dydt c_dydt, t_update_inis u_inis, int jt);
}

RObject rxSymInvCholEnvCalculate(List obj, std::string what, Nullable<NumericVector> theta = R_NilValue);
bool rxIs(const RObject &obj, std::string cls);

List _rxInv;

// These are focei inner options
typedef struct {
  // Integer of ETAs
  unsigned int etaTransN;
  int *etaTrans;
  unsigned int neta;
  unsigned int ntheta;
  double *theta;
  unsigned int thetaTransN;
  int *thetaTrans;
  double scaleTo;
  double epsilon;
  unsigned int maxInnerIterations;
  unsigned int nsim;
  unsigned int nzm;
  unsigned int imp;
} focei_options;

focei_options op_focei;

extern "C" void rxOptionsIniFocei(){
  op_focei.etaTransN = NETAs;
  op_focei.etaTrans = Calloc(NETAs,int);
  op_focei.thetaTrans = Calloc(NTHETAs,int);
  op_focei.thetaTransN = NTHETAs;
  op_focei.theta = Calloc(NTHETAs,double);
}

void foceiThetaN(unsigned int n){
  if (op_focei.thetaTransN < n){
    unsigned int cur = op_focei.thetaTransN;
    while (cur < n){
      cur += NTHETAs;
    }
    Free(op_focei.thetaTrans);
    Free(op_focei.theta);
    op_focei.thetaTrans = Calloc(cur, int);
    op_focei.theta = Calloc(cur, double);
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
  }
}

extern "C" void rxOptionsFreeFocei(){
  Free(op_focei.etaTrans);
  Free(op_focei.thetaTrans);
  Free(op_focei.theta);
}

typedef struct {
  double *eta; // Eta includes the ID number for the patient
  // F and varaibility
  unsigned int nobs;
  unsigned int setup;
  
  double *oldEta;
  
  // Likilihood gradient
  double llik;
  mat lp;// = mat(neta,1);

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
    unsigned int n = op_focei.neta;
    mat L = eye(n,n);
    mat D = mat(n,n,fill::zeros);
    mat H = mat(n,n);
    unsigned int l_n = n*(n+1)/2;
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

CharacterVector rxParams_(const RObject &obj);
List rxModelVars_(const RObject &obj);

void foceiSetupTrans_(CharacterVector pars){
  unsigned int k, j,  ps = pars.size();
  k=ps;
  std::string thetaS;
  std::string etaS;
  std::string cur;
  foceiEtaN(ps);
  foceiThetaN(ps);
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

double likInner(double *eta){
  // id = eta[-1]
  // eta = eta
  double *eprev = &eta[0]-1;
  unsigned int id = (unsigned int)(eprev[0]);
  rx_solving_options_ind ind = rx->subjects[id];
  rx_solving_options *op = rx->op;
  focei_ind fInd = (inds_focei[id]);
  double *par_ptr = ind.par_ptr;
  unsigned int i, j;
  bool recalc = false;
  if (!fInd.setup){
    recalc=true;
    fInd.nobs = ind.n_all_times - ind.ndoses;
    fInd.lp  = arma::mat(op_focei.neta, 1);
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
    std::fill_n(fInd.lp.begin(), fInd.lp.end(), 0.0);
    fInd.llik=0.0;
    double f, err, r, B, fpm, a, rp, c;
    for (j = ind.n_all_times; j--;){
      if (!ind.evid[j]){
	// Observation; Calc LHS.
	inner_calc_lhs((int)id, ind.all_times[j], ind.solve + j * op->neq, ind.lhs);
        f = ind.lhs[0];
	// fInd.f(k, 0) = ind.lhs[0];
	err = f - ind.dv[k];
	// fInd.err(k, 0) = ind.lhs[0] - ind.dv[k]; // pred-dv
        r = ind.lhs[op_focei.neta + 1];
	// fInd.r(k, 0) = ind.lhs[op_focei.neta+1];
        B = 2.0/r;
	// fInd.B(k, 0) = 2.0/ind.lhs[op_focei.neta+1];
	// lhs 0 = F
	// lhs 1-eta = df/deta
	// FIXME faster initiliaitzation via copy or elim
	for (i = op_focei.neta; i--; ){
	  fpm = a = ind.lhs[i + 1]; // Almquist uses different a (see eq #15)
	  rp  = ind.lhs[i + op_focei.neta + 2];
	  c   = rp/r;
	  // lp is eq 12 in Almquist 2015
	  // // .5*apply(eps*fp*B + .5*eps^2*B*c - c, 2, sum) - OMGAinv %*% ETA
	  fInd.lp(i, 0)  += 0.25 * err * err * B * c - 0.5 * c - 0.5 * err * fpm * B;
	}
	// Eq #10
	fInd.llik += -0.5 * err * err/r - 0.5*log(r);
        k--;
      }
    }
    // Now finalize lp
    mat etam = arma::mat(op_focei.neta, 1);
    std::copy(&eta[0], &eta[0] + op_focei.neta, etam.begin()); // fill in etam
    mat omegaInv = getOmegaInv();
    // Finalize eq. #12
    fInd.lp = -(fInd.lp - omegaInv * etam);
    // Partially finalize #10
    mat llikm = -(fInd.llik - 0.5*(etam.t() * omegaInv * etam));
    fInd.llik = llikm(0,0);
    std::copy(&eta[0], &eta[0] + op_focei.neta, &fInd.oldEta[0]);
  }
  return fInd.llik;
}

double *lpInner(double *eta){
  double *eprev = &eta[0]-1;
  unsigned int id = (unsigned int)(eprev[0]);
  focei_ind fInd = (inds_focei[id]);
  likInner(eta);
  return &fInd.oldEta[0];
}

// [[Rcpp::export]]
RObject foceiSetup_(RObject &obj,
		    RObject data,
                    NumericVector theta,
		    RObject rxInv,
		    Nullable<NumericVector> epsilon = R_NilValue,
		    unsigned int maxInnerEvals = 100,
		    Nullable<IntegerVector> nsim = R_NilValue,
		    bool printInner = false){
  if (!rxIs(rxInv, "rxSymInvCholEnv")){
    stop("Omega isn't in the proper format.");
  } else {
    _rxInv = as<List>(_rxInv);
  }
  // if (thetaFixed.size() ==)
  rx = getRxSolve_();
  
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
  List mvi = rxModelVars_(obj);
  rxUpdateInnerFuns(as<SEXP>(mvi["trans"]));
  foceiSetupTrans_(as<CharacterVector>(mvi["params"]));
  return as<RObject>(R_NilValue);
}


