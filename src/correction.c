
#include <math.h>
#include "lsoda.h"
#include "common.h"
#include "blas.h"
#include "lsoda_internal.h"

/* m is the correction count */
int correction(struct lsoda_context_t * ctx, double *y, double pnorm, double *del, double *delp, double told, int *m)
/*
   *corflag = 0 : corrector converged,
              1 : step size to be reduced, redo prediction,
              2 : corrector cannot converge, failure flag.
*/

{
	int             i;
	double          rm, rate, dcon;
	const int neq = ctx->neq;
/*
   Up to maxcor corrector iterations are taken.  A convergence test is
   made on the r.m.s. norm of each correction, weighted by the error
   weight vector _C(ewt).  The sum of the corrections is accumulated in the
   vector _C(acor)[i].  The _C(yh) array is not altered in the corrector loop.
*/

	*m = 0;
	rate = 0.;
	*del = 0.;
	for (i = 1; i <= neq; i++)
		y[i] = _C(yh)[1][i];
	(*ctx->function) (_C(tn), y + 1, _C(savf) + 1, ctx->data);
	_C(nfe)++;
/*
   If indicated, the matrix P = I - _C(h) * _C(el)[1] * J is reevaluated and
   preprocessed before starting the corrector iteration.  _C(ipup) is set
   to 0 as an indicator that this has been done.
*/
	while (1) {
		if (*m == 0) {
			if (_C(ipup) > 0) {
				int ierpj = prja(ctx, y);
				_C(jcur) = 1;
				_C(ipup) = 0;
				_C(rc) = 1.;
				_C(nslp) = _C(nst);
				_C(crate) = 0.7;
				if (!ierpj) {
					return corfailure(ctx, told);
				}
			}
			for (i = 1; i <= neq; i++)
				_C(acor)[i] = 0.;
		}		/* end if ( *m == 0 )   */
		if (_C(miter) == 0) {
/*
   In case of functional iteration, update y directly from
   the result of the last function evaluation.
*/
			for (i = 1; i <= neq; i++) {
				_C(savf)[i] = _C(h) * _C(savf)[i] - _C(yh)[2][i];
				y[i] = _C(savf)[i] - _C(acor)[i];
			}
			*del = vmnorm(neq, y, _C(ewt));
			for (i = 1; i <= neq; i++) {
				y[i] = _C(yh)[1][i] + _C(el)[1] * _C(savf)[i];
				_C(acor)[i] = _C(savf)[i];
			}
		}
		/* end functional iteration   */
		/*
		   In the case of the chord method, compute the corrector error,
		   and solve the linear system with that as right-hand side and
		   P as coefficient matrix.
		 */ 
		else {
			for (i = 1; i <= neq; i++)
				y[i] = _C(h) * _C(savf)[i] - (_C(yh)[2][i] + _C(acor)[i]);
			solsy(ctx, y);
			*del = vmnorm(neq, y, _C(ewt));
			for (i = 1; i <= neq; i++) {
				_C(acor)[i] += y[i];
				y[i] = _C(yh)[1][i] + _C(el)[1] * _C(acor)[i];
			}
		}		/* end chord method   */
/*
   Test for convergence.  If *m > 0, an estimate of the convergence
   rate constant is stored in _C(crate), and this is used in the test.

   We first check for a change of iterates that is the size of
   roundoff error.  If this occurs, the iteration has converged, and a
   new rate estimate is not formed.
   In all other cases, force at least two iterations to estimate a
   local Lipschitz constant estimate for Adams method.
   On convergence, form _C(pdest) = local maximum Lipschitz constant
   estimate.  _C(pdlast) is the most recent nonzero estimate.
*/
		if (*del <= 100. * pnorm * ETA)
			break;
		if (*m != 0 || _C(meth) != 1) {
			if (*m != 0) {
				rm = 1024.0;
				if (*del <= (1024. * *delp))
					rm = *del / *delp;
				rate = fmax(rate, rm);
				_C(crate) = fmax(0.2 * _C(crate), rm);
			}
			double conit = 0.5 / (double) (_C(nq) + 2);
			dcon = *del * fmin(1., 1.5 * _C(crate)) / (_C(tesco)[_C(nq)][2] * conit);
			if (dcon <= 1.) {
				_C(pdest) = fmax(_C(pdest), rate / fabs(_C(h) * _C(el)[1]));
				if (_C(pdest) != 0.)
					_C(pdlast) = _C(pdest);
				break;
			}
		}
/*
   The corrector iteration failed to converge.
   If _C(miter) != 0 and the Jacobian is out of date, prja is called for
   the next try.   Otherwise the _C(yh) array is retracted to its values
   before prediction, and _C(h) is reduced, if possible.  If _C(h) cannot be
   reduced or mxncf failures have occured, exit with corflag = 2.
*/
		(*m)++;
		if (*m == MAXCOR || (*m >= 2 && *del > 2. * *delp)) {
			if (_C(miter) == 0 || _C(jcur) == 1) {
				return corfailure(ctx, told);
			}
			_C(ipup) = _C(miter);
/*
   Restart corrector if Jacobian is recomputed.
*/
			*m = 0;
			rate = 0.;
			*del = 0.;
			for (i = 1; i <= neq; i++)
				y[i] = _C(yh)[1][i];
			(*ctx->function) (_C(tn), y + 1, _C(savf) + 1, ctx->data);
			_C(nfe)++;
		}
/*
   Iterate corrector.
*/
		else {
			*delp = *del;
			(*ctx->function) (_C(tn), y + 1, _C(savf) + 1, ctx->data);
			_C(nfe)++;
		}
	}			/* end while   */
	return 0;
}				/* end correction   */

