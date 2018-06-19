#include "lsoda.h"
#include "common.h"
#include "lsoda_internal.h"

#include <math.h>
#include "blas.h"

/*
   This routine returns from stoda to lsoda.  Hence freevectors() is
   not executed.
*/
#define endstoda() \
{ \
	double          r; \
	int             i; \
 \
	r = 1. / _C(tesco)[_C(nqu)][2]; \
	for (i = 1; i <= neq; i++) \
		_C(acor)[i] *= r; \
	_C(hold) = _C(h); \
}

/*
   The _C(el) vector and related constants are reset
   whenever the order _C(nq) is changed, or at the start of the problem.
*/
#define resetcoeff() \
{ \
	int             i; \
 \
	double el0 = _C(el)[1]; \
	for (i = 1; i <= (_C(nq) + 1); i++) \
		_C(el)[i] = _C(elco)[_C(nq)][i]; \
	_C(rc) = _C(rc) * _C(el)[1] / el0; \
}


int stoda(struct lsoda_context_t * ctx, double *y, int jstart)
{
	int kflag;
	int             i, i1, j, m;
	double          del, delp, dsm, dup, exup, r, rh, told;
	double          pnorm;

	const double hmin = ctx->opt->hmin;
	const int mxords = ctx->opt->mxords;
	const int mxordn = ctx->opt->mxordn;
	const int neq = ctx->neq;

/*
   stoda performs one step of the integration of an initial value
   problem for a system of ordinary differential equations.
   Note.. stoda is independent of the value of the iteration method
   indicator _C(miter), when this is != 0, and hence is independent
   of the type of chord method used, or the Jacobian structure.
   Communication with stoda is done with the following variables:

   jstart = an integer used for input only, with the following
            values and meanings:

               0  perform the first step,
             > 0  take a new step continuing from the last,
              -1  take the next step with a new value of _C(h),
                  neq, _C(meth), _C(miter), and/or matrix parameters.
              -2  take the next step with a new value of _C(h),
                  but with other inputs unchanged.

   kflag = a completion code with the following meanings:

             0  the step was successful,
            -1  the requested error could not be achieved,
            -2  corrector convergence could not be achieved,
            -3  fatal error in prja or solsy.

   _C(miter) = corrector iteration method:

             0  functional iteration,
            >0  a chord method corresponding to jacobian type jt.

*/
	kflag = 0;
	told = _C(tn);
	_C(ncf) = 0;
	delp = 0.;

/*
   On the first call, the order is set to 1, and other variables are
   initialized.  _C(rmax) is the maximum ratio by which _C(h) can be increased
   in a single step.  It is initially 1.e4 to compensate for the small
   initial _C(h), but then is normally equal to 10.  If a filure occurs
   (in corrector convergence or error test), _C(rmax) is set at 2 for
   the next increase.
   cfode is called to get the needed coefficients for both methods.
*/
	int maxord = mxordn;
	if (_C(meth) == 2)
		maxord = mxords;

	if (jstart == 0) {
		_C(nq) = 1;
		_C(ialth) = 2;
		_C(rmax) = 10000.;
		_C(rc) = 0.;
		_C(crate) = 0.7;
		_C(hold) = _C(h);
		_C(nslp) = 0;
		_C(ipup) = _C(miter);
		_C(el)[1] = 1.0;
/*
   Initialize switching parameters.  _C(meth) = 1 is assumed initially.
*/
		_C(icount) = 20;
		_C(irflag) = 0;
		_C(pdest) = 0.;
		_C(pdlast) = 0.;

		cfode(ctx, 1);
		resetcoeff();
	}			/* end if ( jstart == 0 )   */
	/*
	   The following block handles preliminaries needed when jstart = -1.
	   _C(ipup) is set to _C(miter) to force a matrix update.
	   If an order increase is about to be considered ( _C(ialth) = 1 ),
	   _C(ialth) is reset to 2 to postpone consideration one more step.
	   If the caller has changed _C(meth), cfode is called to reset
	   the coefficients of the method.
	   If _C(h) is to be changed, _C(yh) must be rescaled.
	   If _C(h) or _C(meth) is being changed, _C(ialth) is reset to (_C(nq) + 1) = _C(nq) + 1
	   to prevent further changes in _C(h) for that many steps.
	*/
	if (jstart == -1) {
		_C(ipup) = _C(miter);
		if (_C(ialth) == 1)
			_C(ialth) = 2;
		if (_C(meth) != _C(mused)) {
			cfode(ctx, _C(meth));
			_C(ialth) = (_C(nq) + 1);
			resetcoeff();
		}
		if (_C(h) != _C(hold)) {
			rh = _C(h) / _C(hold);
			_C(h) = _C(hold);
			scaleh(ctx, rh);
		}
	}			/* if ( jstart == -1 )   */
	if (jstart == -2) {
		if (_C(h) != _C(hold)) {
			rh = _C(h) / _C(hold);
			_C(h) = _C(hold);
			scaleh(ctx, rh);
		}
	}			/* if ( jstart == -2 )   */
	/*
	   Prediction.
	   This section computes the predicted values by effectively
	   multiplying the _C(yh) array by the pascal triangle matrix.
	   _C(rc) is the ratio of new to old values of the coefficient _C(h) * _C(el)[1].
	   When _C(rc) differs from 1 by more than ccmax, _C(ipup) is set to _C(miter)
	   to force pjac to be called, if a jacobian is involved.
	   In any case, prja is called at least every msbp steps.
	*/
	dsm = 0.0;
	while (1) {
/*
   Before the corrector starts.  _C(jcur) is set to 0
   to signal that the Jacobian involved may need updating later.
*/
		_C(jcur) = 0;
		while (1) {
			if (fabs(_C(rc) - 1.) > CCMAX)
				_C(ipup) = _C(miter);
			if (_C(nst) >= _C(nslp) + MSBP)
				_C(ipup) = _C(miter);
			_C(tn) += _C(h);
			for (j = _C(nq); j >= 1; j--)
				for (i1 = j; i1 <= _C(nq); i1++) {
					for (i = 1; i <= neq; i++)
						_C(yh)[i1][i] += _C(yh)[i1 + 1][i];
				}
			pnorm = vmnorm0(neq, _C(yh)[1], _C(ewt));

			int corflag = correction(ctx, y, pnorm, &del, &delp, told, &m);
			if (corflag == 0)
				break;
			if (corflag == 1) {
				rh = fmax(0.25, hmin / fabs(_C(h)));
				scaleh(ctx, rh);
				continue;
			}
			if (corflag == 2) {
				kflag = -2;
				_C(hold) = _C(h);
				jstart = 1;
				return kflag;
			}
		}		/* end inner while ( corrector loop )   */
/*
   The local error test is done now.
*/
		if (m == 0)
			dsm = del / _C(tesco)[_C(nq)][2];
		if (m > 0)
			dsm = vmnorm0(neq, _C(acor), _C(ewt)) / _C(tesco)[_C(nq)][2];
		if (dsm <= 1.) {
/*
   After a successful step, update the _C(yh) array.
   Decrease _C(icount) by 1, and if it is -1, consider switching methods.
   If a method switch is made, reset various parameters,
   rescale the _C(yh) array, and exit.  If there is no switch,
   consider changing _C(h) if _C(ialth) = 1.  Otherwise decrease _C(ialth) by 1.
   If _C(ialth) is then 1 and _C(nq) < maxord, then _C(acor) is saved for
   use in a possible order increase on the next step.
   If a change in _C(h) is considered, an increase or decrease in order
   by one is considered also.  A change in _C(h) is made only if it is by
   a factor of at least 1.1.  If not, _C(ialth) is set to 3 to prevent
   testing for that many steps.
*/
			kflag = 0;
			_C(nst)++;
			_C(hu) = _C(h);
			_C(nqu) = _C(nq);
			_C(mused) = _C(meth);
			for (j = 1; j <= (_C(nq) + 1); j++) {
				r = _C(el)[j];
				for (i = 1; i <= neq; i++)
					_C(yh)[j][i] += r * _C(acor)[i];
			}
			_C(icount)--;
			if (_C(icount) < 0) {
				methodswitch(ctx, dsm, pnorm, &rh);
				if (_C(meth) != _C(mused)) {
					rh = fmax(rh, hmin / fabs(_C(h)));
					scaleh(ctx, rh);
					_C(rmax) = 10.;
					endstoda();
					break;
				}
			}
/*
   No method switch is being made.  Do the usual step/order selection.
*/
			_C(ialth)--;
			if (_C(ialth) == 0) {
				double rhup = 0.;
				if ((_C(nq) + 1) != maxord + 1) {
					for (i = 1; i <= neq; i++)
						_C(savf)[i] = _C(acor)[i] - _C(yh)[maxord + 1][i];
					dup = vmnorm0(neq, _C(savf), _C(ewt)) / _C(tesco)[_C(nq)][3];
					exup = 1. / (double) ((_C(nq) + 1) + 1);
					rhup = 1. / (1.4 * pow(dup, exup) + 0.0000014);
				}
				int orderflag = orderswitch(ctx, rhup, dsm, &rh, kflag, maxord);
/*
   No change in _C(h) or _C(nq).
*/
				if (orderflag == 0) {
					endstoda();
					break;
				}
/*
   _C(h) is changed, but not _C(nq).
*/
				if (orderflag == 1) {
					rh = fmax(rh, hmin / fabs(_C(h)));
					scaleh(ctx, rh);
					_C(rmax) = 10.;
					endstoda();
					break;
				}
/*
   both _C(nq) and _C(h) are changed.
*/
				if (orderflag == 2) {
					resetcoeff();
					rh = fmax(rh, hmin / fabs(_C(h)));
					scaleh(ctx, rh);
					_C(rmax) = 10.;
					endstoda();
					break;
				}
			}	/* end if ( _C(ialth) == 0 )   */
			if (_C(ialth) > 1 || (_C(nq) + 1) == maxord + 1) {
				endstoda();
				break;
			}
			for (i = 1; i <= neq; i++)
				_C(yh)[maxord + 1][i] = _C(acor)[i];
			endstoda();
			break;
		}
		/* end if ( dsm <= 1. )   */
		/*
		   The error test failed.  kflag keeps track of multiple failures.
		   Restore _C(tn) and the _C(yh) array to their previous values, and prepare
		   to try the step again.  Compute the optimum step size for this or
		   one lower.  After 2 or more failures, _C(h) is forced to decrease
		   by a factor of 0.2 or less.
		 */
		else {
			kflag--;
			_C(tn) = told;
			for (j = _C(nq); j >= 1; j--)
				for (i1 = j; i1 <= _C(nq); i1++) {
					for (i = 1; i <= neq; i++)
						_C(yh)[i1][i] -= _C(yh)[i1 + 1][i];
				}
			_C(rmax) = 2.;
			if (fabs(_C(h)) <= hmin * 1.00001) {
				kflag = -1;
				_C(hold) = _C(h);
				jstart = 1;
				break;
			}
			if (kflag > -3) {
				int orderflag = orderswitch(ctx, 0., dsm, &rh, kflag, maxord);
				if (orderflag == 1 || orderflag == 0) {
					if (orderflag == 0)
						rh = fmin(rh, 0.2);
					rh = fmax(rh, hmin / fabs(_C(h)));
					scaleh(ctx, rh);
				}
				if (orderflag == 2) {
					resetcoeff();
					rh = fmax(rh, hmin / fabs(_C(h)));
					scaleh(ctx, rh);
				}
				continue;
			}
			/* if ( kflag > -3 )   */
			/*
			   Control reaches this section if 3 or more failures have occurred.
			   If 10 failures have occurred, exit with kflag = -1.
			   It is assumed that the derivatives that have accumulated in the
			   _C(yh) array have errors of the wrong order.  Hence the first
			   derivative is recomputed, and the order is set to 1.  Then
			   _C(h) is reduced by a factor of 10, and the step is retried,
			   until it succeeds or _C(h) reaches hmin.
			 */
			else {
				if (kflag == -10) {
					kflag = -1;
					_C(hold) = _C(h);
					jstart = 1;
					break;
				} else {
					rh = 0.1;
					rh = fmax(hmin / fabs(_C(h)), rh);
					_C(h) *= rh;
					for (i = 1; i <= neq; i++)
						y[i] = _C(yh)[1][i];
					(*ctx->function) (_C(tn), y + 1, _C(savf) + 1, ctx->data);
					_C(nfe)++;
					for (i = 1; i <= neq; i++)
						_C(yh)[2][i] = _C(h) * _C(savf)[i];
					_C(ipup) = _C(miter);
					_C(ialth) = 5;
					if (_C(nq) == 1)
						continue;
					_C(nq) = 1;
					resetcoeff();
					continue;
				}
			}	/* end else -- kflag <= -3 */
		}		/* end error failure handling   */
	}			/* end outer while   */

	return kflag;
}				/* end stoda   */
