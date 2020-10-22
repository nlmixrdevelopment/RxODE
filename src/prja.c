#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <math.h>
#include "lsoda.h"
#include "blas.h"
#include "common.h"
#include "lsoda_internal.h"

int prja(struct lsoda_context_t * ctx, double *y)
{
	int             i, ier, j;
	double          fac, hl0, r, r0, yj;
	const int neq = ctx->neq;
/*
   prja is called by stoda to compute and process the matrix
   P = I - _C(h) * _C(el)[1] * J, where J is an approximation to the Jacobian.
   Here J is computed by finite differencin.
   J, scaled by -_C(h) * _C(el)[1], is stored in _C(wm).  Then the norm of J ( the
   matrix norm consistent with the weihted max-norm on vectors given
   by vmnorm ) is computed, and J is overwritten by P.  P is then
   subjected to LU decomposition in preparation for later solution
   of linear systems with p as coefficient matrix.  This is done
   by defa if _C(miter) = 2, and by dgbfa if _C(miter) = 5.
*/
	_C(nje)++;
	hl0 = _C(h) * _C(el)[1];
/*
   If _C(miter) = 2, make neq calls to f to approximate J.
*/
	if (_C(miter) != 2) {
		REprintf("[prja] _C(miter) != 2\n");
		return 0;
	}
	if (_C(miter) == 2) {
		fac = vmnorm0(neq, _C(savf), _C(ewt));
		r0 = 1000. * fabs(_C(h)) * ETA * ((double) neq) * fac;
		if (r0 == 0.)
			r0 = 1.;
		for (j = 1; j <= neq; j++) {
			yj = y[j];
			r = fmax(SQRTETA * fabs(yj), r0 / _C(ewt)[j]);
			y[j] += r;
			fac = -hl0 / r;
			(*ctx->function) (_C(tn), y + 1, _C(acor) + 1, ctx->data);
			for (i = 1; i <= neq; i++)
				_C(wm)[i][j] = (_C(acor)[i] - _C(savf)[i]) * fac;
			y[j] = yj;
		}
		_C(nfe) += neq;
/*
   Compute norm of Jacobian.
*/
		_C(pdnorm) = fnorm0(neq, _C(wm), _C(ewt)) / fabs(hl0);
/*
   Add identity matrix.
*/
		for (i = 1; i <= neq; i++)
			_C(wm)[i][i] += 1.;
/*
   Do LU decomposition on P.
*/
		dgefa0(_C(wm), neq, _C(ipvt), &ier);
		if (ier != 0)
		return 0;
	}
	return 1;
}				/* end prja   */

