#include "lsoda.h"
#include "common.h"
#include "lsoda_internal.h"
#include <math.h>
void scaleh(struct lsoda_context_t * ctx, double rh)
{
	double          r;
	int             j, i;
	const int neq = ctx->neq;
	const double hmxi = ctx->opt->hmxi;
/*
   If _C(h) is being changed, the _C(h) ratio rh is checked against _C(rmax), hmin,
   and hmxi, and the _C(yh) array is rescaled.  _C(ialth) is set to (_C(nq) + 1) = _C(nq) + 1
   to prevent a change of _C(h) for that many steps, unless forced by a
   convergence or error test failure.
*/
	rh = fmin(rh, _C(rmax));
	rh = rh / fmax(1., fabs(_C(h)) * hmxi * rh);
/*
   If _C(meth) = 1, also restrict the new step size by the stability region.
   If this reduces _C(h), set _C(irflag) to 1 so that if there are roundoff
   problems later, we can assume that is the cause of the trouble.
*/
	if (_C(meth) == 1) {
		_C(irflag) = 0;
		double pdh = fmax(fabs(_C(h)) * _C(pdlast), 0.000001);
		if ((rh * pdh * 1.00001) >= sm1[_C(nq)]) {
			rh = sm1[_C(nq)] / pdh;
			_C(irflag) = 1;
		}
	}
	r = 1.;
	for (j = 2; j <= (_C(nq) + 1); j++) {
		r *= rh;
		for (i = 1; i <= neq; i++)
			_C(yh)[j][i] *= r;
	}
	_C(h) *= rh;
	_C(rc) *= rh;
	_C(ialth) = (_C(nq) + 1);

}				/* end scaleh   */
