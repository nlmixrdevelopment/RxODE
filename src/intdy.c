#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <math.h>
#include "lsoda.h"
#include "common.h"
#include "lsoda_internal.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

void RSprintf(const char *format, ...);

int intdy(struct lsoda_context_t * ctx, double t, int k, double *dky)

/*
   Intdy computes interpolated values of the k-th derivative of the
   dependent variable vector y, and stores it in dky.  This routine
   is called within the package with k = 0 and *t = tout, but may
   also be called by the user for any k up to the current order.
   ( See detailed instructions in the usage documentation. )

   The computed values in dky are gotten by interpolation using the
   Nordsieck history array _C(yh).  This array corresponds uniquely to a
   vector-valued polynomial of degree nqcur or less, and dky is set
   to the k-th derivative of this polynomial at t.
   The formula for dky is

             q
   dky[i] = sum c[k][j] * ( t - _C(tn) )^(j-k) * _C(h)^(-j) * _C(yh)[j+1][i]
            j=k

   where c[k][j] = j*(j-1)*...*(j-k+1), q = nqcur, _C(tn) = tcur, _C(h) = hcur.
   The quantities _C(nq) = nqcur, l = _C(nq)+1, neq = neq, _C(tn), and _C(h) are declared
   static globally.  The above sum is done in reverse order.
   *iflag is returned negative if either k or t is out of bounds.
*/

{
	int             i, ic, j, jj, jp1;
	double          c, r, s, tp;

	const int neq = ctx->neq;
	if (k < 0 || k > _C(nq)) {
	  RSprintf(_("[intdy] k = %d illegal\n"), k);
	  return -1;
	}
	tp = _C(tn) - _C(hu) - 100. * ETA * (_C(tn) + _C(hu));
	if ((t - tp) * (t - _C(tn)) > 0.) {
	  RSprintf(_("intdy -- t = %g illegal. t not in interval tcur - _C(hu) to tcur\n"), t);
	  return -2;
	}
	s = (t - _C(tn)) / _C(h);
	ic = 1;
	for (jj = (_C(nq) + 1) - k; jj <= _C(nq); jj++)
		ic *= jj;
	c = (double) ic;
	for (i = 1; i <= neq; i++)
		dky[i] = c * _C(yh)[_C(nq) + 1][i];
	for (j = _C(nq) - 1; j >= k; j--) {
		jp1 = j + 1;
		ic = 1;
		for (jj = jp1 - k; jj <= j; jj++)
			ic *= jj;
		c = (double) ic;
		for (i = 1; i <= neq; i++)
			dky[i] = c * _C(yh)[jp1][i] + s * dky[i];
	}
	if (k == 0)
		return 0;
	r = pow(_C(h), (double) (-k));
	for (i = 1; i <= neq; i++)
		dky[i] *= r;
	return 0;
}				/* end intdy   */

