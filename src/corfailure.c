#include "lsoda.h"
#include "common.h"
#include "lsoda_internal.h"
#include <math.h>

int corfailure(struct lsoda_context_t * ctx, double told)
{
	int             j, i1, i;
	const int neq = ctx->neq;
	const double hmin = ctx->opt->hmin;
	_C(ncf)++;
	_C(rmax) = 2.;
	_C(tn) = told;
	for (j = _C(nq); j >= 1; j--)
		for (i1 = j; i1 <= _C(nq); i1++) {
			for (i = 1; i <= neq; i++)
				_C(yh)[i1][i] -= _C(yh)[i1 + 1][i];
		}
	if (fabs(_C(h)) <= hmin * 1.00001 || _C(ncf) == MXNCF) {
		return 2;
	}
	_C(ipup) = _C(miter);
	return 1;
}

