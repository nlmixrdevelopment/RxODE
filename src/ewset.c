#error this file is deprecated ewset is a marco in lsoda.c
#if 0

#include "lsoda.h"
#include <math.h>
#include <stdio.h>
int ewset(struct lsoda_context_t * ctx, const int neq, double * ewt, const double *rtol, const double *atol, const double *ycur)
{
	int             i;

	for (i = 1; i <= neq; i++)
		ewt[i] = rtol[i] * fabs(ycur[i]) + atol[i];

	for (i = 1; i <= neq; i++) {
		ewt[i] = 1. / ewt[i];
	}
	return 1;
}				/* end ewset   */
#endif
