#include <R.h>
#include <Rinternals.h>
#include "lsoda.h"
#include "common.h"
#include "cfode_static.h"

void cfode_static (struct lsoda_context_t * ctx, int meth)
{
#ifdef CFODE_STATIC
	if (meth == 1) {
		_C(elco) = elco1;
		_C(tesco) = tesco1;

	} else {
		_C(elco) = elco2;
		_C(tesco) = tesco2;
	}
#endif
}
#include <stdio.h>
void printcm12() {
	REprintf("static double cm1[13] = {\n");
	int i;
	for(i = 0; i < 13; i++) {
		REprintf("%a, ", (tesco1)[i][2] *(elco1)[i][i + 1]);
		if((i + 1) % 4 == 0) REprintf("\n  ");
	}
	REprintf("};\n");
	REprintf("static double cm2[13] = {\n");
	for(i = 0; i < 13; i++) {
		REprintf("%a, ", (tesco2)[i][2] *(elco2)[i][i + 1]);
		if((i + 1) % 4 == 0) REprintf("\n  ");
	}
	REprintf("};\n");
	
}
