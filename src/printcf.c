#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include "common.h"
#include "lsoda.h"
int main(void) {
	int i, j;
	cfode(1);
	REprintf("static double tesco1[13][4] = {\n");

	for(i = 0; i < 13; i++) {
		REprintf("{ ");
		for(j = 0; j < 4; j++) {
			REprintf("%a, ", _C(tesco)[i][j]);
		}
		REprintf("}, \n");
	}
	REprintf("} ;\n");

	REprintf("static double elco1[13][14] = {\n");

	for(i = 0; i < 13; i++) {
		REprintf("{ ");
		for(j = 0; j < 14; j++) {
			REprintf("%a, ", _C(elco)[i][j]);
			if((j + 1) % 4 == 0) REprintf("\n  ");
		}
		REprintf("}, \n");
	}
	REprintf("} ;\n");

	cfode(2);

	REprintf("static double tesco2[13][4] = {\n");

	for(i = 0; i < 13; i++) {
		REprintf("{ ");
		for(j = 0; j < 4; j++) {
			REprintf("%a, ", _C(tesco)[i][j]);
		}
		REprintf("}, \n");
	}
	REprintf("} ;\n");

	REprintf("static double elco2[13][14] = {\n");

	for(i = 0; i < 13; i++) {
		REprintf("{ ");
		for(j = 0; j < 14; j++) {
			REprintf("%a, ", _C(elco)[i][j]);
			if((j + 1) % 4 == 0) REprintf("\n  ");
		}
		REprintf("}, \n");
	}
	REprintf("} ;\n");
	REprintf("static double cm1[13] = {\n");
	for(i = 0; i < 13; i++) {
		REprintf("%a, ", _C(tesco1)[i][2] *_C(elco1)[i][i + 1]);
		if((i + 1) % 4 == 0) REprintf("\n  ");
		REprintf("}, \n");
	}
	REprintf("};\n");
	REprintf("static double cm2[13] = {\n");
	for(i = 0; i < 13; i++) {
		REprintf("%a, ", _C(tesco2)[i][2] *_C(elco2)[i][i + 1]);
		if((i + 1) % 4 == 0) REprintf("\n  ");
		REprintf("}, \n");
	}
	REprintf("};\n");
}
