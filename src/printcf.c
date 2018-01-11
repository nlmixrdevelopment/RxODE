#include <R.h>
#include <Rinternals.h>

#include <stdio.h>
#include "common.h"
#include "lsoda.h"
int main(void) {
	int i, j;
	cfode(1);
	Rprintf("static double tesco1[13][4] = {\n");

	for(i = 0; i < 13; i++) {
		Rprintf("{ ");
		for(j = 0; j < 4; j++) {
			Rprintf("%a, ", _C(tesco)[i][j]);
		}
		Rprintf("}, \n");
	}
	Rprintf("} ;\n");

	Rprintf("static double elco1[13][14] = {\n");

	for(i = 0; i < 13; i++) {
		Rprintf("{ ");
		for(j = 0; j < 14; j++) {
			Rprintf("%a, ", _C(elco)[i][j]);
			if((j + 1) % 4 == 0) Rprintf("\n  ");
		}
		Rprintf("}, \n");
	}
	Rprintf("} ;\n");

	cfode(2);

	Rprintf("static double tesco2[13][4] = {\n");

	for(i = 0; i < 13; i++) {
		Rprintf("{ ");
		for(j = 0; j < 4; j++) {
			Rprintf("%a, ", _C(tesco)[i][j]);
		}
		Rprintf("}, \n");
	}
	Rprintf("} ;\n");

	Rprintf("static double elco2[13][14] = {\n");

	for(i = 0; i < 13; i++) {
		Rprintf("{ ");
		for(j = 0; j < 14; j++) {
			Rprintf("%a, ", _C(elco)[i][j]);
			if((j + 1) % 4 == 0) Rprintf("\n  ");
		}
		Rprintf("}, \n");
	}
	Rprintf("} ;\n");
	Rprintf("static double cm1[13] = {\n");
	for(i = 0; i < 13; i++) {
		Rprintf("%a, ", _C(tesco1)[i][2] *_C(elco1)[i][i + 1]);
		if((i + 1) % 4 == 0) Rprintf("\n  ");
		Rprintf("}, \n");
	}
	Rprintf("};\n");
	Rprintf("static double cm2[13] = {\n");
	for(i = 0; i < 13; i++) {
		Rprintf("%a, ", _C(tesco2)[i][2] *_C(elco2)[i][i + 1]);
		if((i + 1) % 4 == 0) Rprintf("\n  ");
		Rprintf("}, \n");
	}
	Rprintf("};\n");
}
