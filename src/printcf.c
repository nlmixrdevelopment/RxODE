
#include <stdio.h>
#include "common.h"
#include "lsoda.h"
int main(void) {
	int i, j;
	cfode(1);
	printf("static double tesco1[13][4] = {\n");

	for(i = 0; i < 13; i++) {
		printf("{ ");
		for(j = 0; j < 4; j++) {
			printf("%a, ", _C(tesco)[i][j]);
		}
		printf("}, \n");
	}
	printf("} ;\n");

	printf("static double elco1[13][14] = {\n");

	for(i = 0; i < 13; i++) {
		printf("{ ");
		for(j = 0; j < 14; j++) {
			printf("%a, ", _C(elco)[i][j]);
			if((j + 1) % 4 == 0) printf("\n  ");
		}
		printf("}, \n");
	}
	printf("} ;\n");

	cfode(2);

	printf("static double tesco2[13][4] = {\n");

	for(i = 0; i < 13; i++) {
		printf("{ ");
		for(j = 0; j < 4; j++) {
			printf("%a, ", _C(tesco)[i][j]);
		}
		printf("}, \n");
	}
	printf("} ;\n");

	printf("static double elco2[13][14] = {\n");

	for(i = 0; i < 13; i++) {
		printf("{ ");
		for(j = 0; j < 14; j++) {
			printf("%a, ", _C(elco)[i][j]);
			if((j + 1) % 4 == 0) printf("\n  ");
		}
		printf("}, \n");
	}
	printf("} ;\n");
	printf("static double cm1[13] = {\n");
	for(i = 0; i < 13; i++) {
		printf("%a, ", _C(tesco1)[i][2] *_C(elco1)[i][i + 1]);
		if((i + 1) % 4 == 0) printf("\n  ");
		printf("}, \n");
	}
	printf("};\n");
	printf("static double cm2[13] = {\n");
	for(i = 0; i < 13; i++) {
		printf("%a, ", _C(tesco2)[i][2] *_C(elco2)[i][i + 1]);
		if((i + 1) % 4 == 0) printf("\n  ");
		printf("}, \n");
	}
	printf("};\n");
}
